// runner.c — 試行ループ(状態攪乱 → ウォームアップ → 計測)
#include "bp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t bp_pat_len = BP_PAT_LEN_DEFAULT;

void bp_print_header(void) {
    printf("mode,control,p1,p2,patlen,trial,branches,misses,cycles,"
           "instructions,fe_stalls,scaled,"
           "miss_per_patbranch,ns_per_patbranch,ipc\n");
}

// 1 系列ぶんの計測（攪乱 → ウォームアップ → 計測 → 出力）
static void measure_one(const bp_config *cfg, uint8_t *pat, uint8_t *noise,
                        int control, uint64_t trial) {
    bp_seed(cfg->seed * 1000003ULL + trial + 1);
    cfg->gen(pat, bp_pat_len, cfg->p1, cfg->p2, control);

    // 予測器状態の攪乱: 前の測定の学習を持ち越さない
    bp_gen_random(noise, bp_pat_len);
    volatile uint64_t s0 = cfg->kernel(noise, bp_pat_len, 4, cfg->ctx); (void)s0;

    // ウォームアップ(計測外): 当該系列を学習させる
    volatile uint64_t s1 = cfg->kernel(pat, bp_pat_len, cfg->warmup, cfg->ctx); (void)s1;

    bp_sample smp;
    bp_ctr_start();
    volatile uint64_t s2 = cfg->kernel(pat, bp_pat_len, cfg->reps, cfg->ctx); (void)s2;
    bp_ctr_stop(&smp);

    double patbr = (double)bp_pat_len * (double)cfg->reps;
    double mpb = (double)smp.misses / patbr;
    double npb = smp.ns / patbr;

    double ipc = smp.cycles ? (double)smp.instructions / (double)smp.cycles : 0.0;

    if (cfg->csv)
        printf("%s,%d,%g,%g,%zu,%llu,%llu,%llu,%llu,%llu,%llu,%d,%.6f,%.4f,%.4f\n",
               cfg->mode, control, cfg->p1, cfg->p2, bp_pat_len,
               (unsigned long long)trial, (unsigned long long)smp.branches,
               (unsigned long long)smp.misses, (unsigned long long)smp.cycles,
               (unsigned long long)smp.instructions,
               (unsigned long long)smp.fe_stalls,
               bp_ctr_scaled(), mpb, npb, ipc);
    else
        printf("%-10s ctl=%d p1=%-6g p2=%-6g trial=%llu  miss/br=%.5f  ns/br=%.4f\n",
               cfg->mode, control, cfg->p1, cfg->p2,
               (unsigned long long)trial, mpb, npb);
    fflush(stdout);
}

int bp_run(const bp_config *cfg) {
    size_t alloc = (bp_pat_len + 63) & ~(size_t)63;
    uint8_t *pat   = aligned_alloc(64, alloc);
    uint8_t *noise = aligned_alloc(64, alloc);
    if (!pat || !noise) { perror("alloc"); return 1; }

    // 第1試行の汚染対策。
    // 割り当て直後のバッファはページフォルトを伴い、キャッシュも冷えているため、
    // 第1試行だけ極端に遅くなる（実測で既定の warmup では吸収しきれず、
    // reps=4 のとき 58 ns/分岐 対 他試行 4 ns/分岐 という外れ値が出た）。
    // バッファを先に触ってページを確保し、計測外でカーネルを一度回しておく。
    memset(pat, 0, alloc);
    memset(noise, 0, alloc);
    bp_gen_random(noise, bp_pat_len);
    volatile uint64_t warm = cfg->kernel(noise, bp_pat_len, 2, cfg->ctx);
    (void)warm;

    for (uint64_t t = 0; t < cfg->trials; t++) {
        if (cfg->pair) {
            // test と control を隣接した時刻に測る。別プロセスで逐次実行すると
            // 周波数ドリフトが差分に系統誤差として乗るため、交互に測って打ち消す。
            measure_one(cfg, pat, noise, 0, t);
            measure_one(cfg, pat, noise, 1, t);
            continue;
        }
        measure_one(cfg, pat, noise, cfg->control, t);
    }
    free(pat); free(noise);
    return 0;
}
