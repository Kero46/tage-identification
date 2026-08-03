// TAGE テーブル構成の同定
//
// 観測は三系統。
//   dual : 二重相関(距離 d1, d2 の XOR)。両者が同一の履歴段で処理される場合に
//          競合が生じることを利用して段の区切りを定める。
//   ctx  : 文脈数 K を制限した相関系列。K を増やしてスラッシングが始まる点から
//          当該テーブルの実効エントリ数を推定する。
//   histd: 履歴長測定と同じ高密度相関。段の検出(距離の細かい掃引)に用いる。
//
// 注意: dual の競合はタグ衝突・インデックス衝突でも生じうる。ctx 掃引と
//       組み合わせて切り分けること(README の対照実験を参照)。
#include "bp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void g_dual (uint8_t *p, size_t n, double a, double b, int c) { bp_gen_dual(p, n, (unsigned)a, (unsigned)b, c); }
static void g_ctx  (uint8_t *p, size_t n, double a, double b, int c) { bp_gen_hist_ctx(p, n, (unsigned)a, (unsigned)b, c); }
static void g_histd(uint8_t *p, size_t n, double a, double b, int c) { (void)b; bp_gen_hist_dense(p, n, (unsigned)a, c); }
// aliasing 実験。S 個の分岐サイトに周期 p の固有パターンを与える。
// S を増やすと衝突が起こる「可能性が高まる」実験であり、観測された劣化を
// aliasing と断定するものではない（コード量も増えるため）。
static void g_alias(uint8_t *p, size_t n, double a, double b, int c) { bp_gen_alias(p, n, (unsigned)a, (unsigned)b, c); }

static void usage(const char *p) {
    fprintf(stderr,
      "使い方: %s --mode <dual|ctx|alias|histd> [オプション]\n"
      "  dual : --param D1  --param2 D2\n"
      "  ctx  : --param D   --param2 K   (K = 文脈数)\n"
      "  alias: --param S   --param2 p   (S = 分岐サイト数 2,4,8,16,32,64 / p = 周期)\n"
      "  histd: --param D   (奇数)\n"
      "  --control / --pair / --patlen / --frontend / --reps / --warmup / --trials / --seed / --csv\n", p);
}

int main(int argc, char **argv) {
    bp_config c;
    memset(&c, 0, sizeof(c));
    c.trials = 7; c.reps = 48; c.warmup = 8; c.seed = 1;
    const char *mode = NULL;
    int want_fe = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--mode")   && i+1 < argc) mode = argv[++i];
        else if (!strcmp(argv[i], "--param")  && i+1 < argc) c.p1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--param2") && i+1 < argc) c.p2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--control"))              c.control = 1;
        else if (!strcmp(argv[i], "--pair"))                 c.pair = 1;
        else if (!strcmp(argv[i], "--patlen")  && i+1 < argc) bp_pat_len = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--reps")   && i+1 < argc) c.reps   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--warmup") && i+1 < argc) c.warmup = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--trials") && i+1 < argc) c.trials = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--seed")   && i+1 < argc) c.seed   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--frontend"))             want_fe = 1;
        else if (!strcmp(argv[i], "--csv"))                  c.csv = 1;
        else { usage(argv[0]); return 2; }
    }
    if (!mode) { usage(argv[0]); return 2; }

    static unsigned nsites;
    if      (!strcmp(mode,"alias")) {
        nsites = (unsigned)c.p1;
        if (nsites != 2 && nsites != 4 && nsites != 8 &&
            nsites != 16 && nsites != 32 && nsites != 64) {
            fprintf(stderr, "エラー: alias の --param は 2,4,8,16,32,64 のいずれか\n");
            return 2;
        }
        c.gen = g_alias; c.kernel = bp_kernel_sites; c.ctx = &nsites;
    }
    else if (!strcmp(mode,"dual"))  c.gen = g_dual;
    else if (!strcmp(mode,"ctx"))   c.gen = g_ctx;
    else if (!strcmp(mode,"histd")) c.gen = g_histd;
    else { usage(argv[0]); return 2; }

    if (!strcmp(mode,"histd") && ((unsigned)c.p1 % 2) == 0) {
        fprintf(stderr, "エラー: histd は D が奇数でなければならない\n");
        return 2;
    }
    c.mode = mode;
    if (!c.kernel) c.kernel = bp_kernel_single;
    bp_ctr_init(want_fe);
    if (!bp_ctr_ok())
        fprintf(stderr, "警告: カウンタ取得不可。時間計測(ns/br)を使うこと。\n");
    if (c.csv) bp_print_header();
    return bp_run(&c);
}
