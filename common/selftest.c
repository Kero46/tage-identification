// selftest.c — パターン生成器の性質検査
//
// 測定結果を解釈する前に、生成器が意図した性質を持つことを確かめる。
// 特に重要なのが「必要履歴長の下限」の検査である。文脈数を制限した系列で、
// 直近の短いビット列から標的が決まってしまうと、狙った履歴長のテーブルでは
// なく短い履歴のテーブルを見ることになる（実際に一度この欠陥を作り込んだ）。
#include "bp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

static int fails = 0;
static void bad(const char *fmt, ...) __attribute__((format(printf, 1, 2)));
static void bad(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "  NG: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    fails++;
}

// ---------------------------------------------------------------- 必要履歴長
// 【検査手法の注意】「直近 h ビットが同じで標的が異なる例が無い」ことを
// 一貫性の根拠にしてはならない。有限サンプルでは h が大きくなると窓が
// 全て相異なるため、決定性ではなく一意性を見てしまい偽陽性になる。
// そこで h を小さく取り、サンプル数を 2^h より十分多く確保したうえで
// 「曖昧さが残ること」を積極的に確認する。
//
// 戻り値: 1 なら曖昧（h ビットでは標的が決まらない）、0 なら曖昧さが見つからない。
static int is_ambiguous_at(const uint8_t *pat, const size_t *tpos, size_t ntp,
                           unsigned h) {
    for (size_t a = 0; a < ntp; a++) {
        if (tpos[a] < h) continue;
        for (size_t b = a + 1; b < ntp; b++) {
            if (tpos[b] < h) continue;
            if (memcmp(&pat[tpos[a] - h], &pat[tpos[b] - h], h) == 0 &&
                pat[tpos[a]] != pat[tpos[b]])
                return 1;       // 同じ窓で標的が違う → h では決まらない
        }
    }
    return 0;
}

// 窓 h ビットで標的が一貫して決まるか（矛盾が無いか）。
// 「決まるべき」ことの確認に使う。偽陽性を避けるため、呼び出し側で
// サンプル数と h の関係に注意すること。
static int is_consistent_at(const uint8_t *pat, const size_t *tpos, size_t ntp,
                            unsigned h) {
    return !is_ambiguous_at(pat, tpos, ntp, h);
}

int bp_selftest(void) {
    size_t n = 1u << 13;
    uint8_t *pat = (uint8_t *)malloc(n);
    size_t  *tp  = (size_t *)malloc(sizeof(size_t) * n);
    if (!pat || !tp) { fprintf(stderr, "alloc 失敗\n"); return 1; }
    fails = 0;

    printf("== 低密度相関 hist: 標的が距離 D の値に一致 ==\n");
    bp_seed(12345);
    for (unsigned D = 2; D <= 64; D++) {
        bp_gen_hist(pat, n, D, 0);
        unsigned blk = D + 1;
        for (size_t i = 0; i + blk <= n; i += blk)
            if (pat[i + D] != pat[i]) { bad("hist D=%u で不一致 (i=%zu)", D, i); break; }
    }

    printf("== 高密度相関 histd: 標的が距離 D の値に一致、参照先は偶数位置 ==\n");
    for (unsigned D = 3; D <= 255; D += 2) {
        bp_gen_hist_dense(pat, n, D, 0);
        for (size_t i = (D | 1); i < n; i += 2) {
            if (pat[i] != pat[i - D]) { bad("histd D=%u で不一致 (i=%zu)", D, i); break; }
            if (((i - D) % 2) != 0)   { bad("histd D=%u で参照先が標的位置", D); break; }
        }
    }

    printf("== 高密度相関 histd: 近傍の履歴だけでは標的が決まらない ==\n");
    // h を小さく固定し、サンプル数を 2^h より十分多く取って曖昧さを確認する。
    {
        const unsigned h = 10;              // 2^10 = 1024 通り
        for (unsigned D = 13; D <= 41; D += 4) {
            bp_gen_hist_dense(pat, n, D, 0);
            size_t ntp = 0;
            for (size_t i = (D | 1); i < n && ntp < 4000; i += 2) tp[ntp++] = i;
            if (!is_ambiguous_at(pat, tp, ntp, h))
                bad("histd D=%u: 直近 %u ビットで標的が決まってしまう", D, h);
        }
    }

    printf("== 文脈数制限 ctx: 共通サフィックス長の履歴では決まらない ==\n");
    // 設計上、近傍 D-j ビットは全文脈で共通なので情報を持たない。
    // よって h = D-j では標的が曖昧でなければならない。
    // 逆に h = D なら識別ビットが全て入るので一貫して決まるべきである。
    {
        unsigned Ds[] = {15, 31, 63};
        unsigned Ks[] = {4, 16, 64};
        for (unsigned di = 0; di < 3; di++) {
            for (unsigned ki = 0; ki < 3; ki++) {
                unsigned D = Ds[di], K = Ks[ki];
                unsigned j = 0;
                while ((1u << j) < K) j++;
                unsigned suflen = D - j;

                bp_gen_hist_ctx(pat, n, D, K, 0);
                size_t ntp = 0;
                for (size_t i = D; i + 1 <= n && ntp < 2000; i += D + 1) tp[ntp++] = i;

                if (!is_ambiguous_at(pat, tp, ntp, suflen))
                    bad("ctx D=%u K=%u: 共通サフィックス %u ビットで標的が"
                        "決まってしまう（識別情報が近傍に漏れている）", D, K, suflen);
                if (!is_consistent_at(pat, tp, ntp, D))
                    bad("ctx D=%u K=%u: 履歴 %u ビットでも標的が決まらない", D, K, D);
            }
        }
    }

    printf("== 分岐間相関 cross: B が d ステップ前の A に一致 ==\n");
    // 偶数位置=A、奇数位置=B。B[i] = A[i-1-2d]
    for (unsigned d = 1; d <= 16; d++) {
        bp_gen_cross(pat, n, d, 0);
        for (size_t i = 1 + 2 * (size_t)d; i < n; i += 2)
            if (pat[i] != pat[(i - 1) - 2 * (size_t)d]) {
                bad("cross d=%u で不一致 (i=%zu)", d, i); break;
            }
    }

    printf("== 分岐間相関 cross: B 自身の履歴には必要ビットが無い ==\n");
    // B の予測に必要なのは A[i-1-2d]。B 自身の過去の結果は A の別の位置の
    // コピーなので、B の履歴だけでは決まらないことを確認する。
    for (unsigned d = 1; d <= 8; d++) {
        bp_gen_cross(pat, n, d, 0);
        size_t ntp = 0;
        for (size_t i = 1 + 2 * (size_t)d; i < n && ntp < 3000; i += 2) tp[ntp++] = i;
        // B の直近 h 個の「自分の結果」だけを並べた列を作って曖昧さを見る
        uint8_t *bhist = (uint8_t *)malloc(ntp);
        size_t *bpos = (size_t *)malloc(sizeof(size_t) * ntp);
        if (!bhist || !bpos) { free(bhist); free(bpos); break; }
        for (size_t k = 0; k < ntp; k++) bhist[k] = pat[tp[k]];
        for (size_t k = 0; k < ntp; k++) bpos[k] = k;
        // bhist 上で「直近 h 個」を窓として曖昧さを確認（h=8）
        if (!is_ambiguous_at(bhist, bpos, ntp, 8))
            bad("cross d=%u: B 自身の履歴 8 個で標的が決まってしまう", d);
        free(bhist); free(bpos);
    }

    printf("== aliasing 用 alias: 各サイトが周期 p を保つ ==\n");
    {
        unsigned Ss[] = {2, 8, 32};
        unsigned p = 8;
        for (unsigned si = 0; si < 3; si++) {
            unsigned S = Ss[si];
            bp_gen_alias(pat, n, S, p, 0);
            size_t iters = n / S;
            int okflag = 1;
            for (unsigned s = 0; s < S && okflag; s++)
                for (size_t i = 0; i + p < iters; i++)
                    if (pat[i * S + s] != pat[(i + p) * S + s]) {
                        bad("alias S=%u: サイト %u が周期 %u を保っていない", S, s, p);
                        okflag = 0; break;
                    }
        }
    }

    printf("== 二重相関 dual: 標的が XOR に一致 ==\n");
    {
        unsigned d1 = 8, d2 = 21, dm = 21, blk = 22;
        bp_gen_dual(pat, n, d1, d2, 0);
        for (size_t i = 0; i + blk <= n; i += blk)
            if (pat[i + dm] != (pat[i + dm - d1] ^ pat[i + dm - d2])) {
                bad("dual で不一致 (i=%zu)", i); break;
            }
    }

    printf("== 周期 period ==\n");
    bp_gen_period(pat, n, 10);
    for (size_t i = 0; i < 200; i++)
        if (pat[i] != (uint8_t)(((i % 10) != 9) ? 1 : 0)) { bad("period 不一致 (i=%zu)", i); break; }

    printf("== 偏り bias ==\n");
    for (double p = 0.5; p <= 0.95; p += 0.15) {
        bp_gen_bias(pat, n, p);
        size_t t = 0;
        for (size_t i = 0; i < n; i++) t += pat[i];
        double f = (double)t / (double)n;
        if (f < p - 0.03 || f > p + 0.03) bad("bias p=%.2f 実測 %.3f", p, f);
    }

    printf("== 対照系列: 標的の相関が消えている ==\n");
    {
        unsigned D = 15;
        bp_gen_hist_dense(pat, n, D, 1);
        size_t hit = 0, tot = 0;
        for (size_t i = (D | 1); i < n; i += 2) { tot++; if (pat[i] == pat[i - D]) hit++; }
        double f = (double)hit / (double)tot;
        if (f < 0.45 || f > 0.55) bad("histd 対照の一致率が 0.5 から離れている (%.3f)", f);
    }

    free(pat); free(tp);
    if (fails == 0) printf("\nselftest: 全項目 PASS\n");
    else            printf("\nselftest: %d 件 FAIL\n", fails);
    return fails ? 1 : 0;
}

int main(void) { return bp_selftest(); }
