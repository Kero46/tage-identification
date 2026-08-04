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
#include <math.h>

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
// 【検査手法の注意（落とし穴 7）】「直近 h ビットが同じで標的が異なる例が無い」
// ことを、そのまま「h ビットで標的が決まる」の根拠にしてはならない。有限
// サンプルでは h が大きいと窓が全て相異なるため、決定性ではなく一意性を見て
// しまう（D=17 で「15 ビットで決まる」と誤判定した）。
//
// 対策として、判定の前に**検出力そのものを測る**。同じ窓を持つ標的ペアが
// 何組観測できたかを数え、それが足りなければ「検査不能」として失敗させる。
// 反例の不在を証拠に使えるのは、反例が現れうる状況を作れているときだけ。
// これを関数の中に固定したので、呼び出し側が h とサンプル数の関係を
// 気にしなくても、同じ誤りが黙って再発しない。
//
// pairs     : 同じ窓を持つ標的ペアの数（= 検査の検出力）
// conflicts : そのうち標的の値が異なるペアの数（= 曖昧さの証拠）
static void window_pairs(const uint8_t *pat, const size_t *tpos, size_t ntp,
                         unsigned h, size_t *pairs, size_t *conflicts) {
    size_t p = 0, c = 0;
    for (size_t a = 0; a < ntp; a++) {
        if (tpos[a] < h) continue;
        for (size_t b = a + 1; b < ntp; b++) {
            if (tpos[b] < h) continue;
            if (memcmp(&pat[tpos[a] - h], &pat[tpos[b] - h], h) != 0) continue;
            p++;
            if (pat[tpos[a]] != pat[tpos[b]]) c++;
        }
    }
    *pairs = p; *conflicts = c;
}

// 検査の成立に必要な同一窓ペアの数。これ未満なら検査不能として失敗させる。
// 窓が重複していない状態では、決定性も曖昧さも判定できない。
#define BP_MIN_PAIRS 200

// 窓 h ビットでは標的が決まらないこと（曖昧さが残ること）を確認する。
static int require_ambiguous(const char *what, const uint8_t *pat,
                             const size_t *tpos, size_t ntp, unsigned h) {
    size_t pairs, conflicts;
    window_pairs(pat, tpos, ntp, h, &pairs, &conflicts);
    if (pairs < BP_MIN_PAIRS) {
        bad("%s: 窓 %u ビットの同一窓ペアが %zu 組しかなく検査不能"
            "（%d 組以上必要。サンプル数か窓幅を見直すこと）", what, h, pairs,
            BP_MIN_PAIRS);
        return 0;
    }
    if (conflicts == 0) {
        bad("%s: 窓 %u ビットで標的が決まってしまう（同一窓 %zu 組すべてで"
            "標的が一致）", what, h, pairs);
        return 0;
    }
    return 1;
}

// 窓 h ビットで標的が一貫して決まることを確認する。
// 反例の不在を根拠にするので、同一窓ペアが十分に観測できていることが前提。
// その前提を関数内で検査するため、呼び出し側で 2^h とサンプル数の関係を
// 見積もる必要はない（見積もりを人手に任せたのが落とし穴 7 の原因だった）。
static int require_determined(const char *what, const uint8_t *pat,
                              const size_t *tpos, size_t ntp, unsigned h) {
    size_t pairs, conflicts;
    window_pairs(pat, tpos, ntp, h, &pairs, &conflicts);
    if (pairs < BP_MIN_PAIRS) {
        bad("%s: 窓 %u ビットの同一窓ペアが %zu 組しかなく検査不能"
            "（%d 組以上必要。窓が重複していないので決定性は判定できない）",
            what, h, pairs, BP_MIN_PAIRS);
        return 0;
    }
    if (conflicts > 0) {
        bad("%s: 窓 %u ビットでも標的が決まらない（同一窓 %zu 組のうち %zu 組で"
            "標的が不一致）", what, h, pairs, conflicts);
        return 0;
    }
    return 1;
}

// ---------------------------------------------------------------- 成立率
// 差分法は「構造が同一で標的の相関だけが異なる二系列」を前提にする。
// 実測で、完全に予測できているパターンでも成立率が変わると実行時間が
// 1.6 倍動くことが分かった（lab-notebook 2026-08-03/04 の診断 6）。
// test と control で成立率がずれていると、差分に予測性能以外の要因が混ざる。
// よって成立率の一致は差分法の前提として検査する。
static double taken_rate(const uint8_t *pat, size_t n) {
    size_t t = 0;
    for (size_t i = 0; i < n; i++) t += pat[i] ? 1u : 0u;
    return (double)t / (double)n;
}

// test と control の成立率が許容差以内で一致することを確認する。
static void check_rate_match(const char *what, const uint8_t *t, const uint8_t *c,
                             size_t n, double tol) {
    double rt = taken_rate(t, n), rc = taken_rate(c, n);
    double d = rt - rc;
    if (d < 0) d = -d;
    printf("   %-18s test %.4f  control %.4f  差 %+.4f  (許容 %.3f)%s\n",
           what, rt, rc, rt - rc, tol, (d > tol) ? "  ← NG" : "");
    if (d > tol)
        bad("%s: test と control の成立率が %.4f ずれている（許容 %.3f）",
            what, d, tol);
}

// 相関が消えていること（一致率が 0.5 付近であること）を確認する。
static void check_half(const char *what, size_t hit, size_t tot) {
    if (tot == 0) { bad("%s: 標本が無く検査不能", what); return; }
    double f = (double)hit / (double)tot;
    if (f < 0.45 || f > 0.55)
        bad("%s 対照: 一致率が 0.5 から離れている (%.3f, 標本 %zu)", what, f, tot);
}

int bp_selftest(void) {
    // サンプル数。ctx は D+1 ごとに 1 標本しか取れないため、D=63 だと
    // 2^13 では標本 128 個・同一窓ペア 106 組しか作れず、決定性の検査が
    // 成立しなかった（require_determined が検査不能で落ちた）。
    // 閾値を下げるのではなく標本を増やして対処する。
    size_t n = 1u << 16;
    uint8_t *pat  = (uint8_t *)malloc(n);
    uint8_t *pat2 = (uint8_t *)malloc(n);          // 対照系列との比較用
    size_t  *tp   = (size_t *)malloc(sizeof(size_t) * n);
    if (!pat || !pat2 || !tp) { fprintf(stderr, "alloc 失敗\n"); return 1; }
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
    // h を小さく固定して曖昧さを確認する。検出力（同一窓ペア数）が足りるかは
    // require_ambiguous が内部で検査する。
    {
        const unsigned h = 10;              // 2^10 = 1024 通り
        for (unsigned D = 13; D <= 41; D += 4) {
            bp_gen_hist_dense(pat, n, D, 0);
            size_t ntp = 0;
            for (size_t i = (D | 1); i < n && ntp < 4000; i += 2) tp[ntp++] = i;
            char what[64]; snprintf(what, sizeof what, "histd D=%u", D);
            require_ambiguous(what, pat, tp, ntp, h);
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
                for (size_t i = D; i + 1 <= n && ntp < 1200; i += D + 1) tp[ntp++] = i;

                char what[64];
                snprintf(what, sizeof what, "ctx D=%u K=%u（識別情報の漏れ）", D, K);
                require_ambiguous(what, pat, tp, ntp, suflen);
                // h = D は窓が長いが、設計上この窓は文脈 k で決まり K 通りしか
                // 無いため、2000 サンプルで同一窓が多数重複する。よって
                // 「反例が無い」ことを決定性の根拠にできる。その前提が本当に
                // 成り立っているか（同一窓ペア数）は require_determined が
                // 実行時に検査するので、成り立たなくなれば検査不能で落ちる。
                snprintf(what, sizeof what, "ctx D=%u K=%u（履歴 D での決定性）", D, K);
                require_determined(what, pat, tp, ntp, D);
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
        char what[64]; snprintf(what, sizeof what, "cross d=%u（B 自身の履歴）", d);
        require_ambiguous(what, bhist, bpos, ntp, 8);
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

    printf("== カーネル同値: 多サイト版が単一版と同じ判定列を実行する ==\n");
    // --sites で多サイトカーネルに切り替えたとき、実行される判定の列が変わって
    // しまうと、履歴の単位を確かめる実験（要素あたりのループ分岐の本数だけを
    // 変える）が成立しない。同じパターンに対する sink が一致することで、
    // サイト分割がパターンの読み位置をずらしていないことを確認する。
    // n は 64 の倍数なので、どの S でも全要素が処理される。
    {
        unsigned Ss[] = {2, 4, 8, 16, 32, 64};
        bp_seed(4242);
        bp_gen_random(pat, n);
        uint64_t ref = bp_kernel_single(pat, n, 1, NULL);
        for (unsigned i = 0; i < 6; i++) {
            uint64_t got = bp_kernel_sites(pat, n, 1, &Ss[i]);
            if (got != ref)
                bad("kernel_sites_%u: 単一版と判定列が違う (%llu 対 %llu)",
                    Ss[i], (unsigned long long)got, (unsigned long long)ref);
        }
        uint64_t cross = bp_kernel_cross(pat, n, 1, NULL);
        if (cross != ref)
            bad("bp_kernel_cross: 単一版と判定列が違う (%llu 対 %llu)",
                (unsigned long long)cross, (unsigned long long)ref);
        // 詰め物つきカーネルも、詰め物の方向に依らず同じ判定列を実行するはず。
        // 詰め物は sink に触らないので sink が変わってはいけない。
        {
            unsigned Ps[] = {0, 1, 2, 4};
            for (unsigned pi = 0; pi < 4; pi++)
                for (unsigned dir = 0; dir <= 1; dir++) {
                    bp_pad_cfg pc; pc.pads = Ps[pi]; pc.dir = dir;
                    uint64_t got = bp_kernel_pad(pat, n, 1, &pc);
                    if (got != ref)
                        bad("kernel_pad_%u(dir=%u): 単一版と判定列が違う"
                            " (%llu 対 %llu)", Ps[pi], dir,
                            (unsigned long long)got, (unsigned long long)ref);
                }
        }

        // ctx 系列でも確認する（この実験で実際に使う系列）
        bp_gen_hist_ctx(pat, n, 66, 2, 0);
        ref = bp_kernel_single(pat, n, 1, NULL);
        for (unsigned i = 0; i < 6; i++) {
            uint64_t got = bp_kernel_sites(pat, n, 1, &Ss[i]);
            if (got != ref)
                bad("kernel_sites_%u（ctx 系列）: 単一版と判定列が違う", Ss[i]);
        }
    }

    printf("== 対照系列 (a): 標的の相関が消えている ==\n");
    {
        size_t hit, tot;

        // hist: 標的がブロック先頭と一致しなくなる
        unsigned D = 16;
        bp_gen_hist(pat, n, D, 1);
        hit = tot = 0;
        for (size_t i = 0; i + D + 1 <= n; i += D + 1) {
            tot++; if (pat[i + D] == pat[i]) hit++;
        }
        check_half("hist D=16", hit, tot);

        // hist_dense: 標的が距離 D の値と一致しなくなる
        for (unsigned Dd = 15; Dd <= 101; Dd += 86) {
            bp_gen_hist_dense(pat, n, Dd, 1);
            hit = tot = 0;
            for (size_t i = (Dd | 1); i < n; i += 2) {
                tot++; if (pat[i] == pat[i - Dd]) hit++;
            }
            char what[64]; snprintf(what, sizeof what, "histd D=%u", Dd);
            check_half(what, hit, tot);
        }

        // dual: 標的が XOR と一致しなくなる
        unsigned d1 = 8, d2 = 21, dm = 21, blk = 22;
        bp_gen_dual(pat, n, d1, d2, 1);
        hit = tot = 0;
        for (size_t i = 0; i + blk <= n; i += blk) {
            tot++;
            if (pat[i + dm] == (uint8_t)(pat[i + dm - d1] ^ pat[i + dm - d2])) hit++;
        }
        check_half("dual 8,21", hit, tot);

        // cross: B が d ステップ前の A と一致しなくなる
        for (unsigned d = 5; d <= 30; d += 25) {
            bp_gen_cross(pat, n, d, 1);
            hit = tot = 0;
            for (size_t i = 1 + 2 * (size_t)d; i < n; i += 2) {
                tot++; if (pat[i] == pat[(i - 1) - 2 * (size_t)d]) hit++;
            }
            char what[64]; snprintf(what, sizeof what, "cross d=%u", d);
            check_half(what, hit, tot);
        }

        // alias: 各サイトの周期性が崩れる
        {
            unsigned S = 8, p = 8;
            bp_gen_alias(pat, n, S, p, 1);
            size_t iters = n / S;
            hit = tot = 0;
            for (unsigned s = 0; s < S; s++)
                for (size_t i = 0; i + p < iters; i++) {
                    tot++;
                    if (pat[i * S + s] == pat[(i + p) * S + s]) hit++;
                }
            check_half("alias S=8 p=8", hit, tot);
        }

        // ctx: 履歴 D ビットでも標的が決まらなくなる（文脈と標的の対応が消える）
        {
            unsigned Dc = 31, K = 16;
            bp_gen_hist_ctx(pat, n, Dc, K, 1);
            size_t ntp = 0;
            for (size_t i = Dc; i + 1 <= n && ntp < 1200; i += Dc + 1) tp[ntp++] = i;
            require_ambiguous("ctx D=31 K=16 対照", pat, tp, ntp, Dc);
        }
    }

    printf("== 対照系列 (b): test と control で成立率が一致 ==\n");
    // 差分法は「構造が同一で標的の相関だけが異なる二系列」を前提にする。
    // 実測で、完全に予測できているパターンでも成立率が変わると実行時間が
    // 1.6 倍動いた（lab-notebook 2026-08-03/04 の診断 6）。成立率がずれていれば
    // 差分に予測性能以外の要因が混ざるため、前提として検査する。
    // 種は runner.c と同じく test と control で同一にする。
    {
        const double tol = 0.01;            // 実測のずれは全て 0.005 以下
        const uint64_t sd = 20260804;

        bp_seed(sd); bp_gen_hist(pat,  n, 16, 0);
        bp_seed(sd); bp_gen_hist(pat2, n, 16, 1);
        check_rate_match("hist D=16", pat, pat2, n, tol);

        bp_seed(sd); bp_gen_hist_dense(pat,  n, 15, 0);
        bp_seed(sd); bp_gen_hist_dense(pat2, n, 15, 1);
        check_rate_match("histd D=15", pat, pat2, n, tol);

        bp_seed(sd); bp_gen_hist_dense(pat,  n, 101, 0);
        bp_seed(sd); bp_gen_hist_dense(pat2, n, 101, 1);
        check_rate_match("histd D=101", pat, pat2, n, tol);

        // ctx の成立率は 0.5 ではない（共通サフィックスの固定値で決まる）。
        // ここで要るのは 0.5 であることではなく test と control の一致である。
        unsigned Ds[] = {15, 31, 63}, Ks[] = {4, 16, 1024};
        for (unsigned i = 0; i < 3; i++) {
            bp_seed(sd); bp_gen_hist_ctx(pat,  n, Ds[i], Ks[i], 0);
            bp_seed(sd); bp_gen_hist_ctx(pat2, n, Ds[i], Ks[i], 1);
            char what[64];
            snprintf(what, sizeof what, "ctx D=%u K=%u", Ds[i], Ks[i]);
            check_rate_match(what, pat, pat2, n, tol);
        }

        bp_seed(sd); bp_gen_dual(pat,  n, 8, 21, 0);
        bp_seed(sd); bp_gen_dual(pat2, n, 8, 21, 1);
        check_rate_match("dual 8,21", pat, pat2, n, tol);

        for (unsigned d = 5; d <= 30; d += 25) {
            bp_seed(sd); bp_gen_cross(pat,  n, d, 0);
            bp_seed(sd); bp_gen_cross(pat2, n, d, 1);
            char what[64]; snprintf(what, sizeof what, "cross d=%u", d);
            check_rate_match(what, pat, pat2, n, tol);
        }
    }

    printf("== 対照系列 (b) の例外: alias は成立率が構造的にずれる ==\n");
    // alias の test 系列の成立率は、周期語 w（S*p ビットの独立乱数）の平均で
    // 決まる。標本が S*p 個しかないため 0.5 から系統的に外れ、control（全長に
    // わたる独立乱数）と一致しない。実測では S=8,p=8 で 0.113、種によっては
    // 0.376 ずれた。これは生成器の欠陥ではなく構成上の帰結である。
    //
    // したがって alias を差分法で解析してはならない。仕様 §4.2 が
    // 「--no-control で test_ns の S 依存を見る」と定めているのはこのためでもある
    // （元の理由は「対照がランダムなので緩やかな劣化が埋もれる」）。
    // ここでは、ずれが有限標本で説明できる範囲に収まることだけを確認する。
    // これを超えるずれが出たら、別の原因（生成器の欠陥）を疑うべきである。
    {
        const uint64_t sd = 20260804;
        unsigned Ss[] = {2, 8, 64, 8}, ps[] = {8, 8, 8, 2};
        for (unsigned i = 0; i < 4; i++) {
            unsigned S = Ss[i], p = ps[i];
            bp_seed(sd); bp_gen_alias(pat,  n, S, p, 0);
            bp_seed(sd); bp_gen_alias(pat2, n, S, p, 1);
            // 独立乱数 S*p 個の平均の標準偏差は 0.5/sqrt(S*p)。その 4 倍まで許す。
            double bound = 4.0 * 0.5 / sqrt((double)S * (double)p) + 0.01;
            char what[64]; snprintf(what, sizeof what, "alias S=%u p=%u", S, p);
            check_rate_match(what, pat, pat2, n, bound);
        }
    }

    free(pat); free(pat2); free(tp);
    if (fails == 0) printf("\nselftest: 全項目 PASS\n");
    else            printf("\nselftest: %d 件 FAIL\n", fails);
    return fails ? 1 : 0;
}

int main(void) { return bp_selftest(); }
