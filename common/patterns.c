// patterns.c — 分岐パターンの生成
//
// 共通の考え方: 標的分岐の予測可能性は「相関あり(test)」と
// 「相関なし(control)」の差分で取り出す。カウンタは分岐ミスの総数しか
// 返さないため、履歴を運ぶフィラー分岐の約50%のミスを相殺する必要がある。
#include "bp.h"
#include <stdlib.h>

// 低密度版: ブロック長 D+1。index 0..D-1 はフィラー(独立乱数)、
// index D が標的で値はブロック先頭(=ちょうど D 個前)に一致。
// 標的が全分岐の 1/(D+1) しかないため、D が大きいと信号が埋もれる。
void bp_gen_hist(uint8_t *pat, size_t n, unsigned D, int control) {
    unsigned blk = D + 1;
    for (size_t i = 0; i < n; i += blk) {
        size_t lim = (i + blk <= n) ? blk : (n - i);
        for (size_t j = 0; j < lim; j++) pat[i + j] = (uint8_t)bp_rnd_bit();
        if (lim == blk) pat[i + D] = control ? (uint8_t)bp_rnd_bit() : pat[i];
    }
}

// 高密度版: 偶数位置=フィラー、奇数位置=標的。D を奇数にすると
// i-D は必ず偶数=フィラーになる。標的が 50% を占めるため信号が D に依存しない。
//
// 短い履歴で当てられないことの保証: 標的 i は filler[i-D] に等しい。
// 窓 i-1..i-h に filler[i-D] が入るのは h >= D のときのみ。
void bp_gen_hist_dense(uint8_t *pat, size_t n, unsigned D, int control) {
    for (size_t i = 0; i < n; i++) pat[i] = (uint8_t)bp_rnd_bit();
    for (size_t i = 1; i < n; i += 2) {
        if (i >= D) pat[i] = control ? (uint8_t)bp_rnd_bit() : pat[i - D];
        else        pat[i] = (uint8_t)bp_rnd_bit();
    }
}

// 文脈数を制限した相関系列（テーブルサイズ測定用）。
//
// 【設計上の要点】素朴に「長さ D の文脈語を K 本用意する」構成は誤りである。
// 文脈語がランダムだと直近の数ビットで語が識別できてしまい、D ビットの履歴を
// 必要としない。実測では D=63, K=1024 のとき直近 27 ビットで語が一意に決まり、
// D=63 のテーブルではなく 27 ビット相当のテーブルを見ていた。
//
// そこで識別情報を「遠く」に置き、近傍は全文脈で共通にする。
//   ブロック = [識別ビット j 個][共通サフィックス D-j 個][標的]
//   標的 = 文脈ごとに固定した乱数ビット tgt[k]
// 識別ビットは標的から見て距離 D..D-j+1 にあるため、標的を当てるには
// 履歴長 >= D が必要。かつ K=2^j 個の文脈を区別するので K エントリを要する。
// 共通サフィックスは全ブロックで同一なので情報を持たず、学習後は予測される。
//
// K は 2 のべきに切り上げて用いる（j = ceil(log2 K)）。D-j >= 1 が必要。
// 標的値は半分ずつ 0/1 に均衡させる（縮退の防止と成立率の固定）。
void bp_gen_hist_ctx(uint8_t *pat, size_t n, unsigned D, unsigned K, int control) {
    if (K < 2) K = 2;
    unsigned j = 0;
    while ((1u << j) < K) j++;
    K = 1u << j;
    if (D <= j) {           // 識別ビットが入らない
        bp_gen_random(pat, n);
        return;
    }
    unsigned suflen = D - j;

    uint8_t *suffix = (uint8_t *)malloc(suflen);
    uint8_t *tgt    = (uint8_t *)malloc(K);
    if (!suffix || !tgt) { free(suffix); free(tgt); return; }
    // 共通サフィックス（全ブロックで同一。固定パターンなので学習可能）
    for (unsigned i = 0; i < suflen; i++) suffix[i] = (uint8_t)bp_rnd_bit();
    // 文脈ごとの標的値。独立乱数にすると K が小さいとき全文脈が同値になる縮退が
    // 起こり（K=4 なら 12.5%）、容量に関係なく予測可能になってしまう。
    // 半分を 1、半分を 0 にして混ぜることで縮退を防ぎ、同時に成立率を K に
    // 依らず 0.5 に固定する（K を変えたとき偏りが変わる交絡も避けられる）。
    for (unsigned k = 0; k < K; k++) tgt[k] = (uint8_t)(k < K / 2 ? 1 : 0);
    for (unsigned k = K; k > 1; k--) {          // Fisher-Yates
        unsigned r = (unsigned)(bp_rnd_u64() % k);
        uint8_t t = tgt[k - 1]; tgt[k - 1] = tgt[r]; tgt[r] = t;
    }

    size_t i = 0;
    while (i + D + 1 <= n) {
        unsigned k = (unsigned)(bp_rnd_u64() % K);
        for (unsigned b = 0; b < j; b++)               // 識別ビット（遠い側）
            pat[i + b] = (uint8_t)((k >> b) & 1u);
        for (unsigned b = 0; b < suflen; b++)          // 共通サフィックス（近い側）
            pat[i + j + b] = suffix[b];
        pat[i + D] = control ? (uint8_t)bp_rnd_bit() : tgt[k];
        i += D + 1;
    }
    while (i < n) pat[i++] = (uint8_t)bp_rnd_bit();
    free(suffix); free(tgt);
}

// 二重相関: 標的 = bit(d1前) XOR bit(d2前)。段の分離に用いる。
void bp_gen_dual(uint8_t *pat, size_t n, unsigned d1, unsigned d2, int control) {
    unsigned dm = (d1 > d2 ? d1 : d2);
    unsigned blk = dm + 1;
    for (size_t i = 0; i < n; i += blk) {
        size_t lim = (i + blk <= n) ? blk : (n - i);
        for (size_t j = 0; j < lim; j++) pat[i + j] = (uint8_t)bp_rnd_bit();
        if (lim == blk)
            pat[i + dm] = control ? (uint8_t)bp_rnd_bit()
                                  : (uint8_t)(pat[i + dm - d1] ^ pat[i + dm - d2]);
    }
}

void bp_gen_period(uint8_t *pat, size_t n, unsigned N) {
    if (N == 0) N = 1;
    for (size_t i = 0; i < n; i++) pat[i] = (uint8_t)(((i % N) != (N - 1)) ? 1 : 0);
}
void bp_gen_bias(uint8_t *pat, size_t n, double p) {
    for (size_t i = 0; i < n; i++) pat[i] = (uint8_t)bp_rnd_biased(p);
}
void bp_gen_random(uint8_t *pat, size_t n) {
    for (size_t i = 0; i < n; i++) pat[i] = (uint8_t)bp_rnd_bit();
}

// ---------------------------------------------------------------- 複数分岐サイト
// 分岐間相関（グローバル履歴の評価）。
//
// 偶数位置が分岐 A、奇数位置が分岐 B の結果に対応する（カーネル側で
// 別の分岐命令として実行される）。B の結果は d ステップ前の A の結果に一致させる。
//
// 【要点】単一の分岐サイトで自分自身の過去と相関させる構成では、ローカル履歴
// だけで学習できてしまうため、グローバル履歴の効果を評価できない。ここでは
// B を予測するのに必要なビットが A の結果であり、B 自身の過去の結果
// （B[i-1]=A[i-3-2d], B[i-2]=A[i-5-2d], ...）には含まれないため、
// グローバル履歴が必要になる。
//
// 【必要なグローバル履歴長】B[i] = A[i-1-2d] であり、A と B が交互に実行される
// ため、必要なビットは分岐列で 2d+1 個前に位置する。したがってパラメータ d に
// 対して必要なグローバル履歴長は 2d+1 分岐であり、d そのものではない。
// 掃引結果を解釈する際はこの換算を用いる。
void bp_gen_cross(uint8_t *pat, size_t n, unsigned d, int control) {
    if (d == 0) d = 1;
    for (size_t i = 0; i + 1 < n; i += 2)
        pat[i] = (uint8_t)bp_rnd_bit();               // A: 独立乱数
    for (size_t i = 1; i < n; i += 2) {
        size_t src = (i - 1) - (size_t)2 * d;         // d ステップ前の A の位置
        if (!control && i >= 1 + (size_t)2 * d) pat[i] = pat[src];
        else                                    pat[i] = (uint8_t)bp_rnd_bit();
    }
}

// aliasing 実験用。S 個のサイトに、それぞれ周期 p の固有パターンを与える。
// 配置: 反復 i のサイト s は pat[i*S + s]。
//
// 各サイトのパターンは単体なら学習可能（周期 p）。S を増やすと予測テーブルに
// 載せるべき対象が増え、エントリ衝突（aliasing）が起こる可能性が高まる。
// ただし S を増やすとコード量も増えるため、フロントエンド側の影響と完全には
// 分離できない。この実験は「aliasing が生じる可能性を高める」ものであり、
// 観測された劣化を aliasing と断定するものではない。
void bp_gen_alias(uint8_t *pat, size_t n, unsigned S, unsigned p, int control) {
    if (S == 0) S = 1;
    if (p == 0) p = 2;
    uint8_t *w = (uint8_t *)malloc((size_t)S * p);
    if (!w) return;
    for (size_t k = 0; k < (size_t)S * p; k++) w[k] = (uint8_t)bp_rnd_bit();

    size_t iters = n / S;
    for (size_t i = 0; i < iters; i++)
        for (unsigned s = 0; s < S; s++)
            pat[i * S + s] = control ? (uint8_t)bp_rnd_bit()
                                     : w[(size_t)s * p + (i % p)];
    for (size_t i = iters * S; i < n; i++) pat[i] = (uint8_t)bp_rnd_bit();
    free(w);
}
