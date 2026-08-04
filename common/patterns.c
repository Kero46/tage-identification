// patterns.c — 分岐パターンの生成
//
// 共通の考え方: 標的分岐の予測可能性は「相関あり(test)」と
// 「相関なし(control)」の差分で取り出す。カウンタは分岐ミスの総数しか
// 返さないため、履歴を運ぶフィラー分岐の約50%のミスを相殺する必要がある。
#include "bp.h"
#include <stdio.h>
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
/* ctx で表現できる最大の文脈数。j <= D-1 が必要なので K <= 2^(D-1)。
 * D=1 では識別ビットを置けないので 0 を返す（表現不能）。 */
unsigned bp_ctx_max_k(unsigned D) {
    if (D < 2) return 0;
    if (D - 1 >= 31) return 1u << 31;      /* オーバーフロー回避 */
    return 1u << (D - 1);
}

void bp_gen_hist_ctx(uint8_t *pat, size_t n, unsigned D, unsigned K, int control) {
    if (K < 2) K = 2;
    unsigned j = 0;
    while ((1u << j) < K) j++;
    K = 1u << j;
    if (D <= j) {
        /* 識別ビットが入らない。**黙って純ランダムを返してはならない。**
         *
         * 以前はここで bp_gen_random() して return していた。呼び出し側から見ると
         * 「ctx を測っているつもりで純ランダムを測っている」状態になり、K 掃引に
         * 偽の劣化が現れる（実測: D=5, K=32 で test_ns が純ランダムの飽和値
         * 3.38 になり、K=16→32 が容量の限界に見えた）。黙って別のものを測るのは
         * 最悪の失敗様式なので、落とし穴 1（分岐の消去）・2（サイトの複製）と
         * 同じ方針で中断させる。
         *
         * 通常は CLI が先に拒否する（bp_ctx_max_k）。ここは将来の呼び出し元が
         * 検証を忘れた場合の最後の防具である。 */
        fprintf(stderr,
                "bp_gen_hist_ctx: D=%u に対し K=%u は表現できません"
                "（識別ビット j=%u >= D）。K <= 2^(D-1) = %u にしてください。\n"
                "  以前はここで黙って純ランダム系列を返していました。"
                "偽の劣化を容量の限界と誤読するため中断します。\n",
                D, K, j, bp_ctx_max_k(D));
        abort();
    }
    unsigned suflen = D - j;

    uint8_t *suffix = (uint8_t *)malloc(suflen);
    uint8_t *tgt    = (uint8_t *)malloc(K);
    if (!suffix || !tgt) { free(suffix); free(tgt); return; }
    // 共通サフィックス（全ブロックで同一。固定パターンなので学習可能）。
    // 独立乱数にすると、系列全体の成立率がこの引き方で決まってしまう。
    // サフィックス長は D-j = D-log2(K) なので、K や D を変えるたびに成立率が
    // 動き、K を掃引したときの絶対時間の基準線がずれる（差分法自体は test と
    // control で成立率が一致するので成立するが、絶対値の比較ができない）。
    // 半分を 1、半分を 0 にして混ぜることで成立率をほぼ 0.5 に安定させる。
    // 全ブロックで同一である限り情報を持たないので、「近傍のビットでは文脈を
    // 識別できない」という設計要件は保たれる（自己検査の曖昧さ判定が担保する）。
    for (unsigned i = 0; i < suflen; i++) suffix[i] = (uint8_t)(i < suflen / 2 ? 1 : 0);
    for (unsigned i = suflen; i > 1; i--) {      // Fisher-Yates
        unsigned r = (unsigned)(bp_rnd_u64() % i);
        uint8_t t = suffix[i - 1]; suffix[i - 1] = suffix[r]; suffix[r] = t;
    }
    // 文脈ごとの標的値。独立乱数にすると K が小さいとき全文脈が同値になる縮退が
    // 起こり（K=4 なら 12.5%）、容量に関係なく予測可能になってしまう。
    // 半分を 1、半分を 0 にして混ぜることで縮退を防ぐ。
    //
    // 【注意】これが固定するのは標的ビットの寄与だけである。系列全体の成立率は
    // ブロックの大半を占める長さ D-j の共通サフィックスが支配するので、標的値を
    // 均衡させても成立率は 0.5 にならない（実測 0.40〜0.72）。しかも
    // j = log2(K) なので K を変えるとサフィックス長が変わり、成立率も動く。
    // そこでサフィックス自体も均衡させる（下記）。
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

// 雑音注入 probe（第2段 (a) の主手法）。梯子の段の位置 L_i を求める。
//
// ブロック構造（左が古い側）:
//   [雑音 m 個][識別ビット 1 個][共通サフィックス D-1 個][標的]
//
// 標的から見た距離:
//   距離 1..D-1   共通サフィックス（全ブロック同一。情報を持たない）
//   距離 D        識別ビット（標的を決める。文脈は 2 通り）
//   距離 D+1..D+m 雑音（独立乱数）
//
// 【原理】距離 D より**古い側**に雑音を置くと、履歴長ごとにテーブルの運命が分かれる。
//   L_i < D   識別ビットが窓の外          → 捉えられない
//   L_i = D   サフィックス + 識別ビットのみ → 使える（文脈 2 通り）
//   L_i > D   雑音を min(L_i-D, m) スロット見る。文脈が 2^(1+その数) に膨れ、
//             C_i を超えると散逸する      → 使えない
//
// つまり雑音は「D より長いテーブルだけを選択的に潰す」道具である。学習の成否は
// 「使えるテーブルが 1 つ以上あるか」で決まる。
//
//   成功(D) ⟺ ∃i : D <= L_i <= D + W_i     W_i ≈ log2(C_i / 2)
//
// 各テーブルは D ∈ [L_i - W_i, L_i] で成功に寄与するので、成功領域は区間の和集合。
// **孤立した成功区間の上端がちょうど L_i** である（D=L_i では窓の下端に入って
// 使えるが、D=L_i+1 では L_i < D となり捉えられない）。
//
// 【m=0 は既存の結果の特殊ケース】雑音がないのでどのテーブルも潰れず、成功条件が
// 「∃ L_i >= D」に退化する。最後に成功する D は最長履歴長そのもので、第1段で測った
// ctx K=2 の 66/67 境界と一致する。**ブロック構造も bp_gen_hist_ctx(D, K=2) と
// 同一**なので、実装の最初の検査は「m=0 で 66/67 が再現すること」である。
//
// 【m は W を超えないと何も潰れない】L_i > D のテーブルが見る雑音は
// min(L_i-D, m) スロットなので、m <= W ではどんなに L_i が長くても潰れず m=0 に
// 退化する。W ≈ 10 なら m ≈ 20 を使う（仕様 §4.1、落とし穴 28）。
//
// 【分解能】隣接する段が分離するのは梯子の間隔が W より広い場合だけ。間隔が W より
// 狭い段は 1 つの成功区間に融合する（落とし穴 29）。
//
// 文脈数は構成上 2 に固定である（識別ビット 1 個）。K は引数に取らない。
void bp_gen_ctx_noise(uint8_t *pat, size_t n, unsigned D, unsigned m, int control) {
    if (D < 1) D = 1;
    unsigned suflen = D - 1;              // 距離 1..D-1
    size_t blk = (size_t)m + D + 1;

    uint8_t *suffix = suflen ? (uint8_t *)malloc(suflen) : NULL;
    if (suflen && !suffix) return;
    // 共通サフィックス。全ブロックで同一なので情報を持たない。成立率を 0.5 付近に
    // 固定するため半分ずつ 0/1 に均衡させる（ctx と同じ理由。§1.9）。
    for (unsigned i = 0; i < suflen; i++) suffix[i] = (uint8_t)(i < suflen / 2 ? 1 : 0);
    for (unsigned i = suflen; i > 1; i--) {          // Fisher-Yates
        unsigned r = (unsigned)(bp_rnd_u64() % i);
        uint8_t t = suffix[i - 1]; suffix[i - 1] = suffix[r]; suffix[r] = t;
    }
    // 文脈ごとの標的値。2 通りなので必ず 0 と 1 を 1 つずつ持たせる（縮退の防止）。
    uint8_t tgt[2];
    tgt[0] = (uint8_t)(bp_rnd_bit());
    tgt[1] = (uint8_t)(1 - tgt[0]);

    size_t i = 0;
    while (i + blk <= n) {
        for (unsigned b = 0; b < m; b++)                    // 雑音（最も古い側）
            pat[i + b] = (uint8_t)bp_rnd_bit();
        unsigned k = (unsigned)(bp_rnd_u64() & 1u);         // 識別ビット
        pat[i + m] = (uint8_t)k;
        for (unsigned b = 0; b < suflen; b++)               // 共通サフィックス
            pat[i + m + 1 + b] = suffix[b];
        pat[i + m + D] = control ? (uint8_t)bp_rnd_bit() : tgt[k];
        i += blk;
    }
    while (i < n) pat[i++] = (uint8_t)bp_rnd_bit();
    free(suffix);
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
// 【相関距離】B[i] = A[i-1-2d] であり、A と B が交互に実行されるため、必要な
// ビットはパターン列で 2d+1 要素前に位置する。したがってパラメータ d に対して
// 相関距離は 2d+1 パターン要素であり、d そのものではない。
//
// **これを履歴長と読んではならない。** 二重に誤りになる。
//   (1) 単位が違う。ハードウェアの履歴が数えるのはパターン要素ではなく成立分岐
//       であり、不成立分岐は履歴を消費しない（仕様 §3.2）
//   (2) cross の転移が測っているものが違う。標的が「距離 L の 1 ビット」で
//       決まる構成なので 2^L 通りの履歴文脈にエントリが必要で、転移は容量側の
//       制約を含む（仕様 §3.4、第2段 (a)）
// 履歴長の下限は ctx の K=2 で文脈数を絞った測定から述べる。
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
