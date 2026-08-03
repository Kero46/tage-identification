// bp.h — 測定基盤の共通インタフェース
#ifndef BP_H
#define BP_H

#include <stddef.h>
#include <stdint.h>

// パターン長の既定値。
// 【重要】パターンは反復して実行されるため、短いと予測器が系列全体を記憶する。
// 実測ではパターン長 1KB のとき純ランダム系列でも約 0.82 ns/分岐（=記憶済み）、
// 16KB で約 4.2 ns/分岐であった。記憶化領域に入ると対照系列も予測可能になり
// 差分法が壊れる。必要な長さは予測器の容量に依存するため機械ごとに校正が必要。
//   env/calibrate_patlen.sh で校正すること。
#ifndef BP_PAT_LEN_DEFAULT
#define BP_PAT_LEN_DEFAULT (1u << 16)
#endif

// 実行時に設定するパターン長（--patlen）。
extern size_t bp_pat_len;

// ------------------------------------------------------------------ 乱数
void     bp_seed(uint64_t s);
uint64_t bp_rnd_u64(void);
int      bp_rnd_bit(void);
int      bp_rnd_biased(double p);

// ------------------------------------------------------------------ パターン生成
// control != 0 のとき、相関を持たない対照系列を生成する。
void bp_gen_hist      (uint8_t *pat, size_t n, unsigned D, int control);
void bp_gen_hist_dense(uint8_t *pat, size_t n, unsigned D, int control);
void bp_gen_dual      (uint8_t *pat, size_t n, unsigned d1, unsigned d2, int control);
void bp_gen_period    (uint8_t *pat, size_t n, unsigned N);
void bp_gen_bias      (uint8_t *pat, size_t n, double p);
void bp_gen_random    (uint8_t *pat, size_t n);
void bp_gen_hist_ctx  (uint8_t *pat, size_t n, unsigned D, unsigned K, int control);

// 複数分岐サイト用。
// cross: 分岐 A（偶数位置）の結果が、d ステップ後の分岐 B（奇数位置）を決める。
//        B 自身の履歴には必要なビットが含まれないため、グローバル履歴を要する。
void bp_gen_cross(uint8_t *pat, size_t n, unsigned d, int control);
// alias: S 個のサイトそれぞれに周期 p の固有パターンを与える。
//        S を増やすとテーブル衝突（aliasing）が起こる可能性が高まる。
void bp_gen_alias(uint8_t *pat, size_t n, unsigned S, unsigned p, int control);

// ------------------------------------------------------------------ 計測
typedef struct {
    uint64_t branches, misses, cycles, instructions, fe_stalls;
    double   ns;
} bp_sample;

// want_frontend != 0 でフロントエンド停止サイクルも取る。
// PMU の物理カウンタ数を超えると多重化されるため、必要なときだけ有効にする。
void bp_ctr_init(int want_frontend);
int  bp_ctr_ok(void);
void bp_ctr_start(void);
void bp_ctr_stop(bp_sample *out);
double bp_now_ns(void);

// ------------------------------------------------------------------ カーネル
typedef uint64_t (*bp_kernel_fn)(const uint8_t *pat, size_t n, uint64_t reps, void *ctx);
uint64_t bp_kernel_single(const uint8_t *pat, size_t n, uint64_t reps, void *ctx);
// 2 サイト（偶数位置=A, 奇数位置=B）。グローバル履歴の評価に用いる。
uint64_t bp_kernel_cross(const uint8_t *pat, size_t n, uint64_t reps, void *ctx);
// S サイト（S は 2,4,8,16,32,64）。ctx に unsigned* で S を渡す。
uint64_t bp_kernel_sites(const uint8_t *pat, size_t n, uint64_t reps, void *ctx);

// ------------------------------------------------------------------ 試行ループ
typedef void (*bp_gen_fn)(uint8_t *pat, size_t n, double p1, double p2, int control);

typedef struct {
    const char  *mode;
    double       p1, p2;
    int          control;   // 単独測定時の対照フラグ
    int          pair;      // 1 なら test/control を同一プロセス内で交互に測る
    uint64_t     trials, reps, warmup, seed;
    int          csv;
    bp_gen_fn    gen;
    bp_kernel_fn kernel;
    void        *ctx;
} bp_config;

void bp_print_header(void);
int  bp_run(const bp_config *cfg);

// パターン生成器の性質検査（必要履歴長の下限を含む）。0 で PASS。
int bp_selftest(void);

// カウンタが多重化で縮尺されたか（1 なら値は信頼できない）
int bp_ctr_scaled(void);

#endif
