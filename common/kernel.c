// kernel.c — 計測対象のループ(単一の条件分岐)
//
// パターンをデータ(配列)で表現し、単一の条件分岐が要素を読んで分岐する。
// 掃引パラメータを変えてもコードサイズ・分岐密度が変わらないため、
// フロントエンド側の交絡が掃引に伴って変動しない。
//
// 分岐先の両方に空 asm を置き、cmov 化と分岐の消去を防ぐ。
// 環境ごとに必ず objdump で条件分岐が残っていることを確認すること。
#include "bp.h"

// 計測カーネルに共通で付ける属性。
//   noinline                        : 呼び出し側に溶けないようにする
//   no-unroll-loops                 : 掃引間でコード形状を一定に保つ
//   reorder-blocks-algorithm=simple : ブロック再配置による分岐サイトの複製を防ぐ。
//     既定の再配置では後続の比較が分岐先ごとに複製され、論理的な 1 サイトが
//     2 つの PC に分かれる（サイト数 S に対し 2S-1 個の比較命令が生成された）。
//     サイト数を制御する実験では致命的なので固定する。
#define BP_KERNEL_ATTR __attribute__((noinline, \
    optimize("no-unroll-loops", "reorder-blocks-algorithm=simple")))

BP_KERNEL_ATTR
uint64_t bp_kernel_single(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    (void)ctx;
    uint64_t sink = 0;
    for (uint64_t r = 0; r < reps; r++) {
        for (size_t i = 0; i < n; i++) {
            if (pat[i]) {
                __asm__ volatile("" : "+r"(sink) : :);
                sink += 1;
            } else {
                __asm__ volatile("" : "+r"(sink) : :);
                sink += 3;
            }
        }
    }
    return sink;
}

// ---------------------------------------------------------------- 複数分岐サイト
// 単一の分岐サイトで自分自身の過去と相関させる構成は、ローカル履歴だけで
// 学習できてしまう。グローバル履歴の効果を評価するには、異なる分岐命令
// （異なる PC）の間の相関を使う必要がある。以下のカーネルは各サイトを
// 別の分岐命令として展開する。
//
// 1 サイトぶんの分岐。両分岐先に空 asm を置き cmov 化を防ぐ。
// 各サイトの直後に合流点を強制する。これを省くとコンパイラが後続の比較を
// 分岐先ごとに複製し（tail duplication）、論理的な 1 サイトが 2 つの PC に
// 分かれてしまう。サイト数を制御した実験では致命的なので必ず入れる。
#define BP_SITE(off)                                        \
    if (pat[base + (off)]) {                                \
        __asm__ volatile("" : "+r"(sink) : :);               \
        sink += 1;                                           \
    } else {                                                \
        __asm__ volatile("" : "+r"(sink) : :);               \
        sink += 3;                                           \
    }                                                       \
    __asm__ volatile("" : : : "memory");

#define BP_S2(b)  BP_SITE(b)     BP_SITE((b) + 1)
#define BP_S4(b)  BP_S2(b)       BP_S2((b) + 2)
#define BP_S8(b)  BP_S4(b)       BP_S4((b) + 4)
#define BP_S16(b) BP_S8(b)       BP_S8((b) + 8)
#define BP_S32(b) BP_S16(b)      BP_S16((b) + 16)
#define BP_S64(b) BP_S32(b)      BP_S32((b) + 32)

// 2 サイト。偶数位置が分岐 A、奇数位置が分岐 B。
// B の予測に必要なビットは A の過去の結果であり、B 自身の過去には含まれない。
BP_KERNEL_ATTR
uint64_t bp_kernel_cross(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    (void)ctx;
    uint64_t sink = 0;
    for (uint64_t r = 0; r < reps; r++)
        for (size_t base = 0; base + 2 <= n; base += 2) { BP_S2(0) }
    return sink;
}

#define BP_DEFINE_SITES(S, BODY)                                              \
    BP_KERNEL_ATTR                                                            \
    static uint64_t kernel_sites_##S(const uint8_t *pat, size_t n,            \
                                     uint64_t reps) {                         \
        uint64_t sink = 0;                                                    \
        for (uint64_t r = 0; r < reps; r++)                                   \
            for (size_t base = 0; base + (S) <= n; base += (S)) { BODY }      \
        return sink;                                                          \
    }

BP_DEFINE_SITES(2,  BP_S2(0))
BP_DEFINE_SITES(4,  BP_S4(0))
BP_DEFINE_SITES(8,  BP_S8(0))
BP_DEFINE_SITES(16, BP_S16(0))
BP_DEFINE_SITES(32, BP_S32(0))
BP_DEFINE_SITES(64, BP_S64(0))

// ctx に unsigned* でサイト数（2,4,8,16,32,64）を渡す。
//
// 注意: サイト数を増やすとコード量も増えるため、フロントエンド側の影響と
// aliasing を完全には分離できない。命令数とフロントエンド停止サイクルを
// 併せて観測し、解釈の際にこの点を明記すること。
uint64_t bp_kernel_sites(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    unsigned S = ctx ? *(unsigned *)ctx : 2u;
    switch (S) {
        case 2:  return kernel_sites_2 (pat, n, reps);
        case 4:  return kernel_sites_4 (pat, n, reps);
        case 8:  return kernel_sites_8 (pat, n, reps);
        case 16: return kernel_sites_16(pat, n, reps);
        case 32: return kernel_sites_32(pat, n, reps);
        case 64: return kernel_sites_64(pat, n, reps);
        default: return kernel_sites_2 (pat, n, reps);
    }
}
