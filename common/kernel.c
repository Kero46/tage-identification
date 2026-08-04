// kernel.c — 計測対象のループ(単一の条件分岐)
//
// パターンをデータ(配列)で表現し、単一の条件分岐が要素を読んで分岐する。
// 掃引パラメータを変えてもコードサイズ・分岐密度が変わらないため、
// フロントエンド側の交絡が掃引に伴って変動しない。
//
// 【分岐の生存をコンパイラに依存せず保証する】
// 以前は「両分岐先に空の asm volatile を置く」方式だったが、これは
// if-conversion を防げない。clang は optimize 属性を無視するため
// （-Wunknown-attributes）、全カーネルのパターン依存分岐が条件選択命令
// （arm64: csinc）に潰され、分岐予測を一切測っていない状態になっていた。
// 属性やヒューリスティクスに頼る設計そのものが脆い。
//
// そこで値の読み出しは C で行い、分岐だけを asm goto にする。asm goto は
// コンパイラから見て不透明な終端であり、その分岐を条件選択に変換できない。
// 生成結果は env/verify_branch.sh（make verify）が毎回検査する。
#include "bp.h"

// ---------------------------------------------------------------- 前提の明示
// 沈黙して壊れるより、対応していない組み合わせでは明示的に失敗させる。
#if defined(__clang__)
#  if __clang_major__ < 9
#    error "asm goto には clang 9 以降が必要"
#  endif
#elif defined(__GNUC__)
#  if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 5)
#    error "asm goto には GCC 4.5 以降が必要"
#  endif
#else
#  error "GCC または clang が必要（分岐の生存を asm goto で保証している）"
#endif

// 値が非零なら lbl へ分岐する。ラベルは被演算子の次の番号なので、
// 入力が 1 つのこの形では %l1 がラベルを指す。
#if defined(__x86_64__)
#  define BP_BR_IF_SET(v, lbl) \
      __asm__ goto ("testb %b0, %b0\n\tjnz %l1" : : "r"((unsigned)(v)) : "cc" : lbl)
#elif defined(__aarch64__)
#  define BP_BR_IF_SET(v, lbl) \
      __asm__ goto ("cbnz %w0, %l1" : : "r"((unsigned)(v)) : : lbl)
#else
#  error "x86-64 と aarch64 のみ対応（分岐命令を明示的に書いているため）"
#endif

// 計測カーネルに共通で付ける属性。
//   noinline                        : 呼び出し側に溶けないようにする（両対応）
//   no-unroll-loops                 : 掃引間でコード形状を一定に保つ
//   reorder-blocks-algorithm=simple : ブロック再配置による分岐サイトの複製を防ぐ。
//     既定の再配置では後続の比較が分岐先ごとに複製され、論理的な 1 サイトが
//     2 つの PC に分かれる（サイト数 S に対し 2S-1 個の比較命令が生成された）。
// optimize 属性は GCC 専用で clang は黙って無視する。clang で警告を出さないよう
// 囲むが、GCC 環境での保護は維持する（外すとブロック複製が戻る）。
// なお分岐そのものの生存は属性ではなく asm goto が保証しており、
// 複製が起きていないことは make verify が「読出 = サイト数」で検査する。
#if defined(__GNUC__) && !defined(__clang__)
#define BP_KERNEL_ATTR __attribute__((noinline, \
    optimize("no-unroll-loops", "reorder-blocks-algorithm=simple")))
#else
#define BP_KERNEL_ATTR __attribute__((noinline))
#endif

// ---------------------------------------------------------------- サイト
// 1 サイトぶんの分岐。分岐は asm goto なので if-conversion されない。
// 直後に合流点を強制する（これを省くとコンパイラが後続の処理を分岐先ごとに
// 複製し、論理的な 1 サイトが 2 つの PC に分かれることがある）。
//
// ラベルは関数スコープなので __COUNTER__ で一意化する。二段の間接展開は
// __COUNTER__ を数値に展開してから貼り付けるため（直接 ## すると展開されない）。
// 複合文にしてあるので、並べるときに区切りの ; は不要（BP_S2 等が連結する）。
// 形の要件が二つある。
//
// (1) 分岐 1 本・無条件分岐 0 本にする。
//     素直に「両アームに sink += k を置いて合流」と書くと、clang は asm goto の
//     ラベル側の辺を cold と見なして taken アームを関数外に飛ばす。すると
//     成立率の高いパターン（period N=10 は 90% 成立）では毎回「飛んで戻る」
//     経路を通り、予測が当たっていても遅くなる（実測 1.24 ns/分岐。
//     仕様 §2.4 の参考値 0.35 に対して 3 倍以上）。
//     そこで分岐を「加算を 1 つ飛ばすか否か」に使う。成立で +1、不成立で +3 と
//     意味は同じまま、無条件分岐と関数外ブロックが消える。
//
// (2) sink をサイトごとにレジスタに materialize させる（"+r"(sink)）。
//     これが無いとコンパイラは各サイトの寄与を独立の値として持ち回り、加算を
//     関数末尾まで遅延させる。サイト数が増えるとレジスタが枯れて退避が入り、
//     1 サイトあたりの命令数が S とともに増えた（実測 6 → 8.5 命令）。
//     コード量の差が掃引に乗るため、サイトあたりの形は S に依らず一定に保つ。
#define BP_SITE_IMPL(v, id)                          \
    {                                                \
        BP_BR_IF_SET((v), bp_s_##id);                \
        sink += 2;                                   \
        __asm__ volatile("" : "+r"(sink) : :);       \
      bp_s_##id: ;                                   \
        sink += 1;                                   \
        __asm__ volatile("" : "+r"(sink) : :);       \
        __asm__ volatile("" : : : "memory");         \
    }
#define BP_SITE_V_(v, id) BP_SITE_IMPL(v, id)
#define BP_SITE_V(v)      BP_SITE_V_(v, __COUNTER__)

BP_KERNEL_ATTR
uint64_t bp_kernel_single(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    (void)ctx;
    uint64_t sink = 0;
    for (uint64_t r = 0; r < reps; r++) {
        for (size_t i = 0; i < n; i++) {
            BP_SITE_V(pat[i])
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
// 各サイトの命令列は S に依らず同一（読み出し 1 個 + 分岐 1 個 + 合流）。
// S を変えて変わるのはサイトの個数とループの刻みだけなので、1 サイトあたりの
// コード量が S に伴って変動することはない。
#define BP_SITE(off) BP_SITE_V(pat[base + (off)])

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

// ---------------------------------------------------------------- 詰め物分岐
// 履歴の単位を切り分けるための道具。1 パターン要素につき、常に同じ方向に分岐する
// 条件分岐を P 本挿入する。値は実行時引数なので分岐は畳まれず（asm goto なので
// そもそも畳めない）、方向は呼び出し側が決める。
//
// 詰め物は完全に予測可能なので**予測ミスを増やさず、履歴スロットだけを消費する**。
// 全条件分岐が履歴を消費するなら、相関距離の境界（要素）× (2+P) が一定になる
// （要素あたりの分岐は サイト 1 本 + 詰め物 P 本 + ループ制御 1 本 = 2+P 本）。
//
// サイト数を変える方法（--sites）との違い: あちらは分岐サイト（PC）の数も同時に
// 変わるため、パス履歴や PC 索引の関与と分離できない。詰め物は PC を P 個増やす
// だけで、標的サイトの構成は変えない。
//
// 【設計上の妥協】P 本を「単一の追加 PC」で出すには内側ループが必要で、その
// ループ制御分岐がもう 1 本増えて 2+P の算術が崩れる。そこで P 本を展開して
// 別々の PC に置く。PC 数は 1+P（P≤4 なので最大 5）で、異常が出た S≥8 の領域には
// 遠い。
#define BP_PAD_IMPL(v, id)                           \
    {                                                \
        BP_BR_IF_SET((v), bp_q_##id);                \
        __asm__ volatile("" : : : "memory");         \
      bp_q_##id: ;                                   \
        __asm__ volatile("" : : : "memory");         \
    }
#define BP_PAD_V_(v, id) BP_PAD_IMPL(v, id)
#define BP_PAD(v)        BP_PAD_V_(v, __COUNTER__)

#define BP_PAD0(v)
#define BP_PAD1(v) BP_PAD(v)
#define BP_PAD2(v) BP_PAD(v) BP_PAD(v)
#define BP_PAD4(v) BP_PAD2(v) BP_PAD2(v)

#define BP_DEFINE_PAD(P)                                                      \
    BP_KERNEL_ATTR                                                            \
    static uint64_t kernel_pad_##P(const uint8_t *pat, size_t n,              \
                                   uint64_t reps, unsigned pv) {              \
        uint64_t sink = 0;                                                    \
        (void)pv;                        /* P=0 では使わない */               \
        for (uint64_t r = 0; r < reps; r++)                                    \
            for (size_t i = 0; i < n; i++) {                                  \
                BP_SITE_V(pat[i])                                             \
                BP_PAD##P(pv)                                                 \
            }                                                                 \
        return sink;                                                          \
    }

BP_DEFINE_PAD(0)
BP_DEFINE_PAD(1)
BP_DEFINE_PAD(2)
BP_DEFINE_PAD(4)

// ctx に bp_pad_cfg* を渡す。pads は 0,1,2,4 のみ。dir は詰め物の分岐方向
// （1 で常に成立、0 で常に不成立）。
uint64_t bp_kernel_pad(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    const bp_pad_cfg *c = (const bp_pad_cfg *)ctx;
    unsigned p  = c ? c->pads : 0u;
    unsigned pv = c ? c->dir  : 1u;
    switch (p) {
        case 1:  return kernel_pad_1(pat, n, reps, pv);
        case 2:  return kernel_pad_2(pat, n, reps, pv);
        case 4:  return kernel_pad_4(pat, n, reps, pv);
        default: return kernel_pad_0(pat, n, reps, pv);
    }
}

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
