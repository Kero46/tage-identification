#!/bin/sh
# 検査自体の回帰検出（env/verify_branch.sh のメタテスト）。
#
# 機械語検査は「壊れているものを NG にする」だけでなく「壊れていないものを
# OK にする」ことも保証しなければ意味がない。片方だけ確かめると、常に NG を
# 返す検査や、常に OK を返す検査を気づかず抱えることになる。
# 実際、旧 verify_branch.sh は cmov しか探しておらず、arm64 では常に「OK」に
# なりうる状態だった（=陰性側の妥当性が未検証だった）。
#
# ここで確認する三点:
#   陽性対照: 本物のカーネルを -O0 でビルドすると分岐が残る → OK になるべき
#   陰性対照: 分岐を持たない（算術で値を使う）カーネル → NG になるべき
#   検査不能: シンボルを strip したバイナリ → NG になるべき（沈黙の見逃し防止）
#
# 使い方: env/verify_selfcheck.sh
DIR=$(dirname "$0")
ROOT="$DIR/.."
CC="${CC:-cc}"
VERIFY="$DIR/verify_branch.sh"

TMPD=$(mktemp -d) || exit 1
trap 'rm -rf "$TMPD"' EXIT
FAIL=0

echo "== 検査自体の回帰検出（verify_branch.sh のメタテスト）=="

# ---------------------------------------------------------------- 陽性対照
# 最適化を切れば if-conversion されないので、全カーネルで分岐が残るはず。
if ! $CC -O0 -g -I "$ROOT/common" -o "$TMPD/pos" \
        "$ROOT/common/kernel.c" "$ROOT/common/rng.c" "$ROOT/common/patterns.c" \
        "$ROOT/common/measure.c" "$ROOT/common/runner.c" \
        "$ROOT/history_length/hist_bench.c" 2>"$TMPD/pos.log"; then
    echo "NG: 陽性対照のビルドに失敗しました（検査不能）"
    tail -5 "$TMPD/pos.log"
    FAIL=1
else
    if sh "$VERIFY" "$TMPD/pos" >"$TMPD/pos.out" 2>&1; then
        echo "OK  陽性対照: -O0 ビルドを OK と判定した"
    else
        echo "NG  陽性対照: 分岐が残っているのに NG と判定した（検査が厳しすぎる）"
        grep -E 'NG:' "$TMPD/pos.out" | head -5
        FAIL=1
    fi
fi

# ---------------------------------------------------------------- 陰性対照
# 分岐を一切持たないカーネル群。パターン値は読むが算術で使うだけなので、
# 「読出 = S だが データ依存分岐 = 0」になる。検査はこれを NG にしなければ
# ならない（今回踏んだ事故と同じ形）。
cat > "$TMPD/plain.c" <<'EOF'
/* 検査の陰性対照。意図的に条件分岐を持たない。測定には使わない。 */
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#define SITE(off) sink += 1 + 2 * (pat[base + (off)] == 0);
#define S2(b)  SITE(b)   SITE((b) + 1)
#define S4(b)  S2(b)     S2((b) + 2)
#define S8(b)  S4(b)     S4((b) + 4)
#define S16(b) S8(b)     S8((b) + 8)
#define S32(b) S16(b)    S16((b) + 16)
#define S64(b) S32(b)    S32((b) + 32)

__attribute__((noinline))
uint64_t bp_kernel_single(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    uint64_t sink = 0; (void)ctx;
    for (uint64_t r = 0; r < reps; r++)
        for (size_t base = 0; base < n; base++) { SITE(0) }
    return sink;
}
__attribute__((noinline))
uint64_t bp_kernel_cross(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    uint64_t sink = 0; (void)ctx;
    for (uint64_t r = 0; r < reps; r++)
        for (size_t base = 0; base + 2 <= n; base += 2) { S2(0) }
    return sink;
}
#define DEF(S, BODY)                                                       \
    __attribute__((noinline))                                              \
    static uint64_t kernel_sites_##S(const uint8_t *pat, size_t n,         \
                                     uint64_t reps) {                      \
        uint64_t sink = 0;                                                 \
        for (uint64_t r = 0; r < reps; r++)                                \
            for (size_t base = 0; base + (S) <= n; base += (S)) { BODY }   \
        return sink;                                                       \
    }
DEF(2,  S2(0))
DEF(4,  S4(0))
DEF(8,  S8(0))
DEF(16, S16(0))
DEF(32, S32(0))
DEF(64, S64(0))
uint64_t bp_kernel_sites(const uint8_t *pat, size_t n, uint64_t reps, void *ctx) {
    unsigned S = ctx ? *(unsigned *)ctx : 2u;
    switch (S) {
        case 4:  return kernel_sites_4 (pat, n, reps);
        case 8:  return kernel_sites_8 (pat, n, reps);
        case 16: return kernel_sites_16(pat, n, reps);
        case 32: return kernel_sites_32(pat, n, reps);
        case 64: return kernel_sites_64(pat, n, reps);
        default: return kernel_sites_2 (pat, n, reps);
    }
}
int main(void) {
    static uint8_t pat[256];
    unsigned S = 8;
    printf("%llu\n", (unsigned long long)(
        bp_kernel_single(pat, 256, 1, 0) + bp_kernel_cross(pat, 256, 1, 0) +
        bp_kernel_sites(pat, 256, 1, &S)));
    return 0;
}
EOF
if ! $CC -O2 -g -o "$TMPD/neg" "$TMPD/plain.c" 2>"$TMPD/neg.log"; then
    echo "NG: 陰性対照のビルドに失敗しました（検査不能）"
    tail -5 "$TMPD/neg.log"
    FAIL=1
else
    if sh "$VERIFY" "$TMPD/neg" >"$TMPD/neg.out" 2>&1; then
        echo "NG  陰性対照: 分岐の無いカーネルを OK と判定した（検査が無効）"
        FAIL=1
    else
        echo "OK  陰性対照: 分岐の無いカーネルを NG と判定した"
    fi
fi

# ---------------------------------------------------------------- 検査不能
# シンボルが無ければ検査できない。これを OK にすると、strip したバイナリが
# 素通りしてしまう。
if [ -f "$TMPD/pos" ] && command -v strip >/dev/null 2>&1; then
    cp "$TMPD/pos" "$TMPD/stripped"
    strip "$TMPD/stripped" 2>/dev/null || strip -x "$TMPD/stripped" 2>/dev/null
    if sh "$VERIFY" "$TMPD/stripped" >"$TMPD/strip.out" 2>&1; then
        echo "NG  検査不能: strip 済みバイナリを OK と判定した（沈黙の見逃し）"
        FAIL=1
    else
        echo "OK  検査不能: strip 済みバイナリを NG と判定した"
    fi
else
    echo "--  検査不能: strip が無いので省略（合否には含めない）"
fi

echo
if [ "$FAIL" = "0" ]; then
    echo "OK: 機械語検査は陽性・陰性の両方向で期待どおり動作しています。"
    exit 0
fi
echo "NG: 機械語検査そのものが壊れています。測定結果を信用しないこと。"
exit 1
