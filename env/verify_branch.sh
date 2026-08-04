#!/bin/sh
# 計測カーネルの機械語検査。測定前に必ず実行する（make verify）。
#
# 検査する不変条件:
#   §1.2 条件分岐が消えていない
#        パターン値を読む条件分岐が条件選択命令（x86: cmov 系 /
#        arm64: csel・csinc 系）に if-conversion されていないこと。
#        分岐が消えていれば分岐予測を測っていない。
#   §1.3 論理的な 1 サイトが 1 つの PC に対応する
#        ブロック再配置で後続の比較が分岐先ごとに複製されると、サイト数 S に
#        対し比較命令が 2S-1 個生成される。サイト数を制御する実験では致命的。
#        以前は「比較命令数 = S」の照合が仕様書の手打ちコマンドのままで
#        どの make ターゲットにも入っていなかったため、ここに組み込む。
#
# 検査対象は全カーネル。bp_kernel_single だけを見ていた時期に、
# bp_kernel_cross（第1段の主実験）と kernel_sites_*（第2段 (d)）の
# if-conversion を見逃した。
#
# 検査不能（objdump が無い・シンボルが無い・命令を解釈できない）は
# 成功にしない。沈黙の見逃しを作らないため、すべて失敗として扱う。
#
# 使い方: env/verify_branch.sh <バイナリ>
BIN="${1:-history_length/hist_bench}"
DIR=$(dirname "$0")
AWKP="$DIR/branch_check.awk"

[ -f "$BIN" ]  || { echo "NG: $BIN がありません。make してください。"; exit 1; }
[ -f "$AWKP" ] || { echo "NG: $AWKP がありません。"; exit 1; }

OBJDUMP=""
for c in objdump gobjdump llvm-objdump; do
    command -v "$c" >/dev/null 2>&1 && { OBJDUMP="$c"; break; }
done
[ -n "$OBJDUMP" ] || {
    echo "NG: objdump が見つかりません。検査不能を成功として扱わないため失敗させます。"
    exit 1
}

TMP=$(mktemp) || exit 1
trap 'rm -f "$TMP"' EXIT
"$OBJDUMP" -d "$BIN" > "$TMP" 2>/dev/null || { echo "NG: $OBJDUMP -d に失敗しました。"; exit 1; }
[ -s "$TMP" ] || { echo "NG: 逆アセンブル結果が空です。"; exit 1; }

# アーキテクチャ判定。判定できなければ失敗（未知の環境で黙って通さない）。
FMT=$(grep -m1 -i 'file format' "$TMP")
case "$FMT" in
    *arm64*|*aarch64*)        ARCH=arm64 ;;
    *x86-64*|*x86_64*|*i386*) ARCH=x86 ;;
    *) case "$(uname -m)" in
           arm64|aarch64)  ARCH=arm64 ;;
           x86_64|amd64)   ARCH=x86 ;;
           *) echo "NG: アーキテクチャを判定できません（$FMT / $(uname -m)）。"; exit 1 ;;
       esac ;;
esac

echo "== $BIN の機械語検査 =="
echo "   objdump=$OBJDUMP  arch=$ARCH"
printf "%-18s %6s %6s %6s %6s %6s  %s\n" カーネル 期待 読出 分岐 条選 条分岐 判定

# カーネル名と期待サイト数。kernel_sites_* は static だが -g 付きビルドでは
# 局所シンボルとして残る。strip すると検査不能になるため strip しないこと。
FAIL=0
for spec in \
    bp_kernel_single:1 \
    kernel_pad_0:1 \
    kernel_pad_1:1 \
    kernel_pad_2:1 \
    kernel_pad_4:1 \
    bp_kernel_cross:2 \
    kernel_sites_2:2 \
    kernel_sites_4:4 \
    kernel_sites_8:8 \
    kernel_sites_16:16 \
    kernel_sites_32:32 \
    kernel_sites_64:64
do
    SYM=${spec%:*}
    EXP=${spec#*:}
    R=$(awk -v SYM="$SYM" -v ARCH="$ARCH" -v EXPECT="$EXP" -f "$AWKP" "$TMP")
    FOUND=$(echo "$R" | sed -n 's/.*found=\([0-9]*\).*/\1/p')
    SRCS=$(echo  "$R" | sed -n 's/.*srcs=\([0-9]*\).*/\1/p')
    SITES=$(echo "$R" | sed -n 's/.*sites=\([0-9]*\).*/\1/p')
    SEL=$(echo   "$R" | sed -n 's/.*sel=\([0-9]*\).*/\1/p')
    COND=$(echo  "$R" | sed -n 's/.*cond=\([0-9]*\).*/\1/p')
    NG=$(echo    "$R" | sed -n 's/.*ng=//p')

    if [ "$NG" = "-" ] && [ "$FOUND" = "1" ]; then
        printf "%-18s %6s %6s %6s %6s %6s  OK\n" "$SYM" "$EXP" "$SRCS" "$SITES" "$SEL" "$COND"
    else
        printf "%-18s %6s %6s %6s %6s %6s  NG: %s\n" "$SYM" "$EXP" "$SRCS" "$SITES" "$SEL" "$COND" "$NG"
        FAIL=1
        # 診断用に当該カーネルの逆アセンブルを出す
        awk -v S="$SYM" '
            /^[0-9a-fA-F]+[ \t]+<.*>:[ \t]*$/ {
                nm=$0; sub(/^[0-9a-fA-F]+[ \t]+</,"",nm); sub(/>:[ \t]*$/,"",nm); sub(/^_/,"",nm)
                f=(nm==S); if(f) print "    --- " nm " ---"; next
            }
            f && /^[ \t]*$/ { f=0 }
            f { print "    " $0 }' "$TMP" | head -40
    fi
done

echo
if [ "$FAIL" = "0" ]; then
    echo "OK: 全カーネルでデータ依存の条件分岐が期待数どおり残っています。"
    exit 0
fi
cat <<'MSG'
NG: 不変条件が破れています。測定しても分岐予測を測っていません。

  読出 = パターン配列のバイト読み出し数（= 論理サイト数のはず）
  分岐 = その値を使う条件分岐の数（レジスタを突き合わせて数えるため
         ループ制御の分岐は含まない）
  条選 = 条件選択命令の数（cmov / csel / csinc 等。0 でなければ if-conversion）

  読出 > 期待なら §1.3 のブロック複製、条選 > 0 なら §1.2 の分岐消去。
  検査を緩めるのではなく、カーネルの実装を直すこと。
MSG
exit 1
