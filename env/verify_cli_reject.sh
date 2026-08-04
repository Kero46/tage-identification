#!/bin/sh
# CLI が「表現できない引数」を拒否することの検査。
#
# 黙って別のものを測るのが最悪の失敗様式なので、CLI の拒否は文書上の約束ではなく
# 機械化する。実際に踏んだもの:
#
#   ctx の K > 2^(D-1)   生成器が黙って純ランダム系列を返していた。K 掃引すると
#                        上限を超えた点で diff が 0 に落ちるため「容量の限界」に
#                        見える（実測: D=5, K=32 で test_ns が純ランダムの飽和値
#                        3.38）。小さい D ほど早く当たるので、梯子の短履歴側の段が
#                        まるごと偽物になる（仕様 §4.1、落とし穴 24）
#   alias の --pair      差分法が成立しないモードで差分を取る（仕様 §4.2）
#   histd の 偶数 D      i-D が奇数（標的位置）になり構成が壊れる（仕様 §3.2）
#
# **陰性対照も置く**（有効な引数を誤って拒否しないこと）。「常に拒否する」検査を
# 抱えたまま気づかないのを防ぐ（旧 verify_branch.sh が実際にその状態だった）。
#
# 使い方: env/verify_cli_reject.sh
TB=table_structure/table_bench
HB=history_length/hist_bench
[ -x "$TB" ] || { echo "NG: $TB がありません。make してください。"; exit 1; }
[ -x "$HB" ] || { echo "NG: $HB がありません。make してください。"; exit 1; }

FAIL=0
SMALL="--patlen 4096 --trials 1 --reps 1 --warmup 0 --csv"

# 拒否されるべき: 終了コードが非零であること
reject() {  # $1=説明 $2...=コマンド
    desc=$1; shift
    if "$@" >/dev/null 2>&1; then
        echo "NG  $desc: 受け付けてしまった（拒否すべき）"
        FAIL=1
    else
        echo "OK  $desc: 拒否した"
    fi
}
# 受け付けられるべき: 終了コードが 0 であること
accept() {  # $1=説明 $2...=コマンド
    desc=$1; shift
    if "$@" >/dev/null 2>&1; then
        echo "OK  $desc: 受け付けた"
    else
        echo "NG  $desc: 拒否してしまった（受け付けるべき）"
        FAIL=1
    fi
}

echo "== CLI が表現できない引数を拒否するか =="

# ---- ctx の K 上限（K <= 2^(D-1)）----
# shellcheck disable=SC2086
reject "ctx D=5 K=32 (> 2^4)"   "$TB" --mode ctx --param 5  --param2 32   $SMALL
reject "ctx D=5 K=64 (> 2^4)"   "$TB" --mode ctx --param 5  --param2 64   $SMALL
reject "ctx D=3 K=8  (> 2^2)"   "$TB" --mode ctx --param 3  --param2 8    $SMALL
reject "ctx D=9 K=512 (> 2^8)"  "$TB" --mode ctx --param 9  --param2 512  $SMALL
reject "ctx D=1 (識別ビット不可)" "$TB" --mode ctx --param 1  --param2 2    $SMALL
# 陰性対照: 有効域の端はちょうど通ること
accept "ctx D=5 K=16 (= 2^4)"   "$TB" --mode ctx --param 5  --param2 16   $SMALL
accept "ctx D=3 K=4  (= 2^2)"   "$TB" --mode ctx --param 3  --param2 4    $SMALL
accept "ctx D=66 K=8192"        "$TB" --mode ctx --param 66 --param2 8192 $SMALL

# ---- alias は差分法に使えない ----
# shellcheck disable=SC2086
reject "alias --pair"    "$TB" --mode alias --param 8 --param2 8 --pair    $SMALL
reject "alias --control" "$TB" --mode alias --param 8 --param2 8 --control $SMALL
accept "alias (対照なし)" "$TB" --mode alias --param 8 --param2 8          $SMALL

# ---- ctxnoise: --param2 は m（雑音の本数）。K ではない ----
# 文脈数は構成上 2 固定なので K を渡す経路は存在しない。代わりに、K のつもりで
# 大きな値を渡した場合（m が過大 → ブロック長が patlen を超える）を拒否する。
# 黙って純ランダムに近い系列を測ることになるので、落とし穴 24 と同じ失敗様式。
# shellcheck disable=SC2086
reject "ctxnoise D=20 m=8192 (K のつもり)" \
    "$TB" --mode ctxnoise --param 20 --param2 8192 $SMALL
reject "ctxnoise D=0"       "$TB" --mode ctxnoise --param 0  --param2 0 $SMALL
accept "ctxnoise D=20 m=0"  "$TB" --mode ctxnoise --param 20 --param2 0 $SMALL
accept "ctxnoise D=20 m=20" "$TB" --mode ctxnoise --param 20 --param2 20 $SMALL

# ---- 掃引軸の名前（--axes）。CSV/.meta.txt の自己記述の出所 ----
# **--param / --param2 の意味はモードごとに違う。** 出所をここ一箇所に保つので、
# 全モードで軸名が取れることを検査する（取れないと解析側が p1/p2 に退化する）。
echo
echo "== 掃引軸の名前（--axes）=="
axis_is() {  # $1=バイナリ $2=モード $3=期待 p1_name $4=期待 p2_name
    out=$("$1" --mode "$2" --axes 2>/dev/null)
    g1=$(printf '%s\n' "$out" | sed -n 's/^p1_name=//p')
    g2=$(printf '%s\n' "$out" | sed -n 's/^p2_name=//p')
    if [ "$g1" = "$3" ] && [ "$g2" = "$4" ]; then
        echo "OK  $2: p1=$g1 p2=$g2"
    else
        echo "NG  $2: p1=$g1 p2=$g2（期待 p1=$3 p2=$4）"
        FAIL=1
    fi
}
axis_is "$TB" dual     D1 D2
axis_is "$TB" ctx      D  K
axis_is "$TB" ctxnoise D  m
axis_is "$TB" alias    S  period
axis_is "$TB" histd    D  unused
axis_is "$HB" cross    d  unused
axis_is "$HB" histd    D  unused
axis_is "$HB" period   N  unused
echo

# ---- histd は D が奇数 ----
# shellcheck disable=SC2086
reject "histd D=16 (偶数)" "$TB" --mode histd --param 16 $SMALL
accept "histd D=15 (奇数)" "$TB" --mode histd --param 15 $SMALL
reject "histd D=16 (第1段)" "$HB" --mode histd --param 16 $SMALL

echo
[ "$FAIL" = "0" ] || {
    echo "NG: CLI が表現できない引数を受け付けています。"
    echo "  黙って別のものを測るのが最悪の失敗様式です。検査を緩めず実装を直すこと。"
    exit 1
}
echo "OK: CLI は表現できない引数を拒否し、有効な引数は受け付けます。"
exit 0
