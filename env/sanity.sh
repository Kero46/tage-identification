#!/bin/sh
# 測定系が分岐予測を実際に測れていることの合否ゲート（make sanity）。
#
# なぜ時間比で判定するのか:
#   条件分岐が条件選択命令に潰される（if-conversion）と、予測可能な系列と
#   純ランダム系列の実行時間差が消える。命令名やアーキテクチャに依存せず
#   検出できるため、objdump による命令検査より頑健な最終防衛線になる。
#   逆アセンブル検査（env/verify_branch.sh）は「なぜ壊れたか」を示し、
#   この比は「壊れているか」を示す。両方を通すこと。
#
# 閾値の根拠:
#   健全な環境（仕様書 §2.4 / lab-notebook の実測）
#     予測可能な周期 約 0.35 ns/分岐 対 純ランダム 約 4.0 ns/分岐 → 比 約 10
#   壊れた環境（macOS/arm64 + clang、2026-08-03 実測。optimize 属性が
#   無視され全サイトが csinc に if-conversion された状態）
#     周期 0.69 対 ランダム 1.01 → 比 1.5
#   閾値 3.0 は健全値の 1/3、壊れた値の 2 倍。両者の間に十分な余裕がある。
#   機械が変わっても「分岐予測が効いていれば数倍以上の差が出る」性質は
#   変わらないため、この水準を機械ごとに調整する必要はない。
MIN_RATIO=3.0

# パターン長。短いと予測器が純ランダム系列まで記憶し、陰性対照が予測可能に
# なって比が縮む（=このゲートが発火する）。それは差分法が壊れている状態なので
# 発火が正しい。校正値を持っているなら PATLEN で渡すこと。
PATLEN="${PATLEN:-65536}"
TRIALS="${TRIALS:-5}"
REPS="${REPS:-16}"

BIN=history_length/hist_bench
DIR=$(dirname "$0")

# 測定条件としての負荷。**合否の判定材料にはしない**（閾値を増やすとゲートが
# 環境に依存してしまう）。ただしコア固定ができない環境では負荷が支配的な誤差要因
# なので、比を読むときの材料として出す。実測で load 20 のとき境界測定が壊れた
# （実験ノート続報 7）。
LOAD=$(uptime 2>/dev/null | sed -n 's/.*load averages*: *//p')
[ -n "$LOAD" ] && echo "load average: ${LOAD}（1/5/15 分。判定には使わない）"
[ -x "$BIN" ] || { echo "NG: $BIN がありません。make してください。"; exit 1; }

# 中央値の ns/分岐。列は名前で解決する（位置で読むと列追加で黙って壊れる）。
# 列が無ければ値が空になり、比が nan になってこのゲートが落ちる。
med_ns() {
    "$BIN" --mode "$1" --param "$2" --patlen "$PATLEN" \
           --trials "$TRIALS" --reps "$REPS" --csv 2>/dev/null \
      | awk -f "$DIR/csvcol.awk" -v COL=ns_per_patbranch | sort -n \
      | awk 'BEGIN{n=0} {v[n++]=$1} END{
              if (n==0) { print "nan"; exit }
              if (n%2) printf "%.4f", v[int(n/2)]
              else     printf "%.4f", (v[n/2-1]+v[n/2])/2 }'
}

echo "== 合否ゲート: 陽性対照（周期）と陰性対照（純ランダム）の時間比 =="
echo "   patlen=$PATLEN trials=$TRIALS reps=$REPS"

POS=$(med_ns period 10)
NEG=$(med_ns random 0)
echo "   陽性対照 period N=10 : $POS ns/分岐"
echo "   陰性対照 random      : $NEG ns/分岐"

RATIO=$(awk -v p="$POS" -v n="$NEG" 'BEGIN{
    if (p=="nan"||n=="nan"||p+0<=0) { print "nan"; exit }
    printf "%.2f", n/p }')
echo "   比 = ${RATIO} （閾値 ${MIN_RATIO} 以上）"
echo

OK=$(awk -v r="$RATIO" -v t="$MIN_RATIO" 'BEGIN{ print (r!="nan" && r+0>=t+0) ? 1 : 0 }')
if [ "$OK" != "1" ]; then
    cat <<MSG
NG: 陽性対照と陰性対照の時間差が小さすぎます（比 ${RATIO} < ${MIN_RATIO}）。
    分岐予測を測れていません。考えられる原因は次の二つです。

    1. 条件分岐が条件選択命令に潰れている（§1.2）
       → make verify で該当カーネルを特定する
    2. パターン長が短く、陰性対照まで記憶されている（§1.6）
       → env/calibrate_patlen.sh で校正し PATLEN= で渡す

    閾値を下げて通すことは不変条件を弱める行為なのでしないこと。
MSG
    exit 1
fi
echo "OK: 比 ${RATIO}。分岐予測の効果が時間に現れています。"

# ノイズ床（帰無実験）。有意判定の基準になるので併せて出す。
echo
echo "== ノイズ床（帰無実験。control 同士なので差は零であるべき）=="
tools/sweep.py "$BIN" --mode histd --p1 15 --patlen "$PATLEN" \
    --trials 7 --reps "$REPS" --null 2>/dev/null \
  || echo "警告: ノイズ床を取得できませんでした（合否には含めません）"
