#!/bin/sh
# パターン長の校正。
#
# パターンは反復実行されるため、短いと予測器が系列全体を記憶してしまう。
# 記憶化領域に入ると「純ランダム」系列まで予測可能になり、対照系列が
# 成立しなくなって差分法が壊れる。
#
# 純ランダム系列の時間がパターン長に対して飽和する点を探し、飽和後の長さを使う。
# 必要な長さは予測器の容量に依存するため、機械ごとに校正する必要がある。
#
# 使い方: taskset -c 2 env/calibrate_patlen.sh
BIN="history_length/hist_bench"
DIR=$(dirname "$0")
[ -x "$BIN" ] || { echo "$BIN がありません。make してください。"; exit 1; }

TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT

echo "測定中..."
for L in 10 11 12 13 14 15 16 17 18 19 20; do
    N=$((1 << L))
    # 列は名前で解決する。以前は位置（12 列目）で読んでおり、カウンタ追加で
    # 列がずれた結果 scaled（常に 0）を校正値として扱っていた。
    R=$("$BIN" --mode random --patlen "$N" --trials 5 --reps 8 --csv 2>/dev/null \
        | awk -f "$DIR/csvcol.awk" -v COL=ns_per_patbranch | sort -n \
        | awk '{a[n++]=$1} END{ if (n==0) printf "nan"; else printf "%.3f", a[int(n/2)]}')
    echo "$L $N $R" >> "$TMP"
done

# 全点を取ってから最大値と比較して判定する（逐次比較ではノイズでぶれる）
echo
echo "パターン長        random の ns/分岐   判定"
awk '{v[NR]=$3; l[NR]=$1; n[NR]=$2; if ($3+0>m) m=$3+0} END {
    if (m <= 0) {
        print "NG: 測定値が取れませんでした（校正不能）。"
        print "    バイナリと CSV の列位置を確認すること。"
        exit 1
    }
    for (i=1; i<=NR; i++) {
        r = v[i]/m
        note = (r < 0.90) ? sprintf("記憶化の疑い(飽和値の%.0f%%)", r*100) : "飽和域"
        printf "2^%-2d (%7d KB)  %8.3f            %s\n", l[i], n[i]/1024, v[i], note
    }
    printf "\n飽和値: %.3f ns/分岐\n", m
}' "$TMP"

cat <<'MSG'

判定のしかた:
  「記憶化の疑い」と出た長さは、予測器が系列を部分的に記憶しているため
  使ってはいけない。飽和域に入った最初の長さ以降を使う。
  ただし長すぎるとパターン配列がキャッシュに収まらず、キャッシュミスが
  ノイズとして混入する。飽和直後の長さを選ぶのが望ましい。

  選んだ値は以降の測定で --patlen に明示的に渡し、実験ノートに記録すること。
MSG
