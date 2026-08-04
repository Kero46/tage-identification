# csvcol.awk — CSV の列を「名前」で取り出す
#
# 位置（cut -d, -f<N>）で読むと、列を追加したときにスクリプトだけが黙って
# 壊れる。実際に起きた: カウンタ（instructions / fe_stalls / scaled）の追加で
# CSV が 12 列から 15 列になり、ns_per_patbranch が 12 → 14 列目に移動した。
# env/calibrate_patlen.sh は 12 列目を読み続けたため、以後ずっと scaled（常に 0）を
# 校正値として扱っていた。エラーも出ないので気づけない。
#
# 使い方: awk -f env/csvcol.awk -v COL=ns_per_patbranch < out.csv
# 列が無ければ何も出力せず終了ステータス 2 で終わる（呼び出し側は値が
# 空になるので、必ず失敗として扱われる形にしておくこと）。
NR == 1 {
    sub(/\r$/, "")
    n = split($0, h, ",")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+|[ \t]+$/, "", h[i])
        if (h[i] == COL) c = i
    }
    if (!c) {
        printf "csvcol: 列 %s が見つかりません（見出し: %s）\n", COL, $0 > "/dev/stderr"
        exit 2
    }
    next
}
{
    sub(/\r$/, "")
    split($0, f, ",")
    print f[c]
}
