#!/bin/sh
# 掃引ドライバ（tools/sweep.py）の出力列のドリフト検出。
#
# ベンチマーク側の CSV 列だけを検査していると、解析側の出力列が変わったときに
# 見逃す。今回の事故（カウンタ追加で 12 → 15 列になり、位置で読んでいた
# スクリプトが黙って壊れた）と同じ形なので、掃引ドライバの出力列も照合する。
#
# cross モードだけは派生列 hist_len_lower を持つ。cross の --param は d だが
# 必要なグローバル履歴長は 2d+1 分岐であり、換算を解析者任せにすると解釈を
# 誤るため。他のモードでは p1 の意味が違うので列そのものを出さない。
# その「モードによって列が変わる」仕様自体を検査対象にする。
#
# 使い方: env/verify_sweep_csv.sh
BASE="p1,p2,test_ns,test_iqr,ctl_ns,ctl_iqr,diff_ns,significant,test_mpb,ctl_mpb,diff_mpb,test_ipc,ctl_ipc,test_fe_per_br,ctl_fe_per_br"
CROSS="p1,hist_len_lower,p2,test_ns,test_iqr,ctl_ns,ctl_iqr,diff_ns,significant,test_mpb,ctl_mpb,diff_mpb,test_ipc,ctl_ipc,test_fe_per_br,ctl_fe_per_br"

SWEEP=tools/sweep.py
HIST=history_length/hist_bench
[ -x "$SWEEP" ] || { echo "NG: $SWEEP に実行権限がありません。"; exit 1; }
[ -x "$HIST" ]  || { echo "NG: $HIST がありません。make してください。"; exit 1; }

TMPD=$(mktemp -d) || exit 1
trap 'rm -rf "$TMPD"' EXIT
FAIL=0
SMALL="--patlen 512 --trials 2 --reps 1 --warmup 1"

echo "== tools/sweep.py の出力列 =="

run() {  # $1=モード $2=p1 $3=出力先
    # shellcheck disable=SC2086
    "$SWEEP" "$HIST" --mode "$1" --p1 "$2" $SMALL -o "$3" >/dev/null 2>&1
}

# 派生列を持たないモード
if ! run histd 15 "$TMPD/histd.csv"; then
    echo "NG: histd の掃引に失敗しました（検査不能）"
    FAIL=1
else
    # csv モジュールの既定行終端は CRLF なので \r を落とす
    H=$(head -1 "$TMPD/histd.csv" | tr -d '\r')
    if [ "$H" != "$BASE" ]; then
        echo "NG: histd の列が期待と違います。"
        echo "  期待: $BASE"
        echo "  実際: $H"
        FAIL=1
    else
        echo "OK  histd: 15 列（派生列なし）"
    fi
fi

# cross のみ派生列を持つ
if ! run cross 7 "$TMPD/cross.csv"; then
    echo "NG: cross の掃引に失敗しました（検査不能）"
    FAIL=1
else
    H=$(head -1 "$TMPD/cross.csv" | tr -d '\r')
    if [ "$H" != "$CROSS" ]; then
        echo "NG: cross の列が期待と違います。"
        echo "  期待: $CROSS"
        echo "  実際: $H"
        FAIL=1
    else
        # 値も確認する。d=7 なら 2d+1 = 15。
        V=$(awk -f env/csvcol.awk -v COL=hist_len_lower < "$TMPD/cross.csv" | head -1)
        P=$(awk -f env/csvcol.awk -v COL=p1 < "$TMPD/cross.csv" | head -1)
        EXP=$(awk -v p="$P" 'BEGIN{printf "%g", 2*p+1}')
        if [ "$(awk -v v="$V" 'BEGIN{printf "%g", v}')" != "$EXP" ]; then
            echo "NG: hist_len_lower が 2d+1 になっていません（p1=$P のとき $V、期待 $EXP）"
            FAIL=1
        else
            echo "OK  cross: 16 列（hist_len_lower = 2d+1 = $EXP を確認）"
        fi
    fi
fi

echo
[ "$FAIL" = "0" ] || {
    echo "NG: 掃引ドライバの出力列が仕様から動いています。"
    echo "  列を変えたなら、この期待値と docs/measurement_guide.md の記述も揃えること。"
    exit 1
}
echo "OK: 掃引ドライバの出力列は期待どおりです。"
exit 0
