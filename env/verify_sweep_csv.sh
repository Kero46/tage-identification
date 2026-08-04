#!/bin/sh
# 掃引ドライバ（tools/sweep.py）の出力列のドリフト検出。
#
# ベンチマーク側の CSV 列だけを検査していると、解析側の出力列が変わったときに
# 見逃す。今回の事故（カウンタ追加で 12 → 15 列になり、位置で読んでいた
# スクリプトが黙って壊れた）と同じ形なので、掃引ドライバの出力列も照合する。
#
# cross モードだけは派生列 corr_dist_elements を持つ。cross の --param は d だが
# 標的が必要とする相関距離は 2d+1 **パターン要素**であり、換算を解析者任せに
# すると解釈を誤るため。他のモードでは p1 の意味が違うので列そのものを出さない。
# その「モードによって列が変わる」仕様自体を検査対象にする。
#
# この列は**相関距離であって履歴長ではない**（仕様 §3.2 / §3.4）。旧名
# hist_len_lower は「履歴長の下限・分岐数」を含意していて二重に誤っていた。
#
# p1_name / p2_name は軸の自己記述。**--param / --param2 の意味はモードごとに違う**
# （dual: D1/D2、ctx: D/K、ctxnoise: D/m、alias: S/period）ので、CSV を単独で
# 読んだときに何の軸なのかが分かるように名前を残す。位置と慣習で意味が決まって
# いる状態は落とし穴 15（CSV を位置で読む）と同種。出所はベンチマークの --axes。
#
# noise_floor / random_sat / anomaly 列は判定の再現に必要なので列構成に含める。
# random_sat は純ランダム飽和値。意図せず純ランダムを測っていないかの検査基準で、
# 生成器の引数が範囲外・パターン長不足・記憶化を同じ signature で捕まえる。
# noise_floor が無いと significant を CSV だけから再計算できず、判定式を
# 変えたときに過去の結果を再評価できない（v1 → v2 で実際に必要になった）。
#
# 使い方: env/verify_sweep_csv.sh
BASE="p1,p1_name,p2,p2_name,test_ns,test_iqr,ctl_ns,ctl_iqr,diff_ns,noise_floor,random_sat,significant,anomaly,test_mpb,ctl_mpb,diff_mpb,test_ipc,ctl_ipc,test_fe_per_br,ctl_fe_per_br"
CROSS="p1,p1_name,corr_dist_elements,p2,p2_name,test_ns,test_iqr,ctl_ns,ctl_iqr,diff_ns,noise_floor,random_sat,significant,anomaly,test_mpb,ctl_mpb,diff_mpb,test_ipc,ctl_ipc,test_fe_per_br,ctl_fe_per_br"

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
        echo "OK  histd: 20 列（派生列なし）"
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
        # 値も確認する。d=7 なら 2d+1 = 15 要素。
        V=$(awk -f env/csvcol.awk -v COL=corr_dist_elements < "$TMPD/cross.csv" | head -1)
        P=$(awk -f env/csvcol.awk -v COL=p1 < "$TMPD/cross.csv" | head -1)
        EXP=$(awk -v p="$P" 'BEGIN{printf "%g", 2*p+1}')
        if [ "$(awk -v v="$V" 'BEGIN{printf "%g", v}')" != "$EXP" ]; then
            echo "NG: corr_dist_elements が 2d+1 になっていません（p1=$P のとき $V、期待 $EXP）"
            FAIL=1
        else
            echo "OK  cross: 21 列（corr_dist_elements = 2d+1 = $EXP 要素を確認）"
        fi
    fi
fi

# 判定式そのものの回帰検出。列構成が正しくても判定式が v1（abs を取る誤った式）に
# 戻れば結論が壊れる。負の差分が yes になる状態を機械的に落とす。
echo
echo "== 有意性の判定式（sweep.py の judge）=="
if python3 - <<'PY'
import sys, os
sys.path.insert(0, "tools")
import importlib.util
spec = importlib.util.spec_from_file_location("sweep", "tools/sweep.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

fail = 0
def check(name, got, want):
    global fail
    if got != want:
        print(f"NG  {name}: {got!r}（期待 {want!r}）")
        fail = 1
    else:
        print(f"OK  {name}")

nf = 0.05
# 負の差分（test が control より遅い）は有意にしてはならない。v1 の誤りの本体。
sig, anom = m.judge(-0.30, nf, 0.01, 0.01)
check("負の差分は有意でない", sig, "no")
check("負の差分は anomaly に残る", anom.startswith("negative_diff("), True)
# 床は超えるがばらつきの内側の負の差分は、本物の疑いと区別する。
_, anom2 = m.judge(-0.08, nf, 0.30, 0.30)
check("ばらつき内側の負の差分は区別される",
      anom2.startswith("negative_diff_within_spread"), True)
# ノイズ床を超え、ばらつきも小さい正の差分は有意。
check("正で床超え・低ばらつきは有意", m.judge(+0.30, nf, 0.01, 0.02)[0], "yes")
# ばらつきが差より大きい点は落とす（trials=7 の cross d=7 の形）。
check("差より IQR が大きい点は有意でない", m.judge(+0.06, nf, 1.21, 0.58)[0], "no")
# ノイズ床以下は有意でない。
check("床以下は有意でない", m.judge(+0.02, nf, 0.01, 0.01)[0], "no")
# NaN は有意でない（測定不能）。
check("NaN は有意でない", m.judge(float("nan"), nf, 0.01, 0.01)[0], "no")
check("判定規則の版が記録されている", bool(m.SIG_RULE), True)

# --- 意図せぬ純ランダムの検出（random_flag）---
# 校正の根拠は sweep.py の RANDOM_TOL_FRAC のコメント。実測値で固定する。
SAT = 3.40
# ctx の K > 2^(D-1) フォールバック: test_ns が飽和値そのもの → 発火すべき
check("純ランダムそのものを検出する",
      bool(m.random_flag(3.379, 0.12, SAT)), True)
# cross の d>=7（正常に学習できない点）: 飽和値から 0.34 離れる → 発火しない
check("学習できない点では発火しない",
      bool(m.random_flag(3.02, 0.03, SAT)), False)
# 学習できている点 → 発火しない
check("学習できている点では発火しない",
      bool(m.random_flag(1.80, 0.03, SAT)), False)
# 飽和値が測れていない（NaN）→ 発火しない（検査無効）
check("飽和値不明なら発火しない",
      bool(m.random_flag(3.40, 0.03, float("nan"))), False)
# IQR が大きい場合は IQR を許容幅に使う
check("IQR が大きければ許容幅に使う",
      bool(m.random_flag(3.00, 0.50, SAT)), True)
sys.exit(fail)
PY
then
    echo "OK: 判定式は符号とばらつきの両方を見ています（規則の版は .meta.txt）。"
else
    echo "NG: 判定式が回帰しています。負の差分を有意と判定する式に戻してはいけません。"
    FAIL=1
fi

echo
[ "$FAIL" = "0" ] || {
    echo "NG: 掃引ドライバの出力列が仕様から動いています。"
    echo "  列を変えたなら、この期待値と docs/measurement_guide.md の記述も揃えること。"
    exit 1
}
echo "OK: 掃引ドライバの出力列は期待どおりです。"
exit 0
