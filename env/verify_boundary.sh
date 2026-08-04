#!/bin/sh
# 境界検出ツール（tools/boundary.py）の自己検査。
#
# 判定ツールそのものを検査対象にする。低密度モード（ctx / dual）の境界は
# `significant` 列では判定できず（仕様の落とし穴 20）、このツールの判定が
# 第2段 (a)(c) の結論を直接支える。したがってツールが壊れたら結論が壊れる。
#
# 検査内容（詳細は tools/boundary.py の selftest）:
#   - 段差がある合成データで検出し、平坦・全域 0 では検出しない
#   - 前提が崩れた場合（反復の食い違い、点数不足、段が複数）を
#     「境界なし」ではなく「判定不能」として区別する
#   - **実測データの回帰検査**。ctx K=2 / --sites 2 / --sites 8 / warmup=512 の
#     実測反復値を埋め込み、人手の読みと同じ結論になることを固定する
#
# 使い方: env/verify_boundary.sh
BIN=tools/boundary.py
[ -x "$BIN" ] || { echo "NG: $BIN に実行権限がありません。"; exit 1; }

echo "== 境界検出ツールの自己検査 =="
if "$BIN" --selftest; then
    exit 0
fi
echo
echo "NG: 境界検出ツールが回帰しています。"
echo "  閾値（DOMINANCE / STEP_MIN_RATIO）を緩めるのではなく実装を直すこと。"
echo "  実測回帰が落ちた場合は、人手の読み（docs/lab-notebook.md）と"
echo "  突き合わせてどちらが誤りかを判断すること。"
exit 1
