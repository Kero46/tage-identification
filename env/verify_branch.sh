#!/bin/sh
# 計測カーネル内の条件分岐が cmov に化けていないかの確認。
# 分岐が消えていれば分岐予測を測っていないため、測定前に必ず実行する。
BIN="${1:-history_length/hist_bench}"
[ -f "$BIN" ] || { echo "$BIN がありません。make してください。"; exit 1; }

echo "== $BIN の bp_kernel_single =="
OUT=$(objdump -d "$BIN" | sed -n '/<bp_kernel_single>:/,/ret/p')
echo "$OUT"
echo
if echo "$OUT" | grep -qiE '\bcmov'; then
    echo "NG: cmov が含まれています。条件分岐が消えているため測定は無効です。"
    exit 1
fi
if echo "$OUT" | grep -qiE '\b(jne|je|jnz|jz|b\.ne|b\.eq)\b'; then
    echo "OK: 条件分岐が残っています。"
else
    echo "警告: 条件分岐を確認できませんでした。手動で確認してください。"
    exit 1
fi
