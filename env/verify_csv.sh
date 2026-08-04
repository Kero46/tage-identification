#!/bin/sh
# CSV 列構成のドリフト検出。
#
# 実際に起きた回帰: カウンタ（instructions / fe_stalls / scaled）を追加した際に
# CSV が 12 列から 15 列になり、位置で列を読んでいた env/calibrate_patlen.sh だけが
# 黙って壊れた（ns_per_patbranch のはずが scaled を読み、校正値が常に 0）。
# 読み手は列名で解決する方式に統一したが、それでも「列名が変わる」「列が消える」
# 変更は検出できないため、見出し行そのものを期待値と突き合わせる。
#
# 期待値の出典は仕様書 §2.3 と common/runner.c の bp_print_header()。
# 列を増減させる場合は、仕様書・bp_print_header・この期待値・tools/sweep.py の
# 参照名を同時に更新すること。ここが落ちたら、まずどれを忘れたかを確認する。
#
# 使い方: env/verify_csv.sh <バイナリ> --mode <モード> [--param <v>]
EXPECT="mode,control,p1,p2,patlen,trial,branches,misses,cycles,instructions,fe_stalls,scaled,miss_per_patbranch,ns_per_patbranch,ipc"

# 解析側（tools/sweep.py・env/*.sh）が名前で参照している列。
# 順序が変わっても壊れないが、消えると壊れるので個別に存在を確かめる。
NEEDED="mode control p1 p2 patlen trial branches misses cycles instructions fe_stalls scaled miss_per_patbranch ns_per_patbranch ipc"

BIN="$1"
[ -n "$BIN" ] || { echo "NG: 使い方: $0 <バイナリ> --mode <モード> [...]"; exit 1; }
shift
[ -x "$BIN" ] || { echo "NG: $BIN がありません。make してください。"; exit 1; }

HDR=$("$BIN" --csv --trials 1 --reps 1 --patlen 512 "$@" 2>/dev/null | head -1)
[ -n "$HDR" ] || { echo "NG: $BIN が CSV 見出しを出力しませんでした（検査不能）。"; exit 1; }

echo "== $BIN の CSV 列構成 =="
if [ "$HDR" != "$EXPECT" ]; then
    echo "NG: 見出しが期待と違います。"
    echo "  期待: $EXPECT"
    echo "  実際: $HDR"
    echo "  列を変えたなら、仕様書 §2.3 / bp_print_header / この期待値 /"
    echo "  tools/sweep.py の参照名をすべて揃えること。"
    exit 1
fi

MISS=""
for c in $NEEDED; do
    case ",$HDR," in
        *",$c,"*) ;;
        *) MISS="$MISS $c" ;;
    esac
done
if [ -n "$MISS" ]; then
    echo "NG: 解析側が参照する列が欠けています:$MISS"
    exit 1
fi

echo "OK: 15 列すべて期待どおり（$(echo "$HDR" | awk -F, '{print NF}') 列）。"
exit 0
