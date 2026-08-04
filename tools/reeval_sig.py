#!/usr/bin/env python3
"""保存済みの掃引 CSV を現行の判定規則（v2）で再評価する。

なぜ必要か: `significant` 列は判定式に依存するので、式を変えた時点で過去の CSV の
判定は旧規則のものになる。v1 は `abs(diff) > noise_floor` で符号を見ておらず、
「test のほうが遅い」点まで有意と判定していた。第1段の結論がその誤判定に
依存していないかを確かめるため、**再測定せずに** 保存済みの
diff_ns / test_iqr / ctl_iqr と、端末ログに残るノイズ床から判定をやり直す。

v1 → v2 で判定が変わりうる向きは一方向しかない。
  - v1 で yes、v2 で no になる: (i) diff が負、(ii) diff <= max(iqr)
  - v1 で no、v2 で yes になることは**ない**（v2 は v1 より狭い）
したがって再評価で境界が動くとしたら「学習できた距離が短くなる」向きだけである。

ノイズ床の入手:
  1. 同名の .log に「ノイズ床（帰無実験）: X」があればそれを使う（確定値）
  2. 無ければ v1 の判定から区間を挟み込む。v1 は abs(diff) > nf なので
     max{|diff| : no} <= nf < min{|diff| : yes}。区間の両端で再評価し、
     判定が区間内で変わらなければ結論はノイズ床の不確かさに依存しない
使い方: tools/reeval_sig.py <csv> [<csv> ...]
"""
import csv
import importlib.util
import os
import re
import sys

_spec = importlib.util.spec_from_file_location(
    "sweep", os.path.join(os.path.dirname(os.path.abspath(__file__)), "sweep.py"))
_sweep = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_sweep)
judge, SIG_RULE = _sweep.judge, _sweep.SIG_RULE


def num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float("nan")


def floor_from_log(path):
    """.log からノイズ床を読む。無ければ None。"""
    log = os.path.splitext(path)[0] + ".log"
    if not os.path.exists(log):
        return None
    with open(log, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = re.search(r"ノイズ床[^:]*:\s*([0-9.]+)", line)
            if m:
                return float(m.group(1))
    return None


def floor_bracket(rows):
    """v1 の判定（abs(diff) > nf）を逆に解いてノイズ床を挟み込む。"""
    lo = 0.0
    hi = float("inf")
    for r in rows:
        d = abs(num(r.get("diff_ns")))
        if d != d:
            continue
        if r.get("significant") == "yes":
            hi = min(hi, d)
        elif r.get("significant") == "no":
            lo = max(lo, d)
    if hi == float("inf"):
        return None
    return lo, hi


def reeval(path):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        print(f"{path}: 空")
        return []
    if not any(r.get("significant") for r in rows):
        print(f"  {os.path.basename(path)}: 対照なしの掃引（判定列が空）。再評価の対象外")
        return []

    stored_nf = "noise_floor" in rows[0] and num(rows[0].get("noise_floor")) == \
        num(rows[0].get("noise_floor"))
    nf = num(rows[0]["noise_floor"]) if stored_nf else floor_from_log(path)
    src = "CSV の noise_floor 列" if stored_nf else "同名 .log"
    bracket = None
    if nf is None:
        bracket = floor_bracket(rows)
        if bracket is None:
            print(f"  {os.path.basename(path)}: ノイズ床が不明で挟み込みもできない。判定不能")
            return []
        # v1 は abs(diff) > nf で yes なので nf ∈ [max{|diff|:no}, min{|diff|:yes})。
        # **上端は区間に含まれない**（そこを代入すると、区間を決めた当の点が
        # diff <= nf になって必ず no に落ちる。最初これで偽の「判定が変わった」を
        # 出した）。代表値は区間に含まれる下端を使い、上端側は開区間の内側を
        # 取って安定性の確認にだけ用いる。
        lo, hi = bracket
        nf = lo
        bracket = (lo, hi - (hi - lo) * 1e-6 if hi > lo else lo)
        src = f"v1 判定からの挟み込み {lo:.4f} <= nf < {hi:.4f}（代表値 {nf:.4f}）"

    changed = []
    anomalies = []
    for r in rows:
        d, ti, ci = (num(r.get("diff_ns")), num(r.get("test_iqr")),
                     num(r.get("ctl_iqr")))
        old = r.get("significant")
        new, anom = judge(d, nf, ti, ci)
        # ノイズ床が区間でしか分からない場合、両端で判定が一致するかも見る
        stable = True
        if bracket:
            hi_new, _ = judge(d, bracket[1], ti, ci)
            stable = (hi_new == new)
        if anom.startswith("negative_diff"):
            anomalies.append((r, d, anom))
        if old != new:
            changed.append((r, d, ti, ci, old, new, anom, stable))
    return [(path, nf, src, rows, changed, anomalies)]


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    print(f"判定規則: {SIG_RULE}"
          "  (diff > noise_floor かつ diff > max(test_iqr, ctl_iqr))")
    print("v1 → v2 で判定が変わる向きは yes → no の一方向のみ。\n")
    total_changed = total_anom = 0
    for path in sys.argv[1:]:
        out = reeval(path)
        for (p, nf, src, rows, changed, anomalies) in out:
            print(f"== {p}")
            print(f"   行数 {len(rows)}  ノイズ床 {nf:.4f}（{src}）")
            if changed:
                total_changed += len(changed)
                print(f"   判定が変わった点 {len(changed)} 件:")
                for (r, d, ti, ci, old, new, anom, stable) in changed:
                    p1 = r.get("p1")
                    p2 = r.get("p2")
                    note = "" if stable else "  ※ノイズ床の区間内で判定が不安定"
                    print(f"     p1={p1:>8} p2={p2:>6} diff={d:+.4f} "
                          f"iqr(t/c)={ti:.4f}/{ci:.4f}  {old} → {new}"
                          f"  [{anom or '理由: 差がばらつき以下'}]{note}")
            else:
                print("   判定が変わった点: なし")
            if anomalies:
                total_anom += len(anomalies)
                print(f"   負の差分の異常 {len(anomalies)} 件"
                      "（test が control よりノイズ床を超えて遅い）:")
                for (r, d, anom) in anomalies:
                    print(f"     p1={r.get('p1'):>8} p2={r.get('p2'):>6} "
                          f"diff={d:+.4f}  {anom}")
            print()
    print(f"合計: 判定が変わった点 {total_changed} 件、"
          f"負の差分の異常 {total_anom} 件")


if __name__ == "__main__":
    main()
