#!/usr/bin/env python3
"""掃引結果の作図。matplotlib が無い環境では端末に簡易グラフを出す。"""
import argparse, csv, sys


def load(path, x, y):
    xs, ys = [], []
    with open(path) as f:
        for r in csv.DictReader(f):
            if r.get(y) in (None, ""):
                continue
            xs.append(float(r[x])); ys.append(float(r[y]))
    return xs, ys


def ascii_plot(xs, ys, w=60, h=18):
    lo, hi = min(ys), max(ys)
    rng = (hi - lo) or 1.0
    print(f"  y: {lo:.4f} .. {hi:.4f}")
    for row in range(h, -1, -1):
        line = ""
        for i in range(len(xs)):
            lvl = round((ys[i] - lo) / rng * h)
            line += "*" if lvl == row else " "
        print(f"  |{line}")
    print("  +" + "-" * len(xs))
    print(f"   x: {xs[0]:g} .. {xs[-1]:g}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("-x", default="p1")
    ap.add_argument("-y", default="diff_ns")
    ap.add_argument("-o", "--out", help="画像として保存(matplotlib 必要)")
    a = ap.parse_args()
    xs, ys = load(a.csv, a.x, a.y)
    if not xs:
        raise SystemExit("データがありません")
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.figure(figsize=(7, 4))
        plt.plot(xs, ys, marker="o", ms=3)
        plt.xlabel(a.x); plt.ylabel(a.y); plt.grid(alpha=.3)
        out = a.out or a.csv.replace(".csv", ".png")
        plt.tight_layout(); plt.savefig(out, dpi=140)
        print(f"{out} を書き出しました")
    except ImportError:
        print("matplotlib が無いため簡易表示します", file=sys.stderr)
        ascii_plot(xs, ys)


if __name__ == "__main__":
    main()
