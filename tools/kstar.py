#!/usr/bin/env python3
"""距離ごとの実効容量 K*(D) を求め、D 方向の梯子（平坦域と段）を探す。

## なぜこれを測るのか（仕様 §4.1）

**距離掃引は原理的に梯子を見せない。** `ctx` を K=2 に絞ると容量が律速しないので
どのテーブルもエントリを保持でき、常に最長履歴のテーブルが提供側になる。D を
変えても提供側が変わらず、D が最長履歴長を超えたところで一度だけ崖が出るだけ。
第1段で観測した「D=66 の単一の崖、中間に段なし」は構成から予想される帰結である。

## K*(D) の期待される形（重要。単調減少する階段）

距離 D の相関を捉えられるのは **履歴長 L_i >= D のテーブルすべて**である。TAGE は
最長一致のテーブルが予測を担うが、そのテーブルが容量で溢れてエントリを保持できな
ければ次に長いテーブルが引き継げる。必要な情報は距離 D にあるので、L_i >= D なら
どのテーブルでも捉えられる。したがって距離 D の文脈を保持できる総容量は

    K*(D) ~= Σ C_i      （L_i >= D を満たす i について）

D を上げると L_i を跨いだ時点でそのテーブルが使用可能集合から脱落し、K* が C_i だけ
落ちる。つまり **K*(D) は単調減少する階段**で

    段の位置 = 梯子の各段 L_i
    段の落差 = そのテーブルの容量 C_i

を直接与える。**(a) と (c) は一体の観測である。**

**平坦域を探してはならない。** K* は単調非増加だが、段の間が平坦になる保証はない
（複数の L_i が近接していれば階段は細かくなる）。探すのは**減少の段**である。

**両端の絶対値も意味を持つ。** 最大の D での K* は最長履歴テーブル単独の容量、
最小の D での K* は使用可能な全テーブルの総容量に対応する。**両端の差が階段の
総落差と整合するかを検算する。**

## 判定が二重に入る

1. **K 方向**（各 D で K* を決める）: `--direction down`。K が大きいほど劣化する
   ので符号は理論から決まる
2. **D 方向**（K*(D) の段を決める）: `--mode raw --direction down`。
   **減少方向の片側**で評価する（C_i > 0 なので K* は増えないはず）。
   増加方向の段が出たら異常として報告する — 上限への衝突、測定ノイズ、
   模型の誤りのいずれかである

**落差を必ず記録する。** 段の位置だけでなく落差が容量の推定値であり、
効果量の併記がここでは物理量そのものになる。

## 落とし穴 26 との関係（必ず patlen 対照を行う）

配列内のブロック数による K の上限は n/(D+1) に比例するので、**上限に当たった K* も
単調減少する。本物の階段と人工物が同じ形をする。** 段が見つかったら必ず patlen を
変えた対照実験を行い、**段の位置と落差の両方が patlen に依存しないこと**を確認する。

## K の上限（三つの制約の最小値）

1. 生成器: K <= 2^(D-1)（CLI が拒否する）
2. 配列内のブロック数: 文脈あたり 10 ブロック確保（上記の交絡）
3. 実務上限

このツールは (1)(2) の範囲外を検出したら警告する。

**それぞれで反復 5 本以上と多重比較の扱いが必要。** そのため入力の種を群に分け、
各群で k=5 の K* を確定させ、群ごとの K* を D 方向の反復として使う。

## 落差の量子化（等比 K 掃引の限界）

K を等比（2 のべき）で掃引すると **K\* は 2 のべきに丸められる**ので、落差も
2 のべき単位に量子化される。合成データで検証したところ、真の C=(4096, 2048, 1024)
に対して検出された落差は (2048, 1024) だった（K\*(13)=7168 が 4096 に、
K\*(29)=3072 が 2048 に丸められたため）。

**段の位置は正しく出るが、落差の絶対値は 2 倍程度の誤差を持つ。** 容量の値を
述べるときはこの丸めを明記し、精度が必要なら段の周辺で K を等差掃引し直す。

使い方:
  tools/kstar.py table_structure/results/kstar_D*_s*.csv --group 5
  tools/kstar.py ... --outdir table_structure/results
"""
import argparse
import collections
import importlib.util
import os
import re
import statistics
import sys

_spec = importlib.util.spec_from_file_location(
    "boundary", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "boundary.py"))
_b = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_b)


def parse_name(path):
    """kstar_D<D>_s<seed>.csv から (D, seed) を取る。"""
    m = re.search(r"_D(\d+)_s(\d+)\.csv$", os.path.basename(path))
    if not m:
        return None
    return int(m.group(1)), int(m.group(2))


def check_k_range(D, ks):
    """K <= 2^(D-1) を確かめる。超えた K は生成器が純ランダムに落ちている。"""
    bad = [k for k in ks if k > 2 ** (D - 1)]
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", nargs="+")
    ap.add_argument("--group", type=int, default=5,
                    help="K* を 1 回確定させるのに使う種の本数（既定 5。"
                         "境界確定の要件 k>=5）")
    ap.add_argument("--outdir", help="K*(D) 系列 CSV の出力先")
    args = ap.parse_args()

    by_d = collections.defaultdict(list)
    for p in args.csv:
        got = parse_name(p)
        if not got:
            print(f"注意: 名前から D/seed を取れないので無視します: {p}",
                  file=sys.stderr)
            continue
        D, seed = got
        by_d[D].append((seed, p))

    if not by_d:
        raise SystemExit("入力がありません（kstar_D<D>_s<seed>.csv の形式が必要）")

    print("=" * 72)
    print("第1段階: K 方向の判定 — 各 D で実効容量 K*(D) を求める")
    print(f"  向き down（K が大きいほど劣化。符号は理論から決まるので片側評価）")
    print(f"  種を {args.group} 本ずつの群に分け、各群で k={args.group} の"
          "境界を確定させる")
    print("=" * 72)

    kstar = collections.defaultdict(list)      # D -> [群ごとの K*]
    for D in sorted(by_d):
        files = sorted(by_d[D])
        # ctx の K 掃引なので軸は p2（K）。**明示する**（自動検出に頼らない）
        ks = sorted(_b.load(files[0][1], "ctx", "p2")[0].keys())
        bad = check_k_range(D, ks)
        print(f"\n--- D={D}  K = {', '.join(f'{k:g}' for k in ks)}")
        if bad:
            print(f"  **警告: K = {', '.join(f'{k:g}' for k in bad)} は "
                  f"2^(D-1)={2**(D-1)} を超えており、生成器が純ランダムに"
                  "フォールバックしています。この点の劣化は容量の限界ではありません")
        ngroup = len(files) // args.group
        if ngroup == 0:
            print(f"  判定不能: 種が {len(files)} 本しかなく "
                  f"k={args.group} の群を作れない")
            continue
        for g in range(ngroup):
            grp = [p for _, p in files[g * args.group:(g + 1) * args.group]]
            r = _b.analyse(grp, "ctx", "down", axis="p2")
            v, pos = _b.report(r, quiet=True)
            seeds = [s for s, _ in files[g * args.group:(g + 1) * args.group]]
            if v == "境界あり":
                kstar[D].append(pos[0])        # 劣化前に保てた最大の K
                eff = r["win"]["effect_ratio"]
                print(f"  群{g + 1} (種 {seeds[0]}..{seeds[-1]}): "
                      f"K* = {pos[0]:g}（劣化開始 {pos[1]:g}）"
                      f"  効果量 {r['win']['effect']:+.3f}"
                      f"（ばらつきの {eff:.1f} 倍）  m·p={r['family']:.3f}")
            else:
                print(f"  群{g + 1} (種 {seeds[0]}..{seeds[-1]}): {v}"
                      f" — {r['reason'][:70]}")

    print("\n" + "=" * 72)
    print("第2段階: D 方向の判定 — K*(D) の減少の段を探す")
    print("  向き down（K*(D) = Σ C_i (L_i>=D) は単調非増加。C_i > 0 なので")
    print("  増える段は出ないはず。片側評価にできる）")
    print("  **平坦域は探さない。** 段の間が平坦になる保証はない")
    print("  **落差 = そのテーブルの容量 C_i**")
    print("=" * 72)

    Ds = sorted(d for d in kstar if kstar[d])
    if not Ds:
        print("\nK*(D) が 1 点も求まらなかったので D 方向の判定はできません。")
        raise SystemExit(1)

    print(f"\n{'D':>6}  {'K* の群ごとの値':<34} 中央")
    for D in Ds:
        vals = "  ".join(f"{v:g}" for v in kstar[D])
        print(f"{D:>6}  {vals:<34} {statistics.median(kstar[D]):g}")

    missing = [d for d in sorted(kstar) if not kstar[d]]
    if missing:
        print(f"\n注意: K* が求まらなかった D: "
              f"{', '.join(str(d) for d in missing)}")

    ngroup = min(len(kstar[D]) for D in Ds)
    if ngroup < 2:
        print(f"\n各 D の K* が {ngroup} 個しかないので D 方向の反復判定ができません。")
        raise SystemExit(1)

    # K*(D) 系列を反復ぶん CSV に書き、boundary.py に raw モードでかける
    outdir = args.outdir or "."
    paths = []
    for g in range(ngroup):
        p = os.path.join(outdir, f"kstar_curve_s{g + 1}.csv")
        with open(p, "w", newline="", encoding="utf-8") as fh:
            fh.write("p1,p2,diff_ns,test_iqr,ctl_iqr,significant\n")
            for D in Ds:
                fh.write(f"{D},0,{kstar[D][g]},0,0,\n")
        with open(p + ".meta.txt", "w", encoding="utf-8") as fh:
            fh.write("モード: raw\n")
            fh.write("内容: K*(D)（距離ごとの実効容量）。diff_ns 列は K* の値\n")
            fh.write(f"由来: ctx の K 方向境界を種 {args.group} 本の群 {g + 1} で確定\n")
        paths.append(p)
    print(f"\nK*(D) 系列を {ngroup} 本書き出しました: "
          f"{', '.join(os.path.basename(p) for p in paths)}")

    print()
    r = _b.analyse(paths, "raw", "down", multi=True, axis="p1")
    verdict, pos = _b.report(r)

    print("\n" + "=" * 72)
    print("梯子の読み取り")
    print("=" * 72)
    # 単調性の確認。K*(D) = Σ C_i (L_i>=D) は単調非増加でなければならない。
    med = [statistics.median(kstar[D]) for D in Ds]
    rises = [(Ds[i], Ds[i + 1], med[i + 1] - med[i])
             for i in range(len(Ds) - 1) if med[i + 1] > med[i]]
    if rises:
        print("  **異常: K*(D) が増加している区間があります**"
              "（模型では単調非増加）:")
        for a, b, dv in rises:
            print(f"    D={a:g} → {b:g} で +{dv:g}")
        print("    原因の候補: 配列ブロック数の上限への衝突、測定ノイズ、模型の誤り")
    if r.get("anomalies"):
        print("  **異常: 増加方向の完全分離があります**"
              "（C_i > 0 なら K* は増えないはず）:")
        for c in r["anomalies"]:
            print(f"    D={c['last']:g}/{c['first']:g} ギャップ {c['gap']:+g}")

    if verdict.startswith("境界あり"):
        steps = r["steps"]
        print(f"\n  減少の段を {len(steps)} 段検出しました"
              "（模型ではテーブル 1 つにつき 1 段）:")
        for st in steps:
            print(f"    梯子の段 L_i は {st['last']:g} < L_i <= {st['first']:g}"
                  f"   **容量 C_i ≈ {st['effect']:.4g}**"
                  f"（{st['effect_med_before']:.4g} → {st['effect_med_after']:.4g}）")
        if verdict.endswith("(暫定)"):
            print("    ※ 反復不足のため暫定。確定には D 方向も k>=5 が必要")
        tot = sum(st["effect"] for st in steps)
        print(f"\n  両端の検算:")
        print(f"    最小 D={Ds[0]:g} の K* = {med[0]:.4g}"
              "（使用可能な全テーブルの総容量に対応）")
        print(f"    最大 D={Ds[-1]:g} の K* = {med[-1]:.4g}"
              "（最長履歴テーブル単独の容量に対応）")
        print(f"    両端の差 {med[0] - med[-1]:.4g}  対  落差の合計 {tot:.4g}")
        if abs((med[0] - med[-1]) - tot) > 0.1 * max(abs(tot), 1):
            print("    → **整合しない。** 掃引点の間に検出できていない段があるか、"
                  "上限に当たっているかを疑うこと")
        else:
            print("    → 整合する（検出した段で両端の差を説明できる）")
        print("\n  **これは実効容量であってテーブルのエントリ数ではありません。**")
        print("  タグ幅・置換方針・有用性ビット・複数テーブルへの割り当てが"
              "混ざります（仕様 §4.1 の交絡）")
        print("  **必ず patlen 対照を行うこと。** 配列ブロック数の上限に当たった")
        print("  K* も単調減少するので、本物の階段と人工物が同じ形をします。")
        print("  段の位置と落差の両方が patlen に依存しないことを確認します。")
    elif verdict == "境界なし":
        print("  **減少の段は見つかりませんでした。** K*(D) は掃引範囲で一定です。")
        print("  結論の候補は二つで、この観測だけでは決められません:")
        print("    (1) 使用可能なテーブルが 1 つで説明できる")
        print("    (2) この手法では段を分離できない")
        print("  **段が無いことは測定の失敗ではありません**（仕様 §4.5）。")
    else:
        print("  判定不能。前提が満たせていないので測り直しが必要です。")
    raise SystemExit(0 if verdict.startswith("境界あり") else 1)


if __name__ == "__main__":
    main()
