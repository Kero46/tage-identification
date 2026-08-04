#!/usr/bin/env python3
"""境界（変化点）の検出。低密度モードの掃引から学習可否の境界を機械的に決める。

## なぜ専用ツールが必要か（仕様の落とし穴 20）

`ctx` のような低密度モードでは **`tools/sweep.py` の `significant` 列を根拠にでき
ない。** 標的が (D+1) 要素につき 1 個しかないため、要素あたりの差分が 1/(D+1) に
希釈される。一方 **ノイズ床は希釈されない** — 帰無実験は同じ系列の総時間の
ばらつきを測っているので、標的密度とは無関係の大きさを持つ。

実測（機械 B、patlen 262144、trials 21、K=2）では D=66 で満額学習しているのに
要素あたりの差分は 0.045 ns しかなく、この構成のノイズ床は 0.12 ns だった。
**信号が床の 1/3 しかないので、単点判定は必ず no になる。**

しかも **diff と床を同じ (D+1) 倍しても比は変わらない。** 素朴な正規化では
解決しない。単点の SNR が実際に低いのである。

## ではなぜ境界が見えるのか

**点を床と比べるのではなく、隣り合う点の反復標本どうしを比べている。**
境界の直前と直後で、独立反復が**完全に分離する**（前の最小 > 後の最大）。
これは順序統計量の性質だけを使う検定で、交換可能性の帰無仮説の下での厳密確率は

    片側（符号が理論から決まっている場合）  1 / C(2k, k)
    両側（符号を予測できない場合）          2 / C(2k, k)
    片側: k=3 → 0.050,  k=5 → 0.0040,  k=7 → 0.00029

`histd` / `ctx` / `cross` の距離掃引では**符号が理論から予測されている**
（学習できている側が速い = 標的あたり差分が大きい）ので片側で評価する。
`--direction any`（K*(D) のような導出量で符号を予測できない場合）は両側になる。

**符号が逆向きの完全分離は境界として数えない。** 学習側のほうが遅いという
ありえない向きなので、分離として数えず**異常として報告する**。

## 多重比較（重要）

**掃引の隣接対をすべて走査するので多重比較が入る。** 掃引点が n 個なら
m = n−1 対を検定し、族全体の誤警報確率の上界は概ね `m × p`（Bonferroni）。

    k=3, m=10  →  m·p ≈ 0.50   ← 偶然の分離が 0.5 個期待される
    k=5, m=10  →  m·p ≈ 0.040

**実測でこれが起きた。** `units_probe.log` の S=2（k=3、m=5）で、真の境界
96/104 のほかに 84/90 が偶然分離した（ギャップは勝者の 4%）。
`DOMINANCE` の補助規則がこれを退けたが、**あれは経験的な後処理であって多重比較の
制御ではない。** 制御は反復数で行う。

したがって **境界の確定には片側 k ≥ 5 を要求する**（`MIN_REPS_BOUNDARY`）。
k < 5 でも検出はするが `暫定` とし、族全体の誤警報を制御できていない旨を出す。

## 効果量を併記する理由

完全分離は順序だけを見るので**尺度に依存しない代わりに大きさを無視する。**
「微小だが一貫した分離」と「150 倍の段差」が同じ「分離あり」になる。
判定の強さを誤解させないため、段差の絶対値と反復ばらつきに対する倍率を必ず併記する。

## ドリフトへの頑健性（設計上の利点）

反復ごとに掃引全体を回すので、熱や周波数のドリフトは**掃引内の全点を同じ方向に
動かす。** 隣接点間の分離はドリフトでは生じにくい（両点が一緒に動く）。
単点をノイズ床と比べる方式はドリフトの影響を直接受けるが、この設計は受けない。

**したがってこのツールが行うのは単点の有意判定ではなく変化点の検出である。**
個々の点について統計的有意性を主張しない。主張するのは「隣接する二つの掃引点で
反復標本が完全に分離する位置がある」ことだけである。

## 手順

1. 標的あたり正規化（`diff_ns × 正規化係数`）。係数はモードごとに生成器の構成から
   決まる（`NORM`）。未知のモードは黙って 1 にせずエラーにする
2. 反復を束ねる。各 p1 について反復の値を集める。**個々の反復に独立して境界を
   要求してはならない** — 反復は「同じ p1 点の測り直し」であって独立した掃引では
   なく、単点では暴れる（実測で warmup=512 の D=67 が −0.09 / −1.77 / +1.28）
3. 隣接する各対で完全分離（`min(前) > max(後)`）を検定し、ギャップを記録する
4. 分離が複数あれば、最大のギャップが次点の `DOMINANCE` 倍を超えることを要求する。
   超えなければ段が複数あるということで**判定不能**
5. 分離が無い場合、束ねた曲線に段があるか（`STEP_MIN_RATIO`）で二つを区別する。
   段があるのに分離しない → **判定不能**（反復が食い違っている）。
   段も無い → **境界なし**
6. 結論は三値。**「判定不能」と「境界なし」を区別する**（落とし穴 20 の本体は
   この混同だった）

使い方:
  tools/boundary.py a.csv b.csv c.csv          # 同一条件の反復。モードは .meta.txt から
  tools/boundary.py --mode ctx a.csv b.csv c.csv
  tools/boundary.py --selftest                 # 合成データで自己検査
"""
import argparse
import csv
import math
import os
import statistics
import sys

# ---------------------------------------------------------------- 正規化係数
# 標的あたりに直すための係数。生成器の構成から決まるのでモードごとに定義する。
# **未知のモードを 1 として通さない。** 黙って希釈されたままの値を比較すると
# 低密度モードで境界を見落とす（落とし穴 20 そのもの）。
#
#   ctx   ブロック長 D+1 に標的 1 個            → (D+1)
#   hist  ブロック長 D+1 に標的 1 個            → (D+1)
#   dual  ブロック長 max(D1,D2)+1 に標的 1 個   → (max+1)
#   histd 奇数位置が標的（密度 50%）            → 2
#   cross 奇数位置が分岐 B（密度 50%）          → 2
#   ctxnoise ブロック長 m+D+1 に標的 1 個     → (D + m + 1)
#   raw   既に正規化済み／導出量（K*(D) など）  → 1
#
# alias は差分法を使わないので対象外（CLI が --pair / --control を拒否する）。
#
# `raw` は第2段 (a) の K*(D) 系列のような**導出量**のためにある。K*(D) は差分では
# なく「K 方向の境界位置」なので標的密度の概念が無い。raw を使うときは
# `--direction any` も検討すること（K*(D) の段の向きは理論から決まらない）。
NORM = {
    "ctx":   lambda p1, p2: p1 + 1.0,
    "hist":  lambda p1, p2: p1 + 1.0,
    "dual":  lambda p1, p2: max(p1, p2) + 1.0,
    "histd": lambda p1, p2: 2.0,
    "cross": lambda p1, p2: 2.0,
    # ctxnoise は [雑音 m][識別 1][サフィックス D-1][標的] なのでブロック長は
    # m+D+1。**(D+1) を使うと雑音ぶんの希釈を過小評価する**（D=66, m=20 で
    # 1/87 対 1/67 の違い。仕様 §4.1 の制約 2）。
    "ctxnoise": lambda p1, p2: p1 + p2 + 1.0,
    "raw":   lambda p1, p2: 1.0,
}

# 軸名 → その軸を持つモード、の逆引きは作らない。軸はベンチマークの --axes が
# 出したものを CSV の p1_name / p2_name 列で受け取る（出所を一箇所に保つ）。

# 分離が複数あるときに勝者が次点の何倍のギャップを要求するか。
#
# 根拠: 完全分離の厳密確率は k=3 で 0.10 なので、掃引点が n 個あれば n-1 対を
# 検定することになり、偶然の分離が期待値 0.1(n-1) 個出る。実測でも
# `units_probe.log` の S=2 で、真の境界 96/104（ギャップ +1.37）のほかに
# 84/90 が偶然分離した（ギャップ +0.05、勝者の 4%）。段が本当に一つなら勝者は
# 桁違いに大きくなるので、2 倍は緩い側に置いた保守的な値である。
#
# **勝者が次点を圧倒しないときは「境界なし」ではなく「判定不能」にする。**
# 段が複数あるのか偶然なのかを、この道具だけでは決められないため。
DOMINANCE = 2.0

# 分離が無いときに「判定不能」と「境界なし」を分ける閾値。束ねた曲線の段差を
# 反復のばらつき（純粋な測定ノイズ）で割った値。
#
# 根拠: これは**区別のための診断値**であって有意性の主張ではない。段差が反復
# ばらつきの 2 倍以上あるのに分離しないなら、反復が食い違っているということなので
# 「境界なし」と言ってはいけない。2 倍という値は、段が無い実データ（`--sites 8`、
# だらだら減衰して段が立たない）で 1.5 以下に収まることから置いた。
STEP_MIN_RATIO = 2.0

# 境界を「確定」と述べるために必要な片側反復数。
#
# 根拠: 完全分離の片側厳密確率は k=3 で 0.050、k=5 で 0.0040。掃引点が n 個なら
# m = n−1 対を検定するので、族全体の誤警報の上界は m·p。典型的な m = 5〜10 では
#   k=3 → m·p = 0.25〜0.50  （偶然の分離が半分の確率で出る。**制御できていない**）
#   k=5 → m·p = 0.020〜0.040
# k=3 で実際に偶然の分離が出た（`units_probe.log` の S=2 の 84/90）。
# **k < 5 の判定は「暫定」とし、確定した境界として扱わない。**
MIN_REPS_BOUNDARY = 5


def sep_pvalue(k, two_sided=False):
    """k 反復が完全分離する厳密確率（交換可能性の帰無仮説の下）。

    2k 個の値を 2 群に分ける組み合わせは C(2k,k) 通りで、そのうち完全分離は
    片側で 1 通り、両側で 2 通り。符号が理論から予測できる場合は片側を使う
    （不必要に保守的にしない）。
    """
    return (2.0 if two_sided else 1.0) / math.comb(2 * k, k)


def family_bound(k, m, two_sided=False):
    """族全体の誤警報確率の上界（Bonferroni）。m は検定した隣接対の数。"""
    return min(1.0, m * sep_pvalue(k, two_sided))


def spread(vals):
    if len(vals) < 2:
        return 0.0
    if len(vals) >= 4:
        q = statistics.quantiles(vals, n=4)
        return q[2] - q[0]
    return max(vals) - min(vals)


def norm_factor(mode, p1, p2):
    if mode not in NORM:
        raise SystemExit(
            f"モード '{mode}' の正規化係数が未定義です。tools/boundary.py の NORM に\n"
            "  生成器の構成から決まる係数を追加してください。"
            "未知のモードを 1 として通すと、\n"
            "  低密度モードで境界を見落とします（仕様の落とし穴 20）。")
    return NORM[mode](p1, p2)


def read_mode(path):
    meta = path + ".meta.txt"
    if not os.path.exists(meta):
        return None
    with open(meta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("モード:"):
                return line.split(":", 1)[1].strip()
    return None


def load(path, mode, axis=None):
    """{掃引軸の値: 正規化後の diff} を返す。

    **軸は明示的に決める。自動検出は最後の手段である。**

    以前は「動いている列」を軸として自動検出していた。K 掃引（p1 が定数）の欠陥を
    直した実装だが、**2 軸が同時に動く掃引では曖昧になる。** `ctxnoise` は D と m の
    両方を振るのでこの曖昧さが実際に効く。**黙って片方を選ぶのが最悪の失敗様式**
    なので、2 軸が動いていたらエラーにする。

    軸の決め方（優先順）:
      1. 引数 axis（"p1" / "p2" / 軸名そのもの。CLI の --axis）
      2. CSV の p1_name / p2_name 列（sweep.py が --axes から記録した自己記述）
         + 実際に動いている列。名前が分かっていれば、どちらが軸かを名前で言える
      3. 動いている列の自動検出（fallback）
    """
    rows, names = [], {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            try:
                p1, p2 = float(r["p1"]), float(r.get("p2") or 0)
                d = float(r["diff_ns"])
            except (KeyError, TypeError, ValueError):
                continue
            if r.get("p1_name"):
                names["p1"] = r["p1_name"]
            if r.get("p2_name"):
                names["p2"] = r["p2_name"]
            rows.append((p1, p2, d))
    if not rows:
        return {}, names

    moving = []
    if len({p1 for p1, _, _ in rows}) > 1:
        moving.append("p1")
    if len({p2 for _, p2, _ in rows}) > 1:
        moving.append("p2")

    which = None
    if axis:
        # 軸名そのもので指定された場合は p1_name / p2_name と照合する
        if axis in ("p1", "p2"):
            which = axis
        elif names.get("p1") == axis:
            which = "p1"
        elif names.get("p2") == axis:
            which = "p2"
        else:
            raise SystemExit(
                f"{path}: 軸 '{axis}' が見つかりません。"
                f"この CSV の軸名は p1={names.get('p1', '?')} / "
                f"p2={names.get('p2', '?')} です。")
    elif len(moving) > 1:
        n1 = names.get("p1", "p1")
        n2 = names.get("p2", "p2")
        raise SystemExit(
            f"{path}: 2 軸が同時に動いています（p1={n1}, p2={n2}）。\n"
            "  どちらを軸にするかを --axis で明示してください"
            f"（例: --axis {n1} または --axis {n2}）。\n"
            "  **黙って片方を選ぶことはしません。** 2 次元の掃引は軸を固定して"
            "分割してから渡してください。")
    elif moving:
        which = moving[0]                  # fallback: 動いている列
    else:
        which = "p1"                       # 単点

    out = {}
    for p1, p2, d in rows:
        key = p2 if which == "p2" else p1
        out[key] = d * norm_factor(mode, p1, p2)
    return out, names


def pooled_step(xs, reps):
    """束ねた曲線の最良の単一分割と、その段差 / 反復ばらつき比を返す（診断用）。

    分割の選択は群内の絶対偏差の和（L1）の最小化で行う。段差の倍率を最大にする
    分割を選んではならない（平坦域を 2 点に切るとばらつきの推定値が偶然小さくなり
    そこが必ず勝つ。実際にそれで真の境界 66/67 に対して 65/66 を返した）。

    基準尺度は**反復間のばらつきの代表値**を使う。平坦域そのものの散らばりを使うと、
    境界前が緩やかに減衰している掃引（実測の `--sites 2` / `--sites 4`）で
    傾きがばらつきとして計上され、段差が埋もれる。
    """
    cons = [statistics.median(reps[x]) for x in xs]
    n = len(cons)
    if n < 2:
        return None
    best = None
    for j in range(1, n):
        lo, hi = cons[:j], cons[j:]
        lm, hm = statistics.median(lo), statistics.median(hi)
        cost = sum(abs(v - lm) for v in lo) + sum(abs(v - hm) for v in hi)
        if best is None or cost < best[0]:
            best = (cost, j, lm - hm)
    _, j, step = best
    noise = statistics.median([spread(reps[x]) for x in xs]) or 1e-12
    return {"j": j, "step": step, "noise": noise, "ratio": step / noise,
            "last": xs[j - 1], "first": xs[j]}


def detect(xs, reps, direction="down", multi=False):
    """境界の検出。

    direction:
      "down" 前 > 後 の分離だけを境界とする（距離掃引。符号は理論から決まる）
      "up"   前 < 後 の分離だけを境界とする
      "any"  両向きを境界とする（符号を予測できない導出量。両側評価）

    multi:
      False 段は 1 つと仮定する。複数の分離があれば勝者が次点を DOMINANCE 倍
            圧倒することを要求し、しなければ判定不能（距離掃引はこちら）
      True  **段が複数あることを前提に、分離をすべて段として報告する。**
            K*(D) = Σ C_i (L_i>=D) はテーブル 1 つにつき 1 段できるので、
            段が複数あるのが正常である。単段を仮定すると、模型どおりの階段が
            「複数あって勝者が圧倒しない」で判定不能になる（合成データで実際に
            そうなった）。DOMINANCE の要求は外し、多重比較は m·p で明示する

    返り値の verdict は 境界あり / 境界なし / 判定不能。
    strength は 確定 / 暫定（反復不足で族全体の誤警報を制御できていない）。
    """
    k = min(len(reps[x]) for x in xs)
    diag = pooled_step(xs, reps)
    two_sided = (direction == "any")
    m = max(0, len(xs) - 1)
    stats = {"k": k, "m": m, "two_sided": two_sided,
             "p": sep_pvalue(k, two_sided),
             "family": family_bound(k, m, two_sided)}

    if len(xs) < 2:
        return {"verdict": "判定不能", "reason": f"掃引点が {len(xs)} 個しかない",
                "cands": [], "anomalies": [], "diag": diag, "strength": "—",
                **stats}

    # 隣接対の完全分離。向きも記録する。
    cands, anomalies = [], []
    for i in range(len(xs) - 1):
        before, after = reps[xs[i]], reps[xs[i + 1]]
        down = min(before) - max(after)          # 前が上（学習側が速い）
        up = min(after) - max(before)            # 後が上（ありえない向き）
        if down > 0:
            hit = {"last": xs[i], "first": xs[i + 1], "gap": down,
                   "dir": "down", "before": before, "after": after}
            if direction in ("down", "any"):
                cands.append(hit)
        elif up > 0:
            hit = {"last": xs[i], "first": xs[i + 1], "gap": up,
                   "dir": "up", "before": before, "after": after}
            if direction in ("up", "any"):
                cands.append(hit)
            else:
                # 学習側のほうが遅い向きの完全分離。設計上ありえないので
                # 境界として数えず異常として残す（黙って捨てると気づけない）。
                anomalies.append(hit)
    cands.sort(key=lambda c: -c["gap"])
    base = {"cands": cands, "anomalies": anomalies, "diag": diag, **stats}

    if not cands:
        # 段があるのに分離しない = 反復の食い違いか向きの異常。
        # 段も無い = 本当に境界なし。
        if anomalies:
            pos = ", ".join(f"{c['last']:g}/{c['first']:g}"
                            for c in anomalies)
            return {"verdict": "判定不能", "strength": "—",
                    "reason": f"想定した向き（{direction}）の分離は無いが、"
                              f"**逆向きの完全分離がある**: {pos}。"
                              "学習側のほうが遅いという設計上ありえない向きなので、"
                              "「境界なし」ではなく生成器か測定系の異常を疑うこと",
                    **base}
        if diag and abs(diag["ratio"]) >= STEP_MIN_RATIO:
            return {"verdict": "判定不能", "strength": "—",
                    "reason": f"束ねた曲線には段差がある"
                              f"（反復ばらつきの {abs(diag['ratio']):.1f} 倍、"
                              f"{diag['last']:g}/{diag['first']:g} 付近）のに、"
                              "どの隣接対でも反復が完全分離しない。反復が食い違っている",
                    **base}
        return {"verdict": "境界なし", "strength": "—",
                "reason": "完全分離する隣接対が無く、束ねた曲線にも段差が無い"
                          + (f"（最大でも反復ばらつきの {abs(diag['ratio']):.1f} 倍）"
                             if diag else ""),
                **base}

    # 効果量。完全分離は順序しか見ないので大きさを併記しないと判定の強さを誤解する。
    # K*(D) では落差そのものがテーブル容量 C_i の推定値になる（物理量）。
    def with_effect(c):
        b_med = statistics.median(c["before"])
        a_med = statistics.median(c["after"])
        noise = max(spread(c["before"]), spread(c["after"])) or 1e-12
        return {**c, "effect": b_med - a_med, "effect_med_before": b_med,
                "effect_med_after": a_med, "effect_noise": noise,
                "effect_ratio": (b_med - a_med) / noise}

    if not multi and len(cands) > 1 and \
            cands[0]["gap"] <= DOMINANCE * cands[1]["gap"]:
        pos = ", ".join(f"{c['last']:g}/{c['first']:g}(ギャップ {c['gap']:+.2f})"
                        for c in cands[:3])
        return {"verdict": "判定不能", "strength": "—",
                "reason": f"完全分離が複数あり勝者が圧倒しない: {pos}。"
                          "段が複数あるのか偶然かを決められない"
                          "（段が複数あるのが正常な観測なら --multi を使う）",
                **base}

    # multi では位置順、単段では落差順に並べる
    steps = [with_effect(c) for c in cands]
    steps.sort(key=(lambda c: c["last"]) if multi else (lambda c: -c["gap"]))
    win = steps[0]
    base = {**base, "steps": steps, "multi": multi}

    if k < MIN_REPS_BOUNDARY:
        return {"verdict": "境界あり", "strength": "暫定", "win": win,
                "reason": f"反復 k={k} < {MIN_REPS_BOUNDARY}。"
                          f"m={m} 対を検定しており族全体の誤警報の上界が "
                          f"{stats['family']:.3f} で、制御できていない",
                **base}
    return {"verdict": "境界あり", "strength": "確定", "reason": "", "win": win,
            **base}


def analyse(paths, mode_arg, direction="down", multi=False, axis=None):
    modes, per_file, axis_names = set(), [], {}
    for p in paths:
        mode = mode_arg or read_mode(p)
        if not mode:
            raise SystemExit(
                f"{p}: モードが分かりません。--mode で渡すか .meta.txt を置いてください。")
        modes.add(mode)
        d, nm = load(p, mode, axis)
        axis_names.update(nm)
        per_file.append((p, d))
    if len(modes) > 1:
        raise SystemExit(f"モードが混在しています: {sorted(modes)}。"
                         "同一条件の反復だけを渡してください。")
    mode = modes.pop()

    sets = [set(d) for _, d in per_file]
    common = set.intersection(*sets)
    if not common:
        raise SystemExit("反復間で共通する p1 がありません。")
    dropped = sorted(set.union(*sets) - common)
    xs = sorted(common)
    reps = {x: [d[x] for _, d in per_file] for x in xs}
    r = detect(xs, reps, direction, multi)
    r.update(mode=mode, paths=[p for p, _ in per_file], xs=xs, reps=reps,
             dropped=dropped, nrep=len(per_file), direction=direction,
             multi=multi, axis=axis, axis_names=axis_names)
    return r


def report(r, quiet=False):
    def say(*a):
        if not quiet:
            print(*a)

    say("判定: 変化点の検出（単点の有意性は主張しない）")
    say("基準: 隣接する掃引点で独立反復が完全分離すること。"
        f"{'両側' if r['two_sided'] else '片側'}厳密確率 "
        f"{'2' if r['two_sided'] else '1'}/C(2k,k) = {r['p']:.4f}（k={r['k']}）")
    say(f"多重比較: 隣接対 m={r['m']} を検定 → "
        f"族全体の誤警報の上界 m·p = {r['family']:.3f}（Bonferroni）")
    say("**significant 列は使いません**"
        "（低密度モードでは信号が希釈され床は希釈されないため、"
        "単点判定は必ず no になる）")
    nm = r.get("axis_names") or {}
    axis_shown = r.get("axis") or "自動検出"
    say(f"\nmode={r['mode']}  反復 {r['nrep']} 本  向き={r['direction']}"
        f"  軸={axis_shown}"
        + (f"（p1={nm.get('p1', '?')} / p2={nm.get('p2', '?')}）" if nm else ""))
    for p in r["paths"]:
        say(f"   {os.path.basename(p)}")
    if r["dropped"]:
        say("   注意: 反復間で共通しない p1 を除外しました: "
            + ", ".join(f"{x:g}" for x in r["dropped"]))
    if r["k"] < MIN_REPS_BOUNDARY:
        say(f"   **警告: 反復 k={r['k']} < {MIN_REPS_BOUNDARY}。**"
            f" 族全体の誤警報の上界が {r['family']:.3f} で制御できていません。"
            f"境界を確定と述べるには {MIN_REPS_BOUNDARY} 本以上必要です")

    say("\n標的あたり（反復間の中央値と生の反復値）:")
    for x in r["xs"]:
        raw = "  ".join(f"{u:+.2f}" for u in r["reps"][x])
        say(f"   p1={x:>6g}  中央 {statistics.median(r['reps'][x]):+.3f}"
            f"   反復: {raw}")

    if r["cands"]:
        say("\n完全分離した隣接対（ギャップ = 前の最小 − 後の最大）:")
        for c in r["cands"]:
            say(f"   {c['last']:>6g}/{c['first']:<6g} ギャップ {c['gap']:+.3f}"
                f"  向き {c['dir']}")
    else:
        say("\n完全分離した隣接対: なし")
    if r.get("anomalies"):
        say("**異常: 向きが逆の完全分離**"
            "（学習側のほうが遅い。設計上ありえないので境界に数えていない）:")
        for c in r["anomalies"]:
            say(f"   {c['last']:>6g}/{c['first']:<6g} ギャップ {c['gap']:+.3f}"
                "  生成器か測定系を疑うこと")
    if r["diag"]:
        d = r["diag"]
        say(f"診断（束ねた曲線の最良分割）: {d['last']:g}/{d['first']:g} "
            f"段差 {d['step']:+.3f} = 反復ばらつき({d['noise']:.3f})の "
            f"{d['ratio']:.1f} 倍")

    if r["verdict"] == "境界あり" and r.get("multi"):
        tag = "" if r["strength"] == "確定" else f"（**{r['strength']}**）"
        say(f"\n判定: 段あり{tag} — {len(r['steps'])} 段を検出")
        say("   位置と落差（落差が物理量の推定値）:")
        for st in r["steps"]:
            say(f"     {st['last']:>8g}/{st['first']:<8g} "
                f"落差 {st['effect']:+12.4g}"
                f"（{st['effect_med_before']:.4g} → {st['effect_med_after']:.4g}）"
                f"  ギャップ {st['gap']:+.4g}"
                f"  ばらつきの {st['effect_ratio']:.1f} 倍")
        tot = sum(st["effect"] for st in r["steps"])
        say(f"   落差の合計 {tot:.4g}")
        if r["strength"] != "確定":
            say(f"   → {r['reason']}")
        return (("境界あり" if r["strength"] == "確定" else "境界あり(暫定)"),
                [(st["last"], st["first"]) for st in r["steps"]])
    if r["verdict"] == "境界あり":
        w = r["win"]
        tag = "" if r["strength"] == "確定" else f"（**{r['strength']}**）"
        say(f"\n判定: 境界あり{tag} — 学習できた最大 {w['last']:g}、"
            f"できなかった最小 {w['first']:g}")
        say(f"   全 {r['nrep']} 反復が分離（前の最小 {min(w['before']):+.2f} > "
            f"後の最大 {max(w['after']):+.2f}、ギャップ {w['gap']:+.3f}）")
        # 効果量。完全分離は順序しか見ないので大きさを併記する。
        say(f"   効果量: 中央値の差 {w['effect']:+.3f}"
            f"（{w['effect_med_before']:+.3f} → {w['effect_med_after']:+.3f}）"
            f" = 反復ばらつき({w['effect_noise']:.3f})の "
            f"{w['effect_ratio']:.1f} 倍")
        if r["strength"] != "確定":
            say(f"   → {r['reason']}")
        return ("境界あり" if r["strength"] == "確定"
                else "境界あり(暫定)"), (w["last"], w["first"])
    if r["verdict"] == "境界なし":
        say(f"\n判定: 境界なし — {r['reason']}")
        say("   → 掃引範囲に境界は無い（検出の前提は満たしている）")
        return "境界なし", None
    say(f"\n判定: 判定不能 — {r['reason']}")
    say("   → **「境界なし」ではない。** 条件を変えて測り直すこと")
    return "判定不能", None


# ---------------------------------------------------------------- 自己検査
def _synth(mode, points, seed, extra=None):
    """合成掃引。points は [(p1, 標的あたりの真の水準)]。

    決定的な擬似乱数でノイズを乗せる（Date/random は使わない）。CSV には希釈された
    diff_ns を入れる（ツール側が正規化することを確かめるため）。
    """
    rows = []
    st = seed * 2654435761 % 2147483647
    for i, (p1, level) in enumerate(points):
        st = (st * 1103515245 + 12345) % 2147483648
        noise = ((st / 2147483648.0) - 0.5) * 0.30
        v = level + noise + ((extra[i] if extra else 0.0))
        rows.append({"p1": p1, "p2": 0,
                     "diff_ns": v / NORM[mode](p1, 0.0),
                     "test_iqr": 0.05, "ctl_iqr": 0.05, "significant": "no"})
    return rows


def _write(path, rows, mode):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                           lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(path + ".meta.txt", "w", encoding="utf-8") as fh:
        fh.write(f"モード: {mode}\n")


def selftest():
    """判定ツール自体を検査対象にする（検査を検査する、の原則）。

    段差がある合成データで検出し、無い合成データで検出せず、前提が崩れた場合を
    「境界なし」と区別して「判定不能」にすることを確かめる。
    **閾値や正規化を変えたらここが落ちるはずである。落ちたら実装を疑うこと。**
    """
    import tempfile
    fail = 0

    def check(name, got, want):
        nonlocal fail
        ok = got == want
        print(f"{'OK ' if ok else 'NG '} {name}"
              f"{'' if ok else f': {got!r}（期待 {want!r}）'}")
        if not ok:
            fail = 1

    with tempfile.TemporaryDirectory() as d:
        def run(name, points, nrep=5, extras=None, direction="down"):
            paths = []
            for i in range(nrep):
                p = os.path.join(d, f"{name}{i}.csv")
                _write(p, _synth("ctx", points, i + 1,
                                 extras[i] if extras else None), "ctx")
                paths.append(p)
            return report(analyse(paths, None, direction), quiet=True)

        # (1) 段差あり。機械 B の ctx K=2 を模す。
        step = [(64, 2.5), (65, 2.5), (66, 2.5), (67, 0.0), (68, 0.0), (69, 0.0)]
        v, pos = run("step", step)
        check("段差ありを検出する", v, "境界あり")
        check("境界の位置が正しい", pos, (66.0, 67.0))

        # (2) 平坦（全域で学習できている）→ 境界なし。
        v, _ = run("flat", [(p, 2.5) for p in range(64, 70)])
        check("平坦を境界と誤検出しない", v, "境界なし")

        # (3) 全域 0（全く学習できていない）→ 境界なし。
        v, _ = run("zero", [(p, 0.0) for p in range(64, 70)])
        check("全域 0 を境界と誤検出しない", v, "境界なし")

        # (4) 境界前が緩やかに減衰する掃引。平坦域そのものの散らばりを基準尺度に
        #     すると傾きがばらつきに計上され段差が埋もれる（実測の --sites 2 /
        #     --sites 4 がこの形。旧実装はこれを取り落とした）。
        v, pos = run("slope", [(84, 2.47), (90, 2.07), (96, 1.86),
                               (104, -0.05), (112, -0.09), (120, 0.05)])
        check("境界前が減衰していても検出する", v, "境界あり")
        check("減衰形でも位置が正しい", pos, (96.0, 104.0))

        # (5) 掃引点が 2 個だけの確認測定（実測の warmup 掃引）。反復が分離すれば
        #     検出できる。平坦域の推定を要求する設計だとここが落ちる。
        v, pos = run("two", [(66, 2.7), (67, 0.0)])
        check("2 点でも反復が分離すれば検出する", v, "境界あり")
        check("2 点でも位置が正しい", pos, (66.0, 67.0))

        # (6) 反復が単点で暴れるが束ねれば段が出る形（実測の warmup=512 が
        #     D=67 で -0.09 / -1.77 / +1.28）。反復ごとに独立判定する設計だと落ちる。
        v, pos = run("noisy", [(64, 2.5), (65, 2.5), (66, 2.7), (67, 0.0), (68, 0.0)],
                     nrep=3,
                     extras=[[0, 0, 1.06, -0.09, 0], [0, 0, -0.02, -1.77, 0],
                             [0, 0, 0.01, 1.28, 0]])
        check("反復が暴れても束ねて検出する", v, "境界あり(暫定)")
        check("暴れた形でも位置が正しい", pos, (66.0, 67.0))

        # (7) 1 本だけ境界が違う → 分離しないが段はある → 判定不能。
        #     **「境界なし」と区別することがこのツールの要点。**
        v, _ = run("disagree", step, nrep=3,
                   extras=[[0] * 6, [0] * 6, [0, 0, 0, 2.5, 0, 0]])
        check("反復が食い違えば判定不能", v, "判定不能")

        # (8) **実測データによる回帰検査。** 合成データは「合意曲線を真の水準として
        #     小さなノイズを乗せる」形になるので、実測より必ず綺麗になる。実際に
        #     `--sites 8` を合成で作ったら段が立ってしまい（反復の散らばりが実測の
        #     1/7 になった）「判定不能」を返した。**実測の反復値をそのまま埋め込む。**
        #     出典は history_length/results/{ctx_boundary_K2,units_probe}.log
        #     （標的あたり ns。人手の読みは docs/lab-notebook.md）。
        def run_real(name, data, direction="down"):
            nrep = len(next(iter(data.values())))
            paths = []
            for rep in range(nrep):
                rows = [{"p1": D, "p2": 2, "diff_ns": data[D][rep] / (D + 1.0),
                         "test_iqr": 0.05, "ctl_iqr": 0.05, "significant": "no"}
                        for D in sorted(data)]
                p = os.path.join(d, f"{name}{rep}.csv")
                _write(p, rows, "ctx")
                paths.append(p)
            return report(analyse(paths, None, direction), quiet=True)

        # ctx K=2（99 成立分岐の根拠）。人手の読み: 66/67
        v, pos = run_real("real_ctxK2", {
            64: [1.90, 2.36, 2.67, 1.88, 2.25],
            66: [2.22, 3.17, 2.28, 2.58, 2.54],
            67: [-0.39, 0.22, 0.11, -0.52, 0.35],
            68: [0.39, 0.72, 0.79, 0.35, -0.03],
            70: [0.13, -0.60, -0.54, 0.20, -0.23]})
        check("実測 ctx K=2 で境界あり", v, "境界あり")
        check("実測 ctx K=2 の位置が人手と一致", pos, (66.0, 67.0))

        # 以下 3 件は **k=3 の実測**なので「暫定」になる。多重比較が制御できて
        # いないことをツールが明示すること自体を検査する（k=5 で再確認が必要）。
        # --sites 8（棚上げ案件）。人手の読み: 明確な境界が出ない
        v, _ = run_real("real_sites8", {
            60: [1.34, 1.11, 1.45], 70: [1.21, 0.63, 1.16],
            80: [0.88, 1.58, 0.61], 90: [0.36, 0.40, 0.66],
            100: [1.08, 0.04, 1.12], 110: [0.30, 0.54, 0.29]})
        check("実測 --sites 8 は境界なし", v, "境界なし")

        # --sites 2（境界前が減衰し、偶然の分離も混じる実測形）。人手: 96/104
        v, pos = run_real("real_sites2", {
            84: [2.32, 2.47, 2.51], 90: [2.27, 1.74, 2.07],
            96: [1.33, 1.86, 2.43], 104: [-0.04, -0.65, -0.05],
            112: [0.12, -0.69, -0.09], 120: [0.05, -0.07, 0.45]})
        check("実測 --sites 2 の位置が人手と一致", pos, (96.0, 104.0))
        check("k=3 の実測は暫定と出す", v, "境界あり(暫定)")

        # warmup=512（単点が暴れる実測形）。人手: 66/67
        v, pos = run_real("real_warmup512",
                          {66: [3.76, 2.68, 2.71], 67: [-0.09, -1.77, 1.28]})
        check("実測 warmup=512 の位置が人手と一致", pos, (66.0, 67.0))
        check("k=3 の warmup も暫定", v, "境界あり(暫定)")

        # (9) 偶然の分離が混じっても勝者が圧倒すれば検出する（実測 --sites 2 の
        #     84/90 が偶然分離した。勝者の 4%）。
        v, pos = run("spur", [(84, 2.43), (90, 2.03), (96, 1.86),
                              (104, -0.05), (112, -0.09), (120, 0.05)])
        check("偶然の分離が混じっても勝者を採る", v, "境界あり")
        check("勝者の位置が正しい", pos, (96.0, 104.0))

        # (10) 同程度の段が二つ → 判定不能（「境界なし」ではない）。
        v, _ = run("twostep", [(60, 3.0), (70, 3.0), (80, 1.5),
                               (90, 1.5), (100, 0.0), (110, 0.0)])
        check("同程度の段が二つなら判定不能", v, "判定不能")

        # (11) K 掃引の形（p1 = D 固定、p2 = K が動く）。軸を p1 に固定した実装だと
        #      全行が 1 キーに潰れて「掃引点が 1 個」になる。第2段 (c) がこの形。
        kpaths = []
        for i in range(5):
            rows = []
            st = (i + 1) * 7919
            for K, level in [(2, 2.5), (4, 2.5), (8, 2.5), (16, 0.0), (32, 0.0)]:
                st = (st * 1103515245 + 12345) % 2147483648
                noise = ((st / 2147483648.0) - 0.5) * 0.30
                rows.append({"p1": 63, "p2": K,
                             "diff_ns": (level + noise) / 64.0,
                             "test_iqr": 0.05, "ctl_iqr": 0.05,
                             "significant": "no"})
            p = os.path.join(d, f"ksweep{i}.csv")
            _write(p, rows, "ctx")
            kpaths.append(p)
        v, pos = report(analyse(kpaths, None), quiet=True)
        check("K 掃引（p2 が動く）を扱える", v, "境界あり")
        check("K 掃引の境界位置が正しい", pos, (8.0, 16.0))

        # (12) 正規化と未知モード
        check("ctx の正規化係数は (D+1)",
              abs(NORM["ctx"](66, 0) - 67.0) < 1e-9, True)
        check("cross の正規化係数は 2", abs(NORM["cross"](6, 0) - 2.0) < 1e-9, True)
        try:
            norm_factor("unknown_mode", 1, 0)
            check("未知モードを拒否する", False, True)
        except SystemExit:
            check("未知モードを拒否する", True, True)

        # (12a) --multi: 段が複数あるのが正常な観測（K*(D)）。単段を仮定すると
        #       模型どおりの階段が「勝者が圧倒しない」で判定不能になる。
        #       実際に合成階段でそうなった。
        # 落差を等しくする（1024 と 1024）。単段仮定の DOMINANCE は
        # 勝者が次点の 2 倍を要求するので、等しい段では必ず落ちる。
        stair = [(13, 4096), (19, 4096), (29, 3072), (43, 3072), (66, 2048)]
        paths = []
        for i in range(5):
            rows = [{"p1": D, "p2": 0, "diff_ns": v + (i - 2) * 0.5,
                     "test_iqr": 0, "ctl_iqr": 0, "significant": ""}
                    for D, v in stair]
            p = os.path.join(d, f"stair{i}.csv")
            _write(p, rows, "raw")
            paths.append(p)
        v, _ = report(analyse(paths, None, "down", multi=False), quiet=True)
        check("多段を単段仮定で見ると判定不能", v, "判定不能")
        rm = analyse(paths, None, "down", multi=True)
        v, pos = report(rm, quiet=True)
        check("--multi なら多段を検出する", v, "境界あり")
        check("段の数が正しい", len(rm["steps"]), 2)
        check("段の位置が正しい", pos, [(19.0, 29.0), (43.0, 66.0)])
        check("落差が容量になる",
              [round(st["effect"]) for st in rm["steps"]], [1024, 1024])
        check("落差の合計が両端の差と一致",
              round(sum(st["effect"] for st in rm["steps"])), 2048)
        # 単調増加する系列を down で見れば段は出ない（異常として拾われる）
        rise = [(13, 1024), (19, 1024), (29, 2048), (43, 2048), (66, 4096)]
        paths = []
        for i in range(5):
            rows = [{"p1": D, "p2": 0, "diff_ns": v + (i - 2) * 0.5,
                     "test_iqr": 0, "ctl_iqr": 0, "significant": ""}
                    for D, v in rise]
            p = os.path.join(d, f"rise{i}.csv")
            _write(p, rows, "raw")
            paths.append(p)
        rr3 = analyse(paths, None, "down", multi=True)
        v, _ = report(rr3, quiet=True)
        check("増加する K*(D) は段にしない", v, "判定不能")
        check("増加を異常として報告する", len(rr3["anomalies"]) >= 1, True)

        # (12a2) 軸の明示。**2 軸が同時に動く CSV では黙って片方を選ばない。**
        #        ctxnoise は D と m の両方を振るのでこの曖昧さが実際に効く。
        two = os.path.join(d, "twoaxis.csv")
        with open(two, "w", newline="", encoding="utf-8") as fh:
            fh.write("p1,p1_name,p2,p2_name,diff_ns,test_iqr,ctl_iqr,significant\n")
            for D, m, v in ((20, 0, 2.5), (20, 8, 2.5), (30, 0, 0.1), (30, 8, 0.1)):
                fh.write(f"{D},D,{m},m,{v},0.05,0.05,no\n")
        with open(two + ".meta.txt", "w", encoding="utf-8") as fh:
            fh.write("モード: ctxnoise\n")
        try:
            load(two, "ctxnoise")
            check("2 軸が動く CSV は軸なしで拒否する", False, True)
        except SystemExit as e:
            check("2 軸が動く CSV は軸なしで拒否する", "2 軸" in str(e), True)
        # 軸名で指定できる（p1_name / p2_name 列と照合）
        got, names = load(two, "ctxnoise", "D")
        check("軸名 D で p1 を選べる", sorted(got.keys()), [20.0, 30.0])
        got, _ = load(two, "ctxnoise", "m")
        check("軸名 m で p2 を選べる", sorted(got.keys()), [0.0, 8.0])
        check("軸名を CSV から読める", (names.get("p1"), names.get("p2")), ("D", "m"))
        try:
            load(two, "ctxnoise", "K")
            check("知らない軸名を拒否する", False, True)
        except SystemExit:
            check("知らない軸名を拒否する", True, True)
        # ctxnoise の正規化はブロック長 m+D+1
        check("ctxnoise の正規化係数は m+D+1",
              abs(NORM["ctxnoise"](66, 20) - 87.0) < 1e-9, True)
        check("ctxnoise の m=0 は D+1 に一致",
              abs(NORM["ctxnoise"](66, 0) - NORM["ctx"](66, 0)) < 1e-9, True)

        # (12b) 片側/両側の厳密確率と多重比較の上界
        check("片側 k=3 の厳密確率", abs(sep_pvalue(3) - 1 / 20) < 1e-12, True)
        check("片側 k=5 の厳密確率", abs(sep_pvalue(5) - 1 / 252) < 1e-12, True)
        check("両側は片側の 2 倍",
              abs(sep_pvalue(5, True) - 2 * sep_pvalue(5)) < 1e-12, True)
        check("族全体の上界は m 倍",
              abs(family_bound(5, 10) - 10 / 252) < 1e-12, True)
        check("族全体の上界は 1 で切る", family_bound(3, 100), 1.0)

        # (12c) 向きが逆の完全分離は境界に数えず異常として報告する。
        rev = [(64, 0.0), (65, 0.0), (66, 0.0), (67, 2.5), (68, 2.5)]
        paths = []
        for i in range(5):
            p = os.path.join(d, f"rev{i}.csv")
            _write(p, _synth("ctx", rev, i + 1), "ctx")
            paths.append(p)
        rr = analyse(paths, None, "down")
        v, _ = report(rr, quiet=True)
        # 逆向きの分離は「境界なし」ではなく「判定不能」にする。設計上ありえない
        # 向きなので、平坦だったのと同じ扱いにしてはいけない。
        check("逆向きの分離は判定不能（境界なしにしない）", v, "判定不能")
        check("逆向きの分離を異常として報告する", len(rr["anomalies"]) >= 1, True)
        check("理由に逆向きと書く", "逆向き" in rr["reason"], True)

        # (12d) --direction any なら逆向きも境界になり、両側評価になる。
        rr2 = analyse(paths, None, "any")
        v2, pos2 = report(rr2, quiet=True)
        check("direction=any は逆向きも境界にする", v2, "境界あり")
        check("direction=any の位置", pos2, (66.0, 67.0))
        check("direction=any は両側評価", rr2["two_sided"], True)

        # (12e) raw モード（K*(D) のような導出量）。正規化しない。
        kstar = [(3, 8192.0), (5, 8192.0), (9, 8192.0),
                 (15, 1024.0), (25, 1024.0), (41, 1024.0)]
        paths = []
        for i in range(5):
            rows = [{"p1": D, "p2": 0, "diff_ns": v + (i - 2) * 0.5,
                     "test_iqr": 0, "ctl_iqr": 0, "significant": ""}
                    for D, v in kstar]
            p = os.path.join(d, f"kstar{i}.csv")
            _write(p, rows, "raw")
            paths.append(p)
        v, pos = report(analyse(paths, None, "any"), quiet=True)
        check("raw モードで K*(D) の段を検出する", v, "境界あり")
        check("K*(D) の段の位置", pos, (9.0, 15.0))
        check("raw の正規化係数は 1", abs(NORM["raw"](66, 0) - 1.0) < 1e-9, True)

        # (13) 厳密確率が反復数で正しく下がること
        check("反復が増えると厳密確率が下がる",
              sep_pvalue(7) < sep_pvalue(5) < sep_pvalue(3), True)

    print()
    if fail:
        print("NG: 境界検出ツールの自己検査が落ちました。"
              "閾値を緩めるのではなく実装を直すこと。")
    else:
        print("OK: 境界検出ツールは段差を検出し、平坦を検出せず、"
              "前提の崩れを判定不能として区別します。")
    return fail


def main():
    ap = argparse.ArgumentParser(
        description="低密度モードの掃引から境界（変化点）を検出する")
    ap.add_argument("csv", nargs="*", help="同一条件・異なる種の反復 CSV")
    ap.add_argument("--mode", help="モード（既定は .meta.txt から読む）")
    ap.add_argument("--direction", choices=["down", "up", "any"], default="down",
                    help="段の向き。down=学習側が高い（既定、符号が理論から決まる"
                         "距離掃引）/ up=逆 / any=両向き（K*(D) など導出量。"
                         "両側評価になる）")
    ap.add_argument("--axis",
                    help="掃引軸を明示する（p1 / p2 / 軸名そのもの。例 --axis D）。"
                         "2 軸が同時に動く CSV では必須。自動検出は 1 軸のみ動いて"
                         "いる場合の fallback")
    ap.add_argument("--multi", action="store_true",
                    help="段が複数あることを前提に、分離をすべて段として報告する"
                         "（K*(D) のように 1 テーブル 1 段できる観測用）")
    ap.add_argument("--selftest", action="store_true",
                    help="合成データと実測回帰で自己検査する")
    args = ap.parse_args()
    if args.selftest:
        raise SystemExit(selftest())
    if not args.csv:
        ap.error("CSV を 1 本以上渡してください（反復一致を見るには 3 本以上）")
    verdict, _ = report(analyse(args.csv, args.mode, args.direction,
                               args.multi, args.axis))
    # 確定のみ 0。暫定は 2 にして、呼び出し側が確定と混ぜないようにする。
    raise SystemExit({"境界あり": 0, "境界あり(暫定)": 2}.get(verdict, 1))


if __name__ == "__main__":
    main()
