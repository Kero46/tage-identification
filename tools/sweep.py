#!/usr/bin/env python3
"""掃引ドライバ。

test（相関あり）と control（相関なし）を **同一プロセス内で交互に** 測り、
その差分から標的分岐の予測可能性を取り出す。差分によって、履歴を運ぶ
フィラー分岐が出すミスの寄与が相殺される。

交互測定にする理由: test と control を別プロセスで逐次実行すると、
周波数ドリフトが差分に系統誤差として乗る。

差が本物かを判断するため、ばらつき（四分位範囲）と、帰無実験
（control 同士の比較。差は零であるべき）によるノイズ床を併せて出す。

判定規則は SIG_RULE で版を付け、ノイズ床を CSV の列と .meta.txt に残す。
判定式を変えたときに、過去の CSV がどの規則で判定されたかを後から決められる
ようにするため（v1 → v2 の変更で実際に必要になった）。

例:
  tools/sweep.py history_length/hist_bench --mode histd --p1 3:255:2 \
      --patlen 65536 -o history_length/results/histd.csv
  tools/sweep.py table_structure/table_bench --mode ctx --p1 31 --p2 2:8192:x2
  tools/sweep.py history_length/hist_bench --mode period --p1 2:1024:x2 --no-control
  tools/sweep.py history_length/hist_bench --mode histd --p1 15 --null
"""
import argparse
import csv
import datetime
import os
import platform
import statistics
import subprocess
import sys

NS = "ns_per_patbranch"
MPB = "miss_per_patbranch"
IPC = "ipc"
FE = "fe_stalls"
INS = "instructions"

# 有意性の判定規則の版。CSV だけを見て「どの規則で判定されたか」を決められる
# ようにするため .meta.txt に残す。式を変えたら必ず版を上げる。
#
#   v1  abs(diff) > noise_floor
#       **誤り。** 学習が起きれば test は control より速くなるので
#       diff = ctl_ns - test_ns は正でなければならない。絶対値を取ると
#       「test のほうが遅い」点まで有意と判定する（実測で cross d=8 の
#       diff = -0.0918 が yes になった）。ばらつきも見ていなかった。
#   v2  diff > noise_floor かつ diff > max(test_iqr, ctl_iqr)
#       現行。下の judge() を参照。
SIG_RULE = "v2"


def judge(diff, nf, test_iqr, ctl_iqr):
    """有意性と異常を判定する。返り値は ("yes"|"no", anomaly 文字列)。

    二つの条件を課す。性質が違うので混同しないこと。

    1. diff > nf（**片側条件。根拠がある**）
       学習が起きれば標的分岐の予測が当たり test 側が速くなるので、
       diff = ctl_ns - test_ns は構造上正でなければならない。符号を要求するのは
       物理的な根拠のある片側条件である。ノイズ床は帰無実験（control 同士）で
       得た「差が零であるべき条件での実測の幅」なので、これを超えることを
       求めるのも根拠がある。
    2. diff > max(test_iqr, ctl_iqr)（**経験的条件。統計的検定ではない**）
       ばらつきが大きい点を落とすためのヒューリスティクスにすぎない。
       四分位範囲は中央値の標準誤差ではないので、これは有意水準を持つ検定では
       なく、「転移域でノイズ床の 5〜24 倍まで IQR が膨らんだ点が yes になる」
       という実測の失敗（trials=7 の cross d=7）に対する経験的な防御である。
       信頼区間として引用してはならない。

    **負の差分は捨てない。** ノイズ床を超える負の差分は、単なるノイズではなく
    test 系列が control より難しくなっていることを示す。設計上そうならない
    （§1.9 で成立率まで揃えてある）ので、生成器か測定系の異常の手がかりになる。
    significant を no にするだけで消してしまうと気づけないため anomaly に残す。
    """
    if diff != diff or nf != nf:          # NaN
        return "no", ""
    spread = max(test_iqr if test_iqr == test_iqr else 0.0,
                 ctl_iqr if ctl_iqr == ctl_iqr else 0.0)
    anomaly = ""
    if -diff > nf:
        # 負の差分をノイズ床だけで拾うと、ばらつきの裾を引いた点まで異常になる
        # （実測で cross d=61 の -0.0373 が床の 1.9 倍で立ったが、IQR 0.07 の
        # ほうが大きく、外部負荷で説明できる）。床を超えたことは記録しつつ、
        # ばらつきも超えたものだけを本物の疑いとして区別する。
        if -diff > spread:
            anomaly = f"negative_diff({diff:+.4f}<-nf={nf:.4f},>iqr={spread:.4f})"
        else:
            anomaly = f"negative_diff_within_spread({diff:+.4f}<-nf={nf:.4f}," \
                      f"<=iqr={spread:.4f})"
    if diff <= nf:
        return "no", anomaly
    if diff <= spread:
        return "no", anomaly or f"noisy(diff={diff:+.4f}<=iqr={spread:.4f})"
    return "yes", anomaly


def parse_range(spec):
    """'3:255:2' → 等差, '2:8192:x2' → 等比, '15' → 単一値"""
    if ":" not in spec:
        return [float(spec)]
    parts = spec.split(":")
    lo, hi = float(parts[0]), float(parts[1])
    step = parts[2] if len(parts) > 2 else "1"
    out = []
    if step.startswith("x"):
        mul, v = float(step[1:]), lo
        while v <= hi:
            out.append(v)
            v *= mul
    else:
        st, v = float(step), lo
        while v <= hi:
            out.append(v)
            v += st
    return out


def invoke(binary, mode, p1, p2, args, pair=False, control=False):
    """バイナリを1回起動し、行のリスト（dict）を返す。"""
    cmd = [binary, "--mode", mode, "--param", f"{p1:g}", "--param2", f"{p2:g}",
           "--trials", str(args.trials), "--reps", str(args.reps),
           "--warmup", str(args.warmup), "--seed", str(args.seed), "--csv"]
    if args.patlen:
        cmd += ["--patlen", str(args.patlen)]
    if args.frontend:
        cmd.append("--frontend")
    if getattr(args, "sites", None):
        cmd += ["--sites", str(args.sites)]
    if getattr(args, "pad", None) is not None:
        cmd += ["--pad", str(args.pad), "--pad-dir", str(args.pad_dir)]
    if pair:
        cmd.append("--pair")
    if control:
        cmd.append("--control")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        raise SystemExit(f"失敗: {' '.join(cmd)}")
    rows = list(csv.DictReader(r.stdout.splitlines()))
    if any(row.get("scaled") == "1" for row in rows):
        print("警告: カウンタが多重化で縮尺されています。misses は信頼できません。",
              file=sys.stderr)
    return rows


def stats(vals):
    """中央値と四分位範囲を返す。"""
    if not vals:
        return float("nan"), float("nan")
    med = statistics.median(vals)
    if len(vals) >= 4:
        q = statistics.quantiles(vals, n=4)
        iqr = q[2] - q[0]
    else:
        iqr = max(vals) - min(vals)
    return med, iqr


def fe_per_branch(rows):
    """フロントエンド停止サイクルを分岐あたりに正規化する。

    サイト数を変える実験ではコード量も変わるため、劣化がフロントエンド由来か
    予測器由来かを切り分ける手がかりになる。--frontend 未指定なら 0。
    """
    vals = []
    for r in rows:
        fe = float(r.get(FE, 0) or 0)
        br = float(r.get("branches", 0) or 0)
        if br > 0:
            vals.append(fe / br)
    return statistics.median(vals) if vals else 0.0


def split_pair(rows, key):
    t = [float(r[key]) for r in rows if r["control"] == "0"]
    c = [float(r[key]) for r in rows if r["control"] == "1"]
    return t, c


# 純ランダム飽和値との区別に使う許容幅。test_ns が「飽和値 − 許容幅」以上なら
# 意図せず純ランダムを測っている可能性を警告する。
#
# 校正: 実測で、意図せぬ純ランダム（ctx の K > 2^(D-1) のフォールバック）は
# test_ns が飽和値そのもの（3.38 対 3.40）になった。一方、正常に「学習できない」
# 点（cross の d>=7）は test_ns が 3.02〜3.06 で飽和値 3.36〜3.40 から 0.34 離れる
# （IQR 0.03 の 10 倍以上）。飽和値の 3% ≈ 0.10 はこの二つの間に十分な余裕をもって
# 入る。IQR がそれより大きいときは IQR を使う。
RANDOM_TOL_FRAC = 0.03


# 測定中に他プロセスが CPU を食っていると反復のばらつきが跳ね上がる。実測で
# load average 3.7 のとき安定していた ctx K=2 の境界測定が、load 20 では
# 「境界なし」になった（標的あたりの値が ±19 まで暴れた。実験ノート続報 7）。
# **コア固定ができない環境では負荷が支配的な誤差要因**なので、測定条件として記録する。
LOAD_WARN = 5.0


def load_avg():
    """1/5/15 分の load average。取得できない環境では None。"""
    try:
        return os.getloadavg()
    except (OSError, AttributeError):
        return None


def fmt_load(la):
    return "取得不可" if la is None else " ".join(f"{v:.2f}" for v in la)


def query_axes(binary, mode):
    """ベンチマークに軸名を問い合わせる（`--axes`）。

    **--param / --param2 の意味はモードごとに違う**（dual: D1/D2、ctx: D/K、
    ctxnoise: D/m、alias: S/period）。位置と慣習で意味が決まっている状態は
    落とし穴 15（CSV を位置で読む）と同種なので、名前を記録して自己記述にする。

    **軸の意味の出所はベンチマーク側（C コード）に一箇所だけ置く。** ここに同じ表を
    持つと必ずずれるので問い合わせる。取得できなければ p1/p2 のままにする。
    """
    r = subprocess.run([binary, "--mode", mode, "--axes"],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return "p1", "p2"
    got = {}
    for line in r.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            got[k.strip()] = v.strip()
    return got.get("p1_name", "p1"), got.get("p2_name", "p2")


def random_saturation(binary, args):
    """この構成での「純ランダム系列」の ns/要素を測る。

    なぜ必要か: 生成器の引数が範囲外・パターン長が不足・記憶化領域など複数の失敗が
    **「気づかないまま純ランダムを測っている」**という同じ signature で現れる。
    実際に `ctx` の K > 2^(D-1) で生成器が黙って純ランダムを返していた（掃引に
    偽の劣化が出て容量の限界に見えた）。飽和値と比べればこの失敗様式一般を捕まえられる。

    `histd --control` は全位置が独立乱数になるので、これが純ランダムの基準線に
    なる（`calibration_sites.txt` と同じ手法）。**同じ patlen / sites / pad で測る**
    ので、それらを変えても基準線が追随する。

    測るのは 1 回だけ（trials 5, reps 8）。掃引全体に対して無視できる費用である。
    """
    cmd = [binary, "--mode", "histd", "--param", "15", "--control",
           "--trials", "5", "--reps", "8", "--csv"]
    if args.patlen:
        cmd += ["--patlen", str(args.patlen)]
    if getattr(args, "sites", None):
        cmd += ["--sites", str(args.sites)]
    if getattr(args, "pad", None) is not None:
        cmd += ["--pad", str(args.pad), "--pad-dir", str(args.pad_dir)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        return float("nan")
    vals = [float(row[NS]) for row in csv.DictReader(r.stdout.splitlines())
            if row.get(NS)]
    return statistics.median(vals) if vals else float("nan")


def random_flag(test_ns, test_iqr, sat):
    """test_ns が純ランダム飽和値と区別できないか。

    **「有意でないこと」の検査ではない。** 学習できていない点で test_ns が対照に
    近づくのは正常である。ここで見ているのは「test 系列そのものが純ランダムに
    なってしまっていないか」であり、生成器の破綻・パターン長不足・記憶化などを
    同じ signature で捕まえる。
    """
    if test_ns != test_ns or sat != sat or sat <= 0:
        return ""
    tol = max(test_iqr if test_iqr == test_iqr else 0.0, RANDOM_TOL_FRAC * sat)
    if test_ns >= sat - tol:
        return (f"test_ns({test_ns:.4f}) が純ランダム飽和値({sat:.4f})と"
                f"区別できない(許容 {tol:.4f})")
    return ""


def noise_floor(binary, mode, p1, p2, args):
    """帰無実験: control 系列を繰り返し測り、試行間の差からノイズ床を求める。

    相関の無い系列を同じ条件で測っているので、試行間の差は零であるべき。
    ここで得られる幅より小さい差は有意と見なせない。
    """
    rows = invoke(binary, mode, p1, p2, args, control=True)
    v = [float(r[NS]) for r in rows]
    if len(v) < 2:
        return float("nan")
    d = [abs(v[i + 1] - v[i]) for i in range(len(v) - 1)]
    return statistics.median(d)


def write_meta(path, args, binary, nf=float("nan"), sat=float("nan"),
               p1_name="p1", p2_name="p2", load_start=None):
    """再現に必要な環境情報と、判定に使った規則・ノイズ床を残す。

    判定規則の版とノイズ床を残す理由: significant 列は判定式に依存するので、
    式を変えた時点で「この CSV はどちらの規則で判定されたのか」が決められないと
    過去の結果と比較できなくなる。実際に v1（abs を取る誤った式）から v2 へ
    変更する必要が生じ、既存の CSV を再評価する作業が発生した。
    """
    meta = path + ".meta.txt"
    cpu = ""
    try:
        with open("/proc/cpuinfo") as fh:
            for line in fh:
                if line.startswith("model name"):
                    cpu = line.split(":", 1)[1].strip()
                    break
    except OSError:
        pass

    def read(p):
        try:
            with open(p) as fh:
                return fh.read().strip()
        except OSError:
            return "?"

    with open(meta, "w") as fh:
        fh.write(f"日時: {datetime.datetime.now().isoformat(timespec='seconds')}\n")
        fh.write(f"CPU: {cpu}\n")
        fh.write(f"プラットフォーム: {platform.platform()}\n")
        fh.write(f"バイナリ: {binary}\n")
        fh.write(f"モード: {args.mode}\n")
        # 軸名。--param / --param2 の意味はモードごとに違うので必ず残す。
        fh.write(f"p1_name={p1_name}  p2_name={p2_name}\n")
        fh.write(f"p1={args.p1} p2={args.p2} patlen={args.patlen or '既定'}\n")
        fh.write(f"trials={args.trials} reps={args.reps} "
                 f"warmup={args.warmup} seed={args.seed}\n")
        fh.write(f"sites={args.sites or '単一'} "
                 f"pad={'なし' if args.pad is None else args.pad}"
                 f"(dir={args.pad_dir})\n")
        # 判定の再現に必要な情報。sig_rule は SIG_RULE の版、noise_floor は
        # 帰無実験で得た値（CSV の noise_floor 列と同じ）。
        fh.write(f"sig_rule={SIG_RULE}"
                 "  (v2: diff > noise_floor かつ diff > max(test_iqr, ctl_iqr)。"
                 "符号は片側条件として根拠あり、IQR 条件は経験的)\n")
        fh.write(f"noise_floor={'測定せず' if nf != nf else f'{nf:.4f}'}"
                 " ns/分岐（帰無実験）\n")
        fh.write(f"random_sat={'測定せず' if sat != sat else f'{sat:.4f}'}"
                 " ns/要素（純ランダム飽和値。意図せぬ純ランダムの検査基準）\n")
        # 測定条件としての負荷。コア固定ができない環境では支配的な誤差要因なので、
        # 測定の前後で記録する（片方だけでは掃引中に負荷が変わった場合を見逃す）。
        fh.write(f"load_avg_start={fmt_load(load_start)}"
                 f"  load_avg_end={fmt_load(load_avg())}"
                 f"（1/5/15 分。{LOAD_WARN} 超で警告。判定には使わない）\n")
        fh.write("ガバナ: "
                 f"{read('/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor')}\n")
        fh.write(f"no_turbo: {read('/sys/devices/system/cpu/intel_pstate/no_turbo')}\n")
        fh.write(f"perf_event_paranoid: {read('/proc/sys/kernel/perf_event_paranoid')}\n")
    print(f"{meta} に環境情報を記録しました")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("binary")
    ap.add_argument("--mode", required=True)
    ap.add_argument("--p1", default="0", help="値 または lo:hi:step / lo:hi:xN")
    ap.add_argument("--p2", default="0")
    ap.add_argument("--no-control", action="store_true",
                    help="対照系列を取らない（period, bias, random 用）")
    ap.add_argument("--null", action="store_true",
                    help="帰無実験のみ実行してノイズ床を出す")
    ap.add_argument("--no-random-check", action="store_true",
                    help="純ランダム飽和値との比較を省く（既定は行う）。"
                         "意図せず純ランダムを測っていないかの検査なので、"
                         "省く理由が無ければ付けないこと")
    ap.add_argument("--patlen", type=int, default=None,
                    help="パターン長（分岐数）。env/calibrate_patlen.sh で校正した値")
    ap.add_argument("--frontend", action="store_true",
                    help="フロントエンド停止サイクルも取る（多重化に注意）。"
                         "サイト数を変える実験では併用を推奨")
    ap.add_argument("--sites", type=int, default=None,
                    help="パターン要素を S 個の分岐サイトに分散して実行する"
                         "（2,4,8,16,32,64）。履歴の単位を確かめる実験用")
    ap.add_argument("--pad", type=int, default=None,
                    help="1 要素あたりの予測可能な詰め物分岐の本数（0,1,2,4）。"
                         "履歴の単位を切り分ける実験用")
    ap.add_argument("--pad-dir", type=int, default=1,
                    help="詰め物の分岐方向（1=常に成立, 0=常に不成立）")
    ap.add_argument("--trials", type=int, default=7)
    ap.add_argument("--reps", type=int, default=48)
    ap.add_argument("--warmup", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-o", "--out")
    args = ap.parse_args()

    if not os.path.exists(args.binary):
        raise SystemExit(f"{args.binary} が見つかりません。make してください。")
    if args.patlen is None:
        print("注意: --patlen が未指定です。既定値を使いますが、"
              "env/calibrate_patlen.sh で校正した値を渡すことを強く推奨します。",
              file=sys.stderr)

    p1s, p2s = parse_range(args.p1), parse_range(args.p2)

    if args.null:
        nf = noise_floor(args.binary, args.mode, p1s[0], p2s[0], args)
        print(f"ノイズ床（帰無実験, ns/分岐）: {nf:.4f}")
        print("これより小さい差は有意と見なせない。")
        return

    # noise_floor を列に持つ理由: これが無いと significant を CSV だけから
    # 再計算できない（ノイズ床が端末出力にしか残らない）。判定式を変えたときに
    # 過去の結果を再評価できなくなるため、判定の入力はすべて CSV に置く。
    fields = ["p1", "p1_name", "p2", "p2_name",
              "test_ns", "test_iqr", "ctl_ns", "ctl_iqr",
              "diff_ns", "noise_floor", "random_sat", "significant", "anomaly",
              "test_mpb", "ctl_mpb", "diff_mpb",
              "test_ipc", "ctl_ipc", "test_fe_per_br", "ctl_fe_per_br"]
    # cross の p1 は d であり、標的が必要とする相関距離は 2d+1 **パターン要素**
    # である（d ではない）。A と B が交互に実行されるため。換算を解析者任せに
    # すると解釈を誤るので派生列として出す。他のモードでは p1 の意味が違い
    # 2*p1+1 に意味が無いため、列そのものを出さない。
    #
    # **これは相関距離であって履歴長ではない。** ハードウェアの履歴が数える単位は
    # パターン要素ではなく成立分岐であり（仕様 §3.2）、さらに cross の転移は
    # 容量側の制約を含む（仕様 §3.4）。旧名 hist_len_lower は「履歴長の下限・
    # 分岐数」を含意していて二重に誤っていたので corr_dist_elements に改めた。
    cross_mode = args.mode == "cross"
    if cross_mode:
        # p1_name の直後に置く（p1 と p1_name を隣接させたまま派生列を足す）
        fields.insert(2, "corr_dist_elements")
    w = fout = None
    if args.out:
        fout = open(args.out, "w", newline="")
        # 行終端を LF に固定する。csv モジュールの既定は CRLF で、位置や文字列で
        # 読む側が \r を落とし忘れると黙って不一致になる（実際に踏んだ）。
        w = csv.DictWriter(fout, fieldnames=fields, lineterminator="\n")
        w.writeheader()

    # 測定条件としての負荷。**停止はさせない**（判定材料にしない）が、記録して
    # 後から結果を捨てられるようにする。
    load_start = load_avg()
    if load_start:
        print(f"load average（開始時）: {fmt_load(load_start)}")
        if load_start[0] > LOAD_WARN:
            print(f"警告: load average が {load_start[0]:.2f} で "
                  f"{LOAD_WARN} を超えています。反復のばらつきが大きくなり、"
                  "境界判定が「境界なし」に落ちることがあります"
                  "（実測: load 20 で ctx K=2 の 66/67 が消えた）。\n"
                  "  他のアプリを止めてから測り直すことを推奨します。"
                  "この警告で測定は止めません。", file=sys.stderr)

    # 軸名。CSV と .meta.txt を自己記述にする（p1/p2 が何なのかを単独で読めるように）。
    p1_name, p2_name = query_axes(args.binary, args.mode)
    print(f"掃引軸: p1={p1_name}  p2={p2_name}（モード {args.mode}）")

    # 純ランダム飽和値の基準線。意図せず純ランダムを測っていないかの検査に使う。
    sat = float("nan")
    if not args.no_random_check:
        sat = random_saturation(args.binary, args)
        if sat != sat:
            print("注意: 純ランダム飽和値を測れませんでした"
                  "（意図せぬ純ランダムの検査は無効）", file=sys.stderr)
        else:
            print(f"純ランダム飽和値（{args.binary} の histd --control）: "
                  f"{sat:.4f} ns/要素"
                  f"  これに近い test_ns は警告します")

    nf = float("nan")
    if not args.no_control:
        nf = noise_floor(args.binary, args.mode, p1s[0], p2s[0], args)
        print(f"ノイズ床（帰無実験）: {nf:.4f} ns/分岐"
              "  （これより小さい差は有意でない）\n")

    hdr = f"{'p1':>8}"
    if cross_mode:
        hdr += f" {'2d+1要素':>8}"
    hdr += f" {'p2':>8} {'test_ns':>9} {'iqr':>7}"
    if args.no_control:
        hdr += f" {'ipc':>6} {'fe/br':>9}"
    if not args.no_control:
        hdr += f" {'ctl_ns':>9} {'iqr':>7} {'diff_ns':>9} {'有意':>5} {'異常':>4}"
    print(hdr)

    for p1 in p1s:
        for p2 in p2s:
            if args.no_control:
                rows = invoke(args.binary, args.mode, p1, p2, args)
                tm, ti = stats([float(r[NS]) for r in rows])
                mm, _ = stats([float(r[MPB]) for r in rows])
                ipc, _ = stats([float(r[IPC]) for r in rows])
                fe = fe_per_branch(rows)
                rflag = random_flag(tm, ti, sat)
                if rflag:
                    print(f"警告: p1={p1:g} p2={p2:g} で {rflag}。"
                          "意図せず純ランダムを測っている可能性があります。",
                          file=sys.stderr)
                lead = f"{p1:>8g}"
                if cross_mode:
                    lead += f" {2 * p1 + 1:>8g}"
                print(f"{lead} {p2:>8g} {tm:>9.4f} {ti:>7.4f} "
                      f"{ipc:>6.3f} {fe:>9.5f}")
                row = dict(p1=p1, p1_name=p1_name, p2=p2, p2_name=p2_name,
                           test_ns=tm, test_iqr=ti, ctl_ns="",
                           ctl_iqr="", diff_ns="", noise_floor="",
                           random_sat=sat,
                           significant="", anomaly="random_like" if rflag else "",
                           test_mpb=mm, ctl_mpb="", diff_mpb="",
                           test_ipc=ipc, ctl_ipc="", test_fe_per_br=fe,
                           ctl_fe_per_br="")
            else:
                rows = invoke(args.binary, args.mode, p1, p2, args, pair=True)
                tv, cv = split_pair(rows, NS)
                tmv, cmv = split_pair(rows, MPB)
                tm, ti = stats(tv)
                cm, ci = stats(cv)
                diff = cm - tm
                sig, anom = judge(diff, nf, ti, ci)
                rflag = random_flag(tm, ti, sat)
                if rflag:
                    print(f"警告: p1={p1:g} p2={p2:g} で {rflag}。"
                          "意図せず純ランダムを測っている可能性があります"
                          "（生成器の引数が範囲外・パターン長不足・記憶化など）。"
                          "有意でないこと自体は異常ではありません。",
                          file=sys.stderr)
                    anom = (anom + "; " if anom else "") + "random_like"
                if anom.startswith("negative_diff("):
                    # 設計上ありえない向きの差で、ばらつきでも説明できない。
                    # 黙って no にせず前に出す。
                    print(f"警告: p1={p1:g} p2={p2:g} で test が control より"
                          f"ノイズ床とばらつきを超えて遅い（{anom}）。"
                          "生成器か測定系の異常を疑うこと。", file=sys.stderr)
                lead = f"{p1:>8g}"
                if cross_mode:
                    lead += f" {2 * p1 + 1:>8g}"
                print(f"{lead} {p2:>8g} {tm:>9.4f} {ti:>7.4f} "
                      f"{cm:>9.4f} {ci:>7.4f} {diff:>+9.4f} {sig:>5} "
                      f"{'!' if anom.startswith('negative_diff') else '':>4}")
                ti_ipc, _ = stats([float(r[IPC]) for r in rows if r["control"] == "0"])
                ci_ipc, _ = stats([float(r[IPC]) for r in rows if r["control"] == "1"])
                row = dict(p1=p1, p1_name=p1_name, p2=p2, p2_name=p2_name,
                           test_ns=tm, test_iqr=ti, ctl_ns=cm,
                           ctl_iqr=ci, diff_ns=diff, noise_floor=nf,
                           random_sat=sat, significant=sig, anomaly=anom,
                           test_mpb=stats(tmv)[0], ctl_mpb=stats(cmv)[0],
                           diff_mpb=stats(cmv)[0] - stats(tmv)[0],
                           test_ipc=ti_ipc, ctl_ipc=ci_ipc,
                           test_fe_per_br=fe_per_branch(
                               [r for r in rows if r["control"] == "0"]),
                           ctl_fe_per_br=fe_per_branch(
                               [r for r in rows if r["control"] == "1"]))
            if cross_mode:
                # 相関距離（パターン要素）。履歴長ではない。cross の p1=d に対し
                # 2d+1 要素。要素と分岐は 1 対 1 でなく、不成立分岐は履歴を
                # 消費しない（仕様 §3.2）。
                row["corr_dist_elements"] = 2 * p1 + 1
            if w:
                w.writerow(row)

    if fout:
        fout.close()
        print(f"\n{args.out} に書き出しました")
        write_meta(args.out, args, args.binary, nf, sat, p1_name, p2_name,
                   load_start)


if __name__ == "__main__":
    main()
