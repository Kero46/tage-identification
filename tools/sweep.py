#!/usr/bin/env python3
"""掃引ドライバ。

test（相関あり）と control（相関なし）を **同一プロセス内で交互に** 測り、
その差分から標的分岐の予測可能性を取り出す。差分によって、履歴を運ぶ
フィラー分岐が出すミスの寄与が相殺される。

交互測定にする理由: test と control を別プロセスで逐次実行すると、
周波数ドリフトが差分に系統誤差として乗る。

差が本物かを判断するため、ばらつき（四分位範囲）と、帰無実験
（control 同士の比較。差は零であるべき）によるノイズ床を併せて出す。

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


def write_meta(path, args, binary):
    """再現に必要な環境情報を残す。"""
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
        fh.write(f"p1={args.p1} p2={args.p2} patlen={args.patlen or '既定'}\n")
        fh.write(f"trials={args.trials} reps={args.reps} "
                 f"warmup={args.warmup} seed={args.seed}\n")
        fh.write(f"sites={args.sites or '単一'} "
                 f"pad={'なし' if args.pad is None else args.pad}"
                 f"(dir={args.pad_dir})\n")
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

    fields = ["p1", "p2", "test_ns", "test_iqr", "ctl_ns", "ctl_iqr",
              "diff_ns", "significant", "test_mpb", "ctl_mpb", "diff_mpb",
              "test_ipc", "ctl_ipc", "test_fe_per_br", "ctl_fe_per_br"]
    # cross の p1 は d であり、必要なグローバル履歴長は 2d+1 分岐である（d ではない）。
    # A と B が交互に実行されるため。換算を解析者任せにすると解釈を誤るので、
    # 派生列として出す。他のモードでは p1 の意味が違い 2*p1+1 に意味が無いため、
    # 列そのものを出さない（意味を持たない数値を並べない）。
    cross_mode = args.mode == "cross"
    if cross_mode:
        fields.insert(1, "hist_len_lower")
    w = fout = None
    if args.out:
        fout = open(args.out, "w", newline="")
        # 行終端を LF に固定する。csv モジュールの既定は CRLF で、位置や文字列で
        # 読む側が \r を落とし忘れると黙って不一致になる（実際に踏んだ）。
        w = csv.DictWriter(fout, fieldnames=fields, lineterminator="\n")
        w.writeheader()

    nf = float("nan")
    if not args.no_control:
        nf = noise_floor(args.binary, args.mode, p1s[0], p2s[0], args)
        print(f"ノイズ床（帰無実験）: {nf:.4f} ns/分岐"
              "  （これより小さい差は有意でない）\n")

    hdr = f"{'p1':>8}"
    if cross_mode:
        hdr += f" {'2d+1':>6}"
    hdr += f" {'p2':>8} {'test_ns':>9} {'iqr':>7}"
    if args.no_control:
        hdr += f" {'ipc':>6} {'fe/br':>9}"
    if not args.no_control:
        hdr += f" {'ctl_ns':>9} {'iqr':>7} {'diff_ns':>9} {'有意':>5}"
    print(hdr)

    for p1 in p1s:
        for p2 in p2s:
            if args.no_control:
                rows = invoke(args.binary, args.mode, p1, p2, args)
                tm, ti = stats([float(r[NS]) for r in rows])
                mm, _ = stats([float(r[MPB]) for r in rows])
                ipc, _ = stats([float(r[IPC]) for r in rows])
                fe = fe_per_branch(rows)
                lead = f"{p1:>8g}"
                if cross_mode:
                    lead += f" {2 * p1 + 1:>6g}"
                print(f"{lead} {p2:>8g} {tm:>9.4f} {ti:>7.4f} "
                      f"{ipc:>6.3f} {fe:>9.5f}")
                row = dict(p1=p1, p2=p2, test_ns=tm, test_iqr=ti, ctl_ns="",
                           ctl_iqr="", diff_ns="", significant="",
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
                sig = "yes" if (diff == diff and nf == nf and abs(diff) > nf) else "no"
                lead = f"{p1:>8g}"
                if cross_mode:
                    lead += f" {2 * p1 + 1:>6g}"
                print(f"{lead} {p2:>8g} {tm:>9.4f} {ti:>7.4f} "
                      f"{cm:>9.4f} {ci:>7.4f} {diff:>+9.4f} {sig:>5}")
                ti_ipc, _ = stats([float(r[IPC]) for r in rows if r["control"] == "0"])
                ci_ipc, _ = stats([float(r[IPC]) for r in rows if r["control"] == "1"])
                row = dict(p1=p1, p2=p2, test_ns=tm, test_iqr=ti, ctl_ns=cm,
                           ctl_iqr=ci, diff_ns=diff, significant=sig,
                           test_mpb=stats(tmv)[0], ctl_mpb=stats(cmv)[0],
                           diff_mpb=stats(cmv)[0] - stats(tmv)[0],
                           test_ipc=ti_ipc, ctl_ipc=ci_ipc,
                           test_fe_per_br=fe_per_branch(
                               [r for r in rows if r["control"] == "0"]),
                           ctl_fe_per_br=fe_per_branch(
                               [r for r in rows if r["control"] == "1"]))
            if cross_mode:
                # 必要なグローバル履歴長の下限（分岐数）。cross の p1=d に対し 2d+1。
                row["hist_len_lower"] = 2 * p1 + 1
            if w:
                w.writerow(row)

    if fout:
        fout.close()
        print(f"\n{args.out} に書き出しました")
        write_meta(args.out, args, args.binary)


if __name__ == "__main__":
    main()
