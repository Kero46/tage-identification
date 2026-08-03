#!/bin/sh
# 測定環境の設定。要 root。測定の質はここで決まる。
# 使い方: sudo env/setup.sh [固定するコア番号(既定 2)]
set -e
CORE="${1:-2}"

echo "== 周波数を性能重視で固定 =="
if command -v cpupower >/dev/null 2>&1; then
    cpupower frequency-set -g performance || echo "  cpupower 失敗(続行)"
else
    for f in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor; do
        [ -w "$f" ] && echo performance > "$f" || true
    done
fi

echo "== turbo/boost を無効化 =="
if [ -w /sys/devices/system/cpu/intel_pstate/no_turbo ]; then
    echo 1 > /sys/devices/system/cpu/intel_pstate/no_turbo
    echo "  intel_pstate: no_turbo=1"
elif [ -w /sys/devices/system/cpu/cpufreq/boost ]; then
    echo 0 > /sys/devices/system/cpu/cpufreq/boost
    echo "  cpufreq: boost=0"
else
    echo "  turbo 制御不可(結果のばらつきに注意)"
fi

echo "== パフォーマンスカウンタの利用許可 =="
sysctl -w kernel.perf_event_paranoid=1

echo "== コア $CORE の SMT 兄弟 =="
SIB="/sys/devices/system/cpu/cpu$CORE/topology/thread_siblings_list"
[ -r "$SIB" ] && echo "  $(cat $SIB)  (兄弟コアに負荷を置かないこと)" || true

echo
echo "完了。測定は taskset で固定して実行すること:"
echo "  taskset -c $CORE tools/sweep.py history_length/hist_bench --mode histd --p1 3:255:2"
