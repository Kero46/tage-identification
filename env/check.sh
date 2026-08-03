#!/bin/sh
# 測定前の環境確認。root 不要。
CORE="${1:-2}"
echo "== CPU =="
grep -m1 'model name' /proc/cpuinfo 2>/dev/null || true
echo "  コア数: $(nproc)"

echo "== 周波数ガバナ =="
G="/sys/devices/system/cpu/cpu$CORE/cpufreq/scaling_governor"
[ -r "$G" ] && echo "  cpu$CORE: $(cat $G)" || echo "  取得不可"

echo "== turbo =="
if [ -r /sys/devices/system/cpu/intel_pstate/no_turbo ]; then
    echo "  no_turbo=$(cat /sys/devices/system/cpu/intel_pstate/no_turbo) (1 が望ましい)"
elif [ -r /sys/devices/system/cpu/cpufreq/boost ]; then
    echo "  boost=$(cat /sys/devices/system/cpu/cpufreq/boost) (0 が望ましい)"
else
    echo "  制御なし"
fi

echo "== perf_event_paranoid =="
P=$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo "?")
echo "  $P (2 以下が必要。カウンタ不可でも時間計測は動作する)"

echo "== SMT =="
[ -r /sys/devices/system/cpu/smt/active ] && echo "  smt active=$(cat /sys/devices/system/cpu/smt/active)" || true

echo "== 負荷 =="
uptime
