// measure.c — パフォーマンスカウンタと時間の計測
//
// カウンタが使えない環境(権限/サンドボックス)でも時間計測は常に有効で、
// これだけで掃引は成立する。課題も時間計測による推定を想定している。
#define _GNU_SOURCE
#include "bp.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <asm/unistd.h>

static int fd_br = -1, fd_miss = -1, fd_cyc = -1, fd_ins = -1, fd_fe = -1, ok = 0;

static long perf_open(uint64_t config, int group) {
    struct perf_event_attr a;
    memset(&a, 0, sizeof(a));
    a.type = PERF_TYPE_HARDWARE;
    a.size = sizeof(a);
    a.config = config;
    a.disabled = (group == -1);
    a.exclude_kernel = 1;
    a.exclude_hv = 1;
    // イベントが多重化されると値が黙って縮尺されるため、リーダを pinned にし、
    // 有効時間と実行時間を読んで縮尺の有無を検出する。
    a.pinned = (group == -1);
    a.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
    return syscall(__NR_perf_event_open, &a, 0, -1, group, 0);
}

// 観測されたミスの要因を切り分けるには、分岐ミスとサイクル数だけでは足りない。
// 命令数と分岐数も取り、必要ならフロントエンド停止サイクルも併せて取る。
void bp_ctr_init(int want_frontend) {
    fd_br = (int)perf_open(PERF_COUNT_HW_BRANCH_INSTRUCTIONS, -1);
    if (fd_br < 0) return;
    fd_miss = (int)perf_open(PERF_COUNT_HW_BRANCH_MISSES, fd_br);
    fd_cyc  = (int)perf_open(PERF_COUNT_HW_CPU_CYCLES, fd_br);
    fd_ins  = (int)perf_open(PERF_COUNT_HW_INSTRUCTIONS, fd_br);
    ok = (fd_miss >= 0 && fd_cyc >= 0 && fd_ins >= 0);
    if (ok && want_frontend) {
        // 物理カウンタ数を超えると多重化されるため任意扱い。
        // 縮尺は scaled フラグで検出できる。
        fd_fe = (int)perf_open(PERF_COUNT_HW_STALLED_CYCLES_FRONTEND, fd_br);
    }
}
int bp_ctr_ok(void) { return ok; }

static double t0_ns;
static int scaled = 0;

// 1イベント読み出し: value, time_enabled, time_running
static uint64_t read_one(int fd, int *is_scaled) {
    uint64_t buf[3] = {0, 0, 0};
    if (read(fd, buf, sizeof(buf)) != (ssize_t)sizeof(buf)) return 0;
    if (buf[2] == 0 || buf[2] < buf[1]) *is_scaled = 1;
    return buf[0];
}

void bp_ctr_start(void) {
    if (ok) {
        ioctl(fd_br, PERF_EVENT_IOC_RESET, PERF_IOC_FLAG_GROUP);
        ioctl(fd_br, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    }
    t0_ns = bp_now_ns();
}
void bp_ctr_stop(bp_sample *o) {
    double t1 = bp_now_ns();
    o->branches = o->misses = o->cycles = 0;
    o->instructions = o->fe_stalls = 0;
    o->ns = t1 - t0_ns;
    if (!ok) return;
    ioctl(fd_br, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
    scaled = 0;
    o->branches     = read_one(fd_br,   &scaled);
    o->misses       = read_one(fd_miss, &scaled);
    o->cycles       = read_one(fd_cyc,  &scaled);
    o->instructions = read_one(fd_ins,  &scaled);
    o->fe_stalls    = (fd_fe >= 0) ? read_one(fd_fe, &scaled) : 0;
}
int bp_ctr_scaled(void) { return scaled; }
#else
void bp_ctr_init(int want_frontend) { (void)want_frontend; }
int  bp_ctr_ok(void) { return 0; }
int  bp_ctr_scaled(void) { return 0; }
static double t0_ns;
void bp_ctr_start(void) { t0_ns = bp_now_ns(); }
void bp_ctr_stop(bp_sample *o) {
    o->branches = o->misses = o->cycles = 0;
    o->instructions = o->fe_stalls = 0;
    o->ns = bp_now_ns() - t0_ns;
}
#endif

double bp_now_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e9 + (double)ts.tv_nsec;
}
