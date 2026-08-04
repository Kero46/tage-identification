// TAGE の履歴長の同定（骨子 1）
//
// 目的: 分岐予測器が用いる履歴長をリバースエンジニアリングして同定する。
//       併せて、以降の測定が乗る基盤を確立する。
//
// 【報告上の注意】実機の予測器は非公開かつ複雑であるため、観測から直接
// 確定できるのは「この距離・この周期までの規則性は学習できた」という事実である。
// 履歴長として述べる際は、それが下限であることを明記する。
// 目的が同定であることと、個々の主張を証拠の強さに合わせることは両立する。
//
// モードの使い分け:
//   cross  分岐間相関。グローバル履歴の効果を評価できる（推奨）
//          パラメータ d に対し、標的が必要とする相関距離は 2d+1 パターン要素
//          である。**これは履歴長ではない。** 要素と分岐は 1 対 1 でなく、
//          不成立分岐は履歴を消費しない（仕様 §3.2）。さらに cross の転移は
//          容量側の制約を含むため履歴長として報告しない（仕様 §3.4）
//   histd  単一サイトの自己相関。ローカル履歴だけでも学習できるため、
//          グローバル履歴の評価には使えない
//   period 周期パターン。ローカル履歴でも学習可能な点に注意
#include "bp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void g_hist   (uint8_t *p, size_t n, double a, double b, int c) { (void)b; bp_gen_hist(p, n, (unsigned)a, c); }
static void g_histd  (uint8_t *p, size_t n, double a, double b, int c) { (void)b; bp_gen_hist_dense(p, n, (unsigned)a, c); }
static void g_period (uint8_t *p, size_t n, double a, double b, int c) { (void)b; (void)c; bp_gen_period(p, n, (unsigned)a); }
static void g_bias   (uint8_t *p, size_t n, double a, double b, int c) { (void)b; (void)c; bp_gen_bias(p, n, a); }
static void g_random (uint8_t *p, size_t n, double a, double b, int c) { (void)a; (void)b; (void)c; bp_gen_random(p, n); }
// 分岐間相関（グローバル履歴の評価）。単一サイトの自己相関はローカル履歴でも
// 学習できるため、グローバル履歴の効果を見るにはこちらを使う。
static void g_cross  (uint8_t *p, size_t n, double a, double b, int c) { (void)b; bp_gen_cross(p, n, (unsigned)a, c); }

// 掃引軸の名前。**--param の意味はモードごとに違う**ので、CSV や .meta.txt を
// 単独で読んだときに何の軸なのかが分かるように名前を持たせる（table_bench と同じ方針）。
// 軸の意味の出所をここ一箇所に保ち、解析側に同じ表を持たせない。
static const char *axis_p1(const char *mode) {
    if (!strcmp(mode, "cross"))  return "d";
    if (!strcmp(mode, "histd"))  return "D";
    if (!strcmp(mode, "hist"))   return "D";
    if (!strcmp(mode, "period")) return "N";
    if (!strcmp(mode, "bias"))   return "p";
    if (!strcmp(mode, "random")) return "unused";
    return "p1";
}
static const char *axis_p2(const char *mode) { (void)mode; return "unused"; }

static void usage(const char *p) {
    fprintf(stderr,
      "使い方: %s --mode <cross|histd|hist|period|bias|random> [オプション]\n"
      "  cross : 分岐間相関。グローバル履歴の評価に用いる（推奨）\n"
      "          --param d に対し相関距離は 2d+1 パターン要素（履歴長ではない）\n"
      "  histd : 単一サイトの自己相関。ローカル履歴でも学習可能な点に注意\n"
      "  --param <v>    D (histd/hist, histd は奇数) / d (cross) / N (period) / p (bias)\n"
      "  --param2 <v>   第2引数。第1段のどのモードも使わない（掃引ドライバが\n"
      "                 一律に渡すため受け付ける。第2段の dual/ctx/alias で使う）\n"
      "  --control      対照系列(相関なし)を生成する\n"
      "  --reps <n>     計測時の反復回数 (既定 48)\n"
      "  --warmup <n>   計測外の学習反復 (既定 8)\n"
      "  --trials <n>   試行回数 (既定 7)\n"
      "  --seed <n>     乱数種 (既定 1)\n"
      "  --pair         test/control を同一プロセス内で交互に測る(推奨)\n"
      "  --patlen <n>   パターン長(分岐数)。短いと予測器が系列を記憶する\n"
      "  --frontend     フロントエンド停止サイクルも取る(多重化に注意)\n"
      "  --axes         そのモードの軸名を表示して終了（p1_name / p2_name）\n"
      "  --csv          CSV 出力\n", p);
}

int main(int argc, char **argv) {
    bp_config c;
    memset(&c, 0, sizeof(c));
    c.trials = 7; c.reps = 48; c.warmup = 8; c.seed = 1;
    const char *mode = NULL;
    int want_fe = 0, want_axes = 0;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--mode")    && i+1 < argc) mode = argv[++i];
        else if (!strcmp(argv[i], "--param")   && i+1 < argc) c.p1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--param2")  && i+1 < argc) c.p2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--control"))               c.control = 1;
        else if (!strcmp(argv[i], "--pair"))                  c.pair = 1;
        else if (!strcmp(argv[i], "--patlen")  && i+1 < argc)  bp_pat_len = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--reps")    && i+1 < argc) c.reps   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--warmup")  && i+1 < argc) c.warmup = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--trials")  && i+1 < argc) c.trials = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--seed")    && i+1 < argc) c.seed   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--axes"))                  want_axes = 1;
        else if (!strcmp(argv[i], "--frontend"))              want_fe = 1;
        else if (!strcmp(argv[i], "--csv"))                   c.csv = 1;
        else { usage(argv[0]); return 2; }
    }
    if (!mode) { usage(argv[0]); return 2; }

    // 掃引軸の名前を問い合わせる経路（table_bench と同じ）。
    if (want_axes) {
        printf("p1_name=%s\n", axis_p1(mode));
        printf("p2_name=%s\n", axis_p2(mode));
        return 0;
    }

    if      (!strcmp(mode,"cross"))  { c.gen = g_cross; c.kernel = bp_kernel_cross; }
    else if (!strcmp(mode,"histd"))  c.gen = g_histd;
    else if (!strcmp(mode,"hist"))   c.gen = g_hist;
    else if (!strcmp(mode,"period")) c.gen = g_period;
    else if (!strcmp(mode,"bias"))   c.gen = g_bias;
    else if (!strcmp(mode,"random")) c.gen = g_random;
    else { usage(argv[0]); return 2; }

    if (!strcmp(mode,"histd") && ((unsigned)c.p1 % 2) == 0) {
        fprintf(stderr, "エラー: histd は D が奇数でなければならない (指定 %g)\n"
                        "  理由: 参照先を必ずフィラー位置(偶数)に保つため\n", c.p1);
        return 2;
    }

    c.mode = mode;
    if (!c.kernel) c.kernel = bp_kernel_single;
    bp_ctr_init(want_fe);
    if (!bp_ctr_ok())
        fprintf(stderr, "警告: カウンタ取得不可。misses は 0。時間計測(ns/br)を使うこと。\n");
    if (c.csv) bp_print_header();
    return bp_run(&c);
}
