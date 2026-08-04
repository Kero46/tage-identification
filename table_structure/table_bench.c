// TAGE テーブル構成の同定
//
// 観測は三系統。
//   dual : 二重相関(距離 d1, d2 の XOR)。両者が同一の履歴段で処理される場合に
//          競合が生じることを利用して段の区切りを定める。
//   ctx  : 文脈数 K を制限した相関系列。K を増やしてスラッシングが始まる点から
//          当該テーブルの実効エントリ数を推定する。
//   histd: 履歴長測定と同じ高密度相関。段の検出(距離の細かい掃引)に用いる。
//
// 注意: dual の競合はタグ衝突・インデックス衝突でも生じうる。ctx 掃引と
//       組み合わせて切り分けること(README の対照実験を参照)。
#include "bp.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void g_dual (uint8_t *p, size_t n, double a, double b, int c) { bp_gen_dual(p, n, (unsigned)a, (unsigned)b, c); }
static void g_ctx  (uint8_t *p, size_t n, double a, double b, int c) { bp_gen_hist_ctx(p, n, (unsigned)a, (unsigned)b, c); }
static void g_histd(uint8_t *p, size_t n, double a, double b, int c) { (void)b; bp_gen_hist_dense(p, n, (unsigned)a, c); }
// aliasing 実験。S 個の分岐サイトに周期 p の固有パターンを与える。
// S を増やすと衝突が起こる「可能性が高まる」実験であり、観測された劣化を
// aliasing と断定するものではない（コード量も増えるため）。
static void g_alias(uint8_t *p, size_t n, double a, double b, int c) { bp_gen_alias(p, n, (unsigned)a, (unsigned)b, c); }

static void usage(const char *p) {
    fprintf(stderr,
      "使い方: %s --mode <dual|ctx|alias|histd> [オプション]\n"
      "  dual : --param D1  --param2 D2\n"
      "  ctx  : --param D   --param2 K   (K = 文脈数)\n"
      "  alias: --param S   --param2 p   (S = 分岐サイト数 2,4,8,16,32,64 / p = 周期)\n"
      "  histd: --param D   (奇数)\n"
      "  --pad <P>    1 要素につき常に同じ方向に分岐する詰め物分岐を P 本挿入する\n"
      "               (0,1,2,4)。予測ミスを増やさず履歴スロットだけを消費するので、\n"
      "               履歴の単位を切り分けられる。--sites との併用は不可\n"
      "  --pad-dir <d> 詰め物の分岐方向 (1=常に成立[既定], 0=常に不成立)\n"
      "  --sites <S>  ctx/dual/histd を S 個の分岐サイトに分散して実行する\n"
      "               (2,4,8,16,32,64)。既定は単一サイト。\n"
      "               目的: 履歴の単位を確かめる。単一サイトのカーネルは\n"
      "               パターン 1 要素につきサイト分岐 1 本とループ制御分岐 1 本を\n"
      "               実行するが、S サイト版は S 要素につきループ分岐 1 本になる。\n"
      "               履歴が条件分岐すべてで更新されるなら、要素で測った相関距離の\n"
      "               限界は S を増やすと伸びる。伸びなければ履歴はサイト分岐だけを\n"
      "               数えている。\n"
      "  --control / --pair / --patlen / --frontend / --reps / --warmup / --trials / --seed / --csv\n", p);
}

int main(int argc, char **argv) {
    bp_config c;
    memset(&c, 0, sizeof(c));
    c.trials = 7; c.reps = 48; c.warmup = 8; c.seed = 1;
    const char *mode = NULL;
    int want_fe = 0;
    unsigned sites = 0;                 // 0 なら単一サイト
    int pads = -1, pad_dir = 1;         // pads<0 なら詰め物なし

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "--mode")   && i+1 < argc) mode = argv[++i];
        else if (!strcmp(argv[i], "--param")  && i+1 < argc) c.p1 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--param2") && i+1 < argc) c.p2 = atof(argv[++i]);
        else if (!strcmp(argv[i], "--control"))              c.control = 1;
        else if (!strcmp(argv[i], "--pair"))                 c.pair = 1;
        else if (!strcmp(argv[i], "--patlen")  && i+1 < argc) bp_pat_len = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--reps")   && i+1 < argc) c.reps   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--warmup") && i+1 < argc) c.warmup = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--trials") && i+1 < argc) c.trials = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--seed")   && i+1 < argc) c.seed   = strtoull(argv[++i],0,10);
        else if (!strcmp(argv[i], "--sites")  && i+1 < argc) sites = (unsigned)atoi(argv[++i]);
        else if (!strcmp(argv[i], "--pad")    && i+1 < argc) pads = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--pad-dir") && i+1 < argc) pad_dir = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--frontend"))             want_fe = 1;
        else if (!strcmp(argv[i], "--csv"))                  c.csv = 1;
        else { usage(argv[0]); return 2; }
    }
    if (!mode) { usage(argv[0]); return 2; }

    static unsigned nsites;
    if      (!strcmp(mode,"alias")) {
        nsites = (unsigned)c.p1;
        if (nsites != 2 && nsites != 4 && nsites != 8 &&
            nsites != 16 && nsites != 32 && nsites != 64) {
            fprintf(stderr, "エラー: alias の --param は 2,4,8,16,32,64 のいずれか\n");
            return 2;
        }
        c.gen = g_alias; c.kernel = bp_kernel_sites; c.ctx = &nsites;
        // alias を差分法で解析してはならないので、対照系列を要求する指定を拒む。
        // 理由は二つある。
        //   1) 対照系列は「学習不可能なランダム」なので、差は「予測可能か否か」に
        //      支配され、aliasing による緩やかな劣化が埋もれる（仕様 §4.2）。
        //   2) test 系列の成立率は周期語 w（S*p ビットの独立乱数）の平均で決まり、
        //      標本が S*p 個しかないため control（全長にわたる乱数）と一致しない。
        //      実測で S=8,p=8 のとき 0.05〜0.11、S=2,p=8 では 0.19 ずれた。
        //      成立率がずれた二系列の差には予測性能以外の要因が混ざる
        //      （成立率で実行時間が 1.6 倍動くため。lab-notebook の診断 6）。
        // 文書だけの約束にしておくと守られないので、ここで機械的に止める。
        if (c.pair || c.control) {
            fprintf(stderr,
                "エラー: alias は差分法で解析できない（--pair / --control は使えない）\n"
                "  対照系列との差は「予測可能か否か」に支配され、aliasing による\n"
                "  緩やかな劣化が埋もれる。加えて test と control で成立率が構造的に\n"
                "  ずれるため、差に予測性能以外の要因が混ざる。\n"
                "  --no-control で test_ns の S 依存を見ること（仕様 §4.2）:\n"
                "    tools/sweep.py table_structure/table_bench --mode alias \\\n"
                "        --p1 2:64:x2 --p2 8 --no-control --frontend\n");
            return 2;
        }
    }
    else if (!strcmp(mode,"dual"))  c.gen = g_dual;
    else if (!strcmp(mode,"ctx"))   c.gen = g_ctx;
    else if (!strcmp(mode,"histd")) c.gen = g_histd;
    else { usage(argv[0]); return 2; }

    if (!strcmp(mode,"histd") && ((unsigned)c.p1 % 2) == 0) {
        fprintf(stderr, "エラー: histd は D が奇数でなければならない\n");
        return 2;
    }
    // --sites: パターン要素を S 個の分岐サイト（異なる PC）に分散して実行する。
    // 生成する系列は変わらない。変わるのは要素あたりのループ制御分岐の本数で、
    // 単一サイトの 1 本から 1/S 本に減る。履歴の単位（条件分岐すべてか、
    // サイト分岐だけか）を判別するための実験用。
    static unsigned nsites_opt;
    if (sites) {
        if (!strcmp(mode, "alias")) {
            fprintf(stderr, "エラー: alias のサイト数は --param で指定する"
                            "（--sites は使えない）\n");
            return 2;
        }
        if (sites != 2 && sites != 4 && sites != 8 &&
            sites != 16 && sites != 32 && sites != 64) {
            fprintf(stderr, "エラー: --sites は 2,4,8,16,32,64 のいずれか\n");
            return 2;
        }
        nsites_opt = sites;
        c.kernel = bp_kernel_sites;
        c.ctx = &nsites_opt;
    }

    // --pad: 予測可能な詰め物分岐を 1 要素あたり P 本挿入する。
    // 分岐サイト（PC）の数を大きく変えずに、要素あたりの条件分岐の本数だけを
    // 変えられるので、--sites より交絡が少ない。
    static bp_pad_cfg padcfg;
    if (pads >= 0) {
        if (sites) {
            fprintf(stderr, "エラー: --pad と --sites は併用できない"
                            "（要素あたりの分岐本数の勘定が混ざる）\n");
            return 2;
        }
        if (pads != 0 && pads != 1 && pads != 2 && pads != 4) {
            fprintf(stderr, "エラー: --pad は 0,1,2,4 のいずれか\n");
            return 2;
        }
        if (pad_dir != 0 && pad_dir != 1) {
            fprintf(stderr, "エラー: --pad-dir は 0 か 1\n");
            return 2;
        }
        padcfg.pads = (unsigned)pads;
        padcfg.dir  = (unsigned)pad_dir;
        c.kernel = bp_kernel_pad;
        c.ctx = &padcfg;
    }

    c.mode = mode;
    if (!c.kernel) c.kernel = bp_kernel_single;
    bp_ctr_init(want_fe);
    if (!bp_ctr_ok())
        fprintf(stderr, "警告: カウンタ取得不可。時間計測(ns/br)を使うこと。\n");
    if (c.csv) bp_print_header();
    return bp_run(&c);
}
