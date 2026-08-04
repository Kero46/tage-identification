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
// 雑音注入 probe。--param D --param2 m。文脈数は構成上 2 固定なので K は取らない。
static void g_ctxnoise(uint8_t *p, size_t n, double a, double b, int c) { bp_gen_ctx_noise(p, n, (unsigned)a, (unsigned)b, c); }

// 掃引軸の名前。**--param / --param2 の意味はモードごとに違う**ので、CSV や
// .meta.txt を単独で読んだときに何の軸なのかが分かるように名前を持たせる。
// 位置と慣習で意味が決まっている状態は落とし穴 15（CSV を位置で読む）と同種。
// ここが単一の出所であり、sweep.py は --axes で問い合わせる（二重管理にしない）。
static const char *axis_p1(const char *mode) {
    if (!strcmp(mode, "dual"))     return "D1";
    if (!strcmp(mode, "ctx"))      return "D";
    if (!strcmp(mode, "ctxnoise")) return "D";
    if (!strcmp(mode, "alias"))    return "S";
    if (!strcmp(mode, "histd"))    return "D";
    return "p1";
}
static const char *axis_p2(const char *mode) {
    if (!strcmp(mode, "dual"))     return "D2";
    if (!strcmp(mode, "ctx"))      return "K";
    if (!strcmp(mode, "ctxnoise")) return "m";
    if (!strcmp(mode, "alias"))    return "period";
    if (!strcmp(mode, "histd"))    return "unused";
    return "p2";
}

static void usage(const char *p) {
    fprintf(stderr,
      "使い方: %s --mode <dual|ctx|ctxnoise|alias|histd> [オプション]\n"
      "\n"
      "  【モードごとの引数】--param / --param2 の意味はモードで違う。誤用しないこと。\n"
      "  モード     --param      --param2        備考\n"
      "  --------   ----------   -------------   --------------------------------\n"
      "  dual       D1 (距離1)   D2 (距離2)      標的 = bit(D1前) XOR bit(D2前)\n"
      "  ctx        D (距離)     K (文脈数)      K <= 2^(D-1) を強制。2 のべきに丸め\n"
      "  ctxnoise   D (距離)     m (雑音の本数)  **文脈数は 2 固定。K は取らない**\n"
      "  alias      S (サイト数) period (周期)   S は 2,4,8,16,32,64 のみ\n"
      "  histd      D (距離、奇数)  未使用       \n"
      "\n"
      "  --axes       そのモードの軸名を表示して終了（p1_name / p2_name）\n"
      "\n"
      "  ctxnoise: 梯子の段の位置 L_i を求める（第2段 (a) の主手法）。\n"
      "     ブロック = [雑音 m][識別ビット 1][共通サフィックス D-1][標的]\n"
      "     距離 D より古い側の雑音が「D より長いテーブルだけ」を選択的に潰す。\n"
      "     成功(D) <=> ∃i: D <= L_i <= D + W_i。成功区間の上端が段の履歴長。\n"
      "     **m は W(~10) を超えないと何も潰れない**（m=0 に退化する）。m~20 を使う。\n"
      "     **m=0 は ctx K=2 と同じ構造**で、第1段の 66/67 境界が再現するはず。\n"
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
    int want_fe = 0, want_axes = 0;
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
        else if (!strcmp(argv[i], "--axes"))                 want_axes = 1;
        else if (!strcmp(argv[i], "--frontend"))             want_fe = 1;
        else if (!strcmp(argv[i], "--csv"))                  c.csv = 1;
        else { usage(argv[0]); return 2; }
    }
    if (!mode) { usage(argv[0]); return 2; }

    // 掃引軸の名前を問い合わせる経路。sweep.py がこれを読んで CSV と .meta.txt に
    // p1_name / p2_name を書く。**軸の意味の出所をここ一箇所に保つ**ため、
    // 解析側に同じ表を持たせない（二重管理は必ずずれる）。
    if (want_axes) {
        printf("p1_name=%s\n", axis_p1(mode));
        printf("p2_name=%s\n", axis_p2(mode));
        return 0;
    }

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
    else if (!strcmp(mode,"ctxnoise")) c.gen = g_ctxnoise;
    else if (!strcmp(mode,"histd")) c.gen = g_histd;
    else { usage(argv[0]); return 2; }

    if (!strcmp(mode,"histd") && ((unsigned)c.p1 % 2) == 0) {
        fprintf(stderr, "エラー: histd は D が奇数でなければならない\n");
        return 2;
    }
    // ctxnoise の引数検証。--param2 は m（雑音の本数）であって K ではない。
    //
    // **文脈数は構成上 2 に固定である**（識別ビット 1 個）。K を渡す経路は存在しない
    // ので「K を拒否する」検査は書けないが、代わりに m の範囲を検証する。
    // 誤って K のつもりで大きな値を渡すと m が過大になり、ブロック長が
    // パターン長を超えて 1 ブロックも生成されなくなる（黙って純ランダムに近い
    // 系列を測ることになる）。これは落とし穴 24 と同じ失敗様式なので拒否する。
    if (!strcmp(mode,"ctxnoise")) {
        unsigned D = (unsigned)c.p1, m = (unsigned)c.p2;
        if (D < 1) {
            fprintf(stderr, "エラー: ctxnoise は D >= 1 が必要（指定 D=%u）\n", D);
            return 2;
        }
        size_t blk = (size_t)m + D + 1;
        size_t blocks = blk ? bp_pat_len / blk : 0;
        // **構造的に不可能な場合はエラー**にする。ブロックが数個も取れないなら、
        // 系列はほぼ末尾の埋め草（純ランダム）になり、黙って別のものを測ることに
        // なる（落とし穴 24 と同じ失敗様式）。K のつもりで大きな値を渡すとここに来る。
        if (blocks < 8) {
            fprintf(stderr,
                "エラー: ctxnoise の D=%u m=%u ではブロック長 %zu に対し patlen %zu が"
                "短すぎる（ブロック %zu 個）\n"
                "  **--param2 は m（雑音の本数）であって K（文脈数）ではない。**\n"
                "  文脈数は構成上 2 固定である。K のつもりで大きな値を渡していないか"
                "確認すること。\n",
                D, m, blk, bp_pat_len, blocks);
            return 2;
        }
        // **測定に足りない場合は警告**にとどめる。閾値で拒否すると小さな patlen の
        // 動作確認ができなくなる（検査スクリプトが実際に落ちた）。
        // 潰したいテーブルの文脈数 2^(W+2) を踏むには文脈あたり 10 ブロック必要で、
        // W~10 なら 20480 ブロック（仕様 §4.1 の制約 4）。
        if (blocks < 1024)
            fprintf(stderr,
                "警告: ctxnoise の D=%u m=%u でブロックが %zu 個しかない"
                "（測定には 20480 個程度必要）。\n"
                "  patlen を上げること（目安 2097152、仕様 §4.1 の制約 4）。\n",
                D, m, blocks);
    }
    // ctx の文脈数の上限。識別ビット j = log2(K) は標的から見て距離 D..D-j+1 に
    // 置くので j <= D-1 が必要、すなわち K <= 2^(D-1)。
    //
    // **以前はこれを検証しておらず、生成器が黙って純ランダム系列を返していた。**
    // K 掃引すると上限を超えた点で diff が 0 に落ちるので「容量の限界」に見える。
    // 実測で D=5, K=32 のとき test_ns が純ランダムの飽和値 3.38 になり、
    // K=16→32 が劣化開始点に見えた（実験ノート 2026-08-04 続報 4）。
    // 小さい D ほど早く当たるため、梯子の短履歴側の段がまるごと偽物になる。
    if (!strcmp(mode,"ctx")) {
        unsigned D = (unsigned)c.p1, K = (unsigned)c.p2;
        if (K < 2) K = 2;
        unsigned kmax = bp_ctx_max_k(D);
        if (kmax == 0 || K > kmax) {
            fprintf(stderr,
                "エラー: ctx は D=%u に対し K <= 2^(D-1) = %u でなければならない"
                "（指定 K=%u）\n"
                "  識別ビット j=log2(K) は標的から見て距離 D..D-j+1 に置くので\n"
                "  j <= D-1 が必要である。範囲外の K では文脈を表現できない。\n"
                "  以前はここを検証せず、生成器が黙って純ランダム系列を返していた。\n"
                "  K 掃引に偽の劣化が現れ、容量の限界と誤読する（仕様 §4.1）。\n"
                "  掃引の上限は min(8192, 2^(D-1)) にすること。\n",
                D, kmax, K);
            return 2;
        }
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
