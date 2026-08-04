# tage-identification

分岐予測器の内部構成を、外部から観測できる情報だけを手がかりに
リバースエンジニアリングして同定する試み。

大学院講義「先進計算機構成論」第11回の課題として実施している。

## 目的

現代 CPU の条件分岐予測器は、構成が公開されていない。本リポジトリでは、
規則性を制御した分岐列を実行し、予測が当たったか外れたかという外部からの
観測だけを手がかりに、その内部構成を同定する。

得られた結果は、構成が既知のシミュレータ（gem5）を真値オラクルとして
手法の妥当性を確認したうえで実機に適用する。実機で観測できるのは予測ミスの
総数のみであり、個々のミスを内部機構に帰属させることはできないため、
実機に対する同定は予測成否のパターンからの推論である。

## 測定対象

| 段 | 対象 | 位置づけ |
|---|---|---|
| 第1段 | TAGE の履歴長 | 本体（必達） |
| 第2段 | TAGE のテーブルのサイズや数 | 本体（必達） |
| 第3段 | 第2段で説明できなかった残差の分析 | 発展（結果次第） |

第1段と第2段で研究は完結する。第3段はその自然な発展であり、内容は第2段の
結果を見てから決める。

実験は次の三本柱で構成する。

- **周期パターン** — 固定周期の系列をどこまで学習できるか
- **相関パターン** — 分岐間の相関をどの距離まで学習できるか
- **aliasing 実験** — 分岐サイト数を増やして衝突が生じる可能性を高める

## ディレクトリ構成

```
.
├── common/               測定基盤（ライブラリ libbp.a）
│   ├── bp.h              共通インタフェース
│   ├── rng.c             再現性のある乱数
│   ├── patterns.c        分岐パターン生成
│   ├── measure.c         パフォーマンスカウンタと時間の計測
│   ├── kernel.c          計測カーネル（単一／2 サイト／多サイト）
│   ├── runner.c          試行ループ
│   └── selftest.c        生成器の性質検査
├── history_length/       第1段：TAGE の履歴長の同定
│   ├── hist_bench.c
│   └── results/
├── table_structure/      第2段：テーブルのサイズや数の同定
│   ├── table_bench.c
│   ├── kernel_index_bits.c  索引ビット特定用（未実装・設計メモ）
│   └── results/
├── gem5/                 真値オラクルとしての検証（未整備）
│   ├── configs/
│   └── results/
├── env/                  測定環境の構築と不変条件の検査
│   ├── setup.sh          周波数固定・turbo 無効・カウンタ許可（要 root）
│   ├── check.sh          測定前の環境確認
│   ├── calibrate_patlen.sh  パターン長の校正
│   ├── sanity.sh         分岐予測を測れていることの合否ゲート
│   ├── verify_branch.sh  データ依存分岐の生存とサイト数の検査
│   ├── branch_check.awk  逆アセンブル解析（レジスタを追って分岐を数える）
│   ├── verify_selfcheck.sh  機械語検査そのものの回帰検出
│   ├── verify_csv.sh     ベンチマークの CSV 列構成のドリフト検出
│   ├── verify_sweep_csv.sh  掃引ドライバの出力列のドリフト検出
│   └── csvcol.awk        CSV 列を名前で解決する（位置読み禁止）
├── tools/                解析
│   ├── sweep.py          掃引ドライバ（交互測定・ノイズ床・メタデータ記録）
│   └── plot.py           作図
└── docs/
    ├── research_plan.md       研究計画書
    ├── implementation_spec.md 実装仕様書
    ├── measurement_guide.md   測定手順書（設計判断と落とし穴）
    ├── stage1_handoff.md      第1段から第2段への引き継ぎ
    └── lab-notebook.md        実験ノート
```

測定対象ごとのディレクトリと、共通基盤・環境スクリプト・解析ツールが
同じ階層に並ぶ。各測定ディレクトリは自分の CLI と結果を持ち、共通基盤を
ライブラリとして使う。

## 必要なもの

**測定は Linux / x86-64 + GCC で行う。**

- Linux — `perf_event_open` を使用。PMU が要る受け入れ基準（alias の `fe/br`・
  `ipc` 併読）はこの環境でしか満たせない
- x86-64 — 分岐命令を明示的に書いているため（`asm goto`）。aarch64 でも動く
- GCC — ブロック複製を防ぐ最適化属性のため。**分岐の生存自体は属性ではなく
  `asm goto` が保証**しており、複製の有無は `make verify` が検査する
- Python 3（解析スクリプト。`matplotlib` は任意）

macOS / arm64 / clang は副次対象。ビルド・実行・全検査は通り、**時間計測のみで
パターン長の校正と第1段の掃引は成立する**（実証済み）。一方で PMU が無いため
alias の受け入れ基準と環境メタデータの要件は満たせない。満たせない基準の一覧は
`docs/implementation_spec.md` の §0.2 にある。

## ビルドと実行

```sh
make            # 全ビルド
make test       # 生成器の性質検査
make verify     # 機械語（分岐の生存・サイト数）と CSV 列構成の検査
make sanity     # 分岐予測を測れていることの合否ゲート（通らなければ非零終了）
make check      # 測定環境の確認
make calibrate  # パターン長の校正
```

`make test && make verify && make sanity` が通らない状態の測定値は使わない。
検査が落ちたら、検査を緩めるのではなく実装を直す。

測定の最小例：

```sh
sudo env/setup.sh 2                      # コア 2 を測定用に設定
taskset -c 2 make calibrate              # PATLEN を決める
taskset -c 2 tools/sweep.py history_length/hist_bench \
    --mode cross --p1 1:64:1 --patlen 65536 \
    -o history_length/results/cross.csv
```

`ar` がランダムアクセス書き込みを必要とするため、共有マウント上ではビルドに
失敗することがある。通常のローカルファイルシステムでビルドすること。

詳しい手順、設計判断、既知の落とし穴は `docs/measurement_guide.md` を参照。

## 実装状況

| 項目 | 状況 |
|---|---|
| 測定基盤（生成器・計測・カーネル・試行ループ） | 実装済 |
| 自己検査、環境スクリプト、掃引ドライバ | 実装済 |
| 第1段の全モード | 実装済 |
| **第1段の測定（機械 B、時間計測のみ）** | **実施済**（`docs/stage1_handoff.md`） |
| 第2段の主要観測（距離掃引・二重相関・文脈数・aliasing） | 実装済・未測定 |
| 第2段：索引ビットの特定 | 未実装 |
| 第2段：学習過程の観測 | 未実装 |
| gem5 による検証（手法の妥当性確認） | **未整備（第1段の保留事項）** |
| Linux / x86-64 + GCC での再現 | 未実施（第1段の保留事項） |

## 前提と限界

- 対象実機の予測器が TAGE 系である保証はない。**それ自体も検証の対象**であり、
  リポジトリ名は想定を表すものにすぎない。
- 実機の予測器構成は非公開であるため、断定ではなく「観測はこの構成と整合する」
  という形で述べる。履歴長やテーブルサイズは下限ないし制約として報告する。
- aliasing 実験は衝突が生じる可能性を高めるものであり、観測された劣化を
  aliasing と断定するものではない。

## ドキュメント

| 文書 | 内容 |
|---|---|
| [`docs/research_plan.md`](docs/research_plan.md) | 研究計画。三段構成と各段の論理 |
| [`docs/implementation_spec.md`](docs/implementation_spec.md) | 実装仕様。API、モード仕様、受け入れ基準 |
| [`docs/measurement_guide.md`](docs/measurement_guide.md) | 測定手順、設計判断、落とし穴 |
| [`docs/stage1_handoff.md`](docs/stage1_handoff.md) | 第1段の結論と第2段への申し送り |
| [`docs/lab-notebook.md`](docs/lab-notebook.md) | 実験ノート。判明した欠陥と修正の記録 |
