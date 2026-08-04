# プロジェクト全体のビルド・確認
DIRS = history_length table_structure

.PHONY: all common $(DIRS) test verify verify-selfcheck sanity check calibrate clean
all: $(DIRS)

common:
	$(MAKE) -C common

$(DIRS): common
	$(MAKE) -C $@

# パターン生成器の性質検査（必要履歴長の下限を含む）。測定前に必ず通すこと。
test:
	$(MAKE) -C common test

# 機械語と出力形式の検査。環境ごとに必須。
#   verify_selfcheck  検査そのものが陽性・陰性の両方向で動くこと（メタテスト）
#   verify_branch     データ依存の条件分岐が全カーネルに期待数だけあること
#   verify_csv        ベンチマークの CSV 列構成が仕様書 §2.3 から動いていないこと
#   verify_sweep_csv  掃引ドライバの出力列（cross のみ派生列を持つ）
verify: all verify-selfcheck
	@env/verify_branch.sh history_length/hist_bench
	@env/verify_branch.sh table_structure/table_bench
	@env/verify_csv.sh history_length/hist_bench --mode random
	@env/verify_csv.sh table_structure/table_bench --mode histd --param 3
	@env/verify_sweep_csv.sh

# 機械語検査の回帰検出。単体でも回せるようにしておく。
verify-selfcheck:
	@env/verify_selfcheck.sh

# 測定系が分岐予測を測れていることの合否ゲート。
# 陽性対照と陰性対照の時間比が閾値を下回ったら非零終了する。
# 校正済みのパターン長があれば PATLEN=<n> make sanity で渡す。
sanity: all
	@env/sanity.sh

# 測定環境の確認(root 不要)
check:
	@env/check.sh

# パターン長の校正（記憶化を避ける長さを決める）
calibrate: all
	@env/calibrate_patlen.sh

clean:
	$(MAKE) -C common clean
	@for d in $(DIRS); do $(MAKE) -C $$d clean; done
