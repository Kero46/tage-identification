# プロジェクト全体のビルド・確認
DIRS = history_length table_structure

.PHONY: all common $(DIRS) test verify sanity check calibrate clean
all: $(DIRS)

common:
	$(MAKE) -C common

$(DIRS): common
	$(MAKE) -C $@

# パターン生成器の性質検査（必要履歴長の下限を含む）。測定前に必ず通すこと。
test:
	$(MAKE) -C common test

# 条件分岐が cmov に化けていないかの確認。環境ごとに必須。
verify: all
	@env/verify_branch.sh history_length/hist_bench
	@env/verify_branch.sh table_structure/table_bench

# 測定系が信号を検出できることの確認（陽性対照・陰性対照・ノイズ床）
sanity: all
	@echo "== 陽性対照: 予測可能な周期パターンは速いはず =="
	@./history_length/hist_bench --mode period --param 10 --trials 3 --reps 16 2>/dev/null
	@echo "== 陰性対照: 純ランダムは遅いはず =="
	@./history_length/hist_bench --mode random --trials 3 --reps 16 2>/dev/null
	@echo "== ノイズ床（帰無実験）=="
	@tools/sweep.py history_length/hist_bench --mode histd --p1 15 \
	    --trials 7 --reps 16 --null 2>/dev/null

# 測定環境の確認(root 不要)
check:
	@env/check.sh

# パターン長の校正（記憶化を避ける長さを決める）
calibrate: all
	@env/calibrate_patlen.sh

clean:
	$(MAKE) -C common clean
	@for d in $(DIRS); do $(MAKE) -C $$d clean; done
