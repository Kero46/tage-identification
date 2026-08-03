CC     ?= gcc
CFLAGS ?= -O2 -g -Wall -Wextra -fno-tree-vectorize -I../common
LIB     = ../common/libbp.a

hist_bench: hist_bench.c $(LIB)
	$(CC) $(CFLAGS) -o $@ hist_bench.c $(LIB)

$(LIB):
	$(MAKE) -C ../common

.PHONY: clean verify
clean:
	rm -f hist_bench
verify: hist_bench
	@../env/verify_branch.sh hist_bench
