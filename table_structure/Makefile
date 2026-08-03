CC     ?= gcc
CFLAGS ?= -O2 -g -Wall -Wextra -fno-tree-vectorize -I../common
LIB     = ../common/libbp.a

table_bench: table_bench.c $(LIB)
	$(CC) $(CFLAGS) -o $@ table_bench.c $(LIB)

$(LIB):
	$(MAKE) -C ../common

.PHONY: clean verify
clean:
	rm -f table_bench
verify: table_bench
	@../env/verify_branch.sh table_bench
