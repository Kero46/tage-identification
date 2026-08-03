CC     ?= gcc
CFLAGS ?= -O2 -g -Wall -Wextra -fno-tree-vectorize -I.
OBJS    = rng.o patterns.o measure.o kernel.o runner.o

libbp.a: $(OBJS)
	ar rcs $@ $(OBJS)

# パターン生成器の性質検査（必要履歴長の下限を含む）
bp_selftest: selftest.c libbp.a
	$(CC) $(CFLAGS) -o $@ selftest.c libbp.a

%.o: %.c bp.h
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONY: clean test
test: bp_selftest
	./bp_selftest

clean:
	rm -f *.o libbp.a bp_selftest
