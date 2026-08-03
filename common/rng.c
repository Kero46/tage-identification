#include "bp.h"
static uint64_t st = 1;
void bp_seed(uint64_t s) { st = s ? s : 1; }
uint64_t bp_rnd_u64(void) {
    uint64_t x = st;
    x ^= x >> 12; x ^= x << 25; x ^= x >> 27;
    st = x;
    return x * 0x2545F4914F6CDD1DULL;
}
int bp_rnd_bit(void) { return (int)(bp_rnd_u64() >> 33) & 1; }
int bp_rnd_biased(double p) {
    return ((double)(bp_rnd_u64() >> 11) / 9007199254740992.0) < p ? 1 : 0;
}
