#pragma once

#include <cmath>
#include <cstdint>

#include "poc_64types.h"

Globals::Globals(int N_in,
                 int bgbit_in,
                 int l_in,
                 int basebit_in,
                 int t_in,
                 int smalln_in,
                 int k_in)
    : k(k_in),
      N(N_in),
      t(t_in),
      smalln(smalln_in),
      bgbit(bgbit_in),
      l(l_in),
      basebit(basebit_in) {
    torusDecompBuf = new uint64_t[N];
    in_key = new int64_t[N];
    lwekey = new int[N];
    tlwekey = new IntPolynomiala(N);
    privKS = nullptr;

    uint64_t offset = 0;
    const uint64_t Bg = UINT64_C(1) << bgbit;
    const int64_t halfBg = Bg / 2;
    for (int p = 0; p < l; ++p) {
        offset += (UINT64_C(1) << (64 - (p + 1) * bgbit)) * halfBg;
    }
    torusDecompOffset = offset;
}

static inline int RoundScaledCoeff(Torus64 x, int64_t scale) {
    const long double v = static_cast<long double>(static_cast<int64_t>(x));
    const long double q = v / static_cast<long double>(scale);
    return static_cast<int>(std::llround(q));
}
