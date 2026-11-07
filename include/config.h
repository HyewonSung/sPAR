#ifndef SPAR_CONFIG_H
#define SPAR_CONFIG_H

// Every value can be overridden at compile time with -D (see build.sh, PRESET=paper).
// Gadgets are given as (log2 of the base, number of levels).

#ifndef REFRESH_HW
#define REFRESH_HW   72     // refresh hasWritten every REFRESH_HW candidates (0 = never)
#endif
#ifndef REFRESH_I
#define REFRESH_I    32     // refresh I every REFRESH_I messages (0 = never)
#endif
#ifndef LWE_N
#define LWE_N        1024   // reduced LWE dimension n
#endif
#ifndef SIGMA_LOG
#define SIGMA_LOG    (-54)  // encryption, evaluation-key and partial-decryption noise, sigma = 2^SIGMA_LOG
#endif
#ifndef SIGMA_KSK_LOG
#define SIGMA_KSK_LOG (-26) // N -> n key-switching key noise
#endif
#ifndef KS_OUT_BASEBIT
#define KS_OUT_BASEBIT 8
#endif
#ifndef KS_OUT_LEVELS
#define KS_OUT_LEVELS 2     // output-KSK keeps the top KS_OUT_LEVELS*KS_OUT_BASEBIT bits
#endif
#ifndef GATE_BG
#define GATE_BG      24     // gate-PBS blind-rotation gadget
#endif
#ifndef GATE_L
#define GATE_L       1
#endif
#ifndef OUT_BG
#define OUT_BG       9      // CB output RGSW gadget (also HomTrace / scheme-switch key base)
#endif
#ifndef OUT_L
#define OUT_L        3
#endif
#ifndef CB_PBS_BG
#define CB_PBS_BG    10     // CB blind-rotation gadget
#endif
#ifndef CB_PBS_L
#define CB_PBS_L     4
#endif
#ifndef CB_AKS_L
#define CB_AKS_L     6      // HomTrace key levels (base 2^OUT_BG)
#endif
#ifndef CB_L_SS
#define CB_L_SS      5      // scheme-switching key levels (base 2^OUT_BG)
#endif
#ifndef PAYLOAD_WIDTH
#define PAYLOAD_WIDTH 8     // payload coefficients per message, 4 bits each; coefficient 0 is the sort key
#endif
#ifndef PARTY_COUNT
#define PARTY_COUNT  4
#endif
#ifndef ETA
#define ETA          4
#endif
#ifndef SEED
#define SEED         0
#endif

// ---- compile-time checks: an invalid combination fails the build instead of the run ----
static_assert(ETA == PARTY_COUNT, "ETA and PARTY_COUNT must be equal (one message per client)");
static_assert(ETA >= 3, "ETA must be at least 3 (each client chooses 3 distinct buckets)");
static_assert(LWE_N > 0 && LWE_N <= 2048, "LWE_N must be in [1, N=2048]");
static_assert(REFRESH_HW >= 0 && REFRESH_I >= 0, "refresh periods must be >= 0");
static_assert(SIGMA_LOG < 0 && SIGMA_KSK_LOG < 0, "noise levels are given as negative log2 values");
static_assert(KS_OUT_BASEBIT > 0 && KS_OUT_BASEBIT < 32 && KS_OUT_LEVELS > 0 &&
              KS_OUT_BASEBIT * KS_OUT_LEVELS <= 64, "N -> n key switching needs basebit < 32 and basebit*levels <= 64");
static_assert(GATE_BG > 0 && GATE_BG < 32 && GATE_L > 0 && GATE_BG * GATE_L <= 64,
              "gate gadget needs 0 < GATE_BG < 32 and GATE_BG*GATE_L <= 64");
static_assert(CB_PBS_BG > 0 && CB_PBS_BG < 32 && CB_PBS_L > 0 && CB_PBS_BG * CB_PBS_L <= 64,
              "CB blind-rotation gadget needs 0 < CB_PBS_BG < 32 and CB_PBS_BG*CB_PBS_L <= 64");
static_assert(OUT_BG > 0 && OUT_BG < 32 && OUT_L > 0 && OUT_BG * OUT_L <= 64,
              "CB output gadget needs 0 < OUT_BG < 32 and OUT_BG*OUT_L <= 64");
static_assert(CB_AKS_L > 0 && OUT_BG * CB_AKS_L <= 64, "HomTrace keys need OUT_BG*CB_AKS_L <= 64");
static_assert(CB_L_SS > 0 && OUT_BG * CB_L_SS <= 64, "scheme-switching key needs OUT_BG*CB_L_SS <= 64");
static_assert(PAYLOAD_WIDTH >= 4 && PAYLOAD_WIDTH <= 2048, "PAYLOAD_WIDTH must be in [4, N=2048]");

#endif
