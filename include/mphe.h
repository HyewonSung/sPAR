#pragma once

#include <vector>

#include "poc_64types.h"

struct Globals;
struct IntPolynomiala;
struct LweSample64;
struct Torus64Polynomial;
struct TLweSample64;

struct MPPartyKey64 {
    IntPolynomiala* control_sk_share;
    IntPolynomiala* data_sk_share;
    TLweSample64* data_pk_share;
    LweSample64* control_lwe_pk_share;

    MPPartyKey64(IntPolynomiala* control_in = nullptr,
                 IntPolynomiala* data_in = nullptr,
                 TLweSample64* data_pk_in = nullptr,
                 LweSample64* control_lwe_pk_in = nullptr)
        : control_sk_share(control_in),
          data_sk_share(data_in),
          data_pk_share(data_pk_in),
          control_lwe_pk_share(control_lwe_pk_in) {}
};

Torus64 EncodeControlBitToTorus64(int bit);

void ClearIntPolynomial64(IntPolynomiala* poly, const Globals* env);
void ClearLweSample64(LweSample64* ct, const Globals* env);
void CopyLweSample64(LweSample64* out, const LweSample64* in, const Globals* env);
void AddLweSample64(LweSample64* out, const LweSample64* in, const Globals* env);
void SubLweSample64(LweSample64* out, const LweSample64* in, const Globals* env);

void ClearTorusPolynomial64(Torus64Polynomial* poly, const Globals* env);
void ClearTLweSample64(TLweSample64* ct, const Globals* env);
void CopyTLweSample64(TLweSample64* out, const TLweSample64* in, const Globals* env);
void AddTLweSample64(TLweSample64* out, const TLweSample64* in, const Globals* env);

void SampleCommonA64(Torus64Polynomial* a, const Globals* env);

void MPAggregateSecretKey64(IntPolynomiala* sk_sum,
                            IntPolynomiala* const* sk_list,
                            int num_users,
                            const Globals* env);

void MPKeyGenShare64(TLweSample64* pk_share,
                     const Torus64Polynomial* common_a,
                     const IntPolynomiala* sk_share,
                     double stdev,
                     const Globals* env);

void MPAggregatePublicKey64(TLweSample64* common_pk,
                            TLweSample64* const* pk_shares,
                            int num_users,
                            const Globals* env);


void MPPartialDecrypt64(Torus64Polynomial* share,
                        const TLweSample64* ct,
                        const IntPolynomiala* sk_share,
                        double stdev,
                        const Globals* env);

void MPCombineShares64(Torus64Polynomial* sum,
                       Torus64Polynomial* const* shares,
                       int num_users,
                       const Globals* env);

void MPFinalPhase64(Torus64Polynomial* phase,
                    const TLweSample64* ct,
                    const Torus64Polynomial* share_sum,
                    const Globals* env);

void MPSetupPartiesLweControl64(std::vector<MPPartyKey64>& parties,
                                LweSample64* control_pk,
                                TLweSample64* data_pk,
                                Torus64* control_common_a,
                                Torus64Polynomial* data_common_a,
                                int num_users,
                                double stdev,
                                const Globals* env);

void MPDestroyParties64(std::vector<MPPartyKey64>& parties);
