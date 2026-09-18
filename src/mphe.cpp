#include "include/mphe.h"

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <vector>

#include "src/64header.h"

namespace {

#ifndef CONTROL_BIT_TORUS_SHIFT64
constexpr int kControlBitTorusShift64 = 59;
#else
constexpr int kControlBitTorusShift64 = CONTROL_BIT_TORUS_SHIFT64;
#endif

Torus64 CenteredControlBitTargetLocal(int bit) {
    const int64_t mu = static_cast<int64_t>(UINT64_C(1) << kControlBitTorusShift64);
    return (bit & 1) ? static_cast<Torus64>(mu) : static_cast<Torus64>(-mu);
}

void AddTorusPolynomial64(Torus64Polynomial* out,
                          const Torus64Polynomial* in,
                          const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        out->coefs[i] += in->coefs[i];
    }
}

void BinarySecretKeyGen64(IntPolynomiala* sk, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        sk->coefs[i] = random_bit();
    }
}

void SampleCommonAVector64(Torus64* a, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        a[i] = random_int64();
    }
}

void MPGenerateSingleBinarySecretShares64(IntPolynomiala* secret_out,
                                          IntPolynomiala* const* share_out, int num_users,
                                          const Globals* env) {
    assert(num_users > 0);
    BinarySecretKeyGen64(secret_out, env);

    for (int u = 0; u < num_users; ++u) {
        ClearIntPolynomial64(share_out[u], env);
    }

    const int half_range = 1 << 23;
    for (int i = 0; i < env->N; ++i) {
        int running_sum = 0;
        for (int u = 0; u + 1 < num_users; ++u) {
            const int share = (std::rand() % (2 * half_range + 1)) - half_range;
            share_out[u]->coefs[i] = share;
            running_sum += share;
        }
        share_out[num_users - 1]->coefs[i] = secret_out->coefs[i] - running_sum;
    }
}

void MPKeyGenControlLweShare64(LweSample64* pk_share,
                               const Torus64* common_a, const IntPolynomiala* sk_share,
                               double stdev, const Globals* env) {
    Torus64 inner = random_gaussian64(0, stdev);

    for (int i = 0; i < env->N; ++i) {
        pk_share->a[i] = common_a[i];
        inner += common_a[i] * sk_share->coefs[i];
    }

    *pk_share->b = inner;
}

void MPAggregateControlLwePublicKey64(LweSample64* common_pk,
                                      LweSample64* const* pk_shares, int num_users,
                                      const Globals* env) {
    assert(num_users > 0);

    for (int i = 0; i < env->N; ++i) {
        common_pk->a[i] = pk_shares[0]->a[i];
    }

    Torus64 b_sum = 0;
    for (int u = 0; u < num_users; ++u) {
        b_sum += *pk_shares[u]->b;
    }
    *common_pk->b = b_sum;
}

}  

Torus64 EncodeControlBitToTorus64(int bit) {
    return CenteredControlBitTargetLocal(bit & 1);
}

void ClearIntPolynomial64(IntPolynomiala* poly, const Globals* env) {
    for (int i = 0; i < env->N; ++i) poly->coefs[i] = 0;
}

void ClearLweSample64(LweSample64* ct, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        ct->a[i] = 0;
    }
    *ct->b = 0;
}

void CopyLweSample64(LweSample64* out, const LweSample64* in, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        out->a[i] = in->a[i];
    }
    *out->b = *in->b;
}

void AddLweSample64(LweSample64* out, const LweSample64* in, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        out->a[i] += in->a[i];
    }
    *out->b += *in->b;
}

void SubLweSample64(LweSample64* out, const LweSample64* in, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        out->a[i] -= in->a[i];
    }
    *out->b -= *in->b;
}

void ClearTorusPolynomial64(Torus64Polynomial* poly, const Globals* env) {
    for (int i = 0; i < env->N; ++i) poly->coefs[i] = 0;
}

void ClearTLweSample64(TLweSample64* ct, const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        for (int i = 0; i < env->N; ++i) {
            ct->a[q].coefs[i] = 0;
        }
    }
}

void CopyTLweSample64(TLweSample64* out, const TLweSample64* in, const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        for (int i = 0; i < env->N; ++i) {
            out->a[q].coefs[i] = in->a[q].coefs[i];
        }
    }
}

void AddTLweSample64(TLweSample64* out, const TLweSample64* in, const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        for (int i = 0; i < env->N; ++i) {
            out->a[q].coefs[i] += in->a[q].coefs[i];
        }
    }
}

void SampleCommonA64(Torus64Polynomial* a, const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        a->coefs[i] = random_int64();
    }
}

void MPAggregateSecretKey64(IntPolynomiala* sk_sum,
                            IntPolynomiala* const* sk_list,int num_users,
                            const Globals* env) {
    ClearIntPolynomial64(sk_sum, env);
    for (int u = 0; u < num_users; ++u) {
        for (int i = 0; i < env->N; ++i) {
            sk_sum->coefs[i] += sk_list[u]->coefs[i];
        }
    }
}

void MPKeyGenShare64(TLweSample64* pk_share,
                     const Torus64Polynomial* common_a,const IntPolynomiala* sk_share,
                     double stdev,const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        pk_share->a[0].coefs[i] = common_a->coefs[i];
        pk_share->b->coefs[i] = random_gaussian64(0, stdev);
    }

    torus64PolynomialMultAddKaratsuba_lvl2(pk_share->b, sk_share, common_a, env);
}

void MPAggregatePublicKey64(TLweSample64* common_pk,
                            TLweSample64* const* pk_shares,int num_users,
                            const Globals* env) {
    assert(num_users > 0);

    for (int i = 0; i < env->N; ++i) {
        common_pk->a[0].coefs[i] = pk_shares[0]->a[0].coefs[i];
        common_pk->b->coefs[i] = 0;
    }

    for (int u = 0; u < num_users; ++u) {
        for (int i = 0; i < env->N; ++i) {
            common_pk->b->coefs[i] += pk_shares[u]->b->coefs[i];
        }
    }
}

void MPPartialDecrypt64(Torus64Polynomial* share,
                        const TLweSample64* ct,const IntPolynomiala* sk_share,
                        double stdev,const Globals* env) {
    torus64PolynomialMultKaratsuba_lvl2(share, sk_share, &ct->a[0], env);

    for (int i = 0; i < env->N; ++i) {
        share->coefs[i] += random_gaussian64(0, stdev);
    }
}

void MPCombineShares64(Torus64Polynomial* sum,
                       Torus64Polynomial* const* shares,int num_users,
                       const Globals* env) {
    ClearTorusPolynomial64(sum, env);
    for (int u = 0; u < num_users; ++u) {
        AddTorusPolynomial64(sum, shares[u], env);
    }
}

void MPFinalPhase64(Torus64Polynomial* phase,
                    const TLweSample64* ct,const Torus64Polynomial* share_sum,
                    const Globals* env) {
    for (int i = 0; i < env->N; ++i) {
        phase->coefs[i] = ct->b->coefs[i] - share_sum->coefs[i];
    }
}


// Each party receives only its own additive secret shares
void MPSetupPartiesLweControl64(std::vector<MPPartyKey64>& parties,
                                LweSample64* control_pk,
                                TLweSample64* data_pk,
                                Torus64* control_common_a,
                                Torus64Polynomial* data_common_a,
                                int num_users,
                                double stdev,
                                const Globals* env) {
    assert(num_users > 0);

    parties.clear();
    parties.reserve(num_users);

    // Sample the public 'a' components shared by all parties
    SampleCommonAVector64(control_common_a, env);
    SampleCommonA64(data_common_a, env);

    IntPolynomiala* control_secret = new IntPolynomiala(env->N);
    IntPolynomiala* data_secret = new IntPolynomiala(env->N);

    std::vector<IntPolynomiala*> control_shares(num_users, nullptr);
    std::vector<IntPolynomiala*> data_shares(num_users, nullptr);
    std::vector<LweSample64*> control_pk_ptrs(num_users, nullptr);
    std::vector<TLweSample64*> data_pk_ptrs(num_users, nullptr);

    for (int u = 0; u < num_users; ++u) {
        control_shares[u] = new IntPolynomiala(env->N);
        data_shares[u] = new IntPolynomiala(env->N);
    }

    // Locally emulate the MPC output: binary secrets are additively shared
    MPGenerateSingleBinarySecretShares64(control_secret, control_shares.data(), num_users, env);
    MPGenerateSingleBinarySecretShares64(data_secret, data_shares.data(), num_users, env);

    for (int u = 0; u < num_users; ++u) {
        LweSample64* control_pk_share = new LweSample64(env->N);
        TLweSample64* data_pk_share = new TLweSample64(env->N);

        // Generate party u's pk contribution from its local secret share
        // These pk shares are needed to construct the global pk from distributed secret shares without reconstructing the global sk.
        MPKeyGenControlLweShare64(control_pk_share, control_common_a, control_shares[u], stdev, env);
        MPKeyGenShare64(data_pk_share, data_common_a, data_shares[u], stdev, env);

        // Store only party u's local secret shares and pk contributions
        parties.emplace_back(
            control_shares[u],
            data_shares[u],
            data_pk_share,
            control_pk_share
        );
        control_pk_ptrs[u] = control_pk_share;
        data_pk_ptrs[u] = data_pk_share;
    }

     // Aggregate all pk contributions into the common pk
    MPAggregateControlLwePublicKey64(control_pk, control_pk_ptrs.data(), num_users, env);
    MPAggregatePublicKey64(data_pk, data_pk_ptrs.data(), num_users, env);

    delete control_secret;
    delete data_secret;
}

void MPDestroyParties64(std::vector<MPPartyKey64>& parties) {
    for (std::size_t i = 0; i < parties.size(); ++i) {
        delete parties[i].control_sk_share;
        delete parties[i].data_sk_share;
        delete parties[i].data_pk_share;
        delete parties[i].control_lwe_pk_share;
        parties[i].control_sk_share = nullptr;
        parties[i].data_sk_share = nullptr;
        parties[i].data_pk_share = nullptr;
        parties[i].control_lwe_pk_share = nullptr;
    }
    parties.clear();
}
