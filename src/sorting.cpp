#include "include/sorting.h"

#include "src/64header.h"
#include "include/mphe.h"


void HomCompareSwap(TLweSample64*& key_a, TLweSample64*& val_a,
                                TLweSample64*& key_b, TLweSample64*& val_b,
                                const ControlGateBootstrappingKey64* cb_bsk,
                                AutoKsKeyFFTa** aks_fft, const GLevCipher64* SSK, int l_ss,
                                const ControlOutputKeySwitchKey64* ksk,
                                const Globals* env_pbs, const Globals* env)
{
    const int N = env->N;
    const int kpl = (env->k + 1) * env->l;
    const bool val_is_key = (val_a == key_a);

    // Extract the two encrypted sort keys as LWE ciphertexts and compute their difference
  // The sign of (key_a - key_b) determines whether the two entries should be swapped
    LweSample64* ka = new LweSample64(N);
    LweSample64* kb = new LweSample64(N);
    SampleExtract_AtZero64(ka, key_a, env);
    SampleExtract_AtZero64(kb, key_b, env);
    LweSample64* diff = new LweSample64(N);
    CopyLweSample64(diff, ka, env);
    SubLweSample64(diff, kb, env);

    // key-switch N -> n
    LweSample64* diff_n_dim = new LweSample64(N);
    KS_controlbit64(diff_n_dim, diff, ksk, env);
    TGswSample64* swap_condition_bsk = new TGswSample64(env->l, N);
    // obtain RGSW swap bit : result of comparison
    CircuitBTS64_ss_mvfbs(swap_condition_bsk, diff_n_dim, cb_bsk, aks_fft, SSK, l_ss, env_pbs, env);
    TGswSampleFFTa* swap_condition_fft = new TGswSampleFFTa(env->l, N);
    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= env->k; ++q)
            TorusPolynomial64_ifft_lvl2(&swap_condition_fft->allsamples[p].a[q],
                                        &swap_condition_bsk->allsamples[p].a[q], env);

    // obliviously reorder the two sort keys based on RGSW swap bit
    TLweSample64* na = new TLweSample64(N);
    TLweSample64* nb = new TLweSample64(N);
    { TLweSample64* t1 = new TLweSample64(N); 
      CopyTLweSample64(t1, key_b, env);
      CMuxFFT(na, swap_condition_fft, key_a, t1, env); delete t1; }
    { TLweSample64* t1 = new TLweSample64(N); 
      CopyTLweSample64(t1, key_a, env);
      CMuxFFT(nb, swap_condition_fft, key_b, t1, env); delete t1; }

    // Move the corresponding payloads using the same RGSW
    if (!val_is_key) {
        TLweSample64* ma = new TLweSample64(N);
        TLweSample64* mb = new TLweSample64(N);
        { TLweSample64* t1 = new TLweSample64(N); 
          CopyTLweSample64(t1, val_b, env);
          CMuxFFT(ma, swap_condition_fft, val_a, t1, env); delete t1; }
        { TLweSample64* t1 = new TLweSample64(N); 
          CopyTLweSample64(t1, val_a, env);
          CMuxFFT(mb, swap_condition_fft, val_b, t1, env); delete t1; }
        CopyTLweSample64(val_a, ma, env); CopyTLweSample64(val_b, mb, env);
        delete ma; delete mb;
    }
    CopyTLweSample64(key_a, na, env);
    CopyTLweSample64(key_b, nb, env);

    delete na; delete nb;
    delete swap_condition_fft; delete swap_condition_bsk;
    delete diff_n_dim; delete diff; delete kb; delete ka;
}

// Oblivious Sorting : the compared positions are fixed and public, while all comparison outcomes and data movements remain encrypted.
// The server cannot learn the ordering before or after the oblivious sorting.
void ObliviousSorting(std::vector<TLweSample64*>& K, std::vector<TLweSample64*>& V,
                               const ControlGateBootstrappingKey64* cb_bsk,
                               AutoKsKeyFFTa** aks_fft, const GLevCipher64* SSK, int l_ss,
                               const ControlOutputKeySwitchKey64* ksk,
                               const Globals* env_pbs, const Globals* env)
{
    HomCompareSwap(K[0],V[0], K[1],V[1], cb_bsk, aks_fft, SSK, l_ss, ksk, env_pbs, env);
    HomCompareSwap(K[1],V[1], K[2],V[2], cb_bsk, aks_fft, SSK, l_ss, ksk, env_pbs, env);
    HomCompareSwap(K[0],V[0], K[1],V[1], cb_bsk, aks_fft, SSK, l_ss, ksk, env_pbs, env);
}


