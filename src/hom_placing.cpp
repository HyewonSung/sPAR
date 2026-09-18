#include "include/config.h"
#include "include/hom_placing.h"

#include <chrono>
#include <cstdint>
#include <vector>

#include "include/bootstrapping.h"
#include "include/mphe.h"
#include "include/spar_helper.h"
#include "src/64header.h"

namespace {

int CeilLog2Pow2Local(int value) {
    int out = 0;
    int cur = 1;
    while (cur < value) {
        cur <<= 1;
        ++out;
    }
    return out;
}

}  

HomPlacingState64* CreateHomPlacingState64(int eta,
                                           int slots,
                                           const BootstrapContext64* boot,
                                           const Globals* env) {
    if (eta <= 0 || slots <= 0 || !boot) return nullptr;

    HomPlacingState64* st = new HomPlacingState64();
    st->eta = eta;
    st->slots = slots;
    st->choice_bitlen = CeilLog2Pow2Local(eta);
    st->L.assign(eta, std::vector<TLweSample64*>(slots, nullptr));
    st->I.assign(eta, std::vector<LweSample64*>(slots, nullptr));

    for (int i = 0; i < eta; ++i) {
        for (int k = 0; k < slots; ++k) {
            st->L[i][k] = new TLweSample64(env->N);
            st->I[i][k] = new LweSample64(env->N);
            MakeDataZero64(st->L[i][k], env);
            ENC_ControlBitConstant64(st->I[i][k], 1, boot, env);
        }
    }

    return st;
}

void DeleteHomPlacingState64(HomPlacingState64* st) {
    if (!st) return;
    for (int i = 0; i < st->eta; ++i) {
        for (int k = 0; k < st->slots; ++k) {
            delete st->L[i][k];
            delete st->I[i][k];
        }
    }
    delete st;
}



// ---- Algorithm 1: homomorphic placing --------------------------------------
//   Phase 1a  (parallel): zI[t] = AND(z, I)               -independent across candidates
//   Phase 1b  (serial)  : h[t] = AND(zI, NOT hasWritten); -update I/hasWritten; refresh them periodically
//   Phase 2   (parallel): write[t] = CircuitBootstrap(h[t])*(payload ciphertext)   -independent
//   Phase 3   (serial)  : L[i][k] += write[t]             -cheap adds, order-independent
void HomPlacingPhased(
    HomPlacingState64* st,
    const TLweSample64* payload_ct,
    const std::vector<std::vector<const LweSample64*>>& z_onehot,
    const BootstrapContext64* boot,
    const ControlGateBootstrappingKey64* cb_bsk,
    AutoKsKeyFFTa** aks_fft,
    const GLevCipher64* SSK, int l_ss,
    const Globals* env_pbs,
    const Globals* env,
    StageTimes* T)
{
    const double t_fn0 = now_ms();
    const int eta = st->eta, slots = st->slots;
    const int num_choices = (int)z_onehot[0].size();
    const int N = env->N;
    const int total_loop = slots * num_choices * eta;
    // linear index t = ((k*num_choices)+d)*eta + i  ->  iterating t in order = the (k,d,i) loop order
    auto K_loop = [&](int t){ return t / (num_choices*eta); };
    auto D_loop = [&](int t){ return (t / eta) % num_choices; };
    auto I_loop = [&](int t){ return t % eta; };

    std::vector<LweSample64*> zI(total_loop), h(total_loop);
    std::vector<TLweSample64*> write(total_loop);
    for (int t=0;t<total_loop;++t){ zI[t]=new LweSample64(N); h[t]=new LweSample64(N); write[t]=new TLweSample64(N); }

    // ---- Phase 1a (PARALLEL): zI = AND(z, I) ----
    const double t_1a = now_ms();
    #pragma omp parallel for schedule(dynamic)
    for (int t=0;t<total_loop;++t)
        PBS_ANDbit64(zI[t], z_onehot[I_loop(t)][D_loop(t)], st->I[I_loop(t)][K_loop(t)], boot, env);
    T->p1a_and_ms = now_ms() - t_1a;

    // ---- Phase 1b (SERIAL): h = AND(zI, NOT hasWritten) ; state update ; refresh them ----
    const double t_1b = now_ms();
    LweSample64* has_written = new LweSample64(N);
    LweSample64* control_zero = new LweSample64(N);
    LweSample64* control_one = new LweSample64(N);
    ENC_ControlBitConstant64(has_written, 0, boot, env);
    ENC_ControlBitConstant64(control_zero, 0, boot, env);
    ENC_ControlBitConstant64(control_one, 1, boot, env);
    for (int t=0;t<total_loop;++t) {
        const int i=I_loop(t), k=K_loop(t);
        LweSample64* not_has_written = new LweSample64(N);
        NOTbit64(not_has_written, has_written, boot, env);       // NOT(hasWritten)
        PBS_ANDbit64(h[t], zI[t], not_has_written, boot, env);       // h = AND(zI, NOT hasWritten)
        delete not_has_written;
        
        // I[i][k] -= h
        LweSample64* I_next = new LweSample64(N);
        CopyLweSample64(I_next, st->I[i][k], env);
        SubLweSample64(I_next, h[t], env);
        AddLweSample64(I_next, control_zero, env); // Re-center the ciphertext to the canonical bit encoding after subtraction.
        CopyLweSample64(st->I[i][k], I_next, env);
        delete I_next;
        
        // hasWritten += h
        LweSample64* hw = new LweSample64(N);
        CopyLweSample64(hw, has_written, env);
        AddLweSample64(hw, h[t], env);
        AddLweSample64(hw, control_one, env); // Re-center the ciphertext to the canonical bit encoding after addition.
        CopyLweSample64(has_written, hw, env);
        delete hw;
        
        // refresh of hasWritten (bootstrap AND with 1) -> resets accumulation
        if (REFRESH_HW > 0 && ((t+1) % REFRESH_HW) == 0) {
            LweSample64* hwr = new LweSample64(N);
            PBS_ANDbit64(hwr, has_written, control_one, boot, env);
            CopyLweSample64(has_written, hwr, env); delete hwr;
        }
    }
    delete has_written; delete control_zero; delete control_one;
    T->p1b_and_ms = now_ms() - t_1b;

    // ---- Phase 2 (PARALLEL): write[t] = CircuitBootstrap(h[t])*(payload ciphertext) ----
    const int kpl = (env->k + 1) * env->l;
    const double t_2 = now_ms();
    double cb_cpu = 0, extprod_cpu = 0;
    #pragma omp parallel for schedule(dynamic) reduction(+:cb_cpu,extprod_cpu)
    for (int t=0;t<total_loop;++t) {
        const double c0 = now_ms();
        TGswSample64* cb = new TGswSample64(env->l, N);
        CircuitBTS64_ss_mvfbs(cb, h[t], cb_bsk, aks_fft, SSK, l_ss, env_pbs, env); //cb=bootstrapped RGSW(h)
        const double c1 = now_ms(); cb_cpu += c1 - c0;
        TGswSampleFFTa* cb_fft = new TGswSampleFFTa(env->l, N);
        for (int p=0;p<kpl;++p) 
            for (int q=0;q<=env->k;++q)
                TorusPolynomial64_ifft_lvl2(&cb_fft->allsamples[p].a[q], &cb->allsamples[p].a[q], env);
        ExternalProductFFT(write[t], cb_fft, payload_ct, env);          // RLWE(h * payload)
        delete cb_fft; delete cb;
        extprod_cpu += now_ms() - c1;
    }
    T->p2_ms = now_ms() - t_2;
    T->p2_cb_cpu_ms = cb_cpu; T->p2_extprod_cpu_ms = extprod_cpu;

    // ---- Phase 3 (SERIAL): accumulate into the bucket array ----
    const double t_3 = now_ms();
    for (int t=0;t<total_loop;++t) 
        AddTLweSample64(st->L[I_loop(t)][K_loop(t)], write[t], env);
    T->p3_add_ms = now_ms() - t_3;

    for (int t=0;t<total_loop;++t)
        {   delete zI[t]; 
            delete h[t]; 
            delete write[t]; }
    T->total_ms = now_ms() - t_fn0;
}


