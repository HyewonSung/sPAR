/*
 ---------------------------------------------------------------------------
                                     BUILD
 ---------------------------------------------------------------------------
 There are two ways of building this file, before running.
 One is for running the file with all eta = 20, 40, 60, 80, 100 and the other one is for just one eta.
 For example,
 - ./build.sh              # eta = 20, 40, 60, 80, 100   (~20 s)
 - ./build.sh 20           # just one                    

 ---------------------------------------------------------------------------
                                RUN
 ---------------------------------------------------------------------------
 After building the project, an executable run-file for each value of eta should be generated. 
 Each runfile can be run as shown below.

 Note that the number of OpenMP threads can be specified using OMP_NUM_THREADS
 The measurements reported in the paper were obtained using 4 OpenMP threads, so to reproduce the results in the paper, run:

 -   OMP_NUM_THREADS=4 ./build/spar_e20      # ~3  minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e40      # ~10 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e60      # ~22 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e80      # ~41 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e100     # ~61 minutes

 Each run ends with "======== PASS ========" if the round was correct.
 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

#include "../src/mphe.cpp"
#include "../src/bootstrapping.cpp"
#include "../src/spar_helper.cpp"
#include "../src/hom_placing.cpp"
#include "../src/sorting.cpp"
#include "include/test_common64.h"

#include "include/config.h"


static uint64_t tabs(Torus64 x) { uint64_t a = (uint64_t)x, b = (uint64_t)(-x); return a < b ? a : b; }
static uint64_t tdist(Torus64 a, Torus64 b) { return tabs((Torus64)((uint64_t)a - (uint64_t)b)); }


int main(int argc, char** argv) {
    using Clock = std::chrono::steady_clock;
    auto ms = [](const Clock::time_point& t){ return std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now()-t).count(); };
    const unsigned seed = (argc > 1) ? (unsigned)std::atoi(argv[1]) : (unsigned)SEED;  // runtime seed override
    std::srand(seed);
    global_random->generator.seed(seed);

    const int    party_count = PARTY_COUNT; // messages
    const int    eta         = ETA; // buckets
    const int    slots       = 3;
    const int    choices     = 3; // choice-of-three
    const int    lwe_n       = LWE_N; // reduced LWE dimension 
    const int    ks_out_bb   = KS_OUT_BASEBIT; // output-KSK decomposition base bits
    const int    ks_out_lv   = KS_OUT_LEVELS;  // -> levels (top ks_out_lv*ks_out_bb bits)

    const int    sigma_big_log = SIGMA_LOG;     // encryption, evaluation keys, partial decryption
    const int    sigma_ksk_log = SIGMA_KSK_LOG; // small-key output-KSK noise
    const double sigma_big   = std::pow(2.0, (double)sigma_big_log);
    const double sigma_ksk   = std::pow(2.0, (double)sigma_ksk_log);

    Globals* env_out  = new Globals(2048, OUT_BG, OUT_L, 8, 256, 2048, 1); // state/output ring (CB output gadget)
    Globals* env_pbs  = new Globals(2048, CB_PBS_BG, CB_PBS_L, 8, 256, 2048, 1); // CB blind-rotation gadget
    Globals* env_gate = new Globals(2048, GATE_BG, GATE_L, 8, 256, 2048, 1); // gate-PBS blind-rotation gadget
    const int N = env_out->N;
    const int t = 16; //plaintext modulus
    const Torus64 Delta = (Torus64)(UINT64_C(1) << 60);
    const int PW = PAYLOAD_WIDTH; 
    auto decode = [&](Torus64 x)->int { uint64_t u=(uint64_t)x+(uint64_t)(Delta/2); return (int)((u>>60)&(uint64_t)(t-1)); };

    printf("======== sPAR: one protocol round ========\n");
    printf("N=%d  n(gate PBS)=%d  parties=%d eta=%d slots=%d choices=%d\n", N, lwe_n, party_count, eta, slots, choices);
    printf("refresh hasWritten every %d cand, I every %d msgs\n", REFRESH_HW, REFRESH_I);
    printf("omp_max_threads=%d\n", omp_get_max_threads());
    printf("sigma=2^%d (encryption, evaluation keys, partial decryption) ; payload=%d coefficients\n",
           sigma_big_log, PW);
    printf("output-KSK: basebit=%d levels=%d ; sigma_ksk=2^%d (abs 2^%d)\n",
           ks_out_bb, ks_out_lv, sigma_ksk_log, sigma_ksk_log + 64);
    printf("gadgets (bg,l): gate (%d,%d)  CB blind rotation (%d,%d)  CB output (%d,%d) ; aks_l=%d l_ss=%d\n",
           GATE_BG, GATE_L, CB_PBS_BG, CB_PBS_L, OUT_BG, OUT_L, CB_AKS_L, CB_L_SS);
    if (eta < choices) { printf("[error] need eta >= choices (%d) for distinct choice-of-three\n", choices); return 2; }

    // ---- setup ----
    auto t0 = Clock::now();
    std::vector<MPPartyKey64> parties;
    LweSample64* c_pk = new LweSample64(N); 
    TLweSample64* d_pk = new TLweSample64(N);
    Torus64* c_a = new Torus64[N]; 
    Torus64Polynomial* d_a = new Torus64Polynomial(N);
    MPSetupPartiesLweControl64(parties, c_pk, d_pk, c_a, d_a, party_count, sigma_big, env_out);  // all eta users share the key 

    IntPolynomiala* s = new IntPolynomiala(N);
    { 
        std::vector<IntPolynomiala*> sh(parties.size()); 
        for (size_t u=0;u<parties.size();++u) 
            sh[u]=parties[u].control_sk_share;
        MPAggregateSecretKey64(s, sh.data(), (int)sh.size(), env_out); 
    }
    for (int i=0;i<N;++i){ 
        env_out->tlwekey->coefs[i]=env_pbs->tlwekey->coefs[i]=env_gate->tlwekey->coefs[i]=s->coefs[i];
        env_out->lwekey[i]=env_pbs->lwekey[i]=env_gate->lwekey[i]=s->coefs[i];
        env_out->in_key[i]=env_pbs->in_key[i]=env_gate->in_key[i]=s->coefs[i]; 
    }

    TLweSample64* common_pk = new TLweSample64(N);
    { 
        Torus64Polynomial* a=new Torus64Polynomial(N); 
        SampleCommonA64(a, env_out);
        std::vector<TLweSample64*> pks(parties.size());
        for (size_t u=0;u<parties.size();++u){ 
            pks[u]=new TLweSample64(N); 
            MPKeyGenShare64(pks[u], a, parties[u].control_sk_share, sigma_big, env_out); 
        }
        MPAggregatePublicKey64(common_pk, pks.data(), (int)pks.size(), env_out);
        for (auto* p:pks) 
            delete p; 
        delete a; 
    }

    //setup for bts
    BootstrapContext64 boot; 
    Config_BTSContext64(&boot); 
    Attach_BTSEnvToContext64(&boot, env_gate);
    ControlGateBootstrappingKey64 gate_bsk;
    Gen_GateBTSKey64(&gate_bsk, parties, lwe_n, sigma_big, env_gate);
    ControlOutputKeySwitchKey64 output_ksk;
    Gen_OutputKSK_Basebit64(&output_ksk, parties, N, lwe_n, ks_out_bb, ks_out_lv, sigma_ksk, env_out);
    Attach_GatePBSKeysToBTSContext64(&boot, &gate_bsk, &output_ksk);

    ControlGateBootstrappingKey64 cb_bsk;
    Gen_GateBTSKey64(&cb_bsk, parties, lwe_n, sigma_big, env_pbs);
    
    // 2024/323 conversion keys: HomTrace automorphism-KS keys + one scheme-switch key RLWE'(s^2)
    const int aks_l = CB_AKS_L, l_ss = CB_L_SS;
    AutoKsKey64**   aks     = AutoKsKeyGenAll(env_out, 1, aks_l, sigma_big);
    AutoKsKeyFFTa** aks_fft = AutoKsKeyAllToFFT(aks, env_out, 1);
    IntPolynomiala* s2 = new IntPolynomiala(N); 
    NegacyclicSquare(s2, s, N);
    GLevCipher64* SSK = new GLevCipher64(l_ss, N); 
    GLevEncryptPoly(SSK, s2, sigma_big, env_out);

    HomPlacingState64* st = CreateHomPlacingState64(eta, slots, &boot, env_out);
    const long setup_ms = ms(t0);
    printf("\n[keygen] setup done in %ld ms\n", setup_ms);

    std::vector<std::vector<int>> pempty(eta, std::vector<int>(slots, 1));
    std::vector<std::vector<std::vector<int>>> pmsg(eta, std::vector<std::vector<int>>(slots, std::vector<int>(PW, 0)));

    // control_one for periodic I refresh across messages
    LweSample64* refresh_one = new LweSample64(N); 
    ENC_ControlBitConstant64(refresh_one, 1, &boot, env_out);

    long total_hp_ms = 0, selector_ms = 0;
    double client_send_ms = 0.0;
    StageTimes agg;
    const double t_e2e0 = now_ms();

    // Client:Send : every client encrypts (payload + 3 queries)
    std::vector<std::vector<int>> choice(party_count, std::vector<int>(choices, 0));
    std::vector<std::vector<int>> symv(party_count, std::vector<int>(PW, 0));
    std::vector<TLweSample64*> payload_ct(party_count, nullptr);
    std::vector<std::vector<TLweSample64*>> query(party_count, std::vector<TLweSample64*>(choices, nullptr));
    const double t_send = now_ms();
    for (int c = 0; c < party_count; ++c) {
        int a[3]; 
        ChoiceofDistinctThree(a, choices, eta);   // 3 distinct buckets
        for (int d=0; d<choices; ++d) 
            choice[c][d] = a[d];

        Torus64Polynomial* payload = new Torus64Polynomial(N);
        for (int j=0;j<N;++j) payload->coefs[j]=0;
        
        
        for (int j=0;j<PW;++j){
            symv[c][j] = (j==0) ? (((c+1)*3+1) % CMP_HALF_RANGE) : (((c+1)*3+j+1) % t);
            payload->coefs[j]=(Torus64)((uint64_t)symv[c][j]*(uint64_t)Delta);
        }
        payload_ct[c] = new TLweSample64(N);
        MP_RLWE_Encrypt64(payload_ct[c], common_pk, payload, sigma_big, env_out);
        for (int d=0; d<choices; ++d) {
            Torus64Polynomial* oh = new Torus64Polynomial(N);
            for (int j=0;j<N;++j) 
                oh->coefs[j]=EncodeControlBitToTorus64(0);
            oh->coefs[a[d]] = EncodeControlBitToTorus64(1);
            query[c][d] = new TLweSample64(N);
            MP_RLWE_Encrypt64(query[c][d], common_pk, oh, sigma_big, env_out);
            delete oh;
        }
        delete payload;
    }
    client_send_ms = now_ms() - t_send;
    printf("\n[client:send] done in %.0f ms\n", client_send_ms);

    // Server:Write : selector extraction + HomPlacing, one message at a time.
    for (int c = 0; c < party_count; ++c) {
        // STATE REFRESH: bootstrap every I every REFRESH_I messages 
        if (REFRESH_I > 0 && c > 0 && (c % REFRESH_I) == 0) {
            for (int i=0;i<eta;++i) 
                for (int k=0;k<slots;++k) {
                    LweSample64* ir = new LweSample64(N);
                    PBS_ANDbit64(ir, st->I[i][k], refresh_one, &boot, env_out);
                    CopyLweSample64(st->I[i][k], ir, env_out);
                    delete ir;
                }
        }
        // Server-side selector extraction (parallel)
        auto tq = Clock::now();
        std::vector<std::vector<LweSample64*>> zown(eta, std::vector<LweSample64*>(choices, nullptr));
        std::vector<std::vector<const LweSample64*>> zview(eta, std::vector<const LweSample64*>(choices, nullptr));
        #pragma omp parallel for collapse(2) schedule(dynamic)
        for (int i=0;i<eta;++i) for (int d=0; d<choices; ++d) {
            LweSample64* z_full = new LweSample64(N);
            ExtractRLWEtoLWE(z_full, query[c][d], i, env_out);
            zown[i][d] = new LweSample64(N);
            KS_controlbit64(zown[i][d], z_full, &output_ksk, env_out);
            zview[i][d] = zown[i][d];
            delete z_full;
        }
        selector_ms += ms(tq);

        auto tp = Clock::now();
        StageTimes T;
        HomPlacingPhased(st, payload_ct[c], zview, &boot, &cb_bsk, aks_fft, SSK, l_ss, env_pbs, env_out, &T);
        long hp = ms(tp); total_hp_ms += hp; agg.add(T);

        bool placed=false;
        for (int k=0;k<slots && !placed;++k) 
            for (int d=0; d<choices && !placed; ++d) {
                int b=choice[c][d]; if (pempty[b][k]) { 
                    pempty[b][k]=0;
                    for (int j=0;j<PW;++j) 
                        pmsg[b][k][j]=symv[c][j]; 
                    placed=true; 
                } 
            }
        if (!placed) printf("[warning] message %d overflows its three buckets (should not happen)\n", c);

        for (int i=0;i<eta;++i) 
            for (int d=0; d<choices; ++d) 
                delete zown[i][d];
        for (int d=0; d<choices; ++d) 
            delete query[c][d];
        delete payload_ct[c];
    }
    printf("[server:write] done  (selector extract %.1f s + HomPlacing %.1f s)\n", selector_ms/1000.0, total_hp_ms/1000.0);

    auto tsort = Clock::now();

    std::vector<double> bucket_ms(eta, 0.0);
    #pragma omp parallel for schedule(dynamic)
    for (int i=0;i<eta;++i) {
        const double b0 = now_ms();
        ObliviousSorting(st->L[i], st->L[i], &cb_bsk, aks_fft, SSK, l_ss, &output_ksk, env_pbs, env_out);
        bucket_ms[i] = now_ms() - b0;
    }
    const long sort_ms = ms(tsort);
    double sort_bucket_ms = 0.0; 
    for (int i=0;i<eta;++i) 
        sort_bucket_ms += bucket_ms[i];
    sort_bucket_ms /= eta;
    printf("[sorting] done in %ld ms  (%d buckets on %d threads)\n", sort_ms, eta, omp_get_max_threads());

    auto tdec = Clock::now();
    double pdec_user_ms = 0.0;   
    double combine_ms   = 0.0;   
    bool ok=true;
    bool order_ok = true; 
    int multi_buckets = 0;
    uint64_t nz_max=0;
    Torus64Polynomial* ph=new Torus64Polynomial(N);
    const int nparties = (int)parties.size();
    std::vector<Torus64Polynomial*> dshare(nparties);
    for (int u=0;u<nparties;++u) 
        dshare[u]=new Torus64Polynomial(N);
    Torus64Polynomial* dsum = new Torus64Polynomial(N);
    for (int i=0;i<eta;++i) {
        std::vector<std::vector<int>> got, want;
        for (int k=0;k<slots;++k) {
       
            const double pd0 = now_ms();
            MPPartialDecrypt64(dshare[0], st->L[i][k], parties[0].control_sk_share, sigma_big, env_out);
            pdec_user_ms += now_ms() - pd0;                 
            for (int u=1;u<nparties;++u)
                MPPartialDecrypt64(dshare[u], st->L[i][k], parties[u].control_sk_share, sigma_big, env_out);
            const double cb0 = now_ms();
            MPCombineShares64(dsum, dshare.data(), nparties, env_out);
            MPFinalPhase64(ph, st->L[i][k], dsum, env_out);
            combine_ms += now_ms() - cb0;                   
            std::vector<int> v(PW);
            for (int j=0;j<PW;++j) 
                v[j]=decode(ph->coefs[j]);
            got.push_back(v);
            for (int j=0;j<N;++j){
                int sym=decode(ph->coefs[j]);
                if (j>=PW && sym!=0)
                    ok=false;
                uint64_t nm = tdist(ph->coefs[j], (Torus64)((uint64_t)sym*(uint64_t)Delta));
                if (nm>nz_max)
                    nz_max=nm; 
            }
            want.push_back(std::vector<int>(pmsg[i][k].begin(), pmsg[i][k].begin()+PW));
        }
        std::vector<std::vector<int>> got_s=got, want_s=want;
        std::sort(got_s.begin(), got_s.end());
        std::sort(want_s.begin(), want_s.end());
        bool bucket_ok = (got_s==want_s);
        ok = ok && bucket_ok;

        bool asc=true;
        for (int k=1;k<slots;++k) 
            if (got[k][0] < got[k-1][0]) 
                asc=false;
        int nonzero=0; 
        for (int k=0;k<slots;++k) 
            if (got[k][0]!=0 || got[k][1]!=0) nonzero++;
        order_ok = order_ok && asc; 
        if (nonzero>=2) multi_buckets++;
    }
    delete ph; 
    delete dsum; 
    for (int u=0;u<nparties;++u) delete dshare[u];
    printf("[partial decryption] done in %ld ms\n", ms(tdec));
    (void)multi_buckets;
    ok = ok && order_ok;


    // noise budget: worst payload noise over all coefficients against the decryption boundary Delta/2
    const double noise_log = nz_max ? log2((double)nz_max) : 0.0;
    printf("\n[noise budget] max |noise| = 2^%.1f, Delta/2 = 2^59 -> margin %.1f bits\n", noise_log, 59.0 - noise_log);
    const long dec_ms = ms(tdec);
    const double e2e_ms = now_ms() - t_e2e0;

    printf("\n================ END-TO-END (eta=%d, %d messages) ================\n", eta, party_count);
    printf("  setup (keygen, one-time)          : %8.1f s\n", setup_ms/1000.0);
    printf("  ------------------------------------------------------\n");
    printf("  client:send (per client, parallel): %8.3f s\n", (client_send_ms/1000.0)/party_count);
    printf("  selector extract (all messages)   : %8.1f s   (server, 3*eta KS/msg, parallel)\n", selector_ms/1000.0);
    printf("  HomPlacing       (all messages)   : %8.1f s\n", total_hp_ms/1000.0);
    printf("      p1a  AND(z,I)        parallel : %8.1f s  (%4.1f%%)\n", agg.p1a_and_ms/1000.0, 100.0*agg.p1a_and_ms/total_hp_ms);
    printf("      p1b  AND(zI,NOT hw)  serial   : %8.1f s  (%4.1f%%)\n", agg.p1b_and_ms/1000.0, 100.0*agg.p1b_and_ms/total_hp_ms);
    printf("      p2   CB + ExtProd    parallel : %8.1f s  (%4.1f%%)\n", agg.p2_ms/1000.0, 100.0*agg.p2_ms/total_hp_ms);
    printf("      p3   accumulate      serial   : %8.1f s  (%4.1f%%)\n", agg.p3_add_ms/1000.0, 100.0*agg.p3_add_ms/total_hp_ms);
    printf("  oblivious sort                    : %8.1f s\n", sort_ms/1000.0);
    printf("  threshold decrypt + verify        : %8.1f s\n", dec_ms/1000.0);
    printf("  ======================================================\n");
    printf("  end-to-end (excluding setup)      : %8.1f s\n", e2e_ms/1000.0);
    printf("  threads                          : %8d\n", omp_get_max_threads());
    printf("\n---- summary (eta=%d) -------------------------------------------------\n", eta);
    printf("  Client:Send      per message          : %.3f s\n", (client_send_ms/1000.0)/party_count);
    printf("  Server:Write     per message          : %.3f s\n", (total_hp_ms/1000.0)/party_count);
    printf("  Partial Dec.     per client            : %.3f s   (3*eta ciphertexts)\n", pdec_user_ms/1000.0);
    printf("  Partial Dec.     server combine        : %.3f s   (all eta shares)\n", combine_ms/1000.0);
    printf("  Sorting          3 slots              : %.3f s\n", sort_bucket_ms/1000.0);

    printf("\n======== %s ========\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
