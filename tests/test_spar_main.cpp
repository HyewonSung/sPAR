/*
 [Implementation code of sPAR]
 - Paper title: sPAR :  (Somewhat) Practical Anonymous Router
 - Authors: Debajyoti Das, Jeoungeun Park and Hyewon Sung
 
 This is the main test code for sPAR.
 It executes one complete round of sPAR for eta clients (one message each).
 
 Main test file proceeds as follows :
 1. MPC-based setup : eta clients jointly generate the common pk and all evaluation keys.
                      Each client keeps local secret key share s_i of master sk s.
 2. Client:Send : Each client encrypts its payload and its three candidate bucket indices (one-hot RLWE queries) under the common pk.
 3. Server:Write : The server extracts the encrypted selector bits from the queries (SampleExtract + key switching) 
                   and runs HomPlacing (Algorithm 1) to insert each message into the least-loaded of its three candidate buckets.
 4. Sorting : Each bucket is sorted obliviously
 5. Decryption : The eta clients partially decrypt the bucket array and the server combines the shares to recover the messages.
 
 If you run this program, you can check a per-step timing breakdown, the end-to-end round latency, and a correctness result(per-bucket multiset match).
 For more details, please refer to README.md.

 ---------------------------------------------------------------------------
                                     BUILD
 ---------------------------------------------------------------------------
 The number of clients is a compile-time constant, so the program is rebuilt for each eta. 
 Set ETA and PARTY_COUNT to the same value.
 
 There are two ways of building this file, before running.
 One is for running the file with all eta = 20, 40, 60, 80, 100 and the other one is for just one eta.
 For example,
 - ./build.sh              # eta = 20, 40, 60, 80, 100   (~20 s)
 - ./build.sh 20           # just one                    

 which is equivalent to:
 -   mkdir -p build

 -   OBJ="spqlios/spqlios-fft-fma.o spqlios/spqlios-ifft-fma.o \
          spqlios/spqlios-fft-impl.o spqlios/fft_processor_spqlios.o \
          spqlios/lagrangehalfc_impl_fma.o"
 
 -   for E in 20 40 60 80 100; do
       g++ -std=c++2a -O3 -w -DUSE_FFT -march=native \
           -DETA=$E -DPARTY_COUNT=$E \
           -Iinclude -I. -Ispqlios -fopenmp -no-pie \
           $OBJ src/global_random.cpp tests/test_spar_main.cpp \
           -o build/spar_e$E
     done
  
 ---------------------------------------------------------------------------
                                RUN
 ---------------------------------------------------------------------------
 After building the project, an executable run-file for each value of eta should be generated. 
 Each runfile can be run as shown below.

 Note that the number of OpenMP threads can be specified using OMP_NUM_THREADS
 The measurements reported in the paper were obtained using 4 OpenMP threads, so to reproduce the results in the paper, run:

 -   OMP_NUM_THREADS=4 ./build/spar_e20      # ~4-5  minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e40      # ~16-17 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e60      # ~36-37 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e80      # ~64-65 minutes
 -   OMP_NUM_THREADS=4 ./build/spar_e100     # ~102-103 minutes
 
 An optional command-line argument overrides the RNG seed:
 -  ./build/spar_e20 12345

 Each run ends with "======== PASS ========" if the round was correct.

 VERBOSITY
 The print_debug flag defaults to FALSE:
 the parameter header, progress of per protocol step, and the three types of summary (failure estimate, END-TO-END, summary).
 
 Build with PRINT_DEBUG=1 to get the per-message and per-bucket detail:
 -   EXTRA_FLAGS=-DPRINT_DEBUG=1 ./build.sh 20

 
 ---------------------------------------------------------------------------
                                PARAMETERS
 ---------------------------------------------------------------------------
 The cryptographic parameters are fixed in the config macro block below and match the table in the paper. 
 Every LWE/RLWE instance arising in the implementation provides at least 117 bits of security (lattice estimator):
 
    GLWE  (N=2048, sigma=2^-54)                     119.7 bit
    small key / output KSK (n=1024, sigma=2^-26)    117.7 bit
    partial decryption share (N=2048, sigma=2^-54)  122.5 bit
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



// print debug flag and is set to FALSE as default
static const bool print_debug = (PRINT_DEBUG != 0);
#define DBG(...) do { if (print_debug) printf(__VA_ARGS__); } while (0)

static int flog2(uint64_t x) { int o = -1; while (x) { ++o; x >>= 1; } return o; }
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
    const int    ks_out_lv   = 64 / ks_out_bb;    // -> levels

    const int    sigma_big_log = -54;             
    const int    sigma_ksk_log = SIGMA_KSK_LOG; // small-key output-KSK noise 
    const double sigma_big   = std::pow(2.0, (double)sigma_big_log);
    const double sigma_ksk   = std::pow(2.0, (double)sigma_ksk_log);
    const int    ks_bb = CB_KS_BB, ks_lv = 64 / CB_KS_BB;  // CB private-KS gadget (ks_bb*ks_lv=64)

    Globals* env_out = new Globals(2048, 9, 3, 8, 256, 2048, 1); // gate/state/output ring
    Globals* env_pbs = new Globals(2048, 9, CB_PBS_L, 8, 256, 2048, 1); // CB blind-rotation gadget
    const int N = env_out->N;
    const int t = 16; //plaintext modulus
    const Torus64 Delta = (Torus64)(UINT64_C(1) << 60);
    const int PW = 8;
    auto decode = [&](Torus64 x)->int { uint64_t u=(uint64_t)x+(uint64_t)(Delta/2); return (int)((u>>60)&(uint64_t)(t-1)); };

    printf("======== sPAR: one protocol round ========\n");
    printf("N=%d  n(gate PBS)=%d  parties=%d eta=%d slots=%d choices=%d\n", N, lwe_n, party_count, eta, slots, choices);
    printf("refresh hasWritten every %d cand, I every %d msgs\n", REFRESH_HW, REFRESH_I);
    printf("omp_max_threads=%d\n", omp_get_max_threads());
    printf("output-KSK: basebit=%d levels=%d ; sigma_ksk=2^%d (abs 2^%d)\n",
           ks_out_bb, ks_out_lv, sigma_ksk_log, sigma_ksk_log + 64);
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
        env_out->tlwekey->coefs[i]=env_pbs->tlwekey->coefs[i]=s->coefs[i];
        env_out->lwekey[i]=env_pbs->lwekey[i]=s->coefs[i];
        env_out->in_key[i]=env_pbs->in_key[i]=s->coefs[i]; 
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

    DBG("\n[keygen] gate bsk + N->n finer output KSK(bb=%d) + CB bsk(l=%d) + HomTrace keys(aks_l=%d) + scheme-switch key(l_ss=%d) ...\n", ks_out_bb, env_pbs->l, CB_AKS_L, CB_L_SS);
    //setup for bts
    BootstrapContext64 boot; 
    Config_BTSContext64(&boot); 
    Attach_BTSEnvToContext64(&boot, env_out);
    ControlGateBootstrappingKey64 gate_bsk;
    Gen_GateBTSKey64(&gate_bsk, parties, lwe_n, sigma_big, env_out);
    ControlOutputKeySwitchKey64 output_ksk;
    Gen_OutputKSK_Basebit64(&output_ksk, parties, N, lwe_n, ks_out_bb, sigma_ksk, env_out);
    Attach_GatePBSKeysToBTSContext64(&boot, &gate_bsk, &output_ksk);

    ControlGateBootstrappingKey64 cb_bsk;
    Gen_GateBTSKey64(&cb_bsk, parties, lwe_n, sigma_big, env_pbs);
    
    // 2024/323 conversion keys: HomTrace automorphism-KS keys + one scheme-switch key RLWE'(s^2)
    const int aks_l = CB_AKS_L, l_ss = CB_L_SS;
    AutoKsKey64**   aks     = AutoKsKeyGenAll(env_out, 1, aks_l);
    AutoKsKeyFFTa** aks_fft = AutoKsKeyAllToFFT(aks, env_out, 1);
    IntPolynomiala* s2 = new IntPolynomiala(N); 
    NegacyclicSquare(s2, s, N);
    GLevCipher64* SSK = new GLevCipher64(l_ss, N); 
    GLevEncryptPoly(SSK, s2, sigma_big, env_out);

    HomPlacingState64* st = CreateHomPlacingState64(eta, slots, &boot, env_out);
    const long setup_ms = ms(t0);
    if (print_debug) 
        printf("  setup done in %ld ms\n", setup_ms);
    else             
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
        
        // Coefficient 0 is the sort key
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
    DBG("\n[Server:Write] starts\n");
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
        DBG("\n[msg %d] choices = {%d,%d,%d} (distinct)\n", c, choice[c][0], choice[c][1], choice[c][2]);

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
        const int ncand = slots*choices*eta;
        DBG("  HomPlacing %ld ms   (candidates=%d, %.1f ms/candidate)\n", hp, ncand, (double)hp/ncand);
        DBG("    p1a  AND(z,I)          parallel : %8.0f ms\n", T.p1a_and_ms);
        DBG("    p1b  AND(zI,NOT hw)    serial   : %8.0f ms\n", T.p1b_and_ms);
        DBG("    p2   CB + ExtProd      parallel : %8.0f ms   [CPU: CB %.0f + ExtProd %.0f]\n", T.p2_ms, T.p2_cb_cpu_ms, T.p2_extprod_cpu_ms);
        DBG("    p3   accumulate into L serial   : %8.0f ms\n", T.p3_add_ms);

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
        if (!placed) DBG(" OVERFLOW (should not happen)\n");

        for (int i=0;i<eta;++i) 
            for (int d=0; d<choices; ++d) 
                delete zown[i][d];
        for (int d=0; d<choices; ++d) 
            delete query[c][d];
        delete payload_ct[c];
    }
    if (!print_debug) printf("[server:write] done  (selector extract %.1f s + HomPlacing %.1f s)\n", selector_ms/1000.0, total_hp_ms/1000.0);

    DBG("\n[sort] oblivious sort of each bucket\n");
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
    if (print_debug) 
        printf("  sort done in %ld ms  (%d buckets on %d threads)\n", sort_ms, eta, omp_get_max_threads());
    else             
        printf("[sorting] done in %ld ms  (%d buckets on %d threads)\n", sort_ms, eta, omp_get_max_threads());
    DBG("  per bucket (3 slots, serial)      : %.3f s\n", sort_bucket_ms/1000.0);

    DBG("\n[verify] MPPartialDecrypt -> Combine -> FinalPhase; compute the noise per-bucket\n");
    auto tdec = Clock::now();
    double pdec_user_ms = 0.0;   
    double combine_ms   = 0.0;   
    bool ok=true;
    bool order_ok = true; 
    int multi_buckets = 0;
    long double nz_sumsq=0.0L; 
    unsigned long nz_cnt=0; 
    uint64_t nz_max=0;   
    Torus64Polynomial* ph=new Torus64Polynomial(N);
    const int nparties = (int)parties.size();
    std::vector<Torus64Polynomial*> dshare(nparties);
    for (int u=0;u<nparties;++u) 
        dshare[u]=new Torus64Polynomial(N);
    Torus64Polynomial* dsum = new Torus64Polynomial(N);
    for (int i=0;i<eta;++i) {
        std::vector<std::vector<int>> got, want;
        std::vector<uint64_t> cell_noise(slots, 0);
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
            uint64_t w=0;
            for (int j=0;j<N;++j){ 
                int sym=decode(ph->coefs[j]);
                if (j>=PW && sym!=0) 
                    ok=false;
                uint64_t nm = tdist(ph->coefs[j], (Torus64)((uint64_t)sym*(uint64_t)Delta));
                w=std::max(w, nm);
  
                long double e=(long double)nm; 
                    nz_sumsq += e*e; 
                    nz_cnt++;
                if (nm>nz_max) 
                    nz_max=nm; 
            }
            cell_noise[k]=w;
            want.push_back(std::vector<int>(pmsg[i][k].begin(), pmsg[i][k].begin()+PW));
        }
        std::vector<std::vector<int>> got_s=got, want_s=want;
        std::sort(got_s.begin(), got_s.end());
        std::sort(want_s.begin(), want_s.end());
        bool bucket_ok = (got_s==want_s);
        ok = ok && bucket_ok;
        uint64_t worst=0; 
        for (int k=0;k<slots;++k) 
            worst=std::max(worst, cell_noise[k]);
        DBG("  bucket %d  L-slots(after sort): ", i);
        if (print_debug) 
            for (int k=0;k<slots;++k) 
                printf("[%d%d%d%d] ", got[k][0],got[k][1],got[k][2],got[k][3]);
        DBG("| answer: ");
        if (print_debug) 
            for (auto& v: want_s) 
                printf("[%d%d%d%d] ", v[0],v[1],v[2],v[3]);

        bool asc=true; 
        for (int k=1;k<slots;++k) 
            if (got[k][0] < got[k-1][0]) 
                asc=false;
        int nonzero=0; 
        for (int k=0;k<slots;++k) 
            if (got[k][0]!=0 || got[k][1]!=0) nonzero++;
        order_ok = order_ok && asc; 
        if (nonzero>=2) multi_buckets++;
        DBG("| MATCH ? %s | noise~2^%d (Delta/2=2^59)\n",bucket_ok?"OK":"FAIL", worst?flog2(worst):-1);
    }
    delete ph; 
    delete dsum; 
    for (int u=0;u<nparties;++u) delete dshare[u];
    if (!print_debug) printf("[partial decryption] done in %ld ms\n", ms(tdec));
    (void)multi_buckets;
    ok = ok && order_ok;


    {
        // Noise budget
        // decryption boundary Delta/2
        const double sigma  = (nz_cnt ? (double)sqrtl(nz_sumsq/(long double)nz_cnt) : 0.0);
        const double margin = (nz_max ? 59.0 - log2((double)nz_max) : 59.0);
        printf("\n---- noise budget (payload) ------------------------------------------\n");
        printf("    samples (coefs)            : %lu\n", nz_cnt);
        printf("    sigma (RMS noise)          : 2^%.1f\n", sigma>0 ? log2(sigma) : 0.0);
        printf("    max |noise| observed       : 2^%d\n", nz_max ? flog2(nz_max) : 0);
        printf("    boundary Delta/2           : 2^59\n");
        printf("    margin                     : %.1f bits\n", margin);
    }
    printf("\n  avg HomPlacing/msg = %ld ms\n", total_hp_ms / party_count);
    const long dec_ms = ms(tdec);
    const double e2e_ms = now_ms() - t_e2e0;
    const int ncand_round = slots*choices*eta*party_count;

    printf("\n================ END-TO-END (eta=%d, %d messages) ================\n", eta, party_count);
    printf("  setup (keygen, one-time)          : %8.1f s\n", setup_ms/1000.0);
    printf("  ------------------------------------------------------\n");
    printf("  client:send (per client, parallel): %8.3f s\n", (client_send_ms/1000.0)/party_count);
    printf("  selector extract (all messages)   : %8.1f s   (server, 3*eta KS/msg, parallel)\n", selector_ms/1000.0);
    printf("  HomPlacing       (all messages)   : %8.1f s\n", total_hp_ms/1000.0);
    printf("      p1a  AND(z,I)        parallel : %8.1f s  (%4.1f%%)\n", agg.p1a_and_ms/1000.0, 100.0*agg.p1a_and_ms/total_hp_ms);
    printf("      p1b  AND(zI,NOT hw)  serial   : %8.1f s  (%4.1f%%)\n", agg.p1b_and_ms/1000.0, 100.0*agg.p1b_and_ms/total_hp_ms);
    printf("      p2   CB + ExtProd    parallel : %8.1f s  (%4.1f%%)   [CPU: CB %.1f s + ExtProd %.1f s]\n",
           agg.p2_ms/1000.0, 100.0*agg.p2_ms/total_hp_ms, agg.p2_cb_cpu_ms/1000.0, agg.p2_extprod_cpu_ms/1000.0);
    printf("      p3   accumulate      serial   : %8.1f s  (%4.1f%%)\n", agg.p3_add_ms/1000.0, 100.0*agg.p3_add_ms/total_hp_ms);
    printf("  oblivious sort                    : %8.1f s\n", sort_ms/1000.0);
    printf("  threshold decrypt + verify        : %8.1f s\n", dec_ms/1000.0);
    printf("  ======================================================\n");
    printf("  end-to-end (excluding setup)      : %8.1f s\n", e2e_ms/1000.0);
    printf("  per message                       : %8.2f s\n", (total_hp_ms/1000.0)/party_count);
    printf("  per candidate                     : %8.1f ms   (%d candidates in round)\n",
           (double)total_hp_ms/ncand_round, ncand_round);
    printf("  threads                           : %8d\n", omp_get_max_threads());
    printf("\n---- summary (eta=%d) -------------------------------------------------\n", eta);
    printf("  Client:Send      per message          : %.3f s\n", (client_send_ms/1000.0)/party_count);
    printf("  Server:Write     per message          : %.3f s\n", (total_hp_ms/1000.0)/party_count);
    printf("  Partial Dec.     per client            : %.3f s   (3*eta ciphertexts)\n", pdec_user_ms/1000.0);
    printf("  Partial Dec.     server combine        : %.3f s   (all eta shares)\n", combine_ms/1000.0);
    printf("  Sorting          3 slots              : %.3f s\n", sort_bucket_ms/1000.0);

    printf("\n======== %s ========\n", ok?"PASS":"FAIL");
    return ok?0:1;
}
