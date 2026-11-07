#include "include/bootstrapping.h"

#include <cassert>
#include <vector>

#include "src/64header.h"

namespace {

enum class ControlBinaryGate64 {
    And,
};

// |x| on the torus: min(x, 2^64 - x)
uint64_t TorusAbs64(Torus64 x) {
    const uint64_t ux = static_cast<uint64_t>(x);
    const uint64_t uy = static_cast<uint64_t>(-x);
    return ux < uy ? ux : uy;
}

// dist(lhs, rhs) = |lhs - rhs| on the torus
uint64_t TorusDist64(Torus64 lhs, Torus64 rhs) {
    const Torus64 diff = static_cast<Torus64>(static_cast<uint64_t>(lhs) - static_cast<uint64_t>(rhs));
    return TorusAbs64(diff);
}


// E(a) + E(b) = (2(a+b) - 2)*mu
// sum_class = a+b = 0,1,2 gives -2mu, 0, +2mu
Torus64 GateInputTarget(int sum_class) {
    const int64_t mu = static_cast<int64_t>(EncodeControlBitToTorus64(1));
    switch (sum_class) {
        case 0:
            return static_cast<Torus64>(-2 * mu);
        case 1:
            return 0;
        case 2:
            return static_cast<Torus64>(2 * mu);
    }
    return 0;
}

// log2(x), x a power of two
int Log2Pow2(int x) {
    assert(x > 0);
    assert((x & (x - 1)) == 0);

    int out = 0;
    while ((1 << out) < x) ++out;
    return out;
}

// slot -> phase: idx -> (idx + 1/2)*2^64/(2N)
Torus64 TorusRepAtIndex(int idx, int N) {
    const int bits = Log2Pow2(2 * N);
    const uint64_t cell = UINT64_C(1) << (64 - bits);
    return static_cast<Torus64>((static_cast<uint64_t>(idx) << (64 - bits)) + (cell >> 1));
}

//inverse of the above
// phase -> slot: x ->  round(x*2N / 2^64) mod 2N
int ModSwitch_2N(Torus64 x, int N, int step = 1) {
    const int slots = 2 * N / step;
    const int bits = Log2Pow2(slots);
    const int shift = 64 - bits;
    const uint64_t ux = static_cast<uint64_t>(x);
    const uint64_t rounding = (shift > 0) ? (UINT64_C(1) << (shift - 1)) : 0;
    return static_cast<int>(((ux + rounding) >> shift) & static_cast<uint64_t>(slots - 1)) * step;
}

// (a, b) = (0, E(bit)) 
void MakeTrivial_ControlBit(LweSample64* out,
                                        int bit,
                                        const Globals* env) {
    ClearLweSample64(out, env);
    *out->b = EncodeControlBitToTorus64(bit & 1);
}

// the gate truth table
// AND(a,b) = 1 iff a+b == 2
int EvalGate_FromClass(ControlBinaryGate64 gate, int sum_class) {
    switch (gate) {
        case ControlBinaryGate64::And:
            return (sum_class == 2) ? 1 : 0;
    }
    return 0;
}

// a+b = argmin_c dist(phase, target(c))
int ClassifyGatePhase(Torus64 phase) {
    int best_class = 0;
    uint64_t best_dist = TorusDist64(
        phase, GateInputTarget(0)
    );
    for (int sum_class = 1; sum_class <= 2; ++sum_class) {
        const uint64_t dist = TorusDist64(
            phase, GateInputTarget(sum_class)
        );
        if (dist < best_dist) {
            best_dist = dist;
            best_class = sum_class;
        }
    }
    return best_class;
}

// out = E(a) + E(b) 
void BuildGateInput(LweSample64* out,
                    const LweSample64* lhs_bit, const LweSample64* rhs_bit,
                    const Globals* env) {
    assert(out != nullptr);
    assert(lhs_bit != nullptr);
    assert(rhs_bit != nullptr);

    CopyLweSample64(out, lhs_bit, env);
    AddLweSample64(out, rhs_bit, env);
}

// s = sum_u share_u
void Recon_ControlSecret(IntPolynomiala* control_secret,
                        const std::vector<MPPartyKey64>& parties,
                        const Globals* env) {
    const int num_users = static_cast<int>(parties.size());
    std::vector<IntPolynomiala*> ptrs(num_users, nullptr);
    for (int u = 0; u < num_users; ++u) {
        ptrs[u] = parties[u].control_sk_share;
    }
    MPAggregateSecretKey64(control_secret, ptrs.data(), num_users, env);
}


// LWE encryption that only uses the first active_dim coordinates of the key (the rest of a is 0)
// cf) output KSK targets the smaller dimension n < N.
// out = (a,b) s.t. b = m + e + <a, s>
void ENC_ControlTorus_smalldim(LweSample64* out,
                               Torus64 message,
                               const IntPolynomiala* control_secret, int active_dim,
                               double stdev, const Globals* env) {
    assert(out != nullptr);
    assert(control_secret != nullptr);
    assert(active_dim >= 0 && active_dim <= env->N);

    Torus64 phase = message + random_gaussian64(0, stdev);
    for (int i = 0; i < active_dim; ++i) {
        out->a[i] = random_int64();
        phase += out->a[i] * control_secret->coefs[i];
    }
    for (int i = active_dim; i < env->N; ++i) {
        out->a[i] = 0;
    }
    *out->b = phase;
}

// out = -in
void NegateTLWE(TLweSample64* out,
                     const TLweSample64* in,
                     const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        for (int j = 0; j < env->N; ++j) {
            out->a[q].coefs[j] = -in->a[q].coefs[j];
        }
    }
}

// out = X^exponent * in  in Z[X]/(X^N + 1)
// cf) X^N = -1
void MulTLWEbyMonomial(TLweSample64* out,
                            const TLweSample64* in,
                            int exponent,
                            const Globals* env) {
    const int twoN = 2 * env->N;
    int exp_mod = exponent % twoN;
    if (exp_mod < 0) exp_mod += twoN;

    if (exp_mod < env->N) {
        shift(out, exp_mod, const_cast<TLweSample64*>(in), env->N);
        return;
    }

    TLweSample64* tmp = new TLweSample64(env->N);
    shift(tmp, exp_mod - env->N, const_cast<TLweSample64*>(in), env->N);
    NegateTLWE(out, tmp, env);
    delete tmp;
}

// phase += constant 
void AddConst_LWEmsg(LweSample64* ct, Torus64 constant) {
    *ct->b += constant;
}

// -x mod 2^64
static inline Torus64 TorusNeg64(Torus64 x) {
    return static_cast<Torus64>(UINT64_C(0) - static_cast<uint64_t>(x));
}


// Build the PBS accumulator for the selected Boolean gate.
void BuildGateAccumulator(TLweSample64* out,
                                            ControlBinaryGate64 gate,
                                            Torus64 phase_offset,
                                            const Globals* env) {
    Torus64Polynomial* poly = new Torus64Polynomial(env->N);
    for (int idx = 0; idx < env->N; ++idx) {
        //slot->phase
        const Torus64 rep = static_cast<Torus64>(static_cast<uint64_t>(TorusRepAtIndex(idx, env->N)) -static_cast<uint64_t>(phase_offset)); 
        const int sum_class = ClassifyGatePhase(rep); //phase of a+b
        const int out_bit = EvalGate_FromClass(gate, sum_class); //AND result
        poly->coefs[idx] = EncodeControlBitToTorus64(out_bit); 
    }
    tLwe64NoiselessTrivial(out, poly, env);
    delete poly;
}

// PBS : blind rotate -> sample extract -> key switch back to n
void RunPBS(LweSample64* out,
                               const TLweSample64* accumulator,
                               const LweSample64* gate_input,
                               const BootstrapContext64* ctx,
                               const Globals* env) {
    assert(out != nullptr);
    assert(accumulator != nullptr);
    assert(gate_input != nullptr);
    assert(ctx != nullptr);
    assert(ctx->gate_bsk != nullptr);
    assert(ctx->output_ksk != nullptr);

    TLweSample64* rotated = new TLweSample64 (env->N);
    LweSample64* extracted = new LweSample64(env->N);

    BlindRotate64( rotated, accumulator, gate_input, ctx->gate_bsk, env);
    SampleExtract_AtZero64(extracted, rotated, env);
    KS_controlbit64(out, extracted, ctx->output_ksk, env);

    delete rotated;
    delete extracted;
}

// decompose one torus value into signed base-2^basebit digits in [-Bg/2, Bg/2), keeping the top levels*basebit bits (rounded)
void Decomp_LWEscalar(int* digits,
                             Torus64 scalar,
                             int basebit,
                             int levels) {
    assert(digits != nullptr);
    assert(basebit > 0 && basebit < 32);
    assert(levels > 0);
    assert(levels * basebit <= 64);

    const uint64_t Bg = UINT64_C(1) << basebit;
    const int prec = levels * basebit;
    uint64_t state = static_cast<uint64_t>(scalar) + ((prec < 64) ? (UINT64_C(1) << (63 - prec)) : 0);
    for (int lev = 0; lev < levels; ++lev) state += (Bg / 2) << (64 - (lev + 1) * basebit);

    for (int lev = 0; lev < levels; ++lev) {
        const int shift = 64 - (lev + 1) * basebit;
        digits[lev] = static_cast<int>((state >> shift) & (Bg - 1)) - static_cast<int>(Bg / 2);
    }
}

}  

void Config_BTSContext64(BootstrapContext64* ctx) {
    assert(ctx != nullptr);
    ctx->gate_bsk = nullptr;
    ctx->output_ksk = nullptr;
    ctx->gate_env = nullptr;
}

void Attach_GatePBSKeysToBTSContext64(
    BootstrapContext64* ctx,
    const ControlGateBootstrappingKey64* gate_bsk,
    const ControlOutputKeySwitchKey64* output_ksk
) {
    assert(ctx != nullptr);
    ctx->gate_bsk = gate_bsk;
    ctx->output_ksk = output_ksk;
}

void Attach_BTSEnvToContext64(BootstrapContext64* ctx,
                                   const Globals* gate_env) {
    assert(ctx != nullptr);
    ctx->gate_env = gate_env;
}

void ENC_ControlBitConstant64(LweSample64* out,
                                       int bit,
                                       const BootstrapContext64* ctx,
                                       const Globals* env) {
    assert(out != nullptr);
    assert(ctx != nullptr);
    (void)ctx;
    MakeTrivial_ControlBit(out, bit, env);
}

void MakeDataZero64(TLweSample64* out, const Globals* env) {
    Torus64Polynomial* mu = new Torus64Polynomial(env->N);
    ClearTorusPolynomial64(mu, env);
    tLwe64NoiselessTrivial(out, mu, env);
    delete mu;
}

// BTS key: bsk[i] = RGSW.Enc(sk = s, msg = s_i) for i < input_dim
//  CMux selectors of blind rotation
void Gen_GateBTSKey64(
    ControlGateBootstrappingKey64* key,
    const std::vector<MPPartyKey64>& parties,
    int input_dim,
    double stdev,
    const Globals* env
) {
    assert(key != nullptr);
    assert(input_dim >= 0 && input_dim <= env->N);

    Del_GateBTSKey64(key);
    key->input_dim = input_dim;
    key->secret_bits.assign(input_dim, nullptr);
    key->secret_bits_fft.assign(input_dim, nullptr);

    IntPolynomiala* control_secret = new IntPolynomiala(env->N);
    Recon_ControlSecret(control_secret, parties, env);

    std::vector<int64_t> saved_tlwe(env->N, 0);
    std::vector<int> saved_lwe(env->N, 0);
    std::vector<int64_t> saved_in(env->N, 0);
    for (int i = 0; i < env->N; ++i) {
        saved_tlwe[i] = env->tlwekey->coefs[i];
        saved_lwe[i] = env->lwekey[i];
        saved_in[i] = env->in_key[i];
        env->tlwekey->coefs[i] = control_secret->coefs[i];
        env->lwekey[i] = control_secret->coefs[i];
        env->in_key[i] = control_secret->coefs[i];
    }

    for (int i = 0; i < input_dim; ++i) {
        key->secret_bits[i] = new TGswSample64(env->l, env->N);

        // Encrypt the i-th control-secret bit as an RGSW ciphertext
        tGsw64Encrypt(key->secret_bits[i], control_secret->coefs[i], stdev, env);
        // Precompute its Fourier representation for faster CMux
        key->secret_bits_fft[i] = new TGswSampleFFTa(env->l, env->N);
        const int kpl = (env->k + 1) * env->l;
        for (int p = 0; p < kpl; ++p) {
            for (int q = 0; q <= env->k; ++q) {
                TorusPolynomial64_ifft_lvl2( &key->secret_bits_fft[i]->allsamples[p].a[q], &key->secret_bits[i]->allsamples[p].a[q], env);
            }
        }
    }

    for (int i = 0; i < env->N; ++i) {
        env->tlwekey->coefs[i] = saved_tlwe[i];
        env->lwekey[i] = saved_lwe[i];
        env->in_key[i] = saved_in[i];
    }

    delete control_secret;
}

void Del_GateBTSKey64(ControlGateBootstrappingKey64* key) {
    assert(key != nullptr);
    for (TGswSample64* entry : key->secret_bits) {
        delete entry;
    }
    for (TGswSampleFFTa* entry_fft : key->secret_bits_fft) {
        delete entry_fft;
    }
    key->secret_bits.clear();
    key->secret_bits_fft.clear();
    key->input_dim = 0;
}

// Output KSK, N -> target_dim:
//   ksk[i][lev] = LWE.Enc(sk = s[0:target_dim], msg = s_i * 2^(64-(lev+1)*basebit))
void Gen_OutputKSK_TargetDim64(
    ControlOutputKeySwitchKey64* key,
    const std::vector<MPPartyKey64>& parties,
    int input_dim,
    int target_dim,
    double stdev,
    const Globals* env
) {
    assert(key != nullptr);
    assert(input_dim >= 0 && input_dim <= env->N);
    assert(target_dim >= 0 && target_dim <= env->N);
    assert(64 % env->basebit == 0);

    Del_OutputKSK64(key);

    key->input_dim = input_dim;
    key->target_dim = target_dim;
    key->levels = 64 / env->basebit;
    key->basebit = env->basebit;
    key->identity = false;
    key->entries.assign(input_dim, std::vector<LweSample64*>(key->levels, nullptr));

    IntPolynomiala* control_secret = new IntPolynomiala(env->N);
    Recon_ControlSecret(control_secret, parties, env);

    for (int i = 0; i < input_dim; ++i) {
        for (int lev = 0; lev < key->levels; ++lev) {
            key->entries[i][lev] = new LweSample64(env->N);
            //gadget
            const Torus64 scale = static_cast<Torus64>(UINT64_C(1) << (64 - (lev + 1) * key->basebit));
            const Torus64 message = control_secret->coefs[i] ? scale : 0;
            ENC_ControlTorus_smalldim(key->entries[i][lev],message,control_secret,target_dim,stdev,env);
        }
    }

    delete control_secret;
}

void Del_OutputKSK64(ControlOutputKeySwitchKey64* key) {
    assert(key != nullptr);
    for (auto& row : key->entries) {
        for (LweSample64* entry : row) {
            delete entry;
        }
    }
    key->entries.clear();
    key->input_dim = 0;
    key->target_dim = 0;
    key->levels = 0;
    key->basebit = 0;
    key->identity = true;
}

// Rotate the test vector by -phase(gate_input)
void BlindRotate64(
    TLweSample64* out_accum,
    const TLweSample64* testvec,
    const LweSample64* gate_input,
    const ControlGateBootstrappingKey64* key,
    const Globals* env,
    int ms_step
) {
    assert(out_accum != nullptr);
    assert(testvec != nullptr);
    assert(gate_input != nullptr);
    assert(key != nullptr);
    assert(key->input_dim >= 0 && key->input_dim <= env->N);
    assert(static_cast<int>(key->secret_bits_fft.size()) == key->input_dim);

    TLweSample64* accum = new TLweSample64(env->N);
    TLweSample64* rotated = new TLweSample64(env->N);
    TLweSample64* next = new TLweSample64(env->N);

    // Initial rotation corresponding to -b
    const int bbar = ModSwitch_2N(TorusNeg64(*gate_input->b), env->N, ms_step);
    MulTLWEbyMonomial(accum, testvec, bbar, env);

    // Conditionally add the rotation corresponding to each a_i s_i
    for (int i = 0; i < key->input_dim; ++i) {
        const int abar = ModSwitch_2N(gate_input->a[i], env->N, ms_step);
        MulTLWEbyMonomial(rotated, accum, abar, env);
        CMuxFFT(next, key->secret_bits_fft[i], accum, rotated, env);
        CopyTLweSample64(accum, next, env);
    }

    CopyTLweSample64(out_accum, accum, env);

    delete accum;
    delete rotated;
    delete next;
}

// Extract coefficient 0 of a TLWE/RLWE ciphertext as an LWE ciphertext
void SampleExtract_AtZero64(LweSample64* out,
                                  const TLweSample64* in_accum,
                                  const Globals* env) {
    assert(out != nullptr);
    assert(in_accum != nullptr);
    assert(env->k == 1);

    ClearLweSample64(out, env);
    out->a[0] = in_accum->a[0].coefs[0];
    for (int i = 1; i < env->N; ++i) {
        out->a[i] = -in_accum->a[0].coefs[env->N - i];
    }
    *out->b = in_accum->b->coefs[0];
}

// Dimension switch N -> n: gadget-decompose each a_i and subtract the matching KSK entries
void KS_controlbit64(
    LweSample64* out,
    const LweSample64* in_extracted,
    const ControlOutputKeySwitchKey64* key,
    const Globals* env
) {
    assert(out != nullptr);
    assert(in_extracted != nullptr);
    assert(key != nullptr);

    assert(key->input_dim >= 0 && key->input_dim <= env->N);
    assert(key->levels > 0);
    assert(static_cast<int>(key->entries.size()) == key->input_dim);

    ClearLweSample64(out, env);
    *out->b = *in_extracted->b;

    std::vector<int> digits(key->levels, 0);

    for (int i = 0; i < key->input_dim; ++i) {
        Decomp_LWEscalar(digits.data(), in_extracted->a[i], key->basebit, key->levels);

        for (int lev = 0; lev < key->levels; ++lev) {
            const int digit = digits[lev];
            if (digit == 0) continue;

            const uint64_t udigit = static_cast<uint64_t>(static_cast<int64_t>(digit));
            for (int j = 0; j < key->target_dim; ++j) {
                out->a[j] = static_cast<Torus64>(
                    static_cast<uint64_t>(out->a[j]) -
                                            udigit * static_cast<uint64_t>(key->entries[i][lev]->a[j])
                );
            }

            *out->b = static_cast<Torus64>(
                static_cast<uint64_t>(*out->b) -
                udigit * static_cast<uint64_t>(*key->entries[i][lev]->b)
            );
        }
    }
}


// NOT gate
void NOTbit64(LweSample64* out,
                       const LweSample64* in_bit,
                       const BootstrapContext64* ctx,
                       const Globals* env) {
    assert(out != nullptr);
    assert(in_bit != nullptr);
    (void)ctx;

    ClearLweSample64(out, env);
    SubLweSample64(out, in_bit, env);
}

// AND-PBS
void PBS_ANDbit64(LweSample64* out,
                       const LweSample64* lhs_bit,
                       const LweSample64* rhs_bit,
                       const BootstrapContext64* ctx,
                       const Globals* env) {
    assert(ctx != nullptr);
    assert(ctx->gate_bsk != nullptr);
    assert(ctx->output_ksk != nullptr);

    const Globals* gate_env = ctx->gate_env ? ctx->gate_env : env;
    LweSample64* gate_input = new LweSample64(gate_env->N);
    BuildGateInput(gate_input, lhs_bit, rhs_bit, gate_env);

    LweSample64* shifted_input = new LweSample64(gate_env->N);
    TLweSample64* accumulator = new TLweSample64(gate_env->N);
    const Torus64 gate_offset = EncodeControlBitToTorus64(1);

    CopyLweSample64(shifted_input, gate_input, gate_env);
    AddConst_LWEmsg(shifted_input, gate_offset);
    //build LUT 
    BuildGateAccumulator(accumulator, ControlBinaryGate64::And, gate_offset, gate_env);

    RunPBS(out, accumulator, shifted_input, ctx, gate_env);

    delete shifted_input;
    delete accumulator;
    delete gate_input;
}

// ============================ Circuit bootstrapping ============================

namespace {
inline Torus64 CB_KSGadgetValue(int lev, int basebit) {
    // Return the gadget value 2^(-basebit*(lev+1)) in 64-bit torus representation
    return static_cast<Torus64>(UINT64_C(1) << (64 - (lev + 1) * basebit));
}

}  

// ==================== circuit bootstrapping (LWE -> RGSW) ====================
// Compute s^2 in Z[X]/(X^N + 1).
void NegacyclicSquare(IntPolynomiala* out, //this used to generate the scheme-switching key
                      const IntPolynomiala* s, int N){
    for (int j=0;j<N;++j) out->coefs[j]=0;
    for (int i=0;i<N;++i){ const long si=(long)s->coefs[i]; if(!si) continue;
        for (int j=0;j<N;++j){ const long v=si*(long)s->coefs[j]; const int kidx=i+j;
            if (kidx<N) 
                out->coefs[kidx]+=(int)v; 
            else
                //X^N = -1 
                out->coefs[kidx-N]-=(int)v; } }
}

// out = <decomp(a_poly), SSK> = RLWE(a_poly * s^2),  for SSK = GLev.Enc(sk = s, msg = s^2)
void GadgetProd_RLWE_ss(TLweSample64* out, const Torus64Polynomial* a_poly,
                                 const GLevCipher64* SSK, int l_ss, const Globals* env){
    const int N=env->N,k=env->k;
    IntPolynomiala* decomp=new_array1<IntPolynomiala>(l_ss,N);
    tGswTorus64PolynomialDecompH_explicit_l(decomp, a_poly, N, l_ss, env->bgbit);
    for (int q=0;q<=k;++q) 
        for (int j=0;j<N;++j) 
            out->a[q].coefs[j]=0;
    for (int p=0;p<l_ss;++p) 
        for (int q=0;q<=k;++q)
            torus64PolynomialMultAddKaratsuba_lvl2(&out->a[q], &decomp[p], &SSK->cts[p]->a[q], env);
    delete_array1<IntPolynomiala>(decomp);
}

// Scheme switch: build the RGSW row from in_row1
// intput : in_row1  = (a, b),  phase b-a*s = m*g
// output : out_row0 =       ,  phase       = -m*g*s
// Using the SSK, which encrypts s^2, first obtain an RLWE encryption (A, B) of a*s^2. 
// Add b to the mask component:
//     B - (A + b)*s
//       = a*s^2 - b*s
//       = -s*(b - a*s)
//       = -m*g*s.
void SchemeSwitchRow_ss(TLweSample64* out_row0, 
                        const TLweSample64* in_row1,
                        const GLevCipher64* SSK, int l_ss, const Globals* env){
    const int N=env->N;
    GadgetProd_RLWE_ss(out_row0, &in_row1->a[0], SSK, l_ss, env);
    for (int j=0;j<N;++j)
        out_row0->a[0].coefs[j]=(Torus64)((uint64_t)out_row0->a[0].coefs[j]+(uint64_t)in_row1->b->coefs[j]);
}

// LWE(bit) -> RGSW(bit), multi-value
void CircuitBTS64_ss_mvfbs(TGswSample64* out, const LweSample64* lwe_in, int in_scale,
                                        const ControlGateBootstrappingKey64* bsk, AutoKsKeyFFTa** aks_fft,
                                        const GLevCipher64* SSK, int l_ss,
                                        const Globals* env_pbs, const Globals* env_out){
    const int N=env_out->N, bg=env_out->bgbit, k=env_out->k, l=env_out->l;
    int p=1; while (p<l) p<<=1;
    Torus64Polynomial* tv=new Torus64Polynomial(N); TLweSample64* testvec=new TLweSample64(N);
    GlweCipher64* rotated=new GlweCipher64(N); GlweCipher64* shifted=new GlweCipher64(N); GlweCipher64* base=new GlweCipher64(N);
    LweSample64* scaled=new LweSample64(N);

    for (int j=0;j<N;++j) scaled->a[j]=(Torus64)((uint64_t)lwe_in->a[j]*(uint64_t)in_scale);
    *scaled->b=(Torus64)((uint64_t)*lwe_in->b*(uint64_t)in_scale);

    // test vector: coefficient j holds half_(j mod p) = 1/(2*Bg^(j mod p + 1)), unused slots are 0
    for (int j=0;j<N;++j){
        const int lev=j%p;
        tv->coefs[j]= lev<l ? (Torus64)(UINT64_C(1)<<(63-(lev+1)*bg)) : 0;
    }
    tLwe64NoiselessTrivial(testvec, tv, env_pbs);
    //BR, phase rounded to multiples of p
    BlindRotate64(rotated, testvec, scaled, bsk, env_pbs, p);
    for (int i=0;i<l;++i){
        const Torus64 half_i=(Torus64)(UINT64_C(1)<<(63-(i+1)*bg));
        // coefficient i -> 0, then shift {-half_i,+half_i} to {0, 1/Bg^(i+1)}
        MulTLWEbyMonomial(shifted, rotated, -i, env_out);
        shifted->b->coefs[0]=(Torus64)((uint64_t)shifted->b->coefs[0]+(uint64_t)half_i);
        // HomTrace -> RLWE(mu*g_out[i])
        RevHomTrace_Alg5_fft(base, shifted, aks_fft, 1, env_out);
        for (int q=0;q<=k;++q) for (int j=0;j<N;++j)
            out->samples[1][i].a[q].coefs[j]=base->a[q].coefs[j];
        SchemeSwitchRow_ss(&out->samples[0][i], &out->samples[1][i], SSK, l_ss, env_out);
    }
    delete tv; delete testvec; delete rotated; delete shifted; delete base; delete scaled;
}

// Generate an output KSK with an explicit decomposition base
// Each entry encrypts a gadget-scaled input-secret coefficient under the reduced control key supported on the first target_dim coordinates.
void Gen_OutputKSK_Basebit64( ControlOutputKeySwitchKey64* key,
                              const std::vector<MPPartyKey64>& parties,
                            int input_dim, int target_dim, int basebit, int levels,
                            double stdev, const Globals* env
) {
    assert(key != nullptr);
    assert(input_dim >= 0 && input_dim <= env->N);
    assert(target_dim >= 0 && target_dim <= env->N);
    assert(basebit > 0 && levels > 0 && levels * basebit <= 64);

    Del_OutputKSK64(key);

    key->input_dim  = input_dim;
    key->target_dim = target_dim;
    key->levels     = levels;
    key->basebit    = basebit;
    key->identity   = false;
    key->entries.assign(input_dim, std::vector<LweSample64*>(key->levels, nullptr));

    IntPolynomiala* control_secret = new IntPolynomiala(env->N);
    Recon_ControlSecret(control_secret, parties, env);

    for (int i = 0; i < input_dim; ++i) {
        for (int lev = 0; lev < key->levels; ++lev) {
            key->entries[i][lev] = new LweSample64(env->N);
            const Torus64 scale = static_cast<Torus64>(
                UINT64_C(1) << (64 - (lev + 1) * key->basebit)
            );
            const Torus64 message = control_secret->coefs[i] ? scale : 0;
            ENC_ControlTorus_smalldim(key->entries[i][lev], message, control_secret, target_dim, stdev, env);
        }
    }
    delete control_secret;
}
