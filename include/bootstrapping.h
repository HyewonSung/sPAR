#pragma once

#include <vector>

#include "poc_64types.h"
#include "mphe.h"

struct ControlGateBootstrappingKey64 {
    int input_dim;
    std::vector<TGswSample64*> secret_bits;
    std::vector<TGswSampleFFTa*> secret_bits_fft;

    ControlGateBootstrappingKey64()
        : input_dim(0) {}
};

struct ControlOutputKeySwitchKey64 {
    int input_dim;
    int target_dim;
    int levels;
    int basebit;
    bool identity;
    std::vector<std::vector<LweSample64*>> entries;

    ControlOutputKeySwitchKey64()
        : input_dim(0), target_dim(0), levels(0), basebit(0), identity(true) {}
};

struct BootstrapContext64 {
    const ControlGateBootstrappingKey64* gate_bsk;
    const ControlOutputKeySwitchKey64* output_ksk;
    const Globals* gate_env;

    BootstrapContext64()
        : gate_bsk(nullptr),
          output_ksk(nullptr),
          gate_env(nullptr) {}
};

void Config_BTSContext64(BootstrapContext64* ctx);

void Attach_GatePBSKeysToBTSContext64(BootstrapContext64* ctx,
                                           const ControlGateBootstrappingKey64* gate_bsk,
                                           const ControlOutputKeySwitchKey64* output_ksk);

void Attach_BTSEnvToContext64(BootstrapContext64* ctx,
                                   const Globals* gate_env);

void ENC_ControlBitConstant64(LweSample64* out,
                                       int bit,
                                       const BootstrapContext64* ctx,
                                       const Globals* env);

void MakeDataZero64(TLweSample64* out,
                    const Globals* env);

void Gen_GateBTSKey64(ControlGateBootstrappingKey64* key,
                                           const std::vector<MPPartyKey64>& parties,
                                           int input_dim,
                                           double stdev,
                                           const Globals* env);

void Del_GateBTSKey64(ControlGateBootstrappingKey64* key);

void Gen_OutputKSK_TargetDim64(
    ControlOutputKeySwitchKey64* key,
    const std::vector<MPPartyKey64>& parties,
    int input_dim,
    int target_dim,
    double stdev,
    const Globals* env
);

void Del_OutputKSK64(ControlOutputKeySwitchKey64* key);

void BlindRotate64(TLweSample64* out_accum,
                                     const TLweSample64* testvec,
                                     const LweSample64* gate_input,
                                     const ControlGateBootstrappingKey64* key,
                                     const Globals* env);

void SampleExtract_AtZero64(LweSample64* out,
                                  const TLweSample64* in_accum,
                                  const Globals* env);

void KS_controlbit64(LweSample64* out,
                                 const LweSample64* in_extracted,
                                 const ControlOutputKeySwitchKey64* key,
                                 const Globals* env);

void NOTbit64(LweSample64* out,
                       const LweSample64* in_bit,
                       const BootstrapContext64* ctx,
                       const Globals* env);

void PBS_ANDbit64(LweSample64* out,
                       const LweSample64* lhs_bit,
                       const LweSample64* rhs_bit,
                       const BootstrapContext64* ctx,
                       const Globals* env);



// ---- circuit bootstrapping (LWE -> RGSW) -----------------------------------
// Multi-value functional bootstrap -> homomorphic trace -> scheme switching.
void NegacyclicSquare(IntPolynomiala* out, const IntPolynomiala* s, int N);
void GadgetProd_RLWE_ss(TLweSample64* out, const Torus64Polynomial* a_poly,
                          const GLevCipher64* SSK, int l_ss, const Globals* env);
void SchemeSwitchRow_ss(TLweSample64* out_row0, const TLweSample64* in_row1,
                        const GLevCipher64* SSK, int l_ss, const Globals* env);
void CircuitBTS64_ss_mvfbs(TGswSample64* out, const LweSample64* lwe_in,
                                 const ControlGateBootstrappingKey64* bsk, AutoKsKeyFFTa** aks_fft,
                                 const GLevCipher64* SSK, int l_ss,
                                 const Globals* env_pbs, const Globals* env_out);

// Output key-switching key N -> n with an explicit gadget base.
void Gen_OutputKSK_Basebit64(
    ControlOutputKeySwitchKey64* ksk, const std::vector<MPPartyKey64>& parties,
    int input_dim, int target_dim, int basebit, double stdev, const Globals* env);
