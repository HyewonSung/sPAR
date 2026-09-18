#pragma once

#include <chrono>
#include <vector>

#include "poc_64types.h"

struct BootstrapContext64;

struct HomPlacingState64 {
    int eta;
    int slots;
    int choice_bitlen;

    std::vector<std::vector<TLweSample64*>> L;
    std::vector<std::vector<LweSample64*>> I;
};

HomPlacingState64* CreateHomPlacingState64(int eta,
                                           int slots,
                                           const BootstrapContext64* boot,
                                           const Globals* env);

void DeleteHomPlacingState64(HomPlacingState64* st);


inline double now_ms() {
    using C = std::chrono::steady_clock;
    return std::chrono::duration<double, std::milli>(C::now().time_since_epoch()).count();
}

// ---- per-phase timing of one HomPlacing call -------------------------------
struct StageTimes {
    double p1a_and_ms   = 0;  // Phase 1a  AND(z, I)                  -- parallel
    double p1b_and_ms   = 0;  // Phase 1b  AND(zI, NOT hasWritten)    -- serial
    double p2_ms        = 0;  // Phase 2   circuit bootstrap + ext. product -- parallel
    double p2_cb_cpu_ms = 0;  //   of which circuit bootstrapping          -- summed over threads
    double p2_extprod_cpu_ms = 0;  //   of which FFT convert + external product -- summed over threads
    double p3_add_ms    = 0;  // Phase 3   accumulate into L          -- serial
    double total_ms     = 0;
    void add(const StageTimes& o) {
        p1a_and_ms+=o.p1a_and_ms; p1b_and_ms+=o.p1b_and_ms; p2_ms+=o.p2_ms;
        p2_cb_cpu_ms+=o.p2_cb_cpu_ms; p2_extprod_cpu_ms+=o.p2_extprod_cpu_ms;
        p3_add_ms+=o.p3_add_ms; total_ms+=o.total_ms;
    }
};


void HomPlacingPhased(HomPlacingState64* st,
                      const TLweSample64* payload_ct,
                      const std::vector<std::vector<const LweSample64*>>& z_onehot,
                      const BootstrapContext64* boot,
                      const ControlGateBootstrappingKey64* cb_bsk,
                      AutoKsKeyFFTa** aks_fft, const GLevCipher64* SSK, int l_ss,
                      const Globals* env_pbs, const Globals* env,
                      StageTimes* T);
