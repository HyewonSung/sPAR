#pragma once

#include <cstdint>
#include <vector>

#include "poc_64types.h"

struct AutoKsKeyFFTa;
struct BootstrapContext64;
struct HomPlacingState64;

struct SortingInsertTimingStats {
    uint64_t candidates = 0;
    uint64_t and_zI_ns = 0;
    uint64_t not_has_written_ns = 0;
    uint64_t and_h_ns = 0;
    uint64_t affine_scale_ns = 0;
    uint64_t pseudo_embed_ns = 0;
    uint64_t rev_trace_ns = 0;
    uint64_t pt_mul_ns = 0;
    uint64_t add_payload_ns = 0;
    uint64_t update_I_ns = 0;
    uint64_t update_has_written_ns = 0;
};

void MP_RLWE_Encrypt64(TLweSample64* ct,
                               const TLweSample64* common_pk,
                               const Torus64Polynomial* torus_msg,
                               double stdev,
                               const Globals* env);

void ExtractRLWEtoLWE(LweSample64* out,
                                const TLweSample64* in,
                                int idx,
                                const Globals* env);


void ChoiceofDistinctThree(int* out, int num, int eta);
