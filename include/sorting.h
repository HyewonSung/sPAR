#ifndef SPAR_SORTING_H
#define SPAR_SORTING_H

#include <vector>
#include "poc_64types.h"
#include "include/config.h"
#include "include/bootstrapping.h"


#define CMP_HALF_RANGE 8


void HomCompareSwap(TLweSample64*& key_a, TLweSample64*& val_a,
                         TLweSample64*& key_b, TLweSample64*& val_b,
                         const ControlGateBootstrappingKey64* cb_bsk,
                         AutoKsKeyFFTa** aks_fft, const GLevCipher64* SSK, int l_ss,
                         const ControlOutputKeySwitchKey64* ksk,
                         const Globals* env_pbs, const Globals* env);


void ObliviousSorting(std::vector<TLweSample64*>& K, std::vector<TLweSample64*>& V,
                        const ControlGateBootstrappingKey64* cb_bsk,
                        AutoKsKeyFFTa** aks_fft, const GLevCipher64* SSK, int l_ss,
                        const ControlOutputKeySwitchKey64* ksk,
                        const Globals* env_pbs, const Globals* env);

#endif
