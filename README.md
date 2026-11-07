# sPAR: (Somewhat) Practical Anonymous Router 


This repository contains a TFHE-based proof-of-concept implementation of sPAR, an FHE-based anonymous router.
It can be used to measure the performance of sPAR.

Note that while preparing this artifact, we found a faster parameter set, which is now the default.
Running the code with the parameters of the paper still yields numbers close to the tables in the paper, and this README describes how to run both the new and the original parameter sets.

The implementation is based on a TFHE backend adapted to 64-bit arithmetic, and all related functions are included in this code.
Therefore, there are no external dependencies beyond a C++ compiler with OpenMP.


## 1. Requirements

| | |
|---|---|
| Compiler | `g++` with C++20 support (tested with 9.4.0, glibc 2.31) |
| CPU | x86-64 with **AVX2 + FMA** (the SPQLIOS FFT uses FMA assembly) |
| OpenMP | required (`-fopenmp`) |
| RAM | ~2 GB (peak RSS is about 1.3–1.4 GB, nearly independent of `eta`) |
| OS | Linux (tested on Ubuntu 20.04.6 LTS) |

Everything is compiled from the sources in this repository, so no installation step is needed.

## 2. Directory layout
This artifact is organized as follows.

```
build.sh                  build script (all eta, or selected ones)
generic_utils.h           RNG helpers and basic scalar types
poc_64types.h             ciphertext / polynomial types 

include/
  config.h                all tunable parameters 
  mphe.h                  multi-party HE: key generation, partial decryption
  bootstrapping.h         gate bootstrapping, key switching, circuit bootstrapping
  sorting.h               oblivious sorting over homomorphic-encrypted values
  hom_placing.h           bucket state (L, I), Algorithm 1 (HomPlacing), per-phase timing
  spar_helper.h           client-side encryption, coefficient extraction
  test_common64.h         Globals (scheme parameters) and test utilities

src/
  64header.h              TFHE core functions: FFT, external product, CMux, automorphisms, homomorphic trace, gadget decomposition
  mphe.cpp                mpc-based key generation, Partial.Dec / Final.Dec
  bootstrapping.cpp       blind rotation, Boolean gates, circuit bootstrapping
  sorting.cpp             oblivious sorting under encryption
  hom_placing.cpp         Algorithm 1 (homomorphic placing)
  spar_helper.cpp         client-side encryption and server-side extraction
  global_random.cpp       global RNG instance

spqlios/                  SPQLIOS negacyclic FFT 

tests/
  test_spar_main.cpp      main experiment file, it runs the protocol one round and verifies the result
```



## 3. Quick start

```bash
./build.sh 20                      # build for eta = 20 (~5 seconds)
OMP_NUM_THREADS=4 ./build/spar_e20 # run one protocol round (~3 minutes)
```

A successful run ends with:

```
======== PASS ========
```


## 4. Building
The user can specify the number of clients, denoted by `eta`, for each experiment. 
If no value is provided, the build script compiles the program for all default settings used in our experiments: `eta = 20, 40, 60, 80, 100`.

For example, the provided script can be used as follows:

```bash
./build.sh                 # build all five: eta = 20, 40, 60, 80, 100
./build.sh 20 100          # build only eta = 20 and eta = 100
PRESET=paper ./build.sh    # parameter set of the submitted paper -> build/spar_e<E>_paper
```


On the reference machine (4 cores) this will take :

| | time |
|---|---:|
| one-time SPQLIOS FFT objects | ~3 s |
| each `eta` | ~3 s |
| `./build.sh` from a clean checkout (all five) | ~20 s |


## 5. Running
After the build completes, an executable file is generated for each selected value of `eta`.

Running an executable performs one complete round of sPAR for the corresponding number of clients.
The number of OpenMP threads can be specified using `OMP_NUM_THREADS`. 
Since the measurements reported in the paper were obtained using 4 threads, use `OMP_NUM_THREADS=4` when reproducing the results in the paper.

For example:
```bash
OMP_NUM_THREADS=4 ./build/spar_e20 [seed]
```

The optional [seed] argument overrides the default RNG seed (0) and can be used to test correctness across independent protocol rounds.

For example, the following command runs the experiment with five different seeds:
```bash
for s in 1 2 3 4 5; do OMP_NUM_THREADS=4 ./build/spar_e20 $s | tail -1; done
```

The same commands can be used for `eta = 40, 60, 80, 100` by replacing spar_e20 with the corresponding executable, such as spar_e40, spar_e60, spar_e80, or spar_e100.



Each run prints the parameter header, the progress of each protocol step, the noise budget, the `END-TO-END` block (Table 5) and the `---- summary ----` block (Table 4), and ends with `PASS` or `FAIL`.



## 6. Reproducing the paper's tables
All measurements reported in the paper were obtained using 4 OpenMP threads on the reference machine (Intel Core i9-12900K, 4 vCPUs under WSL2).
For the full hardware and software configuration of the experimental environment, please refer to the paper.

Depending on the machine configuration and the OpenMP threading environment, the measured runtimes may differ from the values reported in the paper by a few percent.

### 6.1 Run the sweep
The following commands run the experiments for `eta = 20, 40, 60, 80, 100` and print the results corresponding to both Table 4 and Table 5.

On our reference machine, running the full sweep takes approximately **2.3 hours** in total.

```bash
./build.sh                                   # all five binaries
for E in 20 40 60 80 100; do
  OMP_NUM_THREADS=4 ./build/spar_e$E
done
```

Every run must end with `======== PASS ========`.



### 6.2 Table 4 / Table 5
For each value of eta, the values used in Table 4 are printed in the final ---- summary ---- block of the output.


The values used in Table 5 (`Current backend`) are printed with the block beginning with ================ END-TO-END (eta=20, 20 messages) ================ .


Independently of timing, every run verifies the round and ends with the verdict:

```
======== PASS ========
```

`PASS` means that for every bucket, that the messages decrypted from it are exactly the ones the plaintext choice-of-three process would have placed there.


The `[noise budget]` line reports the worst payload noise observed in the round against the decryption boundary `Delta/2 = 2^59`, and the margin, i.e. how many bits of headroom are left. 



## 7. Parameters

All parameters are defined in **`include/config.h`**.
Invalid combinations (e.g. `ETA != PARTY_COUNT`, or a gadget with `B^l > 2^64`) are rejected at compile time by the `static_assert` checks at the end of `config.h`.

The default parameter set differs from the one in the paper.
We found parameters that run faster and made them the default; the original parameter set of the paper can still be built with `PRESET=paper`.

The table below shows how the parameters changed.

| Parameter | Macro | default | `PRESET=paper` |
|---|---|---:|---:|
| LWE dimension `n` | `LWE_N` | 1024 | 1024 |
| encryption / evaluation-key / partial-decryption noise `sigma` | `SIGMA_LOG` | 2^-54 | 2^-54 |
| key-switching noise `sigma_ksk` | `SIGMA_KSK_LOG` | 2^-26 | 2^-26 |
| gate bootstrapping gadget `(B, l)` | `GATE_BG`, `GATE_L` | (2^24, 1) | (2^9, 3) |
| `N -> n` key switching `(B, l)` | `KS_OUT_BASEBIT`, `KS_OUT_LEVELS` | (2^8, 2) | (2^8, 8) |
| circuit-bootstrapping blind rotation `(B, l)` | `CB_PBS_BG`, `CB_PBS_L` | (2^10, 4) | (2^9, 5) |
| circuit-bootstrapping output RGSW `(B, l)` | `OUT_BG`, `OUT_L` | (2^9, 3) | (2^9, 3) |
| homomorphic-trace key levels | `CB_AKS_L` | 6 | 6 |
| scheme-switching key levels | `CB_L_SS` | 5 | 6 |
| refresh period of `hasWritten` (candidates) | `REFRESH_HW` | 72 | 18 |
| refresh period of `I` (messages) | `REFRESH_I` | 32 | 8 |

The plaintext modulus `t = 2^4` and the ciphertext modulus `q = 2^64` are the same in both sets.

## 7.1 Updated Performance Results

The following results were measured with the default parameters on the reference machine (4 OpenMP threads, one run per `eta`, all runs `PASS`).

**Computation time of each step** :

| `eta` | Client:Send<br>per message | Server:Write<br>per message | Partial Dec.<br>`3*eta` ciphertexts | Sorting<br>3 slots |
|---:|---:|---:|---:|---:|
| 20  | 0.002 s | 7.65 s  | 0.014 s | 0.183 s |
| 40  | 0.002 s | 14.95 s | 0.028 s | 0.175 s |
| 60  | 0.002 s | 22.11 s | 0.039 s | 0.168 s |
| 80  | 0.002 s | 30.44 s | 0.055 s | 0.164 s |
| 100 | 0.002 s | 36.32 s | 0.065 s | 0.172 s |

**End-to-end latency per round**, excluding the one-time setup:

| `eta` | end-to-end |
|---:|---:|
| 20  | 155.0 s  |
| 40  | 606.1 s  |
| 60  | 1340.9 s |
| 80  | 2464.3 s |
| 100 | 3675.6 s |


## 8. License

This artifact is released under the **Apache License 2.0**; see `LICENSE`.

We chose Apache 2.0 because the FHE backend is derived from the TFHE library (version 1.0), which is itself Apache 2.0 licensed. 
`LICENSE` also records what comes from TFHE and what was changed:

* `spqlios/` — the SPQLIOS negacyclic FFT from TFHE v1.0. Four files are
  unmodified; three gained 64-bit torus entry points alongside the original
  32-bit ones.
* `src/64header.h` — derived from TFHE v1.0, re-expressed for 64-bit torus
  elements and reduced to the routines sPAR needs.

Everything else is original work of the sPAR authors.
