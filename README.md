# sPAR: (Somewhat) Practical Anonymous Router — Artifact

Proof-of-concept implementation of **sPAR** : an anonymous router built on fully homomorphic encryption (FHE) over a single untrusted server. 

This artifact reproduces the performance numbers reported in the paper (Tables 4 and 5).

The implementation is based on a TFHE backend adapted to 64-bit arithmetic, and all related functions are included in this artifact.
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

spqlios/                  SPQLIOS negacyclic FFT (sources and FMA assembly)

tests/
  test_spar_main.cpp      main experiment file, it runs the protocol one round and verifies the result
```



## 3. Quick start

```bash
./build.sh 20                      # build for eta = 20 (~5 seconds)
OMP_NUM_THREADS=4 ./build/spar_e20 # run one protocol round (~4 minutes)
```

A successful run ends with:

```
======== PASS ========
```


## 4. Building
The main source code is located in `tests/test_spar_main.cpp`.

Before running the experiments, the project must first be built.

The user can specify the number of clients, denoted by `eta`, for each experiment. 
The program is rebuilt according to the specified value of `eta`. 
If no value is provided, the build script compiles the program for all default settings used in our experiments: `eta = 20, 40, 60, 80, 100`.

For example, the provided script can be used as follows:

```bash
./build.sh                 # build all five: eta = 20, 40, 60, 80, 100
./build.sh 20 100          # build only eta = 20 and eta = 100
```


On the reference machine (4 cores) this will take :

| | time |
|---|---:|
| one-time SPQLIOS FFT objects | ~3 s |
| each `eta` | ~3 s |
| `./build.sh` from a clean checkout (all five) | ~20 s |
| `./build.sh 20` with the FFT objects already built | ~3 s |

The equivalent manual command is:

```bash
mkdir -p build

# one-time: SPQLIOS FFT objects
g++ -c -O3 -march=native -Ispqlios -I. spqlios/spqlios-fft-impl.cpp      -o spqlios/spqlios-fft-impl.o
g++ -c -O3 -march=native -Ispqlios -I. spqlios/fft_processor_spqlios.cpp -o spqlios/fft_processor_spqlios.o
g++ -c spqlios/spqlios-fft-fma.s        -o spqlios/spqlios-fft-fma.o
g++ -c spqlios/spqlios-ifft-fma.s       -o spqlios/spqlios-ifft-fma.o
g++ -c spqlios/lagrangehalfc_impl_fma.s -o spqlios/lagrangehalfc_impl_fma.o

# the protocol, for eta = 20
g++ -std=c++2a -O3 -w -DUSE_FFT -march=native \
    -DETA=20 -DPARTY_COUNT=20 \
    -Iinclude -I. -Ispqlios -fopenmp -no-pie \
    spqlios/spqlios-fft-fma.o spqlios/spqlios-ifft-fma.o \
    spqlios/spqlios-fft-impl.o spqlios/fft_processor_spqlios.o \
    spqlios/lagrangehalfc_impl_fma.o \
    src/global_random.cpp tests/test_spar_main.cpp \
    -o build/spar_e20
```

Note that the protocol requires the number of clients to be equal to the number of buckets. 
Therefore, `ETA` (the number of buckets) and `PARTY_COUNT` (the number of clients, with one message per client) must always be set to the same value.


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

Each execution performs one complete protocol round: setup, client encryption, homomorphic placing, oblivious sorting, partial/final decryption, and verification.

### Print debug details
The amount of runtime output is controlled by the `print_debug` flag in `tests/test_spar_main.cpp`, which is disabled by default.

The default output already includes all information required to reproduce the tables in the paper and to verify correctness. 
Enabling debug output additionally prints detailed per-message and per-bucket information.

The flag can be enabled at build time as follows:

```bash
./build.sh 20                                     # default (quiet)
EXTRA_FLAGS=-DPRINT_DEBUG=1 ./build.sh 20     # verbose build
```

The difference between running with and without debug output is shown below:
| | quiet (default) | verbose |
|---|---|---|
| parameter header | yes | yes |
| setup / client:send / sorting timings | yes | yes |
| `[server:write] done`, `[partial decryption] done` | yes | — |
| `[keygen] gate bsk ...`, `[Server:Write] starts`, `[sort] ...` phase headers | — | yes |
| `[msg N] choices = {...}` + per-message `p1a/p1b/p2/p3` breakdown | — | yes |
| per-bucket `L-slots ... MATCH ? OK ... noise~...` | — | yes |
| noise budget | yes | yes |
| `END-TO-END` block (Table 5) | yes | yes |
| `---- summary ----` block (Table 4) | yes | yes |
| `PASS` / `FAIL` | yes | yes |

The quiet log is a fixed ~50 lines whatever `eta` is, because the per-message and
per-bucket lines are exactly what it drops. The verbose log grows by 8 lines per
message, so `eta = 20` comes out at about 216 lines.


## 6. Reproducing the paper's tables
All measurements reported in the paper were obtained using **4 OpenMP threads** on the reference machine (Intel Core i9-12900K, 4 vCPUs under WSL2).
For the full hardware and software configuration of the experimental environment, please refer to the paper.

Both Table 4 and Table 5 are obtained from the same five runs.

Depending on the machine configuration and the OpenMP threading environment, the measured runtimes may differ from the values reported in the paper by a few percent.

### 6.1 Run the sweep
The following commands run the experiments for `eta = 20, 40, 60, 80, 100` and generate the results corresponding to both Table 4 and Table 5.
The output for each value of `eta` is stored in `logs/run_e<E>.log`.

On our reference machine, running the full sweep takes approximately **3.7 hours** in total.

```bash
./build.sh                                   # all five binaries
mkdir -p logs
for E in 20 40 60 80 100; do
  OMP_NUM_THREADS=4 ./build/spar_e$E 2>&1 | tee logs/run_e$E.log
done
```

Every log must end with `======== PASS ========`.



### 6.2 Table 4 
For each value of eta, the values used in Table 4 are printed in the final ---- summary ---- block of the corresponding log.

The following command extracts these blocks from all five runs:

```bash
for E in 20 40 60 80 100; do
  echo "== eta=$E"
  sed -n '/---- summary/,/^$/p' logs/run_e$E.log
done
```


### 6.3 Table 5 — end-to-end latency
In each log, the block beginning with ================ END-TO-END (eta=20, 20 messages) ================ 
reports the end-to-end runtime of one complete protocol round. 
The value reported in the paper corresponds to the end-to-end runtime excluding setup, as shown in (Table 5, `Current backend`).


Only the `Current backend` column is directly measured. 
The **`Zama-projected`** column is estimated by replacing the measured bootstrap costs with Zama's reported GPU latency (0.796 ms per gate bootstrap), while keeping the remaining costs unchanged. 
Circuit-bootstrap costs are scaled by the same factor.


### 6.4 Correctness

Independently of timing, every run verifies the round and ends with the verdict:

```
======== PASS ========
```

`PASS` means that for every bucket, that the messages decrypted from it are exactly the ones the plaintext choice-of-three process would have placed there, and that every payload coefficient decrypts cleanly.

The `noise budget` block reports the worst payload noise observed in the round against the decryption boundary `Delta/2 = 2^59`. 
The margin is how many bits of headroom are left.

A **verbose** build additionally shows the per-bucket evidence:

```
  bucket 0  L-slots(after sort): [0000] [514150] [78910] | answer: [0000] [514150] [78910] | MATCH ? OK | noise~2^51 (Delta/2=2^59)
```

* `MATCH ? OK` — the messages decrypted from the bucket equals the plaintext process would have placed there (`FAIL` otherwise).
* `noise~2^51 (Delta/2=2^59)` — the largest observed noise magnitude in that bucket against the decryption boundary and the gap is the correctness margin.



## 7. Parameters

All parameters are defined in **`include/config.h`** and can be overridden at compile time using `-D`, without modifying the source code.

The default values are those used for the experiments in the paper. 
They were selected so that the noise budget remains sufficient for a complete round with `eta = 100`, while the underlying LWE/RLWE instances meet the target security level.

The noise parameters, gadget decompositions, and refresh periods are closely related. 
Changing one of them may affect correctness or security, so the parameters below are grouped by how safely they can be modified.

### Fixed

The ring dimension `N = 2048` is **not configurable** in the current implementation.

The bundled SPQLIOS FFT backend is specialized for `N = 2048`, and several parts of the implementation assume this value. 
Changing `N` therefore requires modifying the FFT backend as well.

### Safe to change

| Parameter | Constraint |
|---|---|
| `ETA`, `PARTY_COUNT` | **Must be set to the same value** and must be at least 3, since each client chooses 3 distinct candidate buckets. Larger values increase the running time. |
| `SEED` | Can be set to any value. It can also be overridden at run time using `argv[1]`. |
| `PRINT_DEBUG` | `0` or `1`; controls output verbosity only. |

### Change with caution

The following parameters affect the noise growth or security of the scheme. 
They can be changed to explore other parameter settings, but the guarantees and measurements reported in the paper no longer apply.

| Parameter | Effect |
|---|---|
| `REFRESH_HW`, `REFRESH_I` | Increasing them too much may allow the control-bit noise to cross its decision boundary and cause incorrect placement. Smaller values refresh more often and are therefore slower. |
| `LWE_N` | Dimension of the reduced LWE instance. Lowering it reduces the estimated security level and invalidates the 117-bit security setting used in the paper. |
| `SIGMA_KSK_LOG` | Noise parameter of the `N -> n` key-switching key. Lower noise may weaken security, while higher noise reduces the available correctness margin. |
| `CB_PBS_L`, `CB_AKS_L`, `CB_L_SS` | Gadget decomposition levels used in circuit bootstrapping. Changing them affects the circuit-bootstrap noise and running time. |
| `KS_OUT_BASEBIT` | Must **divide 64**. The implementation checks `64 % basebit == 0` and aborts otherwise. |

After changing any of these parameters, verify that the execution still ends with

```
======== PASS ========
```




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
