#!/bin/bash
# Build sPAR for eta in {20,40,60,80,100}.
#   ./build.sh            build all five
#   ./build.sh 20         build only eta = 20
set -e

# EXTRA_FLAGS is appended verbatim, e.g. a verbose build:
#   EXTRA_FLAGS=-DPRINT_DEBUG=1 ./build.sh 20
CXXFLAGS="-std=c++2a -O3 -w -DUSE_FFT -march=native -Iinclude -I. -Ispqlios -fopenmp -no-pie ${EXTRA_FLAGS}"
mkdir -p build

# --- one-time: assemble/compile the SPQLIOS FFT objects -----------------------
if [ ! -f spqlios/fft_processor_spqlios.o ]; then
  echo "[1/2] building SPQLIOS FFT objects"
  g++ -c -O3 -march=native -Ispqlios -I. spqlios/spqlios-fft-impl.cpp    -o spqlios/spqlios-fft-impl.o
  g++ -c -O3 -march=native -Ispqlios -I. spqlios/fft_processor_spqlios.cpp -o spqlios/fft_processor_spqlios.o
  g++ -c spqlios/spqlios-fft-fma.s        -o spqlios/spqlios-fft-fma.o
  g++ -c spqlios/spqlios-ifft-fma.s       -o spqlios/spqlios-ifft-fma.o
  g++ -c spqlios/lagrangehalfc_impl_fma.s -o spqlios/lagrangehalfc_impl_fma.o
fi

OBJ="spqlios/spqlios-fft-fma.o spqlios/spqlios-ifft-fma.o spqlios/spqlios-fft-impl.o \
     spqlios/fft_processor_spqlios.o spqlios/lagrangehalfc_impl_fma.o"

ETAS="${@:-20 40 60 80 100}"
echo "[2/2] building sPAR for eta = $ETAS"
for E in $ETAS; do
  g++ $CXXFLAGS -DETA=$E -DPARTY_COUNT=$E \
      $OBJ src/global_random.cpp tests/test_spar_main.cpp -o build/spar_e$E
  echo "  -> build/spar_e$E"
done

echo
echo "Run with, e.g.:  OMP_NUM_THREADS=4 ./build/spar_e20"
