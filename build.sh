#!/bin/bash
# Build sPAR for eta in {20,40,60,80,100}.
#   ./build.sh                  build all five (default parameters)
#   ./build.sh 20               build only eta = 20
#   PRESET=paper ./build.sh     build with the parameter set of the submitted paper
set -e

# PRESET selects a parameter set; every value can still be overridden with -D in EXTRA_FLAGS.
#   default : include/config.h as shipped              -> build/spar_e<E>
#   paper   : parameters behind the paper's Tables 4-5 -> build/spar_e<E>_paper
case "${PRESET:-default}" in
  default) PRESET_FLAGS=""; SUFFIX="" ;;
  paper)   PRESET_FLAGS="-DGATE_BG=9 -DGATE_L=3 -DKS_OUT_LEVELS=8 -DCB_PBS_BG=9 -DCB_PBS_L=5 \
                         -DCB_L_SS=6 -DREFRESH_HW=18 -DREFRESH_I=8"
           SUFFIX="_paper" ;;
  *) echo "unknown PRESET=$PRESET (expected: default or paper)"; exit 1 ;;
esac

# EXTRA_FLAGS is appended verbatim, e.g. to override one parameter:
#   EXTRA_FLAGS=-DREFRESH_I=16 ./build.sh 20
CXXFLAGS="-std=c++2a -O3 -w -DUSE_FFT -march=native -Iinclude -I. -Ispqlios -fopenmp -no-pie ${PRESET_FLAGS} ${EXTRA_FLAGS}"
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
echo "[2/2] building sPAR (PRESET=${PRESET:-default}) for eta = $ETAS"
for E in $ETAS; do
  g++ $CXXFLAGS -DETA=$E -DPARTY_COUNT=$E \
      $OBJ src/global_random.cpp tests/test_spar_main.cpp -o build/spar_e$E$SUFFIX
  echo "  -> build/spar_e$E$SUFFIX"
done

echo
echo "Run with, e.g.:  OMP_NUM_THREADS=4 ./build/spar_e${ETAS%% *}$SUFFIX"
