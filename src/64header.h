#pragma once

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <sys/time.h>
#include "generic_utils.h"
#include "spqlios/lagrangehalfc_impl.h"
#include "poc_64types.h"
#include <chrono>
#include <random>
#include <vector>
#include <immintrin.h>

extern Random* global_random;





#if defined(USE_FFT)
#ifndef HAVE_GET_FFTP2048






static inline FFT_Processor_Spqlios& get_fftp2048() {
    return fftp2048;
}

static inline void fftp2048_execute_reverse_int(double* out, const int* in) {
    fftp2048.execute_reverse_int(out, in);
}

static inline void fftp2048_execute_direct_torus64(Torus64* out, const double* in) {
    fftp2048.execute_direct_torus64(out, in);
}

static inline void fftp2048_execute_reverse_torus64(double* out, const Torus64* in) {
    fftp2048.execute_reverse_torus64(out, in);
}
#endif
#endif

inline double sample_gaussian(double sigma) {
    static thread_local std::mt19937_64 gen(std::random_device{}());
    static thread_local std::normal_distribution<double> dist(0.0, 1.0);
    return dist(gen) * sigma; 
}


void torus64PolynomialMultNaive_plain_aux(Torus64* __restrict result, const int* __restrict poly1, const Torus64* __restrict poly2, const int N) {
    const int _2Nm1 = 2*N-1;
    Torus64 ri;
    for (int i=0; i<N; i++) {
        ri=0;
        for (int j=0; j<=i; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
    for (int i=N; i<_2Nm1; i++) {
        ri=0;
        for (int j=i-N+1; j<N; j++) ri += poly1[j]*poly2[i-j];
        result[i]=ri;
    }
}




void Karatsuba64_aux(Torus64* R, const int* A, const Torus64* B, const int size, const char* buf){
    const int h = size / 2;
    const int sm1 = size-1;

    
    
    if (h<=4) {
        torus64PolynomialMultNaive_plain_aux(R, A, B, size);
        return;
    }

    
    int* Atemp = (int*) buf; buf += h*sizeof(int);
    Torus64* Btemp = (Torus64*) buf; buf+= h*sizeof(Torus64);
    Torus64* Rtemp = (Torus64*) buf; buf+= size*sizeof(Torus64);
    

    for (int i = 0; i < h; ++i) Atemp[i] = A[i] + A[h+i];
    for (int i = 0; i < h; ++i) Btemp[i] = B[i] + B[h+i];

    
    Karatsuba64_aux(R, A, B, h, buf); 
    Karatsuba64_aux(R+size, A+h, B+h, h, buf); 
    Karatsuba64_aux(Rtemp, Atemp, Btemp, h, buf);
    R[sm1]=0; 
    for (int i = 0; i < sm1; ++i) Rtemp[i] -= R[i] + R[size+i];
    for (int i = 0; i < sm1; ++i) R[h+i] += Rtemp[i];
}





void torus64PolynomialMultKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomiala* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; 

    
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    
    for (int i = 0; i < N2-1; ++i) result->coefs[i] = R[i] - R[N2+i];
    result->coefs[N2-1] = R[N2-1];

    delete[] R;
    delete[] buf;
}

static inline void MultiplyTlweByIntPolynomialLocal(TLweSample64* out,
                                                    const TLweSample64* in,
                                                    const IntPolynomiala* plain,
                                                    const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        torus64PolynomialMultKaratsuba_lvl2(&out->a[q], plain, &in->a[q], env);
    }
    torus64PolynomialMultKaratsuba_lvl2(out->b, plain, in->b, env);
}

static inline bool CanonicalLessTlweLocal(const TLweSample64* lhs,
                                          const TLweSample64* rhs,
                                          const Globals* env) {
    for (int q = 0; q <= env->k; ++q) {
        for (int i = 0; i < env->N; ++i) {
            const uint64_t lv = static_cast<uint64_t>(lhs->a[q].coefs[i]);
            const uint64_t rv = static_cast<uint64_t>(rhs->a[q].coefs[i]);
            if (lv < rv) return true;
            if (lv > rv) return false;
        }
    }
    for (int i = 0; i < env->N; ++i) {
        const uint64_t lv = static_cast<uint64_t>(lhs->b->coefs[i]);
        const uint64_t rv = static_cast<uint64_t>(rhs->b->coefs[i]);
        if (lv < rv) return true;
        if (lv > rv) return false;
    }
    return false;
}





void torus64PolynomialMultAddKaratsuba_lvl2(Torus64Polynomial* result, const IntPolynomiala* poly1, const Torus64Polynomial* poly2, const Globals* env){
    const int N2 = env->N;
    Torus64* R = new Torus64[2*N2-1];
    char* buf = new char[32*N2]; 

    
    Karatsuba64_aux(R, poly1->coefs, poly2->coefs, N2, buf);

    
    for (int i = 0; i < N2-1; ++i) result->coefs[i] += R[i] - R[N2+i];
    result->coefs[N2-1] += R[N2-1];

    delete[] R;
    delete[] buf;
}


void tLwe64EncryptZero(TLweSample64* cipher, const double stdev, const Globals* env){
    const int N = env->N;
    const int k= env->k;
    for (int j = 0; j < N; ++j) cipher->b->coefs[j] = random_gaussian64(0, stdev);

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j) cipher->a[i].coefs[j] = random_int64();
    }
    for (int i = 0; i < k; ++i) torus64PolynomialMultAddKaratsuba_lvl2(cipher->b, &env->tlwekey[i], &cipher->a[i], env);
}

void tLwe64Encrypt(TLweSample64* cipher,const Torus64Polynomial* mess, const double stdev, const Globals* env){
    const int N = env->N;
    

   tLwe64EncryptZero(cipher, stdev, env);

    for (int32_t j = 0; j < N; ++j)
       cipher->b->coefs[j] += mess->coefs[j];
}

void tLwe64EncryptZero_debug(TLweSample64* cipher, const double stdev, const Globals* env){
    const int N = env->N;
    const int k= env->k;

    for (int j = 0; j < N; ++j)
        cipher->b->coefs[j] = random_gaussian64(0, stdev); 

    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < N; ++j)
            cipher->a[i].coefs[j] = random_int64();
    }

    
    for (int i = 0; i < k; ++i)
        torus64PolynomialMultAddKaratsuba_lvl2(cipher->b, &env->tlwekey[i], &cipher->a[i], env);

}

void tLwe64Encrypt_debug(TLweSample64* cipher, const Torus64Polynomial* mess, const double stdev, const Globals* env){
    const int N = env->N;


    tLwe64EncryptZero_debug(cipher, stdev, env);

    for (int32_t j = 0; j < N; ++j)
        cipher->b->coefs[j] += mess->coefs[j];

}



void tLwe64Phase_lvl2(Torus64Polynomial* phase, const TLweSample64* cipher, const Globals* env){
    const int N = env->N;
    const int k= env->k;

    
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -cipher->b->coefs[j];
    }



    for (int i = 0; i < k; ++i) {
        torus64PolynomialMultAddKaratsuba_lvl2(phase, &env->tlwekey[i], &cipher->a[i], env);
    }


    
    for (int j = 0; j < N; ++j) {
        phase->coefs[j] = -phase->coefs[j];
    }

}







/*
LagrangeHalfCPolynomiala* new_LagrangeHalfCPolynomiala_array(int nbelts, int N) {
    return new_array1<LagrangeHalfCPolynomiala>(nbelts,N);
}

void delete_LagrangeHalfCPolynomial_array(int nbelts, LagrangeHalfCPolynomiala* data) {
    delete_array1<LagrangeHalfCPolynomiala>(data);
}
*/

#ifdef USE_FFT
void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomiala* result, const IntPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048_execute_reverse_int(result->values, source->coefs);
}


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomiala* result, const Globals* env) {
    const int N = env->N;
#if defined(__AVX2__)
    const __m256d z = _mm256_setzero_pd();
    int i = 0;
    for (; i + 4 <= N; i += 4) {
        _mm256_storeu_pd(result->values + i, z);
    }
    for (; i < N; ++i) result->values[i] = 0.0;
#else
    for (int i = 0; i < N; ++i) result->values[i] = 0.0;
#endif
}
void LagrangeHalfCPolynomialAddTo_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const Globals* env) {
    const int N = env->N;
#if defined(__AVX2__)
    int i = 0;
    for (; i + 4 <= N; i += 4) {
        __m256d r = _mm256_loadu_pd(result->values + i);
        __m256d x = _mm256_loadu_pd(a->values + i);
        r = _mm256_add_pd(r, x);
        _mm256_storeu_pd(result->values + i, r);
    }
    for (; i < N; ++i) result->values[i] += a->values[i];
#else
    for (int i = 0; i < N; ++i) result->values[i] += a->values[i];
#endif
}
void LagrangeHalfCPolynomialSubTo_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const Globals* env) {
    const int N = env->N;
#if defined(__AVX2__)
    int i = 0;
    for (; i + 4 <= N; i += 4) {
        __m256d r = _mm256_loadu_pd(result->values + i);
        __m256d x = _mm256_loadu_pd(a->values + i);
        r = _mm256_sub_pd(r, x);
        _mm256_storeu_pd(result->values + i, r);
    }
    for (; i < N; ++i) result->values[i] -= a->values[i];
#else
    for (int i = 0; i < N; ++i) result->values[i] -= a->values[i];
#endif
}
void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const LagrangeHalfCPolynomiala* b, const Globals* env) {
    const int Ns2 = env->N/2;
    LagrangeHalfCPolynomialAddMulASM(result->values, a->values, b->values, Ns2);
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048_execute_direct_torus64(result->coefs, source->values);
}

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomiala* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N==2048);
    fftp2048_execute_reverse_torus64(result->values, source->coefs);
}


#else

void IntPolynomial_ifft_lvl2(LagrangeHalfCPolynomiala* result, const IntPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    result->setIntPoly(source, 2048);
}


void LagrangeHalfCPolynomialClear_lvl2(LagrangeHalfCPolynomiala* result, const Globals* env) {
    assert(env->N==2048);
    result->setZeroTorus64Poly(2048);
}

void LagrangeHalfCPolynomialAddMul_lvl2(LagrangeHalfCPolynomiala* result, const LagrangeHalfCPolynomiala* a, const LagrangeHalfCPolynomiala* b, const Globals* env) {
    assert(env->N==2048);
    assert(result->torus64Poly!=0);
assert(a->intPoly!=0);
    assert(b->torus64Poly!=0);
    torus64PolynomialMultAddKaratsuba_lvl2(result->torus64Poly, a->intPoly, b->torus64Poly, env);
}

void TorusPolynomial64_fft_lvl2(Torus64Polynomial* result, const LagrangeHalfCPolynomiala* source, const Globals* env) {
    assert(env->N==2048);
    assert(source->torus64Poly!=0);
    for (int i=0; i<2048; i++) result->coefs[i]=source->torus64Poly->coefs[i];
}

void TorusPolynomial64_ifft_lvl2(LagrangeHalfCPolynomiala* result, const Torus64Polynomial* source, const Globals* env) {
    assert(env->N==2048);
    result->setTorus64Poly(source, 2048);
}
#endif



void tLwe64NoiselessTrivial(TLweSample64* cipher, const Torus64Polynomial* mess, const Globals* env){
    const int N = env->N;
    const int k= env->k;

     for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            cipher->a[i].coefs[j] = 0;
        }
    }
 for (int j = 0; j < N; ++j) {

    cipher->b->coefs[j] = mess->coefs[j];
  }
}

void int_to_bin_digit(unsigned int in, int count, int64_t* out)
{
        unsigned int mask =1U << (count-1);
                int k;
                        for (k=0;k< count; k++){
                                                out[k]=(in & mask) ? 1 : 0;
                                                                       in <<=1;


                                                                                        }

}

void tGswTorus64PolynomialDecompH(IntPolynomiala* result, const Torus64Polynomial* sample, const Globals* env){
            const int N = env->N;
            const int l = env->l;
            const int Bgbit = env->bgbit;
            
            const uint64_t Bg = UINT64_C(1)<<Bgbit;
            const uint64_t mask = Bg-1;
            const int64_t halfBg = Bg/2;
            
            
            static thread_local std::vector<uint64_t> tl_decomp_buf;
            if (static_cast<int>(tl_decomp_buf.size()) < N) tl_decomp_buf.resize(N);
            uint64_t* buf = tl_decomp_buf.data();
            const uint64_t offset = env->torusDecompOffset;

    
    for (int j = 0; j < N; ++j) buf[j]=sample->coefs[j]+offset;

    
    for (int p = 0; p < l; ++p) {
        const int decal = (64-(p+1)*Bgbit);
         int* res_p = result[p].coefs; 
        for (int j = 0; j < N; ++j) {
            uint64_t temp1 = (buf[j] >> decal) & mask;
            res_p[j] = temp1 - halfBg;
        }
    }
}

        void tGsw64DecompH(IntPolynomiala* result, const TLweSample64* sample, const Globals* env){
    const int l = env->l;
    const int k=env->k;
    for (int i = 0; i <= k; ++i) tGswTorus64PolynomialDecompH(result+(i*l), &sample->a[i], env);
}


static void tGswTorus64PolynomialDecompH_explicit_l(
    IntPolynomiala* result, const Torus64Polynomial* sample,
    int N, int l, int Bgbit
) {
    const uint64_t Bg = UINT64_C(1) << Bgbit;
    const uint64_t mask = Bg - 1;
    const int64_t halfBg = (int64_t)(Bg / 2);
    uint64_t offset = 0;
    for (int p = 0; p < l; ++p)
        offset += (UINT64_C(1) << (64 - (p + 1) * Bgbit)) * (uint64_t)halfBg;
    static thread_local std::vector<uint64_t> tl_buf_l;
    if ((int)tl_buf_l.size() < N) tl_buf_l.resize(N);
    uint64_t* buf = tl_buf_l.data();
    for (int j = 0; j < N; ++j) buf[j] = (uint64_t)sample->coefs[j] + offset;
    for (int p = 0; p < l; ++p) {
        const int decal = 64 - (p + 1) * Bgbit;
        int* res_p = result[p].coefs;
        for (int j = 0; j < N; ++j) {
            uint64_t temp1 = (buf[j] >> decal) & mask;
            res_p[j] = (int)(temp1 - (uint64_t)halfBg);
        }
    }
}

        void tGswExternMulToTLwe1(TLweSample64 *accum, const TGswSample64 *sample, const Globals *env) {

    const int32_t N = env->N;
    const int32_t k= env->k;
    const int l = env->l;
    const int32_t kpl = (k+1)*l;
    

 IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
    tGsw64DecompH(decomp, accum, env);
    
 for (int i = 0; i <= k ; ++i) {
        for (int j = 0; j < N; ++j) {
            accum->a[i].coefs[j] = 0;
        }
    }

    for (int32_t i = 0; i < kpl; i++) {
      
      

     for (int j = 0; j <= k; ++j) torus64PolynomialMultAddKaratsuba_lvl2(accum->a+j, &decomp[i], &sample->allsamples[i].a[j], env);
    }

}


        void CMux(TLweSample64 *result, const TGswSample64 *eps, const TLweSample64 *c0, TLweSample64 *c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

         IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N); Torus64Polynomial* phase = new Torus64Polynomial(N);

         TLweSampleFFTa* accFFT = new TLweSampleFFTa(N);
         TGswSampleFFTa* epsFFT= new TGswSampleFFTa(l,N);

         for (int i=0;i<kpl;i++)
                for (int q=0;q<=k;q++)
                  TorusPolynomial64_ifft_lvl2(&epsFFT->allsamples[i].a[q],&eps->allsamples[i].a[q],  env);

 
         for (int q = 0; q <= k; ++q)
         for (int j = 0; j < N; ++j) c1->a[q].coefs[j] -= c0->a[q].coefs[j];

 

        tGsw64DecompH(decomp, c1, env);
        for (int p = 0; p < kpl; ++p) IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);
        
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

        
auto start = std::chrono::high_resolution_clock::now();
         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &epsFFT->allsamples[p].a[q], env);

auto end = std::chrono::high_resolution_clock::now();
std::chrono::duration<double,std::milli> execution_time = end-start;
/* debug disabled */

        
        for (int q = 0; q <= k; ++q) TorusPolynomial64_fft_lvl2(c1->a+q,accFFT->a+q, env);




        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) c1->a[q].coefs[j] += c0->a[q].coefs[j];
for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].coefs[j] = c1->a[q].coefs[j];

}





void ExternalProductFFT(TLweSample64 *result, const TGswSampleFFTa *eps, const TLweSample64 *in, const Globals* env){
    const int N = env->N;
    const int k = env->k;
    const int l = env->l;
    const int kpl = (k+1)*l;

    struct Cache {
        int N=0, k=0, l=0, kpl=0;
        IntPolynomiala* decomp=nullptr;
        LagrangeHalfCPolynomiala* decompFFT=nullptr;
        TLweSampleFFTa* accFFT=nullptr;
    };
    static thread_local Cache cache;

    auto ensure_cache = [&]() {
        if (cache.accFFT && cache.N==N && cache.k==k && cache.l==l) return;
        if (cache.decomp)    delete_array1<IntPolynomiala>(cache.decomp);
        if (cache.decompFFT) delete_array1<LagrangeHalfCPolynomiala>(cache.decompFFT);
        if (cache.accFFT)    delete cache.accFFT;
        cache.N=N; cache.k=k; cache.l=l; cache.kpl=kpl;
        cache.decomp    = new_array1<IntPolynomiala>(kpl, N);
        cache.decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl, N);
        cache.accFFT    = new TLweSampleFFTa(N);
    };
    ensure_cache();


    tGsw64DecompH(cache.decomp, in, env);
    for (int p = 0; p < kpl; ++p)
        IntPolynomial_ifft_lvl2(cache.decompFFT + p, cache.decomp + p, env);


    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialClear_lvl2(cache.accFFT->a + q, env);
    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= k; ++q)
            LagrangeHalfCPolynomialAddMul_lvl2(cache.accFFT->a + q, cache.decompFFT + p, &eps->allsamples[p].a[q], env);


    for (int q = 0; q <= k; ++q)
        TorusPolynomial64_fft_lvl2(result->a + q, cache.accFFT->a + q, env);
}

void CMuxFFT(TLweSample64 *result, const TGswSampleFFTa *eps, const TLweSample64 *c0, TLweSample64 *c1, const Globals* env){
    const int N = env->N;
    const int k = env->k;
    const int l = env->l;
    const int kpl = (k+1)*l;

    
    struct Cache {
        int N=0, k=0, l=0, kpl=0;
        IntPolynomiala* decomp=nullptr;
        LagrangeHalfCPolynomiala* decompFFT=nullptr;
        TLweSampleFFTa* accFFT=nullptr;
    };
    static thread_local Cache cache;

    auto ensure_cache = [&]() {
        if (cache.accFFT && cache.N==N && cache.k==k && cache.l==l) return;
        if (cache.decomp)    delete_array1<IntPolynomiala>(cache.decomp);
        if (cache.decompFFT) delete_array1<LagrangeHalfCPolynomiala>(cache.decompFFT);
        if (cache.accFFT)    delete cache.accFFT;
        cache.N=N; cache.k=k; cache.l=l; cache.kpl=kpl;
        cache.decomp    = new_array1<IntPolynomiala>(kpl, N);
        cache.decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl, N);
        cache.accFFT    = new TLweSampleFFTa(N);
    };
    ensure_cache();

    
    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            c1->a[q].coefs[j] -= c0->a[q].coefs[j];

    
    tGsw64DecompH(cache.decomp, c1, env);
    for (int p = 0; p < kpl; ++p)
        IntPolynomial_ifft_lvl2(cache.decompFFT + p, cache.decomp + p, env);

    
    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialClear_lvl2(cache.accFFT->a + q, env);

    
    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= k; ++q)
            LagrangeHalfCPolynomialAddMul_lvl2(cache.accFFT->a + q, cache.decompFFT + p, &eps->allsamples[p].a[q], env);

    
    for (int q = 0; q <= k; ++q)
        TorusPolynomial64_fft_lvl2(c1->a + q, cache.accFFT->a + q, env);

    
    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            c1->a[q].coefs[j] += c0->a[q].coefs[j];

    
    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            result->a[q].coefs[j] = c1->a[q].coefs[j];
}
void CMuxFFTdb(TLweSampleFFTa *result, const TGswSampleFFTa *eps, const Torus64 c0, const Torus64 c1, const Globals* env){

        const int l=env->l;
        const int N = env->N;
        const int k= env->k;
        const int kpl=(k+1)*l;

        IntPolynomiala* decomp = new_array1<IntPolynomiala>(kpl,N);
         LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl,N);
        TLweSampleFFTa* accFFT= new TLweSampleFFTa(N);
        TLweSampleFFTa* tempFFT= new TLweSampleFFTa(N);
        TLweSample64 *temp = new TLweSample64(N);
        Torus64 cn= 0;
        for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialClear_lvl2(accFFT->a+q, env);

        for (int j = 0; j < N; ++j) {
                temp->a[0].coefs[j]=0;
                temp->a[1].coefs[j] =0;

                                      }
           temp->a[1].coefs[0]=c0;

        for (int q = 0; q <= k; ++q)

        TorusPolynomial64_ifft_lvl2(tempFFT->a+q,temp->a+q, env);




         cn= c1-c0;


         for (int j = 0; j < N; ++j)
                temp->a[1].coefs[0] =cn;


        tGsw64DecompH(decomp, temp , env);


        for (int p = 0; p < kpl; ++p)

        IntPolynomial_ifft_lvl2(decompFFT+p,decomp+p, env);




        


         for (int p = 0; p < kpl; ++p)
          for (int q = 0; q <= k; ++q) LagrangeHalfCPolynomialAddMul_lvl2(accFFT->a+q, decompFFT+p, &eps->allsamples[p].a[q], env);


        

   for (int q = k; q <= k; ++q) LagrangeHalfCPolynomialAddTo_lvl2(accFFT->a+q,tempFFT->a+q, env); 


        for (int q = 0; q <= k; ++q)
          for (int j = 0; j < N; ++j) result->a[q].values[j] = accFFT->a[q].values[j];


        delete_array1<IntPolynomiala>(decomp);
        delete_array1<LagrangeHalfCPolynomiala>(decompFFT);
        delete accFFT;
        delete tempFFT;
        delete temp;


}
void CMuxFFTa(TLweSampleFFTa *result, const TGswSampleFFTa *eps, const TLweSampleFFTa *c0, TLweSampleFFTa *c1, const Globals* env){
    const int N = env->N;
    const int k = env->k;
    const int l = env->l;
    const int kpl = (k+1)*l;

    
    struct Cache {
        int N=0, k=0, l=0, kpl=0;
        IntPolynomiala* decomp=nullptr;
        LagrangeHalfCPolynomiala* decompFFT=nullptr;
        TLweSampleFFTa* accFFT=nullptr;
        TLweSample64*  accTorus=nullptr; 
    };
    static thread_local Cache cache;

    auto ensure_cache = [&]() {
        if (cache.accFFT && cache.N==N && cache.k==k && cache.l==l) return;
        if (cache.decomp)    delete_array1<IntPolynomiala>(cache.decomp);
        if (cache.decompFFT) delete_array1<LagrangeHalfCPolynomiala>(cache.decompFFT);
        if (cache.accFFT)    delete cache.accFFT;
        if (cache.accTorus)  delete cache.accTorus;
        cache.N=N; cache.k=k; cache.l=l; cache.kpl=kpl;
        cache.decomp    = new_array1<IntPolynomiala>(kpl, N);
        cache.decompFFT = new_array1<LagrangeHalfCPolynomiala>(kpl, N);
        cache.accFFT    = new TLweSampleFFTa(N);
        cache.accTorus  = new TLweSample64(N);
    };
    ensure_cache();

    
    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialClear_lvl2(cache.accFFT->a + q, env);

    
    for (int q = 0; q <= k; ++q) {
        LagrangeHalfCPolynomialSubTo_lvl2(c1->a + q, c0->a + q, env);
        
        TorusPolynomial64_fft_lvl2(cache.accTorus->a + q, c1->a + q, env);
    }

    
    tGsw64DecompH(cache.decomp, cache.accTorus, env);
    for (int p = 0; p < kpl; ++p)
        IntPolynomial_ifft_lvl2(cache.decompFFT + p, cache.decomp + p, env);

    
    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= k; ++q)
            LagrangeHalfCPolynomialAddMul_lvl2(cache.accFFT->a + q, cache.decompFFT + p, &eps->allsamples[p].a[q], env);

    
    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialAddTo_lvl2(cache.accFFT->a + q, c0->a + q, env);

    
    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            result->a[q].values[j] = cache.accFFT->a[q].values[j];
}
void CMuxDecompFFT(TLweSample64* c0 ,const TGswSampleFFTa *eps, const  LagrangeHalfCPolynomiala* decompFFT, const Globals* env){
    const int N = env->N;
    const int k = env->k;
    const int l = env->l;
    const int kpl = (k+1)*l;

    
    struct Cache {
        int N=0, k=0, l=0, kpl=0;
        TLweSampleFFTa* accFFTa=nullptr;
        TLweSample64*   tmpTorus=nullptr;
    };
    static thread_local Cache cache;

    auto ensure_cache = [&]() {
        if (cache.accFFTa && cache.N==N && cache.k==k && cache.l==l) return;
        if (cache.accFFTa)  delete cache.accFFTa;
        if (cache.tmpTorus) delete cache.tmpTorus;
        cache.N=N; cache.k=k; cache.l=l; cache.kpl=kpl;
        cache.accFFTa = new TLweSampleFFTa(N);
        cache.tmpTorus = new TLweSample64(N);
    };
    ensure_cache();

    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialClear_lvl2(cache.accFFTa->a + q, env);

    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= k; ++q)
            LagrangeHalfCPolynomialAddMul_lvl2(cache.accFFTa->a + q, decompFFT + p, &eps->allsamples[p].a[q], env);

    
    for (int q = 0; q <= k; ++q)
        TorusPolynomial64_fft_lvl2(cache.tmpTorus->a + q, cache.accFFTa->a + q, env);

    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            c0->a[q].coefs[j] += cache.tmpTorus->a[q].coefs[j];
}
void CMuxDecompFFTa(TLweSampleFFTa *c0 ,const TGswSampleFFTa *eps, const  LagrangeHalfCPolynomiala* decompFFT, const Globals* env){
    const int N = env->N;
    const int k = env->k;
    const int l = env->l;
    const int kpl = (k+1)*l;

    
    struct Cache {
        int N=0, k=0, l=0, kpl=0;
        TLweSampleFFTa* accFFTa=nullptr;
    };
    static thread_local Cache cache;

    auto ensure_cache = [&]() {
        if (cache.accFFTa && cache.N==N && cache.k==k && cache.l==l) return;
        if (cache.accFFTa) delete cache.accFFTa;
        cache.N=N; cache.k=k; cache.l=l; cache.kpl=kpl;
        cache.accFFTa = new TLweSampleFFTa(N);
    };
    ensure_cache();

    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialClear_lvl2(cache.accFFTa->a + q, env);

    for (int p = 0; p < kpl; ++p)
        for (int q = 0; q <= k; ++q)
            LagrangeHalfCPolynomialAddMul_lvl2(cache.accFFTa->a + q, decompFFT + p, &eps->allsamples[p].a[q], env);

    for (int q = 0; q <= k; ++q)
        LagrangeHalfCPolynomialAddTo_lvl2(c0->a + q, cache.accFFTa->a + q, env);
}
void shift(TLweSample64 *result, int j,TLweSample64 *sample,int N){





        for (int i=0; i<N; i++){
          if (i+j<N){
  result->a[0].coefs[j+i] = sample->a[0].coefs[i];
   result->a[1].coefs[j+i] = sample->a[1].coefs[i];
          }
          else
           {
   result->a[0].coefs[(i+j)%N]= -sample->a[0].coefs[i];
   result->a[1].coefs[(i+j)%N] = -sample->a[1].coefs[i];

            }


         }
        }

        void tGsw64Encrypt(TGswSample64* cipher, const int mess, const double stdev, const Globals* env){
            const int l = env->l;
            const int Bgbit = env->bgbit;
            const int k=env->k;
         for (int bloc = 0; bloc <= k; ++bloc) {
         for (int i = 0; i < l; ++i) {
            
            tLwe64EncryptZero(&cipher->samples[bloc][i], stdev, env);
            
            cipher->samples[bloc][i].a[bloc].coefs[0] += mess * (UINT64_C(1) << (64-(i+1)*Bgbit));
                                                }
                                     }
          }


    void tGsw64Encrypt_poly_2(TGswSample64* cipher, const IntPolynomiala* mess, const double stdev, const Globals* env){
    const int N = env->N;
    const int l = env->l;
    const int Bgbit = env->bgbit;
    const int k = env->k;

    
        for (int i = 0; i < l; ++i) {
            tLwe64EncryptZero_debug(&cipher->samples[0][i], stdev, env);
            tLwe64EncryptZero_debug(&cipher->samples[1][i], stdev, env);


            for (int j = 0; j < N; ++j) {
                cipher->samples[0][i].a[0].coefs[j] += mess->coefs[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
                cipher->samples[1][i].a[1].coefs[j] += mess->coefs[j] * (UINT64_C(1) << (64 - (i + 1) * Bgbit));
            }
        }
    
}



Torus64 lwe64Phase_lvl2(const LweSample64* cipher, const Globals* env) {
    const int n = env->N;
    Torus64 res = *cipher->b;
    for (int i = 0; i < n; ++i) {
        res -= cipher->a[i]*env->lwekey[i];
    }
    return res;
}


void KSKGen_RGSW(TGswSample64* ksk, const IntPolynomiala* info_sk, const Globals* env) {
    const int l = env->l;            
    const int Bgbit = env->bgbit;   
    const int N = env->N;           
    const double stdev = pow(2., -55); 

    
    std::vector<uint64_t> gadget_vector(l);
    for (int i = 0; i < l; ++i) {
        gadget_vector[i] = (UINT64_C(1) << (64 - (i + 1) * Bgbit));
    }

    
    for (int i = 0; i < l; ++i) {  
        tLwe64EncryptZero(&ksk->samples[0][i], stdev, env);

        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[0].coefs[j] += info_sk->coefs[j] * gadget_vector[i];
            ksk->samples[0][i].a[1].coefs[j] += info_sk->coefs[j] * gadget_vector[i];
        }
    }

    std::cout << "Key switching key generation completed." << std::endl;
}

void KSKGen_RGSW_2_debug(TGswSample64* ksk, const IntPolynomiala* info_sk, const Globals* env) {
    const int l = env->l;            
    const int Bgbit = env->bgbit;     
    const int N = env->N;             
    const double stdev = pow(2., -55); 

    
    std::vector<uint64_t> gadget_vector(l);
    for (int i = 0; i < l; ++i) {
        gadget_vector[i] = (UINT64_C(1) << (64 - (i + 1) * Bgbit));
    }
    
    for (int i = 0; i < l; ++i) {  

        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].b->coefs[j] = 0; 
        }

        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[0].coefs[j] = random_int64();
        }

        torus64PolynomialMultAddKaratsuba_lvl2(ksk->samples[0][i].b, env->tlwekey, &ksk->samples[0][i].a[0], env);

        for (int j = 0; j < N; ++j) {
            ksk->samples[0][i].a[1].coefs[j] += info_sk->coefs[j] * gadget_vector[i];

        }
    }

}


void unpacking_algorithm4(TGswSample64* result, const TLweSample64** rlweInputs, const TGswSample64* convk, const Globals* env) {
    const int l = env->l;            
    const int N = env->N;            

    for (int k = 0; k <= env->k; ++k) {
        for (int i = 0; i < l; ++i) {
            for (int j = 0; j < N; ++j) {
                result->samples[k][i].a[0].coefs[j] = rlweInputs[i]->a[0].coefs[j];
                result->samples[k][i].a[1].coefs[j] = rlweInputs[i]->a[1].coefs[j];
            }
        }
    }

    
    for (int i = 0; i < l; ++i) {
        
        tGswExternMulToTLwe1(&result->samples[0][i], convk, env);

    }
}


void left_shift_by_one(TLweSample64 *output, TLweSample64 *input, int N){
    for(int i=0; i<N; i++){
        output->a[0].coefs[i] = input->a[0].coefs[i+1];
        output->a[1].coefs[i] = input->a[1].coefs[i+1];
    }
    output->a[0].coefs[N-1] = -input->a[0].coefs[0];
    output->a[1].coefs[N-1] = -input->a[1].coefs[0];
}

void left_shift_by_one_poly(IntPolynomiala *output, IntPolynomiala *input, int N){
    for(int i=0; i<N; i++){
        output->coefs[i] = input->coefs[i+1];
    }

}


inline void tau_d_pos_sign(int j, int d, int N,
                           int& out_index, int& out_sign) {
    long long e = 1LL * j * d;
    int q = (int)(e / N);
    int r = (int)(e % N);
    if (r < 0) { r += N; q -= 1; }

    out_index = r;
    out_sign  = (q & 1) ? -1 : 1;
}



void tau_d_secret(IntPolynomiala* out,
                  const IntPolynomiala* in,
                  int d,
                  const Globals* env)
{
    const int N = env->N;

    for (int i = 0; i < N; ++i)
        out->coefs[i] = 0;

    for (int j = 0; j < N; ++j) {
        int coeff = in->coefs[j];
        if (coeff == 0) continue;

        int idx, sgn;
        tau_d_pos_sign(j, d, N, idx, sgn);

        out->coefs[idx] += sgn * coeff;
    }
}


void tau_d_poly(Torus64Polynomial* out,
                const Torus64Polynomial* in,
                int d,
                const Globals* env)
{
    const int N = env->N;

    for (int i = 0; i < N; ++i)
        out->coefs[i] = 0;

    for (int j = 0; j < N; ++j) {
        Torus64 coeff = in->coefs[j];
        if (coeff == 0) continue;

        int idx, sgn;
        tau_d_pos_sign(j, d, N, idx, sgn);

        out->coefs[idx] += (Torus64)sgn * coeff;
    }
}



void tau_d_cipher(TLweSample64* out,
                  const TLweSample64* in,
                  int d,
                  const Globals* env)
{
    const int k = env->k; 

    for (int i = 0; i <= k; ++i) {
        tau_d_poly(&out->a[i], &in->a[i], d, env);
    }
    tau_d_poly(out->b, in->b, d, env);
}




void glwe_encrypt_poly_scaled_level(
    GlweCipher64* ct,
    const IntPolynomiala* msg, 
    int level,                 
    double stdev,
    const Globals* env
) {
    const int N = env->N;
    const int Bgbit = env->bgbit;


    Torus64Polynomial* tor_msg = new Torus64Polynomial(N);

    
    uint64_t scale = UINT64_C(1) << (64 - (level + 1) * Bgbit);

    for (int j = 0; j < N; ++j) {
        
        tor_msg->coefs[j] = (Torus64)msg->coefs[j] * (Torus64)scale;
    }

    
    tLwe64Encrypt(ct, tor_msg, stdev, env);

    delete tor_msg;
}


void glwe_encrypt_poly_scaled_level_noiseless(
    GlweCipher64* ct,
    const IntPolynomiala* msg, 
    int level,                 
    const Globals* env
) {
    const int N     = env->N;
    const int Bgbit = env->bgbit;

    uint64_t scale = UINT64_C(1) << (64 - (level + 1) * Bgbit);


    
    for (int j = 0; j < N; ++j) {
        ct->a[0].coefs[j] = 0;
        ct->a[1].coefs[j] = 0;          
        ct->b->coefs[j]   = (Torus64)msg->coefs[j] * (Torus64)scale;
    }
}



void glwe_encrypt_poly_scaled_level_evalkey(
    GlweCipher64* ct,
    const IntPolynomiala* msg, 
    int level,                 
    double stdev,              
    const Globals* env
) {
    const int N     = env->N;
    const int Bgbit = env->bgbit;
    const int k     = env->k;   

    const uint64_t scale = UINT64_C(1) << (64 - (level + 1) * Bgbit);

    for (int j = 0; j < N; ++j) {
        ct->a[0].coefs[j] = 0;
        if (k >= 1) ct->a[1].coefs[j] = 0;

        Torus64 m_scaled = (Torus64)msg->coefs[j] * (Torus64)scale;

        

        double e_real = sample_gaussian(stdev);

        long double w = (long double)e_real * (long double)((unsigned long long)1 << 63);
        Torus64 e_torus = (Torus64) llround(w);

        
        ct->b->coefs[j] = m_scaled + e_torus;
    }
}



void GLevEncryptPoly(
    GLevCipher64* glev,
    const IntPolynomiala* msg,   
double stdev,
const Globals* env
) {
    const int l = glev->l;
    for (int i = 0; i < l; ++i) {
        glwe_encrypt_poly_scaled_level_evalkey(
            (*glev)[i],   
            msg,
            i, stdev,
            env
        );
    }
}


AutoKsKey64* AutoKsKeyGen64(const Globals* env, int d, int aks_l = -1) {
    const int N = env->N;
    const int l = (aks_l > 0) ? aks_l : env->l;

    const IntPolynomiala* S = env->tlwekey;

    
    IntPolynomiala* S_d = new IntPolynomiala(N);
    tau_d_secret(S_d, S, d, env);

    
    AutoKsKey64* aks = new AutoKsKey64(d, l, N);

    
double ks_stdev = pow(2.,-55);
GLevEncryptPoly(aks->glev, S_d, ks_stdev, env);

    delete S_d;
    return aks;
}


void GLWE_KS_only(
    TLweSample64* out,                 
    const TLweSample64* in,           
    const AutoKsKey64* aks,           
    const Globals* env
) {
    const int N = env->N;
    const int l = env->l;
    const int k = env->k;  

    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            out->a[q].coefs[j] = 0;

    for (int j = 0; j < N; ++j)
        out->b->coefs[j] = in->b->coefs[j];

    
    IntPolynomiala* decomp = new_array1<IntPolynomiala>(l, N);
    tGswTorus64PolynomialDecompH(decomp, &in->a[0], env);

    Torus64Polynomial* tmp = new Torus64Polynomial(N);

    
    for (int p = 0; p < l; ++p) {
        GlweCipher64* kct = (*aks->glev)[p];  
        
        torus64PolynomialMultKaratsuba_lvl2(tmp, &decomp[p], kct->b, env);

        
        for (int j = 0; j < N; ++j)
            out->b->coefs[j] -= tmp->coefs[j];
    }

    delete_array1<IntPolynomiala>(decomp);
    delete tmp;
}


void EvalAuto(
    GlweCipher64* out,
    const GlweCipher64* in,
    const AutoKsKey64* aks,
    const Globals* env
) {
    const int N = env->N;

    GlweCipher64* ct_tau = new GlweCipher64(N);
    tau_d_cipher(ct_tau, in, aks->d, env);  

    GLWE_KS_only(out, ct_tau, aks, env);

    delete ct_tau;
}

inline Torus64 modswitch_down_Torus64(Torus64 x, int bits) {
    if (bits <= 0) return x;
  
    Torus64 half = (Torus64)1 << (bits - 1);  

    if (x >= 0) {
        return (x + half) >> bits;
    } else {
        Torus64 y = -x;
        y = (y + half) >> bits;
        return -y;
    }
}


inline Torus64 modraise_up_Torus64(Torus64 x, int bits) {
    (void)bits;  
    return x;    
}

void GlweModSwitchDown(
    TLweSample64* out,
    const TLweSample64* in,
    int bits,
    const Globals* env
) {
    const int N = env->N;
    const int k = env->k;

    for (int j = 0; j < N; ++j) {
        for (int q = 0; q <= k; ++q) {
            out->a[q].coefs[j] = modswitch_down_Torus64(in->a[q].coefs[j], bits);
        }
    }

    for (int j = 0; j < N; ++j) {
        out->b->coefs[j] = modswitch_down_Torus64(in->b->coefs[j], bits);
    }
}


void GlweModRaiseUp(
    TLweSample64* out,
    const TLweSample64* in,
    int bits,
    const Globals* env
) {
    (void)bits;  

    const int N = env->N;
    const int k = env->k;

    for (int j = 0; j < N; ++j) {
        for (int q = 0; q <= k; ++q) {
            out->a[q].coefs[j] = in->a[q].coefs[j];
        }
    }

    for (int j = 0; j < N; ++j) {
        out->b->coefs[j] = in->b->coefs[j];
    }
}






inline int int_log2(int x) {
    int r = 0;
    while ((1 << r) < x) ++r;
    return r;
}


AutoKsKey64** AutoKsKeyGenAll(const Globals* env, int n, int aks_l = -1) {
    const int N    = env->N;
    const int logN = int_log2(N);
    const int logn = int_log2(n);

    const int num_k = logN - logn;
    AutoKsKey64** aks_list = new AutoKsKey64*[num_k];

    for (int i = 0; i < num_k; ++i) {
        int k = logn + 1 + i;
        int d = (1 << k) + 1;
        aks_list[i] = AutoKsKeyGen64(env, d, aks_l);
    }
    return aks_list;
}

void DeleteAutoKsKeyAll(AutoKsKey64** aks_list, const Globals* env, int n) {
    if (!aks_list) return;

    const int N    = env->N;
    const int logN = int_log2(N);
    const int logn = int_log2(n);
    const int num_k = logN - logn;

    for (int i = 0; i < num_k; ++i) {
        if (aks_list[i]) {
            delete aks_list[i];
        }
    }
    delete[] aks_list;
}









void RevHomTrace_plaintext(
    IntPolynomiala* out,
    const IntPolynomiala* in,
    int n,
    const Globals* env
) {
    const int N    = env->N;
    const int logN = int_log2(N);
    const int logn = int_log2(n);

    
    IntPolynomiala* Cprime = new IntPolynomiala(N);
    IntPolynomiala* Cbar   = new IntPolynomiala(N);
    IntPolynomiala* Cauto  = new IntPolynomiala(N);

    for (int i = 0; i < N; ++i)
        Cprime->coefs[i] = in->coefs[i];

    
    
    
    
    for (int k = logn + 1; k <= logN; ++k) {
        int d = (1 << k) + 1; 

        
        for (int i = 0; i < N; ++i)
            Cbar->coefs[i] = Cprime->coefs[i];

        
        tau_d_secret(Cauto, Cbar, d, env);

        
        for (int i = 0; i < N; ++i)
            Cprime->coefs[i] = Cbar->coefs[i] + Cauto->coefs[i];
    }

    
    for (int i = 0; i < N; ++i)
        out->coefs[i] = Cprime->coefs[i];

    delete Cprime;
    delete Cbar;
    delete Cauto;
}









void RevHomTrace_Alg5(
    GlweCipher64* out,
    const GlweCipher64* in,
    AutoKsKey64** aks_list,
    int n,
    const Globals* env
){
    const int N = env->N;
    const int k = env->k;

    
    int logN = 0;
    {
        int tmp = N;
        while (tmp > 1) { tmp >>= 1; ++logN; }
    }
    int logn = 0;
    {
        int tmp = n;
        while (tmp > 1) { tmp >>= 1; ++logn; }
    }

    
    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j)
            out->a[q].coefs[j] = in->a[q].coefs[j];
    for (int j = 0; j < N; ++j)
        out->b->coefs[j] = in->b->coefs[j];

    if (n == N) return;

    GlweCipher64* Cbar_low = new GlweCipher64(N);
    GlweCipher64* Cbar     = new GlweCipher64(N);
    GlweCipher64* Cauto    = new GlweCipher64(N);

    
    for (int kk = logn + 1; kk <= logN; ++kk) {
        
        GlweModSwitchDown(Cbar_low, out, /*ell=*/1, env);
        
        GlweModRaiseUp   (Cbar,     Cbar_low, /*ell=*/1, env);

        int idx = kk - (logn + 1);           
        AutoKsKey64* aks = aks_list[idx];    

        
        EvalAuto(Cauto, Cbar, aks, env);

        
        for (int q = 0; q <= k; ++q)
            for (int j = 0; j < N; ++j)
                out->a[q].coefs[j] = Cbar->a[q].coefs[j] + Cauto->a[q].coefs[j];

        for (int j = 0; j < N; ++j)
           out->b->coefs[j] = Cbar->b->coefs[j] + Cauto->b->coefs[j];
    }

    delete Cbar_low;
    delete Cbar;
    delete Cauto;
}


struct AutoKsKeyFFTa {
    int d;
    int l;
    LagrangeHalfCPolynomiala* glev_b_fft;

    AutoKsKeyFFTa(int d_in, int l_in, int N)
        : d(d_in), l(l_in), glev_b_fft(new_array1<LagrangeHalfCPolynomiala>(l_in, N)) {}

    ~AutoKsKeyFFTa() { delete_array1<LagrangeHalfCPolynomiala>(glev_b_fft); }
};

static AutoKsKeyFFTa* AutoKsKeyToFFT(const AutoKsKey64* aks, const Globals* env) {
    const int key_l = aks->glev->l;
    AutoKsKeyFFTa* aks_fft = new AutoKsKeyFFTa(aks->d, key_l, env->N);
    for (int p = 0; p < key_l; ++p)
        TorusPolynomial64_ifft_lvl2(&aks_fft->glev_b_fft[p], (*aks->glev)[p]->b, env);
    return aks_fft;
}

static AutoKsKeyFFTa** AutoKsKeyAllToFFT(AutoKsKey64** aks_list, const Globals* env, int n) {
    int logN = 0; for (int t = env->N; t > 1; t >>= 1) ++logN;
    int logn = 0; for (int t = n;      t > 1; t >>= 1) ++logn;
    const int count = logN - logn;
    AutoKsKeyFFTa** out = new AutoKsKeyFFTa*[count];
    for (int i = 0; i < count; ++i) out[i] = AutoKsKeyToFFT(aks_list[i], env);
    return out;
}

static void DeleteAutoKsKeyFFTAll(AutoKsKeyFFTa** aks_fft_list, const Globals* env, int n) {
    if (!aks_fft_list) return;
    int logN = 0; for (int t = env->N; t > 1; t >>= 1) ++logN;
    int logn = 0; for (int t = n;      t > 1; t >>= 1) ++logn;
    const int count = logN - logn;
    for (int i = 0; i < count; ++i) delete aks_fft_list[i];
    delete[] aks_fft_list;
}

static void GLWE_KS_only_fft(TLweSample64* out, const TLweSample64* in,
                              const AutoKsKeyFFTa* aks_fft, const Globals* env) {
    const int N = env->N;
    const int l = aks_fft->l;
    const int k = env->k;

    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j) out->a[q].coefs[j] = 0;
    for (int j = 0; j < N; ++j) out->b->coefs[j] = in->b->coefs[j];

    IntPolynomiala* decomp = new_array1<IntPolynomiala>(l, N);
    LagrangeHalfCPolynomiala* decompFFT = new_array1<LagrangeHalfCPolynomiala>(l, N);
    LagrangeHalfCPolynomiala accFFT(N);
    Torus64Polynomial acc(N);

    tGswTorus64PolynomialDecompH_explicit_l(decomp, &in->a[0], N, l, env->bgbit);
    for (int p = 0; p < l; ++p)
        IntPolynomial_ifft_lvl2(&decompFFT[p], &decomp[p], env);

    LagrangeHalfCPolynomialClear_lvl2(&accFFT, env);
    for (int p = 0; p < l; ++p)
        LagrangeHalfCPolynomialAddMul_lvl2(&accFFT, &decompFFT[p], &aks_fft->glev_b_fft[p], env);
    TorusPolynomial64_fft_lvl2(&acc, &accFFT, env);

    for (int j = 0; j < N; ++j) out->b->coefs[j] -= acc.coefs[j];

    delete_array1<IntPolynomiala>(decomp);
    delete_array1<LagrangeHalfCPolynomiala>(decompFFT);
}

static void EvalAuto_fft(GlweCipher64* out, const GlweCipher64* in,
                         const AutoKsKeyFFTa* aks_fft, const Globals* env) {
    GlweCipher64* ct_tau = new GlweCipher64(env->N);
    tau_d_cipher(ct_tau, in, aks_fft->d, env);
    GLWE_KS_only_fft(out, ct_tau, aks_fft, env);
    delete ct_tau;
}


void RevHomTrace_Alg5_fft(GlweCipher64* out, const GlweCipher64* in,
                           AutoKsKeyFFTa** aks_fft_list, int n, const Globals* env) {
    const int N = env->N;
    const int k = env->k;
    int logN = 0; { int t = N; while (t > 1) { t >>= 1; ++logN; } }
    int logn = 0; { int t = n; while (t > 1) { t >>= 1; ++logn; } }

    for (int q = 0; q <= k; ++q)
        for (int j = 0; j < N; ++j) out->a[q].coefs[j] = in->a[q].coefs[j];
    for (int j = 0; j < N; ++j) out->b->coefs[j] = in->b->coefs[j];

    if (n == N) return;

    GlweCipher64* Cbar_low = new GlweCipher64(N);
    GlweCipher64* Cbar     = new GlweCipher64(N);
    GlweCipher64* Cauto    = new GlweCipher64(N);

    for (int kk = logn + 1; kk <= logN; ++kk) {
        GlweModSwitchDown(Cbar_low, out, 1, env);
        GlweModRaiseUp   (Cbar,     Cbar_low, 1, env);
        const int idx = kk - (logn + 1);
        EvalAuto_fft(Cauto, Cbar, aks_fft_list[idx], env);
        for (int q = 0; q <= k; ++q)
            for (int j = 0; j < N; ++j)
                out->a[q].coefs[j] = Cbar->a[q].coefs[j] + Cauto->a[q].coefs[j];
        for (int j = 0; j < N; ++j)
            out->b->coefs[j] = Cbar->b->coefs[j] + Cauto->b->coefs[j];
    }

    delete Cbar_low;
    delete Cbar;
    delete Cauto;
}
