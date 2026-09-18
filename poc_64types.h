#ifndef _POC_64TYPESH_
#define _POC_64TYPESH_





#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>



#include "generic_utils.h"




typedef int32_t Torus32; 
typedef int64_t Torus64; 

static const int64_t _two32 = INT64_C(1) << 32; 

inline Torus32 t64tot32(Torus64 x) {
    return int32_t(x/_two32);
}
inline Torus64 t32tot64(Torus32 x) {
    return int64_t(x)*_two32;
}



struct PolynomialParameter32 {
    const int N;
    void* const FFT_PREPROC;
};
struct PolynomialParameter64 {
    const int N;
    void* const FFT_PREPROC;
};



struct Torus32Polynomial {
    Torus32* const coefs;
    
    Torus32Polynomial(int N): coefs(new Torus32[N]) {}
    
    ~Torus32Polynomial() { delete[] coefs; }
};



Torus32Polynomial* new_Torus32Polynomial(int N) {
    Torus32Polynomial* obj = (Torus32Polynomial*) malloc(sizeof(Torus32Polynomial));
    new(obj) Torus32Polynomial(N);
    return obj;
}
Torus32Polynomial* new_Torus32Polynomial_array(int nbelts, int N) {
    Torus32Polynomial* obj = (Torus32Polynomial*) malloc(nbelts*sizeof(Torus32Polynomial));
    for (int i = 0; i < nbelts; i++) new(obj+i) Torus32Polynomial(N);
    return obj;
}

void delete_Torus32Polynomial(Torus32Polynomial* obj) {
    obj->~Torus32Polynomial();
    free(obj);
}
void delete_Torus32Polynomial_array(int nbelts, Torus32Polynomial* obj) {
    for (int i = 0; i < nbelts; i++) (obj+i)->~Torus32Polynomial();
    free(obj);
}





struct Torus64Polynomial {
    Torus64* const coefs;
    
    Torus64Polynomial(int N): coefs(new Torus64[N]) {}
    
    ~Torus64Polynomial() { delete[] coefs; }
};





struct IntPolynomiala {
    int* const coefs;
    
    IntPolynomiala(int N): coefs(new int[N]) {}
    
    ~IntPolynomiala() { delete[] coefs; }
};



#ifdef USE_FFT
struct LagrangeHalfCPolynomiala {
    double* const values;
    LagrangeHalfCPolynomiala(int N): values(new double[N]) {}
    ~LagrangeHalfCPolynomiala() {
        delete[] values;
    }
};

#else

struct LagrangeHalfCPolynomiala {
    IntPolynomiala* intPoly;
    Torus64Polynomial* torus64Poly;
    LagrangeHalfCPolynomiala(int N) {
        intPoly=0;
        torus64Poly=0;
    }
    void clear() {
        if (intPoly) { delete intPoly; intPoly=0; }
        if (torus64Poly) { delete torus64Poly; torus64Poly=0;}
    }
    void setIntPoly(const IntPolynomiala* a, const int N) {
        clear();
        intPoly = new IntPolynomiala(N);
        for (int i=0; i<N; i++) intPoly->coefs[i]=a->coefs[i];
    }
    void setTorus64Poly(const Torus64Polynomial* a, const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=a->coefs[i];
    }
    void setZeroTorus64Poly(const int N) {
        clear();
        torus64Poly = new Torus64Polynomial(N);
        for (int i=0; i<N; i++) torus64Poly->coefs[i]=0;
    }
    ~LagrangeHalfCPolynomiala() {
    }
};

#endif



struct LweSample32 {
    Torus32* const a;
    Torus32* const b; 
    
    LweSample32(int n): a(new Torus32[n+1]), b(&a[n]) {}
    
    ~LweSample32() { delete[] a; }
};






struct LweSample64 {
    Torus64* const a;
    Torus64* const b; 
    
    LweSample64(int n): a(new Torus64[n+1]), b(&a[n]) {}
    
    ~LweSample64() { delete[] a; }
};




/*
struct Torus64Polynomial {
    Torus64* const coefs;
    
    Torus64Polynomial(int N): coefs(new Torus64[N]) {}
    
    ~Torus64Polynomial() { delete[] coefs; }
};
*/


struct TLweSample32 {
    Torus32Polynomial* const a;
    Torus32Polynomial* const b; 
    
    TLweSample32(int N): a(new_array1<Torus32Polynomial>(2,N)), b(&a[1]) {}
    
    ~TLweSample32() { delete_array1<Torus32Polynomial>(a); }
};



struct TLweSample64 {
    Torus64Polynomial* const a;
    Torus64Polynomial* const b; 
    
    TLweSample64(int N): a(new_array1<Torus64Polynomial>(2,N)), b(&a[1]) {}
    
    ~TLweSample64() { delete_array1<Torus64Polynomial>(a); }
};





struct TLweSampleFFTa {
    LagrangeHalfCPolynomiala* const a;
    LagrangeHalfCPolynomiala* const b; 
    
    TLweSampleFFTa(int N): a(new_array1<LagrangeHalfCPolynomiala>(2,N)), b(&a[1]) {}
    
    ~TLweSampleFFTa() { delete_array1<LagrangeHalfCPolynomiala>(a); }
};








struct TGswSample32 {
    TLweSample32** const samples; 
    TLweSample32* const allsamples;
    
    TGswSample32(int l, int N):
        samples(new_array2<TLweSample32>(2,l,N)),
        allsamples(samples[0]) {}
    
    ~TGswSample32() {
        delete_array2<TLweSample32>(samples);
    }
};





struct TGswSample64 {
    TLweSample64** const samples; 
    TLweSample64* const allsamples;
    
    TGswSample64(int l, int N):
        samples(new_array2<TLweSample64>(2,l,N)),
        allsamples(samples[0]) {}
    
    ~TGswSample64() {
        delete_array2<TLweSample64>(samples);
    }
};





struct TGswSampleFFTa {
    TLweSampleFFTa** const samples; 
    TLweSampleFFTa* const allsamples;
    
    TGswSampleFFTa(int l, int N):
        samples(new_array2<TLweSampleFFTa>((2),l,N)),
        allsamples(samples[0]) {}
    
    ~TGswSampleFFTa() {
        delete_array2<TLweSampleFFTa>(samples);
    }
};


using GlweSecretKey = IntPolynomiala;  


struct GlweCipher64 : public TLweSample64 {
    GlweCipher64(int N) : TLweSample64(N) {}
};


struct GLevCipher64 {
    int l;                 
    GlweCipher64** cts;    

    GLevCipher64(int l_in, int N) : l(l_in) {
        cts = new GlweCipher64*[l];
        for (int i = 0; i < l; ++i)
            cts[i] = new GlweCipher64(N);
    }

    ~GLevCipher64() {
        for (int i = 0; i < l; ++i) delete cts[i];
        delete[] cts;
    }

    GlweCipher64* operator[](int i)       { return cts[i]; }
    const GlweCipher64* operator[](int i) const { return cts[i]; }
};




struct AutoKsKey64 {
    int d;               
    GLevCipher64* glev;  

    AutoKsKey64(int d_in, int l, int N) : d(d_in) {
        glev = new GLevCipher64(l, N);
    }
    ~AutoKsKey64() {
        delete glev;
    }
};





/*
class Globals {
    public:
        
        static const int n_lvl0;
        static const int n_lvl1;
        static const int n_lvl2;
        static const int bgbit_lvl1;
        static const int ell_lvl1;
        static const int bgbit_lvl2;
        static const int ell_lvl2;
        static const double bkstdev_lvl2;
        static const double ksstdev_lvl10;
        static const int kslength_lvl10;
        static const int ksbasebit_lvl10;
        static const double ksstdev_lvl21;
        static const int kslength_lvl21;
        static const int ksbasebit_lvl21;

        
        int t_lvl0; 
        int t_lvl1; 
        int N_lvl1; 
        int N_lvl2; 
        
        uint64_t torusDecompOffset;
        
        uint64_t* torusDecompBuf;

        
        int* key_lvl0;
        int* key_lvl1;
        IntPolynomial* Key_lvl1;
        int* key_lvl2;
        IntPolynomial* Key_lvl2;

        
        
        LweSample32*** preKS; 
        
        TGswSample64* bk; 
        TGswSampleFFT* bkFFT; 
        
        TLweSample32**** privKS;

        Globals();
};
*/


class Globals {
    public:
        
        int k;
        int N;
        int t;
        int smalln;
        int bgbit;
        int l;
        int basebit;

       
    uint64_t torusDecompOffset;
        
    uint64_t* torusDecompBuf;
      int* lwekey;
      IntPolynomiala* tlwekey;
        int64_t* in_key;
         
       
        
        TLweSample64**** privKS;
        
       
        
        

        Globals(int N_in = 2048,
                int bgbit_in = 9,
                int l_in = 3,
                int basebit_in = 8,
                int t_in = 256,
                int smalln_in = 2048,
                int k_in = 1);
};


#endif
