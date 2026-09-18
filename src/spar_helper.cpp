#include "include/spar_helper.h"

#include <algorithm>

#include "include/bootstrapping.h"
#include "include/hom_placing.h"
#include "include/mphe.h"
#include "src/64header.h"


void MP_RLWE_Encrypt64(TLweSample64* ct,
                               const TLweSample64* common_pk,
                               const Torus64Polynomial* torus_msg,
                               double stdev,
                               const Globals* env) {

    IntPolynomiala* rand_intpoly = new IntPolynomiala(env->N);
    Torus64Polynomial* e0 = new Torus64Polynomial(env->N);
    Torus64Polynomial* e1 = new Torus64Polynomial(env->N);

    for (int i = 0; i < env->N; ++i) {
        rand_intpoly->coefs[i] = random_bit();
        e0->coefs[i] = random_gaussian64(0, stdev);
        e1->coefs[i] = random_gaussian64(0, stdev);
    }

    torus64PolynomialMultKaratsuba_lvl2(&ct->a[0], rand_intpoly, &common_pk->a[0], env);
    torus64PolynomialMultKaratsuba_lvl2(ct->b, rand_intpoly, common_pk->b, env);

    for (int i = 0; i < env->N; ++i) {
        ct->a[0].coefs[i] += e0->coefs[i];
        ct->b->coefs[i] += e1->coefs[i] + torus_msg->coefs[i];
    }

    delete rand_intpoly;
    delete e0;
    delete e1;
}

void ExtractRLWEtoLWE(LweSample64* out,
                                const TLweSample64* in,
                                int idx,
                                const Globals* env) {
    TLweSample64* rotated = new TLweSample64(env->N);
    MulTLWEbyMonomial(rotated, in, -idx, env);
    SampleExtract_AtZero64(out, rotated, env);
    delete rotated;
}



void ChoiceofDistinctThree(int* out, int num, int eta) {
    for (int d = 0; d < num; ++d) {
        while (true) {
            int cand = std::rand() % eta;
            bool used = false;
            for (int p = 0; p < d; ++p) if (out[p] == cand) { used = true; break; }
            if (!used) { out[d] = cand; break; }
        }
    }
}


