//
//  linalgext.hpp
//  qme-ng
//
//  Created by Anonymous on 10/15/18.
//  Copyright © 2018 Ying Zhou. All rights reserved.
//

#ifndef linalgext_hpp
#define linalgext_hpp

#include <stdio.h>
#include "ap.h"
#include "linalg.h"
typedef alglib::real_1d_array vect;
typedef alglib::integer_1d_array ivect;
typedef std::vector<ivect> vecivect;
typedef alglib::real_2d_array mat;
typedef alglib::integer_2d_array imat;
#define EPSILON 0.001
double sumvec(vect);
int isumvec(ivect);
bool isSincere(ivect);
ivect intizev(vect);
imat intizem(mat);
imat mppo(imat, int);
vect rizev(ivect);
mat rizem(imat);
mat ident(int);
vect identv(int);
bool vecisint(vect);
bool ivecisnonneg(ivect);
imat vec2matc(ivect);
imat vec2matr(ivect);
ivect add(ivect, ivect);
ivect matc2vec(imat, int);
ivect matr2vec(imat, int);
ivect matvecmult(imat, ivect);
imat mult(imat, imat);
inline ivect mult(ivect v, imat M) {
    return matr2vec(mult(vec2matr(v), M), 0);
}
inline ivect mult(imat M, ivect v) {
    return matc2vec(mult(M, vec2matc(v)), 0);
}
bool equals(ivect, ivect);
vect matvecmult(mat, vect);
extern void printVector(ivect);
extern void printVectors(vecivect);
bool belongsTo(ivect, vecivect);
int mFinDeg(imat, ivect);
ivect minvec(imat, imat, ivect);
int sameMSemiOrbit (imat, ivect, ivect);
inline int sameMOrbit (imat M, ivect v1, ivect v2) {
    int d1 = sameMSemiOrbit(M,v1,v2), d2 = sameMSemiOrbit(M,v2,v1);
    if (d1 != -1)
        return d1;
    else if (d2 != -1)
        return -d2;
    return 0;
}
#endif /* linalgext_hpp */
