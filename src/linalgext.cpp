//
//  linalgext.cpp
//  qme-ng
//
//  Created by Anonymous on 10/15/18.
//  Copyright © 2018 Ying Zhou. All rights reserved.
//

#include "linalgext.hpp"
bool belongsTo(ivect array, vecivect vs) {
    vecivect::iterator vsi;
    long i = 0;
    long length = array.length();
    //bool status = false;
    for (vsi = vs.begin(); vsi != vs.end(); vsi++) {
        for (i = 0; i < length; i++) {
            if ((*vsi)[i] != array[i]) {
                break;
            }
        }
        if (i == length) {
            return true;
        }
    }
    return false;
}
//Return ident matrix of size n * n
//Tested
mat ident(int n) {
    int i = 0, j = 0;
    mat M;
    if (n <= 0)
        throw Exception("Error: The number of rows and the number of columns in an (initialized) matrix must be positive.");
    M.setlength(n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j)
                M[i][j] = 1;
            else
                M[i][j] = 0;
        }
    }
    return M;
}
//All 1 vector with length n
//Tested
vect identv(int n) {
    int i = 0;
    vect v;
    if (n <= 0)
        throw Exception("Error: The size of an (initialized) vector must be positive.");
    v.setlength(n);
    for (i = 0; i < n; i++) {
        v[i] = 1;
    }
    return v;
}
//Calc the sum of a real vector
//Tested
double sumvec(vect v) {
    int i = 0;
    double sum = 0;
    int l = v.length();
    if (!l)
        throw Exception("Error: sumvec can not take an empty array as its argument.");
    for (; i < l; i++) {
        sum += v[i];
    };
    return sum;
}
//Calc the sum of an integer vector
//Tested
int isumvec(ivect v) {
    int i = 0;
    int sum = 0;
    int l = v.length();
    if (!l)
        throw Exception("Error: isumvec can not take an empty array as its argument.");
    for (; i < l; i++) {
        sum += v[i];
    };
    return sum;
}
//Check whether a real vec is actually (almost) an int one
//Tested
bool vecisint(vect v) {
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: vecisint can not take an empty array as its argument.");
    for (; i < l; i++) {
        if (v[i] - round(v[i]) > EPSILON || v[i] - round(v[i]) < -EPSILON)
            return false;
    }
    return true;
}
//Check whether an int vec is non negative
//Tested
bool ivecisnonneg(ivect v) {
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: ivecisnonneg can not take an empty array as its argument.");
    for (; i < l; i++) {
        if (v[i] < 0)
            return false;
    }
    return true;
}
//Rounding real vectors
//Tested
ivect intizev(vect v) {
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: intizev can not take an empty array as its argument.");
    ivect iv;
    iv.setlength(l);
    for (; i < l; i++) {
        iv[i] = round(v[i]);
    }
    return iv;
}
//Rounding real matrices
//Tested
imat intizem(mat m) {
    int i = 0, j = 0;
    long r = m.rows(), c = m.cols();
    if (!r || !c)
        throw Exception("Error: intizem can not take an empty array as its argument.");
    imat im;
    im.setlength(r, c);
    for (; i < r; i++) {
        for(j = 0; j < c; j++)
            im[i][j] = round(m[i][j]);
    }
    return im;
}
//Casting int vector to real
//Tested
vect rizev(ivect v) {
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: Can not apply rizev to an empty vector.");
    vect rv;
    rv.setlength(l);
    for (; i < l; i++) {
        rv[i] = v[i];
    }
    return rv;
}
//Casting int matrix to real
//Tested
mat rizem(imat m) {
    int i = 0, j = 0;
    long r = m.rows(), c = m.cols();
    if (!r || !c)
        throw Exception("Error: rizem can not take an empty array as its argument.");
    mat rm;
    rm.setlength(r, c);
    for (; i < r; i++) {
        for(j = 0; j < c; j++)
            rm[i][j] = m[i][j];
    }
    return rm;
}
//Convert an int vector of length n to a n * 1 int matrix
//Tested
imat vec2matc(ivect v) {
    imat M;
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: Can not apply vec2matc to an empty vector.");
    M.setlength(l, 1);
    for (; i < l; i++) {
        M[i][0] = v[i];
    }
    return M;
}
//Convert a vector of length n to a 1 * n matrix
//Tested
imat vec2matr(ivect v) {
    imat M;
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: Can not apply vec2matr to an empty vector.");
    M.setlength(1, l);
    for (; i < l; i++) {
        M[0][i] = v[i];
    }
    return M;
}
//Tested
ivect add(ivect v1, ivect v2) {
    int i = 0;
    long l1 = v1.length();
    long l2 = v2.length();
    ivect rv;
    if (!l1 || !l2)
        throw Exception("Error: Can not add an empty vector to any vector");
    if (l1 != l2) {
        throw Exception("Error: The two vectors have different lengths!");
    }
    rv.setlength(l1);
    for (; i < l1; i++) {
        rv[i] = v1[i] + v2[i];
    }
    return rv;
}
//Take the k-th (starting from 0) column of matrix M
//Tested
ivect matc2vec(imat M, int k) {
    ivect v;
    int i = 0;
    long l1 = M.rows();
    long l2 = M.cols();
    if (!l1 || !l2)
        throw Exception("Error: Can not apply matr2vec to an empty matrix.");
    if (k < 0 || k >= l2)
        throw Exception("Error: The index is out of range.");
    v.setlength(l1);
    for (; i < l1; i++) {
        v[i] = M[i][k];
    }
    return v;
}
//Take the k-th (starting from 0) row of matrix M
//Tested
ivect matr2vec(imat M, int k) {
    ivect v;
    int i = 0;
    long l1 = M.rows();
    long l2 = M.cols();
    if (!l1 || !l2)
        throw Exception("Error: Can not apply matr2vec to an empty matrix.");
    if (k < 0 || k >= l1)
        throw Exception("Error: The index is out of range.");
    v.setlength(l2);
    for (; i < l2; i++) {
        v[i] = M[k][i];
    }
    return v;
}
//Convert an real vector of length n to a n * 1 real matrix
//Tested
mat vec2matc(vect v) {
    mat M;
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: Can not apply vec2matc to an empty vector.");
    M.setlength(l, 1);
    for (; i < l; i++) {
        M[i][0] = v[i];
    }
    return M;
}
//Take the k-th (starting from 0) column of matrix M
//Tested
vect matc2vec(mat M, int k) {
    vect v;
    int i = 0;
    long l1 = M.rows();
    long l2 = M.cols();
    if (!l1 || !l2)
        throw Exception("Error: Can not apply matr2vec to an empty matrix.");
    if (k < 0 || k >= l2)
        throw Exception("Error: The index is out of range.");
    v.setlength(l1);
    for (; i < l1; i++) {
        v[i] = M[i][k];
    }
    return v;
}
//Is an int vector sincere?
//Tested
bool isSincere(ivect v) {
    int i = 0;
    long l = v.length();
    if (!l)
        throw Exception("Error: Can not apply isSincere to an empty array.");
    for (; i < l; i++) {
        if (v[i] == 0)
            return false;
    }
    return true;
}
ivect matvecmult(imat M, ivect v) {
    mat ress;
    imat vv = vec2matc(v);
    long len = M.rows();
    if (M.cols() != v.length()) {
        throw Exception("Error: The matrix and the vector are not compatible.");
    }
    ress.setlength(len, 1);
    alglib::rmatrixgemm(len,1,M.cols(), 1, rizem(M),0,0,0, rizem(vv),0,0,0,0,ress,0,0);
    return matc2vec(intizem(ress), 0);
}

imat mult(imat M1, imat M2) {
    mat ress;
    long len = M1.cols();
    if (M2.rows() != len) {
        throw Exception("Error: The matrices are not compatible.");
    }
    ress.setlength(M1.rows(), M2.cols());
    alglib::rmatrixgemm(M1.rows(),M2.cols(),len, 1, rizem(M1),0,0,0, rizem(M2),0,0,0,0,ress,0,0);
    return intizem(ress);
}

vect matvecmult(mat M, vect v) {
    mat ress;
    mat vv = vec2matc(v);
    long len = M.rows();
    if (M.cols() != v.length()) {
        throw Exception("Error: The matrix and the vector are not compatible.");
    }
    ress.setlength(len, 1);
    alglib::rmatrixgemm(len,1,M.cols(), 1, M,0,0,0, vv,0,0,0,0,ress,0,0);
    return matc2vec(ress, 0);
}

void printVector(ivect vec) {
    long i = 0;
    long n = vec.length();
    for (; i<n; i++) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

void printVectors(vecivect vs) {
    vecivect::iterator vsi;
    for (vsi = vs.begin(); vsi != vs.end(); vsi++) {
        printVector(*vsi);
    }
}

bool equals(ivect v1, ivect v2) {
    long l = v1.length();
    int i = 0;
    if (v2.length() != l)
        return false;
    for (; i < l; i++) {
        if (v1[i] != v2[i])
            return false;
    }
    return true;
}

int mFinDeg(imat M, ivect v) {
    const int MAX = 40 * isumvec(v);
    int i = 0, sum;
    ivect temp = v;
    do {
        temp = matvecmult(M,temp);
        //std::cout << "temp is " << temp.tostring() << " sum is " << isumvec(temp) << std::endl;
        sum = isumvec(temp);
        if (sum < 0) {
            //std::cout << "Success vec: " << temp.tostring() << std::endl;
            return i;
        }
        i++;
        if (equals(temp, v))//tame regulars
            return -i;
    } while(sum <= MAX);
    //std::cout << "Exit vec: " << temp.tostring() << std::endl;
    //std::cout << "Exit sum: " << sum << std::endl;
    return -1;
}
ivect minvec(imat M, imat invM, ivect v) {
    int minsum = isumvec(v);
    const int MAX = 40 * minsum;
    int sum;
    bool sincerity = isSincere(v), tsincerity;
    ivect temp = v;
    ivect best = v;
    do {
        temp = matvecmult(M,temp);
        //std::cout << "temp is " << temp.tostring() << " sum is " << isumvec(temp) << std::endl;
        sum = isumvec(temp);
        tsincerity = isSincere(v);
        if ((sincerity && !tsincerity) || ((sincerity ^ tsincerity) && (minsum > sum))) {
            best = temp;
            minsum = sum;
        }
    } while(sum <= MAX);
    temp = v;
    do {
        temp = matvecmult(invM,temp);
        //std::cout << "temp is " << temp.tostring() << " sum is " << isumvec(temp) << std::endl;
        sum = isumvec(temp);
        tsincerity = isSincere(v);
        if ((sincerity && !tsincerity) || ((sincerity ^ tsincerity) && (minsum > sum))) {
            best = temp;
            minsum = sum;
        }
    } while(sum <= MAX);
    return best;
}
//Only go from v1 to v2
int sameMSemiOrbit (imat M, ivect v1, ivect v2) {
    int sum1 = isumvec(v1), sum2 = isumvec(v2);
    int maxsum = sum1 > sum2 ? sum1 : sum2;
    const int MAX = 40 * maxsum;
    int sum;
    int degdiff = 0;
    ivect temp = v1;
    //std::cout << "temp is " << temp.tostring() << " sum is " << isumvec(temp) << std::endl;
    do {
        degdiff++;
        temp = matvecmult(M,temp);
        sum = isumvec(temp);
        //std::cout << "temp is " << temp.tostring() << " sum is " << sum << std::endl;
        if (equals(temp, v2))
            return degdiff;
    } while(sum <= MAX);
    return -1;
}
//Calculate I+M+M^2+\cdots+M^k
//Tested
imat mppo(imat M, int k) {
    int i = 0;
    int n = M.rows();
    imat temp = intizem(ident(n));
    mat id = ident(n);
    if (M.cols() != n)
        throw Exception("Error: The matrix isn't a square matrix.");
    if (k < 0)
        throw Exception("Error: k must be a nonnegative number.");
    if (k == 0)
        return temp;
    for (; i < k; i++) {
        alglib::rmatrixgemm(n,n,n, 1, rizem(M),0,0,0, rizem(temp),0,0,0,1,id,0,0);
        temp = intizem(id);
        id = ident(n);
    }
    return temp;
}

