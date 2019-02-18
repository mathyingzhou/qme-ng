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
mat ident(int n) {
    int i = 0, j = 0;
    mat M;
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
vect identv(int n) {
    int i = 0;
    vect v;
    v.setlength(n);
    for (i = 0; i < n; i++) {
        v[i] = 1;
    }
    return v;
}
double sumvec(vect v) {
    int i = 0;
    int sum = 0;
    int l = v.length();
    //std::cout << "Sum of vector " << v.tostring(3) << std::endl;
    for (; i < l; i++) {
        sum += v[i];
        //std::cout << v[i] << " " << sum << " ";
    };
    //std::cout << std::endl;
    return sum;
}
int isumvec(ivect v) {
    int i = 0;
    int sum = 0;
    int l = v.length();
    for (; i < l; i++) {
        sum += v[i];
    };
    //std::cout << std::endl;
    return sum;
}
bool vecisint(vect v) {
    int i = 0;
    long l = v.length();
    for (; i < l; i++) {
        if (v[i] - round(v[i]) > EPSILON || v[i] - round(v[i]) < -EPSILON)
            return false;
    }
    return true;
}
bool ivecisnonneg(ivect v) {
    int i = 0;
    long l = v.length();
    for (; i < l; i++) {
        if (v[i] < 0)
            return false;
    }
    return true;
}
ivect intizev(vect v) {
    int i = 0;
    long l = v.length();
    ivect iv;
    iv.setlength(l);
    for (; i < l; i++) {
        iv[i] = round(v[i]);
    }
    return iv;
}
imat intizem(mat m) {
    int i = 0, j = 0;
    long r = m.rows(), c = m.cols();
    imat im;
    im.setlength(r, c);
    for (; i < r; i++) {
        for(j = 0; j < c; j++)
            im[i][j] = round(m[i][j]);
    }
    return im;
}
vect rizev(ivect v) {
    int i = 0;
    long l = v.length();
    vect rv;
    rv.setlength(l);
    for (; i < l; i++) {
        rv[i] = v[i];
    }
    return rv;
}
mat rizem(imat m) {
    int i = 0, j = 0;
    long r = m.rows(), c = m.cols();
    mat rm;
    rm.setlength(r, c);
    for (; i < r; i++) {
        for(j = 0; j < c; j++)
            rm[i][j] = m[i][j];
    }
    return rm;
}
imat vec2matc(ivect v) {
    imat M;
    int i = 0;
    long l = v.length();
    M.setlength(l, 1);
    for (; i < l; i++) {
        M[i][0] = v[i];
    }
    return M;
}
imat vec2matr(ivect v) {
    imat M;
    int i = 0;
    long l = v.length();
    M.setlength(1, l);
    for (; i < l; i++) {
        M[0][i] = v[i];
    }
    return M;
}
ivect add(ivect v1, ivect v2) {
    int i = 0;
    long l = v1.length();
    ivect rv;
    if (v2.length() != l) {
        std::cerr << "The two vectors have different lengths!" << std::endl;
    }
    rv.setlength(l);
    for (; i < l; i++) {
        rv[i] = v1[i] + v2[i];
    }
    return rv;
}
ivect matc2vec(imat M, int k) {
    ivect v;
    int i = 0;
    long l = M.rows();
    v.setlength(l);
    for (; i < l; i++) {
        v[i] = M[i][k];
    }
    return v;
}
ivect matr2vec(imat M, int k) {
    ivect v;
    int i = 0;
    long l = M.rows();
    v.setlength(l);
    for (; i < l; i++) {
        v[i] = M[k][i];
    }
    return v;
}
mat vec2matc(vect v) {
    mat M;
    int i = 0;
    long l = v.length();
    M.setlength(l, 1);
    for (; i < l; i++) {
        M[i][0] = v[i];
    }
    return M;
}
vect matc2vec(mat M, int k) {
    vect v;
    int i = 0;
    long l = M.rows();
    v.setlength(l);
    for (; i < l; i++) {
        v[i] = M[i][k];
    }
    return v;
}
bool isSincere(ivect v) {
    int i = 0;
    long l = v.length();
    for (; i < l; i++) {
        if (v[i] == 0)
            return false;
    }
    return true;
}
ivect matvecmult(imat M, ivect v) {
    ivect def;
    mat ress;
    imat vv = vec2matc(v);
    long len = M.rows();
    if (M.cols() != v.length()) {
        std::cout << "Wrong sizes!" << std::endl;
        return def;
    }
    ress.setlength(len, 1);
    alglib::rmatrixgemm(len,1,M.cols(), 1, rizem(M),0,0,0, rizem(vv),0,0,0,0,ress,0,0);
    return matc2vec(intizem(ress), 0);
}

imat mult(imat M1, imat M2) {
    imat def;
    mat ress;
    long len = M1.cols();
    if (M2.rows() != len) {
        std::cout << "Wrong sizes!" << std::endl;
        return def;
    }
    ress.setlength(M1.rows(), M2.cols());
    alglib::rmatrixgemm(M1.rows(),M2.cols(),len, 1, rizem(M1),0,0,0, rizem(M2),0,0,0,0,ress,0,0);
    return intizem(ress);
}

vect matvecmult(mat M, vect v) {
    vect def;
    mat ress;
    mat vv = vec2matc(v);
    long len = M.rows();
    if (M.cols() != v.length()) {
        std::cout << "Wrong sizes!" << std::endl;
        return def;
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
    int l = v1.length();
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
imat mppo(imat M, int k) {
    int i = 0;
    int n = M.rows();
    imat temp = intizem(ident(n));
    mat id = ident(n);
    if (M.cols() != n)
        return temp;
    if (k <= 0)
        return temp;
    for (; i < k; i++) {
        alglib::rmatrixgemm(n,n,n, 1, rizem(M),0,0,0, rizem(temp),0,0,0,1,id,0,0);
        temp = intizem(id);
        id = ident(n);
    }
    return temp;
}

