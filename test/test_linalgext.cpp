//
//  test_linalgext.cpp
//  test_linalgext
//
//  Created by Anonymous on 3/3/19.
//  Copyright © 2019 Nyaa Studio. All rights reserved.
//
#define BOOST_TEST_MODULE test_linalgext
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include "linalgext.hpp"
#include "Exception.h"
namespace bdata = boost::unit_test::data;
namespace utf = boost::unit_test;
//Disabled by default: Mppo, Rizem
//TODO: Identv & Ident need to throw exceptions
//sumvec
BOOST_AUTO_TEST_SUITE(SumVec);
BOOST_AUTO_TEST_CASE(SumVecTestLen0) {
    vect v;
    BOOST_CHECK_THROW(sumvec(v), Exception);
}
BOOST_AUTO_TEST_CASE(SumVecTestLen0Explicit) {
    vect v;
    v.setlength(0);
    BOOST_CHECK_THROW(sumvec(v), Exception);
}
BOOST_DATA_TEST_CASE(SumVecTestLen1, bdata::xrange<double>(-10, 10, 0.15), rand1) {
    vect v;
    v.setlength(1);
    v[0] = rand1;
    BOOST_CHECK_CLOSE(sumvec(v), rand1, 0.0001f);
}
BOOST_DATA_TEST_CASE(SumVecTestLen2, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    BOOST_CHECK_CLOSE(sumvec(v), rand1 + rand2, 0.0001f);
}
BOOST_DATA_TEST_CASE(SumVecTestLen3, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^  bdata::xrange<double>(-10, 10, 0.1), rand1, rand2, rand3) {
    vect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    BOOST_CHECK_CLOSE(sumvec(v), rand1 + rand2 + rand3, 0.0001f);
}
BOOST_AUTO_TEST_SUITE_END();
//isumvec
BOOST_AUTO_TEST_SUITE(ISumVec);
BOOST_AUTO_TEST_CASE(ISumVecTestLen0) {
    ivect v;
    BOOST_CHECK_THROW(isumvec(v), Exception);
}
BOOST_AUTO_TEST_CASE(ISumVecTestLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(isumvec(v), Exception);
}
BOOST_DATA_TEST_CASE(ISumVecTestLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    BOOST_CHECK_EQUAL(isumvec(v), rand1);
}
BOOST_DATA_TEST_CASE(ISumVecTestLen2, bdata::xrange<int>(-50, 50, 1) ^ bdata::xrange<int>(-50, 50, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    BOOST_CHECK_EQUAL(isumvec(v), rand1 + rand2);
}
BOOST_DATA_TEST_CASE(ISumVecTestLen3, bdata::xrange<int>(-50, 50, 1) ^ bdata::xrange<int>(-50, 50, 1) ^ bdata::xrange<int>(-50, 50, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    BOOST_CHECK_EQUAL(isumvec(v), rand1 + rand2 + rand3);
}
BOOST_AUTO_TEST_SUITE_END();
//isSincere
BOOST_AUTO_TEST_SUITE(Sincere);
BOOST_AUTO_TEST_CASE(SincereLen0) {
    ivect v;
    BOOST_CHECK_THROW(isSincere(v), Exception);
}
BOOST_AUTO_TEST_CASE(SincereLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(isSincere(v), Exception);
}
BOOST_DATA_TEST_CASE(SincereLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    if (rand1)
        BOOST_CHECK(isSincere(v));
    else
        BOOST_CHECK(!isSincere(v));
}
BOOST_DATA_TEST_CASE(SincereLen2, bdata::xrange<int>(-20, 20, 1) * bdata::xrange<int>(-20, 20, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    if (rand1 && rand2)
        BOOST_CHECK(isSincere(v));
    else
        BOOST_CHECK(!isSincere(v));
}
BOOST_DATA_TEST_CASE(SincereLen3, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    if (rand1 && rand2 && rand3)
        BOOST_CHECK(isSincere(v));
    else
        BOOST_CHECK(!isSincere(v));
}
BOOST_AUTO_TEST_SUITE_END();
//Intizev
BOOST_AUTO_TEST_SUITE(Intizev);
BOOST_AUTO_TEST_CASE(IntizevLen0) {
    vect v;
    BOOST_CHECK_THROW(intizev(v), Exception);
}
BOOST_AUTO_TEST_CASE(IntizevLen0Explicit) {
    vect v;
    v.setlength(0);
    BOOST_CHECK_THROW(intizev(v), Exception);
}
BOOST_DATA_TEST_CASE(IntizevTestLen1, bdata::xrange<double>(-10, 10, 0.15), rand1) {
    vect v;
    v.setlength(1);
    v[0] = rand1;
    ivect w = intizev(v);
    BOOST_CHECK_EQUAL(w.length(), 1);
    BOOST_CHECK_CLOSE(w[0], round(rand1), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizevTestLen2, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    ivect w = intizev(v);
    BOOST_CHECK_EQUAL(w.length(), 2);
    BOOST_CHECK_CLOSE(w[0], round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w[1], round(rand2), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizevTestLen3, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^  bdata::xrange<double>(-10, 10, 0.1), rand1, rand2, rand3) {
    vect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    ivect w = intizev(v);
    BOOST_CHECK_EQUAL(w.length(), 3);
    BOOST_CHECK_CLOSE(w[0], round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w[1], round(rand2), 0.0001f);
    BOOST_CHECK_CLOSE(w[2], round(rand3), 0.0001f);
}
BOOST_AUTO_TEST_SUITE_END();
//Intizem
BOOST_AUTO_TEST_SUITE(Intizem);
BOOST_AUTO_TEST_CASE(IntizemLen00) {
    mat v;
    BOOST_CHECK_THROW(intizem(v), Exception);
}
BOOST_DATA_TEST_CASE(IntizemLen0KExplicit, bdata::xrange<int>(0,10,1), rand1) {
    mat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(intizem(v), Exception);
}
BOOST_DATA_TEST_CASE(IntizemLenK0Explicit, bdata::xrange<int>(0,10,1), rand1) {
    mat v;
    v.setlength(rand1, 0);
    BOOST_CHECK_THROW(intizem(v), Exception);
}
BOOST_DATA_TEST_CASE(IntizemTestLen11, bdata::xrange<double>(-10, 10, 0.15), rand1) {
    mat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizemTestLen12, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2) {
    mat v;
    v.setlength(1,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), round(rand2), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizemTestLen21, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2) {
    mat v;
    v.setlength(2,1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), round(rand2), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizemTestLen13, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^  bdata::xrange<double>(-10, 10, 0.1), rand1, rand2, rand3) {
    mat v;
    v.setlength(1,3);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(0,2) = rand3;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 3);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), round(rand2), 0.0001f);
    BOOST_CHECK_CLOSE(w(0,2), round(rand3), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizemTestLen31, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^  bdata::xrange<double>(-10, 10, 0.1), rand1, rand2, rand3) {
    mat v;
    v.setlength(3,1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    v(2,0) = rand3;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 3);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), round(rand2), 0.0001f);
    BOOST_CHECK_CLOSE(w(2,0), round(rand3), 0.0001f);
}
BOOST_DATA_TEST_CASE(IntizemTestLen22, bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::random(bdata::distribution=std::uniform_real_distribution<float>(-10, 10)) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2, rand3, rand4) {
    mat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    imat w = intizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_CLOSE(w(0,0), round(rand1), 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), round(rand2), 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), round(rand3), 0.0001f);
    BOOST_CHECK_CLOSE(w(1,1), round(rand4), 0.0001f);
}
BOOST_AUTO_TEST_SUITE_END();
//mppo
BOOST_AUTO_TEST_SUITE(Mppo, * utf::disabled());
BOOST_DATA_TEST_CASE(mppoLen00, bdata::xrange<int>(-10,10,1), rand1) {
    imat v;
    BOOST_CHECK_THROW(mppo(v,rand1), Exception);
}
BOOST_DATA_TEST_CASE(mppoLen0KExplicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(mppo(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(mppoLenK0Explicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(rand1,0);
    BOOST_CHECK_THROW(mppo(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(mppoTestLen11K0, bdata::xrange<int>(-10, 10, 1), rand1) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    imat w = mppo(v, 0);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen11KN, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-3, 0, 1), rand1, rand2) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    BOOST_CHECK_THROW(mppo(v,rand2), Exception);
}
BOOST_DATA_TEST_CASE(mppoTestLen11K1, bdata::xrange<int>(-10, 10, 1), rand1) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    imat w = mppo(v, 1);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1 + 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen11K2, bdata::xrange<int>(-10, 10, 1), rand1) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    imat w = mppo(v, 2);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1 * rand1 + rand1 + 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen11K3, bdata::xrange<int>(-10, 10, 1), rand1) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    imat w = mppo(v, 3);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1 * rand1 * rand1 + rand1 * rand1 + rand1 + 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen22K0, bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1), rand1, rand2, rand3, rand4) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    imat w = mppo(v, 0);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_EQUAL(w(0,0), 1);
    BOOST_CHECK_EQUAL(w(0,1), 0);
    BOOST_CHECK_EQUAL(w(1,0), 0);
    BOOST_CHECK_EQUAL(w(1,1), 1);
}
BOOST_DATA_TEST_CASE(mppoTestLenUnequal12, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) *  bdata::xrange<int>(-3, 3, 1), rand1, rand2, rand3) {
    imat v;
    v.setlength(1,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    BOOST_CHECK_THROW(mppo(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(mppoTestLenUnequal21, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) *  bdata::xrange<int>(-3, 3, 1), rand1, rand2, rand3) {
    imat v;
    v.setlength(2,1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    BOOST_CHECK_THROW(mppo(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(mppoTestLen22KN, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 0, 1), rand1, rand2, rand3, rand4, rand5) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    BOOST_CHECK_THROW(mppo(v, rand5), Exception);
}
BOOST_DATA_TEST_CASE(mppoTestLen22K1, bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1), rand1, rand2, rand3, rand4) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    imat w = mppo(v, 1);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_EQUAL(w(0,0), rand1 + 1);
    BOOST_CHECK_EQUAL(w(0,1), rand2);
    BOOST_CHECK_EQUAL(w(1,0), rand3);
    BOOST_CHECK_EQUAL(w(1,1), rand4 + 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen22K2, bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1), rand1, rand2, rand3, rand4) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    imat w = mppo(v, 2);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_EQUAL(w(0,0), rand1 * rand1 + rand1 + rand2 * rand3 + 1);
    BOOST_CHECK_EQUAL(w(0,1), rand1 * rand2 + rand2 * rand4 + rand2);
    BOOST_CHECK_EQUAL(w(1,0), rand1 * rand3 + rand3 * rand4 + rand3);
    BOOST_CHECK_EQUAL(w(1,1), rand4 * rand4 + rand4 + rand2 * rand3 + 1);
}
BOOST_DATA_TEST_CASE(mppoTestLen22K3, bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1) * bdata::xrange<int>(-5, 5, 1), rand1, rand2, rand3, rand4) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    imat w = mppo(v, 3);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_EQUAL(w(0,0), rand1 * rand1 * rand1 + 2 * rand1 * rand2 * rand3 + rand2 * rand3 * rand4 + rand1 * rand1 + rand1 + rand2 * rand3 + 1);
    BOOST_CHECK_EQUAL(w(0,1), rand1 * rand1 * rand2 + rand1 * rand2 * rand4 + rand2 * rand2 * rand3 + rand2 * rand4 * rand4 + rand1 * rand2 + rand2 * rand4 + rand2);
    BOOST_CHECK_EQUAL(w(1,0), rand1 * rand1 * rand3 + rand2 * rand3 * rand3 + rand1 * rand3 * rand4 + rand3 * rand4 * rand4 + rand1 * rand3 + rand3 * rand4 + rand3);
    BOOST_CHECK_EQUAL(w(1,1), rand1 * rand2 * rand3 + 2 * rand2 * rand3 * rand4 + rand4 * rand4 * rand4 + rand4 * rand4 + rand4 + rand2 * rand3 + 1);
}
BOOST_AUTO_TEST_SUITE_END();
BOOST_AUTO_TEST_SUITE(Rizev);
BOOST_AUTO_TEST_CASE(RizevLen0) {
    ivect v;
    BOOST_CHECK_THROW(rizev(v), Exception);
}
BOOST_AUTO_TEST_CASE(RizevLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(rizev(v), Exception);
}
BOOST_DATA_TEST_CASE(RizevLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    vect w = rizev(v);
    BOOST_CHECK_EQUAL(w.length(), 1);
    BOOST_CHECK_CLOSE(w[0], rand1, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizevLen2, bdata::xrange<int>(-20, 20, 1) * bdata::xrange<int>(-20, 20, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    vect w = rizev(v);
    BOOST_CHECK_EQUAL(w.length(), 2);
    BOOST_CHECK_CLOSE(w[0], rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w[1], rand2, 0.0001f);
}
BOOST_DATA_TEST_CASE(Rizev3, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    vect w = rizev(v);
    BOOST_CHECK_EQUAL(w.length(), 3);
    BOOST_CHECK_CLOSE(w[0], rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w[1], rand2, 0.0001f);
    BOOST_CHECK_CLOSE(w[2], rand3, 0.0001f);
}
BOOST_AUTO_TEST_SUITE_END();
//Rizem
BOOST_AUTO_TEST_SUITE(Rizem, * utf::disabled());
BOOST_AUTO_TEST_CASE(RizemLen00) {
    imat v;
    BOOST_CHECK_THROW(rizem(v), Exception);
}
BOOST_DATA_TEST_CASE(RizemLen0KExplicit, bdata::xrange<int>(0,10,1), rand1) {
    imat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(rizem(v), Exception);
}
BOOST_DATA_TEST_CASE(RizemLenK0Explicit, bdata::xrange<int>(0,10,1), rand1) {
    imat v;
    v.setlength(rand1, 0);
    BOOST_CHECK_THROW(rizem(v), Exception);
}
BOOST_DATA_TEST_CASE(RizemTestLen11, bdata::xrange<int>(-10, 10, 1), rand1) {
    imat v;
    v.setlength(1,1);
    v(0,0) = rand1;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizemTestLen12, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    imat v;
    v.setlength(1,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), rand2, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizemTestLen21, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    imat v;
    v.setlength(2,1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), rand2, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizemTestLen13, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    imat v;
    v.setlength(1,3);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(0,2) = rand3;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 3);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), rand2, 0.0001f);
    BOOST_CHECK_CLOSE(w(0,2), rand3, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizemTestLen31, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    imat v;
    v.setlength(3,1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    v(2,0) = rand3;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 3);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), rand2, 0.0001f);
    BOOST_CHECK_CLOSE(w(2,0), rand3, 0.0001f);
}
BOOST_DATA_TEST_CASE(RizemTestLen22, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3, rand4) {
    imat v;
    v.setlength(2,2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    mat w = rizem(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_CLOSE(w(0,0), rand1, 0.0001f);
    BOOST_CHECK_CLOSE(w(0,1), rand2, 0.0001f);
    BOOST_CHECK_CLOSE(w(1,0), rand3, 0.0001f);
    BOOST_CHECK_CLOSE(w(1,1), rand4, 0.0001f);
}
BOOST_AUTO_TEST_SUITE_END();
//Identv
BOOST_AUTO_TEST_SUITE(Identv);
BOOST_DATA_TEST_CASE(IdentvAny, bdata::xrange<int>(1, 10, 1), rand1) {
    vect v = identv(rand1);
    int i = 0;
    BOOST_CHECK_EQUAL(v.length(), rand1);
    for (i = 0; i < rand1; i++)
        BOOST_CHECK_CLOSE(v[i], 1, 0.0001f);
}
BOOST_DATA_TEST_CASE(IdentvN, bdata::xrange<int>(-10, 0, 1), rand1) {
    BOOST_CHECK_THROW(identv(rand1), Exception);
}
BOOST_AUTO_TEST_SUITE_END();
//Ident
BOOST_AUTO_TEST_SUITE(Ident);
BOOST_DATA_TEST_CASE(IdentAny, bdata::xrange<int>(1, 10, 1), rand1) {
    mat v = ident(rand1);
    int i = 0, j = 0;
    BOOST_CHECK_EQUAL(v.rows(), rand1);
    BOOST_CHECK_EQUAL(v.cols(), rand1);
    for (i = 0; i < rand1; i++) {
        for (j = 0; j < rand1; j++) {
            if (i == j)
                BOOST_CHECK_CLOSE(v(i,j), 1, 0.0001f);
            else
                BOOST_CHECK_CLOSE(v(i,j), 0, 0.0001f);
        }
    }
}
BOOST_DATA_TEST_CASE(IdentN, bdata::xrange<int>(-10, 0, 1), rand1) {
    BOOST_CHECK_THROW(ident(rand1), Exception);
}
BOOST_AUTO_TEST_SUITE_END();
//vecisint
BOOST_AUTO_TEST_SUITE(Vecisint);
BOOST_AUTO_TEST_CASE(VecisintLen0) {
    vect v;
    BOOST_CHECK_THROW(vecisint(v), Exception);
}
BOOST_AUTO_TEST_CASE(VecisintLen0Explicit) {
    vect v;
    v.setlength(0);
    BOOST_CHECK_THROW(vecisint(v), Exception);
}
BOOST_DATA_TEST_CASE(VecisintLen1Large, bdata::xrange<double>(-10, 10, 0.15), rand1) {
    vect v;
    v.setlength(1);
    v[0] = rand1;
    if (std::abs(v[0] - round(v[0])) <= 0.001)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_DATA_TEST_CASE(VecisintLen1Small, bdata::xrange<double>(-0.0015, 0.0015, 0.0001)*bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    vect v;
    v.setlength(1);
    v[0] = rand1 + rand2;
    double diff = std::abs(v[0] - round(v[0]));
    if (diff <= 0.00095)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else if (diff >= 0.00105)
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_DATA_TEST_CASE(VecisintLen2LL, bdata::xrange<double>(-10, 10, 0.1) ^ bdata::xrange<double>(-10, 10, 0.1), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    double diff0 = std::abs(v[0] - round(v[0]));
    double diff1 = std::abs(v[1] - round(v[1]));
    if (diff0 <= 0.001 && diff1 <= 0.001)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_DATA_TEST_CASE(VecisintLen2SL, bdata::xrange<double>(-0.0015, 0.0015, 0.0001)*bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    double diff0 = std::abs(v[0] - round(v[0]));
    double diff1 = std::abs(v[1] - round(v[1]));
    if (diff0 <= 0.00095 && diff1 <= 0.00095)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else if (diff0 >= 0.00105 || diff1 >= 0.00105)
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_DATA_TEST_CASE(VecisintLen2LS, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<double>(-0.0015, 0.0015, 0.0001), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    double diff0 = std::abs(v[0] - round(v[0]));
    double diff1 = std::abs(v[1] - round(v[1]));
    if (diff0 <= 0.00095 && diff1 <= 0.00095)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else if (diff0 >= 0.00105 || diff1 >= 0.00105)
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_DATA_TEST_CASE(VecisintLen2SS, bdata::xrange<double>(-0.0015, 0.0015, 0.0001)*bdata::xrange<double>(-0.0015, 0.0015, 0.0001), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    double diff0 = std::abs(v[0] - round(v[0]));
    double diff1 = std::abs(v[1] - round(v[1]));
    if (diff0 <= 0.00095 && diff1 <= 0.00095)
        BOOST_CHECK_EQUAL(vecisint(v), true);
    else if (diff0 >= 0.00105 || diff1 >= 0.00105)
        BOOST_CHECK_EQUAL(vecisint(v), false);
}
BOOST_AUTO_TEST_SUITE_END();
//ivecisnonneg
BOOST_AUTO_TEST_SUITE(Ivecisnonneg);
BOOST_AUTO_TEST_CASE(IvecisnonnegLen0) {
    ivect v;
    BOOST_CHECK_THROW(ivecisnonneg(v), Exception);
}
BOOST_AUTO_TEST_CASE(IvecisnonnegLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(ivecisnonneg(v), Exception);
}
BOOST_DATA_TEST_CASE(IvecisnonnegLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    if (rand1 >= 0)
        BOOST_CHECK(ivecisnonneg(v));
    else
        BOOST_CHECK(!ivecisnonneg(v));
}
BOOST_DATA_TEST_CASE(IvecisnonnegLen2, bdata::xrange<int>(-20, 20, 1) * bdata::xrange<int>(-20, 20, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    if (rand1 >= 0 && rand2 >= 0)
        BOOST_CHECK(ivecisnonneg(v));
    else
        BOOST_CHECK(!ivecisnonneg(v));
}
BOOST_DATA_TEST_CASE(IvecisnonnegLen3, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    if (rand1 >= 0 && rand2 >= 0 && rand3 >= 0)
        BOOST_CHECK(ivecisnonneg(v));
    else
        BOOST_CHECK(!ivecisnonneg(v));
}
BOOST_AUTO_TEST_SUITE_END();
//(i)vec2matc
BOOST_AUTO_TEST_SUITE(IVec2matc);
BOOST_AUTO_TEST_CASE(IVec2matcLen0) {
    ivect v;
    BOOST_CHECK_THROW(vec2matc(v), Exception);
}
BOOST_AUTO_TEST_CASE(IVec2matcLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(vec2matc(v), Exception);
}
BOOST_DATA_TEST_CASE(IVec2matcLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    imat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
}
BOOST_DATA_TEST_CASE(IVec2matcLen2, bdata::xrange<int>(-20, 20, 1) * bdata::xrange<int>(-20, 20, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    imat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(1,0), rand2);
}
BOOST_DATA_TEST_CASE(IVec2matcLen3, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    imat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 3);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(1,0), rand2);
    BOOST_CHECK_EQUAL(w(2,0), rand3);
}
BOOST_AUTO_TEST_SUITE_END();
//(i)vec2matr
BOOST_AUTO_TEST_SUITE(IVec2matr);
BOOST_AUTO_TEST_CASE(IVec2matrLen0) {
    ivect v;
    BOOST_CHECK_THROW(vec2matr(v), Exception);
}
BOOST_AUTO_TEST_CASE(IVec2matrLen0Explicit) {
    ivect v;
    v.setlength(0);
    BOOST_CHECK_THROW(vec2matr(v), Exception);
}
BOOST_DATA_TEST_CASE(IVec2matrLen1, bdata::xrange<int>(-50, 50, 1), rand1) {
    ivect v;
    v.setlength(1);
    v[0] = rand1;
    imat w = vec2matr(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
}
BOOST_DATA_TEST_CASE(IVec2matrLen2, bdata::xrange<int>(-20, 20, 1) * bdata::xrange<int>(-20, 20, 1), rand1, rand2) {
    ivect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    imat w = vec2matr(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 2);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(0,1), rand2);
}
BOOST_DATA_TEST_CASE(IVec2matrLen3, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2, rand3) {
    ivect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    imat w = vec2matr(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 3);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(0,1), rand2);
    BOOST_CHECK_EQUAL(w(0,2), rand3);
}
BOOST_AUTO_TEST_SUITE_END();
//(r)vec2imatc
BOOST_AUTO_TEST_SUITE(RVec2matc);
BOOST_AUTO_TEST_CASE(RVec2matcLen0) {
    vect v;
    BOOST_CHECK_THROW(vec2matc(v), Exception);
}
BOOST_AUTO_TEST_CASE(RVec2matcLen0Explicit) {
    vect v;
    v.setlength(0);
    BOOST_CHECK_THROW(vec2matc(v), Exception);
}
BOOST_DATA_TEST_CASE(RVec2matcLen1, bdata::xrange<double>(-2, 2, 0.3), rand1) {
    vect v;
    v.setlength(1);
    v[0] = rand1;
    mat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 1);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
}
BOOST_DATA_TEST_CASE(RVec2matcLen2, bdata::xrange<double>(-2, 2, 0.3) * bdata::xrange<double>(-2, 2, 0.3), rand1, rand2) {
    vect v;
    v.setlength(2);
    v[0] = rand1;
    v[1] = rand2;
    mat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 2);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(1,0), rand2);
}
BOOST_DATA_TEST_CASE(RVec2matcLen3, bdata::xrange<double>(-2, 2, 0.3) * bdata::xrange<double>(-2, 2, 0.3) * bdata::xrange<double>(-2, 2, 0.3), rand1, rand2, rand3) {
    vect v;
    v.setlength(3);
    v[0] = rand1;
    v[1] = rand2;
    v[2] = rand3;
    mat w = vec2matc(v);
    BOOST_CHECK_EQUAL(w.rows(), 3);
    BOOST_CHECK_EQUAL(w.cols(), 1);
    BOOST_CHECK_EQUAL(w(0,0), rand1);
    BOOST_CHECK_EQUAL(w(1,0), rand2);
    BOOST_CHECK_EQUAL(w(2,0), rand3);
}
BOOST_AUTO_TEST_SUITE_END();
//add
BOOST_AUTO_TEST_SUITE(Add);
BOOST_AUTO_TEST_CASE(AddLen001) {
    ivect v1;
    ivect v2;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}

BOOST_AUTO_TEST_CASE(AddLen002) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}

BOOST_AUTO_TEST_CASE(AddLen003) {
    ivect v1;
    ivect v2;
    v2.setlength(0);
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}

BOOST_AUTO_TEST_CASE(AddLen004) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    v2.setlength(0);
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen01, bdata::xrange<int>(-10, 10, 1), rand1) {
    ivect v1;
    ivect v2;
    v2.setlength(1);
    v2[0] = rand1;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen02, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    ivect v1;
    ivect v2;
    v2.setlength(1);
    v2[0] = rand1;
    v2[1] = rand2;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen10, bdata::xrange<int>(-10, 10, 1), rand1) {
    ivect v1;
    ivect v2;
    v2.setlength(1);
    v2[0] = rand1;
    BOOST_CHECK_THROW(add(v2, v1), Exception);
}
BOOST_DATA_TEST_CASE(AddLen20, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    ivect v1;
    ivect v2;
    v2.setlength(1);
    v2[0] = rand1;
    v2[1] = rand2;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen01E, bdata::xrange<int>(-10, 10, 1), rand1) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    v2.setlength(1);
    v2[0] = rand1;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen02E, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    v2.setlength(1);
    v2[0] = rand1;
    v2[1] = rand2;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLen10E, bdata::xrange<int>(-10, 10, 1), rand1) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    v2.setlength(1);
    v2[0] = rand1;
    BOOST_CHECK_THROW(add(v2, v1), Exception);
}
BOOST_DATA_TEST_CASE(AddLen20E, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    ivect v1;
    ivect v2;
    v1.setlength(0);
    v2.setlength(1);
    v2[0] = rand1;
    v2[1] = rand2;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLenIncomp12, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1), rand1, rand2, rand3) {
    ivect v1;
    ivect v2;
    v1.setlength(1);
    v2.setlength(2);
    v1[0] = rand1;
    v2[0] = rand2;
    v2[1] = rand3;
    BOOST_CHECK_THROW(add(v1, v2), Exception);
}
BOOST_DATA_TEST_CASE(AddLenIncomp21, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1), rand1, rand2, rand3) {
    ivect v1;
    ivect v2;
    v1.setlength(1);
    v2.setlength(2);
    v1[0] = rand1;
    v2[0] = rand2;
    v2[1] = rand3;
    BOOST_CHECK_THROW(add(v2, v1), Exception);
}
BOOST_DATA_TEST_CASE(AddLen11, bdata::xrange<int>(-10, 10, 1) * bdata::xrange<int>(-10, 10, 1), rand1, rand2) {
    ivect v1;
    ivect v2;
    v1.setlength(1);
    v2.setlength(1);
    v1[0] = rand1;
    v2[0] = rand2;
    ivect res = add(v1, v2);
    BOOST_CHECK_EQUAL(res[0], rand1 + rand2);
}
BOOST_DATA_TEST_CASE(AddLen22, bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1) * bdata::xrange<int>(-3, 3, 1), rand1, rand2, rand3, rand4) {
    ivect v1;
    ivect v2;
    v1.setlength(2);
    v2.setlength(2);
    v1[0] = rand1;
    v1[1] = rand2;
    v2[0] = rand3;
    v2[1] = rand4;
    ivect res = add(v1, v2);
    BOOST_CHECK_EQUAL(res[0], rand1 + rand3);
    BOOST_CHECK_EQUAL(res[1], rand2 + rand4);
}
BOOST_AUTO_TEST_SUITE_END();
//matc2vec
BOOST_AUTO_TEST_SUITE(Matc2vec);
BOOST_DATA_TEST_CASE(Matc2vec00, bdata::xrange<int>(-10,10,1), rand1) {
    imat v;
    BOOST_CHECK_THROW(matc2vec(v, rand1), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vec0KExplicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vecK0Explicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(rand1, 0);
    BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vec11, bdata::xrange<int>(-10,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(1, 1);
    v(0,0) = rand1;
    if (!rand2) {
        ivect w = matc2vec(v, rand2);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vec21, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    imat v;
    v.setlength(2, 1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    if (!rand3) {
        ivect w = matc2vec(v, rand3);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand2);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vec12, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    imat v;
    v.setlength(1, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    if (!rand3) {
        ivect w = matc2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else if (rand3 == 1) {
        ivect w = matc2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand2);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(Matc2vec22, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3, rand4, rand5)
{
    imat v;
    v.setlength(2, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    if (!rand5) {
        ivect w = matc2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand3);
    }
    else if (rand5 == 1) {
        ivect w = matc2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand2);
        BOOST_CHECK_EQUAL(w[1], rand4);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand5), Exception);
}
BOOST_AUTO_TEST_SUITE_END();
//matr2vec
BOOST_AUTO_TEST_SUITE(Matr2vec);
BOOST_DATA_TEST_CASE(Matr2vec00, bdata::xrange<int>(-10,10,1), rand1) {
    imat v;
    BOOST_CHECK_THROW(matr2vec(v, rand1), Exception);
}
BOOST_DATA_TEST_CASE(Matr2vec0KExplicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(matr2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matr2vecK0Explicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(rand1, 0);
    BOOST_CHECK_THROW(matr2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matr2vec11, bdata::xrange<int>(-10,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    imat v;
    v.setlength(1, 1);
    v(0,0) = rand1;
    if (!rand2) {
        ivect w = matr2vec(v, rand2);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else
        BOOST_CHECK_THROW(matr2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(Matr2ver12, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    imat v;
    v.setlength(1, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    if (!rand3) {
        ivect w = matr2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand2);
    }
    else
        BOOST_CHECK_THROW(matr2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(Matr2vec21, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    imat v;
    v.setlength(2, 1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    if (!rand3) {
        ivect w = matr2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else if (rand3 == 1) {
        ivect w = matr2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand2);
    }
    else
        BOOST_CHECK_THROW(matr2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(Matr2vec22, bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-1,1,1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3, rand4, rand5)
{
    imat v;
    v.setlength(2, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    if (!rand5) {
        ivect w = matr2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand2);
    }
    else if (rand5 == 1) {
        ivect w = matr2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand3);
        BOOST_CHECK_EQUAL(w[1], rand4);
    }
    else
        BOOST_CHECK_THROW(matr2vec(v, rand5), Exception);
}
BOOST_AUTO_TEST_SUITE_END();
//(r)matc2vec
BOOST_AUTO_TEST_SUITE(RMatc2vec);
BOOST_DATA_TEST_CASE(RMatc2vec00, bdata::xrange<int>(-10,10,1), rand1) {
    mat v;
    BOOST_CHECK_THROW(matc2vec(v, rand1), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vec0KExplicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    mat v;
    v.setlength(0, rand1);
    BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vecK0Explicit, bdata::xrange<int>(0,10,1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    mat v;
    v.setlength(rand1, 0);
    BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vec11, bdata::xrange<double>(-1,1,0.1) * bdata::xrange<int>(-10,10,1), rand1, rand2) {
    mat v;
    v.setlength(1, 1);
    v(0,0) = rand1;
    if (!rand2) {
        vect w = matc2vec(v, rand2);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand2), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vec21, bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    mat v;
    v.setlength(2, 1);
    v(0,0) = rand1;
    v(1,0) = rand2;
    if (!rand3) {
        vect w = matc2vec(v, rand3);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand2);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vec12, bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<int>(-10,10,1), rand1, rand2, rand3)
{
    mat v;
    v.setlength(1, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    if (!rand3) {
        vect w = matc2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand1);
    }
    else if (rand3 == 1) {
        vect w = matc2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 1);
        BOOST_CHECK_EQUAL(w[0], rand2);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand3), Exception);
}
BOOST_DATA_TEST_CASE(RMatc2vec22, bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<double>(-0.2,0.2,0.1) * bdata::xrange<int>(-5,5,1), rand1, rand2, rand3, rand4, rand5)
{
    mat v;
    v.setlength(2, 2);
    v(0,0) = rand1;
    v(0,1) = rand2;
    v(1,0) = rand3;
    v(1,1) = rand4;
    if (!rand5) {
        vect w = matc2vec(v, 0);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand1);
        BOOST_CHECK_EQUAL(w[1], rand3);
    }
    else if (rand5 == 1) {
        vect w = matc2vec(v, 1);
        BOOST_CHECK_EQUAL(w.length(), 2);
        BOOST_CHECK_EQUAL(w[0], rand2);
        BOOST_CHECK_EQUAL(w[1], rand4);
    }
    else
        BOOST_CHECK_THROW(matc2vec(v, rand5), Exception);
}
BOOST_AUTO_TEST_SUITE_END();
/*BOOST_AUTO_TEST_CASE(MyTestCaseFalse)
{
    // To simplify this example test, let's suppose we'll test 'float'.
    // Some test are stupid, but all should pass.
    float x = 9.5f;
    
    BOOST_CHECK(x != 4.5f);
    BOOST_CHECK_EQUAL((int)x, 9);
    BOOST_CHECK_CLOSE(x, 9.5f, 0.0001f); // Checks differ no more then 0.0001%
}*/
