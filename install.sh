#!/bin/bash
export QME_NG_PATH=$PWD
cd ..
curl -o gmp.tar.bz2 -k https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2
bzip2 -d gmp.tar.bz2
tar xopf gmp.tar
cd gmp-6.1.2
bash -x ./configure --enable-cxx
make install
make check
cd ..
curl -o boost.tar.bz2 -L -k https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2
bzip2 -d boost.tar.bz2
tar xopf boost.tar
cp -a /boost_1_67_0/. /usr/local/boost
cd boost_1_67_0
bash -x ./bootstrap.sh —-prefix=/usr/local/boost
./b2 install
cd $QME_NG_PATH
make all
