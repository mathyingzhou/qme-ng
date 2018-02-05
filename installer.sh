#!/bin/bash
export QME_NG_PATH = $PWD
cd ..
curl -o gmp.tar.bz2 -k https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2
bzip2 -d gmp.tar.bz2
tar xopf gmp.tar
cd gmp-6.1.2
sudo ./configure —-enable-cxx
sudo make install
sudo make check
cd ..
curl —L -o boost.tar.bz2 -k https://dl.bintray.com/boostorg/release/1.66.0/source/boost_1_66_0.tar.bz2
bzip2 -d boost.tar.bz2
tar xopf boost.tar
cp -a /boost_1_66_0/. /usr/local/boost
cd boost_1_66_0
sudo ./bootstrap.sh —-prefix=/usr/local/boost
sudo ./b2 install
install_name_tool -change libboost_program_options.dylib /usr/local/lib/libboost_program_options.dylib qme-ng
cd $QME_NG_PATH
make all
