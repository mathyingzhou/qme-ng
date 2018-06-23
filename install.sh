#!/bin/bash
echo -e "qme-ng is a quiver mutation software originally written by Mattieu Perotin and improved by Ying Zhou. The installer may takes about 30-45 minutes to finish. Please do not press any key unless you are asked to do so.\n"
echo -e "Now we need to install the Command Line Tools for macOS. You will see a popup window. Please click "install" and wait for the tools to be downloaded and installed.\n"
xcode-select --install
read -p "Please refrain from pressing enter before command line tools get installed. When command line tools have been successfully installed, please press enter to continue."
curl -o qme-ng.zip -L -k  https://github.com/mathyingzhou/qme-ng/archive/master.zip
unzip qme-ng.zip
cd qme-ng-master
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
echo -e "Congratulations! You have successfully installed qme-ng! Please email Ying Zhou at yzhou935@brandeis.edu if you have any questions concerning the software.\n"
