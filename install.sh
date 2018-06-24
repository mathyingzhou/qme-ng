#!/bin/bash
echo "qme-ng is a quiver mutation software originally written by Mattieu Perotin and improved by Ying Zhou. The installer may take about 60-75 minutes to finish because it needs to install all dependencies."
echo "Now we need to install the Command Line Tools for macOS if it hasn't already been installed. If that's the case you will see a popup window. Please click "Install" and wait for the tools to be downloaded and installed. Please refrain from pressing enter before Command Line Tools get installed."
echo "On the other hand if you already have Command Line Tools then you will just see an error message. After that you just need to press enter to continue."
xcode-select --install
read -p "When command line tools have been successfully installed, please press enter to continue."
curl -o qme-ng.zip -L -k  https://github.com/mathyingzhou/qme-ng/archive/master.zip
unzip qme-ng.zip
cd qme-ng-master
export QME_NG_PATH=$PWD
cd ..
curl -o gmp.tar.bz2 -k https://ftp.gnu.org/gnu/gmp/gmp-6.1.2.tar.bz2
bzip2 -d gmp.tar.bz2
tar xopf gmp.tar
cd gmp-6.1.2
bash ./configure --enable-cxx
make install
make check
cd ..
curl -o boost.tar.bz2 -L -k https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.tar.bz2
bzip2 -d boost.tar.bz2
tar xopf boost.tar
cp -a /boost_1_67_0/. /usr/local/boost
cd boost_1_67_0
bash ./bootstrap.sh —-prefix=/usr/local/boost
./b2 install
cd $QME_NG_PATH
make all
echo "Congratulations! You have successfully installed qme-ng! Please email Ying Zhou at yzhou935@brandeis.edu if you have any questions concerning the software. Have fun exploring quivers!"
