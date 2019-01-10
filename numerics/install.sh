#!/usr/bin/env bash

# --------
# Javaplex
# --------

git clone --depth=1 https://github.com/rodrgo/javaplex dependencies/javaplex
ant -buildfile dependencies/javaplex/build.xml

# --------
# PHAT 
# --------

git clone https://bitbucket.org/phat-code/phat.git dependencies/phat
cd dependencies/phat
git reset --hard e9800d103fcdd19f0417e89781761f7f85d8ec9b 
cmake .
make
cd ../..

# --------
# DIPHA
# --------

git clone https://github.com/DIPHA/dipha.git dependencies/dipha
cd dependencies/dipha
git reset --hard 0b874769fbd092c07a12cebc2459adb02117c2fd
cmake .
make
cd ../..

# --------
# Datasets
# --------

for i in dependencies/dipha/test_data/*.complex; do
    echo "$i"
    ./dependencies/dipha/create_phat_filtration "${i}" datasets/tmp.f
    ./dependencies/phat/convert --save-ascii datasets/tmp.f "datasets/$(basename "$i").dat"
    rm datasets/tmp.f
done


