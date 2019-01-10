#!/usr/bin/env bash

for i in dependencies/dipha/test_data/*.complex; do
    echo "$i"
    ./dependencies/dipha/create_phat_filtration "${i}" datasets/tmp.f
    ./dependencies/phat/convert --save-ascii datasets/tmp.f "datasets/$(basename "$i").dat"
    rm datasets/tmp.f
done

