#!/bin/bash

cat benchmark_datasets.m | LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab -nodesktop -nosplash

