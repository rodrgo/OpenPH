#!/bin/bash

LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab -nodesktop -nosplash -nodisplay -r "benchmark_phat;quit;"
