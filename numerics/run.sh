#!/bin/bash

matlab -nodesktop -nosplash -nodisplay -r "benchmark_pms;quit;"
matlab -nodesktop -nosplash -nodisplay -r "benchmark_roc;quit;"
LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab -nodesktop -nosplash -nodisplay -r "benchmark_phat;quit;"
LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6" matlab -nodesktop -nosplash -nodisplay -r "benchmark_datasets;quit;"

