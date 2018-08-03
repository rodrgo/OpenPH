#!/bin/bash

export LD_PRELOAD=$LD_PRELOAD:/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/usr/lib/x86_64-linux-gnu/libprotobuf.so.9
cat /home/rodrigo/workspace/persistent_homology/tda/cuda/pms/test_ph.m | matlab -nodesktop -nosplash

