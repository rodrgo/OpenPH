#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Outputs some information on CUDA-enabled devices on your computer,
including current memory usage.

It's a port of https://gist.github.com/f0k/0d6431e3faa60bffc788f8b4daa029b1
from C to Python with ctypes, so it can run without compiling anything. Note
that this is a direct translation with no attempt to make the code Pythonic.
It's meant as a general demonstration on how to obtain CUDA device information
from Python without resorting to nvidia-smi or a compiled Python extension.

Author: Jan Schlüter

# Modified by RMS
"""

import sys
import ctypes


# Some constants taken from cuda.h
CUDA_SUCCESS = 0
CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16
CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR = 39
CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13
CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE = 36


def ConvertSMVer2Cores(major, minor):
    # Returns the number of CUDA cores per multiprocessor for a given
    # Compute Capability version. There is no way to retrieve that via
    # the API, so it needs to be hard-coded.
    return {(1, 0): 8,
            (1, 1): 8,
            (1, 2): 8,
            (1, 3): 8,
            (2, 0): 32,
            (2, 1): 48,
            }.get((major, minor), 192)  # 3.0 and above


def main():
    libnames = ('libcuda.so', 'libcuda.dylib', 'cuda.dll')
    for libname in libnames:
        try:
            cuda = ctypes.CDLL(libname)
        except OSError:
            continue
        else:
            break
    else:
        raise OSError("could not load any of: " + ' '.join(libnames))

    nGpus = ctypes.c_int()
    name = b' ' * 100
    cc_major = ctypes.c_int()
    cc_minor = ctypes.c_int()
    cores = ctypes.c_int()
    threads_per_core = ctypes.c_int()
    clockrate = ctypes.c_int()
    freeMem = ctypes.c_size_t()
    totalMem = ctypes.c_size_t()

    result = ctypes.c_int()
    device = ctypes.c_int()
    context = ctypes.c_void_p()
    error_str = ctypes.c_char_p()

    result = cuda.cuInit(0)
    if result != CUDA_SUCCESS:
        cuda.cuGetErrorString(result, ctypes.byref(error_str))
        print("cuInit failed with error code %d: %s" % (result, error_str.value.decode()))
        return 1
    result = cuda.cuDeviceGetCount(ctypes.byref(nGpus))
    if result != CUDA_SUCCESS:
        cuda.cuGetErrorString(result, ctypes.byref(error_str))
        print("cuDeviceGetCount failed with error code %d: %s" % (result, error_str.value.decode()))
        return 1
    for i in range(nGpus.value):
        result = cuda.cuDeviceGet(ctypes.byref(device), i)
        if result != CUDA_SUCCESS:
            cuda.cuGetErrorString(result, ctypes.byref(error_str))
            print("cuDeviceGet failed with error code %d: %s" % (result, error_str.value.decode()))
            return 1
        if cuda.cuDeviceComputeCapability(ctypes.byref(cc_major), ctypes.byref(cc_minor), device) == CUDA_SUCCESS:
            print("sm_%d%d" % (cc_major.value, cc_minor.value))
    return 0


if __name__=="__main__":
    sys.exit(main())
