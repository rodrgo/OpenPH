#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Adapted from Jan Schlüter's code in
https://gist.github.com/f0k/0d6431e3faa60bffc788f8b4daa029b1
"""

import sys
import ctypes

def main():

    # Some constants taken from cuda.h
    CUDA_SUCCESS = 0
    CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT = 16
    CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_MULTIPROCESSOR = 39
    CU_DEVICE_ATTRIBUTE_CLOCK_RATE = 13
    CU_DEVICE_ATTRIBUTE_MEMORY_CLOCK_RATE = 36

    # Find libnames
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

    result = ctypes.c_int()
    device = ctypes.c_int()
    context = ctypes.c_void_p()
    error_str = ctypes.c_char_p()

    result = cuda.cuInit(0)
    if result != CUDA_SUCCESS:
        cuda.cuGetErrorString(result, ctypes.byref(error_str))
        print("cuInit failed with error code %d: %s" % (result, error_str.value.decode()), file=sys.stderr)
        return 1
    result = cuda.cuDeviceGetCount(ctypes.byref(nGpus))
    if result != CUDA_SUCCESS:
        cuda.cuGetErrorString(result, ctypes.byref(error_str))
        print("cuDeviceGetCount failed with error code %d: %s" % (result, error_str.value.decode()), file=sys.stderr)
        return 1

    sm_codes = {}
    for i in range(nGpus.value):
        result = cuda.cuDeviceGet(ctypes.byref(device), i)
        if result != CUDA_SUCCESS:
            cuda.cuGetErrorString(result, ctypes.byref(error_str))
            print("cuDeviceGet failed with error code %d: %s" % (result, error_str.value.decode()), file=sys.stderr)
            return 1
        if cuda.cuDeviceComputeCapability(ctypes.byref(cc_major), ctypes.byref(cc_minor), device) == CUDA_SUCCESS:
            sm_codes[i] = "sm_%d%d" % (cc_major.value, cc_minor.value)

    # Get sm_code and print
    sm_code = ''
    if len(list(set(sm_codes.values()))) > 1:
        print("More than one SM found. Using GPU-0", file=sys.stderr)
        sm_code = sm_codes[0]
    print(sm_code)
    return 0

if __name__=="__main__":
    sys.exit(main())
