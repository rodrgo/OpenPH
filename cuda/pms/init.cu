
/*
** This is a single header file for all algorithms.
** The order of the #include is important based on dependencies 
** some of the files.
*/


//#define SAFE


#include <iostream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <cuda.h>
#include <curand.h>
#include <mex.h>
#include <thrust/sort.h>
#include <thrust/count.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/generate.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include "cublas.h"

#include <math.h>
#include "kernels.cu"
#include "functions.cu"
#include "debug_tools.cu"

#include "ph/kernels.cu"
#include "ph/dimensions.cu"
#include "ph/beta.cu"
#include "ph/pms_functions.cu"
#include "ph/standard.cu"
#include "ph/twist.cu"
#include "ph/ph_row.cu"
#include "ph/pms.cu"


