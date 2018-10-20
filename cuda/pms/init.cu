
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
#include "dimensions.cu"
#include "tic_toc.cu"

#include "algorithms/kernels.cu"
#include "algorithms/tracker.cu"
#include "algorithms/standard/standard_parallel.cu"
#include "algorithms/twist/twist_parallel.cu"
#include "algorithms/ph_row/ph_row_parallel.cu"
#include "algorithms/pms/beta.cu"
#include "algorithms/pms/pms_functions.cu"
#include "algorithms/pms/pms.cu"

#include "functions_host.cu"
#include "algorithms/kernels_host.cu"
#include "algorithms/standard/standard.cu"
#include "algorithms/twist/twist.cu"

#include "algorithms/algorithm_factory.cu"

