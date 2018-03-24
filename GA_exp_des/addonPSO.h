#ifndef addonPSO_h
#define addonPSO_h

#include <algorithm>			//min, max
//#include <cmath>				//ceil
#include <fstream>				// ">>", "<<", is_open(), open, close, ignore, seekg
#include <iostream>				//cin, cout
#include <iomanip>
//#include <istream>
//#include <memory>				//allocate
#include <math.h> 				//sqrt
#include <new>					//new
#include <numeric>				//accumulate
#include <stdio.h>				//getchar
#include <string>				//string, .c_str()
//#include <time.h>
#include <vector>				//std::vector<size>

//SVD (LAPACK)
extern "C" void DGESVD(const char *jobu, const char *jobvt, const int *M, const int *N, 
double *A, const int *lda, double *S, double *U, const int *ldu,
double *VT, const int *ldvt, double *work, const int *lwork, int *info);
//General Matrix Mult (BLAS)
extern "C" void DGEMM(const char *transA, const char *transB, const int *M, const int *N, const int *K,
double *alpha, double *A, const int *lda, double *B, const int *ldb,
double *beta, double *C, const int *ldc);

#endif