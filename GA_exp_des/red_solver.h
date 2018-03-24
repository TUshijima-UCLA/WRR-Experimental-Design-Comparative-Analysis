#ifndef red_solver_h
#define red_solver_h

#include <algorithm>			//min, max
#include <fstream>				// ">>", "<<", is_open(), open, close, ignore, seekg
#include <iostream>				//cin, cout
#include <iomanip>
//#include <istream>
//#include <memory>				//allocate
#include <math.h> 				//sqrt, ceil
#include <new>					//new
//#include <numeric>				//accumulate
//#include <stdio.h>				//getchar
//#include <string>				//string, .c_str()
#include <time.h>
#include <vector>				//std::vector<size>

void read_red_solver(std::vector<double> &P, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<double> &well, std::vector<double> &wellrate, int &npc, int &Nn, int &nq, double &dt, double &T);

//LAPACK symmetric solver
extern "C" void dsysv_(const char *uplo, const int *N, const int *nrhs, double *A,
const int *lda, int *ipiv, double *B, const int *ldb, double *work, int *lwork, int *info);
/*
//LAPACK general solver
extern "C" void dgels_(const char *trans, const int *M, const int *N, const int *nrhs,
double *A, const int *lda, double *B, const int *ldb, double *work,
const int *lwork, int *info);
*/

//BLAS general double matrix vector multiplication
extern "C" void dgemv_(const char *trans, int *M, int *N, double *alpha, double *A,
int *lda, double *X, int *incx, double *beta, double *Y, int *incy);

/*
//BLAS general double matrix matrix multiplication (Declared in matrix.h)
extern "C" void dgemm_(const char *transA, const char *transB, int *M, int *N,
int *K, double *alpha, double *A, int *lda, double *B, int *ldb, double *beta,
double *C, int *ldC);
*/
extern std::vector<double> tmes;
extern std::vector<int> optowell;
extern std::vector<double> P;
extern std::vector<double> Ar;
extern std::vector<double> Br;
extern std::vector<double> well;
extern std::vector<double> wellrate;
extern std::vector<double> wellweight;
extern std::vector<int> owellindex;
extern std::vector<int> owellloc;
extern std::vector<int> owellzone;
extern std::vector<double> owellx;
extern std::vector<double> owelly;
extern std::vector<double> I;
//extern std::vector<double> r;
extern std::vector<double> Jdt;
extern int npc;
extern int Nn;
extern int nq;
extern double dt;
extern double T;
extern int ntmes;
extern int nobswell;
extern int ncalls;
extern double trI;
extern double convg;

#endif