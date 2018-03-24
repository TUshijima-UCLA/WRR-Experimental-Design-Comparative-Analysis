#include "matrix.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <string.h>

//General Matrix Vector Multiplication (BLAS)
extern "C" {void dgemv_(const char *trans, const int *M, const int *N, double *alpha, double *A,
const int *lda, double *x, const int *incx, double *beta, double *y, const int *incy);}
//General Matrix Multiply (BLAS)
extern "C" {void dgemm_(const char *transA, const char *transB, const int *M, const int *N, const int *K,
double *alpha, double *A, const int *lda, double *B, const int *ldb,
double *beta, double *C, const int *ldc);}
//SVD (LAPACK)
extern "C" {void dgesvd_(const char *jobu, const char *jobvt, const int *M, const int *N, 
double *A, const int *lda, double *S, double *U, const int *ldu,
double *VT, const int *ldvt, double *work, const int *lwork, int *info);}
//Forms R part of QR Decomp (LAPACK)
extern "C" {void dgeqrf_(int *M, int *N, double *A, int *lda, double *tau, double *work,
int *lwork, int *info);}
//Forms Q part of QR Decomp (LAPACK)
extern "C" {void dorgqr_(int *M, int *N, int *K, double *A, int *lda, double *tau,
double *work, int *lwork, int *info);}
//Norm of a matrix (LAPACK)
extern "C" {double dlange_(char *norm, const int *M, const int *N, double *A, const int *lda, double*work);}
//Eigenvalue Decomposition (LAPACK)
extern "C" {void dsyev_(const char *jobz, const char *uplo, const int *N, double *A,
const int *lda, double *W, double *work, int *lwork, int *info);}

//indexed like normal with C++ (index 0 is the first element)
mat2d::mat2d(int m, int n)
{
r1 = m;
c1 = n;
A.resize(r1*c1,0);
}

mat2d::~mat2d() {A.clear();}
//Get number of rows in matrix (int m = A.rows())
int mat2d::rows() const {return (r1);}
//Get number columns in matrix (int n = A.columns())
int mat2d::columns() const {return (c1);}
//Resize the matrix (A.resize(m,n))
void mat2d::resize(int m, int n)
	{	
		r1 = m;
		c1 = n;
		A.resize(r1*c1);
	}
//Set element m,n to value (A.set_value(m,n,value)) //store by column (first m values are the first column, FORTRAN style)
void mat2d::set_value(int m, int n, double value)
	{
		A[(n)*r1 + (m)] = value;
	}
//Get the value of element m,n (double value = A.get_value(m,n))
double mat2d::get_value(int m, int n) const
	{
		double a;
		a = A[(n)*r1 + (m)];
		return (a);
	}
//Turn full or part of matrix into a vector (std::vector<double> vec = A.get_vec(m1,m2,n1,n2) m2 > =m1, n2 >= n1)
std::vector<double> mat2d::get_vec(int m1, int m2, int n1, int n2) const
	{//1 is first index 2 is second index 
	//indexing brackets selection (if want A(:,1) m1 = 0, m2 = r1, n1 = 0, n2 = 1)
	//if take all elements in the matrix (0,r1,0,c1) vector is in correct LAPACK form
		if (m2 > m1 && n2 > n1 && m2 <= r1 && n2 <= c1)
		{
			std::vector<double> a((m2-m1)*(n2-n1),0);
			int kk = 0;
			for (int ii = n1; ii < n2; ii++)
				{
					for (int jj = m1; jj < m2; jj++)
						{
							a[kk] = A[ii*r1 + jj];
							kk++;
						}
				}
			return (a);
		}
		else
		{
			std::cout << "Vector cannot be returned for sizing reasons";
			std::exit(1);
		}
	}
	
//Print out matrix to screen (A.print())
void mat2d::print() const
	{
		if (r1 == 0 && c1 == 0)
		{
			std::cout << "Matrix is empty" << std::endl;
		}
		else
		{
			for (int ii = 0; ii < r1; ii++)
				{
					for (int jj = 0; jj < c1; jj++)
						{
							std::cout << A[jj*r1 + ii] << '\t';
						}	std::cout << std::endl;
				}
			std::cout << std::endl;
		}
	}
	
//Set one matrix equal to another (A.equals(B) -> A = B)
void mat2d::equals(mat2d &b) //copies the matrix given
	{	
		r1 = b.rows();
		c1 = b.columns();
		A.resize(r1*c1);
		for (int ii = 0; ii < c1; ii++)
			{
				for (int jj = 0; jj < r1; jj++)
					{
						A[ii*r1 + jj] = b.get_value(jj,ii);
					}
			}
	}
	
//Add two matrices (C.add(A,B) -> C = A + B)
void mat2d::add(mat2d &a, mat2d &b)
	{
		r2 = a.rows();
		c2 = a.columns();
		r3 = b.rows();
		c3 = b.columns();
		if (r2 == r3 && c2 == c3)
			{	
				r1 = r2;
				c1 = c2;
				A.resize(r1*c1);
				for (int ii = 0; ii < c1; ii++)
					{
						for (int jj = 0; jj < r1; jj++)
							{
								A[ii*r1 + jj] = a.get_value(jj,ii) + b.get_value(jj,ii);
							}
					}
			}
		else std::cout << "Cannot add two matrices of differing sizes" << std::endl;
	}
//Subtract two matrices (C.subtract(A,B) -> C = A-B)
void mat2d::subtract(mat2d &a, mat2d &b)
	{
		r2 = a.rows();
		c2 = a.columns();
		r3 = b.rows();
		c3 = b.columns();
		if (r2 == r3 && c2 == c3)
			{	
				r1 = r2;
				c1 = c2;
				A.resize(r1*c1);
				for (int ii = 0; ii < c1; ii++)
					{
						for (int jj = 0; jj < r1; jj++)
							{
								A[ii*r1 + jj] = a.get_value(jj,ii) - b.get_value(jj,ii);
							}
					}
			}
		else std::cout << "Cannot subtract two matrices of differing sizes" << std::endl;
	}
//Transpose matrix (mat2d B = A.transpose -> B = A')
mat2d mat2d::transpose()
	{
		mat2d a(c1,r1);
		for (int ii = 0; ii < c1; ii++)
			{
				for (int jj = 0; jj < r1; jj++)
					{
						a.set_value(ii,jj,A[ii*r1+jj]);
					}
			}
		return (a);
	}
	
//Clears matrix (A.clear() -> dim(A) = 1x1, A(0,0) = 0)
void mat2d::clear()
	{
		r1 = 0;
		c1 = 0;
		A.clear();
		//A.resize(r1*c1);
		//A[0] = 0.0;
	}

//Multiply a matrix by a vector (b = A.matvec_mult(alpha,'T',x,beta,y) -> b = alpha*A'*x + beta*y)
std::vector<double> mat2d::matvec_mult(double alpha, char ta, std::vector<double> &x,
double beta, std::vector<double> y) const
{
	char trans;
	std::vector<double> b(0,0);
	if (ta == 'T' || ta == 't') 
	{
		trans = 'T';
		b.resize(c1);
	}
	else if (ta == 'N' || ta == 'n') 
	{
		trans = 'N';
		b.resize(r1);
	}
	else std::cout << "Unknown matrix operation in matrix vecor multiplication" << std::endl;
//This function copies alpha, and beta does not change them on exit
//			int M = a.rows();
//			int N = a.columns();
	int incx = 1;
	int incy = 1;
	
	if (y.size() != b.size())
	{
		//will take a y of non-comforming size and copy the first b.size() elements
		y.resize(b.size());
	}
	double *b2 = &b[0];
	double *y2 = &y[0];
	//double *y2 = new double[y.size()];
	//for (int ii = 0; ii < y.size(); ii++) {y2[ii] = y[ii];}
	double *x2 = &x[0];
	//double *x2 = new double[x.size()];
	//for (int ii = 0; ii < x.size(); ii++) {x2[ii] = x[ii];}
	std::vector<double> A22 = A;
	double *A2 = &A22[0];
	//double *A2 = new double[r1*c1];
	//for (int ii = 0; ii < r1*c1; ii++) {A2[ii] = A[ii];}
	//std::cout << x.size() << '\t' << r1 << std::endl;
	if ((ta == 'T' && x.size() == r1) || (ta == 'N' && x.size() == c1))
	//if (trans == 'T' && x.size() == r1)
	{
		//std::cout << ta << std::endl;
		//::dgemv_(&trans, &r1, &c1, &alpha, A2, &r1, x2, &incx, &beta, y2, &incy);
		::dgemv_(&trans, &r1, &c1, &alpha, A2, &r1, x2, &incx, &beta, b2, &incy);
	}
//			else if (ta = ='N' && x.size() == c1) {::dgemv_(&trans, &r1, &c1, &alpha, A2, &r1, x2, &incx, &beta, y2, &incy);}
	else std::cout << "Cannot multiply non-conforming vector to a matrix" << std::endl;
	
	//std::vector<double> b(y.size(),0);
	//std::cout << b.size() << '\t' << r1 << std::endl;
	//for (int ii = 0; ii < b.size(); ii++) {b[ii] = y2[ii];}
	
	return (b);
	/*
	delete[] y2;
	delete[] x2;
	delete[] A2;
	*/
	
}

	//Multiply two matrices (D.multiply('T',A,'N',B,C,alpha,beta) -> D = alpha*A' * B + beta*C)
	//Be careful what you put in for C, its size can change on exit
	void mat2d::multiply(char ta, mat2d &a, char tb, mat2d &b, mat2d &c, double alpha, double beta)
{
	
	r2 = a.rows();
	c2 = a.columns();
	r3 = b.rows();
	c3 = b.columns();

	char transA, transB;
	int M, N, K, lda, ldb, ldc;
	if (ta == 'T')  {transA = 'T'; M = c2; K = r2; lda = K;}
	else {transA ='N'; M = r2; K = c2; lda = M;}
	if (tb == 'T') {transB = 'T'; N = r2; ldb = N;} 
	else {transB = 'N'; N = c2; ldb = K;}
	ldc = M;
	std::vector<double> av2(r2*c2);
	double *av = &av2[0];
	//av = new double[r2*c2];
	std::vector<double> bv2(r3*c3);
	double *bv = &bv2[0];
	//bv = new double[r3*c3];
	std::vector<double> cv2(ldc*M);
	double *cv = &cv2[0];
	//cv = new double[ldc*M];
	//std::vector<double> av2;
	//std::vector<double> bv2;
	if (beta == 0.0) c.resize(ldc,N);
	int r4 = c.rows();
	int c4 = c.columns();
	//check if dimensions match
	if ((ta == 'N' && tb == 'N' && c2 == r3 && r2 == r4 && c3 == c4) ||
	(ta == 'N' && tb == 'T' && c2 == c3 && r2 == r4 && r3 == c4) ||
	(ta == 'T' && tb == 'N' && r2 == r3 && c2 == r4 && c3 == c4) ||
	(ta == 'N' && tb == 'T' && c2 == c3 && r2 == r4 && r3 == c4))
	{
		av2 = a.get_vec(0,r2,0,c2);
		bv2 = b.get_vec(0,r3,0,c3);
		for (unsigned int ii = 0; ii < av2.size(); ii++) av[ii] = av2[ii];
		for (unsigned int ii = 0; ii < bv2.size(); ii++) bv[ii] = bv2[ii];
		//use :: to call global function in member function section
		::dgemm_(&transA, &transB, &M, &N, &K, &alpha, av, &lda, bv, &ldb, &beta, cv, &ldc);
		r1 = ldc;
		c1 = N;
		A.resize(r1*c1);
		//for (int ii = 0; ii < r1*c1; ii++) A[ii] = cv[ii];
		A = cv2;
	}
	else std::cout << "Cannot multiply, matrices dimensions do not match (also check dimensions of C)" << std::endl;
	/*
	delete[] av;
	delete[] bv;
	delete[] cv;
	*/
}
	//A.svd(U,S,Vt) -> A = USV' (doesn't calculate Vt)
	void mat2d::svd(mat2d &U, mat2d &S, mat2d &Vt) const
{
	U.clear();
	S.clear();
	Vt.clear();
	std::vector<double> A2 = A;	//create vector A2 copying A
	double *av = &A2[0];	//create pointer to first element in A2, this is what we will pass to dgesvd
	//double *av = new double[A.size()];
	//av = new double[A.size()];
	//for (int ii = 0; ii < r1*c1; ii++) {av[ii] = A[ii];}
	char jobU = 'S';
	char jobVt = 'S';
	int M = r1;
	int N = c1;
	
	//std::cout << M << '\t' << N << std::endl;
	
	int lda = M;
	int npc = std::min(M,N);
	std::vector<double> S22(npc);
	double *S2 = &S22[0];
	//double *S2= new double[npc];
	int ldu = M;
	//std::vector<double> U22(ldu*npc);
	std::vector<double> U22(ldu*M);
	double *U2 = &U22[0];
	//double *U2 = new double[ldu*npc];
	//int ldvt = 1;
	int ldvt = N;
	//std::vector<double> Vt22(1);
	std::vector<double> Vt22(ldvt*N);
	double *Vt2 = &Vt22[0];
	//double *Vt2 = new double[1];
	int lwork = std::max(3*npc+std::max(M,N),5*npc);
	//int lwork = -1;
	std::vector<double> work2(lwork);
	double *work = &work2[0];
	//double *work = new double[lwork];
	int info;
	::dgesvd_(&jobU, &jobVt, &M, &N, av, &lda, S2, U2, &ldu, Vt2, &ldvt, work, &lwork, &info);
/*	
		for (int ii = 0; ii < M*npc; ii++) std::cout << U2[ii] << '\t';
		std::cout<<std::endl;
*/
	if (info != 0)
	{
		std::cout << "SVD failed";
		mat2d Out;
		Out.vec_to_mat(A2,r1,c1);
		Out.print();
		std::exit (EXIT_FAILURE);
	}
	S.augment_vec(S22);
	//U.vec_to_mat(U22,ldu,std::min(M,N));
	//Vt.vec_to_mat(Vt22,ldvt,std::min(M,N));
	U.resize(ldu,std::min(M,N));
	//S.resize(npc,1);
	Vt.resize(ldvt,N);
	for (int ii = 0; ii < npc; ii++)
		{
			//S.set_value(ii,0,S22[ii]);
			for (int jj = 0; jj < M; jj++)
				{
					U.set_value(jj,ii,U22[ii*M+jj]);
				}
		}
	for (int ii = 0; ii < N; ii++)
	{
		for (int jj = 0; jj < ldvt; jj++)
		{
			Vt.set_value(jj,ii,Vt22[ii*N+jj]);
		}
	}
/*
	delete[] av;
	delete[] S2;
	delete[] U2;
	delete[] Vt2;
	delete[] work;
*/
}
	void mat2d::save_to_file(std::string name, std::string app, std::string matname) const
// A.save_to_file("outfile.txt","app","Matrix A");
{	
	std::fstream file;
	if (app == "app") 
		{
			file.open(name.c_str(), std::fstream::out | std::fstream::app);
		}
	else 
		{
			file.open(name.c_str(), std::fstream::out);
		}
	if (file.is_open())
	{
		//file.write(matname.c_str());
		file << matname << std::endl;
		for (int ii = 0; ii < r1; ii++)
			{
				for (int jj = 0; jj < c1; jj++)
					{
						file << A[jj*r1 + ii] << '\t';
					}	file << std::endl;
			}
		file << std::endl;
		file.close();
	}
	else std::cout << "Could not open " << name << " to write "  << matname << std::endl;
}
	//A.augment_vec(a) -> A := [A a]
	void mat2d::augment_vec(std::vector<double> &a)
{
	if (a.size() == r1)
		{
			c1++;
			A.resize(r1*c1);
			for (int ii = 0; ii < r1; ii++) A[r1*(c1-1) + ii] = a[ii];
		}
	//If dim(A) = 1 x 1 then A.augment_vec(a) -> A := a
	else if (c1 == 0 && r1 == 0)
		{
			r1 = a.size();
			c1 = 1;
			A.resize(r1*c1);
			for (int ii = 0; ii < r1; ii++) A[ii] = a[ii];
		}
	else 
	{
	std::cout << "Cannot augment a vector whose length does not equal the number of rows of the original matrix" << std::endl;
	std::exit(1);
	}
}
	//A.augment_mat(B) -> A := [A B]
	void mat2d::augment_mat(mat2d &a)
{
	r2 = a.rows();
	c2 = a.columns();
	if (r2 == r1)
		{
			c1 += c2;
			A.resize(r1*c1);
			std::vector<double> av(0,0);
			av = a.get_vec(0,r2,0,c2);
			for (int ii = 0; ii < r2*c2; ii++) A[r1*(c1-c2) + ii] = av[ii];
		}
	else std::cout << "Cannot augment a matrix whose number of rows do not match the number of rows of the original matrix " << std::endl;
}

void mat2d::augment_vecmat(std::vector<double> &a, int r)
{
	if (r == -1) {r2 = r1;}	//if no row dimension is given, assume to be same as original matrix
	else (r2 = r);	//if a row dimesnio is given, set
	if (r1 == 0 && c1 == 0) //if original matrix is empty
	{
		r1 = r2;
		//c1 = 0;
	}	
	c2 = a.size()/r1;	//check how many columns of data exist
	if (a.size() != r2*c2)
	{
		std::cout << "Cannot augment on a matrix in vector form of non-conforming size";
		std::exit(1);
	}
	c1 += c2;
	A.resize(r1*c1);
	for (int ii = 0; ii < r2*c2; ii++) A[r1*(c1-c2) + ii] = a[ii];
//	A.insert(A.end(),a.begin(),a.end());
}
		
	//A.augment_row(a) -> A:= [A ; a]
	void mat2d::augment_row(std::vector<double> &a)
{
	if (a.size() == c1)
		{
			r1++;
			std::vector<double> B(r1*c1,0);
			for (int ii = 0; ii < c1; ii++)
				{
					for (int jj = 0; jj < r1-1; jj++)
						{
							B[ii*r1+jj] = A[ii*(r1-1)+jj];
						}
					B[ii*r1 + r1-1] = a[ii];
				}
			A.resize(r1*c1);
			for (int ii = 0; ii < r1*c1; ii++) A[ii] = B[ii];
		}
	else std::cout << "Cannot augment on a row with a different length than number of columns in the original matrix" << std::endl;
}
	//A.qr(Q,R) -> A = Q*R
	void mat2d::qr(mat2d &q, mat2d &r) const
{
	std::vector<double> A2 = A;
	double *av = &A2[0];
	//av = new double[r1*c1];
	//for (int ii = 0; ii < r1*c1; ii++) av[ii] = A[ii];
	int M = r1;
	int N = c1;
	int lda = M;
	int K = std::min(M,N);
	std::vector<double> tau2(M);
	double *tau = &tau2[0];
	//tau = new double[M];
	int lwork = M;
	std::vector<double> work2(lwork);
	double *work = &work2[0];
	//work = new double[lwork];
	int info;
	::dgeqrf_(&M, &N, av, &lda, tau, work, &lwork, &info);
	std::cout<<"write r\n";
	/*
	A2.resize(c1*c1);
	r.vec_to_mat(A2,c1,c1);
	A2.resize(A.size());
	*/
	
	r.resize(c1,c1);
	for (int ii = 0; ii < c1; ii++)
		{
			for (int jj = 0; jj < ii+1; jj++)
				{
					r.set_value(jj,ii,A2[ii*r1 + jj]);
				}
		}
	::dorgqr_(&M, &N, &K, av, &lda, tau, work, &lwork, &info);
	std::cout<<"write q\n";
	/*
	A2.resize(r1*c1);
	q.vec_to_mat(A2,r1,c1);
	*/
	
	q.resize(r1,c1);
	for (int ii = 0; ii < c1; ii++)
		{
			for (int jj = 0; jj < r1; jj++)
				{
					q.set_value(jj,ii,A2[ii*r1 + jj]);
				}
		}
	
	/*
	delete[] av;
	delete[] tau;
	delete[] work;
	*/
}
	//Calculate norm
	double mat2d::norm(char type) const
{
	std::vector<double> A2 = A;
	double nval = 0.0;
	if (type == 'T' || type == 't')	//Trace
	{
		if (r1 == c1)
		{
			for (int ii = 0; ii < r1; ii++) {nval += A2[ii*r1 + ii];}
		}
		else
		{
			std::cout << "Cannot take trace of non-square matrix" << std::endl;
			return 0;
		}
	}
	else if (type == 'E' || type == 'e')	//min Eigenvalue
	{
		if (r1 == c1)
		{
			char jobz = 'N';
			char uplo2 = 'U';
			int n = r1, ldt = r1,lwork = 3*r1 ,info;
			std::vector<double> W2(n), work2(lwork);
			double  *av = &A2[0], *W = &W2[0], *work = &work2[0];
			::dsyev_(&jobz, &uplo2, &n, av, &ldt, W, work, &lwork, &info);
			nval = W[0];
		}
		else
		{
			std::cout << "Cannot perform eigenvalue decomposition on a non-square matrix" << std::endl;
			return 0;
		}
	}
	else	//'M','O','I','F'
	{
		double *av = &A2[0];
		std::vector<double> work2(r1);
		double *work = &work2[0];
		nval = ::dlange_(&type,&r1,&c1,av,&r1,work);
	}
	return nval;
}

	mat2d mat2d::scale_columns()
{
	std::vector<double> sum(c1,0);
	mat2d scaled(r1,c1);
	//std::vector<double> scaled(r1*c1,0);
	for (int ii = 0; ii < c1; ii++)
	{
//	std::cout << "sum" << std::endl;
		for (int jj = 0; jj < r1; jj++)
		{
			sum[ii] += A[ii*r1 + jj];
		}
//	std::cout << "scale" << std::endl;
		for (int jj = 0; jj < r1; jj++)
		{
			scaled.set_value(jj,ii,A[ii*r1 + jj]/sum[ii]);
		}
	}
	return scaled;
}
	//perform elementwise multiplication or division
	void mat2d::dot(char type,mat2d &a, mat2d &b)
{
	if (a.rows() == b.rows() && a.columns() == b.columns())
	{
		r1 = a.rows();
		c1 = a.columns();
		A.resize(r1*c1);
		std::vector<double> av(0,0);
		std::vector<double> bv(0,0);
		av = a.get_vec(0,r1,0,c1);
		bv = b.get_vec(0,r1,0,c1);
		if (type == 'M' || type == 'm')
		{
			for (int ii = 0; ii < r1*c1; ii++)
			{
				A[ii] = av[ii]*bv[ii];
			}
		}
		else if (type == 'D'|| type == 'd')
		{
			
			for (int ii = 0; ii < r1*c1; ii++)
			{
				if (std::abs(bv[ii]) >= 10e-4)
				{
//					std::cout << "dot divide" << std::endl;
					A[ii] = av[ii]/bv[ii];
				}
				else {A[ii] = av[ii];}
			}
		}
		else std::cout << "Unknown operation in dot operation" << std::endl;
	}
	else std::cout << "Cannot perform dot operation on different sized matrices" << std::endl;
}

void mat2d::erase_rc(int n1, int n2, int dimension)
{
//1 -> columns, 2-> rows
//n1 is first location being removed, n2 is location one after removal
//Ex. If removing columns 1-3 A.erase_rc(0,3,1) (recall, indexed like C++)
if (n2 == -1)
{
	n2 = n1+1;
}
if (n1 > n2)
{
	std::cout << "Cannot erase dimension (n1 > n2)";
	std::exit(1);
}
if (dimension == 1 && n2 > r1)
{
	std::cout << "erase_rc called on row after matrix dimension (n2 > r1)";
	std::exit(1);
}
else if (n1 < 0)
{
	std::cout << "erase_rc called on dimension before matrix dimension (n1 < r1 or n1 < c1)";
	std::exit(1);
}
else if (dimension == 2 && n2 > c1)
{
	std::cout << "erase_rc called on column after matrix dimension (n2 > c1)";
	std::exit(1);
}
else
{
	std::vector<double> B(0,0);
	if (dimension == 2)
	{
		for (int ii = 0; ii < n1; ii++)
		{
			for (int jj = 0; jj < r1; jj++)
			{
				B.push_back(A[ii*r1+jj]);
			}
		}
		for (int ii = n2; ii < c1; ii++)
		{
			for (int jj = 0; jj < r1; jj++)
			{
				B.push_back(A[ii*r1+jj]);
			}
		}
		c1 -= (n2-n1);
		A = B;
	}
	else if (dimension == 1)
	{
		for (int ii = 0; ii < c1; ii++)
		{
			for (int jj = 0; jj < n1; jj++)
			{
				B.push_back(A[ii*r1+jj]);
			}
			for (int jj = n2; jj < r1; jj++)
			{
				B.push_back(A[ii*r1+jj]);
			}
		}
		r1 -= (n2-n1);
		A = B;
	}
	else
	{
		std::cout << "Cannot erase in unknown dimension" << std::endl;
		std::exit(1);
	}
}
}


void mat2d::vec_to_mat(std::vector<double> &a, int r, int c, int dimension)
{
	if (a.size() != r*c)
	{
		std::cout << "Cannot turn this vector into a matrix";
		std::exit(1);
	}
	else if (dimension == 1) //reading in by row
	{
		A.resize(r*c);
		r1 = r;
		c1 = c;
		for (int ii = 0; ii < r1; ii++)
		{
			for (int jj = 0; jj < c1; jj++)
			{
				A[ii+jj*r1] = a[ii*c1+jj];
			}
		} 
	}
	else	//reading in by column
	{
		A = a;
		r1 = r;
		c1 = c;
	}
}

std::vector<double> mat2d::get_diag() const //Get diagonal from matrix
{
	if (r1 == c1)
	{
		std::vector<double> a(r1);
		for (int ii = 0; ii < r1; ii++)
		{
			a[ii] = A[ii*r1 + ii];
		}
		return a;
	}
	else
	{
		std::cout << "Cannot get diagonal from non-square matrix" << std::endl;
		std::exit(1);
	}
}

double mat2d::determinant() const
{
	if (r1 == c1)
	{
		std::vector<double> A2 = A;
		char jobz = 'N';
		char uplo2 = 'U';
		int n = r1, ldt = r1,lwork = 3*r1 ,info;
		std::vector<double> W2(n), work2(lwork);
		double  *av = &A2[0], *W = &W2[0], *work = &work2[0];
		::dsyev_(&jobz, &uplo2, &n, av, &ldt, W, work, &lwork, &info);
		double det = 1.0;
		for (int ii = 0; ii < n; ii++)
		{
			det *= W2[ii];
		}
		return det;
	}
	else
	{
		std::cout << "Cannot calculate determinant of non-square matrix" << std::endl;
		std::exit(1);
	}
}

void mat2d::multiply_const(double m)
{
	for (unsigned int ii = 0; ii < A.size(); ii++)
	{
		A[ii] *= m;
	}
}