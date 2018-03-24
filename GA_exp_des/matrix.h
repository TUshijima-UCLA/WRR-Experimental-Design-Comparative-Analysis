#ifndef matrix_h
#define matrix_h

//#include <cmath>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
/*
//General Matrix Vector Multiplication (BLAS)
extern "C" void dgemv_(const char *trans, const int *M, const int *N, double *alpha, double *A,
const int *lda, double *x, const int *incx, double *beta, double *y, const int *incy);
//General Matrix Multiply (BLAS)
extern "C" void dgemm_(const char *transA, const char *transB, const int *M, const int *N, const int *K,
double *alpha, double *A, const int *lda, double *B, const int *ldb,
double *beta, double *C, const int *ldc);
//SVD (LAPACK)
extern "C" void dgesvd_(const char *jobu, const char *jobvt, const int *M, const int *N, 
double *A, const int *lda, double *S, double *U, const int *ldu,
double *VT, const int *ldvt, double *work, const int *lwork, int *info);
//Forms R part of QR Decomp (LAPACK)
extern "C" void dgeqrf_(int *M, int *N, double *A, int *lda, double *tau, double *work,
int *lwork, int *info);
//Forms Q part of QR Decomp (LAPACK)
extern "C" void dorgqr_(int *M, int *N, int *K, double *A, int *lda, double *tau,
double *work, int *lwork, int *info);
//Norm of a matrix (LAPACK)
extern "C" double dlange_(char *norm, int *M, int *N, double *A, int *lda, double*work);
*/

class mat2d
{
//indexed like normal with C++ (index 0 is the first element)
private:
int r1, r2, r3, c1, c2, c3;
std::vector<double> A;
public:
//Constructor (mat2d A(m,n))
mat2d (int m = 0, int n = 0);
~mat2d();
int rows() const;
int columns() const;
void resize(int m, int n);
//Set element m,n to value (A.set_value(m,n,value)) //store by column (first m values are the first column, FORTRAN style)
void set_value(int m, int n, double value);
//Get the value of element m,n (double value = A.get_value(m,n))
double get_value(int m, int n) const;
//Turn full or part of matrix into a vector (std::vector<double> vec = A.get_vec(m1,m2,n1,n2) m2 > =m1, n2 >= n1)
std::vector<double> get_vec(int m1, int m2, int n1, int n2) const;
//Print out matrix to screen (A.print())
void print() const;
//Set one matrix equal to another (A.equals(B) -> A = B)
void equals(mat2d &b); //copies the matrix given
//Add two matrices (C.add(A,B) -> C = A + B)
void add(mat2d &a, mat2d &b);
//Subtract two matrices (C.subtract(A,B) -> C = A-B)
void subtract(mat2d &a, mat2d &b);
//Transpose matrix (mat2d B = A.transpose() -> B = A')
mat2d transpose();
//Clears matrix (alternate to destructor A.clear() -> dim(A) = 1x1, A(0,0) = 0)
void clear();
//Multiply a matrix by a vector (b = A.matvec_mult(alpha,'T',x,beta,y) -> b = alpha*A'*x + beta*y)
std::vector<double> matvec_mult(double alpha, char ta, std::vector<double> &x,
double beta, std::vector<double> y) const;
//Multiply two matrices (D.multiply('T',A,'N',B,C,alpha,beta) -> D = alpha*A' * B + beta*C)
//Be careful what you put in for C, its size can change on exit
void multiply(char ta, mat2d &a, char tb, mat2d &b, mat2d &c, double alpha, double beta);
//A.svd(U,S,Vt) -> A = USV'
void svd(mat2d &U, mat2d &S, mat2d &Vt) const;
// A.save_to_file("outfile.txt","app","Matrix A");
void save_to_file(std::string name, std::string app, std::string matname) const;
//A.augment_vec(a) -> A := [A a]
void augment_vec(std::vector<double> &a);
//A.augment_mat(B) -> A := [A B]
void augment_mat(mat2d &a);
//A.augment_vecmat(a) -> A:= [A A2] (a is vector form of matrix A2)
void augment_vecmat(std::vector<double> &a, int r = -1);
//A.augment_row(a) -> A:= [A ; a]
void augment_row(std::vector<double> &a);
//A.qr(Q,R) -> A = QR
void qr(mat2d &q, mat2d &r) const;
//Calculates norm of a matrix (norm = I.norm('char'))
//'T':Trace, 'E':min eigenvalue, 'O':max column sum, 'I':max row sum, 'F':Frobenius (sqrt(sum squares))
//'M':max element
double norm(char type) const;
//Scales the rows of a matrix to sum to 1
mat2d scale_columns();
//performs elementwise multiplication or division
//C.dot('M',A,B) -> C = A.*B;
void dot(char type,mat2d &a, mat2d &b);
//Erase elements n1-(n2-1) from the "dimension" of the matrix
//dimension = 2 -> columns; dimension = 1 -> rows
//Ex.   erasing columns 1-3, A.erase_rc(0,3,2)
//Ex2.  erasing row 3, A.erase_rc(2,3,1)
//Ex2a. erasing column 3, A.erase_rc(2)
void erase_rc(int n1, int n2 = -1, int dimension = 2);
//A.vec_to_mat(a,r1,c1) -> A:= A2 (a is vector from of matrix A2) 
//default data stored by column (dimension = 2), dimension = 1 -> data stored by row
void vec_to_mat(std::vector<double> &a, int r, int c, int dimension = 2);
//Get diagonal from matrix
std::vector<double> get_diag() const;
double determinant() const;
void multiply_const(double m);
/*
	//Constructor (mat2d A(m,n))
	inline mat2d(int, int);
	//Destructor (not sure how to use)
	~mat2d() {A.clear();}
	//Get number of rows in matrix (int m = A.rows())
	int rows() {return (r1);}
	//Get number columns in matrix (int n = A.columns())
	int columns() {return (c1);}
	//Resize the matrix (A.resize(m,n))
	void resize(int m, int n)
		{	
			r1 = m;
			c1 = n;
			A.resize(r1*c1);
		}
	//Set element m,n to value (A.set_value(m,n,value)) //store by column (first m values are the first column, FORTRAN style)
	void set_value(int m, int n, double value)
		{
			A[(n)*r1 + (m)] = value;
		}
	//Get the value of element m,n (double value = A.get_value(m,n))
	double get_value(int m, int n)
		{
			double a;
			a = A[(n)*r1 + (m)];
			return (a);
		}
	//Turn full or part of matrix into a vector (std::vector<double> vec = A.get_vec(m1,m2,n1,n2) m2 > =m1, n2 >= n1)
	std::vector<double> get_vec(int m1, int m2, int n1, int n2)
		{//1 is first index 2 is second index 
		//indexing brackets selection (if want A(:,1) m1 = 0, m2 = r1, n1 = 0, n2 = 1)
		//if take all elements in the matrix (0,r1,0,c1) vector is in correct LAPACK form
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
		
	//Print out matrix to screen (A.print())
	void print()
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
		
	//Set one matrix equal to another (A.equals(B) -> A = B)
	void equals(mat2d &b) //copies the matrix given
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
	void add(mat2d &a, mat2d &b)
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
	void subtract(mat2d &a, mat2d &b)
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
	mat2d transpose()
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
	//Clears matrix (alternate to destructor A.destroy() -> dim(A) = 1x1, A(0,0) = 0)
	void destroy()
		{
			r1 = 1;
			c1 = 1;
			A.clear();
		}
	
	//Multiply a matrix by a vector (b = A.matvec_mult(alpha,'T',x,beta,y) -> b = alpha*A'*x + beta*y)
	std::vector<double> matvec_mult(double alpha, char ta, std::vector<double> &x,
	double beta, std::vector<double> y)
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
			
			if (y.size() != x.size())
			{
				//will a y of non-comforming size and copy the first x.size() elements
				y.resize(x.size());
			}
			double *y2 = new double[y.size()];
			double *x2 = new double[x.size()];			
			for (int ii = 0; ii < y.size(); ii++) 
			{
				y2[ii] = y[ii];
				x2[ii] = x[ii];
			}
			
			double *A2 = new double[r1*c1];
			for (int ii = 0; ii < r1*c1; ii++) {A2[ii] = A[ii];}
			//std::cout << x.size() << '\t' << r1 << std::endl;
			if ((ta == 'T' && x.size() == r1) || (ta == 'N' && x.size() == c1))
			//if (trans == 'T' && x.size() == r1)
			{
				//std::cout << ta << std::endl;
				::dgemv_(&trans, &r1, &c1, &alpha, A2, &r1, x2, &incx, &beta, y2, &incy);
			}
//			else if (ta = ='N' && x.size() == c1) {::dgemv_(&trans, &r1, &c1, &alpha, A2, &r1, x2, &incx, &beta, y2, &incy);}
			else std::cout << "Cannot multiply non-conforming vector to a matrix" << std::endl;
			
			//std::vector<double> b(y.size(),0);
			for (int ii = 0; ii < b.size(); ii++) {b[ii] = y2[ii];}
			
			return (b);
			
			delete[] y2;
			delete[] x2;
			delete[] A2;
			
		}

	//Multiply two matrices (D.multiply('T',A,'N',B,C,alpha,beta) -> D = alpha*A' * B + beta*C)
	//Be careful what you put in for C, its size can change on exit
	void multiply(char ta, mat2d &a, char tb, mat2d &b, mat2d &c, double alpha, double beta)
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
			double *av;
			av = new double[r2*c2];
			double *bv;
			bv = new double[r3*c3];
			double *cv;
			cv = new double[ldc*M];
			std::vector<double> av2;
			std::vector<double> bv2;
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
				for (int ii = 0; ii < av2.size(); ii++) av[ii] = av2[ii];
				for (int ii = 0; ii < bv2.size(); ii++) bv[ii] = bv2[ii];
				//use :: to call global function in member function section
				::dgemm_(&transA, &transB, &M, &N, &K, &alpha, av, &lda, bv, &ldb, &beta, cv, &ldc);
				r1 = ldc;
				c1 = N;
				A.resize(r1*c1);
				for (int ii = 0; ii < r1*c1; ii++) A[ii] = cv[ii];
			}
			else std::cout << "Cannot multiply, matrices dimensions do not match (also check dimensions of C)" << std::endl;
			delete[] av;
			delete[] bv;
			delete[] cv;
		}
	void svd(mat2d &U, mat2d &S, mat2d &Vt)
		{
			double *av;
			av = new double[A.size()];
			for (int ii = 0; ii < r1*c1; ii++) av[ii] = A[ii];
			char jobU = 'S';
			char jobVt = 'N';
			int M = r1;
			int N = c1;
			int lda = M;
			double *S2;
			int npc = std::min(M,N);
			S2 = new double[npc];
			int ldu = M;
			double *U2;
			U2 = new double[ldu*npc];
			int ldvt = 1;
			double *Vt2;
			Vt2 = new double[1];
			int lwork = std::max(3*npc+std::max(M,N),5*npc);
			double *work;
			work = new double[lwork];
			int info;
			::dgesvd_(&jobU, &jobVt, &M, &N, av, &lda, S2, U2, &ldu, Vt2, &ldvt, work, &lwork, &info);
			if (info == 0) std::cout << "SVD sucessful" << std::endl;
			else
			{
				std::cout << "SVD failed";
				std::exit (EXIT_FAILURE);
			}
			U.resize(M,npc);
			S.resize(npc,1);
			Vt.resize(1,1);
			for (int ii = 0; ii < npc; ii++)
				{
					S.set_value(ii,0,S2[ii]);
					for (int jj = 0; jj < M; jj++)
						{
							U.set_value(jj,ii,U2[ii*M+jj]);
						}
				}

			delete[] av;
			delete[] S2;
			delete[] U2;
			delete[] Vt2;
			delete[] work;
		}
	void save_to_file(std::string name, std::string app, std::string matname)
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
	void augment_vec(std::vector<double> &a)
		{
			if (a.size() == r1)
				{
					c1++;
					A.resize(r1*c1);
					for (int ii = 0; ii < r1; ii++) A[r1*(c1-1) + ii] = a[ii];
				}
			//If dim(A) = 1 x 1 then A.augment_vec(a) -> A := a
			else if (c1 == 1 && r1 == 1)
				{
					r1 = a.size();
					A.resize(r1*c1);
					for (int ii = 0; ii < r1; ii++) A[ii] = a[ii];
				}
			else std::cout << "Cannot augment a vector whose length does not equal the number of rows of the original matrix" << std::endl;
		}
	//A.augment_mat(B) -> A := [A B]
	void augment_mat(mat2d &a)
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
	//work here	
	//A.augment_row(a) -> A:= [A ; a]
	void augment_row(std::vector<double> &a)
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
	void qr(mat2d &q, mat2d &r)
		{
			double *av;
			av = new double[r1*c1];
			for (int ii = 0; ii < r1*c1; ii++) av[ii] = A[ii];
			int M = r1;
			int N = c1;
			int lda = M;
			int K = std::min(M,N);
			double *tau;
			tau = new double[M];
			int lwork = M;
			double *work;
			work = new double[lwork];
			int info;
			::dgeqrf_(&M, &N, av, &lda, tau, work, &lwork, &info);
			r.resize(c1,c1);
			for (int ii = 0; ii < c1; ii++)
				{
					for (int jj = 0; jj < ii+1; jj++)
						{
							r.set_value(jj,ii,av[ii*r1 + jj]);
						}
				}
			::dorgqr_(&M, &N, &K, av, &lda, tau, work, &lwork, &info);
			q.resize(r1,c1);
			for (int ii = 0; ii < c1; ii++)
				{
					for (int jj = 0; jj < r1; jj++)
						{
							q.set_value(jj,ii,av[ii*r1 + jj]);
						}
				}
			delete[] av;
			delete[] tau;
			delete[] work;
		}
	//Calculate norm
	double norm(char type)
		{
			double *av = new double[r1*c1];
			for (unsigned int ii = 0; ii < r1*c1; ii++) av[ii] = A[ii];
			double *work = new double[r1];
			double nval;
			nval = ::dlange_(&type,&r1,&c1,av,&r1,work);
			
			delete[] av;
			delete[] work;
			return nval;
		}
	mat2d scale_columns()
		{
			std::vector<double> sum(c1,0);
			mat2d scaled(r1,c1);
//			std::vector<double> scaled(r1*c1,0);
			for (int ii = 0; ii < c1; ii++)
			{
//			std::cout << "sum" << std::endl;
				for (int jj = 0; jj < r1; jj++)
				{
					sum[ii] += A[ii*r1 + jj];
				}
//			std::cout << "scale" << std::endl;
				for (int jj = 0; jj < r1; jj++)
				{
					scaled.set_value(jj,ii,A[ii*r1 + jj]/sum[ii]);
				}
			}
			return scaled;
		}
	void dot(char type,mat2d &a, mat2d &b)
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
//							std::cout << "dot divide" << std::endl;
							A[ii] = av[ii]/bv[ii];
						}
						else {A[ii] = av[ii];}
					}
				}
				else std::cout << "Unknown operation in dot operation" << std::endl;
			}
			else std::cout << "Cannot perform dot operation on different sized matrices" << std::endl;
		}
*/
};
/*
mat2d::mat2d(int m, int n)
{
r1 = m;
c1 = n;
A.resize(r1*c1);
}

mat2d::~mat2d() {A.clear();}
*/

//#include "matrix.cpp"
#endif