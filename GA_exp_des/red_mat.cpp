#include <vector>
#include "addon.h"
#include "fort.h"
#include "matrix.h"

extern "C" {
void assf2d_(int &nt, int &irad, int *triang, int *trija, double *lmass,
double *coef1, double *coef2, double *xc, double *xy, double *permx,
double *permy, double *elstor, double *spess, double *area, double *arear,
double *bi, double *ci);
}

void axbsym(const int &N, const int &NTERM, std::vector<int> &TOPOL,
std::vector<int> &JA, std::vector<double> &coef1, std::vector<double> &XVEC,
std::vector<double> &BVEC);
/*
void sys_mat(int &nt, int &irad, int &nterm, std::vector<int> &triang2,
std::vector<int> &trija2, std::vector<double> &lmass2, std::vector<double> &xc2,
std::vector<double> &yc2, std::vector<double> &elstor2, std::vector<double> &spess2,
std::vector<double> &area2, std::vector<double> &arear2, std::vector<double> &bi2,
std::vector<double> &ci2, std::vector<double> &coef12, std::vector<double> &coef22,
std::vector<double> &permc)
{
int *triang = &triang2[0];
//int *triang = new int[triang2.size()];
//for (int ii = 0; ii < triang2.size(); ii++) {triang[ii] = triang2[ii];}
int *trija = &trija2[0];
//int *trija = new int[trija2.size()];
//for (int ii = 0; ii < trija2.size(); ii++) {trija[ii] = trija2[ii];}
double *lmass = &lmass2[0];
//double *lmass = new double[lmass2.size()];
//for (int ii = 0; ii < lmass2.size(); ii++) {lmass[ii] = lmass2[ii];}
double *xc = &xc2[0];
double *yc = &yc2[0];
double *elstor = &elstor2[0];
double *spess = &spess2[0];
double *area = &area2[0];
double *arear = &arear2[0];
double *bi = &bi2[0];
double *ci = &ci2[0];
coef12.resize(nterm, 0.0);
coef22.resize(nterm, 0.0);
double *coef1 = &coef12[0];
double *coef2 = &coef22[0];
double *permx = &permc[0];
double *permy = &permc[0];

std::cout << "assf2d" << std::endl;
assf2d_(nt, irad, triang, trija, lmass, coef1, coef2, xc, yc, permx, permy,
elstor, spess, area, arear, bi, ci);
std::cout << "return from assf2d" << std::endl;

}
*/

void redmat(std::vector<double> &Pv2, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &topol, std::vector<int> &JA, int &nt, int &irad, int &nterm,
std::vector<int> &triang2, std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, std::vector<double> &elstor2, 
std::vector<double> &spess2, std::vector<double> &area2, 
std::vector<double> &arear2, std::vector<double> &bi2, 
std::vector<double> &ci2, int &npc, int &Nn, std::vector<double> &permc)
{
std::vector<double> coef12(0,0), coef22(0,0),
dum1(nterm,0), dum2(Nn,0), dum3(0,0);
//std::cout << "sys_mat" << std::endl;

int *triang = &triang2[0];
int *trija = &trija2[0];
double *lmass = &lmass2[0];
double *xc = &xc2[0];
double *yc = &yc2[0];
double *elstor = &elstor2[0];
double *spess = &spess2[0];
double *area = &area2[0];
double *arear = &arear2[0];
double *bi = &bi2[0];
double *ci = &ci2[0];
coef12.resize(nterm, 0.0);
coef22.resize(nterm, 0.0);
double *coef1 = &coef12[0];
double *coef2 = &coef22[0];
double *permx = &permc[0];
double *permy = &permc[0];




assf2d_(nt, irad, triang, trija, lmass, coef1, coef2, xc, yc, permx, permy,
elstor, spess, area, arear, bi, ci);

Ar.resize(npc*npc);
Br.resize(npc*npc);
//	Nr1 = coef1 * P
std::vector<double> Nr12(Nn*npc);
double *Nr1 = &Nr12[0];
//AP.clear();
//AP.resize(Nn,npc);
//std::cout << "Pv2.size() = " << Pv2.size() << std::endl;
for (int ii = 0; ii < nterm; ii++) {dum1[ii] = coef1[ii];}
//std::cout << "axbsym" << std::endl;
for (int ii = 0; ii < npc; ii++)
{
	for (int jj = 0; jj < Nn; jj++) {dum2[jj] = Pv2[ii*Nn+jj];}
	dum3.clear();
	dum3.resize(Nn);
	axbsym(Nn,nterm,topol,JA,dum1,dum2,dum3);
//	std::cout << ii << std::endl;
	for (int jj = 0; jj < Nn; jj++) 
	{
		Nr12[ii*Nn + jj] = dum3[jj];
//		AP.set_value(jj,ii,dum3[jj]);
	}
}

//	Ar = P'*coef1*P
//std::cout << "Ar" << std::endl;
double *Arn = &Ar[0];
double *Pv = &Pv2[0];
//std::cout << "Pv" << std::endl;
//for (int ii = 0; ii < Pv2.size(); ii++) {Pv[ii] = Pv2[ii];}

char transa = 'T';	//'T' -> op(A) = A', 'N' -> op(A) = A
char transb = 'N';	//same as above
double alpha = 1.0;	//alpha*op(A)*op(B)
double beta = 0.0;	//beta*C
//std::cout << "dgemm" << std::endl;
//degemm_(transa, transb, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
//M =  # rows op(A), N = # col op(B), K = # col op(A), LDA = M, LDB = K, LDC = M
dgemm_(&transa, &transb, &npc, &npc, &Nn, &alpha, Pv, &Nn, Nr1, &Nn,
&beta, Arn, &npc);

std::vector<double> coef2diag(topol.size(), 0.0);

for (unsigned int ii = 0; ii < topol.size(); ii++) {coef2diag[ii] = coef2[topol[ii]-1];}

// Mr1 = coef2*P;
//std::cout << "Mr1" << std::endl;
std::vector<double> Mr12(Nn*npc);
double *Mr1 = &Mr12[0];
//BP.clear();
//BP.resize(Nn,npc);
//std::cout << BP.rows() << '\t' << BP.columns() << std::endl;
for (int ii = 0; ii < npc; ii++)
{
	for(int jj = 0; jj < Nn; jj++)
	{
		Mr12[ii*Nn + jj] = coef2diag[jj]*Pv[ii*Nn + jj];
//		BP.set_value(jj,ii,Mr12[ii*Nn+jj]);
	}
}

//std::cin.get();	//Pause
//Br = P'*coef2*P
double *Brn = &Br[0];
dgemm_(&transa, &transb, &npc, &npc, &Nn, &alpha, Pv, &Nn, Mr1, &Nn,
&beta, Brn, &npc);
//std::cout << "Br" << std::endl;

}