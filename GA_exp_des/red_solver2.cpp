#include "red_solver.h"
#include "matrix.h"
//#include <time.h>
//#include <chrono>

void red_solver(std::vector<double> &P, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &well, std::vector<double> &wellrate, int &npc, int Nn, int &nq,
std::vector<double> &ts, double &dtmult, std::vector<double> &Jdv, int nprt, int lumpwell, 
double dt, mat2d &Ps)
{

//std::chrono::steady_clock::time_point _start(std::chrono::steady_clock::now());

//std::cout << "red solver" << std::endl;
ncalls++;
//double dtstore;
//std::vector<double> P(0,0);
//P = Pt.get_vec(0,Nn,0,npc);

//Write P to pointer
//double *Pn = &P[0];
std::vector<double> An2(Ar.size());
double *An = &An2[0];
int nrhs = 1;
int lwork = std::max(1, npc*npc + std::max(npc*npc, 1));
std::vector<double> work2(lwork);
double *work = &work2[0];
int info;

char uplo = 'L';
std::vector<int> ipiv2(npc);
int *ipiv = &ipiv2[0];

//setup for matrix vector multiplication
char trans = 'N';
//double alpha = 1.0;
double alphadt;
double beta = 0.0;
int incx = 1;
int incy = 1;
//Creating right-hand side
std::vector<double> v2(npc);
double *v = &v2[0];
double *Brn = &Br[0];

//std::vector<double> Y2(Nn,0);
//double *y = &Y2[0];

std::vector<double> qr;
//std::vector<double> rt(npc,0.0), rt1(npc,0.0);	//stores current and previous red solutions
//Looping through the wells
int sets;	//number of sets of wells pumping at a fixed ratio

if (lumpwell == 1) {sets = 1;}
else {sets = nq;}

std::vector<double> b2;
//dtstore.resize(sets);
std::vector<double> Jdi;		//vector to store current column of Jd
//std::cout << "start red model" << std::endl;
//Include all wells
//std::cout << "sets = " << sets << std::endl;
for (int currentwell = 0; currentwell < sets; currentwell++)
{
	b2.clear();
	b2.resize(npc,0.0);
	//std::cout << "b pointer" << std::endl;
	double *b = &b2[0];
	//double dt = ts[nprt*currentwell];
	//double sstime = ts[nprt*currentwell + nprt-1];
	int kk = 0;
	int simstep = 0;
	double simtime = 0.0;
	qr.clear();
	//std::cout << "qr resize" << std::endl;
	qr.resize(npc, 0.0);
	double dtstore1 = dt;
	//std::cout << "set T" << std::endl;
	//for (unsigned int ii = 0; ii < ts.size(); ii++) {std::cout << ts[ii] << std::endl;}
	
	
	double T = ts[nprt*currentwell + nprt-1] + dt;			//T = max sim time
	//std::cout << "T = " << T << std::endl;
	//Use only one well at a time
	if (lumpwell == 1)
	{
		for (int ii = 0; ii < npc; ii++)
		{
			for (int jj = 0; jj < nq; jj++)
			{
				qr[ii] = qr[ii] + P[ii*Nn + well[jj]-1]*wellrate[jj];
			}
		}
	}
	if (lumpwell != 1) 
	{
		
		for (int ii = 0; ii < npc; ii++)
		{
				qr[ii] = qr[ii] + P[ii*Nn + well[currentwell]-1]*wellrate[currentwell+1];
		}
	}
	
	//std::cout << "qr =" << std::endl;
	//for (unsigned int ii = 0; ii < qr.size(); ii++) {std::cout << qr[ii] << '\t';}
	
	//Run reduced model
	while (simtime < T)
	{
		simstep++;
		alphadt = 1/dt;
		for (unsigned int ii = 0; ii < Ar.size(); ii++)
		{
			An2[ii] = -Ar[ii] - Br[ii]/dt;
		}
		dgemv_(&trans, &npc, &npc, &alphadt, Brn, &npc, b , &incx, &beta, v, &incy);
		for (int ii = 0; ii < npc; ii++)
		{
			b2[ii] = qr[ii] - v2[ii];
		}
		dsysv_(&uplo, &npc, &nrhs, An , &npc, ipiv, b, &npc, work, &lwork, &info);
		//rt1 = rt;
		//rt = b2;
		if (kk < nprt && simtime == ts[currentwell*nprt + kk])
		{
			std::vector<double> Jdvi = Ps.matvec_mult(1.0,'N',b2,0.0,b2);
			//std::cout << "Jdv.size(): " << Jdv.size() << "\tJdvi.size(): " << Jdvi.size() << std::endl;
			//std::cout << "Jdvi.insert()" << std::endl;
			
			//store sensitivies*well rate
			Jdv.insert(Jdv.end(),Jdvi.begin(),Jdvi.end());

			//Jdv.insert(Jdv.end(),b2.begin(),b2.end());
			
			//std::cout << Jdv.size() << std::endl;
			
			
			//dgemv_(&trans, &Nn, &npc, &alpha, Pn, &Nn, b, &incx, &beta, y, &incy);
			//Jdi.insert(Jdi.end(),Y2.begin(),Y2.end());
			
			//store reduced solution at observation times
			//r.augment_vec(b2);
			//rv.insert(rv.end(),b2.begin(),b2.end());
			///////////Comment on reduced model output ///////////////
			/*
			rv = [r(k_i,q1); r(k_i,q2) ; ... ; r(k_i,q_n)]
			rv is a vector of dimension nobs*nowell*nq
			*/
			kk++;
		}
		//std::cout << "Jdv.size() = " << Jdv.size() << std::endl;
		dt = dtstore1;
		//grows time step
		dt *= dtmult;
		//limits time step
		if (dt > 2) {dt = 2;}
		simtime += dt;
		dtstore1 = dt;
		//gets exact output times
		if (kk < nprt && simtime > ts[currentwell*nprt + kk])
		{
			dt = simtime - ts[currentwell*nprt + kk];
			simtime = ts[currentwell*nprt + kk];
		}
	}
}
//std::cout << "Jdv.size() = " << Jdv.size() << std::endl;
/*
std::chrono::steady_clock::time_point _end(std::chrono::steady_clock::now());
std::cout << "Time to run reduced model: ";
std::cout << std::chrono::duration_cast<std::chrono::duration<double>>(_end - _start).count();
std::cout << std::endl;
*/
//std::exit(1);
}