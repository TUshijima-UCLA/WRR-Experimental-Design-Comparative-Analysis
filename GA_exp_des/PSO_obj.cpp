//Standard headers
#include <set>
#include <cstdlib>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <ctime>
#include <iomanip>
//#include <cmath>
#include <math.h>
#include <numeric>      // std::accumulate
//Personal headers
#include "owellstruct.h"
#include "matrix.h"
#include "float.h"

//PSO headers
#include <SwarmOps/Tools/Types.h>

//Global variables
extern mat2d P;
extern std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, permc, Pv, wellrate, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;

extern std::vector<double> lb1, ub1, dI,lbi;
extern double multconst, bestfvalue;

extern bool runlongernow;

extern std::vector<double> diagh;
extern std::vector<int> poi;
extern int ntoi,npoi;

extern std::vector<int> topol, JA, triang, trija, well, ozonec;
extern std::vector<owell> obswells;			//all potential obsrvation nodes
extern int ntri, irad, nterm, npc, Nn, nprt, lumpwell, pickgenome, genomeoption, ncalls, 
nowell, nobszones, nq, nz, expdesparam, expdespump;	//nz: total number of hydrologic zones in the aquifer
extern int nopcalls;
extern double dtmult, dt, dH, dQ, convergence;
extern char OptCI, OptCG;
extern float z, zmin;

void redmat(std::vector<double> &Pv2, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &topol, std::vector<int> &JA, int &nt, int &irad, int &nterm,
std::vector<int> &triang2, std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, std::vector<double> &elstor2, 
std::vector<double> &spess2, std::vector<double> &area2, 
std::vector<double> &arear2, std::vector<double> &bi2, 
std::vector<double> &ci2, int &npc, int &Nn, std::vector<double> &permc);

void red_solver(std::vector<double> &P, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &well, std::vector<double> &wellrate, int &npc, int Nn, int &nq,
std::vector<double> &ts, double &dtmult, std::vector<double> &Jdv, int nprt, int lumpwell, 
double dt, mat2d &Ps);

extern clock_t t3, t4;

double Objective(const SO_TElm *x, int n)
{
	/*
	t4 = clock();
	printf("\nTime to return: %g seconds\n", (double)(t4 - t3) / CLOCKS_PER_SEC);
	t3 = t4;

	clock_t t5 = clock();
	*/
	nopcalls++;
	if (nopcalls % 10 == 0)
	{
		std::cout << std::resetiosflags(std::ios_base::scientific) <<
		'\r' << "Current model call number: " << ncalls << " Objective calls: " << nopcalls << " best score: " <<
		bestfvalue << "      ";
	}
	double fvalue = 1E100;
	std::vector<double> particle(n,0.0);
	
	for (int ii = 0; ii < n; ii++)
	{
		particle[ii] = x[ii];
	}

	std::vector<owell> cowell;
	mat2d Ps; //"small" projection matrix (only contain observation well locations)
	std::vector<double> Psv(0,0);	//vector holding Ps in row form		
	int countwell;
	
	if (pickgenome == 1)	//pick the first genome
	{
		countwell = nobszones;
		for (int ii = 0; ii < n; ii++)
		{
			if (particle[ii] == -1) {countwell--;}
		}
		if (countwell > nowell)
		{
			if (genomeoption == 0)
			{
				fvalue = 1E100;
				return fvalue;
			}
			else if (genomeoption == 1)
			{
				std::cout << "No method to fix particle" << std::endl;
				std::exit(1);
			}
		}
		
		if ((OptCI == 'G' || OptCI == 'I') && (countwell < nowell))
		{
			fvalue = 1E100;
			return fvalue;
		}
		else if (genomeoption == 1)
		{
			std::cout << "No method to fix particle" << std::endl;
			std::exit(1);
		}
		int currentzero = 0;
		for (int ii = 0; ii < n; ii++)
		{
			if (particle[ii] != -1)
			{
				std::vector<double> pvi = P.get_vec(obswells[particle[ii] + currentzero].node -1,
				obswells[particle[ii]+currentzero].node,0,npc);
				Psv.insert(Psv.end(),pvi.begin(),pvi.end());
			}
			currentzero += ozonec[ii];	
		}
	}
	else if (pickgenome == 2)
	{
		std::vector<int> owellzones;
		std::vector<owell> cowells;
		for (int ii = 0; ii < n; ii++)	//n = particle.size()
		{
			owellzones.push_back(obswells[particle[ii]].zone);
		}
		std::set<int> check(owellzones.begin(),owellzones.end());	//a set only holds one copy of a value
		if (check.size() < owellzones.size())	//if there are repeated zones 
		{
			fvalue = 1E100;
			return fvalue;
		}
		for (int ii = 0; ii < n; ii++)
		{
			std::vector<double> pvi = P.get_vec(obswells[particle[ii]].node-1, obswells[particle[ii]].node,0,npc);
//			std::cout << genome.gene(ii) << '\t' << obswells[genome.gene(ii)].node <<'\t' << obswells[genome.gene(ii)].zone << '\n';
			Psv.insert(Psv.end(),pvi.begin(),pvi.end());
			countwell = nowell;
		}
		//std::cout << "No method to pick second type of particle" << std::endl;
		//std::exit(1);
	/*
		std::vector<int> owellzones;
		std::vector<owell> cowells;	//obs wells currently being considered
		for (int ii = 0; ii < genome.length(); ii++)
		{
			owellzones.push_back(obswells[genome.gene(ii)].zone);
		}
		std::set<int> check(owellzones.begin(),owellzones.end());	//a set only holds one copy of a value
		if (check.size() < owellzones.size())	//if there are repeated zones 
		{
			fvalue = 0.0;
			return fvalue;
		}
		for (int ii = 0; ii < genome.size(); ii++)
		{
			std::vector<double> pvi = P.get_vec(obswells[genome.gene(ii)].node-1, obswells[genome.gene(ii)].node,0,npc);
			Psv.insert(Psv.end(),pvi.begin(),pvi.end());
			countwell = nowell;
		}
		*/
	}
	else
	{
		std::cout << "Unknown genome" << std::endl;
		std::exit(1);
	}
	if (Psv.size() == 0)
	{
		fvalue = 1E100;
		return fvalue;
	}
	Ps.vec_to_mat(Psv,countwell,npc,1);
	//Get baseline
	std::vector<double> permstore = permc;
	/*
	clock_t t6 = clock();
	printf("\nTime to build sys matrices: %g seconds\n", (double)(t6 - t5) / CLOCKS_PER_SEC);
	*/
	std::vector<double> Ar(npc*npc), Br(npc*npc);
	std::vector<double> Jdvbase, Jdv, Jdiv;	//reduced solutions and columns of the Jacobian in vector form
	mat2d Jd, I, H; //sensitivity, information, and hat matrices
	for (int ii = 0; ii < npc*npc; ii++)
	{
		Ar[ii] = Arbuild[ii];
		Br[ii] = Brbuild[ii];
	}
	red_solver(Pv,Ar,Br,well,wellrate,npc,Nn,nq,obstimes,dtmult,Jdvbase,nprt,lumpwell,dt,Ps);
	/*
	clock_t t7 = clock();
	printf("\nTime base Jd: %g seconds\n", (double)(t7 - t6) / CLOCKS_PER_SEC);
	*/
	/*
	////////////Comments on Jd/////////////////////
	
	Jd =[ds(q1)/dH1 ds(q2)/dH1 ... ds(qn)/dH1 ds(q1)/dH2 ... ds(qn)/dHn2 ds(k1)/dq1 ... ds(k1)/dqn ds(k2)/dq1 ... ds(kn2)/dqn]
	Jd is of dimesnon (nobs*countwell) x  (2*nq*nz) if you're interested in both pumping wells and zones,
	if you're just intersted in pumping wells the columns are just nq, if you're intersted in zones, 
	then the columns are just nz
	*/

	int nobst, nptoi;
	if (OptCI == 'G' || OptCI == 'I')
	{
		nobst = nprt*countwell;	//number of observation points and times
		nptoi = npoi*ntoi;	//number of points and times of interest
		//Jdiv.resize(dI.size()*nobst);
	}

	for (int ii = 0; ii < nz; ii++)
	{
		/*
		clock_t t10 = clock();
		*/
		Jdv.clear();
		for (int jj = 0; jj < npc*npc; jj++)
		{
			Ar[jj] = Arbuild[(ii+1)*npc*npc + jj];
		}
		red_solver(Pv,Ar,Br,well,wellrate,npc,Nn,nq,obstimes,dtmult,Jdv,nprt,lumpwell,dt,Ps);
		//Calculate sensitivies to changes in hydraulic conductivty
		if (expdesparam == 1)
		{
			Jdiv.clear();
			Jdiv.resize(nprt*countwell);
			for (int jj = 0; jj < nprt*countwell; jj++)
			{
				Jdiv[jj] = (Jdvbase[jj] - Jdv[jj])/(permc[ii]*.01);
			}
			Jd.augment_vec(Jdiv);
		}
		/*
		else if (expdespump == 1)	//this part should really go outside the 0:nz loop because we don't care about those values
		{
			for (int jj = 0; jj < nprt*countwell; jj++)
			{
				Jdiv[jj] = (Jdvbase[jj])/dQ;
			}
			Jd.augment_vecmat(Jdiv,(nprt*countwell));
		}
		*/
		/*
		clock_t t11 = clock();
		printf("\nTime for loop: %g seconds\n", (double)(t11 - t10) / CLOCKS_PER_SEC);
		*/
	}
	/*
	clock_t t8 = clock();
	printf("\nTime rest of Jd: %g seconds\n", (double)(t8 - t7) / CLOCKS_PER_SEC);
	*/
	//If requesting T (A-opt) or E (E-opt)
	//For PSO, minimization so taking 1/ for A,E,D (maximize) and normal for G,I (minimize)
	if (OptCI == 'T' || OptCI == 'E')
	{
		I.multiply('T',Jd,'N',Jd,I,1.0,0.0);
		fvalue = std::abs(1.0/I.norm(OptCI));
	}
	else if (OptCI == 'G' || OptCI == 'I')
	{
		mat2d U,S,Vt;
		Jd.svd(U,S,Vt);
		std::vector<double> Sv = S.get_vec(0,S.rows(),0,1);
		//d_x_e is the variance of location x due to design e;
		std::vector<double> d_x_e(nptoi,0.0), xV, x(nz);
		for (int ii = 0; ii < nptoi; ii++)
		{
			for (int jj = 0; jj < nz; jj++)
			{
				x[jj] = dI[jj*nptoi + ii];
			}
			xV = Vt.matvec_mult(1.0,'N',x,0.0,x);
			for (unsigned int jj = 0; jj < xV.size();jj++)
			{
				d_x_e[ii] += (xV[jj]*xV[jj]/(Sv[jj]*Sv[jj]));
			}
		}
		if (OptCI == 'G')
		{
			fvalue = std::abs(*std::minmax_element(d_x_e.begin(),d_x_e.end()).second);
		}
		else
		{
			fvalue = std::abs(std::accumulate(d_x_e.begin(),d_x_e.end(),0.0)/nptoi);
		}
	
	}
	else if (OptCI == 'D')
	{
		I.multiply('T',Jd,'N',Jd,I,1.0,0.0);
		I.multiply_const(multconst);
		fvalue = std::abs(1.0/I.determinant());
	}
	else
	{
		std::cout << "No Optimality Criterion coded to handle OptCI = " << OptCI << " at this time";
		std::exit(1);
	}
	/*
	clock_t t9 = clock();
	printf("\nTime calculate fvalue: %g seconds\n", (double)(t9 - t8) / CLOCKS_PER_SEC);
	*/
	/*
	///////////test code
	//if (ncalls >= 1)
	//if (z > 0)
	{
	std::cout << '\n' << "Test output" << std::endl;
	
	int currentzero = 0;
	for (int ii = 0; ii < genome.length(); ii++)
	{
		std::cout << genome.gene(ii) << '\t';
		
		if (genome.gene(ii) != -1 && pickgenome == 1)
		{
			std::vector<double> pvi = P.get_vec(obswells[genome.gene(ii) + currentzero].node-1,
			obswells[genome.gene(ii) + currentzero].node,0,npc);
	//			std::cout << '\t' << obswells[genome.gene(ii) + currentzero].node << '\t' << obswells[genome.gene(ii) + currentzero].zone;
				std::cout << obswells[genome.gene(ii) + currentzero].node << '\t';
		}
		currentzero += ozonec[ii];
		
	//		std::cout << "currentzero = " << currentzero << std::endl;
	}
		std::cout << std::endl;
	Jd.print();
	Jd.save_to_file("Jd.txt","no","Jd");
	std::cout << std::endl;
	I.print();
	
	std::cout << '\n' << "z = " << z << std::endl;
	std::exit(1);
	}
	
	//std::exit(1);
	//////////end test
	*/

	if (fvalue < bestfvalue){bestfvalue = fvalue;}
	//clock_t t5 = clock();
	//printf("\nTime to call Obj: %g seconds\n", (double)(t5 - t4) / CLOCKS_PER_SEC);
	//std::exit(1);
	return std::abs(fvalue);
	
}