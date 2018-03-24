
#include "matrix.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <random>
#include <ctime>
#include <chrono>

void read_fort(int &ntri, int &irad, int &nterm, std::vector<int> &triang2, 
std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, 
std::vector<double> &elstor2, std::vector<double> &spess2, 
std::vector<double> &area2, std::vector<double> &arear2, 
std::vector<double> &bi2, std::vector<double> &ci2, 
std::vector<int> &ja2, std::vector<int> &topol2);

void read(int &npc, int &Nn, int &nq, std::vector<int> &well, 
std::vector<double> &wellweight, mat2d &P, int &nz);

void redmat(std::vector<double> &Pv2, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &topol, std::vector<int> &JA, int &nt, int &irad, int &nterm,
std::vector<int> &triang2, std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, std::vector<double> &elstor2, 
std::vector<double> &spess2, std::vector<double> &area2, 
std::vector<double> &arear2, std::vector<double> &bi2, 
std::vector<double> &ci2, int &npc, int &Nn, std::vector<double> &permc);

int main(void)
{
mat2d P,sshots;
std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, Pv, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;
std::vector<int> topol, JA, triang, trija, well, ozonec;
//std::vector<owell> obswells;			//all potential obsrvation nodes
int ntri, irad, nterm, npc, Nn, nq, nz;

read_fort(ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,spess,area,arear,bi,ci,JA,topol);
read (npc,Nn,nq,well,wellweight,P,nz);

std::vector<double> permc(nz,1), Ar(npc*npc), Br(npc*npc), dAr(0,0), Arq(npc*npc), Ar0(npc*npc), dArv(npc*npc);
std::vector<double> permstore(0,0);
unsigned seed = 5;
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> distribution(0.1,20.0);
Pv = P.get_vec(0,Nn,0,npc);



redmat(Pv,Ar0,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
spess,area,arear,bi,ci,npc,Nn,permc);



for (int ii = 0; ii < nz; ii++)
{	
	permc[ii] = 2;
	redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
	spess,area,arear,bi,ci,npc,Nn,permc);
	for (int jj = 0; jj < npc*npc; jj++)
	{
		dArv[jj] = Ar[jj] - Ar0[jj];
	}
	dAr.insert(dAr.end(),dArv.begin(),dArv.end());
	permc[ii] = 1;
}

//12/6/14 Adding modifcations to the code to run different interface norms




for (int checkii = 0; checkii < 1; checkii++)
{
	for (int jj = 0; jj < nz; jj++)
	{
		permc[jj] = distribution(generator);
		//permstore[jj + ii*nz] = permc[jj];
		permstore.insert(permstore.end(),permc.begin(),permc.end());
	}
	Arq.clear();
	Arq.resize(npc*npc);
	
	
	std::clock_t timeredstart, timequickstart;
	double durationq, durationred;
	timequickstart = std::clock();
	//Attempting to use chrono to time
	auto begin = std::chrono::high_resolution_clock::now();
	for (int ii = 0; ii < nz; ii++) 
	{
		for (int jj = 0; jj < npc*npc; jj++)
		{
			Arq[jj] += dAr[jj + ii*npc*npc]*permc[ii];
		}
	}
	//std::cout << "pause here" << std::endl;
	//std::cin.get();
	
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() << "ns" << std::endl;
	
	
	//durationq = (std::clock() - timequickstart)/(double)CLOCKS_PER_SEC;
	durationq = (std::clock() - timequickstart);
	std::cout << "Quick build took: " << durationq << std::endl;
	
	timeredstart = std::clock();
	
	redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
	spess,area,arear,bi,ci,npc,Nn,permc);
	
	durationred = (std::clock() - timeredstart)/(double)CLOCKS_PER_SEC;
	std::cout << "Standard build took: " << durationred << std::endl;
	
	double RMSE;
	
	for (int ii = 0; ii < npc*npc; ii++)
	{
		RMSE += (Ar[ii] - Arq[ii])*(Ar[ii] - Arq[ii]);
	}
	
	RMSE = sqrt(RMSE)/(npc*npc);
	
	for (int ii = 0; ii < nz; ii++)
	{
		std::cout << permc[ii] <<  '\t';
	}
	std::cout << RMSE << std::endl;
}

std::ofstream file("quick_red.txt");
if (file.is_open())
{
	for (int ii = 0; ii < npc; ii++)
	{
		for (int jj = 0; jj < npc; jj++)
		{
			file << Ar[jj + ii*npc] << '\t';
		}
		file << '\n';
	}
	
	file << std::endl;
	
	for (int ii = 0; ii < npc; ii++)
	{
		for (int jj = 0; jj < npc; jj++)
		{
			file << Arq[jj + ii*npc] << '\t';
		}
		file << '\n';
	}
}
file.close();
}