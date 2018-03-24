//Standard headers
#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
//#include <cmath>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>
#include <limits>
//#include <time.h>

//My headers
#include "owellstruct.h"
#include "matrix.h"

//Global variables
mat2d P;
std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, permc, Pv, wellrate, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;

std::vector<double> lb1,ub1,dI,lbi,ubi,first_call;	//lb1 is the lower bound of the optimization, lbi is the lower bound of the intial particle
double lb2, ub2, multconst, bestfvalue;

std::vector<double> diagh;
std::vector<int> poi;
int npoi,ntoi;
int nopcalls;

bool runlongernow;

std::vector<int> topol, JA, triang, trija, well, ozonec;
std::vector<owell> obswells;			//all potential obsrvation nodes
int ntri, irad, nterm, npc, Nn, nprt, lumpwell, pickgenome, genomeoption, ncalls, nowell, nobszones,
nq, nz, expdesparam, expdespump;	//nz: total number of hydrologic zones in the aquifer
double dtmult, dt, dH, dQ;
char OptCI, OptCG;	//OptCI = Individual Optimization Criterion (Trace, Eigen, SE) OptCG = Global Opt. Criterion (Robust, Worst Case)
float z, zmin;

//Function prototypes
void read_fort(int &ntri, int &irad, int &nterm, std::vector<int> &triang2, 
std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, 
std::vector<double> &elstor2, std::vector<double> &spess2, 
std::vector<double> &area2, std::vector<double> &arear2, 
std::vector<double> &bi2, std::vector<double> &ci2, 
std::vector<int> &ja2, std::vector<int> &topol2);

void read0(int& nowell, int &npl, int &startk, std::vector<double> &klim, std::vector<double> &obstimes, 
std::vector<int> &obszonelist, std::vector<owell> &obswells, 
std::vector<int> &ozonec, double &dt, double &dtmult, int &lumpwell, double &dH, double &dQ,
int &pickgenome, int &genomeoption, int &nobszones, int &nobstime, int &expdesparam, int &expdespump,
char &OptCI, char &OptCG, double &close, int &restart, int &restartfreq, double &dog, float &nconv,
bool &runlonger);

void read(int &npc, int &Nn, int &nq, std::vector<int> &well, 
std::vector<double> &wellweight, mat2d &P, int &nz);

void read_restart(int &nowell, int &ncalls, std::vector<robustowell> &rowells,
int nz, int &columns);

void read_interesting_points(std::vector<int> &poi, std::vector<double> &toi, int &npoi, int &ntoi);

void write_restart(int nowell, const std::vector<robustowell> &rowells, 
int ncalls, int kleft);

void write_final(int nowell, std::vector<robustowell> &rowells, int ncalls,
int nz, std::vector<owell> &allowell);

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
/*
void wellsame(std::vector<robustowell> &rowells, std::vector<owell> &copt,
double &close, double score, GARealGenome &g, int combindex);
*/

void wellsame(std::vector<robustowell> &rowells, std::vector<owell> &copt,
	double &close, double score, std::vector<double> &particle, int combindex);


void countowells(const std::vector<robustowell> &rowells, std::vector<owell> &allowell,
double close, int pickgenome);

//PSO wrapper
void Optimize_wrapper(std::vector<double> &bestparticle, float nconv);
double Objective_wrapper(std::vector<double> &particle);

int main()
{
	runlongernow = false;
	nopcalls = 0;
	int nl,startk, readnowell, restart, restartfreq;
	double close, dog;
	float nconv;
	bool runlonger;
	std::vector<owell> obswells2;
	//Reduced model parameters
	//mat2d perm, permcheck, permpicked;
	mat2d perm;
	
	//Set the minimum value z can take
	zmin = std::numeric_limits<float>::min();
	
	//Exp design parameters
	std::vector<int> obszonelist, ozonec2;		//list of observation zones, number of nodes in each observation zone
	std::vector<double> klim, toi;			//observation times, hydraulic conducitivity limits; interesting times (G and I opt)
	std::vector<robustowell> rowells,rowellslong;	//stores robust solutions; rowellslong stores the wells if we want to run the optimization longer
	//Read in parameters needed for fortan functions
	read_fort(ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,spess,area,arear,bi,ci,JA,topol);
	//Read in parameters needed for GA search
	
	//std::cout << "read0" << std::endl;
	read0(readnowell,nl,startk, klim, obstimes, obszonelist, obswells2, ozonec2, dt, dtmult, lumpwell, dH, dQ, pickgenome,
	genomeoption, nobszones, nprt, expdesparam, expdespump, OptCI, OptCG, close, restart, restartfreq, dog, nconv,
	runlonger);
	
	//multconst = floor(pow(10,20.0/(readnowell*obstimes.size())));
	multconst = 1.0;
	
	
	
	//Read Points of interest
	read_interesting_points(poi,toi,npoi,ntoi);
	
	//obswells2: (vector) structures containing obs wells parameters (node, zone, etc.), is sorted by zone
	//obszonelist: (vector) observation zones, sorted low to high
	//ozonec2: (vector) count of nodes in each observation zone
	
	
	//Follows is code to overlay a course gird and only pick observation wells on that grid
	//std::vector<int> ozonec(ozonec2.size(),0);
	ozonec.resize(ozonec2.size());
	obswells.clear();
	//std::cout << "max min" << std::endl;
	//find min, max x and y coordinates of each zone
	std::vector<double> minxy(2*obszonelist.size(),10E10), maxxy(2*obszonelist.size(),0.0);
	int zonecounter = 0;
	int currentzero = 0;
	for (unsigned int ii = 0; ii < obszonelist.size(); ii++)
	{
		for (int jj = 0; jj < ozonec2[ii]; jj++)
		{
			if (obswells2[currentzero + jj].x < minxy[2*ii]){minxy[2*ii] = obswells2[currentzero + jj].x;}
			if (obswells2[currentzero + jj].y < minxy[2*ii+1]){minxy[2*ii+1] = obswells2[currentzero + jj].y;}
			if (obswells2[currentzero + jj].x > maxxy[2*ii]){maxxy[2*ii] = obswells2[currentzero + jj].x;}
			if (obswells2[currentzero + jj].y > maxxy[2*ii+1]){maxxy[2*ii+1] = obswells2[currentzero + jj].y;}
			obswells2[currentzero + jj].zonecount = zonecounter;
			//std::cout << obswells2[currentzero + jj].zonecount << '\t';
		}
		currentzero += ozonec2[ii];
		zonecounter++;
	}
	
	//std::cout << std::endl;
	
	//std::cout << "grid 00" << std::endl;
	//find grid "(0,0)" for each zone
	std::vector<double> xy0(2*obszonelist.size(),0);
	for (unsigned int ii = 0; ii < obszonelist.size(); ii++)
	{
		int stepx0 = minxy[2*ii]/dog;
		int stepy0 = minxy[2*ii + 1]/dog;
		xy0[2*ii] = stepx0*dog;
		xy0[2*ii+1] = stepy0*dog;
		//std::cout << xy0[2*ii] << '\t' << xy0[2*ii+1] << std::endl;
	}
	//std::vector<owell> obswells;
	std::vector<owell> obswellstemp;
	
	//assign the closest SAT2D node to course grid node
	currentzero = 0;
	for (unsigned int ii = 0; ii < obszonelist.size(); ii++)
	{
		obswellstemp.clear();
		int countadd = 0;
		obswellstemp.insert(obswellstemp.begin(),obswells2.begin() + currentzero, obswells2.begin() + currentzero + ozonec2[ii]);
		int stepsx = ((maxxy[2*ii] - xy0[2*ii])/dog);
		int stepsy = ((maxxy[2*ii+1] - xy0[2*ii+1])/dog);
		if (stepsx < 1) {stepsx = 1;}
		if (stepsy < 1) {stepsy = 1;}
		//std::cout << stepsx << '\t' << stepsy << std::endl;
		for (int xx = 0; xx < stepsx; xx++)
		{
			double cgridx = xx*dog + xy0[2*ii];
			for (int yy = 0; yy < stepsy; yy++)
			{
				double cgridy = yy*dog + xy0[2*ii+1];
				std::vector<double> disttoogrid;
				double mindist = 10E10;
				owell minowell;
				//minowell.x = 0.0;
				int erase;
				for (unsigned int jj = 0; jj < obswellstemp.size(); jj++)
				{
					double xdist = abs(cgridx - obswellstemp[jj].x);
					double ydist = abs(cgridy - obswellstemp[jj].y);
					double dist = sqrt(xdist*xdist + ydist*ydist);
					if (xdist <= dog && ydist <= dog && dist < mindist)
					{
						minowell = obswellstemp[jj];
						erase = jj;
						mindist = dist;
					}
				}
				if (minowell.zone != 0.0)
				{
					obswells.push_back(minowell);
					countadd++;
					obswellstemp.erase(obswellstemp.begin() +  erase, obswellstemp.begin() + erase + 1);
				}
			}
		}
		currentzero += ozonec2[ii];
		ozonec[ii] = countadd;
	}
	std::cout << "obswells.size() = " << obswells.size() << std::endl;
	//Read in paramaters needed for reduced model
	read(npc,Nn,nq,well,wellweight,P,nz);
	//formulate the upper,middle,lower combinations
	std::vector<double> lmh(nl,0);
	lmh[0] = klim[1];
	lmh[nl-1] = klim[0];
	for (int ii = 1; ii < nl-1; ii++) {lmh[ii] = lmh[0] + ii*(lmh[nl-1] - lmh[0])/(nl-1);}
	//store of all perm combos		
	perm.resize(nz,pow(nl,nz));
	for (int ii = 0; ii < nz; ii++)
	{
		for (int jj = 0; jj < pow(nl,nz)/pow(nl,ii+1); jj++)
		{
			for (int kk = 0; kk < nl; kk++)
			{
				for (int ll = 0; ll < pow(nl,ii); ll++)
				{
					perm.set_value(ii,ll + kk*pow(nl,ii) + jj*pow(nl,ii+1),lmh[kk]);
				}
			}
		}
	}
	int ncomb = perm.columns();
	//Projection matrix in vector form
	Pv = P.get_vec(0,Nn,0,npc);
	//set all wellrate elements to dQ
	wellrate.resize(nq);
	for (int ii = 0; ii < nq; ii++) {wellrate[ii] = dQ;}
	
	/*
	1. Create r_baseline vector, r_i then ds_i = P(r_i - r_baseline), this forms on column of Jd
		a. pass out r from red_solver
	2. Form Jd from columns
	3. I = Jd'*Jd -> norm(I)
	*/
	
	////////////Iterations through the nowell (1:nobszones)
	int gatotalcalls = 1;
	if (readnowell <= 0) //If Number of observation wells == -1 then loop for all wells up to number of observation wells
	{
		gatotalcalls = nobszones;
		nowell = 0;
	}	//otherwise just run for one number of observation wells
	int begincomb = 0;

	std::vector<owell> obswells3 = obswells;


	if (restart == 1)
	{
		int kleft;
		std::cout << "Attempting Restart" << std::endl;
		read_restart(readnowell, ncalls, rowells, nz, kleft);
		gatotalcalls = 1;
		begincomb = ncomb - kleft;
		/////////////fixing x, y values//////////////////
		//std::cout << "setting parameters from restart file" << std::endl;
		for (unsigned int kk = 0; kk < rowells.size(); kk++)
		{
			std::vector<owell> copt;	//current optimal solution
			currentzero = 0;
			for (unsigned int ii = 0; ii < rowells[kk].genome.size(); ii++)
			{
				if (rowells[kk].genome[ii] != -1)
				{
					copt.push_back(obswells[rowells[kk].genome[ii] + currentzero]);
				}
				if (pickgenome == 1)
				{
					currentzero += ozonec[ii];
				}
				//std::cout << "gene " << rowells[kk].genome[ii] << " pickgenome " << pickgenome << std::endl;
				/*
				if (rowells[kk].genome[ii] != -1 && pickgenome == 1)	//Store solution from genome type 1
				{
					copt.push_back(obswells[rowells[kk].genome[ii] + currentzero]);
					//std::cout << genome.gene(ii) << '\t' << copt[copt.size()-1].node << '\t' << copt[copt.size()-1].zone << std::endl;
				}
				else
				{
					copt.push_back(obswells[rowells[kk].genome[ii]]);
					//std::cout << "gene " << rowells[kk].genome[ii] << "\tnode " << obswells[rowells[kk].genome[ii]].node << "\tzone " << obswells[rowells[kk].genome[ii]].zone << std::endl;
				}
				currentzero += ozonec[ii];
				*/
			}
			rowells[kk].optowells = copt;
		}
	}
	std::vector<double> Ar(npc*npc), Br(npc*npc), dAr(0, 0);

	permc.clear();
	permc.resize(nz, 0);
	for (int ii = 0; ii < nz; ii++)
	{
		permc[ii] = 1;
		redmat(Pv, Ar, Br, topol, JA, ntri, irad, nterm, triang, trija, lmass, xc, yc, elstor,
			spess, area, arear, bi, ci, npc, Nn, permc);
		//dAr contains quick reduced matrices
		dAr.insert(dAr.end(), Ar.begin(), Ar.end());
		permc[ii] = 0;
	}
	permc.clear();

	////////////1/6/13 End major modifications

	///////3/31/14 Edit here to activate startk

	std::vector<int> sequence0(ncomb, 0), sequence;
	for (int ii = 0; ii < ncomb; ii++) { sequence0[ii] = ii; }
	//basically cutting the deck so startk is the first k run then the rest run in sequence
	sequence.insert(sequence.end(), sequence0.begin() + startk, sequence0.end());
	sequence.insert(sequence.end(), sequence0.begin(), sequence0.begin() + startk);

	///////3/31/14 End modification

	for (int gacurrentcall = 0; gacurrentcall < gatotalcalls; gacurrentcall++)
	{
		if (restart != 1)
		{
			ncalls = 0;
			Jdv.clear();
			rowells.clear();
		}
		//std::vector<double> bestparticle;

		if (readnowell <= 0) { nowell++; }
		else { nowell = readnowell; }
		if (pickgenome == 1)	//The elements in this genome can only pick from a specific zone (or turn off that zone)
		{
			lb1.resize(nobszones, -1.0);	//-1 to be able to turn off zone
			ub1.resize(nobszones, 0.0);
			for (int ii = 0; ii < nobszones; ii++)
			{
				ub1[ii] = ozonec[ii] - 1;	//-1 to account for C++ counting
			}
		}
		else if (pickgenome == 2)	//The elements in this genome can pick from any potential observation location
		{
			lb1.resize(nowell, 0.0);
			ub1.resize(nowell, (obswells.size() - 1));
		}
		else
		{
			std::cout << "Unknown genome type" << std::endl;
			std::exit(1);
		}

		std::vector<double> bestparticle(lb1.size(),0.0);

		std::vector<double> Piv;	//Pi in vector form
		mat2d Pi;	//Projection matrix for interesting points
		if (OptCI == 'G' || OptCI == 'I')
		{
			if (npoi == 0)
			{
				std::cout << "Setting interesting points as all possible observation points" << '\n';
				poi.resize(obswells.size());
				for (unsigned int ii = 0; ii < obswells.size(); ii++)
				{
					poi[ii] = obswells[ii].node;
				}
				npoi = poi.size();
			}
			if (ntoi == 0)
			{
				std::cout << "Setting interesting times as observation times" << std::endl;
				toi = obstimes;
				ntoi = nprt;
			}
			dI.resize(poi.size()*toi.size()*nz);
			for (unsigned int ii = 0; ii < poi.size(); ii++)
			{
				std::vector<double> piv = P.get_vec(poi[ii] - 1, poi[ii], 0, npc);
				Piv.insert(Piv.end(), piv.begin(), piv.end());
			}
		}
		std::cout << "Create P interesting matrix" << std::endl;
		Pi.vec_to_mat(Piv, poi.size(), npc, 1);

		//start loop over k combinations here
		for (int kk = begincomb; kk < ncomb; kk++)
		{
			Arbuild.clear();
			Brbuild.clear();
			permc = perm.get_vec(0, nz, sequence[kk], sequence[kk] + 1);
			for (int ii = 0; ii < nz; ii++)
			{
				std::cout << permc[ii] << '\t';
			}
			std::cout << std::endl;

			std::vector<double> permstore = permc;
			std::cout << "Building Reduced Matrices" << std::endl;
			std::fill(Ar.begin(), Ar.end(), 0.0);

			for (int ii = 0; ii < nz; ii++)
			{
				for (int jj = 0; jj < npc*npc; jj++)
				{
					Ar[jj] += dAr[jj + ii*npc*npc] * permc[ii];
				}
			}
			Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
			Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());

			//build perturbed reduced matrices
			for (int ii = 0; ii < nz; ii++)
			{
				permc[ii] *= 1.01;
				std::fill(Ar.begin(), Ar.end(), 0.0);

				for (int kk = 0; kk < nz; kk++)
				{
					for (int jj = 0; jj < npc*npc; jj++)
					{
						Ar[jj] += dAr[jj + kk*npc*npc] * permc[kk];
					}
				}
				//////////// 1/6/13 stop modify code here
				Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
				Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());
				permc[ii] = permstore[ii];
			}

			std::cout << "Running PSO" << std::endl;
			//Run PSO
			std::vector<double> di, dibase;	//vectors for changes in interesting points


			if (OptCI == 'G' || OptCI == 'I')
			{
				std::vector<double> Ar(npc*npc), Br(npc*npc);
				for (int ii = 0; ii < npc*npc; ii++)
				{
					Ar[ii] = Arbuild[ii];
					Br[ii] = Brbuild[ii];
				}
				red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, dibase, toi.size(), lumpwell, dt, Pi);

				for (int ii = 0; ii < nz; ii++)
				{
					di.clear();
					for (int jj = 0; jj < npc*npc; jj++)
					{
						Ar[jj] = Arbuild[(ii + 1)*npc*npc + jj];
					}
					red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, di, toi.size(), lumpwell, dt, Pi);
					//Calculate sensitivies to changes in hydraulic conductivty
					if (expdesparam == 1)
					{
						for (unsigned int jj = 0; jj < dibase.size(); jj++)
						{
							dI[ii*dibase.size() + jj] = (dibase[jj] - di[jj]) / (permc[ii] * 0.01);
						}
					}
					else
					{
						std::cout << "Can only run parameter at this point";
						std::exit(1);
					}
				}
			}
			//Run PSO optimization
			//std::vector<double> bestparticle(lb1.size(),0.0);
			//bestparticle.resize(lb1.size(), -1.0);
			std::fill(bestparticle.begin(), bestparticle.end(), -1.0);
			Optimize_wrapper(bestparticle,nconv);
			std::vector<owell> copt;	//current optimal solution
			currentzero = 0;
			for (unsigned int ii = 0; ii < bestparticle.size(); ii++)
			{
				if (pickgenome == 1 && bestparticle[ii] != -1)
				{
					copt.push_back(obswells[bestparticle[ii] + currentzero]);
					std::cout << bestparticle[ii] << '\t' << copt[copt.size() - 1].node << '\t' <<
					copt[copt.size() - 1].zone << std::endl;
				}
				else if (pickgenome == 2)
				{
					copt.push_back(obswells[bestparticle[ii]]);
					std::cout  << bestparticle[ii] << '\t' << copt[copt.size()-1].node << '\t' << 
					copt[copt.size()-1].zone << std::endl;
				}
				currentzero += ozonec[ii];
			}
			wellsame(rowells, copt, close, 0.0, bestparticle, sequence[kk]);

			//if (kk > 0 && kk % restartfreq == 0)
			if (kk % restartfreq == 0)
			{
				std::cout << "Saving Restart Point" << std::endl;
				int kleft = ncomb - (kk + 1);
				write_restart(nowell, rowells, ncalls, kleft);
			}
		}
		//Save final solution before running worst case scenario
		write_restart(nowell, rowells, ncalls, 0);
		
		//edit 1/19/17 add runlonger command
		if (runlonger == true)
		{
			rowellslong.clear();
			
			std::cout << "Running optimization longer" << std::endl;
			lbi.resize(lb1.size(),0);
			ubi.resize(ub1.size(),0);
			for (unsigned int kk = 0; kk < rowells.size(); kk++)
			{
				for (unsigned int ii = 0; ii < lb1.size(); ii++)
				{
					lbi[ii] = rowells[kk].genome[ii];
				}
				ubi = lbi;
				for (unsigned int ii = 0; ii < rowells[kk].paramindex.size(); ii++)
				{
					runlongernow = true;
					first_call.resize(bestparticle.size(), 0);
					Arbuild.clear();
					Brbuild.clear();
					permc = perm.get_vec(0,nz,rowells[kk].paramindex[ii],rowells[kk].paramindex[ii]+1);
					//rowells[kk].paramindex.erase(ii);//not sure how this will affect the iterations
					for (int ii = 0; ii < nz; ii++){std::cout << permc[ii] << '\t';}
					std::cout << std::endl;
					
					std::vector<double> permstore = permc;
					std::cout << "Building Reduced Matrices" << std::endl;
					std::fill(Ar.begin(), Ar.end(), 0.0);
		
					for (int ii = 0; ii < nz; ii++)
					{
						for (int jj = 0; jj < npc*npc; jj++)
						{
							Ar[jj] += dAr[jj + ii*npc*npc] * permc[ii];
						}
					}
					Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
					Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());
		
					//build perturbed reduced matrices
					for (int ii = 0; ii < nz; ii++)
					{
						permc[ii] *= 1.01;
						std::fill(Ar.begin(), Ar.end(), 0.0);
		
						for (int kk = 0; kk < nz; kk++)
						{
							for (int jj = 0; jj < npc*npc; jj++)
							{
								Ar[jj] += dAr[jj + kk*npc*npc] * permc[kk];
							}
						}
						//////////// 1/6/13 stop modify code here
						Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
						Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());
						permc[ii] = permstore[ii];
					}
		
					std::cout << "Running PSO Longer" << std::endl;
					//Run PSO
					std::vector<double> di, dibase;	//vectors for changes in interesting points
					if (OptCI == 'G' || OptCI == 'I')
					{
						std::vector<double> Ar(npc*npc), Br(npc*npc);
						for (int ii = 0; ii < npc*npc; ii++)
						{
							Ar[ii] = Arbuild[ii];
							Br[ii] = Brbuild[ii];
						}
						red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, dibase, toi.size(), lumpwell, dt, Pi);
		
						for (int ii = 0; ii < nz; ii++)
						{
							di.clear();
							for (int jj = 0; jj < npc*npc; jj++)
							{
								Ar[jj] = Arbuild[(ii + 1)*npc*npc + jj];
							}
							red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, di, toi.size(), lumpwell, dt, Pi);
							//Calculate sensitivies to changes in hydraulic conductivty
							if (expdesparam == 1)
							{
								for (unsigned int jj = 0; jj < dibase.size(); jj++)
								{
									dI[ii*dibase.size() + jj] = (dibase[jj] - di[jj]) / (permc[ii] * 0.01);
								}
							}
							else
							{
								std::cout << "Can only run parameter at this point";
								std::exit(1);
							}
						}
					}
					//Run PSO optimization
					std::fill(bestparticle.begin(), bestparticle.end(), -1.0);
					Optimize_wrapper(bestparticle,nconv);
					
					if (std::equal(bestparticle.begin(), bestparticle.end(), &first_call[0])) //check to see if bestparticle matches the first call
					{
						std::cout << "No change for parameter #" << rowells[kk].paramindex[ii] << std::endl;
						bestparticle = lbi;
					}
					else
					{
						std::cout << "Update for parameter #" << rowells[kk].paramindex[ii] << std::endl;
					}

					std::vector<owell> copt;	//current optimal solution
					currentzero = 0;
					for (unsigned int ii = 0; ii < bestparticle.size(); ii++)
					{
						if (pickgenome == 1 && bestparticle[ii] != -1)
						{
							copt.push_back(obswells[bestparticle[ii] + currentzero]);
							std::cout << bestparticle[ii] << '\t' << copt[copt.size() - 1].node << '\t' <<
							copt[copt.size() - 1].zone << std::endl;
						}
						else if (pickgenome == 2)
						{
							copt.push_back(obswells[bestparticle[ii]]);
							std::cout  << bestparticle[ii] << '\t' << copt[copt.size()-1].node << '\t' << 
							copt[copt.size()-1].zone << std::endl;
						}
						currentzero += ozonec[ii];
					}
					wellsame(rowellslong, copt, close, 0.0, bestparticle, rowells[kk].paramindex[ii]);
					if (kk % restartfreq == 0)
					{
						std::cout << "Saving Restart Point" << std::endl;
						int kleft = ncomb - (kk + 1);
						write_restart(nowell, rowellslong, ncalls, 0);
					}
				}
				
			}
			rowells = rowellslong;
			write_restart(nowell, rowells, ncalls, 0);
		}
		//end edit
		
		

		std::vector<owell> allowell(1);
		countowells(rowells, allowell, close, pickgenome);
		//1/26/17 remove checking for high frequency wells or network
		/*
		std::vector<owell> copt;
		if (pickgenome == 1)
		{
			copt.insert(copt.begin(), allowell.begin(), allowell.begin() + nowell);
			std::fill(bestparticle.begin(), bestparticle.end(), -1.0);
			for (unsigned int ii = 0; ii < copt.size(); ii++)
			{
				bestparticle[copt[ii].zonecount] = copt[ii].gene;
				copt[ii].highfreq = true;
			}
			wellsame(rowells, copt, close, 0.0, bestparticle, 0);
		}
		else
		{
			std::vector<int> zonesused;
			copt.push_back(allowell[0]);
			zonesused.push_back(obswells[allowell[0].gene].zone);
			for (unsigned int ii = 1; ii < allowell.size(); ii++)
			{
				int checkzone = allowell[ii].zone;
				if (copt.size() == (unsigned int)nowell){break;}
				if(std::find(zonesused.begin(),zonesused.end(),checkzone) == zonesused.end())
				{
					copt.push_back(allowell[ii]);
					zonesused.push_back(allowell[ii].zone);
				}	
			}
			if (copt.size() == (unsigned int)nowell)
			{
				for (int ii = 0; ii < nowell; ii++)
				{
					bestparticle[ii] = copt[ii].gene;
					copt[ii].highfreq = true;
				}
				wellsame(rowells,copt,close,0.0,bestparticle,0);
			}
		}

		std::sort(rowells.begin(), rowells.end(), by_count_network());	//sort rowells by frequency of network occurance
		rowells[0].highfreqnetwork = true;

		for (unsigned int ii = 0; ii < rowells.size(); ii++)			//move high frequency wells to front of list
		{
			if (rowells[ii].highfreqwells == true)
			{
				robustowell temp = rowells[ii];
				rowells.erase(rowells.begin() + ii, rowells.begin() + ii + 1);
				rowells.insert(rowells.begin(), temp);
				break;
			}
		}
		*/
		//For each perm combination
		std::cout << "Running Worst Case Scenario Test" << std::endl;
		for (int jj = 0; jj < ncomb; jj++)
		{
			std::cout << '\r' << "Parameter Combintaion: " << jj << "                     ";
			//std::cout << jj << std::endl;
			permc = perm.get_vec(0, nz, jj, jj + 1);
			std::vector<double> permstore = permc;
			Arbuild.clear();
			Brbuild.clear();
			//build baseline reduced matrices

			//Ar.clear();
			//Ar.resize(npc*npc);
			std::fill(Ar.begin(), Ar.end(), 0.0);

			//std::vector<double> Ar,Br;
			//std::cout << "crash1" << std::endl;
			////////// 1/6/14 modify code here for quick red matrix
			/*
			redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
			spess,area,arear,bi,ci,npc,Nn,permc);
			*/


			for (int ii = 0; ii < nz; ii++)
			{
				for (int jj = 0; jj < npc*npc; jj++)
				{
					Ar[jj] += dAr[jj + ii*npc*npc] * permc[ii];
				}
			}

			////////// 1/6/14 stop modify code here

			Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
			Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());
			//build perturbed reduced matrices
			for (int kk = 0; kk < nz; kk++)
			{
				permc[kk] *= 1.01;
				std::fill(Ar.begin(), Ar.end(), 0.0);
				for (int ii = 0; ii < nz; ii++)
				{
					for (int jj = 0; jj < npc*npc; jj++)
					{
						Ar[jj] += dAr[jj + ii*npc*npc] * permc[ii];
					}
				}
				Arbuild.insert(Arbuild.end(), Ar.begin(), Ar.end());
				Brbuild.insert(Brbuild.end(), Br.begin(), Br.end());
				permc[kk] = permstore[kk];
			}

			std::vector<double> di, dibase;
			if (OptCI == 'G' || OptCI == 'I')
			{
				std::vector<double> Ar(npc*npc), Br(npc*npc);
				for (int ii = 0; ii < npc*npc; ii++)
				{
					Ar[ii] = Arbuild[ii];
					Br[ii] = Brbuild[ii];
				}
				red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, dibase, toi.size(), lumpwell, dt, Pi);
				for (int ii = 0; ii < nz; ii++)
				{
					di.clear();
					for (int jj = 0; jj < npc*npc; jj++)
					{
						Ar[jj] = Arbuild[(ii + 1)*npc*npc + jj];
					}
					red_solver(Pv, Ar, Br, well, wellrate, npc, Nn, nq, toi, dtmult, di, toi.size(), lumpwell, dt, Pi);
					//Calculate sensitivies to changes in hydraulic conductivty
					if (expdesparam == 1)
					{
						for (unsigned int jj = 0; jj < dibase.size(); jj++)
						{
							//Follwing line is for Jd = [dO/dK; dI/dK] optimized over choice of O but score based on I
							dI[ii*dibase.size() + jj] = (dibase[jj] - di[jj]) / (permc[ii] * 0.01);
						}
					}
					else
					{
						std::cout << "Can only run parameter at this point";
						std::exit(1);
					}
				}
			}
			//check each robust well combination
			for (unsigned int ii = 0; ii < rowells.size(); ii++)
				//for (unsigned int ii = 0; ii < 2; ii++)
			{
				if (jj == 0)
				{
					rowells[ii].maxscore.clear();//vector of all the scores
					rowells[ii].minscore = 10E100;
				}

				for (unsigned int kk = 0; kk < bestparticle.size(); kk++)
				{
					bestparticle[kk] = rowells[ii].genome[kk];
				}
				double fvalue = 1.0/Objective_wrapper(bestparticle); //turn min problem scores into max scores so can make the same comparisons as the GA code
				/*
				if (OptCI == 'D' || OptCI == 'E' || OptCI == 'T')
				{
					fvalue = 1.0 / Objective_wrapper(bestparticle);	//turn it back into the true value
				}
				else	//for G or I
				{
					fvalue = Objective_wrapper(bestparticle);
				}
				*/
				rowells[ii].maxscore.push_back(fvalue);
				if (fvalue <= rowells[ii].minscore)
				{
					rowells[ii].minscore = fvalue;
					rowells[ii].worstcaseK = permc;
				}
			}
		}

		double maxworst = -10E100;
		int maxworstindex;
		for (unsigned int ii = 0; ii < rowells.size(); ii++)
		{
			if (rowells[ii].minscore > maxworst)
			{
				maxworst = rowells[ii].minscore;
				maxworstindex = ii;
			}
		}
		rowells[maxworstindex].robustdesign = true;
		robustowell temp = rowells[maxworstindex];
		rowells.erase(rowells.begin() + maxworstindex, rowells.begin() + maxworstindex + 1);
		rowells.insert(rowells.begin(), temp);
		std::cout << '\n' << std::endl;
		write_final(nowell, rowells, ncalls, nz, allowell);
	}
}