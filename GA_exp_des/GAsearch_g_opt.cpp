//Standard headers
#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <ctime>
//#include <time.h>

//My headers
#include "owellstruct.h"
#include "matrix.h"

//GA headers
#include <ga/ga.h>
#include <ga/std_stream.h>
//Initialize real genome (only do in main header)
#define INSTANTIATE_REAL_GENOME
//Call real GA
#include <ga/GARealGenome.h>
//GA objective function prototype
float Objective(GAGenome &g);

//////////////Check global and local variables, confirm that global variables are not defined locally

//Global variables
mat2d P;
std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, permc, Pv, wellrate, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;

std::vector<double> lb1,ub1;	//vectors for lower and upper bounds of type 1 genomes
double lb2, ub2;				//single values for lower and upper bounds of type 2 genomes

std::vector<double> diagh;

std::vector<int> topol, JA, triang, trija, well, ozonec;
std::vector<owell> obswells;			//all potential obsrvation nodes
int ntri, irad, nterm, npc, Nn, nprt, lumpwell, pickgenome, genomeoption, ncalls, nowell, nobszones,
nq, nz, expdesparam, expdespump;	//nz: total number of hydrologic zones in the aquifer
double dtmult, dt, dH, dQ, convergence;
char OptCI, OptCG;	//OptCI = Individual Optimization Criterion (Trace, Eigen, SE, G, I) OptCG = Global Opt. Criterion (Robust, Worst Case)
float z;
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

void red_solver(mat2d &Pt, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &well, std::vector<double> &wellrate, int &npc, int &Nn, int &nq,
std::vector<double> &ts, mat2d &sshots, mat2d &scaledshots, double &norm, 
double &dtmult, mat2d &r, int nprt, int pickk, int lumpwell, std::vector<double> &dtstore);

void wellsame(std::vector<robustowell> &rowells, std::vector<owell> &copt,
double &close, double score, GARealGenome &g, int combindex);

void countowells(const std::vector<robustowell> &rowells, std::vector<owell> &allowell,
double close);

int main()
{
int nl,startk, readnowell, restart, restartfreq;
double close, dog;
float nconv;
bool runlonger;
std::vector<owell> obswells2;
//Reduced model parameters
//mat2d perm, permcheck, permpicked;
mat2d perm;

//Exp design parameters
std::vector<int> obszonelist, ozonec2;		//list of observation zones, number of nodes in each observation zone
std::vector<double> klim;			//observation times, hydraulic conducitivity limits
std::vector<robustowell> rowells;	//stores robust solutions
//Read in parameters needed for fortan functions
read_fort(ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,spess,area,arear,bi,ci,JA,topol);
//Read in parameters needed for GA search

//std::cout << "read0" << std::endl;
read0(readnowell,nl,startk, klim, obstimes, obszonelist, obswells2, ozonec2, dt, dtmult, lumpwell, dH, dQ, pickgenome,
genomeoption, nobszones, nprt, expdesparam, expdespump, OptCI, OptCG, close, restart, restartfreq, dog, nconv,
runlonger);


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

/*
/////code to print out potential obswell x,y
std::ofstream file("wellxy.txt");
if (file.is_open())
{
	for (unsigned int ii = 0; ii < obswells.size(); ii++)
	{
		file << std::scientific << std::setprecision(13) << obswells[ii].x 
		<< '\t' << obswells[ii].y << std::endl;
	}
	file.close();
}
else
{
	std::cout << "Could not open file to write";
}
*/


//std::exit(1);

//std::cout << "dH = " << dH << std::endl;

/////////test to run "same as" distance funtion
/*
std::vector< std::vector<owell> > optlocs;
std::vector<robustowell> rowells;
std::vector<owell> copt(3), copt2(4);
std::vector<int> counter;

copt[0].x = 10;
copt[0].y = 14;
copt[1].x = 20;
copt[1].y = 5;
copt[2].x = 0;
copt[2].y = 8;
wellsame(rowells,copt,close);
copt[0].x = 0;
copt[0].y = 8;
copt[1].x = 10;
copt[1].y = 14;
copt[2].x = 20;
copt[2].y = 5;
wellsame(rowells,copt,close);
copt[0].x = 1;
copt[0].y = 8;
copt[1].x = 10;
copt[1].y = 14;
copt[2].x = 20;
copt[2].y = 5;
wellsame(rowells,copt,close);
copt[0].x = 0;
copt[0].y = 8;
copt[1].x = 20;
copt[1].y = 14;
copt[2].x = 20;
copt[2].y = 5;
wellsame(rowells,copt,close);
copt[0].x = 1;
copt[0].y = 8;
copt[1].x = 10;
copt[1].y = 14;
copt[2].x = 20;
copt[2].y = 5;
wellsame(rowells,copt,close);
copt2[0].x = 1;
copt2[0].y = 8;
copt2[1].x = 10;
copt2[1].y = 14;
copt2[2].x = 20;
copt2[2].y = 5;
wellsame(rowells,copt2,close);
copt[0].x = 0;
copt[0].y = 8;
copt[1].x = 20;
copt[1].y = 14;
copt[2].x = 20;
copt[2].y = 5;
wellsame(rowells,copt,close);

for (unsigned int ii = 0; ii < rowells.size(); ii++)
{
	std::cout << "Number of bucket " << ii << " encountered is " << rowells[ii].count << std::endl
	<< "Size of bucket is " << rowells[ii].optowells.size() << std::endl;
	for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
	{
		std::cout << rowells[ii].optowells[jj].x << '\t' << rowells[ii].optowells[jj].y << std::endl;
	}
}

std::exit(1);
*/
/////////end test


//std::cout << "expdes = " << expdesparam << '\t' << expdespump << std::endl;



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
//std::vector<int> corresro(ncomb,0);
if (restart == 1) 
{
	int kleft;
	std::cout << "Attempting Restart" << std::endl;
	read_restart(readnowell,ncalls,rowells,nz,kleft);
	gatotalcalls = 1;
	begincomb = ncomb - kleft;
	
	/*
	for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		for (int jj = 0; jj < rowells[ii].count; jj++)
		{
			corresro[rowells[ii].paramindex[jj]] = ii;
		}
	}
	*/
	//for (unsigned int ii = 0; ii < corresro.size(); ii++) {std::cout << corresro[ii] << std::endl;}
	
	
	//std::cout << begincomb << std::endl;
	
	/////////////fixing x, y values//////////////////
	for (unsigned int kk = 0; kk < rowells.size(); kk++)
	{
		std::vector<owell> copt;	//current optimal solution
		currentzero = 0;
		for (unsigned int ii = 0; ii < rowells[kk].genome.size(); ii++)
		{
			if (rowells[kk].genome[ii] != -1)	//Store solution from genome type 1
			{
				copt.push_back(obswells[rowells[kk].genome[ii] + currentzero]);
				//std::cout << genome.gene(ii) << '\t' << copt[copt.size()-1].node << '\t' << copt[copt.size()-1].zone << std::endl;
			}
			currentzero += ozonec[ii];
		}
		rowells[kk].optowells.resize(copt.size());
		for (unsigned int jj = 0; jj < copt.size(); jj++)
		{
			rowells[kk].optowells[jj].x = copt[jj].x;
			rowells[kk].optowells[jj].y = copt[jj].y;
			rowells[kk].optowells[jj].zone = copt[jj].zone;
			rowells[kk].optowells[jj].node = copt[jj].node;
			rowells[kk].optowells[jj].zonecount = copt[jj].zonecount;
		}
	}
	
	//write_restart(nowell,rowells,ncalls,permcheck);
	//std::exit(1);
	
}


//////////// 1/6/13 Modifying code here to impelment fast reduced matrices build
//Comment: Will have to modify code below to accept changes
std::vector<double> Ar(npc*npc), Br(npc*npc), dAr(0,0);

permc.clear();
permc.resize(nz,0);
//Pv = P.get_vec(0,Nn,0,npc);
for (int ii = 0; ii < nz; ii++)
{	
	permc[ii] = 1;
	redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
	spess,area,arear,bi,ci,npc,Nn,permc);
	//dAr contains quick reduced matrices
	dAr.insert(dAr.end(),Ar.begin(),Ar.end());
	permc[ii] = 0;
}
permc.clear();

////////////1/6/13 End major modifications

///////3/31/14 Edit here to activate startk

std::vector<int> sequence0(ncomb,0), sequence;
for (int ii = 0; ii < ncomb; ii++){sequence0[ii] = ii;}
//basically cutting the deck so startk is the first k run then the rest run in sequence
sequence.insert(sequence.end(),sequence0.begin()+startk,sequence0.end());
sequence.insert(sequence.end(),sequence0.begin(),sequence0.begin()+startk);

///////3/31/14 End modification



for (int gacurrentcall = 0; gacurrentcall < gatotalcalls; gacurrentcall++)
{
	if (restart != 1)
	{
		//permcheck.equals(perm);
		ncalls = 0;
		Jdv.clear();
		rowells.clear();
	}
	if (readnowell <= 0) {nowell++;}
	else {nowell = readnowell;}
	
	///////////////Following is GA code///////////////////
	//GA parameters
	//Number required to reach convergence, default is 20
	//float nconv = 100;
	//Percentage of convergence, default is 0.99
	float pconv;
	if (OptCI == 'G' || OptCI == 'I')
	{
		pconv = 0.01;
	}
	else
	{
		pconv = 1;	
	}
	GAParameterList params;
	GASteadyStateGA::registerDefaultParameters(params);
	//params.set(gaNnGenerations, 200);
	//Default popsize is 30
	//params.set(gaNpopulationSize, 50);
	params.set(gaNscoreFrequency, 20);
	params.set(gaNflushFrequency, 20);
	//Default pmutation is 0.01 (1%)	(test to see if default does any different)
	params.set(gaNpMutation, 0.1);	
	params.set(gaNselectScores, (int)GAStatistics::AllScores);
	
	//Set up allele set array
	GARealAlleleSetArray alleles;	//Using an array to size the genome
	if (pickgenome == 1)	//set up first genome (new type)
	{
		//gene lower bound is -1 (gene: off)
		//std::vector<double> lb(nobszones,-1), ub(nobszones,0.0);
		lb1.resize(nobszones,-1);
		ub1.resize(nobszones,0.0);
		//Set gene upper bound to number of nodes in associated zone (-1 to account for C++ counting)
		for (int ii = 0; ii < nobszones; ii++) {ub1[ii] = ozonec[ii] - 1;}
		////Add nobzones gene alleles (# genes = nobszones)
		for (int ii = 0; ii < nobszones; ii++) {alleles.add(lb1[ii],ub1[ii],1,GAAllele::INCLUSIVE,GAAllele::INCLUSIVE);}
	}
	else if (pickgenome == 2)	//set up old type of genome
	{
		lb2 = 0;
		//Upper bound is total number of potential obs nodes
		ub2 = obswells.size() - 1;
		//Add nowell gene alleles (#genes = nowell)
		for (int ii = 0; ii < nowell; ii++) {alleles.add(lb2,ub2,1,GAAllele::INCLUSIVE,GAAllele::INCLUSIVE);}
	}
	
	//Create initial genome
	GARealGenome genome(alleles,Objective);
	//Run Simple GA
	GASimpleGA ga(genome);
	ga.parameters(params);
	//I shouldn't get negative scores but if I do I can switch to Sigma Truncation
	//GANoScaling scaling;
	GASigmaTruncationScaling scaling;
	ga.scaling(scaling);
	//terminate on convergance
	ga.pConvergence(pconv);
	ga.nConvergence(nconv);
	//ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
	//set output file name
	ga.set(gaNscoreFilename, "expdes.dat");
	//maximize OF
	std::cout << "OptCI = " << OptCI << std::endl;
	
	if (OptCI == 'G' || OptCI == 'I')
	{
		ga.minimize();
		std::cout << "minimize" << std::endl;
	}
	else
	{
		ga.maximize();
		std::cout << "maximize" << std::endl;
	}
	
	//ga.maximize();
	
	//mat2d permcheckleft.equals(permcheck);
	//////////////////////Iterations here to search parameter combinations
	//std::cout << begincomb << std::endl;
	/*
	if (runlonger == true) {begincomb = 0;}
	*/
	
	//2/2/15 make some edits here to be able to split up burden between multiple instances
	
	
	
	for (int kk = begincomb; kk < ncomb; kk++)
	{
		Arbuild.clear();
		Brbuild.clear();
		
		//permc = perm.get_vec(0, nz, kk, kk+1); //modified to activate startk
		permc = perm.get_vec(0, nz, sequence[kk], sequence[kk]+1);
		
		
		//std::cout << permc[0] << '\t' << permc[1] << std::endl;
		std::vector<double> permstore = permc;
		//permcheck.erase_rc(0,1);
		std::cout << "Building Reduced Matrices" << std::endl;
		//build baseline reduced matrices
//////////// 1/6/13 Modify code here to accept quick reduced matrices
		
		//Ar.clear();
		//Ar.resize(npc*npc);
		std::fill(Ar.begin(),Ar.end(),0.0);
		
		
		for (int ii = 0; ii < nz; ii++) 
		{
			for (int jj = 0; jj < npc*npc; jj++)
			{
				Ar[jj] += dAr[jj + ii*npc*npc]*permc[ii];
			}
		}
		
		//std::vector<double> Arq;
		//Arq = Ar;
		
		//Br.clear();
		//std::vector<double> Ar,Br;
		
		
		/*
		std::clock_t timeredstart;
		double duration;
		timeredstart = std::clock();
		
		Ar.clear();
		redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
		spess,area,arear,bi,ci,npc,Nn,permc);
		duration = ( std::clock() - timeredstart ) / (double) CLOCKS_PER_SEC;
		std::cout << "Time to build reduced matrices: " << duration << std::endl;
		std::exit(1);
		*/
		/*
		std::cout << Arq.size() << '\t' << Ar.size() << std::endl;
		double RMSE;
	
		for (int ii = 0; ii < npc*npc; ii++)
		{
			RMSE += (Ar[ii] - Arq[ii])*(Ar[ii] - Arq[ii]);
			
		}
		
		RMSE = sqrt(RMSE)/(npc*npc);
		std::cout << RMSE << std::endl;
		
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
					std::cout << jj + ii*npc << std::endl;
				}
				file << '\n';
			}
		}
		file.close();
		
		*/

//////////// 1/6/13 stop modify code here
		Arbuild.insert(Arbuild.end(),Ar.begin(),Ar.end());
		Brbuild.insert(Brbuild.end(),Br.begin(),Br.end());
		
		//build perturbed reduced matrices
		for (int ii = 0; ii < nz; ii++)
		{
			//permc[ii] += dH;
			//write code here for dividing by a percentage (1%) of permc
			permc[ii] *= 1.01;
			
			
			
			///////// 1/6/13 Modify code here to accept quick reduced matrices
			
			//Ar.clear();
			//Ar.resize(npc*npc);
			std::fill(Ar.begin(),Ar.end(),0.0);
			
			for (int kk = 0; kk < nz; kk++) 
			{
				for (int jj = 0; jj < npc*npc; jj++)
				{
					Ar[jj] += dAr[jj + kk*npc*npc]*permc[kk];
				}
			}
			
			/*
			redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
			spess,area,arear,bi,ci,npc,Nn,permc);
			*/
			//////////// 1/6/13 stop modify code here
			Arbuild.insert(Arbuild.end(),Ar.begin(),Ar.end());
			Brbuild.insert(Brbuild.end(),Br.begin(),Br.end());
			permc[ii] = permstore[ii];
		}
		//permpicked.augment_vec(permc);	
		
		std::cout << "Running GA" << std::endl;
		//Run GA
		unsigned int seed = 1E4;
		ga.initialize(seed);
		int count;
		
		//////////something here for runlonger
		/*
		if (runlonger == true)
		{
			std::cout << "Run Longer" << std::endl;
			for (int ii = 0; ii < genome.length(); ii++)
			{
				genome.gene(ii,rowells[corresro[kk]].genome[ii]);
				std::cout << genome.gene(ii) << std::endl;
			}
			ncalls = 0;
		}
		*/
		
		
		while (!ga.done())
		{
		count++;
		ga.step();
		if (OptCI == 'G' || OptCI == 'I')
		{
			std::cout << std::resetiosflags(std::ios_base::scientific) << 
			'\r' << "Current model call number: " << ncalls << " best score: " << 
			ga.statistics().minEver() << " Current Convergence: " << ga.convergence() << "                         " << std::endl;
			ga.flushScores();
		}
		else
		{
			std::cout << std::resetiosflags(std::ios_base::scientific) << 
			'\r' << "Current model call number: " << ncalls << " best score: " << 
			ga.statistics().maxEver() << " Current Convergence: " << ga.convergence() << "                         " << std::endl;
			ga.flushScores();
		}
		std::cout << std::endl;
		
		
		/*
		convergence = ga.convergence();
		if (convergence == 0 && ncalls > 1000)
		{
			for (unsigned int yy = 0; yy < diagh.size(); yy++)
			{
				std::cout << diagh[yy] << '\t';
			}
			std::cout << std::endl;
			genome = ga.statistics().bestIndividual();
			
				//Storing unique solutions, cataloging scores of all occurances of non-unique solutions
			std::vector<owell> copt;	//current optimal solution
			currentzero = 0;
			std::cout << pickgenome << '\t' << genome.length() << std::endl;
		for (int ii = 0; ii < genome.length(); ii++)
		{
			std::cout << genome.gene(ii) << '\t';
			if (pickgenome == 1 && genome.gene(ii) != -1)	//Store solution from genome type 1
			{
				copt.push_back(obswells[genome.gene(ii) + currentzero]);
				std::cout  << genome.gene(ii) << '\t' << copt[copt.size()-1].node << '\t' << 
				copt[copt.size()-1].zone << std::endl;
			}
			else if (pickgenome == 2) {copt.push_back(obswells[genome.gene(ii)]);}	//Store solution from genome type 2
			currentzero += ozonec[ii];
		}
		*/
		/*
		//std::cout << "wellsame" << std::endl;
		if (OptCI == 'G' || OptCI == 'I')
		{
			wellsame(rowells,copt,close,ga.statistics().minEver(),genome, kk);
		}
		else
		{
			wellsame(rowells,copt,close,ga.statistics().maxEver(),genome, kk);	//Check to see if the current solution is unique
		}
		int kleft = ncomb - (kk+1);
		write_restart(nowell,rowells,ncalls,kleft);
		
		//std::exit(1);
		*/
		}
		
		
		//std::cout << std::endl;
		genome = ga.statistics().bestIndividual();
		//Storing unique solutions, cataloging scores of all occurances of non-unique solutions
		std::vector<owell> copt;	//current optimal solution
		currentzero = 0;
		for (int ii = 0; ii < genome.length(); ii++)
		{
			if (pickgenome == 1 && genome.gene(ii) != -1)	//Store solution from genome type 1
			{
				copt.push_back(obswells[genome.gene(ii) + currentzero]);
				std::cout  << genome.gene(ii) << '\t' << copt[copt.size()-1].node << '\t' << 
				copt[copt.size()-1].zone << std::endl;
			}
			else if (pickgenome == 2) {copt.push_back(obswells[genome.gene(ii)]);}	//Store solution from genome type 2
			currentzero += ozonec[ii];
		}
		//std::cout << "wellsame" << std::endl;
		if (OptCI == 'G' || OptCI == 'I')
		{
			wellsame(rowells,copt,close,ga.statistics().minEver(),genome, kk);
		}
		else
		{
			wellsame(rowells,copt,close,ga.statistics().maxEver(),genome, kk);	//Check to see if the current solution is unique
		}
		
		//if (kk > 0 && kk % restartfreq == 0)
		if (kk % restartfreq == 0)
		{
			std::cout << "Saving Restart Point" << std::endl;
			int kleft = ncomb - (kk+1);
			write_restart(nowell,rowells,ncalls,kleft);
			//std::exit(1);
		}
	
	}
	
	
	//Save final solution before running worst case scenario
	write_restart(nowell,rowells,ncalls,0);
	//std::cout << "crash allowell" << std::endl;
	std::vector<owell> allowell(1);
	//count number of occurances of an observation well
	countowells(rowells,allowell,close);
	std::vector<owell> copt;
	//store top nowell as the highest frequency occurances
	copt.insert(copt.begin(),allowell.begin(),allowell.begin()+nowell);
	
	//for (unsigned int ii = 0; ii < copt.size(); ii++) {std::cout << copt[ii].gene << '\t';}
	//std::cout << std::endl;
	//assign the genome
	for (int ii = 0; ii < genome.length(); ii++)
	{
		//genome.gene(0,-1);
		genome.gene(ii,-1);
	}
	for (unsigned int ii = 0; ii < copt.size(); ii++)
	{
		genome.gene(copt[ii].zonecount,copt[ii].gene);
		copt[ii].highfreq = true;
	}
	
	//for (int ii = 0; ii < genome.length(); ii++){std::cout << genome.gene(ii) << '\t';}
	
	//std::cout << '\n' << "check score" << std::endl;
	//z = Objective(genome);
	//z = 10E10;
	
	//std::exit(1);
	
	
	//check to see if this one has occured before
	//std::cout << "crash see if occured before" << std::endl;
	wellsame(rowells,copt,close,0.0,genome,0);
	
	
	//std::cout << "crash sort by network count" << std::endl;
	std::sort(rowells.begin(),rowells.end(),by_count_network());	//sort rowells by frequency of network occurance
	rowells[0].highfreqnetwork = true;
	//std::cout << "crash move owell count to front" << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)			//move high frequency wells to front of list
	{
		if (rowells[ii].highfreqwells == true)
		{
			robustowell temp = rowells[ii];
			rowells.erase(rowells.begin() + ii, rowells.begin() + ii + 1);
			rowells.insert(rowells.begin(),temp);
			break;
		}
	}
	
	

	//std::cout << "crash worst case" << std::endl;
	//Running Worst Case Scenario for all combinations in rowells
	//for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		//For each perm combination
		std::cout << "Running Worst Case Scenario Test" << std::endl;
		for (int jj = 0; jj < perm.columns(); jj++)
		{
			std::cout << '\r' << "Parameter Combintaion: " << jj << "                     ";
			//std::cout << jj << std::endl;
			permc = perm.get_vec(0, nz, jj, jj+1);
			std::vector<double> permstore = permc;
			Arbuild.clear();
			Brbuild.clear();
			//build baseline reduced matrices
			
			//Ar.clear();
			//Ar.resize(npc*npc);
			std::fill(Ar.begin(),Ar.end(),0.0);

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
					Ar[jj] += dAr[jj + ii*npc*npc]*permc[ii];
				}
			}
			
			////////// 1/6/14 stop modify code here
			
			Arbuild.insert(Arbuild.end(),Ar.begin(),Ar.end());
			Brbuild.insert(Brbuild.end(),Br.begin(),Br.end());
			//build perturbed reduced matrices
			for (int kk = 0; kk < nz; kk++)
			{
				//std::cout << kk << std::endl;
				permc[kk] += dH;
				////////// 1/6/14 modify code here for quick red matrix
				/*
				redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
				spess,area,arear,bi,ci,npc,Nn,permc);
				*/
				//Ar.clear();
				//Ar.resize(npc*npc);
				std::fill(Ar.begin(),Ar.end(),0.0);
				for (int ii = 0; ii < nz; ii++) 
				{
					for (int jj = 0; jj < npc*npc; jj++)
					{
						Ar[jj] += dAr[jj + ii*npc*npc]*permc[ii];
					}
				}
				
				////////// 1/6/14 stop modify code here
				Arbuild.insert(Arbuild.end(),Ar.begin(),Ar.end());
				Brbuild.insert(Brbuild.end(),Br.begin(),Br.end());
				permc[kk] = permstore[kk];
			}
			//check each robust well combination
			for (unsigned int ii = 0; ii < rowells.size(); ii++)
			{
				if (jj == 0) 
				{
					rowells[ii].maxscore.clear();
					rowells[ii].minscore = 10E100;
				}
				for (int kk = 0; kk < genome.length(); kk++)
				{
					genome.gene(kk,rowells[ii].genome[kk]);
				}
				//std::cout << "crash2" << std::endl;
				z = Objective(genome);
				rowells[ii].maxscore.push_back(z);	//conversion from float to double
				if (z <= rowells[ii].minscore)
				{
					rowells[ii].minscore = z;
					rowells[ii].worstcaseK = permc;
				}
			}
		}
		double maxworst = -10E100;
		int maxworstindex;
		for (unsigned int ii = 0; ii < rowells.size(); ii++)
		{
			if(rowells[ii].minscore > maxworst)
			{
				maxworst = rowells[ii].minscore;
				maxworstindex = ii;
			}
		}
		rowells[maxworstindex].robustdesign = true;
		robustowell temp = rowells[maxworstindex];
		rowells.erase(rowells.begin() + maxworstindex, rowells.begin() + maxworstindex + 1);
		rowells.insert(rowells.begin(),temp);	
	}
	std::cout << '\n' << std::endl;
	/*
	std::cout << "crash sort by network count" << std::endl;
	std::sort(rowells.begin(),rowells.end(),by_count_network());	//sort rowells by frequency of network occurance
	std::cout << "crash move owell count to front" << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)			//move high frequency wells to front of list
	{
		if (rowells[ii].highfreqwells == true)
		{
			robustowell temp = rowells[ii];
			rowells.erase(rowells.begin() + ii, rowells.begin() + ii + 1);
			rowells.insert(rowells.begin(),temp);
		}
	}
	*/
	//std::cout << "crash write final" << std::endl;
	write_final(nowell, rowells, ncalls, nz, allowell);
	/*
	//Print out final solutions
	std::ostringstream convert;
	convert << nowell;
	//std::string snowell = convert.str();
	std::string filename = ("Robust solutions for "+convert.str()+" observation wells.txt").c_str();
	
	std::ofstream file (filename);
	if (file.is_open())
	{
		for (unsigned int ii = 0; ii < rowells.size(); ii++)
		{
			//Write to file
			file << "Number of bucket " << ii << " encountered is " << rowells[ii].count << std::endl
			<< "Size of bucket is " << rowells[ii].optowells.size() << std::endl
			<< "Worst Case Scenario score is " << rowells[ii].minscore << std::endl
			<< "Worst Case Hydraulic Conductivity is" << std::endl;
			for (int jj = 0; jj < nz; jj++) {file << rowells[ii].worstcaseK[jj] << '\t';}
			file << std::endl;
			for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
			{
				file << "x = " << rowells[ii].optowells[jj].x << '\t' << "y = " << rowells[ii].optowells[jj].y << std::endl;
			}
			
			//Print to screen
			std::cout << "Number of bucket " << ii << " encountered is " << rowells[ii].count << std::endl
			<< "Size of bucket is " << rowells[ii].optowells.size() << std::endl
			<< "Worst Case Scenario score is " << rowells[ii].minscore << std::endl
			<< "Worst Case Hydraulic Conductivity is" << std::endl;
			for (int jj = 0; jj < nz; jj++) {std::cout << rowells[ii].worstcaseK[jj] << '\t';}
			std::cout << std::endl;
			for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
			{
				std::cout << "x = " << rowells[ii].optowells[jj].x << '\t' << "y = " << rowells[ii].optowells[jj].y << std::endl;
			}
		}
		//Print out total number of reduced model calls
		file << '\n' << "Total number of reduced model calls is " << ncalls << std::endl;
		std::cout << '\n' << "Total number of reduced model calls is " << ncalls << '\n' << std::endl;
		
	}
	else
	{
		std::cout << "Cannot open file to write solutions";
	}
	*/
	
}




/*
/////////////Figuring out code
Two genomes, one like old one, new one where each genome has assigned zone to search, new genome has two options, 
(1) throw out or (2) repair bad individuals
GARealAlleleSetArray alleles
objective...
*/
/*
int nozone = obszonelist.size();
std::vector< std::vector<int> > genestore(nobs, std::vector<int>(0,0));
*/

/*
Use GAArrayAlleleSet and GA1DArrayAlleleGenome
*/
/*
std::cout << std::endl;
int t = clock();
float time = t/CLOCKS_PER_SEC;
//std::string tunit = " seconds";

if (time >= 3600)
{
time /= 3600;
std::cout << "The GA took " << time << " hours to run" << std::endl;
}	else if (time >= 60)
{
time /= 60;
std::cout << "The GA took " << time << " minutes to run" << std::endl;
}	else
{
std::cout << "The GA took " << time << " seconds to run" << std::endl;
}
*/

}