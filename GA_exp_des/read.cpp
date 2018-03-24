#include "addon.h"
#include "matrix.h"
#include "fort.h"
#include <cstdlib>
#include "owellstruct.h"
#include <algorithm>

void read_fort(int &nt, int &irad, int &nterm, std::vector<int> &triang2, 
std::vector<int> &trija2, std::vector<double> &lmass2,
std::vector<double> &xc2, std::vector<double> &yc2, 
std::vector<double> &elstor2, std::vector<double> &spess2, 
std::vector<double> &area2, std::vector<double> &arear2, 
std::vector<double> &bi2, std::vector<double> &ci2, 
std::vector<int> &ja2, std::vector<int> &topol2)
{
int iterm;						//If you change any of the following values make sure to change SAT2D.H
int ntrmax = 57888;				//Max # triangles
int nmax = 29197;				//Max # nodes
int npmax = 504;				//Max # Dirichlet nodes
int nqmax = 25;					//Max # Neumann nodes
int nrmax = 25;					//Max # of partial output nodes
int maxprt = 200;				//Max # of detailed output
//int n1max = 20;					//Max # of element connections to a node
//int ntpmax = n1max*nmax;
int maxzon = 15;				//Max # of material zones
//int maxtrm = n1max*nmax;

triang2.resize(ntrmax*4);

int *triang = &triang2[0];
std::vector<double> x2(nmax);
double *x = &x2[0];
std::vector<double> y2(nmax);
double *y = &y2[0];
std::vector<double> permx2(maxzon);
double *permx = &permx2[0];
std::vector<double> permy2(maxzon);
double *permy = &permy2[0];
elstor2.resize(maxzon);
double *elstor = &elstor2[0];
spess2.resize(maxzon);
double *spess = &spess2[0];
std::vector<int> contp2(npmax);
int *contp = &contp2[0];
std::vector<int> contq2(nqmax);
int *contq = &contq2[0];
std::vector<int> contr2(nrmax);
int *contr = &contr2[0];
std::vector<double> timprt2(maxprt), ptimep2(ntrmax);
double *timprt = &timprt2[0];
double *ptimep = &ptimep2[0];

double tetaf, deltat, dtmax, tmax, dtmaga, dtmagm, tolcg;
int itmxcg, lump, iprt1, iprt, indp, nq, nzone, n1, nr, nprt, n, np;
//std::cout << "openio" << std::endl;
openio_(iterm);
//std::cout << "datin" << std::endl;
datin_(iterm, triang, x, y, permx, permy, elstor, spess, contp, contq,
contr, timprt, ptimep, tetaf, deltat, dtmax, tmax, dtmaga, dtmagm, itmxcg, tolcg,
lump, irad, iprt1, iprt, indp, nq, nzone, n1, nr, nprt, n, np, nt);

std::vector<int> ip32(9);
ip32[0] = 1;
ip32[1] = 2;
ip32[2] = 3;
ip32[3] = 2;
ip32[4] = 3;
ip32[5] = 1;
ip32[6] = 3;
ip32[7] = 1;
ip32[8] = 2;
int *ip3 = &ip32[0];
bi2.resize(nt*3);
ci2.resize(nt*3);
double *bi = &bi2[0];
double *ci = &ci2[0];
std::vector<double> areanod2(n);
double *areanod = &areanod2[0];
area2.resize(nt);
arear2.resize(nt);
double *area = &area2[0];
double *arear = &arear2[0];
std::vector<int> iarea2(nt);
int *iarea = &iarea2[0];
//std::cout << "areabas" << std::endl;
arebas_(n, nt, triang, ip3, bi, ci, areanod, area, arear, iarea, x, y, iterm);

ja2.resize(n1*n);
int *ja = &ja2[0];
topol2.resize(n+1);
int *topol = &topol2[0];
trija2.resize(3*3*ntrmax);
int *trija = &trija2[0];
int imax = 2147483647, ndz;
std::vector<int> ier2(7);
int *ier = &ier2[0];
strpic_(n, nterm, triang, ja, topol, nt, n1, imax, iterm);
chkpic_(n, nterm, topol, ja, ndz, ier);
tripic_(nt, triang, ja, topol, trija);

lmass2.resize(3*3);

double *lmass = &lmass2[0];
//std::cout << "locmas" << std::endl;
locmas_(lmass, lump);
xc2.resize(ntrmax);
yc2.resize(ntrmax);
double *xc = &xc2[0];
double *yc = &yc2[0];
//std::cout << "nodelt x" << std::endl;
nodelt_(nt, triang, x, xc);
//std::cout << "nodelt y" << std::endl;
nodelt_(nt, triang, y, yc);
//std::cout << "resizing" << std::endl;
elstor2.resize(nzone);
spess2.resize(nzone);
trija2.resize(3*3*nt);
triang2.resize(4*nt);
xc2.resize(n);
yc2.resize(n);
//std::cout << "end of read_fort" << std::endl;
}

void read0(int& nowell, int &npl, int &startk, std::vector<double> &klim, std::vector<double> &obstimes, 
std::vector<int> &obszonelist, std::vector<owell> &obswells, 
std::vector<int> &ozonec, double &dt, double &dtmult, int &lumpwell, double &dH, double &dQ,
int &pickgenome, int &genomeoption, int &nobszones, int &nobstime, int &expdesparam, int &expdespump,
char &OptCI, char &OptCG, double &close, int &restart, int &restartfreq, double &dog, float &nconv,
bool &runlonger)
//std::vector<int> &obswellloc, std::vector<double> &xy, std::vector<int> &obszones)
{
std::ifstream file("GA_search_params.txt");
if (file.is_open())
{
	file.seekg(0);
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file >> nowell;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file >> npl;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	klim.resize(2);
	file >> startk >> klim[0] >> klim[1];
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	int runlongertemp;
	file >> dt >> dtmult >> dH >> dQ >> close >> restart >> restartfreq >> dog >> nconv >> runlongertemp;
	if (runlongertemp == 1) {runlonger = true;}
	else {runlonger = false;}
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file >> lumpwell >> pickgenome >> genomeoption >> expdesparam >> expdespump >> OptCI >> OptCG;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	//int nobstime;
	file >> nobstime;
	obstimes.resize(nobstime);
	for (int ii = 0; ii < nobstime; ii++) {file >> obstimes[ii];}
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	//int nobszones;
	file >> nobszones;
	obszonelist.resize(nobszones);
	for (int ii = 0; ii < nobszones;ii++) {file >> obszonelist[ii];}
	std::sort(obszonelist.begin(),obszonelist.end());
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	int nobsloc;
	file >> nobsloc;
	//std::cout << nobsloc << std::endl;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	//load up vector of custom structures (working here on using structures to hold obswell data)
	obswells.resize(nobsloc);
	std::vector<int> countzones(nobsloc,0);
	int count, zone, noden;
	double x, y;
	for (int ii = 0; ii < nobsloc; ii++)
	{
		file >> count >> noden >> zone >> x >> y;
		obswells[ii].node = noden;
		obswells[ii].zone = zone;
		countzones[ii] = zone;
		obswells[ii].x = x;
		obswells[ii].y = y;
	}
	//std::cout << obswells[0].x << std::endl;
	//std::cout << "check read0 " << obswells[obswells.size()-1].y << std::endl;
	//sort by zone
	std::sort(obswells.begin(), obswells.end(), by_zone());
	
	//count number of nodes in each obs zone
	ozonec.resize(nobszones);
	for (int ii = 0; ii < nobszones; ii++)
	{
		ozonec[ii] = std::count(countzones.begin(), countzones.end(), obszonelist[ii]);
	}
	/*
	std::cout << obswells[0].node << '\t' << obswells[0].zone << std::endl;
	std::cout << obswells[1].node << '\t' << obswells[1].zone << std::endl;
	std::sort(obswells.begin(), obswells.begin()+2, by_zone());
	std::cout << obswells[0].node << '\t' << obswells[0].zone << std::endl;
	std::cout << obswells[1].node << '\t' << obswells[1].zone << std::endl;
	*/
	/*
	obswellloc.resize(nobsloc);
	xy.resize(2*nobsloc);
	obszones.resize(nobsloc);
	int count;
	for (int ii = 0; ii < nobsloc; ii++)
	{
		file >> count >> obswellloc[ii] >> obszones[ii] >> xy[ii] >> xy[ii+1];
	}
	*/
	file.close();
}
else
{
std::cout<<"Cannot open GA_search_params.txt to read" << std::endl;
std::exit(1);
}
}


void read(int &npc, int &Nn, int &nq, std::vector<int> &well, 
std::vector<double> &wellweight, mat2d &P, int &nz)
{
std::ifstream file("red_model_param.txt");
if (file.is_open())
{
	file.seekg(0);
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	//std::cout << "various" << std::endl;
	file >> npc >> Nn >> nq >> nz;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	well.resize(nq);
	//std::cout << "wells" << std::endl;
	for (int ii = 0; ii < nq; ii++)
	{
		file >> well[ii];
	}
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	wellweight.resize(2*nq);
	//std::cout << "well weight" << std::endl;
	for (int ii = 0; ii < nq*2; ii++)
	{
		file >> wellweight[ii];
	}
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	std::vector<double> Ptemp(Nn*npc);
	std::cout << "Reading in P matrix" << std::endl;
	for (int ii = 0; ii < Nn*npc; ii++)
	{
		file >> Ptemp[ii];
		//std::cout << ii << '\t';
	}
	//std::cout << "check read " << Ptemp[3] << '\t' << Ptemp[Ptemp.size() - 20];
	P.vec_to_mat(Ptemp,Nn,npc,1);
	//P.print();
	//std::cin.get();
	file.close();
	
	//for (unsigned int ii = 0; ii < wellweight.size(); ii++) {std::cout << wellweight[ii] << '\t';}
	//std::exit(1);
	
}
else {
std::cout << "Could not open red_model_param.txt to read" << std::endl;
std::exit(1);
}
}

void read_restart(int &nowell, int &ncalls, std::vector<robustowell> &rowells,
int nz, int &columns)
{
std::ifstream file("Restart_Exp_Des.txt");
//std::cout << "Trying to open restart file" << std::endl;
if (file.is_open())
{
	//std::cout << "Opened restart file" << std::endl;
	int rowellssize, genomesize;
	file.seekg(0);
	file >> rowellssize >> nowell >> columns >> genomesize >> ncalls;
	rowells.resize(rowellssize);
	//std::cout << "crash here1" << std::endl;
	for (int ii = 0; ii < rowellssize; ii++)
	{
		//double x, y;
		file >> rowells[ii].count;
		rowells[ii].paramindex.resize(rowells[ii].count);
		for (int jj = 0; jj < rowells[ii].count; jj++)
		{
			file >> rowells[ii].paramindex[jj];
		}
		//std::cout << "ii = " << ii << '\t' << rowells[ii].paramindex[0] << std::endl;
		/*
		rowells[ii].optowells.resize(nowell);
		for (int jj = 0; jj < nowell; jj++)
		{
			file >> x >> y;
			rowells[ii].optowells[jj].x = x;
			rowells[ii].optowells[jj].y = y;
		}
		*/
		//rowells[ii].genome.resize(genomesize);
		int checkgenomesize;
		file >> checkgenomesize;
		rowells[ii].genome.resize(checkgenomesize);
		for (int jj = 0; jj < checkgenomesize; jj++)
		{
			file >> rowells[ii].genome[jj];
		}
		rowells[ii].minscore = 10E100;
		//std::cout << rowells[ii].genome[genomesize-1] << std::endl;
	}
	//file.ignore(500, '\n');
	//file.ignore(500, '\n');
	/*
	std::vector<double> permcheckv(nz*columns);
	//double permchecki;
	for (int ii = 0; ii < nz*columns; ii++)
	{
		file >> permcheckv[ii];
	}
	permcheck.vec_to_mat(permcheckv,nz,columns,1);
	*/
	//std::cout << "close file" << std::endl;
	file.close();
}
else
{
	std::cout << "Could not open restart file to read" << std::endl;
}
}

void read_run_param(std::vector<double> &ts, double &m, std::vector<double> &permc, int &nz)
{
std::cout << "read_run_param" << std::endl;
std::ifstream file("Model Run params.txt");
if (file.is_open())
{
	file.seekg(0);
	file.ignore(500, '\n');
	int otimes;
	file >> otimes;
	std::cout << "otimes: " << otimes << std::endl;
	ts.resize(otimes);
	for (int ii = 0; ii < otimes; ii++)
	{
		file >> ts[ii];
		std::cout << ts[ii] << std::endl;
	}
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	file >> m;
	file.ignore(500, '\n');
	file.ignore(500, '\n');
	permc.resize(nz);
	for (int ii = 0; ii < nz; ii++)
	{
		file >> permc[ii];
	}
	file.close();
}
else {std::cout << "Could not open model run parameter file to read" << std::endl;}
}

void read_interesting_points(std::vector<int> &poi, std::vector<double> &toi, int &npoi, int &ntoi)
{
std::ifstream file("Points of interest.txt");
if (file.is_open())
{
	file.seekg(0);
	file >> npoi;
	poi.clear();
	poi.resize(npoi);
	file.ignore(500, '\n');
	for (int ii = 0; ii < npoi; ii++)
	{
		file >> poi[ii];
		//std::cout << poi[ii];
	}
	file >> ntoi;
	toi.clear();
	toi.resize(ntoi);
	file.ignore(500, '\n');
	for (int ii = 0; ii < ntoi; ii++)
	{
		file >> toi[ii];
	}
	file.close();
}
else {std::cout << "Could not open interesting points file to read" << std::endl;}
}