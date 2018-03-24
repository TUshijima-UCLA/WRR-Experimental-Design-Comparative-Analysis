//Standard headers
#include <set>
#include <cstdlib>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <ctime>
#include <cmath>
//Personal headers
#include "owellstruct.h"
#include "matrix.h"
#include "float.h"
//GA headers
#include <ga/ga.h>
#include <ga/std_stream.h>
#include <ga/GARealGenome.h>
//Global variables
extern mat2d P;
extern std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, permc, Pv, wellrate, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;

extern std::vector<double> lb1, ub1, dI;
extern double multconst;

extern std::vector<double> diagh;
extern std::vector<int> poi;
extern int ntoi,npoi;

extern std::vector<int> topol, JA, triang, trija, well, ozonec;
extern std::vector<owell> obswells;			//all potential obsrvation nodes
extern int ntri, irad, nterm, npc, Nn, nprt, lumpwell, pickgenome, genomeoption, ncalls, 
nowell, nobszones, nq, nz, expdesparam, expdespump;	//nz: total number of hydrologic zones in the aquifer
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

//Recursive function to randomly repair a bad #1 genome (genomeoption == 1)
//genome, current zero, #owell allowed, #owell currently on
void fixgenome(GARealGenome &g, int cz, int noa, int no)
{
//Base case 1: If the number of allowable obs wells equals the number of active obs well
if (noa == no) {return;}
//Base case 2: If there are no more allowable obs wells, turn all remaining genes off
if (noa == 0)
{
	for (cz = cz; cz < g.length(); cz++) {g.gene(cz,-1);}
	no = 0;
	return;
}
//fixed seed for fix genome function
unsigned seed = 5;
//"random" seed for fix genome function
//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

//set up random number generator
std::default_random_engine generator (seed);
//use a uniform distribution on [0,noa];
std::uniform_int_distribution<int> distribution(0,noa);
//random owell on to turn off
int E = distribution(generator);

int ecount = 0;		//counter till turn off
while (ecount < E)	//count number of active obs wells passed
{
	if (g.gene(cz) != -1) {ecount++;}
	cz++;
}
while (g.gene(cz) == -1) {cz++;}	//move to next active obs well
g.gene(cz, -1);	//turn the Eth active obs well off
noa -= E;		//How many allowable obs wells in remaining genes
no -= (E+1);	//How many active obs wells left in remaining genes
cz++;			//current gene being considered
fixgenome(g,cz,noa,no);
}

void fixgenometoofew(GARealGenome &g, int cz, int noa, int no, int np)
{
//Base case 1: If I don't need to turn any on (noa) then exit
if (noa == 0) {
//std::cout << "Base case 1" << std::endl;
return;
}

//fixed seed for fix genome function
unsigned seed = 1E4;
//unsigned seed = 5;
//"random" seed for fix genome function
//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
//set up random number generator
std::default_random_engine generator (seed);

//Base case 2: If the number to turn on equals the number remaining, turn all remaining on
if (noa == (g.length() - cz))
{
	for (cz = cz; cz < g.length(); cz++) 
	{
		//use a uniform distribution on [0,noa];
		std::uniform_int_distribution<int> distribution(0,ub1[cz]);
		int E = distribution(generator);
		g.gene(cz,E);
	}
	no = noa;
	//std::cout << "Base case 2" << std::endl;
	return;
}
//If neither base case, pick random gene to turn on
int ubstep = g.length() - cz - noa - (no-np);
std::uniform_int_distribution<int> distribution(0,ubstep);
//random owell on to turn on
int E = distribution(generator);
//std::cout << "ubstep " << ubstep << "\tstep " << E << "\tcz " << cz << "\tnoa " << noa << "\tno " << no << "\tnp " << np << std::endl;
//int ecount = 0;		//counter till turn on
/*
while (ecount < E)	//count number of active obs wells passed
{
	if (g.gene(cz) == -1) {ecount++;}
	else {np++;}
	cz++;
}
*/
for (int ii = cz; ii < cz+E; ii++)
{
	if (g.gene(ii) != -1) {np++;}
}
cz += E;
while (g.gene(cz) != -1) {cz++;np++;}	//move to next inactive obs well
std::uniform_int_distribution<int> distribution2(0,ub1[cz]);
int E2 = distribution2(generator);
g.gene(cz,E2);	//turn the Eth inactive obs well on
noa--;
no++;
np++;
cz++;
/*
for (unsigned int yy = 0; yy < ub1.size(); yy++)
{
	//std::cout << ub1[yy] << '\t';
	std::cout << g.gene(yy) << '\t';
}
std::cout << std::endl;
std::cout << "E2 " << E2 << "\tcz " << cz << "\tnoa " << noa << "\tno " << no <<std::endl;
*/

fixgenometoofew(g,cz,noa,no,np);

}



float Objective(GAGenome &g)
{
//float z = 0.0;	//score value
z = 0.0;
GARealGenome &genome = (GARealGenome&)g;
//Get genome into observation wells, into Projection matrix
std::vector<owell> cowell;
mat2d Ps; //"small" projection matrix (only contain observation well locations)
std::vector<double> Psv(0,0);	//vector holding Ps in row form
//2 genomes, one an array w/ each element taking on one of the allowable zones, two ways of handeling, throw out bad or repair bad
//the other like the one is the last GA
//All options pass out the same thing, some small projection matrix
//First option is hard w/ 1 owell, easier for almost all owells; second option is easy w/ 1 owell, hard w/ almost all owells
int countwell;

//	std::cout << "start of obj fn" << std::endl;

if (pickgenome == 1)	//pick the first genome
{
	/*
	std::cout << "Setting genome" << std::endl;
	////////test code by setting gene
	
	genome.gene(0,-1);
	genome.gene(1,-1);
	genome.gene(2,-1);
	genome.gene(3,25);
	*/
	/*
	genome.gene(0,-1);
	genome.gene(1,0);
	genome.gene(2,5);
	genome.gene(3,25);
	*/
	////////////////////////////////
	
	/*
	genome.gene(0, -1);
	genome.gene(1, -1);
	genome.gene(2, -1);
	genome.gene(3,  4);
	genome.gene(4, -1);
	genome.gene(5, 34);
	genome.gene(6,  1);
	genome.gene(7, 10);
	genome.gene(8, -1);
	genome.gene(9, 16);
	genome.gene(10,39);
	genome.gene(11, 4);
	*/	
	/*
	genome.gene(0, -1);
	genome.gene(1, -1);
	genome.gene(2, 0);
	genome.gene(3, 3);
	genome.gene(4, -1);
	genome.gene(5, -1);
	genome.gene(6,  -1);
	genome.gene(7,  18);
	genome.gene(8, -1);
	genome.gene(9, 7);
	genome.gene(10,38);
	genome.gene(11, 0);
	*/
	/*
	genome.gene(0, -1);
	genome.gene(1, -1);
	genome.gene(2, 13);
	genome.gene(3, -1);
	genome.gene(4, 5);
	genome.gene(5, -1);
	*/
	/*
	countwell = nobszones;
	for (int ii = 0; ii < genome.length(); ii++) 
	{
		if (genome.gene(ii) == -1) {countwell--;}
//		std::cout << genome.gene(ii) << '\t';
	}
	if (countwell < nowell)
	{
		std::cout << countwell << std::endl;
		fixgenometoofew(genome,0,nowell - countwell,countwell,0);
		countwell = nowell;	
	}
	
	for (unsigned int yy = 0; yy < ub1.size(); yy++)
	{
		std::cout << ub1[yy] << '\t';
		//std::cout << genome.gene(yy) << '\t';
	}
	std::cout << std::endl;

	std::exit(1);
	*/
	
	
	countwell = nobszones;
	for (int ii = 0; ii < genome.length(); ii++) 
	{
		if (genome.gene(ii) == -1) {countwell--;}
//		std::cout << genome.gene(ii) << '\t';
	}
		
//	std::cout << "countwell = " << countwell << std::endl;
	if (countwell > nowell)
	{
//		std::cout << "genomeoption = " << genomeoption << std::endl;
		//throw out bad
		if (genomeoption == 0)
		{
			z = 0.0;
			return z;
		}
		//repair the bad randomly
		else if (genomeoption == 1)	
		{
			fixgenome(genome,0,nowell,countwell);
			countwell = nowell;
		}
	}
	
	if ((OptCI == 'G' || OptCI == 'I') && (countwell < nowell))
	{
		if (genomeoption == 0)
		{
			//z = zmin;
			z = 0.0;
			return z;
		}
		else if (genomeoption == 1)	
		{
			fixgenometoofew(genome,0,nowell - countwell,countwell,0);
			countwell = nowell;
		}
	}
	
	
//	P.print();
	//create Ps from this genome type
	int currentzero = 0;
	for (int ii = 0; ii < genome.length(); ii++)
	{
//		std::cout << genome.gene(ii) << '\t';
		if (genome.gene(ii) != -1)
		{
			std::vector<double> pvi = P.get_vec(obswells[genome.gene(ii) + currentzero].node-1,
			obswells[genome.gene(ii) + currentzero].node,0,npc);
//			std::cout << '\t' << obswells[genome.gene(ii) + currentzero].node << '\t' << obswells[genome.gene(ii) + currentzero].zone;
//			std::cout << obswells[genome.gene(ii) + currentzero].node << '\t';
			Psv.insert(Psv.end(),pvi.begin(),pvi.end());
		}
//		std::cout << std::endl;
		currentzero += ozonec[ii];
//		std::cout << "currentzero = " << currentzero << std::endl;
	}
//	std::cout << std::endl;
}
else if (pickgenome == 2)	//pick second genome, currently no method to repair bad genes
{
	std::vector<int> owellzones;
	std::vector<owell> cowells;	//obs wells currently being considered
	//std::cout << "check if repeated zones" << std::endl;
	for (int ii = 0; ii < genome.length(); ii++)
	{
		owellzones.push_back(obswells[genome.gene(ii)].zone);
		//std::cout << "gene " << genome.gene(ii) << "\tzone " << owellzones[ii] << std::endl;
	}
	//std::cout << std::endl;
	std::set<int> check(owellzones.begin(),owellzones.end());	//a set only holds one copy of a value
	if (check.size() < owellzones.size())	//if there are repeated zones 
	{
		/*
		//std::cout << "reject for mutliple zones" << std::endl;
		if (OptCI == 'G' || OptCI == 'I')
		{
			//z = zmin;
			z = 0.0;
		}
		else
		{
			z = 0.0;
		}
		return z;	//test by taking out return call here
		*/
		//changing code here because all have min of zero with new formulation
		z = 0.0;
		return z;
	}
	//create Ps from this genome type
	for (int ii = 0; ii < genome.size(); ii++)
	{
		std::vector<double> pvi = P.get_vec(obswells[genome.gene(ii)].node-1, obswells[genome.gene(ii)].node,0,npc);
//		std::cout << genome.gene(ii) << '\t' << obswells[genome.gene(ii)].node <<'\t' << obswells[genome.gene(ii)].zone << '\n';
		Psv.insert(Psv.end(),pvi.begin(),pvi.end());
		countwell = nowell;
	}
}
else
{
std::cout << "Unknown genome" << std::endl;
std::exit(1);
}
/*
///////////////Test code to test outputs
Psv.clear();
for (int ii = 0; ii < Nn; ii++)
{
	std::vector<double> pvi = P.get_vec(ii,ii+1,0,npc);
	Psv.insert(Psv.end(),pvi.begin(),pvi.end());
	countwell = Nn;
}
//////////////end test
*/

if (Psv.size() == 0) 
{
	/*
	//std::cout << "empty" << std::endl;
	if (OptCI == 'G' || OptCI == 'I')
	{
		//z = zmin;
		z = 0.0;
		return z;
	}
	else
	{
		z = 0.0;
		return z;
	}
	*/
	z = 0.0;
	return z;
}
//	std::cout << "Create Ps" << std::endl;
Ps.vec_to_mat(Psv,countwell,npc,1);
//////////Testing/////////////////

/*
int currentzero = 0;
for (int ii = 0; ii < genome.length(); ii++)
{
	std::cout << genome.gene(ii) << '\t' << obswells[genome.gene(ii) + currentzero].node << '\n';
	currentzero += ozonec[ii];
}
std::cout << std::endl;
std::cin.get();
*/
//Ps.print();
//std::exit(1);	//exit
//std::cin.get();	//pause
//std::cout << "setting test perm values" << std::endl;
//permc[0] = 15; permc[1] = 5;	//set permc to test output
//std::cout << obstimes.size() << '\t' << "nprt = " << nprt << std::endl;
//////////////end testing////////////

//std::cout << "countwell = " << countwell << '\t' << "nprt = " << nprt << std::endl;
//Get baseline
std::vector<double> permstore = permc;

std::vector<double> Ar(npc*npc), Br(npc*npc);
/*
//Pv is the vector form of P
//std::cout << "redmat base" << std::endl;
redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
spess,area,arear,bi,ci,npc,Nn,permc);

//std::cout << "red_solver base" << std::endl;
*/
std::vector<double> Jdvbase, Jdv, Jdiv;	//reduced solutions and columns of the Jacobian in vector form
mat2d Jd, I, H; //sensitivity, information, and hat matrices

//Get out baseline reduced matrices
//std::cout << "Arbuild.size() = " << Arbuild.size() << std::endl << "Brbuild.size() = " << Brbuild.size();
for (int ii = 0; ii < npc*npc; ii++)
{
	Ar[ii] = Arbuild[ii];
	Br[ii] = Brbuild[ii];
}
/*
std::clock_t start;
double duration;
	
start = std::clock();
*/
red_solver(Pv,Ar,Br,well,wellrate,npc,Nn,nq,obstimes,dtmult,Jdvbase,nprt,lumpwell,dt,Ps);
/*
duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
std::cout << "Time to run reduced model is: " << duration << std::endl;
std::exit(1);
*/
/*
///////////testcode to test outputs
mat2d s;
s.vec_to_mat(Jdvbase,countwell,nprt);
std::cout << "saving...";
s.save_to_file("test output.txt","no","test output");
//std::exit(1);
/////////end test
*/
/*
for (int ii = 0; ii < Ar.size(); ii++) {std::cout << Ar[ii] << '\t';}
std::cout << '\n' << std::endl;
for (int ii = 0; ii < Br.size(); ii++) {std::cout << Br[ii] << '\t';}
*/
//std::cout << Jdvbase.size() << std::endl;
//////////testing sensitivies
//for (unsigned int ii = 0; ii < Jdvbase.size(); ii++) {std::cout << Jdvbase[ii] << std::endl;}

/////////end testing

/*
////////////Comments on Jd/////////////////////

Jd =[ds(q1)/dH1 ds(q2)/dH1 ... ds(qn)/dH1 ds(q1)/dH2 ... ds(qn)/dHn2 ds(k1)/dq1 ... ds(k1)/dqn ds(k2)/dq1 ... ds(kn2)/dqn]
Jd is of dimesnon (nobs*countwell) x  (2*nq*nz) if you're interested in both pumping wells and zones,
if you're just intersted in pumping wells the columns are just nq, if you're intersted in zones, 
then the columns are just nz
*/
/*
int sets;
if (lumpwell == 1) {sets = 1;}
else {sets = nq;}
*/
int nobst, nptoi;
if (OptCI == 'G' || OptCI == 'I')
{
	nobst = nprt*countwell;	//number of observation points and times
	nptoi = npoi*ntoi;	//number of points and times of interest
	//std::cout << "nptoi = " << nptoi << "\tnobst = " << nobst << std::endl;
	//Jdiv.resize(npoi*ntoi*nz*nobst);
	Jdiv.resize(dI.size()*nobst);
} 

for (int ii = 0; ii < nz; ii++)
{
	/*
	std::clock_t start;
	double duration;
	start = std::clock();
	*/
	Jdv.clear();
	for (int jj = 0; jj < npc*npc; jj++)
	{
		Ar[jj] = Arbuild[(ii+1)*npc*npc + jj];
	}
	red_solver(Pv,Ar,Br,well,wellrate,npc,Nn,nq,obstimes,dtmult,Jdv,nprt,lumpwell,dt,Ps);
	//Calculate sensitivies to changes in hydraulic conductivty
	//std::cout << "expdesparam " << expdesparam << '\t' << "expdespump " << expdespump << std::endl;
	if (expdesparam == 1)
	{
		Jdiv.clear();
		Jdiv.resize(nprt*countwell);
		for (int jj = 0; jj < nprt*countwell; jj++)
		{
		//1/7/14 To change the sensitivity to take 0.01% need to replace dH with permc[ii]*.01
		//This should work because permc is fixed through the GA simulation
			Jdiv[jj] = (Jdvbase[jj] - Jdv[jj])/(permc[ii]*.01);
		}
		/*
		//This part should be outside the expdesparam loop
		if (OptCI == 'G' || OptCI == 'I')
		{
			Jdiv.insert(Jdiv.end(),dI.begin() + ii*nptoi,dI.begin() + (ii+1)*nptoi);
		//	Jd.augment_vecmat(Jdiv,(nprt*countwell+nptoi));
		}
		*/
		/*
		else
		{
			Jd.augment_vecmat(Jdiv,(nprt*countwell));
		}
		*/
		Jd.augment_vec(Jdiv);
		//Jd.print();
	}
	///////////////////////////
	/*
	NEED TO EDIT CODE AT THIS POINT FOR NEXT PROJECT
	1. Some wells need to be set known and others need to be set for exp des
	2. Could sort wells (first n_knonw are known and used for pump test) remaining are unknonw and used in exp design
	*/
	///////////////////////////
	else if (expdespump == 1)	//this part should really go outside the 0:nz loop because we don't care about those values
	{
		for (int jj = 0; jj < nprt*countwell; jj++)
		{
			Jdiv[jj] = (Jdvbase[jj])/dQ;
		}
		Jd.augment_vecmat(Jdiv,(nprt*countwell));
	}
	/*
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	std::cout << "Time to loop is: " << duration << std::endl;
	*/
}
//std::exit(1);

//Jd.print();


//If requesting T (A-opt) or E (E-opt)
if (OptCI == 'T' || OptCI == 'E')
{
	I.multiply('T',Jd,'N',Jd,I,1.0,0.0);
	z = std::abs(I.norm(OptCI));
}
else if (OptCI == 'G' || OptCI == 'I')
{
	//mat2d Q,R;
	//std::cout << "QR(Jd)\n";
	//Jd.qr(Q,R);
	mat2d U,S,Vt;
	Jd.svd(U,S,Vt);
	std::vector<double> Sv = S.get_vec(0,S.rows(),0,1);
	/*
	Jd.save_to_file("Jd.txt","no","Jd");
	std::cout << "\nU" << std::endl;
	U.print();
	std::cout << "\nVt" << std::endl;
	Vt.print();
	std::cout << "\nS" << std::endl;
	S.print();
	*/
	//d_x_e is the variance of location x due to design e;
	std::vector<double> d_x_e(nptoi,0.0), xV, x(nz);
	for (int ii = 0; ii < nptoi; ii++)
	{
		for (int jj = 0; jj < nz; jj++)
		{
			x[jj] = dI[jj*nptoi + ii];
		}
		xV = Vt.matvec_mult(1.0,'N',x,0.0,x);
		/*
		std::cout << "x\n";
		for (int jj = 0; jj < nz; jj++)
		{
			std::cout << x[jj] << '\t';
		}
		std::cout << std::endl;
		std::cout << "xV\n";
		for (unsigned int jj = 0; jj < xV.size(); jj++)
		{
			std::cout << xV[jj] << '\t';
		}
		std::cout << std::endl;
		*/
		
		for (unsigned int jj = 0; jj < xV.size();jj++)
		{
			d_x_e[ii] += (xV[jj]*xV[jj]/(Sv[jj]*Sv[jj]));
			/*
			if (!(d_x_e[ii] >= DBL_MAX && d_x_e[ii] <= -DBL_MAX))
			{
				std::cout << "exit on\t" << jj << std::endl;
				std::exit(1);
			}
			*/
		}
		
		//std::cout << "d_x_e = " << d_x_e[ii] << std::endl;
	}
//std::cout << "\nfind z" << std::endl;
	if (OptCI == 'G')
	{

		//z = -(*std::minmax_element(d_x_e.begin(),d_x_e.end()).second);
		z = 1/(*std::minmax_element(d_x_e.begin(),d_x_e.end()).second);
	}
	else
	{
		//z = -(std::accumulate(d_x_e.begin(),d_x_e.end(),0.0)/nptoi);
		z = 1/(std::accumulate(d_x_e.begin(),d_x_e.end(),0.0)/nptoi);
	}

/*
mat2d U, S, Vt;
Jd.svd(U,S,Vt);
H.multiply('N',U,'T',U,H,1.0,0.0);
diagh = H.get_diag();
	if (OptCI == 'G')
	{
		//minimization is having some isseues of convergance changing to maiximize always
		z = -(*std::minmax_element(diagh.end()-nptoi,diagh.end()).second);
	}
	else
	{
		z = -(std::accumulate(diagh.end()-nptoi,diagh.end(),0.0)/nptoi);
	}
*/
}
else if (OptCI == 'D')
{
	
	I.multiply('T',Jd,'N',Jd,I,1.0,0.0);
	/*
	std::vector<double> scale(I.rows()*I.columns(),10E10);
	mat2d scalemat, I2;
	scalemat.vec_to_mat(scale,I.rows(),I.columns());
	I2.dot('M',I,scalemat);
	*/
	//mat2d I2;
	//I2.equals(I);
	//double multconst = pow(10.0E25,(1/float(I.rows())));
	//std::cout << multconst << '\t' << I.rows() << std::endl;
	//I2.multiply_const(multconst);
	//I.multiply_const(10E10);
	I.multiply_const(multconst);
	//I2.print();
	//I2.save_to_file("I2.txt","no","I2");
	//z = std::abs(I2.determinant());
	//std::cout << z << std::endl;
	//std::exit(1);
	z = std::abs(I.determinant());
	//std::exit(1);
}

//Put something here to handle S (for SE)
//Need to add G and I optimality, both work off H(hat) matrix
//if X = USV^T then H = UU^T
//G is max(diag(H)), I is mean(diag(H))

else
{
std::cout << "No Optimality Criterion coded to handle OptCI = " << OptCI << " at this time";
std::exit(1);
}

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

return z;
}