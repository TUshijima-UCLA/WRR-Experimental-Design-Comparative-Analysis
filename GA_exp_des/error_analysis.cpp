#include "matrix.h"
//#include "matrix.cpp"
#include <iostream>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <random>

int ncalls = 0;


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

void red_solver(std::vector<double> &P, std::vector<double> &Ar, std::vector<double> &Br,
std::vector<int> &well, std::vector<double> &wellrate, int &npc, int Nn, int &nq,
std::vector<double> &ts, double &dtmult, std::vector<double> &Jdv, int nprt, int lumpwell, 
double dt, mat2d &Ps);

void read_names(std::string &pname, std::string &potname,
std::string &matname, int &nz, int &Nn, std::vector<double> &material)
{
std::string gname;
std::ifstream file("sat2d.fnames");
if (file.is_open())
{
	file.seekg(0);
	file.ignore(500, '\n');		//ignore first line	(ITERM)
	file.ignore(500, '/');
	file >> pname;				//get parameter file name (unit 5)
	file.ignore(500, '\n');	
	file.ignore(500, '/');
	file >> gname;				//grid info
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '/');
	file >> matname;
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '\n');	
	file.ignore(500, '/');
	file >> potname;				//get potential head file name (unit 40)
	//std::cout << pname << '\t' << gname << '\t' << matname << '\t'<< potname << std::endl;
	file.close();
}
else
{std::cout << "Could not open fnames file to read_1";
//exitIO = 1;
std::exit (EXIT_FAILURE);
}
std::ifstream file1(("input\\"+gname).c_str());
if (file1.is_open())
{
	file1.seekg(0);
	file1 >> nz; //get # hydraulogic zones
//	std::cout << nz << std::endl;
	file1.ignore(500, '\n');
	file1 >> Nn;
	file1.close();
}
else
{std::cout << "Could not open grid file to read_1";
std::exit (EXIT_FAILURE);
}

std::ifstream file8(("input\\store_"+matname).c_str());
//std::ifstream file8("store_material_real_3.sat2d");
if (file8.is_open())
{
	double temp1, temp2;
	material.resize(2*nz);
//	std::cout << material.size() << std::endl;
	for (int ii = 0; ii < nz; ii++)
	{
		file8 >> temp1;
		file8 >> temp2;
		file8 >> material[2*ii];
		file8 >> material[2*ii + 1];
		file8.ignore(500, '\n');
	}
	file8.close();
}
else
{
	std::cout << "Could not open material file to read_1";
	std::exit (EXIT_FAILURE);
}

//std::cout << nz << '\t' << Nn << '\t' << material[0] << '\t' << material[1] << '\t' << material[2] << '\t' << material[3] << std::endl;
//std::exit(1);

}

void write_material(std::string &matname, std::vector<double> permc, std::vector<double> &material)
{

std::ofstream file(("input\\"+matname).c_str());
if (file.is_open())
{
	int nz = permc.size();
	for (int ii = 0; ii < nz; ii++)
	{
		file << std::setw(13) << permc[ii] << std::setw(13) << permc[ii] 
		<< std::setw(13) << material[2*ii] 
		<< std::setw(13) << material[2*ii + 1] << std::endl;
	}
	file << '\n' << std::setw(13) << "Perm x" << std::setw(13) << "Perm y"
	<< std::setw(13) << "Elstor" << std::setw(13) << "Depth" << std::endl;
}
else
{
	std::cout << "Cannot open material (hydraulic condutivity) file to write" << std::endl;
	std::exit (EXIT_FAILURE); //Exit upon failure of I/O process
}
//std::exit(1);
}


void read_head(std::string &potname, mat2d &sshots, int &Nn, int nt, int count)
{
std::ifstream file(("output\\"+potname).c_str());
if (file.is_open())
{	
//	std::cout << "pot open" << std::endl;
	double temp;
	mat2d tempmat(Nn,nt);
	file.seekg(0);
	file.ignore(500, '\n');
	for (int ii = 0; ii < ceil(Nn/5.0); ii++)
	{
		file.ignore(500, '\n');
	}
	for (int jj = 0; jj < nt; jj++)
	{	file.ignore(500, '\n');
		for (int ii = 0; ii < Nn; ii++)
		{
			//file >> temp;
//			std::cout << temp << std::endl;
//			sshots.push_back(temp);
//			std::cout << sshots[ii] << std::endl;
//			std::cout << sshots.size() << std::endl;
//			file >> sshots[(count-1)*Nn*nt + jj*Nn + ii];
			file >> temp;
//			std::cout << temp << '\t';
			tempmat.set_value(ii,jj,temp);
		}
		file.ignore(500, '\n');
	}
	file.close();
//	std::cout << "copy tempmat" << std::endl;
	if (count == 0) {sshots.equals(tempmat);}
	else {sshots.augment_mat(tempmat);}
}
else
{std::cout << "Could not open potential head file to read_3";
std::exit (EXIT_FAILURE);
}
}

void write_results(std::vector<double> &RMSE, std::vector<double> &perror, mat2d &red, mat2d &full,
std::vector<double> &permstore, int nz)
{
std::ofstream file("Result_error.txt");
if (file.is_open())
{
std::cout << "write error" << std::endl;
	for (unsigned int ii = 0; ii < RMSE.size(); ii++)
	{
		file << std::setprecision(4) << RMSE[ii] << '\t' << perror[ii] << '\t';
		for (int jj = 0; jj < nz; jj++)
		{
			file << std::setprecision(6) << permstore[jj+ii*nz] << '\t';
		}
		file << std::endl;
	}
	file.close();
}
else
{
	std::cout << "Could not open Result_error file to write";
	std::exit(1);
}

std::cout << "write red" << std::endl;
red.save_to_file("Result_error.txt","app","Reduced");
/*
std::cout << "write full" << std::endl;
full.save_to_file("Result_error.txt","app","Full");
*/
}


int main()
{
mat2d P,sshots;
std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, Pv, 
Jdv, obstimes, wellweight, Arbuild, Brbuild, material;
std::vector<int> topol, JA, triang, trija, well, ozonec;
//std::vector<owell> obswells;			//all potential obsrvation nodes
int ntri, irad, nterm, npc, Nn, nprt = 1, nq, nz;	
std::string pname, potname, matname;

read_fort(ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,spess,area,arear,bi,ci,JA,topol);
read_names(pname,potname,matname,nz,Nn,material);
read (npc,Nn,nq,well,wellweight,P,nz);


std::vector<double> permc(nz,0);
/*model specific parameters permstore(number of k's to check error)
wellrate (well rate), ts (time to take error measurement), dtmult, dt
*/
int nkcheck = 1;
std::vector<double> wellrate(nq,0), ts(1,50);
for (int ii = 0; ii < nq; ii++) {wellrate[ii] = -wellweight[ii];}
double dtmult = 1.1, dt = 0.1;

std::vector<double> permstore(nz*nkcheck,0);
unsigned seed = 5;
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> distribution(0.1,20.0);
//Projection matrix in vector form
Pv = P.get_vec(0,Nn,0,npc);
std::string program = "sat2d.exe";
std::vector<double> RMSE (nkcheck,0), perror(nkcheck,0);
mat2d sshotssave, redshotssave;

// 12-7-14 edit to test error in full and reduced sat2d model
//nkcheck = 1;

for (int ii = 0; ii < nkcheck; ii++)
{
	std::cout << "Parameter Number: " << ii << std::endl;
	///////////Testing/////////////
	/*
	permc[0] = 15;
	permc[1] = 5;
	*/
	///////////////////////////////	
	
	for (int jj = 0; jj < nz; jj++)
	{
		permc[jj] = distribution(generator);
		permstore[jj + ii*nz] = permc[jj];
	}
	/*
	//randomly generated hydraulic conductivity values
	
	permc[0] = 13.188;
	permc[1] = 5.936;
	permc[2] = 1.986;
	permc[3] = 7.947;
	permc[4] = 13.502;
	permc[5] = 11.982;
	permc[6] = 3.540;
	*/
	permc[0] = 20;
	permc[1] = 20;
	permc[2] = 0.1;

	write_material(matname,permc,material);
	system(program.c_str());
	
	sshots.clear();
	
	read_head(potname, sshots, Nn, 1, 0);
	
	//std::cout << sshots.rows() << '\t' << sshots.columns() << std::endl;
	
	
	std::vector<double> Ar(npc*npc), Br(npc*npc), redshotsv;
	
	redmat(Pv,Ar,Br,topol,JA,ntri,irad,nterm,triang,trija,lmass,xc,yc,elstor,
	spess,area,arear,bi,ci,npc,Nn,permc);
	
	red_solver(Pv,Ar,Br,well,wellrate,npc,Nn,nq,ts,dtmult,redshotsv,nprt,1,dt,P);
	/*
	mat2d hbase;
	hbase.vec_to_mat(redshotsv,Nn,1,2);
	hbase.save_to_file("redshotsv.txt","no","redshotsv");
	*/


	std::vector<double> sshotsv = sshots.get_vec(0,Nn,0,1);
	std::vector<double> diffv(Nn,0);
	for (int kk = 0; kk < Nn; kk++)
	{
		diffv[kk] = sshotsv[kk] - redshotsv[kk];
		RMSE[ii] += diffv[kk]*diffv[kk];
		if (std::abs(sshotsv[kk]) > 1E-2)
		{
			perror[ii] += diffv[kk]*diffv[kk]/(sshotsv[kk]*sshotsv[kk]);
		}
		//std::cout << sshotsv[kk] << std::endl;
	}
	
	RMSE[ii] = sqrt(RMSE[ii])/(Nn);
	perror[ii] = sqrt(perror[ii])/(Nn);
	
	if (ii == 0){sshotssave.equals(sshots);}
	else{sshotssave.augment_mat(sshots);}
	redshotssave.augment_vec(redshotsv);
	
}

// 12-7-14 edit to test RMSE of full sat2d model
//std::cout << "RMSE = " << RMSE[0] << std::endl;
//std::exit(1);

std::cout << "write results" << std::endl;
write_results(RMSE,perror,redshotssave,sshotssave, permstore, nz);

}