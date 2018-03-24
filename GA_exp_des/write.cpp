//Standard headers
#include <vector>
#include <iostream>
#include <iomanip>      // std::setprecision
#include <string>
#include <sstream>

//My headers
#include "owellstruct.h"
#include "matrix.h"

extern char OptCI;

void write_final(int nowell, std::vector<robustowell> &rowells, int ncalls,
int nz, std::vector<owell> &allowell)
{
//Print out final solutions
std::ostringstream convert;
convert << nowell;
//std::string snowell = convert.str();
std::string filename = ("Robust solutions for "+convert.str()+" observation wells_"+OptCI+".txt").c_str();

std::ofstream file (filename);
if (file.is_open())
{
	file << "Total number of buckets:" << '\n' << rowells.size() << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		//Write to file
		file << "Number of bucket " << ii << " encountered: " << rowells[ii].count << std::endl;
		if (rowells[ii].highfreqwells == true) {file << "Highest frequency observation wells" << std::endl;}
		if (rowells[ii].highfreqnetwork == true) {file << "Highest frequency observation network" << std::endl;}
		if (rowells[ii].robustdesign == true) {file << "Robust Design" << std::endl;}
		//std::cout << rowells[ii].highfreqwells << '\t' << rowells[ii].highfreqnetwork <<'\t' << rowells[ii].robustdesign << std::endl;
		file << "Size of bucket: " << rowells[ii].optowells.size() << std::endl;
	
		file << std::resetiosflags(std::ios_base::scientific) << "All Scores:" << std::endl;
		for (unsigned int jj = 0; jj < rowells[ii].maxscore.size(); jj++)
		{
			file << rowells[ii].maxscore[jj] << '\t';
		}
		
		file << '\n' << "Worst Case Scenario score: " << rowells[ii].minscore << std::endl
		<< "Worst Case Hydraulic Conductivity:" << std::endl;
		
		for (int jj = 0; jj < nz; jj++) {file  << 
		rowells[ii].worstcaseK[jj] << '\t';}
		file << std::endl;
		for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
		{
			file << std::scientific << std::setprecision(12) << "x = " 
			<< rowells[ii].optowells[jj].x << '\t' << "y = " 
			<< rowells[ii].optowells[jj].y << '\n';
		}
		file << std::endl;
		//Print to screen
		std::cout << "Number of bucket " << ii << " encountered is " << rowells[ii].count << std::endl;
		if (ii == 0) {std::cout << "Highest frequency observation wells" << std::endl;}
		std::cout << "Size of bucket is " << rowells[ii].optowells.size() << std::endl;
		
		std::cout << std::resetiosflags(std::ios_base::scientific) << "All Scores:" << std::endl;
		for (unsigned int jj = 0; jj < rowells[ii].maxscore.size(); jj++)
		{
			std::cout << rowells[ii].maxscore[jj] << '\t';
		}
		
		std::cout << '\n' << "Worst Case Scenario score is " << rowells[ii].minscore << std::endl
		<< "Worst Case Hydraulic Conductivity is" << std::endl;
		
		for (int jj = 0; jj < nz; jj++) {std::cout <<
		rowells[ii].worstcaseK[jj] << '\t';}
		std::cout << std::endl;
		for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
		{
			std::cout << std::scientific << std::setprecision(12) << "x = " 
			<< rowells[ii].optowells[jj].x << '\t' << "y = " 
			<< rowells[ii].optowells[jj].y << '\n';
		}
		std::cout << std::endl;
		/*
		std::cout << "Number of bucket " << ii << " encountered: " << rowells[ii].count << std::endl
		<< "Size of bucket: " << rowells[ii].optowells.size() << std::endl
		<< "Worst Case Scenario score: " << rowells[ii].minscore << std::endl
		<< "Worst Case Hydraulic Conductivity:" << std::endl;
		
		
		
		
		std::cout << '\n' << "Observation well coordinates:" << std::endl;
		for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
		{
			std::cout << "x = " << rowells[ii].optowells[jj].x << '\t' << "y = " << rowells[ii].optowells[jj].y << std::endl;
		}
		*/
	}
	
	
	//Print out total number of reduced model calls
	file << '\n' << "Total number of reduced model calls: " << ncalls << '\n' << std::endl;
	std::cout << '\n' << "Total number of reduced model calls is " << ncalls << '\n' << '\n' << std::endl;
	
	//Print all wells to file
	std::cout << "List of all observation wells (node, zone, # occurances, x, y)" << std::endl;
	file << "List of all observation wells (node, zone, # occurances, x, y):" << std::endl;
	for (unsigned int ii = 0; ii < allowell.size(); ii++)
	{
		std::cout << allowell[ii].node << '\t' << allowell[ii].zone << '\t' << 
		allowell[ii].count << '\t' << allowell[ii].x << '\t' << allowell[ii].y << std::endl; 
		file << allowell[ii].node << '\t' << allowell[ii].zone << '\t' << 
		allowell[ii].count << '\t' << allowell[ii].x << '\t' << allowell[ii].y << std::endl;
	}
	//std::cout << "close file" << std::endl;
	file.close();
}
else
{
	std::cout << "Cannot open file to write solutions";
}

//Print out Restart file with only robust solution

std::string filename2 = ("Robust Restart for " + convert.str() + " observation wells_"+OptCI+".txt").c_str();

std::ofstream file2(filename2);
if (file2.is_open())
{
	file2 << 1 << '\t' << nowell << '\t' << 0
		<< '\t' << rowells[0].genome.size() << '\t' << ncalls << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		if (rowells[ii].robustdesign == true)
		{
			file2 << rowells[ii].count << std::endl;
			for (int jj = 0; jj < rowells[ii].count; jj++)
			{
				file2 << rowells[ii].paramindex[jj] << '\t';
			}
			file2 << std::endl;
			file2 << rowells[ii].genome.size() << '\n';
			for (unsigned int jj = 0; jj < rowells[ii].genome.size(); jj++)
			{
				file2 << rowells[ii].genome[jj] << '\n';
			}
			break;
		}
	}
	file2.close();
}
else
{
	std::cout << "Cannot open file to write robust restart";
}

//Save Final Restart file
std::string filename3 = ("Restart for " + convert.str() + " observation wells_" + OptCI + ".txt").c_str();
std::ofstream file3(filename3);
if (file3.is_open())
{
	file3 << rowells.size() << '\t' << nowell << '\t' << 0
		<< '\t' << rowells[0].genome.size() << '\t' << ncalls << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		//Write to file
		file3 << rowells[ii].count << std::endl;
		for (int jj = 0; jj < rowells[ii].count; jj++)
		{
			file3 << rowells[ii].paramindex[jj] << '\t';
		}
		file3 << std::endl;
		/*
		for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
		{
		file << std::scientific << std::setprecision(12) << rowells[ii].optowells[jj].x << '\t' << rowells[ii].optowells[jj].y << std::endl;
		}
		*/
		file3 << rowells[ii].genome.size() << '\n';
		for (unsigned int jj = 0; jj < rowells[ii].genome.size(); jj++)
		{
			file3 << rowells[ii].genome[jj] << '\n';
		}

	}
	file3.close();
	//permcheck.save_to_file("Restart_Exp_Des.txt","app","Permcheck");
}
else
{
	std::cout << "Could not open file to save final restart data" << std::endl;
}

}

void write_restart(int nowell, const std::vector<robustowell> &rowells, 
int ncalls, int kleft)
{
std::ofstream file("Restart_Exp_Des.txt");
if (file.is_open())
{
	file << rowells.size() << '\t' << nowell << '\t' << kleft
	<< '\t' << rowells[0].genome.size() << '\t' << ncalls << std::endl;
	for (unsigned int ii = 0; ii < rowells.size(); ii++)
	{
		//Write to file
		file << rowells[ii].count << std::endl;
		for (int jj = 0; jj < rowells[ii].count; jj++)
		{
			file << rowells[ii].paramindex[jj] << '\t';
		}
		file << std::endl;
		/*
		for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
		{
			file << std::scientific << std::setprecision(12) << rowells[ii].optowells[jj].x << '\t' << rowells[ii].optowells[jj].y << std::endl;
		}
		*/
		file << rowells[ii].genome.size() << '\n';
		for (unsigned int jj = 0; jj < rowells[ii].genome.size(); jj++)
		{
			file << rowells[ii].genome[jj] << '\n';
		}
		
	}
	file.close();
	//permcheck.save_to_file("Restart_Exp_Des.txt","app","Permcheck");
}
else
{
	std::cout << "Could not open file to save restart data" << std::endl;
}
}