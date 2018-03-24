//Standard headers
#include <vector>
#include <iostream>
#include <cmath>

//My headers
#include "owellstruct.h"
#include "matrix.h"

//GA headers
#include <ga/ga.h>
#include <ga/std_stream.h>
//Call real GA
#include <ga/GARealGenome.h>

void wellsame(std::vector<robustowell> &rowells, std::vector<owell> &copt,
double &close, double score, GARealGenome &g, int combindex)
/*
void wellsame(std::vector<robustowell> &rowells, std::vector<owell> &copt,
	double &close, double score, std::vector<double> &particle, int combindex)
*/
{
double d2;
//double close = 2;
bool isclose = false;
unsigned int ii = 0;	//ii searches the list of unique optimal solutions

while (ii < rowells.size())
{
	//std::cout << "ii = " << ii << std::endl;
	unsigned int jj = 0;	//jj searches the nodes of a stored unique solution
	unsigned int kk = 0;	//kk searches the nodes of the solution currently under consideration
	std::vector<owell> temp = rowells[ii].optowells;
	//std::cout << "temp.size() = " << temp.size() << '\t' << "copt.size() = " << copt.size() << std::endl;
	if (temp.size() == copt.size())
	{
		while (jj < temp.size())
		{
			//std::cout << "temp.size() = " << temp.size() << std::endl;
			//std::cout << "jj = " << jj << std::endl;
			isclose = false;
			//std::cout << "kk = " << kk << std::endl;
			d2 = sqrt((copt[kk].x - temp[jj].x)*(copt[kk].x - temp[jj].x) + (copt[kk].y - temp[jj].y)*(copt[kk].y - temp[jj].y));
			if (d2 <= close && temp.size() > 0)	//if it's close erase mark off jj and move to next kk
			{
				isclose = true;
				temp.erase(temp.begin() + jj);
				kk++;
				jj = 0;
			}
			else {jj++;}	//move on to next jj
		}
	}
	else {isclose = false;}
	if (isclose == true)	//if each kk has a corresponding "close" jj, the solution is "close" -> cataloge occurance
	{
		
		rowells[ii].maxscore.push_back(score);
		if (copt[0].highfreq == true)
		{
			rowells[ii].highfreqwells = true;
		}
		else
		{
			rowells[ii].count++;	//only count GA solutions
			rowells[ii].paramindex.push_back(combindex);	//record param combinations that gave you this network
		}
		ii = rowells.size();
	}
	else {ii++;}	//if it's not close to the ii unique solution, move on to ii+1 unique solution
}

if (isclose == false)	//if it's not close to any unique solution, store as new unique solution 
{
	robustowell temp2;
	temp2.optowells = copt;					//store optimal solution
	if (copt[0].highfreq == false)			//only count GA solutions
	{
		temp2.count = 1;					//count once
		temp2.paramindex.push_back(combindex);
	}
	else {temp2.count = 0;}
	temp2.maxscore.push_back(score);		//store score
	temp2.minscore = score;					//set minscore
	temp2.highfreqwells = copt[0].highfreq;	//check to see if it's the highest frequency for well occurance
	
	for (int ii = 0; ii < g.length(); ii++)	//store genome
	{
		temp2.genome.push_back(g.gene(ii));
	}
	
	/*
	for (unsigned int ii = 0; ii < particle.size(); ii++)
	{
		temp2.genome.push_back(particle[ii]);
	}
	*/
	//temp2.genome = &g;
	//temp2.worstgenome = &g;
	rowells.push_back(temp2);
}

}

void countowells(const std::vector<robustowell> &rowells, std::vector<owell> &allowell,
double close, int pickgenome)
{
allowell[0] = rowells[0].optowells[0];
allowell[0].count = 0;
//std::cout << allowell[0].zonecount << std::endl;
if (pickgenome == 1)
{
	allowell[0].gene = rowells[0].genome[allowell[0].zonecount];
}
else
{
	allowell[0].gene = rowells[0].genome[0];
}
//std::cout << "zone = " << allowell[0].zone << std::endl;
for (unsigned int ii = 0; ii < rowells.size(); ii++)
{
	/*
	for (unsigned int jj = 0; jj < rowells[ii].genome.size(); jj++)
	{
		std::cout << rowells[ii].genome[jj] << '\t';
	}
	std::cout << '\n' << std::endl;
	*/
	//std::vector<owell> tempowell = rowells[ii].optowells;
	for (unsigned int jj = 0; jj < rowells[ii].optowells.size(); jj++)
	{
		for (unsigned int kk = 0; kk < allowell.size(); kk++)
		{
			double dist = sqrt((allowell[kk].x - rowells[ii].optowells[jj].x)*(allowell[kk].x - rowells[ii].optowells[jj].x) + 
			(allowell[kk].y - rowells[ii].optowells[jj].y)*(allowell[kk].y - rowells[ii].optowells[jj].y));
			
			if (dist <= close)
			{
				allowell[kk].count += rowells[ii].count;
				break;	//if you found one it matches stop and move on
			}
			else if (dist > close && kk == (allowell.size()-1))
			{
				allowell.push_back(rowells[ii].optowells[jj]);
				allowell[kk+1].count = rowells[ii].count;
				if (pickgenome == 1)
				{
					allowell[kk+1].gene = rowells[ii].genome[allowell[kk+1].zonecount];
				}
				else
				{
					allowell[kk+1].gene = rowells[ii].genome[jj];
				}
				//allowell[kk+1].zone = rowells[ii].optowells[jj].zone;
				//std::cout << allowell[kk+1].gene << '\t' << allowell[kk+1].zone << std::endl;
				break;	//break so the "push_back" doesn't enter an infinite loop
			}
		}
	}
}
/*
for (unsigned int ii = 0; ii < allowell.size(); ii++) {std::cout << allowell[ii].gene << '\t';}
std::cout << std::endl;
*/
//returns all observation wells sorted by occurance;
std::sort(allowell.begin(),allowell.end(),by_count_owell());
/*
for (unsigned int ii = 0; ii < allowell.size(); ii++) {std::cout << allowell[ii].gene << '\t';}
std::cout << std::endl;
*/
}