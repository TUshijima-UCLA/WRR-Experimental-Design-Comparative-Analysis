#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>
#inclue<cstdlib>
#include<string>
std::vector<double> ph(0,0);
std::vector<int> wells(20,0);
wells[0] = 21833;   
wells[1] =24682;
wells[2] =21593;
wells[3] =16568;
wells[4] =25227;
wells[5] =3648;
wells[6] =9882;
wells[7] =12038;
wells[8] =21705;
wells[9] =11399;
wells[10] =12081;
wells[11] =14724;
wells[12] =11253;
wells[13] =2136;
wells[14] =6491;
wells[15] =10787;
wells[16] =17303;
wells[17] =8710;
wells[18] =7119;
wells[19] =12497;

int count = 0;
int times = 4;		//number of times output
int Nn = 29197;	//number of nodes
int nq = 20;

void read(void)
{
double nr = ceil(Nn/5);

std::ifstream file("output\sat2d.pot");
file.ignore(500, '\n');
for (int ii = 0; ii < nr+1; ii++)
{
	file.ignore(500, '\n');
}
for (int ii = 0; ii < times; ii++)
{
	file.ignore(500, '\n');
	for (int jj = 0; jj < Nn; jj++)
	{
		file >> ph[count*Nn*times + ii*Nn + jj];
	}file.ignore(500, '\n');
}
file.close();
/*
std::ofstream file2("output\satdrawtimes.txt");
for (int ii = 0; ii < times; ii++)
{
	for (int jj = 0; jj < Nn; jj++)
	{
		file2 << ph[ii*Nn + jj] << '\t';
	}	file2 << '\n';
}
*/
}

void write(void)
{
std::ofstream file("input\neu_read.sat2d");
file << "1     nq" << '\n' << wells[count] << '\n' << '\n' << "0.0  TIMEIN"
<< '\n' << "-1000" << '\n' << '\n' << "1.0d+40  TIMEIN" << '\n' << "-1000" << '\n';
}

void write2(void)
{
std::ofstream file2("output\satdrawtimes.txt")
for (int ii = 0; ii < nq; ii++)
{
	for (int jj = 0; jj < times; jj++)
	{
		for (int kk = 0; kk < Nn; kk++)
		{
			file2 << ph[ii*times*Nn + jj*Nn + kk] << '\t';
		}	file2 << std::endl;
	}
}
}

int main(void)
{
std::string program = "sat2d.exe";
ph.resize(Nn*nq*times)
for (int ii = 0; ii < nq; ii++)
{
	write();
	system(program.c_str());
	read();
	count++;
}
}