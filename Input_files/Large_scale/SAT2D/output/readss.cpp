#include<iostream>
#include<fstream>
#include<vector>
#include<cmath>

int main(void)
{
int Nn = 29197;
double nr = ceil(Nn/5);
std::vector<double> ph(Nn,0);
std::ifstream file("sat2d.pot");
file.ignore(500, '\n');
for (int ii = 0; ii < nr+1; ii++)
{
	file.ignore(500, '\n');
}
file.ignore(500, '\n');
for (int ii = 0; ii < Nn; ii++)
{
	file >> ph[ii];
}
file.close();

std::ofstream file2("satdraw.txt",std::ios::app);
{
	file2 << '\n';
	for (int ii = 0; ii < Nn; ii++)
	{
		file2 << ph[ii] << '\t';
	}
}
}