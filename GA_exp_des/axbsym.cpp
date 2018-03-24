#include "addon.h"

void axbsym(const int &N, const int &NTERM, std::vector<int> &TOPOL,
std::vector<int> &JA, std::vector<double> &coef1, std::vector<double> &XVEC,
std::vector<double> &BVEC)
{
	//BVEC.clear();
	//BVEC.resize(N);
/*	for (int ii = 0; ii < BVEC.size(); ii++)
	{
		BVEC[ii] = 0.0;
	}*/
	for (int k = 0; k < N; k++)
	{
		int M = TOPOL[k]-1;
		//std::cout << M << '\t';
		int MM = TOPOL[k+1]-2;
		//std::cout << MM << '\t';
		BVEC[k] = BVEC[k] + coef1[M]*XVEC[JA[M]-1];
		for (int i = M+1; i < MM+1; i++)
		{
			BVEC[k] = BVEC[k] + coef1[i]*XVEC[JA[i]-1];
			BVEC[JA[i]-1] = BVEC[JA[i]-1] + coef1[i]*XVEC[k];
		}
	}
}