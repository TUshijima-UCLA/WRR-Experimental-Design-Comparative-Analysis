//Standard headers
#include <set>
#include <cstdlib>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <ctime>
#include <time.h>
//#include <cmath>
#include <math.h>
//Personal headers
#include "owellstruct.h"
#include "matrix.h"
#include "float.h"

#include <assert.h>

//PSO headers
#include <SwarmOps/Optimize.h>					/* Convenient function for optimizing benchmark problems. */
#include <RandomOps/Random.h>					/* Pseudo-random number generator, for seeding. */

#include <SwarmOps/Methods/Methods.h>			/* Optimization method ID-handles. */
#include <SwarmOps/Tools/Types.h>
#include <SwarmOps/Tools/Vector.h>				/* Vector operations, such as print. */

//Global variables
extern mat2d P;
extern std::vector<double> lmass, xc, yc, elstor, spess, area, arear, bi, ci, permc, Pv, wellrate, 
Jdv, obstimes, wellweight, Arbuild, Brbuild;

extern std::vector<double> lb1, ub1, dI, lbi, first_call;
extern double multconst, bestfvalue;

extern bool runlongernow;

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

double Objective(const SO_TElm *x, int n);

clock_t t3, t4;

struct ExpDesContext
{
	SO_TDim	n;	///Problem dimension
};

SO_TFitness ExpDes(const SO_TElm *x, void *context, const SO_TFitness fitnessLimit)
{
	struct ExpDesContext const* c = (struct ExpDesContext const*) context;
	double n = c->n;
	double fvalue;

	if (runlongernow == true)
	{
		SO_TElm *x2 = SO_NewVector(n);
		//first_call.resize(n);
		for (int ii = 0; ii < n; ii++)
		{
			x2[ii] = lbi[ii];
			first_call[ii] = x[ii];
			fvalue = Objective(x2, n);
		}
		runlongernow = false;
	}
	else
	{
		fvalue = Objective(x, n);
	}

	return fvalue;
};

SO_TDim ExpDesGradient(const SO_TElm *x, SO_TElm *v, void *context)
{
	//Problem is non-differentiable
	return 0;
};

SO_TElm* Optimize(int nowell, int pickgenome, double* lb, double* ub, int particle_size, float nconv)
{
	const size_t kMethodId = SO_kMethodPSOINT;
	const size_t kNumRuns = 50;
	const SO_TDim kDim = particle_size;
	int dimf = nconv;
	const size_t kDimFactor = dimf;
//#define kNumIterations (kDimFactor*kDim)
#define kNumIterations (kDimFactor)
	//std::cout << kNumIterations << std::endl;
	const char* kTraceFilename = "FitnessTrace.txt";

	struct ExpDesContext context;
	context.n = particle_size;		//particle_size = nowell if genome 2 or nobszones if genome 1
	
	struct SO_Results res;

	/* Seed the pseudo-random number generator. */
	RO_RandSeedClock(43455783);
	SO_TElm* kLowerBound = SO_NewVector(particle_size);
	SO_TElm* kUpperBound = SO_NewVector(particle_size);

	//std::cout << "set boundaries" << std::endl;
	for (int ii = 0; ii < particle_size; ii++)
	{
		kLowerBound[ii] = lb[ii];
		//std::cout << lb[ii] << std::endl;
		kUpperBound[ii] = ub[ii];
	}

	SO_TElm* kLowerInit = SO_NewVector(particle_size);
	SO_TElm* kUpperInit = SO_NewVector(particle_size);
	/*
	if (runlongernow == true)
	{
		std::cout << "set initiaion limits" << std::endl;
		for (int ii = 0; ii < particle_size; ii++)
		{
			
			//kLowerInit[ii] = ubi[ii];
			kLowerInit[ii] = kLowerBound[ii];
			std::cout << kLowerInit[ii] << std::endl;
			//kUpperInit[ii] = ubi[ii];
			kUpperInit[ii] = kLowerBound[ii];
			std::cout << kUpperInit[ii] << std::endl;
		}
	}
	else
	{
		kLowerInit = kLowerBound;
		kUpperInit = kUpperBound;
		//kUpperInit = kLowerBound;
	}
	*/

	kLowerInit = kLowerBound;
	kUpperInit = kUpperBound;


	bestfvalue = 1E100;
	printf("Method: %s\n", SO_kMethodName[kMethodId]);
	t3 = clock();
	res = SO_Optimize(kMethodId, kNumRuns, kNumIterations, 0, ExpDes, ExpDesGradient, (void*)&context,
		kDim, kLowerInit, kUpperInit, kLowerBound, kUpperBound, kTraceFilename);
	/*
	res = SO_Optimize(kMethodID, kNumRuns, kNumIterations, 0, MyProblem, MyGradient, (void*)&context, kDim,
	kLowerInit, kUpperInit, kLowerBound, kUpperBound, kTraceFilename);
	*/
	std::cout << std::endl;
	printf("Fitness average: %g\n", res.stat.fitnessAvg);
	printf("Fitness std.dev: %g\n", res.stat.fitnessStdDev);
	printf("Best fitness: %g\n", res.best.fitness);
	printf("Best solution: ");
	SO_PrintVector(res.best.x, res.best.dim);
	std::cout << std::endl;
	SO_TElm* bestx = res.best.x;
	//bestscore = res.best.fitness;
	/*
	SO_FreeResults(&res);
	SO_FreeVector(kLowerBound);
	SO_FreeVector(kUpperBound);
	SO_FreeVector(kLowerInit);
	SO_FreeVector(kUpperInit);
	*/
	return bestx;
};


void Optimize_wrapper(std::vector<double> &bestparticle, float nconv)
{

	int particle_size = bestparticle.size();
	
	double *lbp = &lb1[0];
	double *ubp = &ub1[0];
	
	
	SO_TElm* bestx = Optimize(nowell,pickgenome,lbp,ubp,particle_size,nconv);
	//std::cout << "Writing best particle" << std::endl;
	for (int ii = 0; ii < particle_size; ii++)
	{
		bestparticle[ii] = bestx[ii];
	}
}

double Objective_wrapper(std::vector<double> &particle)
{
int particle_size = particle.size();
SO_TElm* particle2 = SO_NewVector(particle_size);
for (int ii = 0; ii < particle_size; ii++)
{
	particle2[ii] = particle[ii];
}
double fvalue = Objective(particle2, particle_size);

return fvalue;
}