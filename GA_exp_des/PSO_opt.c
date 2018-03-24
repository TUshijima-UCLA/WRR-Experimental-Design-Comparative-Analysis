#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>

#include <SwarmOps/Optimize.h>					/* Convenient function for optimizing benchmark problems. */
#include <RandomOps/Random.h>					/* Pseudo-random number generator, for seeding. */

#include <SwarmOps/Methods/Methods.h>			/* Optimization method ID-handles. */
#include <SwarmOps/Tools/Vector.h>				/* Vector operations, such as print. */

#include <assert.h>

//Personal headers
#include "owellstruct.h"
#include "matrix.h"
#include "float.h"


double Objective(const SO_TElm *x, int n);

/*
struct MyContext
{
	SO_TElm a;
	SO_TElm b;
};

SO_TFitness MyProblem(const SO_TElm *x, void *context, const SO_TFitness fitnessLimit)
{
	struct MyContext const* c = (struct MyContext const*) context;
	SO_TElm a = c->a; ///from struct c get a and set as a
	SO_TElm b = c->b;

	return a * pow(x[0], 4) - b * pow(x[1], 3);
};

SO_TDim MyGradient(const SO_TElm *x, SO_TElm *v, void *context)
{

	struct MyContext const* c = (struct MyContext const*) context;

	SO_TElm a = c->a;
	SO_TElm b = c->b;

	v[0] = a * 4 * pow(x[0], 3);
	v[1] = -b * 3 * x[1];

	return 0;
};
*/


struct ExpDesContext
{
	SO_TDim	n;	///Problem dimension
};

/*
struct ExpDesContext MakeExpDesContext(SO_TDim n)
{
struct ExpDesContext c;

c.n = n;

return c;
}
*/

SO_TFitness ExpDes(const SO_TElm *x, void *context, const SO_TFitness fitnessLimit)
{
	struct ExpDesContext const* c = (struct ExpDesContext const*) context;
	double n = c-> n;
	double fvalue;
	fvalue = Objective(x,n);	
	
	return fvalue;
}

SO_TDim ExpDesGradient(const SO_TElm *x, SO_TElm *v, void *context)
{
	//Problem is non-differentiable
	return 0;
}

SO_TElm* Optimize(int nowell, int pickgenome, double* lb, double* ub, int particle_size)
{
	const size_t kMethodID = SO_kMethodPSO;
	const size_t kNumRuns = 50;
	const SO_TDim kDim = particle_size;
	const size_t kDimFactor = 200;
	#define kNumIterations (kDimFactor*kDim)
	const char* kTraceFilename = "FitnessTrace.txt";
	
	struct ExpDesContext context;
	context.n = particle_size;		//particle_size = nowell if genome 2 or nobszones if genome 1
	struct SO_Results res;
	
	/* Seed the pseudo-random number generator. */
	RO_RandSeedClock(43455783);
	SO_TElm* kLowerBound = SO_NewVector(particle_size);
	SO_TElm* kUpperBound = SO_NewVector(particle_size);
	
	for (int ii = 0; ii < particle_size; ii++)
	{
		kLowerBound[ii] = lb[ii];
		kUpperBound[ii] = ub[ii];
	}
	
	SO_TElm* kLowerInit = kLowerBound;
	SO_TElm* kUpperInit = kUpperBound;
	
	res = SO_Optimize(kMethodID, kNumRuns, kNumIterations, 0, ExpDes, ExpDesGradient, (void*) &context,
						kDim, kLowerInit, kUpperInit, kLowerBound, kUpperBound, kTraceFilename);
	/*
	res = SO_Optimize(kMethodID, kNumRuns, kNumIterations, 0, MyProblem, MyGradient, (void*)&context, kDim,
		kLowerInit, kUpperInit, kLowerBound, kUpperBound, kTraceFilename);
	*/
	printf("Fitness average: %g\n", res.stat.fitnessAvg);
	printf("Fitness std.dev: %g\n", res.stat.fitnessStdDev);
	printf("Best fitness: %g\n", res.best.fitness);
	printf("Best solution: ");
	SO_PrintVector(res.best.x, res.best.dim);
	SO_TElm* bestx = res.best.x;
	//bestscore = res.best.fitness;
	SO_FreeResults(&res);
	SO_FreeVector(kLowerBound);
	SO_FreeVector(kUpperBound);
	SO_FreeVector(kLowerInit);
	SO_FreeVector(kUpperInit);
	
	return bestx;
};