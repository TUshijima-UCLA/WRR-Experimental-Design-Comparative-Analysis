#ifndef fort_h
#define fort_h


#include <iostream>
#include <fstream>
#include <memory>

extern "C" {
void datin_(int &iterm, int *triang, double *x, double *y, double *permx,
double *permy, double *elstor, double *spess, int *contp, int *conq, int *contr,
double *timprt, double *ptimep, double &tetaf, double &deltat, double &dtmax,
double &tmax, double &dtmaga, double &dtmagm, int &itmxcg, double &tolcg,
int &lump, int &irad, int &iprt1, int &iprt, int &indp, int &nq, int &nzone,
int &n1, int &nr, int &nprt, int &n, int &np, int &nt);
}

extern "C" {
void openio_(int &iterm);
}

extern "C" {
void closio_(int &iterm);
}
/*
extern "C" {
void assf2d_(int &nt, int &irad, int *triang, int *trija, double *lmass,
double *coef1, double *coef2, double *xc, double *xy, double *permx,
double *permy, double *elstor, double *spess, double *area, double *arear,
double *bi, double *ci);
}
*/
extern "C"{
void arebas_(int &n, int &nt, int *triang, int *ip3, double *bi,
double *ci, double *areanod, double *area, double *arear, int *iarea,
double *x, double *y, int &iterm);
}

extern "C"{
void init0r_(int &num, double *rvec);
}

extern "C"{
void areas_(int *ip3, int *triang, double *x, double *y, double &are);
}

extern "C"{
void basis2_(int *ip3, int *triang, double *x, double *y, double *bi,
double *ci);
}

extern "C"{
void strpic_(int &n, int &nterm, int *triang, int *ja, int *topol,
int &nt, int &n1, int &imax, int &iterm);
}

extern "C"{
void chkpic_(int &n, int &nterm, int *topol, int *ja, int &ndz, int *ier);
}

extern "C"{
void tripic_(int &nt, int *triang, int *ja, int *topol, int *trija);
//void tripic_(int &nt, std::unique_ptr<int> *triang, int *ja, int *topol, int *trija);
}

extern "C"{
void locmas_(double *lmass, int &lump);
//void locmas_(std::unique_ptr<double> &lmass, int &lump);
}

extern "C"{
void nodelt_(int &nt, int *triang, double *vnod, double *velt);
}

extern "C"{
void init0i_(int &num, int *ivec);
}
#endif