#ifndef _PART_FUNC_H_
#define _PART_FUNC_H_
#include "part_init.h"
#include "part_parm.h"
#include "fdtd_vmap.h"
int getI(double D);
Particle part_init(int N);
Particle part_add(Particle pHead, Particle pAdd);
Particle part_del(Particle pHead, Particle pDel);
double part_density(Grid *g, int sp, double x, double y);
double emit_flux(Grid *g, int s, double vo, double x, double t);
double gauss();
double one();
#endif