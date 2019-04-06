#include "fdtd_vmap.h"
#include "fdtd_aloc.h"
#include "fdtd_phyC.h"
#include <math.h>
void gridInit(Grid *g) {
double stdy;
Lamd=1.0E-6;
Lx=20.0*Lamd;
Ly=20.0*Lamd;
SizeX = 400; // x size of domain
SizeY = 400; // y size of domain

Ul=2*M_PI/Lamd;
Ut=2*M_PI*C/Lamd;

cpmlxl=0;
cpmlxr=0;
cpmlyu=0;
cpmlyd=0;

stdy=1/sqrt(4.0);
dt=Lx / SizeX / C * stdy;

MaxTime = 20000; // duration of simulation

ALLOC_2D(g->hx, SizeX, SizeY - 1, double);
ALLOC_2D(g->hxold, SizeX, SizeY - 1, double);
ALLOC_2D(g->chxh, SizeX, SizeY - 1, double);
ALLOC_2D(g->chxe, SizeX, SizeY - 1, double);
ALLOC_2D(g->mx, SizeX, SizeY - 1, double);
ALLOC_2D(g->peyz, SizeX, SizeY - 1, double);
ALLOC_1D(g->peyz_b,SizeY - 1,double);
ALLOC_1D(g->peyz_c,SizeY - 1,double);
ALLOC_1D(g->peyz_k,SizeY - 1,double);

ALLOC_2D(g->hy, SizeX - 1, SizeY, double);
ALLOC_2D(g->hyold, SizeX - 1, SizeY, double);
ALLOC_2D(g->chyh, SizeX - 1, SizeY, double);
ALLOC_2D(g->chye, SizeX - 1, SizeY, double);
ALLOC_2D(g->my, SizeX - 1, SizeY, double);
ALLOC_2D(g->pexz, SizeX - 1, SizeY, double);
ALLOC_1D(g->pexz_b, SizeX - 1,double);
ALLOC_1D(g->pexz_c, SizeX - 1,double);
ALLOC_1D(g->pexz_k, SizeX - 1,double);

ALLOC_2D(g->hz, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->hzold, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->chzh, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->chze, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->mz, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->pexy, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->peyx, SizeX - 1, SizeY - 1, double);
ALLOC_2D(g->density, SizeX - 1, SizeY - 1, double);

ALLOC_2D(g->ex, SizeX - 1, SizeY, double);
ALLOC_2D(g->cexe, SizeX - 1, SizeY, double);
ALLOC_2D(g->cexh, SizeX - 1, SizeY, double);
ALLOC_2D(g->jx, SizeX - 1, SizeY, double);
ALLOC_2D(g->phyz, SizeX - 1, SizeY, double);
ALLOC_1D(g->phyz_b,SizeY,double);
ALLOC_1D(g->phyz_c,SizeY,double);
ALLOC_1D(g->phyz_k,SizeY,double);

ALLOC_2D(g->ey, SizeX, SizeY - 1, double);
ALLOC_2D(g->ceye, SizeX, SizeY - 1, double);
ALLOC_2D(g->ceyh, SizeX, SizeY - 1, double);
ALLOC_2D(g->jy, SizeX, SizeY - 1, double);
ALLOC_2D(g->phxz, SizeX, SizeY - 1, double);
ALLOC_1D(g->phxz_b,SizeX,double);
ALLOC_1D(g->phxz_c,SizeX,double);
ALLOC_1D(g->phxz_k,SizeX,double);

ALLOC_2D(g->ez, SizeX, SizeY, double);
ALLOC_2D(g->ceze, SizeX, SizeY, double);
ALLOC_2D(g->cezh, SizeX, SizeY, double);
ALLOC_2D(g->jz, SizeX, SizeY, double);
ALLOC_2D(g->phxy, SizeX, SizeY, double);
ALLOC_2D(g->phyx, SizeX, SizeY, double);

return;
}