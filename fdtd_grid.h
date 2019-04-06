#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H

struct Grid {
double *hx, *hxold, *chxh, *chxe, *mx, *peyz, *peyz_b, *peyz_c, *peyz_k;
double *hy, *hyold, *chyh, *chye, *my, *pexz, *pexz_b, *pexz_c, *pexz_k;
double *hz, *hzold, *chzh, *chze, *mz, *peyx, *pexy;
double *ex, *cexe, *cexh, *jx, *phyz, *phyz_b, *phyz_c, *phyz_k;
double *ey, *ceye, *ceyh, *jy, *phxz, *phxz_b, *phxz_c, *phxz_k;
double *ez, *ceze, *cezh, *jz, *phyx, *phxy;
double *density;
double lx, ly, Dt;
double lamd;
double ul, ut;

int sizeX, sizeY;
int time, maxTime;
int pmlxl, pmlxr, pmlyu, pmlyd;

};
typedef struct Grid Grid;
#endif