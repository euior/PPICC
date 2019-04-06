#ifndef _FDTD_FUNC_H
#define _FDTD_FUNC_H
#include "fdtd_grid.h"
#include "part_init.h"
/* Function prototypes */
void gridInit(Grid *g);
void geomInit(Grid *g);
void snapshot2d(Grid *g);
void updateE2d(Grid *g);
void updateH2d(Grid *g);
void updateSrc(Grid *g);
void hardSrc(Grid *g);

/* Functions for cPML*/
void CPML_coeffs(Grid *g);
void aux_cpmlh(Grid *g);
void aux_cpmle(Grid *g);

/* Funcitons for particles*/


Particle part_disp(Grid *g);
Particle part_push(Particle ptc, Grid *g);
Particle part_check(Particle pt, Grid *g);
Particle check_pperiodic2(Particle pt, Grid *g);
Particle part_emit(Particle pp, Grid *g);


void part_snap(Particle ph, Grid *g);
void check_Jperiodic2(Grid *g);
void check_fperiodic2(Grid *g);
/* Connection */
void fdtd_getJ2001(Grid *g, Particle pat);
void fdtd_getD(Grid *g, Particle pat);
#endif