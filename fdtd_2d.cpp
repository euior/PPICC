#include "fdtd_aloc.h"
#include "fdtd_vmap.h"
#include "fdtd_func.h"

int main()
{
Grid *g;

Particle p;

ALLOC_1D(g,1,Grid);//allocate memeory for Grid

gridInit(g);//Initialize the grid

geomInit(g);//Initialize the gemotry

CPML_coeffs(g);//Initialize the cPML

p=part_disp(g);

/* do time stepping */

for (Time = 0; Time < MaxTime; Time++) {

hardSrc(g);   // update hard source, such a laser beam

check_fperiodic2(g); //periodic boundary condition ||periodic yes

aux_cpmlh(g); //cPML phsi_e

updateH2d(g); // update magnetic field

p=part_push(p, g);

updateSrc(g); // update source

fdtd_getJ2001(g, p); //Get the current to push EM field

p=check_pperiodic2(p, g); //periodic boundary condition ||periodic yes

p=part_check(p, g);  // ||periodic no

//p=part_emit(p, g);//particle emitter, from x axis

check_Jperiodic2(g); //Rearragement current ||periodic yes

//----old version magnetic field split algorithm----//

//aux_cpmlh(g); //cPML phsi_e

//updateH2d(g); // update magnetic field

//aux_cpmle(g); //cPML phsi_h

//updateE2d(g); // update electric field

//----old version magnetic field split algorithm----//

aux_cpmle(g); //cPML phsi_h

updateE2d(g); // update electric field

//fdtd_getD(g, p);

snapshot2d(g); 

} 
return 0;

}