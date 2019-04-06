#include "fdtd_vmap.h"
#include "fdtd_phyC.h"

void geomInit(Grid *g) {

/*a description for the geometry of the 2D object */
int mm, nn;
double mium, sigm, epsj, sigj;
double Dt;

Dt=dt;

for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {
mium=MU0;
sigm=0;

Chxh(mm, nn) = (2*mium - sigm*Dt)/(2*mium + sigm*Dt);
Chxe(mm, nn) = 2*Dt/(2*mium + sigm*Dt);

}

for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY; nn++) {
mium=MU0;
sigm=0;

Chyh(mm, nn) = (2*mium - sigm*Dt)/(2*mium + sigm*Dt);
Chye(mm, nn) = 2*Dt/(2*mium + sigm*Dt);


}

for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {

mium=MU0;
sigm=0;

Chzh(mm, nn) = (2*mium - sigm*Dt)/(2*mium + sigm*Dt);
Chze(mm, nn) = 2*Dt/(2*mium + sigm*Dt);


}
for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY; nn++) {

epsj=EP0;
sigj=0;

Cexe(mm, nn) = (2*epsj - sigj*Dt)/(2*epsj + sigj*Dt);
Cexh(mm, nn) = 2*Dt/(2*epsj + sigj*Dt);


}

for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {

epsj=EP0;
sigj=0;

Ceye(mm, nn) = (2*epsj - sigj*Dt)/(2*epsj + sigj*Dt);
Ceyh(mm, nn) = 2*Dt/(2*epsj + sigj*Dt);

}

for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY; nn++) {

epsj=EP0;
sigj=0;

Ceze(mm, nn) = (2*epsj - sigj*Dt)/(2*epsj + sigj*Dt);
Cezh(mm, nn) = 2*Dt/(2*epsj + sigj*Dt);

}


return;


}