#include "fdtd_vmap.h"
/* update the source */
void updateSrc(Grid *g) {

int mm, nn;


for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {
Mx(mm, nn)=0;
}

for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY; nn++) {
My(mm, nn)=0;
}

for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {

Mz(mm, nn)=0;
Density(mm, nn)=0;
}
for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY; nn++) {

Jx(mm, nn)=0;

}

for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY - 1; nn++) {

Jy(mm, nn)=0;

}

for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY; nn++) {

Jz(mm, nn)=0;

}


return;
}
