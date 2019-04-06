#include "fdtd_vmap.h"
void aux_cpmlh(Grid *g){
int mm, nn;
double dx, dy;

dx = Lx / SizeX;
dy = Ly / SizeY;


for (mm = 0; mm < SizeX; mm++)
for (nn = 0; nn < SizeY - 1; nn++)
Peyz(mm, nn) = Peyz_b(nn) * Peyz(mm, nn) + Peyz_c(nn) * (Ez(mm, nn + 1) - Ez(mm, nn)) / dy;
for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY; nn++)
Pexz(mm, nn) = Pexz_b(mm) * Pexz(mm, nn) + Pexz_c(mm) * (Ez(mm + 1, nn) - Ez(mm, nn)) / dx;
for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY - 1; nn++)
{
Peyx(mm, nn) = Peyz_b(nn) * Peyx(mm, nn) + Peyz_c(nn) * (Ex(mm, nn + 1) - Ex(mm, nn)) / dy;
Pexy(mm, nn) = Pexz_b(mm) * Pexy(mm, nn) + Pexz_c(mm) * (Ey(mm + 1, nn) - Ey(mm, nn)) / dx;
}

return;

}


void aux_cpmle(Grid *g){
int mm, nn;
double dx, dy;

dx = Lx / SizeX;
dy = Ly / SizeY;



for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 1; nn < SizeY - 1; nn++)
Phyz(mm, nn) = Phyz_b(nn) * Phyz(mm, nn) + Phyz_c(nn) * (Hz(mm, nn) - Hz(mm, nn - 1)) / dy;
for (mm = 1; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY - 1; nn++)
Phxz(mm, nn) = Phxz_b(mm) * Phxz(mm, nn) + Phxz_c(mm) * (Hz(mm, nn) - Hz(mm - 1, nn)) / dx;
for (mm = 1; mm < SizeX - 1; mm++)
for (nn = 1; nn < SizeY - 1; nn++)
{
Phyx(mm, nn) = Phyz_b(nn) * Phyx(mm, nn) + Phyz_c(nn) * (Hx(mm, nn) - Hx(mm, nn - 1)) / dy;
Phxy(mm, nn) = Phxz_b(mm) * Phxy(mm, nn) + Phxz_c(mm) * (Hy(mm, nn) - Hy(mm - 1, nn)) / dx;
}

return;

}