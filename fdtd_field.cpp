#include "fdtd_vmap.h"
/* update magnetic field */
void updateH2d(Grid *g) {
int mm, nn;
double dx, dy;

dx=Lx/SizeX;
dy=Ly/SizeY;



for (mm = 0; mm < SizeX; mm++){
	for (nn = 0; nn < SizeY - 1; nn++){
		Hxold(mm, nn) = Hx(mm, nn);
		Hx(mm, nn) = Chxh(mm, nn) * Hx(mm, nn)- Chxe(mm, nn) / Peyz_k(nn) / dy * (Ez(mm, nn + 1) - Ez(mm, nn)) - Chxe(mm, nn) * (Mx(mm, nn) + Peyz(mm, nn));
	}
}
for (mm = 0; mm < SizeX - 1; mm++){
	for (nn = 0; nn < SizeY; nn++){
		Hyold(mm, nn) = Hy(mm, nn);
		Hy(mm, nn) = Chyh(mm, nn) * Hy(mm, nn)+ Chye(mm, nn) / Pexz_k(mm) / dx * (Ez(mm + 1, nn) - Ez(mm, nn)) - Chye(mm, nn) * (My(mm, nn) - Pexz(mm, nn));
	}
}
for (mm = 0; mm < SizeX - 1; mm++){
	for (nn = 0; nn < SizeY - 1; nn++){
		Hzold(mm, nn) = Hz(mm, nn);		
		Hz(mm, nn) = Chzh(mm, nn) * Hz(mm, nn)+ Chze(mm, nn) / Peyz_k(nn) / dy * (Ex(mm, nn+1) - Ex(mm, nn)) - Chze(mm, nn) / Pexz_k(mm) / dx * (Ey(mm + 1, nn) - Ey(mm, nn)) - Chze(mm, nn) * (Mz(mm, nn) + Pexy(mm, nn) -Peyx(mm, nn));
	}
}

return;
}
/* update electric field */
void updateE2d(Grid *g) {
int mm, nn;
double dx, dy;

dx=Lx/SizeX;
dy=Ly/SizeY;

for (mm = 0; mm < SizeX - 1; mm++)
for (nn = 1; nn < SizeY - 1; nn++)
Ex(mm, nn) = Cexe(mm, nn) * Ex(mm, nn)+ Cexh(mm, nn) / Phyz_k(nn) / dy * (Hz(mm, nn) - Hz(mm, nn - 1)) - Cexh(mm, nn) * (Jx(mm, nn) - Phyz(mm, nn));
for (mm = 1; mm < SizeX - 1; mm++)
for (nn = 0; nn < SizeY - 1; nn++)
Ey(mm, nn) = Ceye(mm, nn) * Ey(mm, nn)- Ceyh(mm, nn) / Phxz_k(mm) / dx * (Hz(mm, nn) - Hz(mm - 1, nn)) - Ceyh(mm, nn) * (Jy(mm, nn) + Phxz(mm, nn));
for (mm = 1; mm < SizeX - 1; mm++)
for (nn = 1; nn < SizeY - 1; nn++)
Ez(mm, nn) = Ceze(mm, nn) * Ez(mm, nn)+ Cezh(mm, nn) / Phxz_k(mm) / dx * (Hy(mm, nn) - Hy(mm - 1, nn)) - Cezh(mm, nn) / Phyz_k(nn) / dy * (Hx(mm, nn) - Hx(mm, nn - 1)) - Cezh(mm, nn) * (Jz(mm, nn) - Phxy(mm, nn) + Phyx(mm, nn));

return;
}
