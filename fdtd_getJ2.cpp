#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void fdtd_getJ2(Grid *g, Particle pat){
int ii, jj;
double dx, dy, xx, yy;
double wx1, wx2, wx3, wy1, wy2, wy3;
double gamma;
dx=Lx/SizeX;
dy=Ly/SizeY;
while(pat){
	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);
	ii=floor((pat->x - dx) / dx + 1 / 2);
	jj=floor((pat->y - dy / 2) / dy + 1 / 2);
	xx=(pat->x - dx) / dx - ii;
	yy=(pat->y - dy / 2) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	Jx(ii - 1, jj - 1)+=pat->q/dx/dy* wx1 * wy1 * pat->px/gamma*C;
	Jx(ii - 1, jj)+=pat->q/dx/dy* wx1 * wy2 * pat->px/gamma*C;
	Jx(ii - 1, jj + 1)+=pat->q/dx/dy* wx1 * wy3 * pat->px/gamma*C;
	Jx(ii, jj - 1)+=pat->q/dx/dy* wx2 * wy1 * pat->px/gamma*C;
	Jx(ii, jj)+=pat->q/dx/dy* wx2 * wy2 * pat->px/gamma*C;
	Jx(ii, jj + 1)+=pat->q/dx/dy* wx2 * wy3 * pat->px/gamma*C;
	Jx(ii + 1, jj - 1)+=pat->q/dx/dy* wx3 * wy1 * pat->px/gamma*C;
	Jx(ii + 1, jj)+=pat->q/dx/dy* wx3 * wy2 * pat->px/gamma*C;
	Jx(ii + 1, jj + 1)+=pat->q/dx/dy* wx3 * wy3 * pat->px/gamma*C;



	ii=floor((pat->x - dx / 2) / dx + 1 / 2);
	jj=floor((pat->y - dy) / dy + 1 / 2);
	xx=(pat->x - dx / 2) / dx - ii;
	yy=(pat->y - dy) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	Jy(ii - 1, jj - 1)+=pat->q/dx/dy* wx1 * wy1 * pat->py/gamma*C;
	Jy(ii - 1, jj)+=pat->q/dx/dy* wx1 * wy2 * pat->py/gamma*C;
	Jy(ii - 1, jj + 1)+=pat->q/dx/dy* wx1 * wy3 * pat->py/gamma*C;
	Jy(ii, jj - 1)+=pat->q/dx/dy* wx2 * wy1 * pat->py/gamma*C;
	Jy(ii, jj)+=pat->q/dx/dy* wx2 * wy2 * pat->py/gamma*C;
	Jy(ii, jj + 1)+=pat->q/dx/dy* wx2 * wy3 * pat->py/gamma*C;
	Jy(ii + 1, jj - 1)+=pat->q/dx/dy* wx3 * wy1 * pat->py/gamma*C;
	Jy(ii + 1, jj)+=pat->q/dx/dy* wx3 * wy2 * pat->py/gamma*C;
	Jy(ii + 1, jj + 1)+=pat->q/dx/dy* wx3 * wy3 * pat->py/gamma*C;


	ii=floor((pat->x - dx / 2) / dx + 1 / 2);
	jj=floor((pat->y - dy / 2) / dy + 1 / 2);
	xx=(pat->x - dx / 2) / dx - ii;
	yy=(pat->y - dy / 2) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	Jz(ii - 1, jj - 1)+=pat->q/dx/dy* wx1 * wy1 * pat->pz/gamma*C;
	Jz(ii - 1, jj)+=pat->q/dx/dy* wx1 * wy2 * pat->pz/gamma*C;
	Jz(ii - 1, jj + 1)+=pat->q/dx/dy* wx1 * wy3 * pat->pz/gamma*C;
	Jz(ii, jj - 1)+=pat->q/dx/dy* wx2 * wy1 * pat->pz/gamma*C;
	Jz(ii, jj)+=pat->q/dx/dy* wx2 * wy2 * pat->pz/gamma*C;
	Jz(ii, jj + 1)+=pat->q/dx/dy* wx2 * wy3 * pat->pz/gamma*C;
	Jz(ii + 1, jj - 1)+=pat->q/dx/dy* wx3 * wy1 * pat->pz/gamma*C;
	Jz(ii + 1, jj)+=pat->q/dx/dy* wx3 * wy2 * pat->pz/gamma*C;
	Jz(ii + 1, jj + 1)+=pat->q/dx/dy* wx3 * wy3 * pat->pz/gamma*C;

	pat=pat->pNext;

}

}