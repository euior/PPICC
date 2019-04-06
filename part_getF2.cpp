#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void part_getF2(Grid *g, Particle pat, double F[6]){
int ii, jj;
double dx, dy, xx, yy;
double wx1, wx2, wx3, wy1, wy2, wy3;
dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;

	ii=floor((pat->x - dx) / dx + 0.5);
	jj=floor((pat->y - dy / 2) / dy + 0.5);
	xx=(pat->x - dx) / dx - ii;
	yy=(pat->y - dy / 2) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	F[0]= Ex(ii - 1,jj - 1)* wx1 * wy1 + Ex(ii - 1,jj) * wx1 * wy2 + Ex(ii - 1,jj+1) * wx1 * wy3 + \
		  Ex(ii, jj - 1) * wx2 * wy1 + Ex(ii,jj) * wx2 * wy2 + Ex(ii,jj + 1) * wx2 * wy3 + \
		  Ex(ii + 1,jj - 1) * wx3 * wy1 + Ex(ii + 1,jj) * wx3 * wy2 + Ex(ii + 1, jj + 1) * wx3 * wy3;
	F[4]= (Hy(ii - 1,jj - 1) + Hyold(ii - 1,jj - 1)) / 2 * wx1 * wy1 + (Hy(ii - 1,jj) + Hyold(ii - 1,jj)) / 2 * wx1 * wy2 + (Hy(ii - 1,jj + 1) + Hyold(ii - 1,jj + 1)) / 2 * wx1 * wy3 + \
		  (Hy(ii, jj - 1) + Hyold(ii, jj - 1)) / 2 * wx2 * wy1 + (Hy(ii,jj) + Hyold(ii,jj)) / 2 * wx2 * wy2 + (Hy(ii,jj + 1) + Hyold(ii,jj + 1)) / 2 * wx2 * wy3 + \
		  (Hy(ii + 1,jj - 1) + Hyold(ii + 1,jj - 1)) / 2 * wx3 * wy1 + (Hy(ii + 1,jj) + Hyold(ii + 1,jj)) / 2 * wx3 * wy2 + (Hy(ii + 1, jj + 1) + Hyold(ii + 1, jj + 1)) / 2 * wx3 * wy3;    
	
	
	ii=floor((pat->x - dx / 2) / dx + 0.5);
	jj=floor((pat->y - dy) / dy + 0.5);
	xx=(pat->x - dx / 2) / dx - ii;
	yy=(pat->y - dy) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	F[1]= Ey(ii - 1,jj - 1)* wx1 * wy1 + Ey(ii - 1,jj) * wx1 * wy2 + Ey(ii - 1,jj+1) * wx1 * wy3 + \
		  Ey(ii, jj - 1) * wx2 * wy1 + Ey(ii,jj) * wx2 * wy2 + Ey(ii,jj + 1) * wx2 * wy3 + \
		  Ey(ii + 1,jj - 1) * wx3 * wy1 + Ey(ii + 1,jj) * wx3 * wy2 + Ey(ii + 1, jj + 1) * wx3 * wy3;
	F[3]= (Hx(ii - 1,jj - 1) + Hxold(ii - 1,jj - 1)) / 2 * wx1 * wy1 + (Hx(ii - 1,jj) + Hxold(ii - 1,jj)) / 2 * wx1 * wy2 + (Hx(ii - 1,jj + 1) + Hxold(ii - 1,jj + 1)) / 2 * wx1 * wy3 + \
		  (Hx(ii, jj - 1) + Hxold(ii, jj - 1)) / 2 * wx2 * wy1 + (Hx(ii,jj) + Hxold(ii,jj)) / 2 * wx2 * wy2 + (Hx(ii,jj + 1) + Hxold(ii,jj + 1)) / 2 * wx2 * wy3 + \
		  (Hx(ii + 1,jj - 1) + Hxold(ii + 1,jj - 1)) / 2 * wx3 * wy1 + (Hx(ii + 1,jj) + Hxold(ii + 1,jj)) / 2 * wx3 * wy2 + (Hx(ii + 1, jj + 1) + Hxold(ii + 1, jj + 1)) / 2 * wx3 * wy3;

	ii=floor((pat->x - dx / 2) / dx + 0.5);
	jj=floor((pat->y - dy / 2) / dy + 0.5);
	xx=(pat->x - dx / 2) / dx - ii;
	yy=(pat->y - dy / 2) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	F[2]= Ez(ii - 1,jj - 1)* wx1 * wy1 + Ez(ii - 1,jj) * wx1 * wy2 + Ez(ii - 1,jj+1) * wx1 * wy3 + \
		  Ez(ii, jj - 1) * wx2 * wy1 + Ez(ii,jj) * wx2 * wy2 + Ez(ii,jj + 1) * wx2 * wy3 + \
		  Ez(ii + 1,jj - 1) * wx3 * wy1 + Ez(ii + 1,jj) * wx3 * wy2 + Ez(ii + 1, jj + 1) * wx3 * wy3;

	
	ii=floor((pat->x - dx) / dx + 0.5);
	jj=floor((pat->y - dy) / dy + 0.5);
	xx=(pat->x - dx) / dx - ii;
	yy=(pat->y - dy) / dy - jj;
	
	wx1 = 0.5 * (xx - 0.5) * (xx - 0.5);
	wx2 = 0.75 - xx * xx;
	wx3 = 1 - wx1 - wx2;
	wy1 = 0.5 * (yy - 0.5) * (yy - 0.5);
	wy2 = 0.75 - yy * yy;
	wy3 = 1 - wy1 - wy2;

	F[5]= (Hz(ii - 1,jj - 1) + Hzold(ii - 1,jj - 1)) / 2 * wx1 * wy1 + (Hz(ii - 1,jj) + Hzold(ii - 1,jj)) / 2 * wx1 * wy2 + (Hz(ii - 1,jj + 1) + Hzold(ii - 1,jj + 1)) / 2 * wx1 * wy3 + \
		  (Hz(ii, jj - 1) + Hzold(ii, jj - 1)) / 2 * wx2 * wy1 + (Hz(ii,jj) + Hzold(ii,jj)) / 2 * wx2 * wy2 + (Hz(ii,jj + 1) + Hzold(ii,jj + 1)) / 2 * wx2 * wy3 + \
		  (Hz(ii + 1,jj - 1) + Hzold(ii + 1,jj - 1)) / 2 * wx3 * wy1 + (Hz(ii + 1,jj) + Hzold(ii + 1,jj)) / 2 * wx3 * wy2 + (Hz(ii + 1, jj + 1) + Hzold(ii + 1, jj + 1)) / 2 * wx3 * wy3;


}