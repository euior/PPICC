#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void fdtd_getCJ1(Grid *g, Particle pat){


int i1, i2, j1, j2;
double dx, dy, x1, y1, x2, y2;
double Fx1, Fy1, Fx2, Fy2;
double wip1, wjp1;
double wim2, wio2, wip2, wjm2, wjo2, wjp2;
int ii, jj;
double xx, yy;

double gamma;




dx=Lx/SizeX;
dy=Ly/SizeY;
while(pat){
	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);
	x1=pat->x/dx - 0.5;
	y1=pat->y/dy - 0.5;
	x2=(pat->x + pat->px/gamma*dt*C)/dx - 0.5;
	y2=(pat->y + pat->py/gamma*dt*C)/dy - 0.5;

	i1=floor(x1 + 0.5);
	i2=floor(x2 + 0.5);
	j1=floor(y1 + 0.5);
	j2=floor(y2 + 0.5);

 	wip1=(x1 + x2)/2.0 - i1;
    wjp1=(y1 + y2)/2.0 - j1;

	Fx1 = pat->q*pat->px/gamma*C*(0.5 - wip1);
	Fy1 = pat->q*pat->py/gamma*C*(0.5 - wjp1);
	Fx2 = pat->q*pat->px/gamma*C*(0.5 + wip1);
	Fy2 = pat->q*pat->py/gamma*C*(0.5 + wjp1);
	
	wim2 = 0.5*(0.5 - wip1)*(0.5 - wip1);
	wio2 = 0.75 - wip1*wip1;
	wip2 = 0.5*(0.5 + wip1)*(0.5 + wip1);

	wjm2 = 0.5*(0.5 - wjp1)*(0.5 - wjp1);
	wjo2 = 0.75 - wjp1*wjp1;
	wjp2 = 0.5*(0.5 + wjp1)*(0.5 + wjp1);

	Jx(i1 - 1, j1 - 1)+= Fx1*wjm2/dx/dy;
	Jx(i1 - 1, j1)+= Fx1*wjo2/dx/dy;
	Jx(i1 - 1, j1 + 1)+= Fx1*wjp2/dx/dy;

	Jx(i1, j1 - 1)+= Fx2*wjm2/dx/dy;
	Jx(i1, j1)+= Fx2*wjo2/dx/dy;
	Jx(i1, j1 + 1)+= Fx2*wjp2/dx/dy;

	Jy(i1 - 1, j1 - 1)+= Fy1*wim2/dx/dy;
	Jy(i1, j1 - 1)+= Fy1*wio2/dx/dy;
	Jy(i1 + 1, j1 - 1)+= Fy1*wip2/dx/dy;

	Jy(i1 - 1, j1)+= Fy2*wim2/dx/dy;
	Jy(i1, j1)+= Fy2*wio2/dx/dy;
	Jy(i1 + 1, j1)+= Fy2*wip2/dx/dy;



	ii=floor((pat->x - dx / 2) / dx);
	jj=floor((pat->y - dy / 2) / dy);
	xx=pat->x - dx / 2 - ii*dx;
	yy=pat->y - dy / 2 - jj*dy;
	
	Jz(ii, jj)+=pat->q/dx/dy*(dx - xx) * (dy - yy) / dx / dy * pat->pz/gamma*C;
    Jz(ii+1, jj)+=pat->q/dx/dy* xx * (dy - yy) / dx / dy * pat->pz/gamma*C;
	Jz(ii, jj+1)+=pat->q/dx/dy* yy * (dx - xx) / dx / dy * pat->pz/gamma*C;
    Jz(ii+1, jj+1)+=pat->q/dx/dy* yy * xx / dx / dy * pat->pz/gamma*C;

	pat=pat->pNext;

}

}