#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void fdtd_getCJ(Grid *g, Particle pat){


int j1, j2, k1, k2;
double dx, dy, xold, yold, xnew, ynew, xr, yr;
double Fx1, Fy1, Fx2, Fy2;
double xx11, xx12, xx13, xx21, xx22, xx23, yy11, yy12, yy13, yy21, yy22, yy23;
double wx1, wx2, wx3, wx4, wx5, wx6, wy1, wy2, wy3, wy4, wy5, wy6;
int j11, j22, k11, k22;

int ii, jj;
double xx, yy;

double gamma;




dx=Lx/SizeX;
dy=Ly/SizeY;
while(pat){
	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);
	xold=pat->x/dx;
	yold=pat->y/dy;
	xnew=(pat->x + pat->px/gamma*dt*C)/dx;
	ynew=(pat->y + pat->py/gamma*dt*C)/dy;

	j1=floor(xold);
	j2=floor(xnew);
	k1=floor(yold);
	k2=floor(ynew);

	if(j1==j2) xr = (xold + xnew - 1.0)/2.0;
	else xr = (j1 + j2)/2.0;
	if(k1==k2) yr = (yold + ynew - 1.0)/2.0;
	else yr = (k1 + k2)/2.0;

	Fx1 = pat->q*(xr - xold + 0.5)*dx/dt;
	Fy1 = pat->q*(yr - yold + 0.5)*dy/dt;
	Fx2 = pat->q*pat->px/gamma*C - Fx1;
	Fy2 = pat->q*pat->py/gamma*C - Fy1;
	xx11 = (xold - 0.5 + xr)/2.0 - j1 + 1;
	xx12 = fabs((xold - 0.5 + xr)/2.0 - j1);
	xx13 = j1 + 1 - (xold - 0.5 + xr)/2.0;
	xx21 = (xnew - 0.5 + xr)/2.0 - j2 + 1;
	xx22 = fabs((xnew - 0.5 + xr)/2.0 - j2);
	xx23 = j2 + 1 - (xnew - 0.5 + xr)/2.0;
	yy11 = (yold - 0.5 + yr)/2.0 - k1 + 1;
	yy12 = fabs((yold - 0.5 + yr)/2.0 - k1);
	yy13 = k1 + 1 - (yold - 0.5 + yr)/2.0;
	yy21 = (ynew - 0.5 + yr)/2.0 - k2 + 1;
	yy22 = fabs((ynew - 0.5 + yr)/2.0 - k2);
	yy23 = k2 + 1 - (ynew - 0.5 + yr)/2.0;

	wx1 = 0.125*(2*xx11 - 3.0)*(2*xx11 - 3.0);
	wx2 = 0.75 - xx12*xx12;
	wx3 = 0.125*(2*xx13 - 3.0)*(2*xx13 - 3.0);
	wx4 = 0.125*(2*xx21 - 3.0)*(2*xx21 - 3.0);
	wx5 = 0.75 - xx22*xx22;
	wx6 = 0.125*(2*xx23 - 3.0)*(2*xx23 - 3.0);
	wy1 = 0.125*(2*yy11 - 3.0)*(2*yy11 - 3.0);
	wy2 = 0.75 - yy12*yy12;
	wy3 = 0.125*(2*yy13 - 3.0)*(2*yy13 - 3.0);
	wy4 = 0.125*(2*yy21 - 3.0)*(2*yy21 - 3.0);
	wy5 = 0.75 - yy22*yy22;
	wy6 = 0.125*(2*yy23 - 3.0)*(2*yy23 - 3.0);

    j11 = floor(xold - 0.5);
	j22 = floor(xnew - 0.5);
	k11 = floor(yold - 0.5);
	k22 = floor(ynew - 0.5);

	Jx(j11, k1 - 1)+= Fx1*wy1/dx/dy;
	Jx(j22, k2 - 1)+= Fx2*wy4/dx/dy;
	Jx(j11, k1)+= Fx1*wy2/dx/dy;
	Jx(j22, k2)+= Fx2*wy5/dx/dy;
	Jx(j11, k1 + 1)+= Fx1*wy3/dx/dy;
	Jx(j22, k2 + 1)+= Fx2*wy6/dx/dy;

	Jy(j1 - 1, k11)+= Fy1*wx1/dx/dy;
	Jy(j2 - 1, k22)+= Fy2*wx4/dx/dy;
	Jy(j1, k11)+= Fy1*wx2/dx/dy;
	Jy(j2, k22)+= Fy2*wx5/dx/dy;
	Jy(j1 + 1, k11)+= Fy1*wx3/dx/dy;
	Jy(j2 + 1, k22)+= Fy2*wx6/dx/dy;
	




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