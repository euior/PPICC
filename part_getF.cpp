#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void part_getF(Grid *g, Particle pat, double F[6]){
int ii, jj;
double dx, dy, xx, yy;
dx=Lx/SizeX;
dy=Ly/SizeY;

	ii=floor((pat->x - dx) / dx);
	jj=floor((pat->y - dy / 2) / dy);
	xx=pat->x - dx - ii*dx;
	yy=pat->y - dy / 2 - jj*dy;

	F[0]=Ex(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Ex(ii+1,jj) * xx * (dy - yy)/dx/dy + Ex(ii,jj+1) * (dx - xx) * yy /dx/dy + Ex(ii+1, jj+1) * xx * yy /dx/dy;
	F[4]=Hy(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Hy(ii+1,jj) * xx * (dy - yy)/dx/dy + Hy(ii,jj+1) * (dx - xx) * yy /dx/dy + Hy(ii+1, jj+1) * xx * yy /dx/dy;
    
	
	
	ii=floor((pat->x - dx / 2) / dx);
	jj=floor((pat->y - dy) / dy);
	xx=pat->x - dx / 2 - ii*dx;
	yy=pat->y - dy - jj*dy;
	F[1]=Ey(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Ey(ii+1,jj) * xx * (dy - yy)/dx/dy + Ey(ii,jj+1) * (dx - xx) * yy /dx/dy + Ey(ii+1, jj+1) * xx * yy /dx/dy;
	F[3]=Hx(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Hx(ii+1,jj) * xx * (dy - yy)/dx/dy + Hx(ii,jj+1) * (dx - xx) * yy /dx/dy + Hx(ii+1, jj+1) * xx * yy /dx/dy;

	ii=floor((pat->x - dx / 2) / dx);
	jj=floor((pat->y - dy / 2) / dy);
	xx=pat->x - dx / 2 - ii*dx;
	yy=pat->y - dy / 2 - jj*dy;
	F[2]=Ez(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Ez(ii+1,jj) * xx * (dy - yy)/dx/dy + Ez(ii,jj+1) * (dx - xx) * yy /dx/dy + Ez(ii+1, jj+1) * xx * yy /dx/dy;
	
	ii=floor((pat->x - dx) / dx);
	jj=floor((pat->y - dy) / dy);
	xx=pat->x - dx - ii*dx;
	yy=pat->y - dy - jj*dy;
	F[5]=Hz(ii,jj)*(dx - xx) * (dy - yy) / dx / dy + Hz(ii+1,jj) * xx * (dy - yy)/dx/dy + Hz(ii,jj+1) * (dx - xx) * yy /dx/dy + Hz(ii+1, jj+1) * xx * yy /dx/dy;


}