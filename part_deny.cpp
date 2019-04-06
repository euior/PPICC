#include "fdtd_vmap.h"
double part_density(Grid *g, int sp, double x, double y){
	double dx, dy, nc;
	nc=1.1e27 / (Lamd / 1e-6) / (Lamd / 1e-6);
	dx=Lx/SizeX;
	dy=Ly/SizeY;

if( (y>30*dy) && (y<=60*dy) ) return  1.0 * 0.1 * nc * (y - 30 * dy) / 30 / dy;
else if(y>=60*dy&&y<=380*dy) return 1.0 * 0.1 * nc;
else return 0.0;

/*
if(x>2*dx&&x<(SizeX-2)*dx){	
if( (y>2*dy) && (y<(SizeY-2)*dy) )
{
	if(sp==0) return 0.1*nc;
	else return 0.9*nc;
}   else return 0.0;
}   else return 0.0;
*/
//return 0.0;

}