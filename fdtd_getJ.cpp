#include "fdtd_vmap.h"
#include "part_init.h"
#include "fdtd_phyC.h"
void fdtd_getJ(Grid *g, Particle pat){
int ii, jj;
double dx, dy, xx, yy;
double gamma;
dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;
while(pat){
	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);
	ii=floor((pat->x - dx) / dx);
	jj=floor((pat->y - dy / 2) / dy);
	xx=pat->x - dx - ii*dx;
	yy=pat->y - dy / 2 - jj*dy;
	
	Jx(ii, jj)+=pat->q/dx/dy*(dx - xx) * (dy - yy) / dx / dy * pat->px/gamma*C;
    Jx(ii+1, jj)+=pat->q/dx/dy* xx * (dy - yy) / dx / dy * pat->px/gamma*C;
	Jx(ii, jj+1)+=pat->q/dx/dy* yy * (dx - xx) / dx / dy * pat->px/gamma*C;
    Jx(ii+1, jj+1)+=pat->q/dx/dy* yy * xx / dx / dy * pat->px/gamma*C;

	ii=floor((pat->x - dx / 2) / dx);
	jj=floor((pat->y - dy) / dy);
	xx=pat->x - dx / 2 - ii*dx;
	yy=pat->y - dy - jj*dy;
	
	Jy(ii, jj)+=pat->q/dx/dy*(dx - xx) * (dy - yy) / dx / dy * pat->py/gamma*C;
    Jy(ii+1, jj)+=pat->q/dx/dy* xx * (dy - yy) / dx / dy * pat->py/gamma*C;
	Jy(ii, jj+1)+=pat->q/dx/dy* yy * (dx - xx) / dx / dy * pat->py/gamma*C;
    Jy(ii+1, jj+1)+=pat->q/dx/dy* yy * xx / dx / dy * pat->py/gamma*C;



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


void fdtd_getJ2001(Grid *g, Particle pat){
int i, j, ii, jj, iii, jjj;
double s0[5][2], s1[5][2], w[5][5][3];
double J[4];
double dx, dy, ddx, ddy, xx, yy;
double gamma;

ddx=Lx/SizeX;
ddy=Ly/SizeY;

dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;


while(pat){

	for(ii=0;ii<5;ii++){
	s1[ii][0]=0;
	s1[ii][1]=0;
	}

	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);

	ii=floor((pat->xold - dx / 2) / dx + 0.5);
	jj=floor((pat->yold - dy / 2) / dy + 0.5);
	xx=(pat->xold - dx / 2) / dx - ii;
	yy=(pat->yold - dy / 2) / dy - jj;

	i = ii;
	j = jj;


	s0[0][0] = 0.0;
	s0[1][0] = 0.5 * (xx - 0.5) * (xx - 0.5);
	s0[2][0] = 0.75 - xx * xx;
	s0[3][0] = 1 - s0[1][0] - s0[2][0];
	s0[4][0] = 0.0;

	s0[0][1] = 0.0;
	s0[1][1] = 0.5 * (yy - 0.5) * (yy - 0.5);
	s0[2][1] = 0.75 - yy * yy;
	s0[3][1] = 1 - s0[1][1] - s0[2][1];
	s0[4][1] = 0.0;

	iii=floor((pat->x - dx / 2) / dx + 0.5);
	jjj=floor((pat->y - dy / 2) / dy + 0.5);
	xx=(pat->x - dx / 2) / dx - iii;
	yy=(pat->y - dy / 2) / dy - jjj;

	iii = floor(iii*1.0 - ii + 0.5);
	jjj = floor(jjj*1.0 - jj + 0.5);

	s1[1+iii][0] = 0.5 * (xx - 0.5) * (xx - 0.5);
	s1[2+iii][0] = 0.75 - xx * xx;
	s1[3+iii][0] = 1 - s1[1+iii][0] - s1[2+iii][0];

	s1[1+jjj][1] = 0.5 * (yy - 0.5) * (yy - 0.5);
	s1[2+jjj][1] = 0.75 - yy * yy;
	s1[3+jjj][1] = 1 - s1[1+jjj][1] - s1[2+jjj][1];

	for(ii = 0;ii < 5;ii++){
	s1[ii][0] = s1[ii][0] - s0[ii][0];
	s1[ii][1] = s1[ii][1] - s0[ii][1];
	}

	for(ii = 0;ii < 5;ii++){
		for(jj = 0;jj < 5;jj++){
		w[ii][jj][0] = s1[ii][0] * (s0[jj][1] + 0.5 * s1[jj][1]);
		w[ii][jj][1] = s1[jj][1] * (s0[ii][0] + 0.5 * s1[ii][0]);
		w[ii][jj][2] = s0[ii][0] * s0[jj][1] + 0.5 * s1[ii][0] * s0[jj][1] + 0.5 * s0[ii][0] * s1[jj][1] + 1.0 / 3 * s1[ii][0] * s1[jj][1];		
		}
	}

	for(jj =-2;jj <= 2;jj++){
		J[0] = - pat->q / ddx / ddy * ddx / dt * w[0][jj+2][0];
		J[1] = J[0] - pat->q / ddx / ddy * ddx / dt * w[1][jj+2][0];
		J[2] = J[1] - pat->q / ddx / ddy * ddx / dt * w[2][jj+2][0];
		J[3] = J[2] - pat->q / ddx / ddy * ddx / dt * w[3][jj+2][0];

		Jx(i-2, j + jj)+= J[0];
		Jx(i-1, j + jj)+= J[1];
		Jx(i+0, j + jj)+= J[2];
		Jx(i+1, j + jj)+= J[3];

		J[0] = - pat->q / ddx / ddy * ddy / dt * w[jj+2][0][1];
		J[1] = J[0] - pat->q / ddx / ddy * ddy / dt * w[jj+2][1][1];
		J[2] = J[1] - pat->q / ddx / ddy * ddy / dt * w[jj+2][2][1];
		J[3] = J[2] - pat->q / ddx / ddy * ddy / dt * w[jj+2][3][1];

		Jy(i + jj, j-2)+= J[0];
		Jy(i + jj, j-1)+= J[1];
		Jy(i + jj, j+0)+= J[2];
		Jy(i + jj, j+1)+= J[3];

		for(ii =-2;ii<=2;ii++){
			Jz(i + ii, j + jj)+= pat->q / ddx / ddy * pat->pz / gamma * C * w[ii+2][jj+2][2];		
		}
	}
	pat=pat->pNext;

}

}



void fdtd_getD(Grid *g, Particle pat){
int ii, jj;
double dx, dy, xx, yy;
double wx1, wx2, wx3, wy1, wy2, wy3;
double gamma;
dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;
while(pat){
	if(pat->ID==8899166){
	gamma=sqrt(1.0 + pat->px*pat->px + pat->py*pat->py + pat->pz*pat->pz);


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

	Density(ii - 1, jj - 1)+= pat->m / ME * wx1 * wy1;
	Density(ii - 1, jj)+= pat->m / ME * wx1 * wy2;
	Density(ii - 1,jj + 1)+= pat->m / ME * wx1 * wy3;
	Density(ii, jj - 1)+= pat->m / ME * wx2 * wy1;
	Density(ii, jj)+= pat->m / ME * wx2 * wy2;
	Density(ii,jj + 1)+= pat->m / ME * wx2 * wy3;
	Density(ii + 1,jj - 1)+= pat->m / ME * wx3 * wy1;
	Density(ii + 1,jj)+= pat->m / ME * wx3 * wy2;
	Density(ii + 1, jj + 1)+= pat->m / ME * wx3 * wy3;
	}

	pat=pat->pNext;

}

}