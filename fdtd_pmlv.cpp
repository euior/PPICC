#include "fdtd_vmap.h"
#include "fdtd_phyC.h"

void CPML_coeffs(Grid *g){

void calc_CPML_coeffs(double depth, double PML_depth, double sigmamax, double kappamax, double alphamax, double m, double ma, double *ptr);

int mm, ncpml[4];
double sigmx, sigmy;
double alpmx;
double alpmy;
double kapmx=5.0;//1.0
double kapmy=5.0;//1.0
double dx, dy, dep, dep_pml;
double kbc[4];
double mox, moax;
double moy, moay;
dx=Lx/SizeX;
dy=Ly/SizeY;
ncpml[0]=cpmlxl;
ncpml[1]=cpmlxr;
ncpml[2]=cpmlyd;
ncpml[3]=cpmlyu;

kbc[3]=dt;

mox=3.0;
moy=3.0;
moax=1.0;
moay=1.0;

alpmx=0.05*1.0; //2*M_PI*EP0*C/Lamd/10;
alpmy=0.05*1.0; //2*M_PI*EP0*C/Lamd/10;

sigmx= (mox + 1) * 0.8 / 377 / dx;//EP0/2/dt;//
sigmy= (moy + 1) * 0.8 / 377 / dy;//EP0/2/dt;//

for(mm=0;mm<SizeY;mm++){

Phyz_k(mm) = 1.0;
if(mm<ncpml[2]){

	dep_pml = ncpml[2] * dy;
	dep = dep_pml - (mm + 0.5) * dy;
	calc_CPML_coeffs(dep,dep_pml,sigmy,kapmy,alpmy,moy,moay,kbc);
    Phyz_k(mm)=kbc[0];
	Phyz_b(mm)=kbc[1];
	Phyz_c(mm)=kbc[2];
}
if(mm>SizeY-1-ncpml[3]){

	dep_pml = ncpml[3] * dy;
	dep = dep_pml - (SizeY - mm - 0.5) * dy;
	calc_CPML_coeffs(dep,dep_pml,sigmy,kapmy,alpmy,moy,moay,kbc);
    Phyz_k(mm)=kbc[0];
	Phyz_b(mm)=kbc[1];
	Phyz_c(mm)=kbc[2];
}

}

for(mm=0;mm<SizeX;mm++){

Phxz_k(mm) = 1.0;
if(mm<ncpml[0]){

	dep_pml = ncpml[0] * dx;
	dep = dep_pml - (mm + 0.5) * dx;
	calc_CPML_coeffs(dep,dep_pml,sigmx,kapmx,alpmx,mox,moax,kbc);
    Phxz_k(mm)=kbc[0];
	Phxz_b(mm)=kbc[1];
	Phxz_c(mm)=kbc[2];
}
if(mm>SizeX-1-ncpml[1]){

	dep_pml = ncpml[1] * dx;
	dep = dep_pml - (SizeX - mm - 0.5) * dx;
	calc_CPML_coeffs(dep,dep_pml,sigmx,kapmx,alpmx,mox,moax,kbc);
    Phxz_k(mm)=kbc[0];
	Phxz_b(mm)=kbc[1];
	Phxz_c(mm)=kbc[2];
}

}
for(mm=0;mm<SizeY-1;mm++){

Peyz_k(mm) = 1.0;
if(mm<ncpml[2]){

	dep_pml = ncpml[2] * dy;
	dep = dep_pml - (mm + 1.0) * dy;
	calc_CPML_coeffs(dep,dep_pml,sigmy,kapmy,alpmy,moy,moay,kbc);
    Peyz_k(mm)=kbc[0];
    Peyz_b(mm)=kbc[1];
    Peyz_c(mm)=kbc[2];

}
if(mm>SizeY-2-ncpml[3]){

	dep_pml = ncpml[3] * dy;
	dep = dep_pml - (SizeY - mm - 1.0) * dy;
	calc_CPML_coeffs(dep,dep_pml,sigmy,kapmy,alpmy,moy,moay,kbc);
    Peyz_k(mm)=kbc[0];
    Peyz_b(mm)=kbc[1];
    Peyz_c(mm)=kbc[2];


}

}
for(mm=0;mm<SizeX-1;mm++){

Pexz_k(mm) = 1.0;
if(mm<ncpml[0]){

	dep_pml = ncpml[0] * dx;
	dep = dep_pml - (mm + 1.0) * dx;
	calc_CPML_coeffs(dep,dep_pml,sigmx,kapmx,alpmx,mox,moax,kbc);
    Pexz_k(mm)=kbc[0];
    Pexz_b(mm)=kbc[1];
    Pexz_c(mm)=kbc[2];

}
if(mm>SizeX-2-ncpml[1]){

	dep_pml = ncpml[1] * dx;
	dep = dep_pml - (SizeX - mm - 1.0) * dx;
	calc_CPML_coeffs(dep,dep_pml,sigmx,kapmx,alpmx,mox,moax,kbc);
    Pexz_k(mm)=kbc[0];
    Pexz_b(mm)=kbc[1];
    Pexz_c(mm)=kbc[2];
}

}

return;
}



void calc_CPML_coeffs(double depth, double PML_depth, double sigmamax, double kappamax, double alphamax, double m, double ma, double *ptr)
{
    double sigma = sigmamax*pow((depth/PML_depth), m);
    double alpha = alphamax*pow(((PML_depth-depth)/PML_depth), ma);
    ptr[0] = 1.0 + (kappamax - 1.0)*pow((depth/PML_depth), m);
    ptr[1] = exp( -(sigma/ptr[0] + alpha)*ptr[3]/EP0 );
    ptr[2] = sigma*(ptr[1] - 1.0)/(sigma*ptr[0] + ptr[0]*ptr[0]*alpha);
	return;
}
