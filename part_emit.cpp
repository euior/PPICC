#include "fdtd_phyC.h"
#include "part_func.h"
Particle part_emit(Particle pp, Grid *g){
int ii, jj, ss;
int mm, nn, ptca;
int id;
double dx, dy, flux, v0;
Particle part;

srand(time(NULL));
dx=Lx/SizeX;
dy=Ly/SizeY;

part = pp;

v0=1.5;
for(ss=0;ss<1;ss++){
	ptca=getI(sqrt(PTC1));
	
	for(ii=cpmlxl+2;ii<SizeX-cpmlxr-2;ii++){
		for(jj=42;jj<43;jj++){
			for(mm=0; mm<ptca; mm++){
				for(nn=0;nn<ptca;nn++){
					flux=emit_flux(g, ss, v0, ii*dx, Time * dt);
					if(flux>EPSL0){
				Particle pNew=(Particle)malloc(sizeof(particle));
				if(SIG){
					pNew->x=(ii * dx + dx / ptca * (mm + 0.5))*Ul;
					pNew->y=(jj * dy + dy / ptca * (nn * v0 / sqrt(v0 * v0 + 1) * C * dt / dy + 0.5))*Ul;
					pNew->z=0.0*Ul;
					pNew->xold=pNew->x;
					pNew->yold=pNew->y;
					pNew->zold=pNew->z;
				}else{
					pNew->x=(ii * dx + dx*one())*Ul;
					pNew->y=(jj * dy + dy*one())*Ul;
					pNew->z=0.0*Ul;
					pNew->xold=pNew->x;
					pNew->yold=pNew->y;
					pNew->zold=pNew->z;
				}
					pNew->m=M1*flux/ptca/ptca;
					pNew->q=Q1*flux/ptca/ptca;
					pNew->px=gauss()*0.625/15.0*1.0;
					pNew->py=gauss()*0.625 + v0;
					pNew->pz=gauss()*0.625/15.0*1.0;
					pNew->ID=8899166;
					part=part_add(part, pNew);
					}

				}
			}

		}

	}

}
return part;
}

double emit_flux(Grid *g, int s, double vo, double x, double t){

double sigmax, perd, nc;
double t0, t1, x0;
nc = 1.1e27 / (Lamd / 1e-6) / (Lamd / 1e-6);
perd=Lamd/C;
sigmax=10.6*Lamd;
t1=perd*31;

t0=perd;
x0=10.0 * Lamd;

if(t<=t0) return 0.01 * nc * exp(-(x-x0)*(x-x0)/2/sigmax/sigmax) * (t / perd) * Lx / SizeX * dt * vo / sqrt(vo*vo + 1)*C;
else if(t<=(t0+t1)) return 0.01 * nc * exp(-(x-x0)*(x-x0)/2/sigmax/sigmax) * Lx / SizeX * dt * vo / sqrt(vo*vo + 1) * C;
}



