#include "fdtd_phyC.h"
#include "part_func.h"
Particle part_disp(Grid *g){
int ii, jj, ss;
int mm, nn, ptca;
int id;
double dx, dy, density;
srand(time(NULL));
dx=Lx/SizeX;
dy=Ly/SizeY;
Particle part=(Particle)malloc(sizeof(particle));
part=NULL;

id=0;
for(ss=0;ss<nSpecies;ss++){
	ptca=getI(sqrt(PTC1));
	
	for(ii=cpmlxl+2;ii<SizeX-cpmlxr-2;ii++){
		for(jj=cpmlyd+2;jj<SizeY-cpmlyu-2;jj++){
			for(mm=0; mm<ptca; mm++){
				for(nn=0;nn<ptca;nn++){
					density=part_density(g, ss, ii*dx, jj*dy);
					if(density>EPSL0){
				Particle pNew=(Particle)malloc(sizeof(particle));
				if(SIG){
					pNew->x=(ii * dx + dx / ptca * (mm + 0.5))*Ul;
					pNew->y=(jj * dy + dy / ptca * (nn + 0.5))*Ul;
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
					pNew->m=M1*density*dx*dy/ptca/ptca*(1-ss)+1836.0*M1*density*dx*dy/ptca/ptca*ss;
					pNew->q=Q1*density*dx*dy/ptca/ptca*(0.5-ss)*2.0;
					pNew->px=gauss()*(1-ss)*0.1398*0;
					pNew->py=gauss()*(1-ss)*0.1398*0;
					pNew->pz=gauss()*(1-ss)*0.1398*0;				
					pNew->ID=ss;
					part=part_add(part, pNew);
					id++;
					}

				}
			}

		}

	}

}
return part;
}