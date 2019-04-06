#include "part_init.h"
#include "fdtd_vmap.h"
#include "part_func.h"
Particle part_check(Particle pt, Grid *g){
double dx, dy;
Particle pH, pM;
pH=pt;

dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;

while(pH){
	if(pH->x<=2*dx||pH->x>(SizeX-2)*dx||pH->y<=2*dy||pH->y>(SizeY-2)*dy) 
	{
		if(pH!=pt){
		pM->pNext=pH->pNext;
		pH=pH->pNext;
		continue;
		}
		else{ pH=pH->pNext; 
		
		      pt=pH;
			  
			  continue;}
	}
	pM=pH;	
	pH=pH->pNext;
}

return pt;

}

Particle check_periodic(Particle pt, Grid *g){
double dx, dy;
int i;
Particle pH;
pH=pt;

//field periodic rearragement
for(i=0;i<SizeX;i++){
	Ez(i,0)=Ez(i,SizeY-2);
	Ez(i,SizeY-1)=Ez(i,1);
}
for(i=0;i<SizeY;i++){
	Ez(0,i)=Ez(SizeX-2,i);
	Ez(SizeX-1,i)=Ez(1,i);
}
for(i=0;i<SizeX-1;i++){
	Ex(i,0)=Ex(i,SizeY-2);
	Ex(i,SizeY-1)=Ex(i,1);
}
for(i=0;i<SizeY-1;i++){
	Ey(0,i)=Ey(SizeX-2,i);
	Ey(SizeX-1,i)=Ey(1,i);
}

//particle periodic rearragement

dx=Lx/SizeX;
dy=Ly/SizeY;

while(pH){
	pH->x=pH->x-floor((pH->x-dx)/(SizeX-2)/dx)*(SizeX-2)*dx;
	pH->y=pH->y-floor((pH->y-dy)/(SizeY-2)/dy)*(SizeY-2)*dy;
	pH=pH->pNext;
}

return pt;

}



void check_Jperiodic(Grid *g){
int i;

//field periodic rearragement
for(i=0;i<SizeX;i++){
	Jz(i, 1)+=Jz(i, SizeY-1);
	Jz(i, SizeY-2)+=Jz(i, 0);
}
for(i=0;i<SizeY;i++){
	Jz(1, i)+=Jz(SizeX-1, i);
	Jz(SizeX-2, i)+=Jz(0, i);
}
for(i=0;i<SizeX-1;i++){
	Jx(i, 1)+=Jx(i, SizeY-1);
	Jx(i, SizeY-2)+=Jx(i, 0);
}
for(i=0;i<SizeY-1;i++){
	Jy(1, i)+=Jy(SizeX-1, i);
	Jy(SizeX-2, i)+=Jy(0, i);
}
for(i=0;i<SizeY;i++){
	Jx(0, i)+=Jx(SizeX-2, i);
	Jx(SizeX-2, i)=Jx(0, i);
}
for(i=0;i<SizeX;i++){
	Jy(i, 0)+=Jy(i, SizeY-2);
	Jy(i, SizeY-2)=Jy(i, 0);
}




}






void check_fperiodic2(Grid *g){
int i;

//field periodic rearragement
if(PeroidY){
for(i=0;i<SizeX;i++){
	Ez(i,0)=Ez((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-4);
	Ez(i,SizeY-1)=Ez((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),3);
	Ez(i,1)=Ez((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-3);
	Ez(i,SizeY-2)=Ez((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),2);
}
for(i=0;i<SizeX-1;i++){
	Ex(i,0)=Ex((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-4);
	Ex(i,SizeY-1)=Ex((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),3);
	Ex(i,1)=Ex((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-3);
	Ex(i,SizeY-2)=Ex((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),2);
}

for(i=0;i<SizeX;i++){
	Ey(i,0)=Ey((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-4);
	Ey(i,SizeY-2)=Ey((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),2);
	Ey(i,1)=Ey((2+(i+SizeX-6)%(SizeX-4))*PeroidX+i*(1-PeroidX),SizeY-3);
}
}
if(PeroidX){

for(i=0;i<SizeY;i++){
	Ez(0,i)=Ez(SizeX-4,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ez(SizeX-1,i)=Ez(3,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ez(1,i)=Ez(SizeX-3,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ez(SizeX-2,i)=Ez(2,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
}
for(i=0;i<SizeY-1;i++){
	Ey(0,i)=Ey(SizeX-4,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ey(SizeX-1,i)=Ey(3,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ey(1,i)=Ey(SizeX-3,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ey(SizeX-2,i)=Ey(2,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
}
for(i=0;i<SizeY;i++){
	Ex(0,i)=Ex(SizeX-4,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ex(SizeX-2,i)=Ex(2,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
	Ex(1,i)=Ex(SizeX-3,(2+(i+SizeY-6)%(SizeY-4))*PeroidY+i*(1-PeroidY));
}

}


}

Particle check_pperiodic2(Particle pt, Grid *g) {
double dx, dy;
Particle pH;
pH=pt;

//particle periodic rearragement

dx=Lx/SizeX*Ul;
dy=Ly/SizeY*Ul;
if(PeroidX||PeroidY){
while(pH){
	if(PeroidX) pH->x=pH->x-floor((pH->x-2*dx)/(SizeX-4)/dx)*(SizeX-4)*dx;
	if(PeroidY) pH->y=pH->y-floor((pH->y-2*dy)/(SizeY-4)/dy)*(SizeY-4)*dy;
	pH=pH->pNext;
}
}
return pt;
}




void check_Jperiodic2(Grid *g){
int i,j;

//field periodic rearragement
if(PeroidY){
for(i=0;i<SizeX;i++){
	Jz(i, 3)+=Jz(i, SizeY-1);
	Jz(i, SizeY-4)+=Jz(i, 0);
	Jz(i, 2)+=Jz(i, SizeY-2);
	Jz(i, SizeY-3)+=Jz(i, 1);
}

for(i=0;i<SizeX-1;i++){
	Jx(i, 3)+=Jx(i, SizeY-1);
	Jx(i, SizeY-4)+=Jx(i, 0);
	Jx(i, 2)+=Jx(i, SizeY-2);
	Jx(i, SizeY-3)+=Jx(i, 1);
}
for(i=0;i<SizeX;i++){
	Jy(i, 2)+=Jy(i, SizeY-2);
	Jy(i, SizeY-4)+=Jy(i, 0);
}

for(i=0;i<SizeX;i++){
	Jy(i, 1)+=Jy(i, SizeY-3);
	Jy(i, SizeY-3)=Jy(i, 1);
}
}

if(PeroidX){

for(i=0;i<SizeY;i++){
	Jz(3, i)+=Jz(SizeX-1, i);
	Jz(SizeX-4, i)+=Jz(0, i);
	Jz(2, i)+=Jz(SizeX-2, i);
	Jz(SizeX-3, i)+=Jz(1, i);
}

for(i=0;i<SizeY-1;i++){
	Jy(3, i)+=Jy(SizeX-1, i);
	Jy(SizeX-4, i)+=Jy(0, i);
	Jy(2, i)+=Jy(SizeX-2, i);
	Jy(SizeX-3, i)+=Jy(1, i);
}
for(i=0;i<SizeY;i++){
	Jx(2, i)+=Jx(SizeX-2, i);
	Jx(SizeX-4, i)+=Jx(0, i);
}
for(i=0;i<SizeY;i++){
	Jx(1, i)+=Jx(SizeX-3, i);
	Jx(SizeX-3, i)=Jx(1, i);
}

}


}