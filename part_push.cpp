#include "part_init.h"
#include "fdtd_vmap.h"
#include "fdtd_phyC.h"
#include "part_push.h"

Particle part_push(Particle ptc, Grid *g){
double P[9], F[6];
double Dt;

Particle pHd;

pHd=ptc;

Dt=dt*Ut;



while(pHd){

part_getF2(g, pHd, F);

F[0]=F[0]*E*Lamd/ME/C/C/2/M_PI;
F[1]=F[1]*E*Lamd/ME/C/C/2/M_PI;
F[2]=F[2]*E*Lamd/ME/C/C/2/M_PI;
F[3]=MU0*F[3]*E*Lamd/ME/C/2/M_PI;
F[4]=MU0*F[4]*E*Lamd/ME/C/2/M_PI;
F[5]=MU0*F[5]*E*Lamd/ME/C/2/M_PI;



pHd->xold=pHd->x;
pHd->yold=pHd->y;
pHd->zold=pHd->z;


P[0]=pHd->x;
P[1]=pHd->y;
P[2]=pHd->z;
P[3]=pHd->px;
P[4]=pHd->py;
P[5]=pHd->pz;
P[6]=pHd->ax;
P[7]=pHd->ay;
P[8]=pHd->az;




part_Boris(P, F, Dt, pHd->q/pHd->m/E*ME);



pHd->x=P[0];
pHd->y=P[1];
//pHd->z=P[2];
pHd->px=P[3];
pHd->py=P[4];
pHd->pz=P[5];
pHd->ax=P[6];
pHd->ay=P[7];
pHd->az=P[8];

if(pHd->x!=pHd->xold||pHd->y!=pHd->yold||pHd->z!=pHd->zold)
Dt=Dt;


pHd=pHd->pNext;

}

return ptc;

}



void part_Boris(double rpp[9],double emm[6],double del_t, double q_over_m)
{
      double beta,gam_o,gam_n,gam,px_minus,py_minus,pz_minus,px_plus,py_plus,pz_plus;
	  double tx,ty,tz,sx,sy,sz;
	  double vx,vy,vz;
	  beta = 0.5*q_over_m * del_t;

      gam_o=sqrt(1+rpp[3]*rpp[3]+rpp[4]*rpp[4]+rpp[5]*rpp[5]);
	  vx=rpp[3]/gam_o;
	  vy=rpp[4]/gam_o;
	  vz=rpp[5]/gam_o;
	  
	  px_minus = rpp[3] + beta*emm[0];
      py_minus = rpp[4] + beta*emm[1];
      pz_minus = rpp[5] + beta*emm[2];
      gam_n = sqrt(1.0 + px_minus*px_minus + py_minus*py_minus+ pz_minus*pz_minus);
      tx = beta*emm[3]/gam_n;
      ty = beta*emm[4]/gam_n;
      tz = beta*emm[5]/gam_n;
      sx = 2.0*tx/(1.0 + tx*tx+ty*ty+tz*tz);
      sy = 2.0*ty/(1.0 + tx*tx+ty*ty+tz*tz);
      sz = 2.0*tz/(1.0 + tx*tx+ty*ty+tz*tz);
      px_plus = (1.0-sy*ty-sz*tz)*px_minus+       (tx*sy+sz)*py_minus+       (tx*sz-sy)*pz_minus;
      py_plus =        (ty*sx-sz)*px_minus+(1.0-sz*tz-tx*sx)*py_minus+       (ty*sz+sx)*pz_minus;
      pz_plus =        (tz*sx+sy)*px_minus+       (tz*sy-sx)*py_minus+(1.0-tx*sx-sy*ty)*pz_minus;
      rpp[3] = px_plus + beta*emm[0];
      rpp[4] = py_plus + beta*emm[1];
      rpp[5] = pz_plus + beta*emm[2];

      gam = sqrt(1.0+rpp[3]*rpp[3]+rpp[4]*rpp[4]+rpp[5]*rpp[5]);


      rpp[0] = rpp[0] + rpp[3]/gam *del_t;
      rpp[1] = rpp[1] + rpp[4]/gam *del_t;
      rpp[2] = rpp[2] + rpp[5]/gam *del_t;

	  rpp[6] = (rpp[3]/gam - vx)/del_t;
	  rpp[7] = (rpp[4]/gam - vy)/del_t;
	  rpp[8] = (rpp[5]/gam - vz)/del_t;
}
