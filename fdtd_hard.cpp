#include "fdtd_vmap.h"
#include "fdtd_phyC.h"
void hardSrc(Grid *g) {	

int kk;
void laser(Grid *g, int mm, int nn);
if(Time>0){
for(kk=0;kk<SizeY;kk++)
{
laser(g, 0, kk);
}
}else return;
/*
double omig;
omig=C*2*M_PI/Lamd;
Hz(SizeX/2, SizeY/2)=10*exp(-(Time-20.0)*(Time-20.0)/10.0/10.0);
*/
return;

}
void laser(Grid *g, int mm, int nn) {
double xfocus, yfocus;
double Omg, K, omg, w0, w, p0, z0;
double phi_Gb, phi_Ge, phi_Reb, phi_Ree, phi_Rbe, phi_Rbb;
double phi_Ey, phi_Ez, phi_By, phi_Bz, phi_0;
double dx, dy, xb, xe, yb, ye;
double E0, t, tp, t0;
dx=Lx/SizeX;
dy=Ly/SizeY;

phi_0 = 0;
xfocus = 8 * Lamd;
yfocus = 0;
w0 = 1.5 * Lamd;
p0 =-1;
E0 =1.0 * 0.5e13;

t0=2*Lamd/C;
tp=3*t0;




K=2*M_PI/Lamd;
omg=K*C;

t=Time*dt;

xb=mm*dx - xfocus;
xe=(mm-0.5)*dx-xfocus;
yb=(nn - SizeY / 2)*dy-yfocus;
ye=(nn - 0.5 - SizeY / 2)*dy-yfocus;

z0=M_PI*w0*w0/Lamd;

phi_Gb=atan(xb/z0);
phi_Ge=atan(xe/z0);

phi_Rbe=0.5*ye*ye*xb/(xb*xb+z0*z0)*K;
phi_Ree=0.5*ye*ye*xe/(xe*xe+z0*z0)*K;
phi_Reb=0.5*yb*yb*xe/(xe*xe+z0*z0)*K;
phi_Rbb=0.5*yb*yb*xb/(xb*xb+z0*z0)*K;

phi_Ey = phi_Ge - phi_Reb + omg*t - K*xe + phi_0;
phi_Ez = phi_Ge - phi_Ree + omg*t - K*xe + phi_0;
phi_By = phi_Gb - phi_Rbe + omg*(t + dt/2) - K*xb + phi_0;
phi_Bz = phi_Gb - phi_Rbb + omg*(t + dt/2) - K*xb + phi_0;





if(p0==0) E0 = E0 * 1.1412;

if(nn==SizeY-1) {
Ez(mm, nn)= E0*(1-p0)/2*sqrt(1/sqrt(1.0+xe*xe/z0/z0))*exp(-ye*ye/w0/w0/(1.0+xe*xe/z0/z0))*sin(phi_Ez)*exp(-(t-tp)*(t-tp)/2/t0/t0);
Hy(mm, nn)= E0*(1-p0)/2*sqrt(1/sqrt(1.0+xb*xb/z0/z0))*exp(-ye*ye/w0/w0/(1.0+xb*xb/z0/z0))*sin(phi_By)*exp(-(t-tp+dt/2)*(t-tp+dt/2)/2/t0/t0)*MU0/C;
}
else {
Ey(mm, nn)= E0*(1+p0)/2*sqrt(1/sqrt(1.0+xe*xe/z0/z0))*exp(-yb*yb/w0/w0/(1.0+xe*xe/z0/z0))*cos(phi_Ey)*exp(-(t-tp)*(t-tp)/2/t0/t0);
Ez(mm, nn)= E0*(1-p0)/2*sqrt(1/sqrt(1.0+xe*xe/z0/z0))*exp(-ye*ye/w0/w0/(1.0+xe*xe/z0/z0))*sin(phi_Ez)*exp(-(t-tp)*(t-tp)/2/t0/t0);
Hy(mm, nn)= E0*(1-p0)/2*sqrt(1/sqrt(1.0+xb*xb/z0/z0))*exp(-ye*ye/w0/w0/(1.0+xb*xb/z0/z0))*sin(phi_By)*exp(-(t-tp+dt/2)*(t-tp+dt/2)/2/t0/t0)*MU0/C;
Hz(mm, nn)=-E0*(1+p0)/2*sqrt(1/sqrt(1.0+xb*xb/z0/z0))*exp(-yb*yb/w0/w0/(1.0+xb*xb/z0/z0))*cos(phi_Bz)*exp(-(t-tp+dt/2)*(t-tp+dt/2)/2/t0/t0)*MU0/C;
}

return;
}
