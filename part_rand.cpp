#include "fdtd_phyC.h"
double gauss()
{
     static double V1, V2, S;
     static int phase = 0;
     double X;

     if(phase == 0)
     {
        do{
              double U1 = (double)rand() / RAND_MAX;
              double U2 = (double)rand() / RAND_MAX;
               
              V1 = 2 * U1 - 1;
              V2 = 2 * U2 - 1;
              S  = V1 * V1 + V2 * V2;
          }while(S>=1||S==0);
     
          X = V1 * sqrt (-2 * log(S) / S);
     }
     else
     {
          X = V2 * sqrt(-2 * log(S) / S);
     }
     phase = 1 - phase;
     return X;
}
double one()
{
return (double)rand()/RAND_MAX;
}


int getI(double D)
{
    double NewNum=floor(D);
 return ( (D-NewNum)-0.5 <= 0 )? NewNum:NewNum+1; 
}


double besselj0(double x) {

  double ax, z;
  double xx, y, ans, ans1, ans2;
  if((ax=fabs(x))<8.0) {  // calculate directly
    y = x*x;
    ans1 = 57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
    ans2 = 57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
    ans = ans1/ans2;
  }
  else {
    z = 8.0/ax;
    y = z*z;
    xx = ax-0.785398164;
    ans1 = 1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
    ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)));
    ans = sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
  }
  return ans;
}