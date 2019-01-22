#include <math.h>
#include "fun_class.h"

extern double rad_0,RAD_0,theta_0,psi_0,cof1,cof2,cof3,cof4,q,h,del,momen_p,u_0;

double *fun2_2(double rad,double theta£©
 { double *p;
   p=dif_fun(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=dif_fun(rad+h*f[0],theta+h*f[1]);
   return p;
}

double *solv_2(double rad, double theta)
 { double *p;
   double f1[2],f2[2];
   p=dif_fun(rad,theta);
   f1[0]=*p,f1[1]=*(p+1);
   p=fun2_2(rad,theta);
   f2[0]=*p,f2[1]=*(p+1);
   double rad1=rad,theta1=theta;
   rad=rad1+(h/2.0)*(f1[0]+f2[0]);
   theta=theta1+(h/2.0)*(f1[1]+f2[1]);
   double f5[2];
   f5[0]=rad;
   f5[1]=theta;
   p=f5;
   return(p);
 }
