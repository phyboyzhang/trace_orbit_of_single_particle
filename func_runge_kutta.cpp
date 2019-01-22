#include <math.h>
#include <iostream>
using namespace std;

extern double rad_0,RAD_0,theta_0,psi_0,cof1,cof2,cof3,cof4,q,h,del,momen_p,u_0;


double U(double rad,double theta)
 { double UU=(momen_p-0.5*cof1*pow(rad,2)/q)*sqrt(pow(rad,2)+pow(RAD_0*q,2))
            /(cof2*RAD_0*q*(RAD_0+rad*cos(theta)));
   return(UU);
 }

double momen(double rad,double theta)
 { double mom=cof1*pow(rad,2)*0.5*(RAD_0+rad*cos(theta))/(q*(RAD_0+rad*cos(theta)))
	     +cof2*RAD_0*q*U(rad,theta)*(RAD_0+rad*cos(theta))/sqrt(pow(rad,2)+pow(RAD_0*q,2));
   return mom;
 }

double momen_0(double rad,double theta)
 { double mom=cof1*pow(rad,2)*0.5*(RAD_0+rad*cos(theta))/(q*(RAD_0+rad*cos(theta)))
	     +cof2*RAD_0*q*u_0*(RAD_0+rad*cos(theta))/sqrt(pow(rad,2)+pow(RAD_0*q,2));
   return mom;
 }

double *funrk(double rad, double theta)
 { double funcss[3];
   funcss[0]=-cof1*log((RAD_0+rad*cos(theta))/RAD_0)*0.5*RAD_0*sin(theta)
	        +cof1*0.5*RAD_0*rad*sin(theta)*cos(theta)/(RAD_0+rad*cos(theta));
   funcss[1]=-cof1*log((RAD_0+rad*cos(theta))/RAD_0)*0.5*RAD_0*rad*cos(theta)
	       -cof1*RAD_0*pow(rad,2)*pow(sin(theta),2)*0.5/(RAD_0+rad*cos(theta))
		   +cof2*U(rad,theta)*pow(rad,2)/sqrt(pow(rad,2)+pow(RAD_0*q,2));
   funcss[2]=cof3*sqrt(pow(rad,2)/pow(q,2)+pow(RAD_0,2))/(RAD_0+rad*cos(theta))
	       +cof4*pow(U(rad,theta),2);
   double *p;
   p=funcss;
   return(p);
 }

double *fundiff(double rad,double theta)
 // fun[3][2]保存P_r,p_theta,H的对r,theta的导数。
 { double fun[3][2],fun1[3],fun2[3]; 
   double *p; 
   p=funrk(rad,theta);
   for(int i=0;i<=2;i++)
	   fun1[i]=*(p+i);
   p=funrk(rad+del,theta);
   for(int i=0;i<=2;i++)
	   fun2[i]=*(p+i);
   for(int i=0;i<=2;i++)
	   fun[i][0]=(fun2[i]-fun1[i])/del;
   p=funrk(rad,theta+del);
   for(int i=0;i<=2;i++)
	   fun2[i]=*(p+i);
   for(int i=0;i<=2;i++)
	   fun[i][1]=(fun2[i]-fun1[i])/del;
// k1[2]是右边的驱动项。
   double k1[2];
   k1[0]=fun[2][1]/(fun[0][1]-fun[1][0]);
   k1[1]=fun[2][0]/(fun[1][0]-fun[0][1]);
   p=k1;
   return(p);
 }

double *fun_2(double rad,double theta)
 { //double k2[2];
   double *p;
   p=fundiff(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=fundiff(rad+0.5*h*f[0],theta+0.5*h*f[1]);
   return(p);
 }

double *fun_3(double rad,double theta)
 { double *p;
   p=fun_2(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=fundiff(rad+0.5*h*f[0],rad+0.5*h*f[1]);
   return(p);
 }

double *fun_4(double rad,double theta)
 { double *p;
   p=fun_3(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=fundiff(rad+h*f[0],rad+h*f[1]);
   return(p);
 }	

double *solv(double rad, double theta)
 { double *p;
   double f1[2],f2[2],f3[2],f4[2];
   p=fundiff(rad,theta);
   f1[0]=*p,f1[1]=*(p+1);
   p=fun_2(rad,theta);
   f2[0]=*p,f2[1]=*(p+1);
   p=fun_3(rad,theta);
   f3[0]=*p,f3[1]=*(p+1);
   p=fun_4(rad,theta);
   f4[0]=*p,f4[1]=*(p+1);
   double rad1=rad,theta1=theta;
   rad=rad1+(h/6.0)*(f1[0]+2.0*f2[0]+2.0*f3[0]+f4[0]);
   theta=theta1+(h/6.0)*(f1[1]+2.0*f2[1]+2.0*f3[1]+f4[1]);
   double f5[2];
   f5[0]=rad;
   f5[1]=theta;
   p=f5;
   return(p);
 }