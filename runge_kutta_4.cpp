#include <math.h>

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


double *dif_fun(double rad,double theta)
 { double func[3][2];
   func[0][0]=-cof1*0.5*RAD_0*sin(theta)*(RAD_0/(RAD_0+rad*cos(theta)))*cos(theta)/RAD_0
	         +cof1*RAD_0*sin(theta)*cos(theta)*0.5/(RAD_0+rad*cos(theta))
			 +cof1*RAD_0*rad*sin(theta)*cos(theta)*0.5*(-cos(theta))/pow((RAD_0+rad*cos(theta)),2);
   
   func[0][1]=-cof1*RAD_0*0.5*log((RAD_0+rad*cos(theta))/(RAD_0))*cos(theta)
	         -cof1*RAD_0*sin(theta)*0.5*((RAD_0)/(RAD_0+rad*cos(theta)))*(-rad)*sin(theta)/RAD_0
             +cof1*RAD_0*0.5*rad*pow(cos(theta),2)/(RAD_0+rad*cos(theta))
			 +cof1*RAD_0*rad*(-1.0)*pow(sin(theta),2)*0.5/(RAD_0+rad*sin(theta))
			 +cof1*RAD_0*rad*sin(theta)*cos(theta)*0.5*rad*sin(theta)/pow(RAD_0+rad*cos(theta),2);

   func[1][0]=-cof1*log((RAD_0+rad*cos(theta))/(RAD_0))*0.5*RAD_0*cos(theta)
	         -cof1*RAD_0*rad*cos(theta)*0.5*((RAD_0)/(RAD_0+rad*cos(theta)))*cos(theta)/RAD_0
			 -cof1*RAD_0*pow(sin(theta),2)*rad/(RAD_0+rad*cos(theta))
			 -cof1*RAD_0*pow(rad*sin(theta),2)*0.5*(-cos(theta))/pow(RAD_0+rad*cos(theta),2)
			 +2.0*cof2*rad*(momen_p-cof1*pow(rad,2)/(2.0*q))/(cof2*q*RAD_0*(RAD_0+rad*cos(theta)))
             +cof2*pow(rad,2)*((-2.0*cof1*rad)/(2.0*q))/(cof2*q*RAD_0*(RAD_0+rad*cos(theta)))
			 +cof2*pow(rad,2)*(momen_p-cof1*pow(rad,2)/(2.0*q))*(-cos(theta))
			 /(cof2*q*RAD_0*pow(RAD_0+rad*cos(theta),2));
            
   func[1][1]=-cof1*log((RAD_0+rad*cos(theta))/(RAD_0))*0.5*RAD_0*rad*(-sin(theta))
	          -cof1*RAD_0*rad*cos(theta)*0.5*(RAD_0/(RAD_0+rad*cos(theta)))*(rad*sin(theta)/RAD_0)
			  -cof1*RAD_0*pow(rad,2)*sin(theta)*cos(theta)/(RAD_0+rad*cos(theta))
			  -RAD_0*pow(rad,2)*pow(sin(theta),2)*rad*sin(theta)/(2.0*pow(RAD_0+rad*cos(theta),2))
			  +cof2*pow(rad,2)*(momen_p-cof1*pow(rad,2)/(2.0*q))*rad*sin(theta)
			  /(cof2*q*RAD_0*pow(RAD_0+rad*cos(theta),2));

   func[2][0]=0.5*cof3*2.0*(rad/pow(q,2))/((RAD_0+rad*cos(theta))*sqrt(pow(rad/q,2)+pow(RAD_0,2)))
	         +cof3*sqrt(pow(rad/q,2)+pow(RAD_0,2))*(-cos(theta))/pow(RAD_0+rad*cos(theta),2)
			 +cof4*2*(momen_p-cof1*pow(rad,2)*0.5/q)*(-cof1*rad/q)*(pow(rad,2)+pow(RAD_0*q,2))
			 /pow(cof2*q*RAD_0*(RAD_0+rad*cos(theta)),2)
			 +cof4*pow(momen_p-cof1*pow(rad,2)*0.5/q,2)*(2.0*rad)/pow(cof2*q*RAD_0*(RAD_0+rad*cos(theta)),2)
			 +cof4*pow(momen_p-cof1*pow(rad,2)*0.5/q,2)*(pow(rad,2)+pow(RAD_0*q,2))*(-2.0*cos(theta))
			 /(pow(cof2*q*RAD_0,2)*pow(RAD_0+rad*cos(theta),3));

   func[2][1]=cof3*sqrt(pow(rad/q,2)+pow(RAD_0,2))*rad*sin(theta)/pow(RAD_0+rad*cos(theta),2)
	         +cof4*pow(momen_p-cof1*pow(rad,2)*0.5/q,2)*(pow(rad,2)+pow(RAD_0*q,2))*2.0*rad*sin(theta)
			 /(pow(cof2*q*RAD_0,2)*pow(RAD_0+rad*cos(theta),3));
   
   double fun_k[2];
   fun_k[0]=func[2][1]/(func[0][1]-func[1][0]);
   fun_k[1]=func[2][0]/(func[1][0]-func[0][1]);
   double *p;
   p=fun_k;
   return p;
 }

double *fun_2(double rad,double theta)
 { double *p;
   p=dif_fun(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=dif_fun(rad+0.5*h*f[0],theta+0.5*h*f[1]);
   return p;
 }

double *fun_3(double rad,double theta)
 { double *p;
   p=fun_2(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=dif_fun(rad+0.5*h*f[0],theta+0.5*h*f[1]);
   return p;
 }

double *fun_4(double rad,double theta)
 { double *p;
   p=fun_3(rad,theta);
   double f[2];
   f[0]=*p,f[1]=*(p+1);
   p=dif_fun(rad+h*f[0],theta+h*f[1]);
   return p;
 }

double *solv(double rad, double theta)
 { double *p;
   double f1[2],f2[2],f3[2],f4[2];
   p=dif_fun(rad,theta);
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