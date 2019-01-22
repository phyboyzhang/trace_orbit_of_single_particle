#include <math.h>
#include "fun_class.h"

extern double rad_0,RAD_0,theta_0,psi_0,cof1,cof2,cof3,cof4,q,h,del;

double *dif_matr(double rad,double theta,double psi, double u)
 { double dif[4][4];

   dif[0][0]=- (0.5*cof1*rad*RAD_0*pow(cos(theta),2)*sin(theta))/pow((RAD_0 + rad*cos(theta)),2);
   
   dif[0][1]=(0.5*cof1*rad*RAD_0*pow(cos(theta),2))/(RAD_0+rad*cos(theta)) 
            -0.5*cof1*RAD_0*cos(theta)*log((RAD_0+rad*cos(theta))/RAD_0) + 
             (0.5*cof1*pow(rad,2)*RAD_0*cos(theta)*
             pow(sin(theta),2))/pow(RAD_0 + rad*cos(theta),2);
   
   dif[0][2]=0;
   
   dif[0][3]=0;

   dif[1][0]=-(cof2*pow(rad,3)*u)/pow(pow(rad,2) + pow(q*RAD_0,2),1.5) + (2.0*cof2*rad*u)/
   sqrt(pow(rad,2) + pow(q*RAD_0,2)) - 
   (0.5*cof1*rad*RAD_0*pow(cos(theta),2))/(RAD_0 + rad*cos(theta)) - 
   0.5*cof1*RAD_0*cos(theta)*log((RAD_0 + rad*cos(theta))/RAD_0) + 
   (0.5*cof1*pow(rad,2)*RAD_0*cos(theta)*
    pow(sin(theta),2))/pow(RAD_0 + rad*cos(theta),2) - 
   (1.0*cof1*rad*RAD_0*pow(sin(theta),2))/(RAD_0 + rad*cos(theta));
   
   dif[1][1]=-0.5*cof1*pow(rad,2)*RAD_0*cos(theta)*sin(theta)/(RAD_0 + rad*cos(theta)) + 
   0.5*cof1*rad*RAD_0*log((RAD_0 + rad*cos(theta))/RAD_0)*sin(theta) - 
   (0.5*cof1*pow(rad,3)*RAD_0*pow(sin(theta),3))/pow((RAD_0 + rad*cos(theta)),2);
   
   dif[1][2]=0;
   
   dif[1][3]=(cof2*pow(rad,2))/sqrt(pow(rad,2) + pow(q*RAD_0,2));

   dif[2][0]=(RAD_0 + rad*
             cos(theta))*(-((cof2*q*rad*RAD_0*u)/pow((pow(rad,2) + pow(q*RAD_0,2)),1.5)) - 
             (0.5*cof1*pow(rad,2)*cos(theta))/(q*pow((RAD_0 + rad*cos(theta)),2)) + 
             (1.0*cof1*rad)/(q*(RAD_0 + rad*cos(theta)))) + 
              cos(theta)*((cof2*q*RAD_0*u)/sqrt(pow(rad,2) + pow(q*RAD_0,2)) + 
             (0.5*cof1*pow(rad,2))/(q*(RAD_0 + rad*cos(theta))));

   dif[2][1]=-cof2*q*rad*RAD_0*u*sin(theta)/sqrt(pow(rad,2)+pow(q*RAD_0,2));
	   
	  /* (0.5*cof1*pow(rad,3)*sin(theta))/(q*(RAD_0 + rad*cos(theta))) - 
             rad*((cof2*q*RAD_0*u)/sqrt(pow(rad,2) + pow(q*RAD_0,2)) + (0.5*cof1*pow(rad,2))/
             (q*(RAD_0 + rad*cos(theta))))*sin(theta);*/

   dif[2][2]=0;

   dif[2][3]=(cof2*q*RAD_0*(RAD_0 + rad*cos(theta)))/sqrt(pow(rad,2) + pow(q*RAD_0,2));

   dif[3][0]=-(cof3*sqrt(pow(rad/q,2) + pow(RAD_0,2))*
              cos(theta))/pow((RAD_0 + rad*cos(theta)),2)+ 
             (cof3*rad)/(pow(q,2)*sqrt(pow(rad/q,2) + pow(RAD_0,2))
			 *(RAD_0 + rad*cos(theta)));

   dif[3][1]=(cof3*rad*sqrt(pow(rad/q,2) + pow(RAD_0,2))*sin(theta))/pow((RAD_0 + rad*cos(theta)),2);
   
   dif[3][2]=0;

   dif[3][3]=2.0*cof4*u;
   
   double fun[16];
   for(int i=0;i<=3;i++)
	 for(int j=0;j<=3;j++)
	   fun[4*i+j]=dif[i][j];
   double *p;
   p=fun;
   return p;
}

double *gb_solve(double cor1[4],double cor2[4])
 { double *p;
   double a[4][4];
   p=dif_matr(cor2[0],cor2[1],cor2[2],cor2[3]);
   for(int i=0;i<=3;i++)
     { for(int j=0;j<=3;j++)
	     a[i][j]=*(p+i*4+j);
     }
   double tran[4][4];
   tran[0][0]=0,tran[1][1]=0,tran[2][2]=0,tran[3][3]=0;
   tran[0][1]=a[1][0]-a[0][1];
   tran[0][2]=a[2][0]-a[0][2];
   tran[1][2]=a[2][1]-a[1][2];
   tran[0][3]=-a[0][3];
   tran[1][3]=-a[1][3];
   tran[2][3]=-a[2][3];
   tran[1][0]=-tran[0][1];
   tran[2][0]=-tran[0][2];
   tran[2][1]=-tran[1][2];
   tran[3][0]=-tran[0][3];
   tran[3][1]=-tran[1][3];
   tran[3][2]=-tran[2][3];

   double hdiff[4];
   for(int i=0;i<=3;i++)
   {hdiff[i]=4*h*a[3][i];}
   double tr[16];
   for(int i=0;i<=3;i++)
	 for(int j=0;j<=3;j++)
	  tr[i*4+j]=tran[i][j];
   agaus(tr,hdiff,4);  //解就是返回的hdiff
   for(int i=0;i<=3;i++)
   cor2[i]=hdiff[i]+cor1[i];

 //  double *p1;
   p=cor2;
   return p;
 }
