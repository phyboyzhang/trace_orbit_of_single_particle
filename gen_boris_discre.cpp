#include <math.h>
#include "fun_class.h"
#include <iostream>
//#include "parameters.h"
using namespace std;

extern double rad_0,RAD_0,theta_0,psi_0,cof1,cof2,cof3,cof4,q,h,del;

/*double coff3()
 { double cof=sqrt(pow(rad_0,2)/pow(q,2)+pow(RAD_0,2))/(RAD_0+rad_0*cos(theta_0));
   return cof;
 }*/

double *func(double rad,double theta,double psi,double u)
 { double funcs[4];
   funcs[0]=-cof1*log((RAD_0+rad*cos(theta))/RAD_0)*0.5*RAD_0*sin(theta)
	      +cof1*0.5*RAD_0*rad*sin(theta)*cos(theta)/(RAD_0+rad*cos(theta));
   funcs[1]=-cof1*log((RAD_0+rad*cos(theta))/RAD_0)*0.5*RAD_0*rad*cos(theta)
	      -cof1*RAD_0*pow(rad,2)*pow(sin(theta),2)*0.5/(RAD_0+rad*cos(theta))
		  +cof2*u*pow(rad,2)/sqrt(pow(rad,2)+pow(RAD_0*q,2));
   funcs[2]=cof1*pow(rad,2)*0.5*(RAD_0+rad*cos(theta))/(q*(RAD_0+rad*cos(theta)))
	      +cof2*RAD_0*q*u*(RAD_0+rad*cos(theta))/sqrt(pow(rad,2)+pow(RAD_0*q,2));
   funcs[3]=cof3*sqrt(pow(rad,2)/pow(q,2)+pow(RAD_0,2))/(RAD_0+rad*cos(theta))
	      +cof4*pow(u,2);
  double *p;
   p=funcs;
   double exp[4];
   for(int i=0;i<4;i++)
	 exp[i]=*(p+i);
   return p;
 }



double *A(double rad,double theta, double psi, double u)
 { double a[4][4];
   double *p1;
   double *p2;
   p1=func(rad,theta,psi,u);
   double exp1[4];
   for(int i=0;i<4;i++)
	 exp1[i]=*(p1+i);
   p2=func(rad+del,theta,psi,u);
   double exp2[4];
   for(int i=0;i<4;i++)
	 exp2[i]=*(p1+i);
   a[0][0]=(exp2[0]-exp1[0])/del;
   a[1][0]=(exp2[1]-exp1[1])/del;
   a[2][0]=(exp2[2]-exp1[2])/del;
   a[3][0]=(exp2[3]-exp1[3])/del;
   p2=func(rad,theta+del,psi,u);
   for(int i=0;i<4;i++)
	 exp2[i]=*(p1+i);
   a[0][1]=(exp2[0]-exp1[0])/del;
   a[1][1]=(exp2[1]-exp1[1])/del;
   a[2][1]=(exp2[2]-exp1[2])/del;
   a[3][1]=(exp2[3]-exp1[3])/del;
   p2=func(rad,theta,psi+del,u);
   for(int i=0;i<4;i++)
	 exp2[i]=*(p1+i);
   a[0][2]=(exp2[0]-exp1[0])/del;
   a[1][2]=(exp2[1]-exp1[1])/del;
   a[2][2]=(exp2[2]-exp1[2])/del;
   a[3][2]=(exp2[3]-exp1[3])/del;
   p2=func(rad,theta,psi,u+del);
   for(int i=0;i<4;i++)
	 exp2[i]=*(p1+i);
   a[0][3]=(exp2[0]-exp1[0])/del;
   a[1][3]=(exp2[1]-exp1[1])/del;
   a[2][3]=(exp2[2]-exp1[2])/del;
   a[3][3]=(exp2[3]-exp1[3])/del;


  double b[16];
  for(int i=0;i<=3;i++)
   { for(int j=0;j<=3;j++)
     b[i*4+j]=a[i][j]; 
   }
//  double *p;
  p1=b;
  return(p1);
 }
// 求k+1步的解(L,theta,phi,U)
double *solve(double cor1[4],double cor2[4])
 { double *p;
   double a[4][4];
   p=A(cor2[0],cor2[1],cor2[2],cor2[3]);
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
