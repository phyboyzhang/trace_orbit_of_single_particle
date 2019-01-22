#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "fun_class.h"
using namespace std;

//  v=1e+5
double  rad_0=1,RAD_0=4,theta_0=0,psi_0=0,q=2.0,u_0=sqrt(0.01),v_0=1e+7,h=0.02,epsilon=1e-9; /*h=0.02 T=300*/
double  cof1,cof2,cof3,cof4;
double  del=1e-6;
double  momen_p;



int main()
{ double asp=rad_0/RAD_0;
	cof1=((1.0e-10)/3.0)*sqrt(asp)*v_0/(q*(RAD_0+rad_0));
	cof2=(9.1/4.8)*(1e-18)*pow(v_0,2)*sqrt(asp)/(q*(RAD_0+rad_0));
	cof3=0.5*(9.1/4.8)*(1e-18)*pow(v_0,2)*(RAD_0+rad_0*cos(theta_0))*(1-pow(u_0,2))
		/sqrt(pow(rad_0/q,2)+pow(RAD_0,2));
	cof4=0.5*(9.1/4.8)*(1e-18)*pow(v_0,2);
	
  double *p;
  double cor1[4]={rad_0,0.0,0.0,u_0},cor2[4];

  momen_p=momen_0(rad_0,theta_0);
  p=solv(cor1[0],cor1[1]);
  cor2[0]=*p,cor2[1]=*(p+1);
    cor2[3]=(U(cor2[0],cor2[1])+u_0)/2;
  cor2[2]=cor2[3]*h/(RAD_0+rad_0);

  double com1[4],com2[4];
  for(int i=0;i<=3;i++)
	com1[i]=cor1[i],com2[i]=cor2[i];

  int n=68000;
  double cor3[4],cor4[4];

  ofstream outfile;
  // general boris
//  outfile.open("gb_orbit.txt",ios::out);
  outfile.open("orbit.txt",ios::out);
  for(int i=0;i<=-1;i++)
   { for(int j=0;j<=3;j++)
	 cor3[j]=cor2[j],cor4[j]=cor1[j];
	 p=solve(cor1,cor2);
	 for(int k=0;k<=3;k++)
	 cor2[k]=*(p+k);
	 for(int j=0;j<=3;j++)
	 cor1[j]=cor3[j];
	 double x,y;
/*	 x=(RAD_0+cor2[0]*cos(cor2[1]))*cos(cor2[2]);%
	 y=(RAD_0+cor2[0]*cos(cor2[1]))*sin(cor2[2]);
	 z=cor2[0]*sin(cor2[1]);*/

     x=RAD_0+cor2[0]*cos(cor2[1]);
	 y=cor2[0]*sin(cor2[1]);

	 double Eva=(1-pow(u_0,2))*(sqrt(pow(cor2[0]/q,2)+pow(RAD_0,2))/(RAD_0+cor2[0]*cos(cor2[1])))
		       /(sqrt(pow(rad_0/q,2)+pow(RAD_0,2))/(RAD_0+rad_0*cos(theta_0)))+pow(cor2[3],2);
		      
	 outfile<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)
		 <<setprecision(5);
	 outfile<<setw(20)<<x<<setw(20)<<y<<setw(20)<<Eva<<setw(20)<<i<<endl;
	 
/*	 p=gb_solve(cor4,cor3);
	 for(int j=0;j<=3;j++)
		cor3[j]=*(p+j);*/
   }
  outfile.close(); 

// vsa
  for(int i=0;i<=3;i++)
	cor1[i]=com1[i],cor2[i]=com2[i];

  outfile.open("vsaorbit.txt",ios::out);
  for(int i=0;i<=n;i++)
   { double cor3[4];
     for(int j=0;j<=3;j++)
		cor3[j]=cor2[j];
	 NRSOLVE(cor1,cor2);
     double x=RAD_0+cor2[0]*cos(cor2[1]);
	 double y=cor2[0]*sin(cor2[1]);

	 double Eva=(1-pow(u_0,2))*(sqrt(pow(cor2[0]/q,2)+pow(RAD_0,2))/(RAD_0+cor2[0]*cos(cor2[1])))
		       /(sqrt(pow(rad_0/q,2)+pow(RAD_0,2))/(RAD_0+rad_0*cos(theta_0)))+pow(cor2[3],2);;

	 outfile<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)
		 <<setiosflags(ios::fixed)<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)
		 <<setiosflags(ios::fixed)<<setprecision(5);

     outfile<<setw(20)<<x<<setw(20)<<y<<setw(20)<<cor2[0]<<setw(20)<<cor2[1]<<setw(20)<<cor2[2]
	 <<setw(20)<<cor2[3]<<setw(20)<<Eva<<setw(20)<<i<<endl;

	 for(int j=0;j<=3;j++)
		cor1[j]=cor3[j];
   }
  outfile.close();   

// RK4 

    for(int i=0;i<=3;i++)
	cor1[i]=com1[i],cor2[i]=com2[i];

  outfile.open("rkorbit.txt",ios::out);
  double rad=rad_0,theta=theta_0;
  for(int i=0;i<=-1;i++)
   { p=solv(rad,theta);
     rad=*p,theta=*(p+1);
     double x=RAD_0+rad*cos(theta);
	 double y=rad*sin(theta);
	 double vpar=U(rad,theta);
	 double mom=momen(rad,theta);
	 double Er=(1-pow(u_0,2))*sqrt(pow(rad/q,2)+pow(RAD_0,2))*pow(RAD_0+rad*cos(theta),-1)
		      *pow(RAD_0+rad_0*cos(theta_0),1)/sqrt(pow(rad_0/q,2)+pow(RAD_0,2))+pow(vpar,2);
     outfile<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)<<setiosflags(ios::fixed)
		    <<setiosflags(ios::fixed)<<setiosflags(ios::fixed)<<setprecision(5);
	 outfile<<setw(15)<<x<<setw(15)<<y<<setw(15)<<mom<<setw(15)<<Er<<setw(15)<<i<<endl;
   }
  outfile.close(); 

  return 0;
}