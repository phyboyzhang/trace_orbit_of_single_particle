#include <math.h>
#include "fun_class.h"

extern double rad_0,RAD_0,theta_0,psi_0,cof1,cof2,cof3,cof4,q,h,del,epsilon;

//这个函数求得微分系数
double *difmatrix(double cor[4],double cor1[4])
 { double avercor[4];
   avercor[0]=(cor[0]+cor1[0])/2.0;
   avercor[1]=(cor[1]+cor1[1])/2.0;
   avercor[2]=(cor[2]+cor1[2])/2.0;
   avercor[3]=(cor[3]+cor1[3])/2.0;
   double *p;
   p=func(avercor[0],avercor[1],avercor[2],avercor[3]);
   double fun1[4], fun2[4];
   for(int i=0;i<=3;i++)
	fun1[i]=*(p+i);
   double diffun_k[4][4];

   for(int i=0;i<=3;i++)
   { avercor[i]=(cor[i]+cor1[i])/2.0+del/2.0;
     for(int j=0;j<=3;j++)
	 { if (j!=i)
	   avercor[j]=(cor[j]+cor1[j])/2.0;
	 }
     p=func(avercor[0],avercor[1],avercor[2],avercor[3]);
	 for(int j=0;j<=3;j++)
	 { fun2[j]=*(p+j);
	   diffun_k[j][i]=(fun2[j]-fun1[j])/del;
	 }
   }
   double diff[16];
   for(int i=0;i<=3;i++)
   { for(int j=0;j<=3;j++)
      diff[4*i+j]=diffun_k[i][j];
   }
   p=diff;
   return p;
 }

//这个函数求得Euler-Lagrangian方程
double *eul_lag(double cor[4],double cor1[4],double cor2[4])
 { double *p,*p1;
   double cor3[4],cor4[4],fun1[4],fun2[4];
   for(int i=0;i<=3;i++)
    { cor3[i]=(cor[i]+cor1[i])/2.0;
      cor4[i]=(cor1[i]+cor2[i])/2.0;
    }
   p=func(cor3[0],cor3[1],cor3[2],cor3[3]);
   for(int i=0;i<=3;i++)
     fun1[i]=*(p+i);
   p1=func(cor4[0],cor4[1],cor4[2],cor4[3]);
   for(int i=0;i<=3;i++)
     fun2[i]=*(p1+i);
   p=difmatrix(cor,cor1);
   double diff_k[4][4],diff_k1[4][4];
   for(int i=0;i<=3;i++)
    { for(int j=0;j<=3;j++)
		diff_k[i][j]=*(p+4*i+j);       
    }
   p1=difmatrix(cor1,cor2);
   for(int i=0;i<=3;i++)
    { for(int j=0;j<=3;j++)
        diff_k1[i][j]=*(p1+4*i+j);      
    }
   double elequ[4];
   for (int j=0;j<=2;j++)
    { double sum=0;
      for(int i=0;i<=2;i++)
        sum=sum+diff_k[i][j]*(cor1[i]-cor[i])+diff_k1[i][j]*(cor2[i]-cor1[i]);
	  elequ[j]=sum-h*diff_k[3][j]-h*diff_k1[3][j]+fun1[j]-fun2[j];
    }
   elequ[3]=0;
   for(int i=0;i<=2;i++)
    { elequ[3]=elequ[3]+diff_k[i][3]*(cor1[i]-cor[i])+diff_k1[i][3]*(cor2[i]-cor1[i]);
    }
   elequ[3]=elequ[3]-h*diff_k[3][3]-h*diff_k1[3][3];
   p=elequ;
   return p;
 }

//euler-lagrangian方程的导数
double *NRdiff(double cor[4],double cor1[4],double cor2[4])
 { double *p;
   p=eul_lag(cor,cor1,cor2);
   double elequ_0[4];
   for(int i=0;i<=3;i++)
	 elequ_0[i]=*(p+i);
   double cor3[4];
   for (int k=0;k<=3;k++)
	 cor3[k]=cor2[k];
   double sol_matr[4][4]; //   euler-lagrangian equations 的 jiacob矩阵
   for(int i=0;i<=3;i++)
    { cor2[i]=cor2[i]+del;
      for(int j=0;j<=3;j++)
	  { if(j!=i)
	    cor2[j]=cor3[j];
	  }
	  p=eul_lag(cor,cor1,cor2);
	  double elequ_1[4];
	  for (int n=0;n<=3;n++)
	   {elequ_1[n]=*(p+n);
	    sol_matr[n][i]=(elequ_1[n]-elequ_0[n])/del;
	   }
    }
    double sol[16];
	for (int i=0;i<=3;i++)
	 for(int j=0;j<=3;j++)
       sol[4*i+j]=sol_matr[i][j];
	p=sol;
	return p;
}
// newton-raphson求解方程

void NRSOLVE(double cor[4],double cor1[4])
 { double *p;
   double cor2[4];
   for(int i=0;i<=3;i++)
    cor2[i]=cor1[i];
   double sum=1;
   while (sum>=epsilon)
    { sum=0;
	  double difmar[16];
      p=NRdiff(cor,cor1,cor2);
	  for(int i=0;i<=15;i++)
		difmar[i]=*(p+i);
	  p=eul_lag(cor,cor1,cor2);
	  double fun[4];
	  for(int i=0;i<=3;i++)
		{fun[i]=*(p+i);
//	    sum=sum+fabs(fun[i]);
	    }
      agaus(difmar,fun,4);
	  for(int i=0;i<=3;i++)
      cor2[i]=cor2[i]-fun[i];
	  p=eul_lag(cor,cor1,cor2);
	  for(int i=0;i<=3;i++)
		sum=sum+fabs(*(p+i));
    } 
   for(int i=0;i<=3;i++)
	 cor1[i]=cor2[i];
 }

 
