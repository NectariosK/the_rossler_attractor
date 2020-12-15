/*NECTARIOS_KISIIGHA_CS2_PROJECT*/
/*STUDENT NO: 288701*/

/*The Rossler attractor*/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include "rk4.h"
#include "rossler.h" /* rossler(double t, double X[], double dXdt[]) */
#define N 3 /*dimension of state variable*/
#define h 0.01 /*Runge-Kutta step*/
#define T 10000 /*Numbers of caculation of Runge-Kutta*/

int main()
{
  int t, i;
  /* Dynamically allocating arrays for 'double X0[N],' and 'double X1[N]'. */
  double *X0 = vector(N),*X1 = vector(N);
  /* Initial setting */
  X0[0] = 10.0;
  X0[1] = 20.0;
  X0[2] = 30.0;
  /*Main part*/
  for(t=0; t<= T-1; t++)
    {
      for(i=0; i<=N-1; i++)
	{
	  printf("%f", X0[i]);
	  if(i == N - 1) putchar('\n');
	  else putchar(' ');
	}
      /*** Putting the function 'rossler' and the state variable 'X0' at time
      ** 'h*t' to get the state variable 'X1' at time 'h*(t+1)'.*/
      rk4(h,N, rossler, h*t, X0, X1);
      copy_vector(N, X1, X0);
	}
  return 0;
}
/*releasing the definitions */
#undef dxdt
#undef dydt
#undef dzdt
#undef x
#undef y
#undef z