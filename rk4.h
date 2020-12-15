/* Runge-Kutta routine (rk4.h)
** This header file is for resolving a rank-1 ordinary differential equations.
** The number of state variables of each equation is 'N'.
** This time, the Runge-Kutta method was used.
** Since using computer simulation cannot resolve analytically
** the differential equation 'dx(t) / dt = F(x(t), t),'
** the time 't' should be discretized.
** The order 4 Runge-Kutta method is as follows:
** d1 = hF(x(t), t),
** d2 = hF(x(t) + d1 / 2, t + h / 2),
** d3 = hF(x(t) + d2 / 2, t + h / 2),
** d4 = hF(x(t) + d3, t + h),
** x(t + h) = x(t) + d1 / 6 + d2 / 3 + d3 / 3 + d4 / 6,
** where 'h' is the size of discretized time, and is the Runge-Kutta step.*/
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>

double *vector(int N) /* Allocating vector region */
{
  /*** Dynamically allocating an array using "malloc".
  ** The size of the array is 'N'.
  */
  return (double *)malloc(N * sizeof(double));
}

void free_vector(double *v) /* Releasing vector region */
{
  /* Releasing the region of the pointer 'v'. */
  free(v);
}

void copy_vector(int N, double a[], double b[]) /* Copying vectors */
{
  int i;
  /* Copying the vector 'a' to the vector 'b' */
  for(i = 0; i <= N - 1; i++) b[i] = a[i];
}
/*** Putting the Runge-Kutta step 'h,' the rank of the differential equation N,
** the function 'dXdt' which describes the differential equation,
** and the vector 'X0' which contains the state variable at time 't'
** to get the state variable 'X1' at time 't + h'.
** 'dXdt,' which type is 'foo(double t, double X[], double dXdt[])' should be
** defined by user.*/

 void rk4(double h, int N, void(*dXdt)(double t, double X[], double dXdt[]),double t, double X0[], double X[]) /* Runge-Kutta */
{
  int i;
  /*** Dynamically allocating arrays for 'double d1[N],' 'double d2[N],'
  ** and 'double d3[N]'.  */
  double *d1 = vector(N),*d2 = vector(N),*d3 = vector(N);

   /** Dynamically allocating arrays for 'double Xa[N],' and 'double X[N]'.*/
  double *Xa = vector(N),*dX = vector(N);
  /* d1 = hF(x(t), t) */
  dXdt(t, X0, dX);
  for(i = 0; i <= N - 1; i++)
    {
      d1[i] = h * dX[i];
      Xa[i] = X0[i] + 0.5 * d1[i];
    }
  /* d2 = hF(x(t) + d1 / 2, t + h / 2) */
  dXdt((t + (0.5*h)),Xa,dX);
  for(i = 0; i <= N - 1; i++)
    {
      d2[i] = h * dX[i];
      Xa[i] = X0[i] + 0.5 * d2[i];
    }
  /* d3 = hF(x(t) + d2 / 2, t + h / 2) */ 
  dXdt(t + 0.5 * h, Xa, dX);
  for(i = 0; i <= N-1; i++)
    {
      d3[i] = h * dX[i];
      Xa[i] = X0[i] + d3[i];
    }
  /* x(t + h) = x(t) + d1 / 6 + d2 / 3 + d3 / 3 */
  dXdt(t + h, Xa, dX);
  for(i = 0; i <= N - 1; i++)
    X[i] = X0[i] + (d1[i] + d2[i] * 2 + d3[i] * 2 + h * dX[i]) / 6.0;
  /* Releasing the regions of arrays */
  free_vector(d1); free_vector(d2); free_vector(d3);
  free_vector(Xa); free_vector(dX);
}