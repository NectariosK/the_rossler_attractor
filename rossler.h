/*Rossler.h */
/* The Rossler equations (rossel.h) */
/*
** This header file describes function type of the Rossler equations
** for the Runge-Kutta routine (rk4.h).
** The Rossler equations derived by simplifying of convection rolls arising
** in the atmosphere, and described as ordinary differential equations
** that have 3 variables:
** dx / dt = -y - z,
** dy / dt = x + ay,
** dz / dt = b + zx - zc,
** ,where a, b, and c are parameters, and they are set to
** as follows:
** a = b = 0.2, and c = 5.7.
** The Rossler equations exhibit sensitive dependence on initial conditions.*/


/* The notation of functions and variables are consistent with 'rk4.h'. */
#define _CRT_SECURE_NO_WARNINGS
#define dxdt dXdt[0]
#define dydt dXdt[1]
#define dzdt dXdt[2]
#define x X[0]
#define y X[1]
#define z X[2]

/* The Rossler equations */
void rossler(double t, double X[], double dXdt[])
{
  double a=0.2,b=0.2,c=5.7;
  dxdt = -y - z;
  dydt = x + a*y;
  dzdt = b + z*x - z*c;
}