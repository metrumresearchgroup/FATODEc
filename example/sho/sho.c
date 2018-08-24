#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <fatodec.h>

static void f(int*, double*, double[], double[]);

int main(int argc, char *argv[])
{

  double tin = 0.0;
  double tout = 10.0;
  int nvar = 2;
  double   fy[2] = { 0.0, 2.0 };
  double rtol[2] = { 1.E-6, 1.E-6 };
  double atol[2] = { 1.E-10, 1.E-10 };
  int icntrl_u[20];
  double rcntrl_u[20];
  int istatus_u[20];
  double rstatus_u[20];
  int ierr_u;

  icntrl_u[0] = 0;
  icntrl_u[1] = 0;
  icntrl_u[2] = 4;
  icntrl_u[3] = 1000;

  rcntrl_u[0] = 1.E-4;
  rcntrl_u[1] = 0.1;
  rcntrl_u[2] = 0.1;
  rcntrl_u[3] = 0.2;
  rcntrl_u[4] = 6.0;
  rcntrl_u[5] = 0.1;
  rcntrl_u[6] = 0.9;
  rcntrl_u[9] = 1.0;
  rcntrl_u[10] = 1.2;

  integrate_fatode_fwd_erk( &tin, &tout, &nvar, fy, rtol, atol, f, icntrl_u, rcntrl_u, istatus_u, rstatus_u, &ierr_u );    

  printf("No. of function calls                          :  %d\n", istatus_u[0]);
  printf("No. of jacobian calls                          :  %d\n", istatus_u[1]);
  printf("No. of steps                                   :  %d\n", istatus_u[2]);
  printf("No. of accepted steps                          :  %d\n", istatus_u[3]);
  printf("No. of rejected steps (except at the beginning :  %d\n", istatus_u[4]);
  printf("No. of LU decompositions                       :  %d\n", istatus_u[5]);
  printf("No. of forward/backward substitutions          :  %d\n", istatus_u[6]);
  printf("No. of singular matrix decompositions          :  %d\n", istatus_u[7]);
                                                         
  printf("T at exit                                      :  %f\n", rstatus_u[0]);
  printf("H at exit                                      :  %f\n", rstatus_u[1]);
  printf("H new                                          :  %f\n", rstatus_u[2]);

  return 0;
}

static void f(int* n, double* t, double y[], double fy[]) {
  double theta = 0.15;
  fy[0] = y[1];
  fy[1] = -y[0] - theta * y[1];
} 
