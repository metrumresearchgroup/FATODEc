#ifndef FATODEC_H
#define FATODEC_H

extern void integrate_fatode_fwd_erk( double* tin, double* tout,
                                      int* nvar, double var[],
                                      double rtol[], double atol[],
                                      void (*fun) (int*, double*, double[], double[]),
                                      int icntrl_u[],
                                      double rcntrl_u[],
                                      int istatus_u[],
                                      double rstatus_u[],
                                      int* ierr_u );

#endif
