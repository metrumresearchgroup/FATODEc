#ifndef FATODEC_HPP
#define FATODEC_HPP

extern "C" {

  void integrate_fatode_fwd_erk( double* tin, double* tout,
                                 int* nvar, double var[],
                                 double rtol[], double atol[],
                                 void (*fun) (int*, double*, double[], double[]),
                                 int icntrl_u[],
                                 double rcntrl_u[],
                                 int istatus_u[],
                                 double rstatus_u[],
                                 int* ierr_u );

  void integrate_fatode_fwd_rk( double* tin, double* tout,
                                int* nvar, int* nnzero, double var[],
                                double rtol[], double atol[],
                                void (*fun) (int*, double*, double[], double[]),
                                void (*jac) (int*, double*, double[], double[]),
                                int icntrl_u[],
                                double rcntrl_u[],
                                int istatus_u[],
                                double rstatus_u[],
                                int* ierr_u );

  void integrate_fatode_fwd_ros( double* tin, double* tout,
                                 int* nvar, int* nnzero, double var[],
                                 double rtol[], double atol[],
                                 void (*fun) (int*, double*, double[], double[]),
                                 void (*jac) (int*, double*, double[], double[]),
                                 int icntrl_u[],
                                 double rcntrl_u[],
                                 int istatus_u[],
                                 double rstatus_u[],
                                 int* ierr_u );

  void integrate_fatode_fwd_sdirk( double* tin, double* tout,
                                   int* nvar, int* nnzero, double var[],
                                   double rtol[], double atol[],
                                   void (*fun) (int*, double*, double[], double[]),
                                   void (*jac) (int*, double*, double[], double[]),
                                   int icntrl_u[],
                                   double rcntrl_u[],
                                   int istatus_u[],
                                   double rstatus_u[],
                                   int* ierr_u );

}

#endif
