#ifndef FATODEC_HPP
#define FATODEC_HPP

#include <vector>

extern "C" {

  void integrate_fatode_fwd_erk( double* tin, double* tout,
                                 int* nvar, double var[],
                                 double rtol[], double atol[],
                                 void (*fun) (int*, double*, double[], double[]),
                                 int icntrl_u[], double rcntrl_u[],
                                 int istatus_u[], double rstatus_u[],
                                 int* ierr_u );

  void integrate_fatode_fwd_rk( double* tin, double* tout,
                                int* nvar, int* nnzero, double var[],
                                double rtol[], double atol[],
                                void (*fun) (int*, double*, double[], double[]),
                                void (*jac) (int*, double*, double[], double[]),
                                int icntrl_u[], double rcntrl_u[],
                                int istatus_u[], double rstatus_u[],
                                int* ierr_u );

  void integrate_fatode_fwd_ros( double* tin, double* tout,
                                 int* nvar, int* nnzero, double var[],
                                 double rtol[], double atol[],
                                 void (*fun) (int*, double*, double[], double[]),
                                 void (*jac) (int*, double*, double[], double[]),
                                 int icntrl_u[], double rcntrl_u[],
                                 int istatus_u[], double rstatus_u[],
                                 int* ierr_u );

  void integrate_fatode_fwd_sdirk( double* tin, double* tout,
                                   int* nvar, int* nnzero, double var[],
                                   double rtol[], double atol[],
                                   void (*fun) (int*, double*, double[], double[]),
                                   void (*jac) (int*, double*, double[], double[]),
                                   int icntrl_u[], double rcntrl_u[],
                                   int istatus_u[], double rstatus_u[],
                                   int* ierr_u );

  void integrate_fatode_tlm_erk( double* tin, double* tout,
                                 int* nvar, int* ntlm, double var[], double var_tlm[],
                                 double rtol_tlm[], double atol_tlm[],
                                 double rtol[], double atol[],
                                 void (*fun) (int*, double*, double[], double[]),
                                 void (*jac) (int*, double*, double[], double[]),                                      
                                 int icntrl_u[], double rcntrl_u[],
                                 int istatus_u[], double rstatus_u[],
                                 int* ierr_u );

  // void integrate_fatode_tlm_rk ( double* tin, double* tout,
  //                                int* nvar, int* ntlm, int* nnzero, double var[], double var_tlm[],
  //                                double rtol_tlm[], double atol_tlm[],
  //                                double rtol[], double atol[],
  //                                void (*fun) (int*, double*, double[], double[]),
  //                                void (*jac) (int*, double*, double[], double[]),                                      
  //                                int icntrl_u[], double rcntrl_u[],
  //                                int istatus_u[], double rstatus_u[],
  //                                int* ierr_u );

  void integrate_fatode_tlm_ros ( double* tin, double* tout,
                                  int* nvar, int* ntlm, int* nnzero, double var[], double var_tlm[],
                                  double rtol_tlm[], double atol_tlm[],
                                  double rtol[], double atol[],
                                  void (*fun) (int*, double*, double[], double[]),
                                  void (*jac) (int*, double*, double[], double[]),
                                  void (*hess_vec) (int*, double*, double[], double[], double[], double[]),
                                  int icntrl_u[], double rcntrl_u[],
                                  int istatus_u[], double rstatus_u[],
                                  int* ierr_u );

  void integrate_fatode_tlm_sdirk ( double* tin, double* tout,
                                    int* nvar, int* ntlm, int* nnzero, double var[], double var_tlm[],
                                    double rtol_tlm[], double atol_tlm[],
                                    double rtol[], double atol[],
                                    void (*fun) (int*, double*, double[], double[]),
                                    void (*jac) (int*, double*, double[], double[]),                                      
                                    int icntrl_u[], double rcntrl_u[],
                                    int istatus_u[], double rstatus_u[],
                                    int* ierr_u );

  // void integrate_fatode_adj_erk ( double* tin, double* tout,
  //                                 int* nvar, int* np, int* nadj, int* nnz,
  //                                 double var[], double lambda[],
  //                                 double rtol[], double atol[],
  //                                 void (*fun) (int*, double*, double[], double[]),
  //                                 void (*jac) (int*, double*, double[], double[]),                                      
  //                                 void (*adjinit) (int*, int*, int*, double*, double[], double[], double[]),
  //                                 int icntrl_u[], double rcntrl_u[],
  //                                 int istatus_u[], double rstatus_u[],
  //                                 int* ierr_u,
  //                                 double mu[],
  //                                 void (*jacp) (int*, int*, double*, double[], double[]),
  //                                 void (*drdy) (int*, int*, int*, double*, double[], double[]),
  //                                 void (*drdp) (int*, int*, int*, double*, double[], double[]),
  //                                 void (*qfun) (int*, int*, double*, double[], double[]),
  //                                 double q[] );

  // void integrate_fatode_adj_rk  ( double* tin, double* tout,
  //                                 int* nvar, int* np, int* nadj, int* nnz,
  //                                 double var[], double lambda[],
  //                                 double rtol_adj[], double atol_adj[],
  //                                 double rtol[], double atol[],
  //                                 void (*fun) (int*, double*, double[], double[]),
  //                                 void (*jac) (int*, double*, double[], double[]),                                      
  //                                 void (*adjinit) (int*, int*, int*, double*, double[], double[], double[]),
  //                                 int icntrl_u[], double rcntrl_u[],
  //                                 int istatus_u[], double rstatus_u[],
  //                                 int* ierr_u,
  //                                 double mu[],
  //                                 void (*jacp) (int*, int*, double*, double[], double[]),
  //                                 void (*drdy) (int*, int*, int*, double*, double[], double[]),
  //                                 void (*drdp) (int*, int*, int*, double*, double[], double[]),
  //                                 void (*qfun) (int*, int*, double*, double[], double[]),
  //                                 double q[] );

  // void integrate_fatode_adj_ros  ( double* tin, double* tout,
  //                                  int* nvar, int* np, int* nadj, int* nnz,
  //                                  double var[], double lambda[],
  //                                  double rtol_adj[], double atol_adj[],
  //                                  double rtol[], double atol[],
  //                                  void (*fun) (int*, double*, double[], double[]),
  //                                  void (*jac) (int*, double*, double[], double[]),                                      
  //                                  void (*adjinit) (int*, int*, int*, double*, double[], double[], double[]),
  //                                  void (*hess_vec) (int*, double*, double[], double[], double[], double[]),
  //                                  int icntrl_u[], double rcntrl_u[],
  //                                  int istatus_u[], double rstatus_u[],
  //                                  int* ierr_u,
  //                                  double mu[],
  //                                  void (*jacp) (int*, int*, double*, double[], double[]),
  //                                  void (*drdy) (int*, int*, int*, double*, double[], double[]),
  //                                  void (*drdp) (int*, int*, int*, double*, double[], double[]),
  //                                  void (*hesstr_vec_f_py) (int*, int*, double*, double[], double[], double[], double[]),
  //                                  void (*hesstr_vec_r_py) (int*, int*, double*, double[], double[], double[], double[]),
  //                                  void (*hesstr_vec_r) (int*, int*, double*, double[], double[], double[], double[]),
  //                                  void (*qfun) (int*, int*, double*, double[], double[]),
  //                                  double q[] );

  // void integrate_fatode_adj_sdirk ( double* tin, double* tout,
  //                                   int* nvar, int* np, int* nadj, int* nnz,
  //                                   double var[], double lambda[],
  //                                   double rtol_adj[], double atol_adj[],
  //                                   double rtol[], double atol[],
  //                                   void (*fun) (int*, double*, double[], double[]),
  //                                   void (*jac) (int*, double*, double[], double[]),            
  //                                   void (*adjinit) (int*, int*, int*, double*, double[], double[], double[]),
  //                                   int icntrl_u[], double rcntrl_u[],
  //                                   int istatus_u[], double rstatus_u[],
  //                                   int* ierr_u,
  //                                   double mu[],
  //                                   void (*jacp) (int*, int*, double*, double[], double[]),
  //                                   void (*drdy) (int*, int*, int*, double*, double[], double[]),
  //                                   void (*drdp) (int*, int*, int*, double*, double[], double[]),
  //                                   void (*qfun) (int*, int*, double*, double[], double[]),
  //                                   double q[] );

}

void integrate_ode_fwd_erk( double tin, double tout,
                            int nvar,
                            std::vector<double>& var,
                            std::vector<double>& rtol,
                            std::vector<double>& atol,
                            void (*fun) (int*, double*, double[], double[]),
                            std::vector<int>& icntrl_u,
                            std::vector<double>& rcntrl_u,
                            std::vector<int>& istatus_u,
                            std::vector<double>& rstatus_u,
                            int ierr_u );

#endif
