#include <fatodec.hpp>

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
                            int ierr_u ) {
  integrate_fatode_fwd_erk( &tin, &tout,
                            &nvar, var.data(),
                            rtol.data(), atol.data(), fun,
                            icntrl_u.data(), rcntrl_u.data(),
                            istatus_u.data(), rstatus_u.data(),
                            &ierr_u );
}
