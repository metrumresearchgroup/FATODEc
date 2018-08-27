#ifndef FATODE_INTEGRATE_ODE_FWD_ERK_HPP
#define FATODE_INTEGRATE_ODE_FWD_ERK_HPP

#include <vector>
#include <fatode_fun.hpp>
#include <fatode_ode.hpp>
#include <fatodec.hpp>

namespace fatode_cc {

  template<typename F>
  void integrate_ode_fwd_erk( double tin, double tout,
                              int nvar,
                              std::vector<double>& var,
                              std::vector<double>& rtol,
                              std::vector<double>& atol,
                              F& f,
                              std::vector<int>& icntrl_u,
                              std::vector<double>& rcntrl_u,
                              std::vector<int>& istatus_u,
                              std::vector<double>& rstatus_u,
                              int ierr_u,
                              std::vector<double>& theta,
                              std::vector<double>& x_r,
                              std::vector<int>& x_i,
                              std::ostream* msgs ) {
    FatOde<F> ode(f, theta, x_r, x_i, msgs);
    void* user_data = static_cast<void*>(&ode);

    integrate_fatode_fwd_erk_cc( &tin, &tout,
                                 &nvar, var.data(),
                                 rtol.data(), atol.data(),
                                 fatode_fun<FatOde<F> >(),
                                 icntrl_u.data(), rcntrl_u.data(),
                                 istatus_u.data(), rstatus_u.data(),
                                 &ierr_u,
                                 user_data );
    
  }
}

#endif
