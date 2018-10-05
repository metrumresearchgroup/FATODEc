#ifndef FATODE_INTEGRATE_ODE_TLM_ROS_HPP
#define FATODE_INTEGRATE_ODE_TLM_ROS_HPP

#include <vector>
#include <fatode_fun.hpp>
#include <fatode_ode.hpp>
#include <fatodec.hpp>

namespace fatode_cc {

  template<typename F, typename Fj, typename Fh>
  void integrate_ode_tlm_ros( double tin, double tout,
                              int n,
                              int ntlm,
                              int nnzero,
                              std::vector<double>& var,
                              std::vector<double>& var_tlm,
                              std::vector<double>& rtol_tlm,
                              std::vector<double>& atol_tlm,
                              std::vector<double>& rtol,
                              std::vector<double>& atol,
                              F& f,
                              Fj& fj,
                              Fh& fh,
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
    FatOde<Fj> odejac(fj, theta, x_r, x_i, msgs);
    FatOde<Fh> odehess(fh, theta, x_r, x_i, msgs);
    void* fy_user_data = static_cast<void*>(&ode);
    void* fjac_user_data = static_cast<void*>(&odejac);
    void* fhess_user_data = static_cast<void*>(&odehess);

    integrate_fatode_tlm_ros_cc( &tin, &tout, &n, &ntlm, &nnzero,
                                 var.data(), var_tlm.data(),
                                 rtol_tlm.data(), atol_tlm.data(),
                                 rtol.data(), atol.data(),
                                 fatode_fun<FatOde<F> >(),
                                 fatode_jac<FatOde<Fj> >(),
                                 fatode_hess<FatOde<Fh> >(),
                                 icntrl_u.data(), rcntrl_u.data(),
                                 istatus_u.data(), rstatus_u.data(),
                                 &ierr_u,
                                 fy_user_data,
                                 fjac_user_data,
                                 fhess_user_data);
  }
}

#endif
