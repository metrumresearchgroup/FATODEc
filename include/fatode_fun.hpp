#ifndef FATODE_CC_FUN_HPP
#define FATODE_CC_FUN_HPP

#include <iostream>

/* convert a C++ styple functor to FATODE's RHS function pointer.
 * @tparam Ode ODE system that contains RHS functor with a stan functor signature 
 *           vector = f(t, y, theta, x_r, x_i, std::ostream*)
 *           and @c theta, @c x_r, @c x_i
 */
namespace fatode_cc {

  /*
   * return a lambda that can be passed to Fortran function
   * @c xxx_cc to calculate RHS vector.
   */
  template<typename Ode>
  auto fatode_fun() {
    return [](int* n, double* t, double y[], double fy[],
              void* user_data) {
      Ode* ode = static_cast<Ode*>(user_data);
      std::vector<double> yv(y, y + *n);
      std::vector<double> fyv(*n);
      fyv = ode -> f(*t, yv, ode -> theta, ode -> x_r, ode -> x_i, ode -> msgs);
      for (int i = 0; i < *n; ++i) {
        fy[i] = fyv[i];
      }
    };
  }
  
  /*
   * return a lambda that can be passed to Fortran function
   * @c xxx_cc to calculate Jacobi matrix.
   */
  template<typename OdeJac>
  auto fatode_jac() {
    return [](int* n, double* t, double y[], double fjac[],
              void* user_data) {
      OdeJac* jac = static_cast<OdeJac*>(user_data);
      std::vector<double> yv(y, y + *n);
      std::vector<double> fjacv((*n) * (*n));
      fjacv = jac -> f(*t, yv, jac -> theta, jac -> x_r, jac -> x_i, jac -> msgs);
      for (int i = 0; i < (*n) * (*n); ++i) {
        fjac[i] = fjacv[i];
      }
    };
  }

}


#endif
