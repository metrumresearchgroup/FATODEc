#ifndef FATODE_ODE_HPP
#define FATODE_ODE_HPP

#include <vector>

namespace fatode_cc {

template<typename F>
struct FatOde
{
  const F& f;
  const std::vector<double>& theta;
  const std::vector<double>& x_r;
  const std::vector<int>& x_i;
  std::ostream* msgs;

  FatOde( const F& f_in,
          const std::vector<double>& theta_in,
          const std::vector<double>& x_r_in,
          const std::vector<int>& x_i_in,
          std::ostream* msgs_in) :
    f(f_in), theta(theta_in), x_r(x_r_in), x_i(x_i_in), msgs(msgs_in)
  {}
};

}

#endif
