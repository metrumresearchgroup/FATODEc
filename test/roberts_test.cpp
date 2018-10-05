#include <fatodec.hpp>
#include <fatode_cc.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

static void froberts(int* n, double* t, double y[], double fy[]) {
  const std::vector<double> c{0.04, 1.0E4, 3.0E7};
  fy[0] = -c[0]*y[0] +c[1]*y[1]*y[2];
  fy[2] = c[2]*y[1]*y[1];
  fy[1] = -fy[0] - fy[2];
} 

static void froberts_jac(int* n, double* t, double y[], double fjac[]) {
  const std::vector<double> c{0.04, 1.0E4, 3.0E7};
  // Column major for Fortran
  fjac[0] = -c[0];              // df0/dy0
  fjac[1] = c[0];               // df1/dy0
  fjac[2] = 0.0;                // df2/dy0
  fjac[3] = c[1]*y[2];         // df0/dy1
  fjac[4] = -c[1]*y[2] - 2*c[2]*y[1]; // df1/dy1
  fjac[5] = 2*c[2]*y[1];               // df2/dy1
  fjac[6] = c[1]*y[1];                  // df0/dy2
  fjac[7] = -c[1]*y[1];                 // df1/dy2
  fjac[8] = 0.0;                        // df2/dy2
}

// static void froberts_hess(int* n, double* t, double y[], double u[], double v[], double hv[]) {
//   hv[0] = 0.0;
//   hv[1] = 0.0;
//   hv[2] = 0.0;
//   hv[3] = 0.0;
// }

struct RobertsFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> fy(3);
    const std::vector<double>& c(theta);
    const std::vector<double>& y(y_in);
    fy[0] = -c[0]*y[0] +c[1]*y[1]*y[2];
    fy[2] = c[2]*y[1]*y[1];
    fy[1] = -fy[0] - fy[2];
    return fy;
  }
};

struct RobertsJacobiFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> fjac(9);
    const std::vector<double>& c(theta);
    const std::vector<double>& y(y_in);
    fjac[0] = -c[0];                    // df0/dy0
    fjac[1] = c[0];                     // df1/dy0
    fjac[2] = 0.0;                      // df2/dy0
    fjac[3] = c[1]*y[2];                // df0/dy1
    fjac[4] = -c[1]*y[2] - 2*c[2]*y[1]; // df1/dy1
    fjac[5] = 2*c[2]*y[1];              // df2/dy1
    fjac[6] = c[1]*y[1];                // df0/dy2
    fjac[7] = -c[1]*y[1];               // df1/dy2
    fjac[8] = 0.0;                      // df2/dy2
    return fjac;
  }
};

/*
 * single harmonic oscillator test.
 * In TLM methods the sensitivity is calclated w.r.t. the
 * initial condition.
 */
struct FATODEBindingTest_roberts : public testing::Test {
  double tin;
  double tout;
  int n;
  int nnzero;
  int ntlm;
  std::vector<double> init;
  std::vector<double> fy;
  std::vector<double> init_tlm;
  std::vector<double> y_tlm;
  std::vector<double> rtol_tlm;
  std::vector<double> atol_tlm;
  std::vector<double> rtol;
  std::vector<double> atol;
  std::vector<int> icntrl_u;
  std::vector<double> rcntrl_u;
  std::vector<int> istatus_u;
  std::vector<double> rstatus_u;
  int ierr_u;
  RobertsFunctor f;
  RobertsJacobiFunctor fj;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
  }
  FATODEBindingTest_roberts() :
    tin(0.0),
    tout(4.0),
    n(3),
    nnzero(0),
    ntlm(n),
    init{1.0, 0.0, 0.0},
    fy(init),
    init_tlm{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0},
    y_tlm(init_tlm),
    rtol_tlm(n * ntlm, 1.E-6),
    atol_tlm(n * ntlm, 1.E-10),
    rtol{1.E-6, 1.E-6, 1.E-6},
    atol{1.E-8, 1.E-14, 1.E-6},
    icntrl_u(20, 0),
    rcntrl_u(20, 0.0),
    istatus_u(20),
    rstatus_u(20),
    ierr_u(0),
    theta{0.04, 1.0E4, 3.0E7},
    msgs(nullptr) {
      icntrl_u[3] = 100000;     // max nb. of steps
      icntrl_u[11] = 1;         // use TLM error
      rcntrl_u[2]=1.0e-3;       // init step size

      SetUp();
    }
};

TEST_F(FATODEBindingTest_roberts, FWD_ERK) {
  const std::vector<double> fy_sol{0.9055186786 , 2.240474479e-05 , 0.09445891667};

  integrate_fatode_fwd_erk( &tin, &tout, &n, fy.data(), rtol.data(), atol.data(), froberts, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_erk( tin, tout, n, fy, rtol, atol, f, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                         theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_roberts, FWD_rk) {
  nnzero = 9;

  const std::vector<double> fy_sol{0.9055186786 , 2.240475769e-05 , 0.09445891667};
  integrate_fatode_fwd_rk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), froberts, froberts_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_rk( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_roberts, FWD_ros) {
  nnzero = 9;

  const std::vector<double> fy_sol{0.9055186786 , 2.240466623e-05 , 0.09445891667};
  integrate_fatode_fwd_ros( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), froberts, froberts_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  // std::cout.precision(10);
  // std::cout << "taki test: " << fy[1] << "\n";

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_ros( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_roberts, FWD_sdirk) {
  nnzero = 9;

  const std::vector<double> fy_sol{0.9055186786 , 2.240476437e-05 , 0.09445891667};
  integrate_fatode_fwd_sdirk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), froberts, froberts_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_sdirk( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}
