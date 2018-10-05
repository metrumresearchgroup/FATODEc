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
  fjac[0] = -c[0];                    // df0/dy0
  fjac[1] = c[0];                     // df1/dy0
  fjac[2] = 0.0;                      // df2/dy0
  fjac[3] = c[1]*y[2];                // df0/dy1
  fjac[4] = -c[1]*y[2] - 2*c[2]*y[1]; // df1/dy1
  fjac[5] = 2*c[2]*y[1];              // df2/dy1
  fjac[6] = c[1]*y[1];                // df0/dy2
  fjac[7] = -c[1]*y[1];               // df1/dy2
  fjac[8] = 0.0;                      // df2/dy2
}

/*
 * hv = (d^2 fy/dydy * v) * u = (H * v) * u = sum(H_ijk * v_k * u_j)
 * where
 * H_ijk = d^2 fy_i/(dy_j dy_k)
 */
static void froberts_hess(int* n, double* t, double y[], double u[], double v[], double hv[]) {
  const std::vector<double> c{0.04, 1.0E4, 3.0E7};
  const double c1 = c[0], c2 = c[1], c3 = c[2];
  const double u1 = u[0], u2 = u[1], u3 = u[2];
  const double v1 = v[0], v2 = v[1], v3 = v[2];
  hv[0] = c2 * u2 * v3 + c2 * u3 * v2;
  hv[1] = -2 * c3 * u2 * v2 - c2 * u2 * v3 - c2 * u3 * v2;
  hv[2] = 2 * c3 * u2 * v2;
}

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

struct RobertsHessianFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& u_in,
             const std::vector<double>& v_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> hv(3);
    const std::vector<double>& c(theta);
    const std::vector<double>& u(u_in);
    const std::vector<double>& v(v_in);
    const double c1 = c[0], c2 = c[1], c3 = c[2];
    const double u1 = u[0], u2 = u[1], u3 = u[2];
    const double v1 = v[0], v2 = v[1], v3 = v[2];
    hv[0] = c2 * u2 * v3 + c2 * u3 * v2;
    hv[1] = -2 * c3 * u2 * v2 - c2 * u2 * v3 - c2 * u3 * v2;
    hv[2] = 2 * c3 * u2 * v2;
    return hv;
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
  RobertsHessianFunctor fh;
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

TEST_F(FATODEBindingTest_roberts, TLM_erk) {
  nnzero = 9;
  icntrl_u[2] = 6;              // Dopri853 scheme
  icntrl_u[3] = 0;              // default nb. of steps
  icntrl_u[4] = 1;              // user-supplied func for Jacobian calc

  const std::vector<double> fy_sol {0.9055186786, 2.240475687e-05, 0.09445891666};
  const std::vector<double> y_tlm_sol {0.9203321395,
      8.286429025e-06,
      0.07965957411,
      0.5438957265,
      -3.514859954e-05,
      0.4561394221,
      0.5439342307,
      -3.514415674e-05,
      0.4561009134};

  integrate_fatode_tlm_erk( &tin, &tout, &n, &ntlm,
                            fy.data(), y_tlm.data(),
                            rtol_tlm.data(), atol_tlm.data(),
                            rtol.data(), atol.data(),
                            froberts, froberts_jac,
                            icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i)     EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
  for (int i = 0; i < n * n; ++i) EXPECT_FLOAT_EQ(y_tlm[i], y_tlm_sol[i]);

  // reset init condition
  fy = init;
  y_tlm = init_tlm;
  fatode_cc::integrate_ode_tlm_erk( tin, tout, n, ntlm, fy, y_tlm,
                                    rtol_tlm, atol_tlm, rtol, atol,
                                    f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                    theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i)     EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
  for (int i = 0; i < n * n; ++i) EXPECT_FLOAT_EQ(y_tlm[i], y_tlm_sol[i]);
}

TEST_F(FATODEBindingTest_roberts, TLM_ros) {
  nnzero = 9;
  icntrl_u[3] = 0;

  const std::vector<double> fy_sol {0.9055186772, 2.240475658e-05, 0.09445891807};
  const std::vector<double> y_tlm_sol {0.9203321368,
      8.286430121e-06 ,
      0.07965957679   ,
      0.5438957522    ,
      -3.514859061e-05,
      0.4561393964    ,
      0.5439342564    ,
      -3.514414781e-05,
      0.4561008878    };

  integrate_fatode_tlm_ros( &tin, &tout, &n, &ntlm, &nnzero,
                            fy.data(), y_tlm.data(),
                            rtol_tlm.data(), atol_tlm.data(),
                            rtol.data(), atol.data(),
                            froberts, froberts_jac, froberts_hess,
                            icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i)     EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
  for (int i = 0; i < n * n; ++i) EXPECT_FLOAT_EQ(y_tlm[i], y_tlm_sol[i]);

  // reset init condition
  fy = init;
  y_tlm = init_tlm;
  fatode_cc::integrate_ode_tlm_ros( tin, tout, n, ntlm, nnzero, fy, y_tlm,
                                    rtol_tlm, atol_tlm, rtol, atol,
                                    f, fj, fh, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                    theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i)     EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
  for (int i = 0; i < n * n; ++i) EXPECT_FLOAT_EQ(y_tlm[i], y_tlm_sol[i]);
}
