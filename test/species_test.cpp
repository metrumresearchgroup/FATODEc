#include <fatodec.hpp>
#include <fatode_cc.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

static void fspc(int* n, double* t, double y[], double fy[]) {
  // c1 = [r1, k1, b1]
  // c2 = [r2, k2, b2]
  const std::vector<double> c1{1.0, 1.0, 1.0}; // 
  const std::vector<double> c2{0.75, 0.75, 0.6666666667};
  fy[0] = c1[0] * (1.0 - y[0]/c1[1] - c1[2] * y[1]) * y[0];
  fy[1] = c2[0] * (1.0 - y[1]/c2[1] - c2[2] * y[0]) * y[1];
} 

static void fspc_jac(int* n, double* t, double y[], double fjac[]) {
  const std::vector<double> c1{1.0, 1.0, 1.0};
  const std::vector<double> c2{0.75, 0.75, 0.6666666667};
  // Column major for Fortran
  fjac[0] = -2.0 * c1[0]/c1[1] * y[0]; // df0/dy0
  fjac[1] = -c2[0] * c2[2] * y[1];     // df1/dy0
  fjac[2] = -c1[0] * c1[2] * y[0];     // df0/dy1
  fjac[3] = -2.0 * c2[0]/c2[1] * y[1]; // df1/dy1
}

static void fspc_hess(int* n, double* t, double y[], double u[], double v[], double hv[]) {
  const std::vector<double> c1{1.0, 1.0, 1.0};
  const std::vector<double> c2{0.75, 0.75, 0.6666666667};
  hv[0] = -2.0 * c1[0]/c1[1] * u[0] * v[0] - c1[0] * c1[2] * u[0] * v[1]
    - c1[0] * c1[2] * u[1] * v[0];
  hv[1] = -2.0 * c2[0]/c2[1] * u[1] * v[1] - c2[0] * c2[2] * u[0] * v[1]
    -c2[0] * c2[2] * u[1] * v[0];
}

struct SpcFunctor {
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> fy(2);
    const std::vector<double>& y(y_in);
    const std::vector<double> c1{theta[0], theta[1], theta[2]};
    const std::vector<double> c2{theta[3], theta[4], theta[5]};
    fy[0] = c1[0] * (1.0 - y[0]/c1[1] - c1[2] * y[1]) * y[0];
    fy[1] = c2[0] * (1.0 - y[1]/c2[1] - c2[2] * y[0]) * y[1];
    return fy;
  }
};

struct SpcJacobiFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> fjac(4);
    const std::vector<double> c1{theta[0], theta[1], theta[2]};
    const std::vector<double> c2{theta[3], theta[4], theta[5]};
    const std::vector<double>& y(y_in);
    fjac[0] = -2.0 * c1[0]/c1[1] * y[0]; // df0/dy0
    fjac[1] = -c2[0] * c2[2] * y[1];     // df1/dy0
    fjac[2] = -c1[0] * c1[2] * y[0];     // df0/dy1
    fjac[3] = -2.0 * c2[0]/c2[1] * y[1]; // df1/dy1
    return fjac;
  }
};

struct SpcHessianFunctor {
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& u_in,
             const std::vector<double>& v_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> hv(2);
    const std::vector<double> c1{theta[0], theta[1], theta[2]};
    const std::vector<double> c2{theta[3], theta[4], theta[5]};
    const std::vector<double>& u(u_in);
    const std::vector<double>& v(v_in);
    hv[0] = -2.0 * c1[0]/c1[1] * u[0] * v[0] - c1[0] * c1[2] * u[0] * v[1]
      - c1[0] * c1[2] * u[1] * v[0];
    hv[1] = -2.0 * c2[0]/c2[1] * u[1] * v[1] - c2[0] * c2[2] * u[0] * v[1]
      -c2[0] * c2[2] * u[1] * v[0];
    return hv;
  }
};

/*
 * single harmonic oscillator test.
 * In TLM methods the sensitivity is calclated w.r.t. the
 * initial condition.
 */
struct FATODEBindingTest_spc : public testing::Test {
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
  SpcFunctor f;
  SpcJacobiFunctor fj;
  SpcHessianFunctor fh;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
  }
  FATODEBindingTest_spc() :
    tin(0.0),
    tout(2.0),
    n(2),
    nnzero(0),
    ntlm(n),
    init{7.8, 4.5},
    fy(init),
    init_tlm{1.0, 0.0, 0.0, 1.0},
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
    theta{1.0, 1.0, 1.0, 0.75, 0.75, 0.6666666667},
    msgs(nullptr) {
      icntrl_u[3] = 10000;      // max nb. of steps
      icntrl_u[11] = 1;         // use TLM error
      rcntrl_u[2]=1.0e-3;       // init step size

      SetUp();
    }
};

TEST_F(FATODEBindingTest_spc, FWD_ERK) {
  const std::vector<double> fy_sol{0.5438183918, 0.6387706587};

  integrate_fatode_fwd_erk( &tin, &tout, &n, fy.data(), rtol.data(), atol.data(), fspc, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  std::cout.precision(10);
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_erk( tin, tout, n, fy, rtol, atol, f, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                         theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_spc, FWD_rk) {
  const std::vector<double> fy_sol{0.5438179172, 0.6387706465};
  integrate_fatode_fwd_rk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fspc, fspc_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_rk( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_spc, FWD_ros) {
  const std::vector<double> fy_sol{0.5438096963, 0.6387921035};

  integrate_fatode_fwd_ros( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fspc, fspc_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_ros( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_spc, FWD_sdirk) {
  const std::vector<double> fy_sol{0.5438180256, 0.6387717331};

  integrate_fatode_fwd_sdirk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fspc, fspc_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);

  fy = init;                    // reset init condition
  fatode_cc::integrate_ode_fwd_sdirk( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  for (int i = 0; i < n; ++i) EXPECT_FLOAT_EQ(fy[i], fy_sol[i]);
}

TEST_F(FATODEBindingTest_spc, TLM_erk) {
  nnzero = 4;
  icntrl_u[2] = 6;              // Dopri853 scheme
  icntrl_u[3] = 0;              // default nb. of steps
  icntrl_u[4] = 1;              // user-supplied func for Jacobian calc

  const std::vector<double> fy_sol {0.5438178805, 0.6387706098};
  const std::vector<double> y_tlm_sol {0.02072591983 , -0.01255852118 , -0.044087096 , 0.03100945304};

  integrate_fatode_tlm_erk( &tin, &tout, &n, &ntlm,
                            fy.data(), y_tlm.data(),
                            rtol_tlm.data(), atol_tlm.data(),
                            rtol.data(), atol.data(),
                            fspc, fspc_jac,
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

TEST_F(FATODEBindingTest_spc, TLM_ros) {
  nnzero = 4;
  icntrl_u[3] = 0;

  const std::vector<double> fy_sol {0.5438158804, 0.6387788157};
  const std::vector<double> y_tlm_sol {0.02072822873, -0.01255963875, -0.04409122583, 0.03101172573};

  integrate_fatode_tlm_ros( &tin, &tout, &n, &ntlm, &nnzero,
                            fy.data(), y_tlm.data(),
                            rtol_tlm.data(), atol_tlm.data(),
                            rtol.data(), atol.data(),
                            fspc, fspc_jac, fspc_hess,
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
