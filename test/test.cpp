#include <fatodec.hpp>
#include <fatode_cc.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>


static void fsho(int* n, double* t, double y[], double fy[]) {
  double theta = 0.15;
  fy[0] = y[1];
  fy[1] = -y[0] - theta * y[1];
} 

static void fsho_jac(int* n, double* t, double y[], double fjac[]) {
  double theta = 0.15;
  fjac[0] = 0.0;
  fjac[1] = -1.0;
  fjac[2] = 1.0;
  fjac[3] = -theta;
}

struct ShoFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> res;
    res.push_back(y_in.at(1));
    res.push_back(-y_in.at(0) - theta.at(0) * y_in.at(1));
    return res;
  }
};

struct ShoJacobiFunctor{
  inline std::vector<double>
  operator()(const double& t_in, const std::vector<double>& y_in,
             const std::vector<double>& theta, const std::vector<double>& x_r,
             const std::vector<int>& x_i, std::ostream* msgs) const {
    std::vector<double> res {0.0, -1.0, 1.0, -theta.at(0)};
    return res;
  }
};

struct FATOdeBindingTest : public testing::Test {
  double tin;
  double tout;
  int n;
  int nnzero;
  std::vector<double> fy;
  std::vector<double> rtol;
  std::vector<double> atol;
  std::vector<int> icntrl_u;
  std::vector<double> rcntrl_u;
  std::vector<int> istatus_u;
  std::vector<double> rstatus_u;
  int ierr_u;
  ShoFunctor f;
  ShoJacobiFunctor fj;
  std::vector<double> theta;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;

  void SetUp() {
  }
  FATOdeBindingTest() :
    tin(0.0),
    tout(5.52),
    n(2),
    nnzero(0),
    fy{0.0, 1.0},
    rtol{1.E-6, 1.E-6},
    atol{1.E-10, 1.E-10},
    icntrl_u(20),
    rcntrl_u(20),
    istatus_u(20),
    rstatus_u(20),
    ierr_u(0),
    theta{0.15},
    msgs(nullptr) {
      // control defaults
      icntrl_u[0] = 0;
      icntrl_u[1] = 0;
      icntrl_u[2] = 4;
      icntrl_u[3] = 1000;

      rcntrl_u[0] = 0.0;
      rcntrl_u[1] = 0.1;
      rcntrl_u[2] = 0.1;
      rcntrl_u[3] = 0.2;
      rcntrl_u[4] = 6.0;
      rcntrl_u[5] = 0.1;
      rcntrl_u[6] = 0.9;
      rcntrl_u[7] = 0.001;
      rcntrl_u[8] = 0.03;
      rcntrl_u[9] = 1.0;
      rcntrl_u[10] = 1.2;      

      SetUp();
    }
};

TEST_F(FATOdeBindingTest, FWD_ERK) {
  integrate_fatode_fwd_erk( &tin, &tout, &n, fy.data(), rtol.data(), atol.data(), fsho, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  const std::vector<double> fy_sol {-0.4655835298, 0.5054222626};
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);

  fy[0] = 0.0;                  // reset init condition
  fy[1] = 1.0;
  fatode_cc::integrate_ode_fwd_erk( tin, tout, n, fy, rtol, atol, f, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                         theta, x_r, x_i, msgs );
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);
}

TEST_F(FATOdeBindingTest, FWD_rk) {
  icntrl_u[0] = 0;
  icntrl_u[1] = 0;
  icntrl_u[2] = 4;
  icntrl_u[3] = 1000;
  icntrl_u[4] = 0;
  icntrl_u[5] = 0;
  icntrl_u[9] = 0;
  icntrl_u[10] = 0;

  integrate_fatode_fwd_rk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fsho, fsho_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  const std::vector<double> fy_sol {-0.465583533, 0.5054222670};
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);

  fy[0] = 0.0;                  // reset init condition
  fy[1] = 1.0;
  fatode_cc::integrate_ode_fwd_rk( tin, tout, n, nnzero, fy, rtol, atol, f, fj, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,
                                   theta, x_r, x_i, msgs );
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);  
}

TEST_F(FATOdeBindingTest, FWD_ros) {
  icntrl_u[0] = 0;
  icntrl_u[1] = 0;
  icntrl_u[2] = 4;
  icntrl_u[3] = 1000;
  icntrl_u[4] = 0;
  icntrl_u[5] = 0;
  icntrl_u[9] = 0;
  icntrl_u[10] = 0;

  const std::vector<double> fy_sol {-0.46558318, 0.50542152};
  integrate_fatode_fwd_ros( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fsho, fsho_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);
}

TEST_F(FATOdeBindingTest, FWD_sdirk) {
  icntrl_u[0] = 0;
  icntrl_u[1] = 0;
  icntrl_u[2] = 4;
  icntrl_u[3] = 1000;
  icntrl_u[4] = 0;
  icntrl_u[5] = 0;
  icntrl_u[9] = 0;
  icntrl_u[10] = 0;

  const std::vector<double> fy_sol {-0.46558353, 0.50542212};
  integrate_fatode_fwd_sdirk( &tin, &tout, &n, &nnzero, fy.data(), rtol.data(), atol.data(), fsho, fsho_jac, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);
}
