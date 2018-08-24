#include <fatodec.hpp>
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

struct FATOdeBindingTest : public testing::Test {
  double tin;
  double tout;
  int n;
  std::vector<double> fy;
  std::vector<double> rtol;
  std::vector<double> atol;
  std::vector<int> icntrl_u;
  std::vector<double> rcntrl_u;
  std::vector<int> istatus_u;
  std::vector<double> rstatus_u;
  int ierr_u;

  void SetUp() {
  }
  FATOdeBindingTest() :
    tin(0.0),
    tout(10.0),
    n(2),
    fy{0.0, 1.0},
    rtol{1.E-6, 1.E-6},
    atol{1.E-10, 1.E-10},
    icntrl_u(20),
    rcntrl_u(20),
    istatus_u(20),
    rstatus_u(20),
    ierr_u(0) {
      // control defaults
      icntrl_u[0] = 0;
      icntrl_u[1] = 0;
      icntrl_u[2] = 4;
      icntrl_u[3] = 1000;

      rcntrl_u[0] = 1.E-4;
      rcntrl_u[1] = 0.1;
      rcntrl_u[2] = 0.1;
      rcntrl_u[3] = 0.2;
      rcntrl_u[4] = 6.0;
      rcntrl_u[5] = 0.1;
      rcntrl_u[6] = 0.9;
      rcntrl_u[9] = 1.0;
      rcntrl_u[10] = 1.2;      

      SetUp();
    }
};

TEST_F(FATOdeBindingTest, FWD_ERK) {
  integrate_fatode_fwd_erk( &tin, &tout, &n, fy.data(), rtol.data(), atol.data(), fsho, icntrl_u.data(), rcntrl_u.data(), istatus_u.data(), rstatus_u.data(), &ierr_u );
  std::cout.precision(17);
  const std::vector<double> fy_sol {-0.2464078806, -0.3849482586};
  EXPECT_FLOAT_EQ(fy[0], fy_sol[0]);
  EXPECT_FLOAT_EQ(fy[1], fy_sol[1]);
}
