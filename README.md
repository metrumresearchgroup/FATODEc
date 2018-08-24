
# Table of Contents

1.  [FATODEc](#org32123a0)
2.  [Build](#orgc11d2eb)
3.  [Quick start](#orge37a877)
4.  [Test](#orgf454046)
5.  [Use](#orgdccc08d)
    1.  [Forward integration solvers(FWD)](#orgc4451a2)
    2.  [Tangential linear model(TLM)](#org2e2d078)
6.  [TO-DO](#orga83fc5c)
    1.  [Tangential linear model(TLM)](#org0644f6e)
    2.  [Adjoint sensitivity solvers(ADJ)](#org230c3ef)


<a id="org32123a0"></a>

# FATODEc

C/C++ bindings of [FATODE](http://people.cs.vt.edu/asandu/Software/FATODE/index.html) solver library. Currently doesn't
support sparse matrix solvers, as the application is mainly
for PKPD modeling.


<a id="orgc11d2eb"></a>

# Build

To make `libfatode.a`, in `FATODEc` run

    make all


<a id="orge37a877"></a>

# Quick start

One can find the following example in `example/sho`. Consult
[FATODE user guide](http://people.cs.vt.edu/%7Easandu/Software/FATODE/FATODE_user_guide.pdf) for details.

    #include <stdio.h>
    #include <stdlib.h>
    #include <math.h>
    
    #include <fatodec.h>
    
    static void f(int*, double*, double[], double[]);
    
    int main(int argc, char *argv[])
    {
    
      double tin = 0.0;
      double tout = 10.0;
      int nvar = 2;
      double   fy[2] = { 0.0, 2.0 };
      double rtol[2] = { 1.E-6, 1.E-6 };
      double atol[2] = { 1.E-10, 1.E-10 };
      int icntrl_u[20];
      double rcntrl_u[20];
      int istatus_u[20];
      double rstatus_u[20];
      int ierr_u;
    
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
    
      integrate_fatode_fwd_erk( &tin, &tout, &nvar, fy, rtol, atol, f, icntrl_u, rcntrl_u, istatus_u, rstatus_u, &ierr_u );    
    
      printf("No. of function calls                          :  %d\n", istatus_u[0]);
      printf("No. of jacobian calls                          :  %d\n", istatus_u[1]);
      printf("No. of steps                                   :  %d\n", istatus_u[2]);
      printf("No. of accepted steps                          :  %d\n", istatus_u[3]);
      printf("No. of rejected steps (except at the beginning :  %d\n", istatus_u[4]);
      printf("No. of LU decompositions                       :  %d\n", istatus_u[5]);
      printf("No. of forward/backward substitutions          :  %d\n", istatus_u[6]);
      printf("No. of singular matrix decompositions          :  %d\n", istatus_u[7]);
    
      printf("T at exit                                      :  %f\n", rstatus_u[0]);
      printf("H at exit                                      :  %f\n", rstatus_u[1]);
      printf("H new                                          :  %f\n", rstatus_u[2]);
    
      return 0;
    }
    
    static void f(int* n, double* t, double y[], double fy[]) {
      double theta = 0.15;
      fy[0] = y[1];
      fy[1] = -y[0] - theta * y[1];
    }

Sample output

    No. of function calls                          :  728
    No. of jacobian calls                          :  0
    No. of steps                                   :  104
    No. of accepted steps                          :  103
    No. of rejected steps (except at the beginning :  0
    No. of LU decompositions                       :  0
    No. of forward/backward substitutions          :  0
    No. of singular matrix decompositions          :  0
    T at exit                                      :  10.000000
    H at exit                                      :  0.020000
    H new                                          :  0.100000


<a id="orgf454046"></a>

# Test

In `test`

    make test; ./test


<a id="orgdccc08d"></a>

# Use

Use C header `FATODEc/include/fatodec.h` or C++ header `FATODEc/include/fatodec.hpp`.


<a id="orgc4451a2"></a>

## Forward integration solvers(FWD)

-   Explicit Runge-Kutta `integrate_fatode_fwd_erk`
-   Fully Implicit Runge-Kutta methods `integrate_fatode_fwd_rk`
-   Rosenbrock methods `integrate_fatode_fwd_ros`
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_fwd_sdirk`


<a id="org2e2d078"></a>

## Tangential linear model(TLM)

-   Explicit Runge-Kutta `integrate_fatode_tlm_erk`
-   Rosenbrock methods `integrate_fatode_tlm_ros`
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_tlm_sdirk`


<a id="orga83fc5c"></a>

# TO-DO


<a id="org0644f6e"></a>

## Tangential linear model(TLM)

-   Fully Implicit Runge-Kutta methods `integrate_fatode_tlm_rk`


<a id="org230c3ef"></a>

## Adjoint sensitivity solvers(ADJ)

-   Explicit Runge-Kutta `integrate_fatode_adj_erk`
-   Fully Implicit Runge-Kutta methods `integrate_fatode_adj_rk`
-   Rosenbrock methods `integrate_fatode_adj_ros`
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_adj_sdirk`

