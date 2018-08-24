
# Table of Contents

1.  [FATODEc](#orgf5ce894)
2.  [Build](#org886ad74)
3.  [Quick start](#org9da5f64)
4.  [Test](#org32badbf)
5.  [Use](#orgdd80e94)
    1.  [Forward integration solvers](#orgba409f3)
        1.  [Explicit Runge-Kutta](#org18af9a5)
        2.  [Fully Implicit Runge-Kutta methods](#org4d83f24)
        3.  [Rosenbrock methods](#orge3b84a6)
        4.  [Singly Diagonally Implicit Runge-Kutta methods](#orgbda3931)


<a id="orgf5ce894"></a>

# FATODEc

C/C++ bindings of [FATODE](http://people.cs.vt.edu/asandu/Software/FATODE/index.html) solver library.


<a id="org886ad74"></a>

# Build

To make `libfatode.a`, in `FATODEc` run

    make all


<a id="org9da5f64"></a>

# Quick start

One can find the following example in `example/sho`.

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


<a id="org32badbf"></a>

# Test

To make `test/test`, in `FATODEc/test` run

    make test

To run `test`, in `FATODEc/test` run

    ./test


<a id="orgdd80e94"></a>

# Use

Use C header `FATODEc/include/fatodec.h` or C++ header `FATODEc/include/fatodec.hpp`.


<a id="orgba409f3"></a>

## Forward integration solvers


<a id="org18af9a5"></a>

### Explicit Runge-Kutta

    void integrate_fatode_fwd_erk( double* tin, double* tout,
                                   int* nvar, double var[],
                                   double rtol[], double atol[],
                                   void (*fun) (int*, double*, double[], double[]),
                                   int icntrl_u[],
                                   double rcntrl_u[],
                                   int istatus_u[],
                                   double rstatus_u[],
                                   int* ierr_u );


<a id="org4d83f24"></a>

### Fully Implicit Runge-Kutta methods

    void integrate_fatode_fwd_rk( double* tin, double* tout,
                                  int* nvar, int* nnzero, double var[],
                                  double rtol[], double atol[],
                                  void (*fun) (int*, double*, double[], double[]),
                                  void (*jac) (int*, double*, double[], double[]),
                                  int icntrl_u[],
                                  double rcntrl_u[],
                                  int istatus_u[],
                                  double rstatus_u[],
                                  int* ierr_u );


<a id="orge3b84a6"></a>

### Rosenbrock methods

    void integrate_fatode_fwd_ros( double* tin, double* tout,
                                   int* nvar, int* nnzero, double var[],
                                   double rtol[], double atol[],
                                   void (*fun) (int*, double*, double[], double[]),
                                   void (*jac) (int*, double*, double[], double[]),
                                   int icntrl_u[],
                                   double rcntrl_u[],
                                   int istatus_u[],
                                   double rstatus_u[],
                                   int* ierr_u );


<a id="orgbda3931"></a>

### Singly Diagonally Implicit Runge-Kutta methods

    void integrate_fatode_fwd_sdirk( double* tin, double* tout,
                                     int* nvar, int* nnzero, double var[],
                                     double rtol[], double atol[],
                                     void (*fun) (int*, double*, double[], double[]),
                                     void (*jac) (int*, double*, double[], double[]),
                                     int icntrl_u[],
                                     double rcntrl_u[],
                                     int istatus_u[],
                                     double rstatus_u[],
                                     int* ierr_u );

