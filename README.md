- [FATODEc](#org3b282bc)
- [Build](#org202635a)
- [Quick start](#org2430e5f)
- [Test](#orgf855a13)
- [Use](#orgdd4cdbd)
  - [Forward integration solvers(FWD)](#org4b3f200)
  - [Tangential linear model(TLM)](#org557773b)
- [TO-DO](#org0ec929b)
  - [Tangential linear model(TLM)](#org74fd95c)
  - [Adjoint sensitivity solvers(ADJ)](#orgce40f87)


<a id="org3b282bc"></a>

# FATODEc

C/C++ bindings of [FATODE](http://people.cs.vt.edu/asandu/Software/FATODE/index.html) solver library. Currently doesn't support sparse matrix solvers, as the application is mainly for PKPD modeling.


<a id="org202635a"></a>

# Build

To make `libfatode.a`, in `FATODEc` run

```bash
make all
```


<a id="org2430e5f"></a>

# Quick start

One can find the following example in `example/sho`. Consult [FATODE user guide](http://people.cs.vt.edu/%7Easandu/Software/FATODE/FATODE_user_guide.pdf) for details.

```c
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
```

Sample output

```text
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
```


<a id="orgf855a13"></a>

# Test

In `test`

```
make test; ./test
```


<a id="orgdd4cdbd"></a>

# Use

Use C header `FATODEc/include/fatodec.h` or C++ header `FATODEc/include/fatode_cc.hpp`. The C++ interface is in namespace `fatode_cc`. Currently neither SUPER-LU or UMFPACK linear solver is supported.


<a id="org4b3f200"></a>

## Forward integration solvers(FWD)

-   Explicit Runge-Kutta `integrate_fatode_fwd_erk` (C & C++)
-   Fully Implicit Runge-Kutta methods `integrate_fatode_fwd_rk` (C & C++)
-   Rosenbrock methods `integrate_fatode_fwd_ros` (C & C++)
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_fwd_sdirk` (C & C++)


<a id="org557773b"></a>

## Tangential linear model(TLM)

-   Explicit Runge-Kutta `integrate_fatode_tlm_erk` (C & C++)
-   Rosenbrock methods `integrate_fatode_tlm_ros` (C & C++)


<a id="org0ec929b"></a>

# TO-DO


<a id="org74fd95c"></a>

## Tangential linear model(TLM)

-   Fully Implicit Runge-Kutta methods `integrate_fatode_tlm_rk`
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_tlm_sdirk`


<a id="orgce40f87"></a>

## Adjoint sensitivity solvers(ADJ)

-   Explicit Runge-Kutta `integrate_fatode_adj_erk`
-   Fully Implicit Runge-Kutta methods `integrate_fatode_adj_rk`
-   Rosenbrock methods `integrate_fatode_adj_ros`
-   Singly Diagonally Implicit Runge-Kutta methods `integrate_fatode_adj_sdirk`
