! torsten's interface to fatode library
module fatode_c

  use iso_c_binding

  implicit none

  integer, parameter:: dp=kind(0.d0)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd erk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_erk( tin, tout, nvar, var, rtol, atol, fun, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use erk_f90_integrator
  implicit none

  integer, intent(in) :: nvar
  real(dp), intent(inout) :: var(nvar)
  real(dp), intent(in) :: rtol(nvar)
  real(dp), intent(in) :: atol(nvar)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time

  ! control and output
  integer,  intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,  intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,  intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  call integrate( tin, tout, nvar, var, rtol, atol, fun, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_fwd_erk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd erk module wrapper for C++
!
! func_cc has signature
! void func_cc(int* n, double* t, double y[], double fy[], 
!              double theta[], double x_r[], int x_i[] )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_erk_cc( tin, tout, nvar, var, rtol, atol, fun_cc, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u, &
     & user_data ) bind(c)

  use iso_c_binding
  use erk_f90_integrator
  implicit none

  integer, intent(in) :: nvar
  real(dp), intent(inout) :: var(nvar)
  real(dp), intent(in) :: rtol(nvar)
  real(dp), intent(in) :: atol(nvar)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time
  
  ! control and output
  integer,  intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,  intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,  intent(out) :: ierr_u

  ! c++ functor inputs
  type(c_ptr), intent(inout) :: user_data;

  interface
     subroutine func_cc(n, t, y, fy, user_data) bind (c)
       use iso_c_binding
       integer, parameter:: dp=kind(0.d0)

       integer,  intent(in) :: n
       real(dp), intent(in) :: t, y(n)
       real(dp), intent(inout) :: fy(n)
       type(c_ptr), intent(inout) :: user_data;
     end subroutine func_cc
  end interface
  procedure(func_cc) :: fun_cc

  call integrate( tin, tout, nvar, var, rtol, atol, fun, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

contains
  subroutine fun(n, t, y, fy) bind (c)
    use iso_c_binding
    integer, intent(in) :: n
    double precision, intent(in) :: t, y(n)
    double precision, intent(inout) :: fy(n)

    call fun_cc(n, t, y, fy, user_data)

  end subroutine fun

end subroutine integrate_fatode_fwd_erk_cc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd rk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_rk( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use rk_f90_integrator
  implicit none

  integer, intent(in) :: n, nnzero
  real(dp), intent(inout) :: var(n)
  real(dp), intent(in) :: rtol(n)
  real(dp), intent(in) :: atol(n)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time

  ! control and output
  integer,  intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,  intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,  intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in)    :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  call integrate( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_fwd_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd rk module wrapper for C++
!
! func_cc has signature
! void func_cc(int* n, double* t, double y[], double fy[],
!              double theta[], double x_r[], int x_i[] )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_rk_cc( tin, tout, n, nnzero, var, rtol, atol, fun_cc, jac_cc, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u,&
     & fy_user_data, fjac_user_data ) bind(c)

  use iso_c_binding
  use rk_f90_integrator
  implicit none

  integer, intent(in) :: n, nnzero
  real(dp), intent(inout) :: var(n)
  real(dp), intent(in) :: rtol(n)
  real(dp), intent(in) :: atol(n)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time

  ! control and output
  integer,  intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,  intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,  intent(out) :: ierr_u

  ! c++ functor inputs
  type(c_ptr), intent(inout) :: fy_user_data;
  type(c_ptr), intent(inout) :: fjac_user_data;

  interface
     subroutine func_cc(n, t, y, fy, user_data) bind (c)
       use iso_c_binding
       integer, parameter:: dp=kind(0.d0)

       integer,  intent(in) :: n
       real(dp), intent(in) :: t, y(n)
       real(dp), intent(inout) :: fy(n)
       type(c_ptr), intent(inout) :: user_data;
     end subroutine func_cc
  end interface

  interface
     subroutine fjac_cc(n, t, y, fjac, user_data) bind (c)
       use iso_c_binding
       integer, parameter:: dp=kind(0.d0)

       integer,  intent(in) :: n
       real(dp), intent(in) :: t, y(n)
       real(dp), intent(inout) :: fjac(n, n)
       type(c_ptr), intent(inout) :: user_data;
     end subroutine fjac_cc
  end interface

  procedure(func_cc) :: fun_cc
  procedure(fjac_cc) :: jac_cc

  call integrate( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

contains
  subroutine fun(n, t, y, fy) bind (c)
    use iso_c_binding
    integer, intent(in) :: n
    double precision, intent(in) :: t, y(n)
    double precision, intent(inout) :: fy(n)
    call fun_cc(n, t, y, fy, fy_user_data)
  end subroutine fun

  subroutine jac(n, t, y, fjac) bind (c)
    use iso_c_binding
    integer, intent(in) :: n
    double precision, intent(in) :: t, y(n)
    double precision, intent(inout) :: fjac(n, n)
    call jac_cc(n, t, y, fjac, fjac_user_data)
  end subroutine jac

end subroutine integrate_fatode_fwd_rk_cc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd ros module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_ros( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use ros_f90_integrator
  implicit none

  integer, intent(in) :: n, nnzero
  real(dp), intent(inout) :: var(n)
  real(dp), intent(in) :: rtol(n)
  real(dp), intent(in) :: atol(n)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time

  ! control and output
  integer,          intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,          intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,          intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  call integrate( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_fwd_ros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fwd sdirk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_fwd_sdirk( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
     icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use sdirk_f90_integrator
  implicit none

  integer, intent(in) :: n, nnzero
  real(dp), intent(inout) :: var(n)
  real(dp), intent(in) :: rtol(n)
  real(dp), intent(in) :: atol(n)
  real(dp), intent(in) :: tin  ! start time
  real(dp), intent(in) :: tout ! end time

  ! control and output
  integer,          intent(in)  :: icntrl_u(20)
  real(dp), intent(in)  :: rcntrl_u(20)
  integer,          intent(out) :: istatus_u(20)
  real(dp), intent(out) :: rstatus_u(20)
  integer,          intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  call integrate( tin, tout, n, nnzero, var, rtol, atol, fun, jac, &
       icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_fwd_sdirk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tlm erk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_tlm_erk( tin, tout, nvar, ntlm, var, var_tlm, atol_tlm, &
     rtol_tlm, atol, rtol,&
     fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use erk_tlm_f90_integrator
  implicit none

   integer, intent(in) :: nvar,ntlm
   real(dp), intent(inout) :: var(nvar),var_tlm(nvar,ntlm),&
                           rtol_tlm(nvar,ntlm), atol_tlm(nvar,ntlm)
   real(dp), intent(in) :: rtol(nvar)
   real(dp), intent(in) :: atol(nvar)
   real(dp), intent(in) :: tin  ! start time
   real(dp), intent(in) :: tout ! end time

   integer,          intent(in)  :: icntrl_u(20)
   real(dp), intent(in)  :: rcntrl_u(20)
   integer,          intent(out) :: istatus_u(20)
   real(dp), intent(out) :: rstatus_u(20)
   integer,          intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  call integrate_tlm(tin, tout, nvar, ntlm, var, var_tlm, atol_tlm, &
       rtol_tlm, atol, rtol,&
       fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_tlm_erk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tlm rk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine integrate_fatode_tlm_rk( tin, tout, n, ntlm, nnzero, y, y_tlm, atol_tlm, rtol_tlm, &
!  atol, rtol, fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

!   use iso_c_binding
!   use rk_tlm_f90_integrator
!   implicit none

!   integer :: n, ntlm, nnzero
!    real(dp), intent(inout) :: y(n), y_tlm(n,ntlm),&
!                            rtol_tlm(n,ntlm), atol_tlm(n,ntlm)
!    real(dp), intent(in) :: rtol(n)
!    real(dp), intent(in) :: atol(n)
!    real(dp), intent(in) :: tin  ! start time
!    real(dp), intent(in) :: tout ! end time

!    integer,          intent(in)  :: icntrl_u(20)
!    real(dp), intent(in)  :: rcntrl_u(20)
!    integer,          intent(out) :: istatus_u(20)
!    real(dp), intent(out) :: rstatus_u(20)
!    integer,          intent(out) :: ierr_u

!   interface
!      subroutine fun(n, t, y, fy) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fy(n)
!      end subroutine fun
!   end interface

!   interface
!      subroutine jac(n, t, y, fjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fjac(n, n)
!      end subroutine jac
!   end interface

!   call integrate_tlm( n, ntlm, nnzero, y, y_tlm, tin, tout, atol_tlm, rtol_tlm, &
!        atol, rtol, fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

! end subroutine integrate_fatode_tlm_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tlm ros module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_tlm_ros( tin, tout, n, ntlm, nnzero, y, y_tlm, atol_tlm, rtol_tlm, &
     atol, rtol, fun, jac, hess_vec, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use ros_tlm_f90_integrator
  implicit none

  integer :: n, ntlm, nnzero
   real(dp), intent(inout) :: y(n), y_tlm(n,ntlm),&
                           rtol_tlm(n,ntlm), atol_tlm(n,ntlm)
   real(dp), intent(in) :: rtol(n)
   real(dp), intent(in) :: atol(n)
   real(dp), intent(in) :: tin  ! start time
   real(dp), intent(in) :: tout ! end time

   integer,          intent(in)  :: icntrl_u(20)
   real(dp), intent(in)  :: rcntrl_u(20)
   integer,          intent(out) :: istatus_u(20)
   real(dp), intent(out) :: rstatus_u(20)
   integer,          intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  interface
     subroutine hess_vec(n,t,y,u,v,hv) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in)    :: t,y(n),u(n),v(n)
       double precision, intent(inout) :: hv(n)
     end subroutine hess_vec
  end interface

  call integrate_tlm( n, ntlm, nnzero, y, y_tlm, tin, tout, atol_tlm, rtol_tlm, &
       atol, rtol, fun, jac, hess_vec, icntrl_u, rcntrl_u, istatus_u, rstatus_u )

end subroutine integrate_fatode_tlm_ros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! tlm sdirk module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine integrate_fatode_tlm_sdirk( tin, tout, n, ntlm, nnzero, y, y_tlm, atol_tlm, rtol_tlm, &
 atol, rtol, fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u ) bind(c)

  use iso_c_binding
  use sdirk_tlm_f90_integrator
  implicit none

  integer :: n, ntlm, nnzero
   real(dp), intent(inout) :: y(n), y_tlm(n,ntlm),&
                           rtol_tlm(n,ntlm), atol_tlm(n,ntlm)
   real(dp), intent(in) :: rtol(n)
   real(dp), intent(in) :: atol(n)
   real(dp), intent(in) :: tin  ! start time
   real(dp), intent(in) :: tout ! end time

   integer,          intent(in)  :: icntrl_u(20)
   real(dp), intent(in)  :: rcntrl_u(20)
   integer,          intent(out) :: istatus_u(20)
   real(dp), intent(out) :: rstatus_u(20)
   integer,          intent(out) :: ierr_u

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  interface
     subroutine jac(n, t, y, fjac) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fjac(n, n)
     end subroutine jac
  end interface

  call integrate_tlm( n, ntlm, nnzero, y, y_tlm, tin, tout, atol_tlm, rtol_tlm, &
       atol, rtol, fun, jac, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u )

end subroutine integrate_fatode_tlm_sdirk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! adj erk module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine integrate_fatode_adj_erk( tin, tout, nvar, np, nadj, nnz, y, lambda, &
!      atol, rtol, fun, jac, adjinit,  icntrl_u, rcntrl_u, istatus_u, &
!      rstatus_u, ierr_u, mu, jacp, drdy, drdp, qfun, q ) bind(c)

!   use iso_c_binding
!   use erk_adj_f90_integrator
!   implicit none

!    integer, intent(in) :: nvar,np,nadj
!    integer, intent(in) :: nnz

!    real(dp), intent(inout) :: y(nvar)

!    real(dp), intent(inout) :: q(nadj)
!    real(dp), intent(inout)  :: lambda(nvar,nadj)
!    real(dp), intent(inout) ::  mu(np, nadj)
!    real(dp), intent(in)  ::  atol(nvar), rtol(nvar)
!    real(dp), intent(in) :: tin  ! start time
!    real(dp), intent(in) :: tout ! end time

!    integer,          intent(in)  :: icntrl_u(20)
!    real(dp), intent(in)  :: rcntrl_u(20)
!    integer,          intent(out) :: istatus_u(20)
!    real(dp), intent(out) :: rstatus_u(20)
!    integer,          intent(out) :: ierr_u

!    double precision :: rcntrl(20), rstatus(20), t1, t2
!    integer       :: icntrl(20), istatus(20), ierr

!   interface
!      subroutine fun(n, t, y, fy) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fy(n)
!      end subroutine fun
!   end interface

!   interface
!      subroutine jac(n, t, y, fjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fjac(n, n)
!      end subroutine jac
!   end interface

!   interface
!      subroutine adjinit(n,np,nadj,t,y,lambda,mu) bind (c)
!        use iso_c_binding
!        integer,intent(in) :: n,np,nadj
!        double precision,intent(in) :: t,y(n),lambda(n,nadj)
!        double precision, optional,intent(inout) :: mu(np,nadj)
!      end subroutine adjinit
!   end interface

!   interface
!      subroutine drdp(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine drdp
!   end interface
  
!   interface
!      subroutine drdy(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine drdy
!   end interface

!   interface
!      subroutine jacp(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine jacp
!   end interface

!   interface
!      subroutine qfun(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine qfun
!   end interface

!   call integrate_adj( nvar, np, nadj, nnz, y, lambda, tin, tout, &
!        atol, rtol, fun, jac, adjinit,  icntrl_u, rcntrl_u, istatus_u, &
!        rstatus_u, ierr_u, mu, jacp, drdy, drdp, qfun, q )

! end subroutine integrate_fatode_adj_erk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! adj rk module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine integrate_fatode_adj_rk( tin, tout, nvar, np, nadj, nnz, y, lambda, &
!      atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit,  icntrl_u, rcntrl_u, istatus_u, &
!      rstatus_u, ierr_u, mu, jacp, drdy, drdp, qfun, q ) bind(c)

!   use iso_c_binding
!   use rk_adj_f90_integrator
!   implicit none

!    integer, intent(in) :: nvar,np,nadj
!    integer, intent(in) :: nnz

!    real(dp), intent(inout) :: y(nvar)

!    real(dp), intent(inout) :: q(nadj)
!    real(dp), intent(inout)  :: lambda(nvar,nadj)
!    real(dp), intent(inout) ::  mu(np, nadj)
!    real(dp), intent(in) :: atol_adj(nvar,nadj), rtol_adj(nvar,nadj)
!    real(dp), intent(in)  ::  atol(nvar), rtol(nvar)
!    real(dp), intent(in) :: tin  ! start time
!    real(dp), intent(in) :: tout ! end time

!    integer,          intent(in)  :: icntrl_u(20)
!    real(dp), intent(in)  :: rcntrl_u(20)
!    integer,          intent(out) :: istatus_u(20)
!    real(dp), intent(out) :: rstatus_u(20)
!    integer,          intent(out) :: ierr_u

!    double precision :: rcntrl(20), rstatus(20), t1, t2
!    integer       :: icntrl(20), istatus(20), ierr

!   interface
!      subroutine fun(n, t, y, fy) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fy(n)
!      end subroutine fun
!   end interface

!   interface
!      subroutine jac(n, t, y, fjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fjac(n, n)
!      end subroutine jac
!   end interface

!   interface
!      subroutine adjinit(n,np,nadj,t,y,lambda,mu) bind (c)
!        use iso_c_binding
!        integer,intent(in) :: n,np,nadj
!        double precision,intent(in) :: t,y(n),lambda(n,nadj)
!        double precision, optional,intent(inout) :: mu(np,nadj)
!      end subroutine adjinit
!   end interface

!   interface
!      subroutine drdp(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine drdp
!   end interface
  
!   interface
!      subroutine drdy(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine drdy
!   end interface

!   interface
!      subroutine jacp(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine jacp
!   end interface

!   interface
!      subroutine qfun(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine qfun
!   end interface

!   call integrate_adj(nvar, np, nadj, nnz, y, lambda, tin, tout,  &
!        atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit, icntrl_u, &
!        rcntrl_u, istatus_u, rstatus_u, ierr_u, mu, jacp, drdy, drdp, &
!        qfun, q )

! end subroutine integrate_fatode_adj_rk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! adj ros module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine integrate_fatode_adj_ros( tin, tout, nvar, np, nadj, nnz, y, lambda, &
!      atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit, hesstr_vec,&
!      icntrl_u, rcntrl_u, istatus_u, &
!      rstatus_u, ierr_u, mu, jacp, drdy, drdp, hesstr_vec_f_py, hesstr_vec_r_py, hesstr_vec_r, qfun, q ) bind(c)

!   use iso_c_binding
!   use ros_adj_f90_integrator
!   implicit none

!    integer, intent(in) :: nvar,np,nadj
!    integer, intent(in) :: nnz

!    real(dp), intent(inout) :: y(nvar)

!    real(dp), intent(inout) :: q(nadj)
!    real(dp), intent(inout)  :: lambda(nvar,nadj)
!    real(dp), intent(inout) ::  mu(np, nadj)
!    real(dp), intent(in) :: atol_adj(nvar,nadj), rtol_adj(nvar,nadj)
!    real(dp), intent(in)  ::  atol(nvar), rtol(nvar)
!    real(dp), intent(in) :: tin  ! start time
!    real(dp), intent(in) :: tout ! end time

!    integer,          intent(in)  :: icntrl_u(20)
!    real(dp), intent(in)  :: rcntrl_u(20)
!    integer,          intent(out) :: istatus_u(20)
!    real(dp), intent(out) :: rstatus_u(20)
!    integer,          intent(out) :: ierr_u

!    double precision :: rcntrl(20), rstatus(20), t1, t2
!    integer       :: icntrl(20), istatus(20), ierr

!   interface
!      subroutine fun(n, t, y, fy) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fy(n)
!      end subroutine fun
!   end interface

!   interface
!      subroutine jac(n, t, y, fjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fjac(n, n)
!      end subroutine jac
!   end interface

!   interface
!      subroutine adjinit(n,np,nadj,t,y,lambda,mu) bind (c)
!        use iso_c_binding
!        integer,intent(in) :: n,np,nadj
!        double precision,intent(in) :: t,y(n),lambda(n,nadj)
!        double precision, optional,intent(inout) :: mu(np,nadj)
!      end subroutine adjinit
!   end interface

!   interface
!      subroutine drdp(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine drdp
!   end interface
  
!   interface
!      subroutine drdy(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine drdy
!   end interface

!   interface
!      subroutine jacp(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine jacp
!   end interface

!   interface
!      subroutine qfun(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine qfun
!   end interface

!   interface
!      subroutine hesstr_vec(n,t,y,u,v,hv) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in)    :: t,y(n),u(n),v(n)
!        double precision, intent(inout) :: hv(n)
!      end subroutine
!   end interface

!   interface
!      subroutine hesstr_vec_f_py(ny,np,t,y,u,k,htvg) bind (c)
!        use iso_c_binding
!        ! htvg = (f_py x k)^t * u = (d(f_p^t * u)/dy) * k
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvg(np)
!      end subroutine
!   end interface
    
!   interface
!      subroutine hesstr_vec_r_py(ny,np,t,y,u,k,htvr) bind (c)
!        use iso_c_binding
!        ! htvr = (f_py x k)^t * u = (d(f_p^t * u)/dy) * k
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvr(np)
!      end subroutine
!   end interface

!   interface
!      subroutine hesstr_vec_r(ny,np,t,y,u,k,htvr) bind (c)
!        use iso_c_binding
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvr(np)
!      end subroutine
!   end interface

!   call integrate_adj( nvar, np, nadj, nnz,  y, lambda, mu, tin, tout,  &
!        atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit, hesstr_vec,     &
!        jacp, drdy, drdp, hesstr_vec_f_py, hesstr_vec_r_py, hesstr_vec_r,  &
!        icntrl_u, rcntrl_u, istatus_u, rstatus_u, q, qfun ) 

! end subroutine integrate_fatode_adj_ros

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! adj sdirk module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine integrate_fatode_adj_sdirk( tin, tout, nvar, np, nadj, nnz, y, lambda, &
!      atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit,  icntrl_u, rcntrl_u, istatus_u, &
!      rstatus_u, ierr_u, mu, jacp, drdy, drdp, qfun, q ) bind(c)

!   use iso_c_binding
!   use sdirk_adj_f90_integrator
!   implicit none

!    integer, intent(in) :: nvar,np,nadj
!    integer, intent(in) :: nnz

!    real(dp), intent(inout) :: y(nvar)

!    real(dp), intent(inout) :: q(nadj)
!    real(dp), intent(inout)  :: lambda(nvar,nadj)
!    real(dp), intent(inout) ::  mu(np, nadj)
!    real(dp), intent(in) :: atol_adj(nvar,nadj), rtol_adj(nvar,nadj)
!    real(dp), intent(in)  ::  atol(nvar), rtol(nvar)
!    real(dp), intent(in) :: tin  ! start time
!    real(dp), intent(in) :: tout ! end time

!    integer,          intent(in)  :: icntrl_u(20)
!    real(dp), intent(in)  :: rcntrl_u(20)
!    integer,          intent(out) :: istatus_u(20)
!    real(dp), intent(out) :: rstatus_u(20)
!    integer,          intent(out) :: ierr_u

!    double precision :: rcntrl(20), rstatus(20), t1, t2
!    integer       :: icntrl(20), istatus(20), ierr

!   interface
!      subroutine fun(n, t, y, fy) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fy(n)
!      end subroutine fun
!   end interface

!   interface
!      subroutine jac(n, t, y, fjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n
!        double precision, intent(in) :: t, y(n)
!        double precision, intent(inout) :: fjac(n, n)
!      end subroutine jac
!   end interface

!   interface
!      subroutine adjinit(n,np,nadj,t,y,lambda,mu) bind (c)
!        use iso_c_binding
!        integer,intent(in) :: n,np,nadj
!        double precision,intent(in) :: t,y(n),lambda(n,nadj)
!        double precision, optional,intent(inout) :: mu(np,nadj)
!      end subroutine adjinit
!   end interface

!   interface
!      subroutine drdp(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine drdp
!   end interface
  
!   interface
!      subroutine drdy(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine drdy
!   end interface

!   interface
!      subroutine jacp(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine jacp
!   end interface

!   interface
!      subroutine qfun(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine qfun
!   end interface

!   call integrate_adj( nvar, np, nadj, nnz, y, lambda, mu, tin, tout, &
!        atol_adj, rtol_adj, atol, rtol, fun, jac, adjinit, jacp, drdy, &
!        drdp, icntrl_u, rcntrl_u, istatus_u, rstatus_u, ierr_u, q, qfun )

! end subroutine integrate_fatode_adj_sdirk
end module fatode_c
