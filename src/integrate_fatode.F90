! torsten's interface to FATODE library

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FWD ERK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_fwd_erk( TIN, TOUT, NVAR, VAR, RTOL, ATOL, FUN, &
     ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U ) bind(c)

  use iso_c_binding
  use ERK_f90_Integrator
  implicit none

  INTEGER, INTENT(IN) :: NVAR
  DOUBLE PRECISION, INTENT(INOUT) :: VAR(NVAR)
  DOUBLE PRECISION, INTENT(IN) :: RTOL(NVAR)
  DOUBLE PRECISION, INTENT(IN) :: ATOL(NVAR)
  DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
  DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

  ! control and output
  INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
  DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
  INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
  DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
  INTEGER,          INTENT(OUT) :: Ierr_U

  interface
     subroutine fun(n, t, y, fy) bind (c)
       use iso_c_binding
       integer, intent(in) :: n
       double precision, intent(in) :: t, y(n)
       double precision, intent(inout) :: fy(n)
     end subroutine fun
  end interface

  call integrate( TIN, TOUT, NVAR, VAR, RTOL, ATOL, FUN, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

end SUBROUTINE integrate_fatode_fwd_erk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FWD RK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_fwd_rk( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
     ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U ) bind(c)

  use iso_c_binding
  use RK_f90_Integrator
  implicit none

  INTEGER, INTENT(IN) :: N, NNZERO
  DOUBLE PRECISION, INTENT(INOUT) :: VAR(N)
  DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
  DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
  DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
  DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

  ! control and output
  INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
  DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
  INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
  DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
  INTEGER,          INTENT(OUT) :: Ierr_U

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

  call integrate( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

end SUBROUTINE integrate_fatode_fwd_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FWD ROS module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_fwd_ros( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
     ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U ) bind(c)

  use iso_c_binding
  use ROS_f90_Integrator
  implicit none

  INTEGER, INTENT(IN) :: N, NNZERO
  DOUBLE PRECISION, INTENT(INOUT) :: VAR(N)
  DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
  DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
  DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
  DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

  ! control and output
  INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
  DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
  INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
  DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
  INTEGER,          INTENT(OUT) :: Ierr_U

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

  call integrate( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

end SUBROUTINE integrate_fatode_fwd_ros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FWD SDIRK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_fwd_sdirk( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
     ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U ) bind(c)

  use iso_c_binding
  use SDIRK_f90_Integrator
  implicit none

  INTEGER, INTENT(IN) :: N, NNZERO
  DOUBLE PRECISION, INTENT(INOUT) :: VAR(N)
  DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
  DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
  DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
  DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

  ! control and output
  INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
  DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
  INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
  DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
  INTEGER,          INTENT(OUT) :: Ierr_U

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

  call integrate( TIN, TOUT, N, NNZERO, VAR, RTOL, ATOL, FUN, JAC, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

end SUBROUTINE integrate_fatode_fwd_sdirk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TLM ERK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_tlm_erk( TIN, TOUT, NVAR, NTLM, VAR, VAR_TLM, ATOL_TLM, &
     RTOL_TLM, ATOL, RTOL,&
     FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U ) bind(c)

  use iso_c_binding
  use ERK_TLM_f90_Integrator
  implicit none

   INTEGER, INTENT(IN) :: NVAR,NTLM
   DOUBLE PRECISION, INTENT(INOUT) :: VAR(NVAR),VAR_TLM(NVAR,NTLM),&
                           RTOL_TLM(NVAR,NTLM), ATOL_TLM(NVAR,NTLM)
   DOUBLE PRECISION, INTENT(IN) :: RTOL(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: ATOL(NVAR)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

   INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
   INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
   INTEGER,          INTENT(OUT) :: IERR_U

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

  call INTEGRATE_TLM(TIN, TOUT, NVAR, NTLM, VAR, VAR_TLM, ATOL_TLM, &
       RTOL_TLM, ATOL, RTOL,&
       FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

end SUBROUTINE integrate_fatode_tlm_erk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TLM RK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE integrate_fatode_tlm_rk( TIN, TOUT, N, NTLM, NNZERO, Y, Y_TLM, ATOL_TLM, RTOL_TLM, &
!  ATOL, RTOL, FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U ) bind(c)

!   use iso_c_binding
!   use RK_TLM_f90_Integrator
!   implicit none

!   INTEGER :: N, NTLM, NNZERO
!    DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Y_TLM(N,NTLM),&
!                            RTOL_TLM(N,NTLM), ATOL_TLM(N,NTLM)
!    DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
!    DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
!    DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
!    DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

!    INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
!    DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
!    INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
!    DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
!    INTEGER,          INTENT(OUT) :: IERR_U

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

!   call INTEGRATE_TLM( N, NTLM, NNZERO, Y, Y_TLM, TIN, TOUT, ATOL_TLM, RTOL_TLM, &
!        ATOL, RTOL, FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

! end SUBROUTINE integrate_fatode_tlm_rk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TLM ROS module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_tlm_ros( TIN, TOUT, N, NTLM, NNZERO, Y, Y_TLM, ATOL_TLM, RTOL_TLM, &
     ATOL, RTOL, FUN, JAC, HESS_VEC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U ) bind(c)

  use iso_c_binding
  use ROS_TLM_f90_Integrator
  implicit none

  INTEGER :: N, NTLM, NNZERO
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Y_TLM(N,NTLM),&
                           RTOL_TLM(N,NTLM), ATOL_TLM(N,NTLM)
   DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
   DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

   INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
   INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
   INTEGER,          INTENT(OUT) :: IERR_U

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

  call INTEGRATE_TLM( N, NTLM, NNZERO, Y, Y_TLM, TIN, TOUT, ATOL_TLM, RTOL_TLM, &
       ATOL, RTOL, FUN, JAC, HESS_VEC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U )

end SUBROUTINE integrate_fatode_tlm_ros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TLM SDIRK module wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE integrate_fatode_tlm_sdirk( TIN, TOUT, N, NTLM, NNZERO, Y, Y_TLM, ATOL_TLM, RTOL_TLM, &
 ATOL, RTOL, FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U ) bind(c)

  use iso_c_binding
  use SDIRK_TLM_f90_Integrator
  implicit none

  INTEGER :: N, NTLM, NNZERO
   DOUBLE PRECISION, INTENT(INOUT) :: Y(N), Y_TLM(N,NTLM),&
                           RTOL_TLM(N,NTLM), ATOL_TLM(N,NTLM)
   DOUBLE PRECISION, INTENT(IN) :: RTOL(N)
   DOUBLE PRECISION, INTENT(IN) :: ATOL(N)
   DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
   DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

   INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
   DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
   INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
   DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
   INTEGER,          INTENT(OUT) :: IERR_U

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

  call INTEGRATE_TLM( N, NTLM, NNZERO, Y, Y_TLM, TIN, TOUT, ATOL_TLM, RTOL_TLM, &
       ATOL, RTOL, FUN, JAC, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U )

end SUBROUTINE integrate_fatode_tlm_sdirk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ADJ ERK module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE integrate_fatode_adj_erk( TIN, TOUT, NVAR, NP, NADJ, NNZ, Y, Lambda, &
!      ATOL, RTOL, FUN, JAC, AdjInit,  ICNTRL_U, RCNTRL_U, ISTATUS_U, &
!      RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, QFUN, Q ) bind(c)

!   use iso_c_binding
!   use ERK_ADJ_f90_Integrator
!   implicit none

!    INTEGER, INTENT(IN) :: NVAR,NP,NADJ
!    INTEGER, INTENT(IN) :: NNZ

!    DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)

!    DOUBLE PRECISION, INTENT(INOUT) :: Q(NADJ)
!    DOUBLE PRECISION, INTENT(INOUT)  :: Lambda(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(INOUT) ::  Mu(NP, NADJ)
!    DOUBLE PRECISION, INTENT(IN)  ::  ATOL(NVAR), RTOL(NVAR)
!    DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
!    DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

!    INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
!    DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
!    INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
!    DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
!    INTEGER,          INTENT(OUT) :: Ierr_U

!    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
!    INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

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
!      subroutine DRDP(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine DRDP
!   end interface
  
!   interface
!      subroutine DRDY(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine DRDY
!   end interface

!   interface
!      subroutine JACP(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine JACP
!   end interface

!   interface
!      subroutine QFUN(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine QFUN
!   end interface

!   call INTEGRATE_ADJ( NVAR, NP, NADJ, NNZ, Y, Lambda, TIN, TOUT, &
!        ATOL, RTOL, FUN, JAC, AdjInit,  ICNTRL_U, RCNTRL_U, ISTATUS_U, &
!        RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, QFUN, Q )

! end SUBROUTINE integrate_fatode_adj_erk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ADJ RK module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE integrate_fatode_adj_rk( TIN, TOUT, NVAR, NP, NADJ, NNZ, Y, Lambda, &
!      ATOL_ADJ, RTOL_ADJ, ATOL, RTOL, FUN, JAC, AdjInit,  ICNTRL_U, RCNTRL_U, ISTATUS_U, &
!      RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, QFUN, Q ) bind(c)

!   use iso_c_binding
!   use RK_ADJ_f90_Integrator
!   implicit none

!    INTEGER, INTENT(IN) :: NVAR,NP,NADJ
!    INTEGER, INTENT(IN) :: NNZ

!    DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)

!    DOUBLE PRECISION, INTENT(INOUT) :: Q(NADJ)
!    DOUBLE PRECISION, INTENT(INOUT)  :: Lambda(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(INOUT) ::  Mu(NP, NADJ)
!    DOUBLE PRECISION, INTENT(IN) :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(IN)  ::  ATOL(NVAR), RTOL(NVAR)
!    DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
!    DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

!    INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
!    DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
!    INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
!    DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
!    INTEGER,          INTENT(OUT) :: Ierr_U

!    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
!    INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

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
!      subroutine DRDP(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine DRDP
!   end interface
  
!   interface
!      subroutine DRDY(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine DRDY
!   end interface

!   interface
!      subroutine JACP(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine JACP
!   end interface

!   interface
!      subroutine QFUN(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine QFUN
!   end interface

!   call INTEGRATE_ADJ(NVAR, NP, NADJ, NNZ, Y, Lambda, TIN, TOUT,  &
!        ATOL_adj, RTOL_adj, ATOL, RTOL, FUN, JAC, AdjInit, ICNTRL_U, &
!        RCNTRL_U, ISTATUS_U, RSTATUS_U, IERR_U, Mu, JACP, DRDY, DRDP, &
!        QFUN, Q )

! end SUBROUTINE integrate_fatode_adj_rk

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ADJ ROS module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE integrate_fatode_adj_ros( TIN, TOUT, NVAR, NP, NADJ, NNZ, Y, Lambda, &
!      ATOL_ADJ, RTOL_ADJ, ATOL, RTOL, FUN, JAC, AdjInit, HESSTR_VEC,&
!      ICNTRL_U, RCNTRL_U, ISTATUS_U, &
!      RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, HESSTR_VEC_F_PY, HESSTR_VEC_R_PY, HESSTR_VEC_R, QFUN, Q ) bind(c)

!   use iso_c_binding
!   use ROS_ADJ_f90_Integrator
!   implicit none

!    INTEGER, INTENT(IN) :: NVAR,NP,NADJ
!    INTEGER, INTENT(IN) :: NNZ

!    DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)

!    DOUBLE PRECISION, INTENT(INOUT) :: Q(NADJ)
!    DOUBLE PRECISION, INTENT(INOUT)  :: Lambda(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(INOUT) ::  Mu(NP, NADJ)
!    DOUBLE PRECISION, INTENT(IN) :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(IN)  ::  ATOL(NVAR), RTOL(NVAR)
!    DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
!    DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

!    INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
!    DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
!    INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
!    DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
!    INTEGER,          INTENT(OUT) :: Ierr_U

!    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
!    INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

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
!      subroutine DRDP(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine DRDP
!   end interface
  
!   interface
!      subroutine DRDY(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine DRDY
!   end interface

!   interface
!      subroutine JACP(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine JACP
!   end interface

!   interface
!      subroutine QFUN(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine QFUN
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
!      subroutine HESSTR_VEC_F_PY(ny,np,t,y,u,k,htvg) bind (c)
!        use iso_c_binding
!        ! htvg = (f_py x k)^T * u = (d(f_p^T * u)/dy) * k
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvg(np)
!      end subroutine
!   end interface
    
!   interface
!      subroutine HESSTR_VEC_R_PY(ny,np,t,y,u,k,htvr) bind (c)
!        use iso_c_binding
!        ! htvr = (f_py x k)^T * u = (d(f_p^T * u)/dy) * k
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvr(np)
!      end SUBROUTINE
!   end interface

!   interface
!      subroutine HESSTR_VEC_R(ny,np,t,y,u,k,htvr) bind (c)
!        use iso_c_binding
!        integer :: ny,np
!        double precision :: t,y(ny),u(ny),k(ny),htvr(np)
!      end SUBROUTINE
!   end interface

!   call INTEGRATE_ADJ( NVAR, NP, NADJ, NNZ,  Y, Lambda, Mu, TIN, TOUT,  &
!        ATOL_adj, RTOL_adj, ATOL, RTOL, FUN, JAC, ADJINIT, HESSTR_VEC,     &
!        JACP, DRDY, DRDP, HESSTR_VEC_F_PY, HESSTR_VEC_R_PY, HESSTR_VEC_R,  &
!        ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Q, QFUN ) 

! end SUBROUTINE integrate_fatode_adj_ros

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ADJ SDIRK module wrapper
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE integrate_fatode_adj_sdirk( TIN, TOUT, NVAR, NP, NADJ, NNZ, Y, Lambda, &
!      ATOL_ADJ, RTOL_ADJ, ATOL, RTOL, FUN, JAC, AdjInit,  ICNTRL_U, RCNTRL_U, ISTATUS_U, &
!      RSTATUS_U, Ierr_U, Mu, JACP, DRDY, DRDP, QFUN, Q ) bind(c)

!   use iso_c_binding
!   use SDIRK_ADJ_f90_Integrator
!   implicit none

!    INTEGER, INTENT(IN) :: NVAR,NP,NADJ
!    INTEGER, INTENT(IN) :: NNZ

!    DOUBLE PRECISION, INTENT(INOUT) :: Y(NVAR)

!    DOUBLE PRECISION, INTENT(INOUT) :: Q(NADJ)
!    DOUBLE PRECISION, INTENT(INOUT)  :: Lambda(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(INOUT) ::  Mu(NP, NADJ)
!    DOUBLE PRECISION, INTENT(IN) :: ATOL_adj(NVAR,NADJ), RTOL_adj(NVAR,NADJ)
!    DOUBLE PRECISION, INTENT(IN)  ::  ATOL(NVAR), RTOL(NVAR)
!    DOUBLE PRECISION, INTENT(IN) :: TIN  ! Start Time
!    DOUBLE PRECISION, INTENT(IN) :: TOUT ! End Time

!    INTEGER,          INTENT(IN)  :: ICNTRL_U(20)
!    DOUBLE PRECISION, INTENT(IN)  :: RCNTRL_U(20)
!    INTEGER,          INTENT(OUT) :: ISTATUS_U(20)
!    DOUBLE PRECISION, INTENT(OUT) :: RSTATUS_U(20)
!    INTEGER,          INTENT(OUT) :: Ierr_U

!    DOUBLE PRECISION :: RCNTRL(20), RSTATUS(20), T1, T2
!    INTEGER       :: ICNTRL(20), ISTATUS(20), Ierr

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
!      subroutine DRDP(nadj,n,nr,t,y,rp) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: rp(nr,nadj)
!      end subroutine DRDP
!   end interface
  
!   interface
!      subroutine DRDY(nadj,n,nr,t,y,ry) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: nadj,n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: ry(nr,nadj)
!      end subroutine DRDY
!   end interface

!   interface
!      subroutine JACP(n,np,t,y,fpjac) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,np
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: fpjac(n,np)
!      end subroutine JACP
!   end interface

!   interface
!      subroutine QFUN(n,nr,t,y,r) bind (c)
!        use iso_c_binding
!        integer, intent(in) :: n,nr
!        double precision, intent(in)    :: t,y(n)
!        double precision, intent(inout) :: r(nr)
!      end subroutine QFUN
!   end interface

!   call INTEGRATE_ADJ( NVAR, NP, NADJ, NNZ, Y, Lambda, Mu, TIN, TOUT, &
!        ATOL_adj, RTOL_adj, ATOL, RTOL, FUN, JAC, ADJINIT, JACP, DRDY, &
!        DRDP, ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U, Q, QFUN )

! end SUBROUTINE integrate_fatode_adj_sdirk
