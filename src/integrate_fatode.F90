! torsten's interface to FATODE library

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
