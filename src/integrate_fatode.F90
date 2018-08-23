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
  ! Optional input parameters and statistics
  INTEGER,       INTENT(IN),  OPTIONAL :: ICNTRL_U(20)
  DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: RCNTRL_U(20)
  INTEGER,       INTENT(OUT), OPTIONAL :: ISTATUS_U(20)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: RSTATUS_U(20)
  INTEGER,       INTENT(OUT), OPTIONAL :: Ierr_U

  interface
     subroutine fun(nvar, t, yy, fy) bind (c)
       use iso_c_binding
       integer nvar
       double precision :: fy(nvar), yy(nvar), t
     end subroutine fun
  end interface

  call integrate( TIN, TOUT, NVAR, VAR, RTOL, ATOL, FUN, &
       ICNTRL_U, RCNTRL_U, ISTATUS_U, RSTATUS_U, Ierr_U )

end SUBROUTINE integrate_fatode_fwd_erk
