$ Declare
   SUBROUTINE ReadStressRanges (nr, ntot, sfile, ncc)
!
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)           :: nr        ! Period number
   DOUBLE PRECISION, INTENT(IN)  :: ntot      ! Number of cycles in that period
   CHARACTER, INTENT(IN)         :: sfile*100 ! File with stress range histogram for that period
   INTEGER, INTENT(OUT)          :: ncc       ! Number of rows in that file
!
! Declarations
   DOUBLE PRECISION              :: factn, nfile
   INTEGER                       :: i, j
!
! Initialise stress spectrum file
   DO i = 1, 3
      DO j = 1, 100
		 Spect (i,j,nr) = 0.D0
	  END DO
   END DO
!
! Stress ranges
   nfile = 0.D0
   OPEN(8,file=sfile,status='old')
   DO j = 1, 101
      READ (8,*,END=99) Spect(1,j,nr), Spect(2,j,nr), Spect(3,j,nr)
	  nfile = nfile + Spect(1,j,nr)
	  IF (j.EQ.101) STOP 'Error in subroutine ReadStressRanges: Max. no of rows in stress spectrum file = 100'
   END DO
99 CONTINUE
   ncc = j-1
   factn = ntot / nfile
   write(*,*) factn, ntot, nfile
   DO j = 1, ncc
      Spect(1,j,nr) = Spect(1,j,nr) * factn
   END DO
   CLOSE(8)
!
! Einde subroutine
   RETURN
   END