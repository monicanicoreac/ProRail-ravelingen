$ Declare
   SUBROUTINE Ncrack (ntot, ULS, SMult, Mat, Geo, ad, Csr, Z) 
!                    
! Subroutine voor het bepalen van het aantal kritische scheuren, en het wegschrijven van de resultaten
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)  :: ntot(ninsp+1)  ! total number of cycles between inspections
   DOUBLE PRECISION, INTENT(IN)  :: ULS(4)		   ! load variables in ULS
   DOUBLE PRECISION, INTENT(IN)  :: SMult(nloc)    ! stress multiplication factor
   DOUBLE PRECISION, INTENT(IN)  :: Mat(14,nloc)   ! material variables
   DOUBLE PRECISION, INTENT(IN)  :: Geo(8,nloc)    ! geometry variables
   DOUBLE PRECISION, INTENT(IN)  :: ad(ninsp,nloc) ! detectable crack size
   DOUBLE PRECISION, INTENT(IN)  :: Csr(nmaxCsr,nloc) ! stress multiplication factor per stress range
   DOUBLE PRECISION, INTENT(OUT) :: Z              ! limit state
!
! Declarations
   DOUBLE PRECISION, ALLOCATABLE :: nm(:)          ! Cycles to failure
   LOGICAL, POINTER              :: det(:)         ! TRUE if the crack is detected
   INTEGER, POINTER              :: ndet(:)		   ! Number of detected cracks
   DOUBLE PRECISION              :: nkrit, nti, nres
   INTEGER	                     :: i
   INTEGER                       :: loc
   LOGICAL                       :: hit
   LOGICAL                       :: fail           ! TRUE if failure calculated before last inspection
!
! Initialisatie
   ALLOCATE (nm(nloc))
   ALLOCATE (det(ninsp))
   ALLOCATE (ndet(ninsp))
   nti = 0.D0
   DO i = 1, ninsp
      ndet(i) = 0
	  nti = nti + ntot(i)
   END DO
   nkrit = nti + ntot(ninsp+1)
   hit = .true.
!
! Calculation routine
   DO loc = 1, nloc
      CALL FracMap (ntot, loc, ULS, SMult(loc), Mat, Geo, ad, Csr, nm(loc), det, fail)
      IF (fail) hit = .false.        ! failure has occurred before last inspection
      DO i = 1, ninsp
	     IF (det(i)) THEN ! Crack was detected and is repaired. It is assumed that it has sufficient life.
		    nm(loc) = nkrit + 1.D0
		    ndet(i) = ndet(i) + 1
         END IF
	  END DO
   END DO
!
! Write in sequence of failure
   CALL Sort (nloc, nm)
   nres = nm(1)-nti
   DO i = 1, ninsp
      IF (ndet(i).NE.nd(i)) hit = .false. ! MC run does not match inspection result
   END DO
   IF (hit) THEN     ! MC run matches inspection result, failure has not occured before last inspection
      IF (nres.LT.0.D0) THEN
	     WRITE (*,*) fail, nres, nm(1), nti, ndet, nkrit 
		 STOP 'fout'
	  END IF
	  Z = nres/ntot(ninsp+1) - 1.D0
	  WRITE(10,'(2ES13.3, 2ES18.3, I4, 10I5)') nm(1), nres, nm(1)/nkrit, nres/ntot(ninsp+1), (ndet(i),i=1,ninsp)
!!	  CALL WriteData (12, 14, nkrit, ULS, SMult, Mat, Geo, ad, Csr, nm, ndet)
   ELSE ! MC run does not match inspection result
	  Z = 1.D0
	  WRITE(11,'(2ES13.3, 2ES18.3, I4, 10I5)') nm(1), nres, nm(1)/nkrit, nres/ntot(ninsp+1), (ndet(i),i=1,ninsp)
!!	  CALL WriteData (13, 15, nkrit, ULS, SMult, Mat, Geo, ad, Csr, nm, ndet)
   END IF
!
! Einde routine
   RETURN
   END