$ Declare
   SUBROUTINE FracMAp (ntot, loc, ULS, SMult, Mat, Geo, ad, Csr, nm, det, fail)
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(IN)  :: ntot(ninsp+1)  ! total number of cycles between inspections
   INTEGER, INTENT(IN)           :: loc            ! location number
   DOUBLE PRECISION, INTENT(IN)  :: ULS(4)         ! CrackGrowthParameters
   DOUBLE PRECISION, INTENT(IN)  :: SMult          ! stress multiplication factor
   DOUBLE PRECISION, INTENT(IN)  :: Mat(14,nloc)   ! CrackGrowthParameters
   DOUBLE PRECISION, INTENT(IN)  :: Geo(8,nloc)    ! geometry variables
   DOUBLE PRECISION, INTENT(IN)  :: ad(ninsp,nloc) ! detectable crack size
   DOUBLE PRECISION, INTENT(IN)  :: Csr(nmaxCsr,nloc) ! stress multiplication factor per stress range
   DOUBLE PRECISION, INTENT(OUT) :: nm
   LOGICAL, INTENT(OUT)			 :: det(ninsp)     ! TRUE if crack is detected
   LOGICAL, INTENT(OUT)			 :: fail           ! TRUE if failure occurs before all cycles are passed
!
! Declarations
   DOUBLE PRECISION, ALLOCATABLE :: GeoInOut(:)
   INTEGER                       :: i, ifin
   LOGICAL                       :: tdet
   LOGICAL                       :: dummy1, dummy2
!
! Initialisation
   ALLOCATE (GeoInOut(3))
   GeoInOut(1) = Geo(1,loc)
   GeoInOut(2) = Geo(2,loc)
   GeoInOut(3) = 0.D0   ! 1.D0 = through thickness flaw, 0.D0 = surface flaw (only relevant for 3D geometries)
   nm = 0.D0
   DO i = 2, ninsp
      det(i) = .false.
   END DO
!
! Calculate life if inspections are not carried out
   IF (ninsp.EQ.0) THEN
      CALL life (1, ntot(1), loc, ULS, SMult, Mat, Geo, GeoInOut, 1.D6, Csr, nm, dummy1, dummy2)
	  fail = .false.
   ELSE
!
! Calculate life if inspections are carried out
      CALL life (1, ntot(1), loc, ULS, SMult, Mat, Geo, GeoInOut, ad(1,loc), Csr, nm, det(1), fail)
	  tdet = det(1)
      IF (ninsp.GT.1)  THEN
         DO i = 2, ninsp
	        IF ((.not.tdet).AND.(.not.fail)) CALL life (i, ntot(i), loc, ULS, SMult, Mat, Geo, GeoInOut, ad(i,loc), Csr, nm, det(i), fail)
			IF (det(i)) tdet = .true.
         END DO
      END IF
      ifin = ninsp+1
      IF ((.not.tdet).AND.(.not.fail)) CALL life (ifin, ntot(ifin), loc, ULS, SMult, Mat, Geo, GeoInOut, 1.D6, Csr, nm, dummy1, dummy2)
   END IF
!
! End subroutine
   RETURN
   END