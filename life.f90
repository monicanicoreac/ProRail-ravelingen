$ Declare
   SUBROUTINE life (iper, ntot, loc, ULS, SMult, Mat, Geo, GeoInOut, ad, Csr, nm, det, fail)
!
! Calculation of the life
!
! Parameters
   USE PublicPar
   IMPLICIT NONE
!
   INTEGER, INTENT(IN)             :: iper         ! period, i.e. 3rd dimension of stress range histogram
   DOUBLE PRECISION, INTENT(IN)    :: ntot         ! number of cycles between inspections or after last insp.
   INTEGER, INTENT(IN)             :: loc          ! location number
   DOUBLE PRECISION, INTENT(IN)    :: ULS(4)       ! Load variables in ultimate limit state
   DOUBLE PRECISION, INTENT(IN)    :: SMult        ! stress multiplication factor
   DOUBLE PRECISION, INTENT(IN)    :: Mat(14,nloc) ! material variables
   DOUBLE PRECISION, INTENT(IN)    :: Geo(8,nloc)  ! geometry variables
   DOUBLE PRECISION, INTENT(INOUT) :: GeoInOut(3)  ! geometry variables
   DOUBLE PRECISION, INTENT(IN)    :: ad           ! detectable crack size
   DOUBLE PRECISION, INTENT(IN)    :: Csr(nmaxCsr,nloc) ! stress multiplication factor per stress range
   DOUBLE PRECISION, INTENT(INOUT) :: nm           ! number of cycles (to failure)
   LOGICAL, INTENT(OUT)			   :: det          ! TRUE = crack is detected
   LOGICAL, INTENT(OUT)            :: fail         ! TRUE = failure has occurred
!
! Declarations
   DOUBLE PRECISION    :: a, c
   DOUBLE PRECISION    :: a0, c0
   DOUBLE PRECISION    :: amax
   DOUBLE PRECISION    :: W, B, tw, L, phi, rho
   DOUBLE PRECISION    :: fy, fu
   DOUBLE PRECISION    :: K1C, FAD
   DOUBLE PRECISION    :: m1, m2
   DOUBLE PRECISION    :: DK0
   DOUBLE PRECISION    :: A1, A2
   DOUBLE PRECISION    :: DKtr, mstar, Astar
   DOUBLE PRECISION    :: DKa, DKc
   DOUBLE PRECISION    :: Yma, Yba, Ymc, Ybc
   DOUBLE PRECISION    :: Spm, Spb, Ssm, Ssb, Temp
   DOUBLE PRECISION    :: Sref, Lr, Kr, Kra, Krc, fKrLr
   INTEGER 			   :: tt
   DOUBLE PRECISION    :: n, dn
   DOUBLE PRECISION    :: dadn, dadn0, sumdadn
   DOUBLE PRECISION    :: dcdn, dcdn0, sumdcdn
   DOUBLE PRECISION    :: nc
   INTEGER             :: j
   INTEGER             :: iC
   DOUBLE PRECISION    :: ma
   DOUBLE PRECISION    :: Sfact
   DOUBLE PRECISION    :: DSm, DSb
   INTEGER             :: err
!
! Constants
   DOUBLE PRECISION    :: pi, thetac
   pi = 3.14159265358979D0
   thetac =  pi/2.D0
!
! ULS load
   Spm   = ULS(1)
   Spb   = ULS(2)
   Ssm   = ULS(3)
   Ssb   = ULS(4)
!
! Geometry
   IF ((FlawType.EQ.1).OR.(FlawType.EQ.3).OR.(FlawType.EQ.4).OR.(FlawType.EQ.6).OR.(FlawType.EQ.9)) THEN
      tt = 1
   ELSE
      tt = NINT(GeoInOut(3))
   END IF
   a0   = GeoInOut(1)
   c0   = GeoInOut(2)
   a    = a0
   c    = c0
   W    = Geo(3,loc)
   B    = Geo(4,loc)
   tw   = Geo(5,loc)
   L    = Geo(6,loc)
   phi  = Geo(7,loc)
   rho  = Geo(8,loc)
   amax = B * 0.96D0
!
! Material variables
   mstar = 22
   fy  = Mat(1,loc)
   fu  = Mat(2,loc)
   K1C = Mat(3,loc)
   FAD = Mat(4,loc)
   m1  = Mat(5,loc)
   m2  = Mat(6,loc)
   DK0 = Mat(7,loc)
!
! Initialisation
   ma      = 5.D-2
   n       = 0.D0
   dn      = 0.D0
   dadn0   = 0.D0
   iC      = 0
   det     = .false.
   fail    = .true.
   IF (tt.EQ.0) THEN
!
! Calculate crack grow
	  DO WHILE (a.LT.amax)
	     IF (a*(1.D0+ma).GT.amax) ma = (amax - a) / a
! 
! Calculate geometry factors Ya
         IF (LoadType.EQ.1) THEN
		    CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, thetac, Yma, err )
		    CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, 0.0D0,  Ymc, err )
		    Yba = 0.D0
			Ybc = 0.D0
         ELSEIF (LoadType.EQ.2) THEN
		    CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, thetac, Yba, err )
		    CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, 0.0D0,  Ybc, err )
		    Yma = 0.D0
			Ymc = 0.D0
		 ELSE
		    CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, thetac, Yma, err )
		    CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, thetac, Yba, err )
		    CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, 0.0D0,  Ymc, err )
		    CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, 0.0D0,  Ybc, err )
		 END IF
!!		 CALL StressFact ( DetailType, c, Sfact)
         Sfact = 1.D0
!
! Determine crack growth parameters
	     IF (a.LT.((amax - a0)/5.D0*1.D0+a0)) THEN
		    A1 = Mat(8,loc) * Mat(10,loc)
		    A2 = Mat(9,loc) * Mat(10,loc)
	     ELSE IF (a.LT.((amax - a0)/5.D0*2.D0+a0)) THEN
		    A1 = Mat(8,loc) * Mat(11,loc)
		    A2 = Mat(9,loc) * Mat(11,loc)
	     ELSE IF (a.LT.((amax - a0)/5.D0*3.D0+a0)) THEN
		    A1 = Mat(8,loc) * Mat(12,loc)
		    A2 = Mat(9,loc) * Mat(12,loc)
	     ELSE IF (a.LT.((amax - a0)/5.D0*4.D0+a0)) THEN
		    A1 = Mat(8,loc) * Mat(13,loc)
		    A2 = Mat(9,loc) * Mat(13,loc)
	     ELSE
		    A1 = Mat(8,loc) * Mat(14,loc)
		    A2 = Mat(9,loc) * Mat(14,loc)
	     END IF
	     DKtr = (A2/A1)**(1.D0/(M1-M2))
!         IF (DKtr.GE.DK0) THEN
!            Astar = A1a / (DK0**(mstar-m1))
!         ELSE
!            Astar = A2a / (DK0**(mstar-m2))
!         END IF
!
! Calculate da/dn and dc/dn
		 sumdadn = 0
		 sumdcdn = 0
		 DO j = 1, nccmax
		    iC = iC + 1
		    nc = Spect(1,j,iper)
			DSm = Spect(2,j,iper)
            DSb = Spect(3,j,iper)
			DKa = Smult * Csr(iC, loc) * DSQRT(pi * a) * (DSm * Yma * (Bnom/B)**p + DSb * Yba * (Bnom/B)**2.D0)
			DKc = Sfact * SMult * Csr(iC, loc) * DSQRT(pi * a) * (DSm * Ymc * (Bnom/B)**p + DSb * Ybc * (Bnom/B)**2.D0)
			CALL CrackGrowthLaw (A1, m1, A2, m2, DK0, K1C, fy, fu, (Spm+Spb), DKa, a, dadn)
			CALL CrackGrowthLaw (A1, m1, A2, m2, DK0, K1C, fy, fu, (Spm+Spb), DKc, a, dcdn)
            sumdadn = sumdadn + dadn * nc
            sumdcdn = sumdcdn + dcdn * nc
			IF (iC.EQ.nmaxCsr) iC = 0
		 END DO
!
! Update n, a and c
		 dadn = sumdadn / ntot
		 dcdn = sumdcdn / ntot
!
! Procedure if all stress intensity ranges (a, c, all cycles) are lower than the threshold
		 IF ((dadn .EQ. 0.D0) .AND. (dcdn .EQ. 0.D0)) THEN
            n = ntot
            fail = .false.
			GOTO 99
!
! Procedure with a as controlling variable
		 ELSE IF (dadn.GT.0.02D0*dcdn) THEN
			dn = a * ma / (dadn + dadn0) * 2.D0 
			IF (a.NE.a0) n = n + dn
!
! Check if the calculated number of cycles is larger than the applied number of cycles
			IF (n.GE.ntot) THEN
			   a = a + (ntot - n + dn) * (dadn + dadn0) / 2.D0
			   c = c + (ntot - n + dn) * (dcdn + dcdn0) / 2.D0
			   n = ntot
			   fail = .false.
		       GOTO 99
			END IF
			dadn0 = dadn
			dcdn0 = dcdn
			a = a * (1 + ma)
			c = c + dn * (dcdn + dcdn0) / 2.D0 
!
! Evaluation of fracture
		    CALL SIGrefC ( a, c, B, W, Spm, Spb, Sref)
		    Lr = Sref / fy
		    CALL FactorKr ( a, c, B, W, thetac, Yma, Yba, Spm, Spb, Ssm, Ssb, K1C, Lr, Kra )
		    CALL FactorKr ( a, c, B, W, 0.0D0, Yma, Yba, Spm, Spb, Ssm, Ssb, K1C, Lr, Krc )
		    Kr = MAX (Kra, Krc)
		    fKrLr = FAD - DSQRT(Kr*Kr+Lr*Lr)
!			WRITE(16,'(2F10.3, 2E15.3, 3F10.3)') a, c, n, ntot, Kr, Lr, fKrLr
		    IF (fKrLr .LT. 0.D0) GOTO 99
!
! Check if the calculated crack size is larger than the detectable crack size
			IF ((ttd.EQ.0).AND.(a.GE.ad)) THEN
			   det = .true.
		       fail = .false.
			   GOTO 99
			END IF
		 ELSE
!
! Procedure with c as controlling variable
			dn = c * ma / (dcdn + dcdn0) * 2.D0
			IF (a.NE.a0) n = n + dn
!
! Check if the calculated number of cycles is larger than the applied number of cycles
			IF (n.GE.ntot) THEN
			   a = a + (ntot - n + dn) * (dadn + dadn0) / 2.D0
			   c = c + (ntot - n + dn) * (dcdn + dcdn0) / 2.D0
			   n = ntot
			   fail = .false.
		       GOTO 99
			END IF
		    dadn0 = dadn
			dcdn0 = dcdn
			a = a + dn * (dadn + dadn0) / 2.D0
			c = c * (1 + ma)
!
! Evaluation of fracture
		    CALL SIGrefC ( a, c, B, W, Spm, Spb, Sref)
		    Lr = Sref / fy
		    CALL FactorKr ( a, c, B, W, thetac, Yma, Yba, Spm, Spb, Ssm, Ssb, K1C, Lr, Kra)
		    CALL FactorKr ( a, c, B, W, 0.0D0, Yma, Yba, Spm, Spb, Ssm, Ssb, K1C, Lr, Krc)
		    Kr = MAX (Kra, Krc)
		    fKrLr = FAD - DSQRT(Kr*Kr+Lr*Lr)
!			WRITE(16,'(2F10.3, 2E15.3, 3F10.3)') a, c, n, ntot, Kr, Lr, fKrLr
		    IF (fKrLr .LT. 0.D0) GOTO 99
!
! Check if the calculated crack size is larger than the detectable crack size
			IF ((ttd.EQ.0).AND.(a.GE.ad)) THEN
				det = .true.
		        fail = .false.
				GOTO 99
			END IF
		 END IF
	  END DO
!
! Prepare for through thickness crack calculation
	  a = c
	  a0 = c
	  dadn0 = dcdn
   END IF
   dn = 0.D0
   tt = 1
   c = a
!
! Calculate crack grow through thickness
   DO WHILE (a.LT.cmax)
	  IF (a*(1.D0+ma).GT.cmax) ma = (cmax - a) / a
! 
! Calculate geometry factors Ya
      IF (LoadType.EQ.1) THEN
	     CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, thetac, Yma, err )
	     Yba = 0.D0
      ELSEIF (LoadType.EQ.2) THEN
	     CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, thetac, Yba, err )
	     Yma = 0.D0
	  ELSE
	     CALL SICFactorY ( 1, a, c, B, W, L, tw, phi, rho, thetac, Yma, err )
	     CALL SICFactorY ( 2, a, c, B, W, L, tw, phi, rho, thetac, Yba, err )
	  END IF
!
! Calculate da/dn and dc/dn
	  sumdadn = 0
!
! Determine crack growth parameters
	  IF (a.LT.((cmax - c)/5.D0*1.D0+c)) THEN
		 A1 = Mat(8,loc) * Mat(10,loc)
		 A2 = Mat(9,loc) * Mat(10,loc)
	  ELSE IF (a.LT.((cmax - c)/5.D0*2.D0+c)) THEN
		 A1 = Mat(8,loc) * Mat(11,loc)
		 A2 = Mat(9,loc) * Mat(11,loc)
	  ELSE IF (a.LT.((cmax - c)/5.D0*3.D0+c)) THEN
		 A1 = Mat(8,loc) * Mat(12,loc)
		 A2 = Mat(9,loc) * Mat(12,loc)
	  ELSE IF (a.LT.((cmax - c)/5.D0*4.D0+c)) THEN
		 A1 = Mat(8,loc) * Mat(13,loc)
		 A2 = Mat(9,loc) * Mat(13,loc)
	  ELSE
		 A1 = Mat(8,loc) * Mat(14,loc)
		 A2 = Mat(9,loc) * Mat(14,loc)
	  END IF
	  DKtr = (A2/A1)**(1.D0/(M1-M2))
!      IF (DKtr.GE.DK0) THEN
!         Astar = A1a / (DK0**(mstar-m1))
!      ELSE
!         Astar = A2a / (DK0**(mstar-m2))
!      END IF
	  DO j = 1, nccmax
	     iC = iC + 1
	     nc = Spect(1,j,iper)
		 DSm = Spect(2,j,iper)
		 DSb = Spect(3,j,iper)
		 DKa = SMult * Csr(iC, loc) * DSQRT(pi * a) * (DSm * Yma * (Bnom/B)**p + DSb * Yba * (Bnom/B)**2.D0)
		 CALL CrackGrowthLaw (A1, m1, A2, m2, DK0, K1C, fy, fu, (Spm+Spb), DKa, a, dadn)
         sumdadn = sumdadn + dadn * nc
!!	     WRITE(17,'(4E15.3, I5, 4E15.3)') a, nc, DSm, SMult, iC, Csr(iC,loc), DKa, dadN
		 IF (iC.EQ.nmaxCsr) iC = 0
	  END DO
!
! Update n, a and c
	  dadn = sumdadn / ntot
!!	  WRITE(17,'(E15.3)') dadn
!
! Procedure if all stress intensity ranges are lower than the threshold
	  IF ((dadn .EQ. 0.D0)) THEN
	     n = ntot
	     fail = .false.
		 GOTO 99
	  END IF
!
! Continue the normal procedure
	  dn = a * ma / (dadn + dadn0) * 2.D0
	  IF (a.NE.a0) n = n + dn
!
! Check if the calculated number of cycles is larger than the applied number of cycles
	  IF (n.GE.ntot) THEN
		 a = a + (ntot - n + dn) * (dadn + dadn0) / 2.D0
		 n = ntot
		 fail = .false.
!	     WRITE(16,'(I5, 4E15.3, 3F10.3)') iper, a, ad, n, ntot
		 GOTO 99
	  END IF
	  dadn0 = dadn
      a = a * (1.D0 + ma)
!
! Evaluation of fracture
      CALL SIGrefC ( a, c, B, W, Spm, Spb, Sref)
      Lr = Sref / fy
      CALL FactorKr ( a, c, B, W, thetac, Yma, Yba, Spm, Spb, Ssm, Ssb, K1C, Lr, Kr)
	  fKrLr = FAD - DSQRT(Kr*Kr+Lr*Lr)
!	  WRITE(16,'(I5, 4E15.3, 3F10.3)') iper, a, ad, n, ntot, Kr, Lr, fKrLr
      IF (fKrLr .LT. 0.D0) GOTO 99
!
! Check if the calculated crack size is larger than the detectable crack size
	  IF (a.GE.ad) THEN
		 det = .true.
		 fail = .false.
		 GOTO 99
	  END IF
   END DO
!
! Total number of cycles of all periods
99 CONTINUE
   nm = nm + n
!
! Upade geometry
   GeoInOut(1) = a
   GeoInOut(2) = c
   GeoInOut(3) = tt * 1.D0
!
! Einde subroutine
   RETURN
   END