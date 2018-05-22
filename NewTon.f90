$ Declare
   PROGRAM NewTon
!
! Computerprogramma NewTon voor de berekening van de faalkans voor
! verschillende faalmechanismen.
! 
! TNO Bouw Juni 2005 SNH
! ----------------------------------------------------------------------
!
! Declaraties
   USE PublicPar
   IMPLICIT NONE
!
   DOUBLE PRECISION, POINTER     :: alfav(:), x(:), CM(:,:)
   DOUBLE PRECISION              :: bester, Pfster, betax
   INTEGER                       :: mechnr, bermth, nstoch, conv
   INTEGER                       :: mstoch, maxmch, maxxup
   INTEGER                       :: sw(40), swu(3), ferr, werr
!
! Parameters
   DOUBLE PRECISION              :: B, B2, W, tw, phi, rho
   DOUBLE PRECISION              :: D
   INTEGER                       :: i, j, k, joff
   DOUBLE PRECISION              :: ma0, Va0 
   DOUBLE PRECISION              :: aRa0c0
   DOUBLE PRECISION              :: fy, fu, K1C
   DOUBLE PRECISION              :: mSpm, sSpm, Spb, Ssm, Ssb, MT, ST
   DOUBLE PRECISION              :: Rfyfu, Rao, Raoco, Rb, Rtw, Rphi, Rrho, Rdaf, Rc, Rkc, Rt0, Rfad, Rko, Ra1, Ra2, Rad
   DOUBLE PRECISION              :: V1, V2, V3, V4, sC, mC, VC, sCdep, VCdep, mDaf, sDaf
   DOUBLE PRECISION              :: mlog10A1, slog10A1, mlog10A2, slog10A2
   DOUBLE PRECISION              :: mDK0, VDK0
   DOUBLE PRECISION              :: mT27J, sT27J
   DOUBLE PRECISION              :: mfy, Vfy
   DOUBLE PRECISION              :: mfu, Vfu
   DOUBLE PRECISION              :: m1, m2
   DOUBLE PRECISION              :: PoD1, PoD2, PoD3
   DOUBLE PRECISION, ALLOCATABLE :: ntot(:)      ! Number of cycles between inspections and after last inspection
   CHARACTER, ALLOCATABLE        :: sfile*100(:) ! Names of stress files
   INTEGER                       :: ncc
!
! Number of potential crack locations (number of details with same nominal loading and geometry)
   nloc   =48
!
! Geometry
   LoadType     = 1
   FlawType     = 9
   WeldGeometry = 0
   Bnom      = 6.D0    ! Plate thickness
   B2        = 0.D0    ! Thickness of 2nd plate or length of longitudinal stiffener (if any)
   W         = 590.D0  ! Plate width
   tw        = 5.D0    ! Weld throat
   phi       = 45.D0   ! Weld toe angle
   rho       = 0.5D0   ! Weld toe radius
   cmax      = 200.D0  ! (Semi) crack length at which failure is assumed if not already fractured
   p         = 1.0D0   ! Thickness exponent
   ma0       = 0.15D0  ! Mean of the initial crack depth
   Va0       = 0.66D0  ! Coefficient of variation of the initial crack depth
   aRa0c0    = 0.62D0  ! Aspect ratio of the initial crack size 
   alpha     = 1.5D0   ! Fitnet: 1.0 for plane stress to 3.0 for plane strain. Usual psteels ranging from 1.5 to 2.0
!
! Inspections
   ninsp = 2     ! number of inspections
   ALLOCATE (ntot(ninsp+1))
   ALLOCATE (sfile(ninsp+1))
   ALLOCATE (nd(ninsp))
   ttd       = 1       ! 1 = through thickness crack can be detected, 0 = surface crack can be detected
   ! PoD is dependent on the working conditions: easy access for visual inspection, normal working conditions fro ACFM/MPI, no differentiation for UT)
   PoD1       = 43.76D0  ! scale parameter lambda: 38.76D0 + 5.D0 for visual inspection (between easy and moderate  access) and 0.3638D0 + 0.8D0for ACFM (ground welds)
   PoD2       = 0.45D0   ! shape parameter alfa:    0.45D0 for visual inspection (between easy and moderate  access) and 0.48D0 for ACFM (ground welds)
   PoD3       = 5.D0     ! scale parameter a*: 5.D0 for visual and 0.8D0 for ACFM
   nd(1)     = 0      ! number of detected cracks during 1st inspection
   nd(2)     = 0       ! number of detected cracks during 2nd inspection
!
! Information on stress range histograms. To be provided (ninsp + 1) times
   R         = -3.0D0       ! Stress ratio (=smin/smax)
   ntot(1)   = 6.118D7      ! Number of cycles at 1st inspection  
   ntot(2)   = 7.761D6      ! Number of cycles between 1st and 2nd inspection  
   ntot(3)   = 2.394D7      ! Number of cycles after last inspection  
   sfile(1)  = 'spectrum111a.txt' ! File name of stress histogram before 1st inspection
   sfile(2)  = 'spectrum111b.txt' ! File name of stress histogram between 1st and 2nd inspection
   sfile(3)  = 'spectrum111c.txt' ! File name of stress histogram after last inspection
!
! Material properties
   fynom     = 235.D0
   mT27J     = 15.D0
   sT27J     = 15.D0
   law       = 1       ! Crack growth law 1 = Forman Mettu, 2 = BS 7910 air
!
! ULS loading
   mSpm      = 69.3D0   ! primary tensile stress (scale parameter for Weibull distribution)
   sSpm      = 6.46D0  ! primary tensile stress (form parameter for Weibull distribution)
   Spb       = 0.D0	   ! primary bending stress
   Ssm       = 0.5*fynom  ! secondary bending stress 
   MT        = 10.D0   ! mean temperature
   ST        = 8.D0	   ! standard deviation temperature
!
! Model uncertainty and dynamic amplification
   V1 = DSQRT(0.12D0*0.12D0 + 0.12D0*0.12D0) ! coefficient of variation for the uncertainty factor for the stress range (different per stress range)
   V2 = 0.15D0 ! coefficient of variation for the accuracy of the model prediction over locations (same for all stress ranges)
   V3 = 0.2D0  ! coefficient of variation for the uncertainty factor for the SIF (measurements/model)
   nmaxCsr = 9 ! max. no. lines related to V1
   sC = DSQRT (DLOG(V2*V2+1.D0) + DLOG(V3*V3+1.D0))
   mC = DLOG(1.D0)-(DLOG(V2*V2+1.D0)/2.D0) + DLOG(1.D0)-(DLOG(V3*V3+1.D0)/2.D0)
   VC = sC/mC
   sCdep = DSQRT (DLOG(V3*V3+1.D0))
   VCdep = sCdep/mC
   mDaf = 1.042D0 * 1.0992D0
   sDaf = 0.D0
!
! Correlation values
   Rfyfu = 0.75D0
   Rao = 0.D0
   Raoco = 0.D0
   Rb = 0.D0
   Rtw = 0.5D0
   Rphi = 0.7D0
   Rrho = 0.2D0
   Rdaf = 0.8D0
   Rc = 1.D0 - ((VCdep*VCdep) / (VC*VC))
   Rkc = 0.0D0
   Rt0 = 0.D0
   Rfad = 0.D0
   Ra1 = 0.65D0
   Ra2 = 0.45D0
   Rko = 0.75D0
   Rad = 0.7D0 ! For visual inspection
!   Rad = 0.3D0 ! For ACFM
!
! Determine the material strength properties
   IF (fynom.EQ.235D0) THEN
      D = 1.5D0
   ELSE
      D = 1.1D0
   END IF
   Vfy       = 0.07D0
   Vfu       = 0.04D0
   mfy       = fynom / (1.D0 - 2.D0 * Vfy)
   mfu       = mfy * D
!
! Determine the crack growth parameters
   IF (law.EQ.1) THEN
      m1       = 3.D0
      mlog10A1 = -12.447D0
	  IF (R.GE.0.5D0) THEN
	     slog10A1 = 0.171D0
	  ELSE
	     slog10A1 = 0.115D0
	  END IF
	  m2       = m1
	  mlog10A2 = mlog10A1
	  slog10A2 = slog10A1
	  mDK0     = 158.D0
	  VDK0     = 0.1D0
   ELSEIF ((law.EQ.2).AND.(R.GE.0.5D0)) THEN
      m1       = 5.1
      mlog10A1 = -17.32D0
	  slog10A1 = 0.320D0
      m2       = 2.88
	  mlog10A2 = -12.24D0
	  slog10A2 = 0.171D0
	  mDK0     = 63.D0
	  VDK0     = 0.1D0
   ELSEIF ((law.EQ.2).AND.(R.LT.0.5D0)) THEN
      m1       = 8.16
      mlog10A1 = -25.92D0
	  slog10A1 = 0.279D0
      m2       = 2.88
	  mlog10A2 = -12.41D0
	  slog10A2 = 0.115D0
	  mDK0     = MIN (170.D0, 170.D0-214.D0*R)
	  VDK0     = 0.1D0
   END IF
!
! Read stress range histograms in subsequent periods
   nccmax = 0.D0
   ALLOCATE (Spect(3,100,ninsp+1))
   DO j = 1, ninsp+1
      CALL ReadStressRanges (j, ntot(j), sfile(j), ncc)
	  WRITE (*,*) ntot(j), sfile(j), ncc
      IF (ncc.GT.nccmax) nccmax = ncc
   END DO
!
! Initialise
   mstoch = 12 + ninsp + nloc * (20+ninsp+nmaxCsr)
   maxmch = 1
   maxxup = 0
   ALLOCATE (alfav(mstoch))
   ALLOCATE (x(mstoch))
   ALLOCATE (CM(mstoch,mstoch))
   OPEN (8,file='outPfAlpha.dat',status='unknown')
   OPEN (9,file='outPf.dat',status='unknown')
   OPEN (10,file='outHit.dat',status='unknown')
   OPEN (11,file='outFalse.dat',status='unknown')
   OPEN (12,file='outAvHit.dat',status='unknown')
   OPEN (13,file='outAvFalse.dat',status='unknown')
   OPEN (14,file='outAllHit.dat',status='unknown')
   OPEN (15,file='outAllFalse.dat',status='unknown')
   WRITE (10,*)'Calculations that match the inspection results.'
   WRITE (10,*)'_______________________________________________________________________________'
   WRITE (10,*)' # cycles to | as previous, | failure # cycles/ | as previous, |detected cracks'
   WRITE (10,*)'   failure   | period after | applied # cycles  | period after |in subsequent'
   WRITE (10,*)'             |  last insp.  |                   |  last insp.  |inspections'
   WRITE (10,*)'_____________|______________|___________________|______________|_______________'
   WRITE (11,*)'Calculations that do not match the inspection results.'
   WRITE (11,*)'_______________________________________________________________________________'
   WRITE (11,*)' # cycles to | as previous, | failure # cycles/ | as previous, |detected cracks'
   WRITE (11,*)'   failure   | period after | applied # cycles  | period after |in subsequent'
   WRITE (11,*)'             |  last insp.  |                   |  last insp.  |inspections'
   WRITE (11,*)'_____________|______________|___________________|______________|_______________'
   DO i = 12, 15
      WRITE (i,*) 'life     ,life/krit,ndet,  Spm    ,  Spb    ,  Ssm    ,  Ssb    ,  SMult  ,   fy    ,   fu    ,   K1C   ,  FAD    ,   m1    ,   m2    ,   DK0   ,   A1    ,   A2    ,   VA1   ,   VA2   ,   VA3   ,   VA4   ,   VA5   ,   a0    ,   c0    ,   W     ,   B     ,   tw    ,   L     ,   phi   ,   rho   , ad'
   END DO
!
! Set correlation matrix
   CALL VariablesCorrMatrix (Rfyfu, Rao, Raoco, Rb, Rtw, Rphi, Rrho, Rdaf, Rc, Rkc, Rt0, Rfad, Rko, Ra1, Ra2, Rad, mstoch, CM)
!
! Initialise probabilistic database
   CALL PrbIni (mstoch, maxmch, maxxup)
!
! IO initialiseren
!      swu = 0
!	   sw  = 0
!      CALL SetIo (1, 1, 1, swu, sw)
!
! Set numerical control values
   CALL SetNd (1, 50, 1, .15D0, 0.01D0, 0.01D0, 0.3D0, 2, 0, 1000000, 100, -5.D0, 5.D0)
!  CALL SetNd (1, 50, 1, .15D0, 0.01D0, 0.01D0, 0.3D0, 2, 0, 1, 1, -5.D0, 5.D0)
!
! Mechanism and calculation number
   mechnr = 1
   bermth = 3 ! 1 is FORM, 3 = crude MC
!
! Set correlation among all the variables
   CALL SetCor(CM, mstoch, ferr, werr)
!
! Initialise random variables
   CALL SetSto (1, 'W',      1, 0, W, 0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (2, 'm1',     1, 0, m1,  0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (3, 'm2',     1, 0, m2,  0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (4, 'Spm',    1, 9, mSpm, sSpm, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (5, 'Spb',    1, 0, Spb, 0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (6, 'Ssm',    1, 0, Ssm, 0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (7, 'Ssb',    1, 0, Ssb, 0.D0, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (8, 'Temp',   1, 2, mT, sT, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (9, 'fy',     1, 18, mfy, mfy*Vfy, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (10,'fu',     1, 18, mfu, mfy*Vfu, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   CALL SetSto (11,'T27J',   1, 2, mT27J, sT27J, 0.D0, 0.D0, ferr)
   IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   DO j = 1, ninsp+1
      CALL SetSto (11+j, 'ntot', 1, 0, ntot(j), 0.D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
   END DO
   DO j = 1, nloc
      joff = (j-1) * (20+ninsp+nmaxCsr) + ninsp + 13
      CALL SetSto (joff+0,  'a0',     1, 18, ma0,  Va0*ma0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+1,  'aRa0c0', 1, 18, aRa0c0,  0.4D0*aRa0c0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+2,  'B',      1,  2, Bnom, Bnom*0.03D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+3,  'B2',     1,  2, B2, B2*0.03D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+4,  'tw',     1,  2, tw, tw*0.15D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+5,  'phi',    1,  2, phi, phi*0.10D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+6,  'rho',    1, 18, rho, rho*0.15D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+7,  'daf',    1, 18, mDaf, sDaf, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+8,  'LnCunc',   1, 2, mC, sC, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+9,  'Kmatnor',1,  2, 0.D0, 1.D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+10, 'T0',     1,  2, 18.D0, 15.D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+11, 'FAD',    1, 18, 1.7D0, 0.35D0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+12, 'DK0',    1, 18, mDK0,  VDK0*mDK0, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+13, 'A1',     1,  2, mlog10A1,  slog10A1, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      CALL SetSto (joff+14, 'A2',     1,  2, mlog10A2,  slog10A2, 0.D0, 0.D0, ferr)
      IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      DO i = 1, 5
	     CALL SetSto (joff+14+i, 'VA', 1, 2, 1.D0, 0.1D0, 0.D0, 0.D0, ferr)
         IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
	  END DO
	  DO i = 1, ninsp
         CALL SetSto (joff+19+i, 'ad', 1, 9, PoD1, PoD2, PoD3, 0.D0, ferr)
         IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
	  END DO
	  DO i = 1, nmaxCsr
	     CALL SetSto (joff+19+ninsp+i, 'Csr', 1, 18, 1.D0, V1*1.D0, 0.D0, 0.D0, ferr)
         IF (ferr .NE. 0) WRITE(*,'('' Fatale fout :'',I5)') ferr
      END DO
   END DO
!
! Choice of random variables per mechanism      
   DO i = 1, mstoch
      CALL SetStM(i, 1)
   END DO
!
! Execute calculation per mechanism
   CALL BerMch (mechnr, bermth, 'Mech  1', nstoch, alfav, &
                   bester, Pfster, x, conv, ferr, werr)
   IF (werr .NE. 0) THEN
       WRITE(*,'('' Waarschuwing:'',I5)') werr
       WRITE(8,'('' Waarschuwing:'',I5)') werr
       WRITE(9,'('' Waarschuwing:'',I5)') werr
   END IF
   IF (ferr .NE. 0) THEN
       WRITE(*,'('' Fatale fout :'',I5)') ferr
       WRITE(8,'('' Fatale fout :'',I5)') ferr
       WRITE(9,'('' Fatale fout :'',I5)') ferr
   END IF
!
! Write results
   WRITE(*,'('' Beta = ''ES10.2)') bester
   WRITE(*,'('' Pf   = ''ES10.2)') Pfster
   DO i = 1, 20
      WRITE(*,'(I3,F8.3,ES10.2)')i, alfav(i), x(i)
   END DO
   WRITE(9,'(2ES10.3)') bester, Pfster
   WRITE(8,'('' Beta = ''ES10.2)') bester
   WRITE(8,'('' Pf   = ''ES10.2)') Pfster
   DO i = 1, 20
      WRITE(8,'(I3,F10.5,ES15.7)')i, alfav(i), x(i)
   END DO
   CLOSE(8)
   CLOSE(9)
   CLOSE(10)
   CLOSE(11)
   CLOSE(12)
   CLOSE(13)
   CLOSE(14)
   CLOSE(15)
!
! End of program 
   END