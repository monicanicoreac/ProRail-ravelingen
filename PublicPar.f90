MODULE PublicPar
!
! Stress histogram parameters
   DOUBLE PRECISION, ALLOCATABLE :: Spect(:,:,:) ! Stress histograms
   INTEGER                       :: nccmax       ! size of 2nd column of stress histograms (number of rows in stress range input files)
   DOUBLE PRECISION              :: R            ! Stress ratio
   INTEGER                       :: nmaxCsr      ! Maximum number of stress range uncertainty values
!
! Inspection parameters
   INTEGER                       :: ninsp        ! number of inspections
   INTEGER                       :: ttd          ! 0 = surface flaws can be detected
					    		                 ! 1 = only through thickness flaws can be detected
   INTEGER, ALLOCATABLE          :: nd(:)        ! number of detected cracks during each inspection
!
! Number of detail locations
   INTEGER                       :: nloc         ! number of potential crack locations
!
! Material parameters
   DOUBLE PRECISION              :: fynom        ! nominal yield stress
   INTEGER                       :: law          ! Crack growth equation
                                                 ! 1 = Forman Mettu
												 ! 2 = BS 7910:2013 in air
!
! Geometry parameters
   INTEGER                       :: LoadType     ! 1 = Membrane loading
                                                 ! 2 = Bending loading
												 ! 3 = Both membrane and bending
   INTEGER                       :: FlawType     ! 1 = Through-thickness flaws in plates
                                                 ! 2 = Surface flaws in plates
												 ! 3 = Through-width flaws in plates
									             ! 4 = semi-circular surface flaw in round bar
												 ! 5 = Corner flaws in plates
												 ! 6 = Semi-circular surface flaw in bolts, solution 2
												 ! 7 = Edge flaws in plates
												 ! 8 = Through-thickness deck plate crack (orthotropic bridge deck)
												 ! 9 = Cope longitudinal girger (raveling) ProRail
   INTEGER                       :: WeldGeometry ! 0 = No weld
                                                 ! 1 = According to fig. M.23
                                                 ! 2 = According to fig. M.24
   DOUBLE PRECISION              :: Bnom         ! nominal plate thickness
   DOUBLE PRECISION              :: p            ! exponent for the dependency of stress on plate thickness (1 for membrane, 2 for bending)
   DOUBLE PRECISION              :: cmax         ! Semi crack length at which end-of-life is assumed if not already fractured
   DOUBLE PRECISION              :: alpha        ! Constraint constant
!   INTEGER, PARAMETER :: ndof = 6
!
END MODULE PublicPar