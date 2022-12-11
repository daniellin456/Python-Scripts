MODULE CONSTANTS

IMPLICIT NONE

REAL( KIND = 8 ), ALLOCATABLE, DIMENSION(:,:) :: Dens_Table, Ener_Table, DensEint2Pres_Table
REAL( KIND = 8 ), ALLOCATABLE, DIMENSION(:,:) :: S_Table, DensEint2Cs_Table, DensEint2Temp_Table
REAL( KIND = 8 ), ALLOCATABLE, DIMENSION(:,:) :: DensTemp2Eint_Table, DensPres2Eint_Table

REAL( KIND = 8 ), PARAMETER :: kB      = 1.3806200d-16
REAL( KIND = 8 ), PARAMETER :: mH      = 1.6600000d-24
REAL( KIND = 8 ), PARAMETER :: mu_gas  = 2.31
REAL( KIND = 8 ), PARAMETER :: gamma   = 1.66666667
REAL( KIND = 8 ), PARAMETER :: Grav    = 6.67e8
REAL( KIND = 8 ), PARAMETER :: scale_l = 3.08e18
REAL( KIND = 8 ), PARAMETER :: scale_d = mu_gas*mH
REAL( KIND = 8 ), PARAMETER :: scale_t = 1.0 / sqrt(Grav*scale_d)
REAL( KIND = 8 ), PARAMETER :: scale_v = scale_l / scale_t
REAL( KIND = 8 ) :: RhoMin
REAL( KIND = 8 ) :: RhoMax
REAL( KIND = 8 ) :: Emax
REAL( KIND = 8 ) :: Emin
REAL( KIND = 8 ) :: yHe
REAL( KIND = 8 ) :: Tmax = 1.0d5
REAL( KIND = 8 ) :: Tmin = 3.0d0
REAL( KIND = 8 ) :: Pmax = 1e15
REAL( KIND = 8 ) :: Pmin = 1e8

INTEGER  :: nRho
INTEGER  :: nEnergy
INTEGER  :: nTemp
INTEGER  :: nPres

CONTAINS

    
END MODULE CONSTANTS