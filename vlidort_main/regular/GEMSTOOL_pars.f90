   MODULE GEMSTOOL_PARS_m

!  10/18/16. Added SIF control (new sub-structure GEMSTOOL_SurfaceLeaving)
!            Expanded to include Water-leaving inputs (10/25/16 - see below)

!  10/24/16. Added new substructure for using BRDF options - this is the standard
!            BRDF multi-kernel VLIDORT system, should not be used for Ocean Color.

!  10/25/16. Added Control for Water-leaving inputs in sub-structure GEMSTOOL_SurfaceLeaving
!            This also includes use of an optional BRDF Cox-Munk glitter function.

!  10/26/16. Added new substructure for using Instrument band data. Options
!             limited to Himawari-8 at the moment, and only for UVN tool

   IMPLICIT NONE

!  Real number type definitions

   INTEGER, PARAMETER :: GEMSTOOL_SPKIND = SELECTED_REAL_KIND(6)
   INTEGER, PARAMETER :: GEMSTOOL_DPKIND = SELECTED_REAL_KIND(15)
   INTEGER, PARAMETER :: GTPK = GEMSTOOL_DPKIND

!  Version number
!  ==============

   CHARACTER(Len=3), PARAMETER :: GEMSTOOL_VERSION_NUMBER = '2.1'

!  Basic dimensions
!  ================

!  wavelengths

  ! INTEGER, PARAMETER :: GT_maxwav     = 300000
   INTEGER, PARAMETER :: GT_maxwav     = 30000
!   INTEGER, PARAMETER :: GT_maxwav     = 501     UVN Testing 315-335
!   INTEGER, PARAMETER :: GT_maxwav     = 2501
!   INTEGER, PARAMETER :: GT_maxwav     = 301       ! Original value
   INTEGER, PARAMETER :: GT_maxwavConv = 2         ! new 10/27/16

!  10/26/16. Satellite data. Maximum No. of instruments, instrument wavelengths
!   Inst_maxlambdas = maximum channel wavelengths you think you will need.

   INTEGER, PARAMETER :: GT_maxInstruments = 1
   integer, parameter :: GT_maxInstlambdas = 1 !5001

!  Maximum number of geometries, angles

   INTEGER, PARAMETER :: GT_maxGeometries = 1!11 ! Formerly, hard-wired in the "Closure" type structure
   INTEGER, PARAMETER :: GT_maxVzas = 1!11
   INTEGER, PARAMETER :: GT_maxSzas = 1!11
!   INTEGER, PARAMETER :: GT_maxGeometries = 1 ! Formerly, hard-wired in the "Closure" type structure
!   INTEGER, PARAMETER :: GT_maxVzas = 1
!   INTEGER, PARAMETER :: GT_maxSzas = 1

!  Geophysics dimensioning
!  -----------------------

!  Gases and Layering (Increase layers for New Aerosol Loading

!   integer, parameter :: GT_maxlayers      = 25     ! VLIDORT Value
   integer, parameter :: GT_maxlayers      = 99    ! VLIDORT Value
   integer, parameter :: GT_maxgases       = 5
   integer, parameter :: GT_maxfine        = 2      ! At least 20    May not need this

!  10/19/16.
!     Number of aerosol modes (Local parameter setting here - could be changed later)
!     2 = Bimodal. Should be the same as in the AerProperties Modules.

   integer, parameter :: GT_maxAerRegimes  = 2
   integer, parameter :: GT_maxAerModes    = 2

!  Aerosol and Cloud Moment Dimensions
!   ---> These should be the same as each other
!   ---> These should be the same as MAXMOMENTS_INPUT in Vlidort_pars.f90

   integer, parameter :: GT_maxmoments  = 2!500    ! VLIDORT Value
   integer, parameter :: GT_maxaermoms  = 2!500
   integer, parameter :: GT_maxcldmoms  = 2!500

!  Surface

   integer, parameter :: GT_maxClosureBands = 20 ! Formerly, hard-wired in the "Closure" type structure                
   integer, parameter :: GT_maxAlbCoeffs    = 5  ! Formerly, hard-wired in the "Closure" type structure                
   integer, parameter :: GT_maxKernels      = 4  ! Maximum number of BRDF Kernels
   integer, parameter :: GT_maxKernelpars   = 3  ! Maximum number of BRDF Kernel parameters

!  Messages

   integer, parameter :: GT_maxmessages    = 25

!  Linearization settings.
!     GT_maxaerwfs = GT_maxAerRegimes * (3 + 6*GT_maxAerModes)

!  3/10/20. Rob addition. Define GT_maxRaywfs. Set this to 1

   integer, parameter :: GT_maxsurfacewfs = 1!4
   integer, parameter :: GT_maxgaswfs     = 1
   integer, parameter :: GT_maxRaywfs     = 1
   integer, parameter :: GT_maxaerwfs     = 1!3

!  3/10/20. Rob change. Add to GT_maxRaywfs to GT_maxatmoswfs

   integer, parameter :: GT_maxatmoswfs   = GT_maxaerwfs + GT_maxgaswfs + GT_maxRaywfs
!   integer, parameter :: GT_maxatmoswfs   = GT_maxaerwfs + GT_maxgaswfs

!  RTM Settings. MUST BE SAME AS VLIDORT

   integer, parameter :: GT_maxstreams  = 16
   integer, parameter :: GT_max2streams = 2 * GT_maxstreams

!  PCA settings. T_maxbins is increased to 100 for V4 binning (J. Bak)

   integer, parameter :: GT_maxbins = 100
!   integer, parameter :: GT_maxbins = 9
   integer, parameter :: GT_maxeofs = 5
   integer, parameter :: GT_maxeofs2p1 = 2 * GT_maxeofs + 1

!  Constants
!  =========

!  Numbers

   REAL(GTPK), parameter :: GTZERO = 0.0_GTPK
   REAL(GTPK), parameter :: GTONE  = 1.0_GTPK
   REAL(GTPK), parameter :: GTTWO  = 2.0_GTPK
   REAL(GTPK), parameter :: GT100  = 100.0_GTPK

!  Critical size for aerosol moment

   real(GTPK)   , parameter      :: momsize_cutoff = 0.001_GTPK

!  CO2 mixing ratio. @@@@@@ Updated to 400, 

!   real(GTPK)   , PARAMETER      :: CO2_PPMV_MIXRATIO = 385.0d0
!   real(GTPK)   , PARAMETER      :: CO2_PPMV_MIXRATIO = 395.0_GTPK   ! 2/28/13
   real(GTPK)   , PARAMETER      :: CO2_PPMV_MIXRATIO = 405.0_GTPK   ! 2018 Value

!  Limiting single scatter albedo
    !	This should be pre-set to coincide with the internal limit in
    !	the LIDORT models. Suggested value = 0.999999

   REAL(GTPK), parameter :: omega_lim = 0.999999_GTPK
   REAL(GTPK), parameter :: omega_smallnum = 1.0E-15_GTPK

!  Gas-law constants
!    (Loschmidt's number (particles/cm3), STP values

    real(GTPK)   , PARAMETER :: RHO_STANDARD = 2.68675e+19_GTPK
    real(GTPK)   , PARAMETER :: PZERO = 1013.25_GTPK
    real(GTPK)   , PARAMETER :: TZERO = 273.15_GTPK
    real(GTPK)   , PARAMETER :: xo2  = 0.209476_GTPK, kb = 1.381e-23_GTPK
!  End of file.

   END MODULE GEMSTOOL_PARS_m

