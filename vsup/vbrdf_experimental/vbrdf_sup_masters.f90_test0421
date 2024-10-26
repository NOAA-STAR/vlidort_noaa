
! ###############################################################
! #                                                             #
! #                       VLIDORT_2p8p3                         #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Authors :     Robert. J. D. Spurr (1)                      #
! #                Matt Christi                                 #
! #                                                             #
! #  Address (1) : RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                                                             #
! #  Tel:          (617) 492 1183                               #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  This Version :   VLIDORT_2p8p3                             #
! #  Release Date :   31 March 2021                             #
! #                                                             #
! #  Previous VLIDORT Versions under Standard GPL 3.0:          #
! #  ------------------------------------------------           #
! #                                                             #
! #      2.7   F90, released        August 2014                 #
! #      2.8   F90, released        May    2017                 #
! #      2.8.1 F90, released        August 2019                 # 
! #      2.8.2 F90, limited release May    2020                 # 
! #                                                             #
! #  Features Summary of Recent VLIDORT Versions:               #
! #  -------------------------------------------                #
! #                                                             #
! #      NEW: TOTAL COLUMN JACOBIANS         (2.4)              #
! #      NEW: BPDF Land-surface KERNELS      (2.4R)             #
! #      NEW: Thermal Emission Treatment     (2.4RT)            #
! #      Consolidated BRDF treatment         (2.4RTC)           #
! #      f77/f90 Release                     (2.5)              #
! #      External SS / New I/O Structures    (2.6)              #
! #                                                             #
! #      SURFACE-LEAVING / BRDF-SCALING      (2.7)              #
! #      TAYLOR Series / OMP THREADSAFE      (2.7)              #
! #      New Water-Leaving Treatment         (2.8)              #
! #      LBBF & BRDF-Telescoping, enabled    (2.8)              #
! #      Several Performance Enhancements    (2.8)              #
! #      Water-leaving coupled code          (2.8.1)            #
! #      Planetary problem, media properties (2.8.1)            #
! #      Doublet geometry post-processing    (2.8.2)            #
! #      Reduction zeroing, dynamic memory   (2.8.2)            #
! #                                                             #
! #  Features Summary of This VLIDORT Version                   #
! #  ----------------------------------------                   #
! #                                                             #
! #   2.8.3, released 31 March 2021.                            #
! #     ==> Green's function RT solutions (Nstokes = 1 or 3)    #
! #     ==> Sphericity Corrections using MS source terms        #
! #     ==> BRDF upgrades, including new snow reflectance       #
! #     ==> SLEAVE Upgrades, extended water-leaving treatment   #
! #                                                             #
! ###############################################################

! ###################################################################
! #                                                                 #
! # This is Version 2.8.3 of the VLIDORT_2p8 software library.      #
! # This library comes with the Standard GNU General Public License,#
! # Version 3.0, 29 June 2007. Please read this license carefully.  #
! #                                                                 #
! #      VLIDORT Copyright (c) 2003-2021.                           #
! #          Robert Spurr, RT Solutions, Inc.                       #
! #          9 Channing Street, Cambridge, MA 02138, USA.           #
! #                                                                 #
! # This file is part of VLIDORT_2p8p3 ( Version 2.8.3 )            #
! #                                                                 #
! # VLIDORT_2p8p3 is free software: you can redistribute it         #
! # and/or modify it under the terms of the Standard GNU GPL        #
! # (General Public License) as published by the Free Software      #
! # Foundation, either version 3.0 of the License, or any           #
! # later version.                                                  #
! #                                                                 #
! # VLIDORT_2p8p3 is distributed in the hope that it will be        #
! # useful, but WITHOUT ANY WARRANTY; without even the implied      #
! # warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR         #
! # PURPOSE. See the Standard GNU General Public License (GPL)      #
! # for more details.                                               #
! #                                                                 #
! # You should have received a copy of the Standard GNU General     #
! # Public License (GPL) Version 3.0, along with the VLIDORT_2p8p3  #
! # code package. If not, see <http://www.gnu.org/licenses/>.       #
! #                                                                 #
! ###################################################################

! ###############################################################
! #                                                             #
! # Subroutines in this Module                                  #
! #                                                             #
! #            VBRDF_INPUTMASTER                                #
! #            VBRDF_MAINMASTER                                 #
! #                                                             #
! ###############################################################

!  5/5/22. Experimental Version for NASA-LARC PCRTM team
!    ** VBRDF_INPUTMASTER is not changed
!    ** VBRDF_MAINMASTER has major re-engineering
!        -- now only does kernel-scaling and Fourier components
!        -- Assumes Kernel cores already computed and stoked in "core" type structure
!        -- Also uses a bookkeeping type structure ("book") with saved geometrical quantities

      MODULE vbrdf_sup_masters_m

      PRIVATE
      PUBLIC :: VBRDF_INPUTMASTER, VBRDF_MAINMASTER

      CONTAINS

      SUBROUTINE VBRDF_INPUTMASTER ( &
        FILNAM,               & !Input
        VBRDF_Sup_In,         & !Outputs
        VBRDF_Sup_InputStatus ) !Outputs

!  Input routine for BRDF program

!  Version 2.6 notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VBRDF Upgrades for Version 2.7
!  ------------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 02 April 2014, Revised, 14-15 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the VBRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the VBRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the VSLEAVE
!  supplement, we have given the VBRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the VSLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the VSLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for VBRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrades for Version 2.8
!  ------------------------

!  Presence of 2 new kernels, 

!  RTK_HOTSPOT : This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.
!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)
!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)
!   -- Ross-Thick Hotspot Kernel from

!  MODFRESNEL. This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.
!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

!  4 parameters now allowed in Kernels

!  Patch Overhaul for Version 2.8.1. RobFix 11/8/19
!    - Fourier routine has one less argument
!    - Introduce reflectivity Mask for proper identification of Matrix elements

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel. First introduced to VLIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  1/31/21. Version 2.8.3. Add DOUBLET_GEOMETRY option, throughout code

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Scalar and vector amplitudes given separately, controlled with flag
!    -- remove calculation of USER_BRDF_F_0. Flag for this.
!    -- reduced kernel for SNOW BRDF

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Basic model is the ART (Asymptotic radiative transfer) model of Kokhanovsky and Breon
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      USE VLIDORT_PARS_m

      USE VBRDF_FINDPAR_m

      USE vbrdf_sup_inputs_def_m
      USE vbrdf_sup_outputs_def_m

!  Implicit none

      IMPLICIT NONE

!  Arguments
!  ---------

      CHARACTER (LEN=*), INTENT(IN) :: FILNAM

      TYPE(VBRDF_Sup_inputs), INTENT(OUT) :: VBRDF_Sup_In

      TYPE(VBRDF_Input_Exception_Handling), INTENT(OUT) :: VBRDF_Sup_InputStatus

!  Local variables
!  ---------------

!  Stream angle flag

      LOGICAL ::          DO_USER_STREAMS

!  BRDF surface flag
!    ---> Really should be true here

      LOGICAL ::          DO_BRDF_SURFACE

!  Surface emission

      LOGICAL ::          DO_SURFACE_EMISSION

!  Solar sources + Geometry flags
!    -- 1/31/21. Version 2.8.3. Doublet geometry option added

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Number of Stokes components

      INTEGER ::          NSTOKES

!  Number and index-list and names of bidirectional functions

      INTEGER ::            N_BRDF_KERNELS
      INTEGER ::            WHICH_BRDF ( MAX_BRDF_KERNELS )
      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER ::          N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      DOUBLE PRECISION :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  Lambertian Surface control

      LOGICAL ::          LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors
!    2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

      LOGICAL          :: DO_SEPARATE_FACTORS ( MAX_BRDF_KERNELS ) 
      DOUBLE PRECISION :: BRDF_FACTORS        ( MAX_BRDF_KERNELS ) 
      DOUBLE PRECISION :: BRDF_VECTOR_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!  Rob Fix. 9/27/14. Added variable output

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Number of azimuth quadrature streams for BRDF

      INTEGER ::          NSTREAMS_BRDF

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL ::          DO_SHADOW_EFFECT

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL   :: DO_EXACT
!      LOGICAL   :: DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents
!     -- default should be to use direct-bounce for this term

      LOGICAL ::          DO_FOURIER_USERSUN

!  Multiple reflectance correction for Glitter kernels

      LOGICAL ::          DO_MSRCORR
      INTEGER ::          MSRCORR_ORDER
      LOGICAL ::          DO_MSRCORR_DBONLY ! Rob Fix 9/25/14, Variable renamed
      INTEGER ::          MSRCORR_NMUQUAD
      INTEGER ::          MSRCORR_NPHIQUAD

!  Number of discrete ordinate streams

      INTEGER ::          NSTREAMS

!  Local angle control

      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Angles

      DOUBLE PRECISION :: BEAM_SZAS (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS (MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES (MAX_USER_STREAMS)

!  Local Observational Geometry control and angles

      INTEGER ::          N_USER_OBSGEOMS
      DOUBLE PRECISION :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  Local Doublet Geometry control and angles
!mick mod 1/5/2021 - added these doublet vars

      INTEGER ::          N_USER_DOUBLETS
      DOUBLE PRECISION :: USER_DOUBLETS (MAX_USER_STREAMS,2)

!  New Cox-Munk Glint reflectance options (bypasses the usual Kernel system)
!  -------------------------------------------------------------------------

!  Overall flags for this option. New GCM Flag, 7/4/15

      LOGICAL   :: DO_NewCMGLINT, DO_NewGCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input wavelength in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Exception handling
!  ------------------

!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) :: ACTIONS ( 0:MAX_MESSAGES )

!  local variables
!  ===============

      CHARACTER (LEN=9), PARAMETER :: PREFIX = 'BRDFSUP -'

      LOGICAL ::            ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER ::            I, K, L, FILUNIT, NM

!  MODIS-style brdfs, some checks on the input

      LOGICAL            :: DO_ROSS, DO_LI, DO_LAMB, DO_CM, DO_RTS

!  Check list of Kernel names
!  Rob Extension 12/2/14. BPDF Kernels (replace BPDF2009)
!    DO_NewGCMGLINT flag, add new Name (#16). 7.4.15

!  New for  Version 2.8. RossThick-Hotspot and Modified-Fresnel kernels

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- Kokhanovsky and Breon, IEEE GeoScience & Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- First introduced to VLIDORT, 18 November 2020. Now 19 Kernels !!!

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter

      CHARACTER (LEN=10) :: BRDF_CHECK_NAMES ( MAXBRDF_IDX )

      BRDF_CHECK_NAMES = (/ &
                           'Lambertian', &
                           'Ross-thin ', &
                           'Ross-thick', &
                           'Li-sparse ', &
                           'Li-dense  ', &
                           'Hapke     ', &
                           'Roujean   ', &
                           'Rahman    ', &
                           'Cox-Munk  ', &
                           'GissCoxMnk', &
                           'GCMcomplex', &
                           'BPDF-Soil ', &
                           'BPDF-Vegn ', &
                           'BPDF-NDVI ', &
                           'NewCMGlint', &
                           'NewGCMGlit', &
                           'RtkHotSpot', &
                           'ModFresnel', &
                           'SnowModel ', &
                           'SnowRTS   ' /)

!  Initialize Exception handling

      STATUS = VLIDORT_SUCCESS

      MESSAGES(1:MAX_MESSAGES) = ' '
      ACTIONS (1:MAX_MESSAGES) = ' '

      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Read of VLIDORT Input file'
      ACTIONS(0)      = 'No Action required for this Task'

!  Local error handling initialization

      ERROR  = .FALSE.
      NM     = NMESSAGES

!  Open file

      FILUNIT = VLIDORT_INUNIT
      OPEN(VLIDORT_INUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

!  Initialize control and angles
!  =============================

!  1/31/21. Version 2.8.3. Add doublet geometry option
!mick mod 1/5/2021 - modified some initializations to make more
!                    vector-like as in the VSLEAVE sup
!                  - initialize N_USER_DOUBLETS & USER_DOUBLETS

!  Control flags

      DO_USER_STREAMS     = .FALSE.
      DO_SOLAR_SOURCES    = .FALSE.  !@@ New line
      DO_USER_OBSGEOMS    = .FALSE.  !@@ New line
      DO_DOUBLET_GEOMETRY = .FALSE.  !@@ New line 1/31/21. Version 2.8.3
      DO_SURFACE_EMISSION = .FALSE.

!  Integer control

      NSTOKES  = 0
      NSTREAMS = 0

!  Conventional geometry inputs

      NBEAMS         = 0
      N_USER_STREAMS = 0
      N_USER_RELAZMS = 0
      BEAM_SZAS      = ZERO
      USER_ANGLES    = ZERO
      USER_RELAZMS   = ZERO

! Observational Geometry

      N_USER_OBSGEOMS = 0
      USER_OBSGEOMS   = ZERO

! Doublet geometry

      N_USER_DOUBLETS = 0
      USER_DOUBLETS   = ZERO

!  Initialize Surface stuff
!  ========================

      NSTREAMS_BRDF  = 0
      N_BRDF_KERNELS = 0

      WHICH_BRDF = 0
      BRDF_NAMES = '          '
!mick fix 2/20/2016 - initialize N_BRDF_PARAMETERS
      N_BRDF_PARAMETERS = 0

!  2/25/22. Version 2.8.5, Initialize vector factors, and control for this option.

      DO K = 1, MAX_BRDF_KERNELS
        LAMBERTIAN_KERNEL_FLAG(K) = .FALSE.
        DO_SEPARATE_FACTORS(K) = .false.
        BRDF_FACTORS(K)        = ZERO
        BRDF_VECTOR_FACTORS(K) = ZERO
        DO L = 1, MAX_BRDF_PARAMETERS
          BRDF_PARAMETERS(K,L) = ZERO
        ENDDO
      ENDDO

      DO_SHADOW_EFFECT    = .FALSE.

!  Rob Fix 9/25/14. Two variables replaced
!      DO_EXACT             = .FALSE.       !@@  New line
!      DO_EXACTONLY         = .FALSE.
      DO_DBONLY            = .FALSE.

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents

      DO_FOURIER_USERSUN = .false.

!  Multiple-scattering Glitter options

      DO_MSRCORR           = .FALSE.
      DO_MSRCORR_DBONLY    = .FALSE.
      MSRCORR_ORDER        = 0
      MSRCORR_NMUQUAD      = 0
      MSRCORR_NPHIQUAD     = 0

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!   Rob Fix 9/27/14; added output option flag

      DO_WSABSA_OUTPUT = .false.
      DO_WSA_SCALING   = .false.
      DO_BSA_SCALING   = .false.
      WSA_VALUE        = zero
      BSA_VALUE        = zero

!  NewCM options

      DO_NewCMGLINT  = .false.
      DO_NewGCMGLINT = .false.        ! New 7/4/15
      SALINITY   = zero
      WAVELENGTH = zero
      WINDSPEED  = zero
      WINDDIR    = zero
      DO_GlintShadow   = .false.
      DO_FoamOption    = .false.
      DO_FacetIsotropy = .false.

!  Number of Stokes components
!  ===========================

      PAR_STR = 'Number of Stokes vector components'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) NSTOKES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Read Angle stuff
!  ================

!  Basic control for solar sources

      PAR_STR = 'Use solar sources?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_SOLAR_SOURCES
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  User-defined Stream angle

      PAR_STR = 'Use user-defined viewing zenith angles?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
          READ (FILUNIT,*,ERR=998) DO_USER_STREAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Discrete ordinates

      PAR_STR = 'Number of half-space streams'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS
      CALL FINDPAR_ERROR (ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  All numbers are now checked against maximum dimensions

      IF ( NSTREAMS .GT. MAXSTREAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of half-space streams > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  Observational & Doublet Geometry
!  ================================
!mick mod 1/5/2021 - modifed this section similar to VSLEAVE supplement

      IF ( DO_SOLAR_SOURCES .and. DO_USER_STREAMS ) THEN

!  !@@ New flag, Observational Geometry

         PAR_STR = 'Do Observation Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_USER_OBSGEOMS
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21. Version 2.8.3. Add doublet geometry option

         PAR_STR = 'Do Doublet Geometry?'
         IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_DOUBLET_GEOMETRY
         CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  1/31/21. Version 2.8.3. Safety check on geometry options

        IF ( DO_USER_OBSGEOMS .and. DO_DOUBLET_GEOMETRY) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Not allowed to have both Observation and Doublet Geometry options'
           ACTIONS(NM)  = 'Re-set input to one or the other of these flags'
           STATUS       = VLIDORT_SERIOUS
           NMESSAGES    = NM
           GO TO 764
        ENDIF

!  End Obsgeom/Doublet control flags

      ENDIF

!  !@@ Observational Geometry control
!     ---- check not exceeding dimensioned number
!mick mod 1/5/2021 - modifed this section similar to VSLEAVE supplement

      IF ( DO_USER_OBSGEOMS ) THEN

!  Number of Observation Geometry inputs

        PAR_STR = 'Number of Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_OBSGEOMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( N_USER_OBSGEOMS .GT. MAX_USER_OBSGEOMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of ObsGeometry inputs > Maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_OBSGEOMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  Observational Geometry inputs

        PAR_STR = 'Observation Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_OBSGEOMS
             READ (FILUNIT,*,ERR=998) USER_OBSGEOMS(I,1), USER_OBSGEOMS(I,2), USER_OBSGEOMS(I,3)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of NBEAMS, N_USER_STREAMS, N_USER_RELAZMS

        NBEAMS          = N_USER_OBSGEOMS
        N_USER_STREAMS  = N_USER_OBSGEOMS
        N_USER_RELAZMS  = N_USER_OBSGEOMS
        DO_USER_STREAMS = .TRUE.

!  Automatic setting of BEAM_SZAS, USER_VZANGLES, and USER_RELAZMS

         BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
         USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
         USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)

!  Finish - go to control point

         GO TO 5667

!  End Observational geometry clause

      ENDIF

!  Solar beam input (Lattice or Doublet)
!  =====================================

!  Number of Solar zenith angles

      PAR_STR = 'Number of solar zenith angles'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NBEAMS
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

      IF ( NBEAMS .GT. MAXBEAMS ) THEN
        NM = NM + 1
        MESSAGES(NM) = 'Number of solar zenith angles > maximum dimension'
        ACTIONS(NM)  = 'Re-set input value or increase MAXBEAMS dimension'
        STATUS = VLIDORT_SERIOUS
        NMESSAGES = NM
        GO TO 764
      ENDIF

!  TOA solar zenith angle inputs

      PAR_STR = 'Solar zenith angles (degrees)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
        DO I = 1, NBEAMS
          READ (FILUNIT,*,ERR=998) BEAM_SZAS(I)
        ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Doublet Geometry input
!  ======================
!mick fix 1/5/2021 - added this doublet section (same as VSLEAVE supplement)
!  1/31/21. Version 2.8.3. Section completely new

      IF ( DO_DOUBLET_GEOMETRY ) THEN

!  Number of Doublet Geometries

        PAR_STR = 'Number of Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_DOUBLETS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_DOUBLETS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of doublet geometry inputs > maximum dimension'
          ACTIONS(NM)  = 'Re-set input or increase MAX_USER_STREAMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  Read Doublet Geometry values

        PAR_STR = 'Doublet Geometry inputs'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           DO I = 1, N_USER_DOUBLETS
             READ (FILUNIT,*,ERR=998)USER_DOUBLETS(I,1), USER_DOUBLETS(I,2)
           ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Automatic setting of N_USER_STREAMS, N_USER_RELAZMS

        N_USER_STREAMS  = N_USER_DOUBLETS
        N_USER_RELAZMS  = N_USER_DOUBLETS

!  Automatic setting of user angles

        USER_ANGLES  (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,1)
        USER_RELAZMS (1:N_USER_DOUBLETS) = USER_DOUBLETS(1:N_USER_DOUBLETS,2)

!  Finish - go to control point

        GO TO 5667

!  End Doublet geometry clause

      ENDIF

!  Lattice Geometry controls
!  =========================
!mick fix 1/5/2021 - replaced old lattice angle sections with this rewritten
!                    section (same as VSLEAVE supplement)

!  Make explicit the control here

      IF ( DO_USER_STREAMS .and. .not.DO_USER_OBSGEOMS .and. .not.DO_DOUBLET_GEOMETRY ) then

!  Number of User defined viewing zenith angles

        PAR_STR = 'Number of user-defined viewing zenith angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_USER_STREAMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_STREAMS .GT. MAX_USER_STREAMS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of viewing zenith angles > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_USER_STREAMS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  User-defined viewing zenith angles (should be positive)

        PAR_STR = 'User-defined viewing zenith angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_STREAMS
            READ (FILUNIT,*,ERR=998) USER_ANGLES(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Number of azimuth angles

        PAR_STR = 'Number of user-defined relative azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) N_USER_RELAZMS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check not exceeding dimensioned number

        IF ( N_USER_RELAZMS .GT. MAX_USER_RELAZMS ) THEN
          NM = NM + 1
          MESSAGES(NM) =  'Number of relative azimuth angles > maximum dimension'
          ACTIONS(NM)  =  'Re-set input value or increase MAX_USER_RELAZMS dimension'
          STATUS       = VLIDORT_SERIOUS
          NMESSAGES    = NM
          GO TO 764
        ENDIF

!  User defined Azimuth angles

        PAR_STR = 'User-defined relative azimuth angles (degrees)'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_USER_RELAZMS
            READ (FILUNIT,*,ERR=998) USER_RELAZMS(I)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  End lattice input geometry clause

      ENDIF

!  !@@ Continuation point for Skipping the Lattice-input or Doublet geometry angles

5667  continue

!  Surface stuff
!  =============

!  BRDF input
!  ----------

!  Basic flag

      PAR_STR = 'Do BRDF surface?'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) DO_BRDF_SURFACE
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  NewCM Flags. New for Version 2.7. 7/4/15  Extension.

      IF ( DO_BRDF_SURFACE ) THEN
        PAR_STR = 'Do NewCM Ocean BRDF reflectance?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_NewCMGLINT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        PAR_STR = 'Do NewGCM Ocean BRDF reflectance?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
              READ (FILUNIT,*,ERR=998) DO_NewGCMGLINT
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  7/4/15. Cannot have both options

      if ( DO_NewCMGLINT .and. DO_NewGCMGLINT ) then
         NM = NM + 1
         MESSAGES(NM) = 'Bad input: cannot have both NewCM and NewGCM options set'
         ACTIONS(NM)  = 'Fix the input'
         STATUS = VLIDORT_SERIOUS
         NMESSAGES = NM
         GO TO 764
      ENDIF

!  Surface emission flag. Not for NewCM or NewGCM flags, update 7/4/15

      IF ( .not.DO_NewCMGLINT .or. .not.DO_NewGCMGLINT ) then
        PAR_STR = 'Do surface emission?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_SURFACE_EMISSION
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      ENDIF

!  Only if set

      IF ( DO_BRDF_SURFACE ) THEN

!  Get the NewCM parameters (either, updated 7/4/15)
!  ------------------------

        IF ( DO_NewCMGLINT .or. DO_NewGCMGLINT ) then

!  Flags

          PAR_STR = 'Do NewCM glint shadowing?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_GlintShadow
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM whitecap (foam) reflectance?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FoamOption
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'Do NewCM facet isotropy?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) DO_FacetIsotropy
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Wavelength, Salinity

          PAR_STR = 'NewCM Wavelength [Microns]?'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) WAVELENGTH
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

          PAR_STR = 'NewCM Ocean water salinity [ppt]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) SALINITY
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Windspeed

          PAR_STR = 'NewCM Windspeed in [m/s]'
          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
             READ (FILUNIT,*,ERR=998) WINDSPEED
          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Wind-directions, only required for non DO_FacetIsotropy
!   Multiple SZAS not allowed with Facet ANISOTROPY. Upgrade 1/5/16.
!      LOGIC changed 1/5/16 to read wind directions.....

          IF ( .not. Do_FacetIsotropy ) then
            if ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            else
              PAR_STR = 'NewCM Wind directions (degrees) relative to sun positions'
              IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
                DO I = 1, NBEAMS
                  READ (FILUNIT,*,ERR=998) WINDDIR(I)
                ENDDO
              ENDIF
              CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
            endif
          ENDIF

        ENDIF

!  NewCM input: Set kernel defaults
!     One kernel, no surface emission, no scaling
!    7/4/15, relax Scalar-only condition

!    2/25/22. Version 2.8.5, Only 1 kernel, so don't need vector factors (already initialized)

        IF ( DO_NewCMGLINT ) then
           DO_NewGCMGLINT         = .false.
           N_BRDF_KERNELS         = 1
           BRDF_NAMES(1)          = 'NewCMGlint'
           WHICH_BRDF(1)          = NewCMGLINT_IDX
           BRDF_FACTORS(1)        = one
           DO_SEPARATE_FACTORS(1) = .false.
           BRDF_VECTOR_FACTORS(1) = ZERO
           N_BRDF_PARAMETERS(1)   = 2
           BRDF_PARAMETERS(1,1)   = WINDSPEED
           BRDF_PARAMETERS(1,2)   = SALINITY
        endif

!  7/4/15, Introduce Vector condition

        IF ( DO_NewGCMGLINT ) then
           DO_NewCMGLINT          = .false.
           N_BRDF_KERNELS         = 1
           BRDF_NAMES(1)          = 'NewGCMGlit'
           WHICH_BRDF(1)          = NewGCMGLINT_IDX
           BRDF_FACTORS(1)        = one
           DO_SEPARATE_FACTORS(1) = .false.
           BRDF_VECTOR_FACTORS(1) = ZERO
           N_BRDF_PARAMETERS(1)   = 2
           BRDF_PARAMETERS(1,1)   = WINDSPEED
           BRDF_PARAMETERS(1,2)   = SALINITY
        endif

!  Skip next section if doing the one of the NewCM kernels

        IF ( DO_NewCMGLINT .or. DO_NewGCMGLINT  ) go to 656

!  Basic KERNEL BRDF inputs
!  ------------------------

!  Number of kernels

        PAR_STR = 'Number of BRDF kernels'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
           READ (FILUNIT,*,ERR=998) N_BRDF_KERNELS
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Check Dimension. Rob Fix 3/17/15, No longer 3

        IF ( N_BRDF_KERNELS .GT. MAX_BRDF_KERNELS ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of BRDF Kernels > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAX_BRDF_KERNELS dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  Main kernel input
!   Rob Fix 9/25/14. Version 2.7, F8.4 format for the BRDF_FACTORS
!   Rob Fix 2/22/16, Version 2.8, Now 4 parameters per kernel, change formatting!

!  2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

        PAR_STR = 'Kernel names, indices, amplitude control and factors, # parameters, parameters'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
          DO I = 1, N_BRDF_KERNELS
            READ (FILUNIT,56,ERR=998) &
                BRDF_NAMES(I), WHICH_BRDF(I), DO_SEPARATE_FACTORS(I), BRDF_FACTORS(I), BRDF_VECTOR_FACTORS(I), &
               N_BRDF_PARAMETERS(I),(BRDF_PARAMETERS(I,K),K=1,4)
          ENDDO
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

! 56     FORMAT( A10, I2, F8.4, I2, 3F12.6 ) ! Version 2.7
! 56     FORMAT( A10, I2, F8.4, I2, 4F12.6 ) ! Version 2.8.3
 56     FORMAT( A10, I3, L2, F8.4, F8.4, I2, 4F12.6)

!  Check Kernel indices are within bounds. Check BRDF name is on accepted list

        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).GT.MAXBRDF_IDX.OR.WHICH_BRDF(K).LE.0) THEN
            NM = NM + 1
            MESSAGES(NM) = 'Bad input: BRDF Index not on list of indices'
            ACTIONS(NM)  = 'Re-set input value: Look in VLIDORT_PARS for correct index'
            STATUS = VLIDORT_SERIOUS
            NMESSAGES = NM
            GO TO 764
          ELSE
            IF ( BRDF_NAMES(K).NE.BRDF_CHECK_NAMES(WHICH_BRDF(K)) ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: BRDF kernel name not one of accepted list'
              ACTIONS(NM)  = 'Re-set input value: Look in VLIDORT_PARS for correct name'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
            ENDIF
          ENDIF
        ENDDO

!  **************************************************************
!  Rob Fix 9/25/14. New section on use of MODIS-type kernels
!   2/22/16. Version 2.8 Also includes Ross-Thick Hotspot Alternative (RTKHOTSPOT)
!   2/28/22. Version 2.8.5. Add Checks on use of Kernel 20 (do_RTS)

        DO_ROSS = .false. ; do_Li = .false. ; do_Lamb = .false. ; do_CM = .false. ; do_RTS = .false.
        DO K = 1, N_BRDF_KERNELS
          IF ( WHICH_BRDF(K).eq.ROSSTHIN_IDX.or.WHICH_BRDF(K).eq.ROSSTHICK_IDX &
                                            .or.WHICH_BRDF(K).eq.RTKHOTSPOT_IDX ) DO_ROSS = .true.
          IF ( WHICH_BRDF(K).eq.LISPARSE_IDX.or.WHICH_BRDF(K).eq.LIDENSE_IDX    ) DO_LI   = .true.
          IF ( WHICH_BRDF(K).eq.LAMBERTIAN_IDX ) DO_Lamb = .true.
          IF ( WHICH_BRDF(K).eq.COXMUNK_IDX )    DO_CM   = .true.
          IF ( WHICH_BRDF(K).eq.SNOWRTS_IDX )    DO_RTS  = .true.
        ENDDO

        IF ( ( DO_ROSS .and. .not. DO_Lamb ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross-type and no Lambertian type kernel NOT ALLOWED'
          ACTIONS(NM)  = 'Re-set input: Ross/Lambertian kernels must be together for MODIS-based BRDFs'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

        IF ( DO_RTS .and. .not. DO_ROSS ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: SnowRTS-type and no Ross-type kernel NOT ALLOWED'
          ACTIONS(NM)  = 'Re-set input: Ross/RTSSnow kernels must be together for MODIS-based BRDFs'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

        IF ( .not.DO_ROSS .and. DO_LI ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Li-type and no Ross-type kernel NOT ALLOWED'
          ACTIONS(NM)  = 'Re-set input: Ross/Li kernels must be together for MODIS-based BRDFs'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

        IF ( ( DO_ROSS .and. DO_LI ) .and. (.not. DO_Lamb .and. .not. do_CM ) ) then
          NM = NM + 1
          MESSAGES(NM) = 'Bad input: Ross/Li kernels must go with Lambertian or Cox-Munk'
          ACTIONS(NM)  = 'Re-set input kernels to include Lambertian or Cox-Munk kernel'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  END section. Rob Fix 9/25/14. Check use of MODIS-type kernels
!  **************************************************************

!  2/25/22. Version 2.8.5. 3/28/22. Expanded error handling to cover all cases.
!    -- Scalar kernels (1-9,15,17,19,20) have no separate vector amplitudes. 
!    -- Vector kernels (10-14,16,18) MAY have separate vector amplitudes - check validity

        DO K = 1, N_BRDF_KERNELS
           IF ( WHICH_BRDF(K).le.COXMUNK_IDX    .or. WHICH_BRDF(K).eq.RTKHOTSPOT_IDX .or. &
                WHICH_BRDF(K).eq.NewCMGLINT_IDX .or. WHICH_BRDF(K).eq.SNOWBRDF_IDX   .or. &
                WHICH_BRDF(K).eq.SNOWRTS_IDX ) then
              IF ( DO_SEPARATE_FACTORS(k) ) then
                NM = NM + 1
                MESSAGES(NM) = BRDF_NAMES(K)//' is a SCALAR kernel: DO_SEPARATE_FACTORS flag should be OFF!'
                ACTIONS(NM)  = 'Turn off input kernel flag DO_SEPARATE_FACTORS'
                STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ;   GO TO 764
              ELSE if ( BRDF_VECTOR_FACTORS(K) .ne. ZERO ) THEN
                NM = NM + 1
                MESSAGES(NM) = BRDF_NAMES(K)//' is a SCALAR kernel: BRDF_VECTOR_FACTORS should be 0.0!'
                ACTIONS(NM)  = 'Set input kernel factor BRDF_VECTOR_FACTORS to Zero'
                STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ;   GO TO 764
              ENDIF
           ELSE 
              IF ( DO_SEPARATE_FACTORS(k) ) then
                 if ( BRDF_VECTOR_FACTORS(K) .eq. ZERO ) THEN
                    NM = NM + 1
                    MESSAGES(NM) = BRDF_NAMES(K)//' is VECTOR kernel with Separate factor ==> BRDF_VECTOR_FACTORS should NOT = 0.0!'
                    ACTIONS(NM)  = 'Set input kernel factor BRDF_VECTOR_FACTORS to a positive non-zero value'    
                    STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ;   GO TO 764
                 else if ( BRDF_FACTORS(K) .ne. ZERO ) THEN
                    NM = NM + 1
                    MESSAGES(NM) = BRDF_NAMES(K)//' is VECTOR kernel with Separate factor ==> BRDF_FACTORS should = 0.0!'
                    ACTIONS(NM)  = 'Set input kernel factor BRDF_FACTORS to zero'    
                    STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ;   GO TO 764
                 ENDIF
              ELSE
                 if ( BRDF_VECTOR_FACTORS(K) .eq. ZERO ) THEN
                    NM = NM + 1
                    MESSAGES(NM) = BRDF_NAMES(K)//' is VECTOR kernel with NO Separate factor ==> BRDF_VECTOR_FACTORS should = 0.0!'
                    ACTIONS(NM)  = 'Set input kernel factor BRDF_VECTOR_FACTORS to zero'    
                    STATUS = VLIDORT_SERIOUS ; NMESSAGES = NM ;   GO TO 764
                 ENDIF
              ENDIF 
           ENDIF
        ENDDO

!  Set the Lambertian kernel flags

        DO I = 1, N_BRDF_KERNELS
          IF ( BRDF_NAMES(I) .EQ. 'Lambertian' ) THEN
            LAMBERTIAN_KERNEL_FLAG(I) = .true.
          ENDIF
        ENDDO

!  Shadowing input (for Cox-Munk types)

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do shadow effect for glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_SHADOW_EFFECT
           ENDIF
          ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  Continuation point for skipping Regular kernel input

656     continue

!  General inputs, all types
!  -------------------------

!  Number of BRDF azimuth streams, check this value

        PAR_STR = 'Number of BRDF azimuth angles'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) &
            READ (FILUNIT,*,ERR=998) NSTREAMS_BRDF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

        IF ( NSTREAMS_BRDF .GT. MAXSTREAMS_BRDF ) THEN
          NM = NM + 1
          MESSAGES(NM) = 'Number of  BRDF streams > maximum dimension'
          ACTIONS(NM)  = 'Re-set input value or increase MAXSTREAMS_BRDF dimension'
          STATUS = VLIDORT_SERIOUS
          NMESSAGES = NM
          GO TO 764
        ENDIF

!  **********************************************************************
!  Rob Fix 9/25/14. New Variable and more explicit character string

! FORMER CODE
!  !@@ Overall-Exact flag
!        PAR_STR = 'Do Overall-Exact kernels?'
!        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!           READ (FILUNIT,*,ERR=998)DO_EXACT
!        ENDIF
!        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!  Exact only flag. Only if above is set (!@@)
!        IF ( DO_EXACT ) THEN
!          PAR_STR = 'Do Exact-only (no Fourier) kernels?'
!          IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
!             READ (FILUNIT,*,ERR=998)DO_EXACTONLY
!          ENDIF
!          CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
!        ENDIF
! END FORMER CODE

!  New flag for turning on the Direct-Bounce only flag ("exact" BRDF)
!   Normally this would be .FALSE., if set, then no multiple-scatter BRDFs will be done

        PAR_STR = 'Do direct-bounce only (no multiple-scatter contributions to BRDF)?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_DBONLY
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents
!     -- default should be to use direct-bounce for this term

        IF ( DO_USER_STREAMS ) then
           PAR_STR = 'Do include truncated Fourier solution for direct-bounce BRDF?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)DO_FOURIER_USERSUN
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
        ENDIF

!   END Replacement Section Rob Fix 9/25/14
! ********************************************************************************

!  Multiple reflectance correction (for Cox-Munk types)
!  ----------------------------------------------------

!  General flag

        DO I = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(I) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(I) .EQ. 'GissCoxMnk' .OR. &
              BRDF_NAMES(I) .EQ. 'GCMcomplex' ) THEN
           PAR_STR = 'Do multiple reflectance for all glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR
           ENDIF
         ENDIF
        ENDDO
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Specific MSRCORR inputs
!   Rob Fix 9/25/14. Variable name changed, the wording made more explicit

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Do multiple reflectance for just the direct-bounce glitter kernels?'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
            READ (FILUNIT,*,ERR=998)DO_MSRCORR_DBONLY
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR scattering order

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering order for glitter kernels'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_ORDER
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  MSRCORR quadrature orders

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Polar quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NMUQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

        IF ( DO_MSRCORR ) THEN
           PAR_STR = 'Multiple reflectance scattering; Azimuth quadrature order'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR )) THEN
            READ (FILUNIT,*,ERR=998)MSRCORR_NPHIQUAD
           ENDIF
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Check MSCORR dimensions

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_NMUQUAD .gt. max_msrs_muquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR polar quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_muquad in vlidort_pars'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
           IF ( MSRCORR_NPHIQUAD .gt. max_msrs_phiquad ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR azimuth quadrature No. > Dimensioning'
              ACTIONS(NM)  = 'Increase value of max_msrs_phiquad in vlidort_pars'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check on MSRCORR order

        IF ( DO_MSRCORR ) THEN
           IF ( MSRCORR_ORDER .EQ.0 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: MSR is on, but scattering order = 0'
              ACTIONS(NM)  = 'Turn off MSRCORR flags and proceed with warning'
              DO_MSRCORR = .false. ; DO_MSRCORR_DBONLY = .false.
              STATUS = VLIDORT_WARNING
              NMESSAGES = NM
           ENDIF
        ENDIF

!  White-Sky and Black-Sky Albedo scalings. New for Version 2.7
!  ============================================================

!  Skip this section for one of the NewCM kernels

        if ( DO_NewCMGLINT .or. DO_NewGCMGLINT ) go to 646

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  Output option flag.
!  -------------------

!   Rob Fix 9/27/14; added output option flag.
!   If you are doing WSA or BSA scaling, this should be set automatically

        PAR_STR = 'Do white-sky and black-sky albedo output?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSABSA_OUTPUT
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  White-Sky inputs
!  ----------------

!  White-sky Albedo scaling

        PAR_STR = 'Do white-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_WSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  WSA value. This could be extracted from a data set.....

        IF ( DO_WSA_SCALING  ) THEN
           PAR_STR = 'White-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)WSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check WSA value

        IF ( DO_WSA_SCALING  ) THEN
           IF ( WSA_VALUE .le.zero .or. WSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: White-sky albedo value not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Black-Sky inputs
!  ----------------

!  Black-sky Albedo scaling.

        PAR_STR = 'Do black-sky albedo scaling?'
        IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
           READ (FILUNIT,*,ERR=998)DO_BSA_SCALING
        ENDIF
        CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)

!  Cannot have BSA and WSA together

        IF ( DO_BSA_SCALING .and. DO_WSA_SCALING ) THEN
           NM = NM + 1
           MESSAGES(NM) = 'Bad input: Cannot apply both Black-sky albedo and White-sky albedo scalings!'
           ACTIONS(NM)  = 'Make a choice of which one you want! '
           STATUS = VLIDORT_SERIOUS
           NMESSAGES = NM
           GOTO 764
        ENDIF

!  BSA value. This could be extracted from a data set.....
!    WARNING: ONLY ALLOWED ONE VALUE HERE...................

        IF ( DO_BSA_SCALING  ) THEN
           PAR_STR = 'Black-sky albedo value'
           IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) THEN
              READ (FILUNIT,*,ERR=998)BSA_VALUE
           ENDIF
           CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS)
        ENDIF

!  Check BSA value

        IF ( DO_BSA_SCALING  ) THEN
           IF ( BSA_VALUE .le.zero .or. BSA_VALUE .gt. one ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Black-sky albedo value is not in the range [0,1]'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check that Solar source flag is on for Black sky albedo

        IF ( DO_BSA_SCALING  ) THEN
           IF ( .not. DO_SOLAR_SOURCES ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Cannot have Black-sky albedo if Solar_sources not turned on'
              ACTIONS(NM)  = 'Fix the input (turn on DO_SOLAR_SOURCES)'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Check that there is only one beam for Black sky albedo
!   Revision 12 august 2014. Same check applies to the output option

        IF ( DO_BSA_SCALING  ) THEN
           IF ( NBEAMS.gt.1 ) THEN
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: Cannot have Black-sky albedo with more than 1 solar angle'
              ACTIONS(NM)  = 'Fix the input (set NBEAMS = 1)'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF
        ENDIF

!  Continuation point for avoiding the WSA/BSA albedos

646     continue

!  Checking NewCM Kernels. New for Version 2.7. 7/4/15 upgrade
!    Single kernel, Solar Sources only, No Surface Emission. 
!      No MSR (Multiple-surface reflections), No Scaling

        if ( DO_NewCMGLINT .or. DO_NewGCMGLINT ) then

           IF ( WINDSPEED .le.zero ) then
              NM = NM + 1
              MESSAGES(NM) = 'Bad input: NewCM windspeed value is negative'
              ACTIONS(NM)  = 'Fix the input'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF

!  Rob Fix 9/27/14
!   Added check here on number of Solar beams ( = 1 for FacetIsotropy = .false. )

           IF ( .not. DO_FacetIsotropy .and. NBEAMS .gt. 1 ) then
              NM = NM + 1
              MESSAGES(NM) = 'Facet Anisotropy (wind direction) not Allowed for NBEAMS > 1'
              ACTIONS(NM)  = 'Either re-set NBEAMS = 1 or re-set Do_FacetIsotropy = .true.'
              STATUS = VLIDORT_SERIOUS
              NMESSAGES = NM
              GO TO 764
           ENDIF

        endif

!  End BRDF surface clause

      ENDIF

!  Successful finish

      CLOSE(FILUNIT)

!mick fix
      NMESSAGES = NM

!  Copy Control inputs
!    -- 1/31/21. Version 2.8.3. Add doublet geometry option

      VBRDF_Sup_In%BS_DO_BRDF_SURFACE     = DO_BRDF_SURFACE
      VBRDF_Sup_In%BS_DO_USER_STREAMS     = DO_USER_STREAMS
      VBRDF_Sup_In%BS_DO_SOLAR_SOURCES    = DO_SOLAR_SOURCES    !@@
      VBRDF_Sup_In%BS_DO_USER_OBSGEOMS    = DO_USER_OBSGEOMS    !@@
      VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY = DO_DOUBLET_GEOMETRY !@@
      VBRDF_Sup_In%BS_DO_SURFACE_EMISSION = DO_SURFACE_EMISSION

!  Copy Geometry results
!mick fix 1/5/2021 - added N_USER_DOUBLETS & USER_DOUBLETS

      VBRDF_Sup_In%BS_NSTOKES           = NSTOKES
      VBRDF_Sup_In%BS_NSTREAMS          = NSTREAMS

      VBRDF_Sup_In%BS_NBEAMS            = NBEAMS
      VBRDF_Sup_In%BS_BEAM_SZAS         = BEAM_SZAS
      VBRDF_Sup_In%BS_N_USER_STREAMS    = N_USER_STREAMS
      VBRDF_Sup_In%BS_USER_ANGLES_INPUT = USER_ANGLES
      VBRDF_Sup_In%BS_N_USER_RELAZMS    = N_USER_RELAZMS
      VBRDF_Sup_In%BS_USER_RELAZMS      = USER_RELAZMS

      VBRDF_Sup_In%BS_N_USER_OBSGEOMS   = N_USER_OBSGEOMS !@@
      VBRDF_Sup_In%BS_USER_OBSGEOMS     = USER_OBSGEOMS   !@@
      VBRDF_Sup_In%BS_N_USER_DOUBLETS   = N_USER_DOUBLETS
      VBRDF_Sup_In%BS_USER_DOUBLETS     = USER_DOUBLETS

!  Copy BRDF inputs

      VBRDF_Sup_In%BS_NSTREAMS_BRDF          = NSTREAMS_BRDF
      VBRDF_Sup_In%BS_N_BRDF_KERNELS         = N_BRDF_KERNELS

!    2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

      VBRDF_Sup_In%BS_WHICH_BRDF             = WHICH_BRDF
      VBRDF_Sup_In%BS_BRDF_NAMES             = BRDF_NAMES
      VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = LAMBERTIAN_KERNEL_FLAG

      VBRDF_Sup_In%BS_DO_SEPARATE_FACTORS    = DO_SEPARATE_FACTORS
      VBRDF_Sup_In%BS_BRDF_FACTORS           = BRDF_FACTORS
      VBRDF_Sup_In%BS_BRDF_VECTOR_FACTORS    = BRDF_VECTOR_FACTORS
      VBRDF_Sup_In%BS_N_BRDF_PARAMETERS      = N_BRDF_PARAMETERS
      VBRDF_Sup_In%BS_BRDF_PARAMETERS        = BRDF_PARAMETERS

!  Rob Fix 9/25/14. Two variables replaced
!      VBRDF_Sup_In%BS_DO_EXACT               = DO_EXACT         !@@
!      VBRDF_Sup_In%BS_DO_EXACTONLY           = DO_EXACTONLY
      VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY   = DO_DBONLY

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents

      VBRDF_Sup_In%BS_DO_FOURIER_USERSUN     = DO_FOURIER_USERSUN

!  Shadow effect flag (only for old style Cox-Munk type kernels)

      VBRDF_Sup_In%BS_DO_SHADOW_EFFECT       = DO_SHADOW_EFFECT

!  Multiple-scattering Glitter options

      VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR           = DO_MSRCORR
      VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY    = DO_MSRCORR_DBONLY  !  Rob Fix 9/25/14, name 
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER        = MSRCORR_ORDER
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD      = MSRCORR_NMUQUAD
      VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD     = MSRCORR_NPHIQUAD

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!  Rob Fix 9/27/14. Output variable added

      VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT = DO_WSABSA_OUTPUT
      VBRDF_Sup_In%BS_DO_WSA_SCALING   = DO_WSA_SCALING
      VBRDF_Sup_In%BS_DO_BSA_SCALING   = DO_BSA_SCALING
      VBRDF_Sup_In%BS_WSA_VALUE        = WSA_VALUE
      VBRDF_Sup_In%BS_BSA_VALUE        = BSA_VALUE

!  NewCM options. Add  vector term, 7/4/15

      VBRDF_Sup_In%BS_DO_NewCMGLINT    = DO_NewCMGLINT
      VBRDF_Sup_In%BS_DO_NewGCMGLINT   = DO_NewGCMGLINT

      VBRDF_Sup_In%BS_SALINITY         = SALINITY
      VBRDF_Sup_In%BS_WAVELENGTH       = WAVELENGTH

      VBRDF_Sup_In%BS_WINDSPEED        = WINDSPEED
      VBRDF_Sup_In%BS_WINDDIR          = WINDDIR

      VBRDF_Sup_In%BS_DO_GlintShadow   = DO_GlintShadow
      VBRDF_Sup_In%BS_DO_FoamOption    = DO_FoamOption
      VBRDF_Sup_In%BS_DO_FacetIsotropy = DO_FacetIsotropy

!  Exception handling

      VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Normal return

      RETURN

!  Open file error

300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//trim(adjustl(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
      GO TO 764

!  Line read error - abort immediately

998   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'read failure for '//trim(adjustl(PAR_STR))
      ACTIONS(NMESSAGES)  = 'Re-set: Entry is incorrect in input file'
      CLOSE(FILUNIT)

!  Final error copying

764   CONTINUE

      VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD = STATUS
      VBRDF_Sup_InputStatus%BS_NINPUTMESSAGES   = NMESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTMESSAGES    = MESSAGES
      VBRDF_Sup_InputStatus%BS_INPUTACTIONS     = ACTIONS

!  Finish

      RETURN
      END SUBROUTINE VBRDF_INPUTMASTER

!

      SUBROUTINE VBRDF_MAINMASTER ( &
        DO_DEBUG_RESTORATION,     & ! Inputs
        NMOMENTS_INPUT,           & ! Inputs
        VBRDF_Sup_In, Core, Book, & ! Inputs
        VBRDF_Sup_Out,            & ! Outputs
        VBRDF_Sup_OutputStatus )   ! Output Status

!  Prepares the bidirectional reflectance functions necessary for VLIDORT.

!  Version 2.6 notes
!  -----------------

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012. 
!       Observation-Geometry input control.       DO_USER_OBSGEOMS
!       Observation-Geometry input control.       N_USER_OBSGEOMS
!       User-defined Observation Geometry angles. USER_OBSGEOMS
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  VBRDF Upgrades for Version 2.7
!  ------------------------------

!  A. White-sky and Black-sky scaling options
!  ==========================================

!  WSA and BSA scaling options.
!   first introduced 02 April 2014, Revised, 14-15 April 2014
!      WSA = White-sky albedo. BSA = Black-sky albedo.

!  These options are mutually exclusive. If either is set, the VBRDF code
!  will automatically perform an albedo calculation (either WS or BS) for
!  the (1,1) component of the complete 3-kernel BRDF, then normalize the
!  entire BRDF with this albedo before scaling up with the externally chosen
!  WSA or BSA (from input).

!  Additional exception handling has been introduced to make sure that
!  the spherical or planar albedos for a complete 3-kernel BRDF are in
!  the [0,1] range. Otherwise, the WSA/BSA scaling makes no sense. It is
!  still the case that some of the MODIS-tpe kernels give NEGATIVE albedos.

!  The albedo scaling process has been linearized for all existing surface
!  property Jacobians available to the VBRDF linearized supplement. In
!  addition, it is also possible to derive single Jacobians of the BRDFs
!  with respect to the WS or BS albedo - this is separate from (and orthogonal
!  to) the usual kernel derivatives.

!  B. Alternative Cox-Munk Glint Reflectance
!  =========================================

!  In conjunction with new Water-leaving code developed for the VSLEAVE
!  supplement, we have given the VBRDF supplement a new option to 
!  return the (scalar) Cox-Munk glint reflectance, based on code originally
!  written for the 6S code.

!  Developed and tested by R. Spurr, 21-29  April 2014
!  Based in part on Modified-6S code by A. Sayer (NASA-GSFC).
!  Validated against Modified-6S OCEABRDF.F code, 24-28 April 2014.

!  The new glint option depends on Windspeed/direction, with refractive
!  indices now computed using salinity and wavelength. There is now an 
!  optional correction for (Foam) Whitecaps (Foam). These choices come from
!  the 6S formulation.

!  Need to make sure that the wind input information is the same as that
!  use for glint calculations in the VSLEAVE supplement when the Glitter
!  kernels are in use. Also, the Foam correction applied here in the
!  surface-leaving code should also be applied in the VSLEAVE system..

!  Choosing this new glint option bypasses the normal kernel inputs (and
!  also the WSA/BSA inputs). Instead, a single-kernel glint reflectance
!  is calculated with amplitude 1.0, based on a separate set of dedicated
!  inputs. This option does not apply for surface emission - the solar
!  sources flag must be turned on. There is only 1 Jacobian available - 
!  with respect to the windspeed. 

!  Note that the use of the facet isotropy flag is recommended for
!  multi-beam runs. This is because the wind-direction is a function of
!  the solar angle, and including this wind-direction in the glint
!  calculations for VBRDF will only work if there is just one SZA. This
!  condition is checked for internally.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!  Upgrades for Version 2.8
!  ------------------------

!  Presence of 2 new kernels, 

!  RTK_HOTSPOT : This is the Old Ross-Thick kernel with a hot-spot modification
!     Derives from the following references.
!  F. M. Breon, F. Maignan, M. Leroy and I. Grant, 
!    "Analysis of hot spot directional siganatures measured from space",
!      J. Geophys. Res., 107, D16, 4282, (2002)
!  E. Vermote C. Justice, and F. M. Breon,
!    "Towards a generalized approach for correction of the BRDF effect in MODIS reflectances"
!      IEEE Trans. Geo. Rem. Sens., 10.1109/TGRS.2008.2005997 (2008)
!   -- Ross-Thick Hotspot Kernel from

!  MODFRESNEL. This is a "Modified Fresnel" kernel developed for polarized reflectances
!   Taken from the following reference.
!  P. Litvinov, O. Hasekamp and B. Cairns,
!    "Models for surface reflection of radiance and polarized radiance: Comparison
!     with airborne multi-angle photopolarimetric measurements and implications for
!     modeling top-of-atmopshere measurements,
!       Rem. Sens. Env., 115, 781-792, (2011).

!  4 parameters now allowed in Kernels

!  Patch Overhaul for Version 2.8.1. RobFix 11/8/19
!    - Fourier routine has one less argument
!    - Introduce reflectivity Mask for proper identification of Matrix elements

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel. First introduced to VLIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  1/31/21. Version 2.8.3. Add doublet geometry option

!  7/28/21. Version 2.8.3. Some Changes
!    -- Half-Range azimuth integration introduced, now the default option. Controlled by parameter setting
!    -- Adjustment of +/- signs in Fourier-component sine-series settings. Validates against TestBed
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/BAX/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Scalar and vector amplitudes given separately, controlled with flag
!    -- remove calculation of USER_BRDF_F_0. Flag for this.
!    -- reduced kernel for SNOW BRDF

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    ** VBRDF_MAINMASTER has major re-engineering
!        -- now only does kernel-scaling and Fourier components
!        -- Assumes Kernel cores already computed and stored in "core" type structure
!        -- Also uses a bookkeeping type structure ("book") with saved geometrical quantities

! #####################################################################
! #####################################################################

      USE vlidort_pars_m

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    ** VBRDF_MAINMASTER has major re-engineering
!        -- now only does kernel-scaling and Fourier components
!        -- Module vbrdf_sup_core_def_m has stored "core" and  "book" type structures
!        -- Module vbrdf_sup_kernels_experimental_m   has the scaling routines (that depend on wavelength)
!        -- Module vbrdf_sup_makerfourier_routines_m  applies scaling to kernel cores, then does Fourier terms

      USE vbrdf_sup_inputs_def_m
      USE vbrdf_sup_core_def_m
      USE vbrdf_sup_outputs_def_m
      USE vbrdf_sup_kernels_experimental_m, only : vbrdf_water_refracindex, vbrdf_whitecap_reflectance
      USE vbrdf_sup_makerfourier_routines_m

!  Implicit none

      IMPLICIT NONE

!  Inputs
!  ------

!  Debug flag for restoration

      LOGICAL, INTENT(IN) ::              DO_DEBUG_RESTORATION

!  Input number of moments (only used for restoration debug)

      INTEGER, INTENT(IN) ::              NMOMENTS_INPUT

!  Input structure
!  ---------------

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    ** New "Core" and "Book" type structures

      TYPE(VBRDF_Sup_Inputs), INTENT(IN)   :: VBRDF_Sup_In
      TYPE(VBRDF_Sup_Core)  , INTENT(IN)   :: Core
      TYPE(VBRDF_Sup_Book)  , INTENT(IN)   :: Book

!  Output structure
!  ----------------

      TYPE(VBRDF_Sup_Outputs), INTENT(OUT) :: VBRDF_Sup_Out

!  Exception handling introduced 02 April 2014 for Version 2.7

      TYPE(VBRDF_Output_Exception_Handling), INTENT(OUT) :: VBRDF_Sup_OutputStatus

!  VBRDF local and proxy variables
!  +++++++++++++++++++++++++++++++

!  User stream Control

      LOGICAL ::          DO_USER_STREAMS

!  Surface emission

      LOGICAL ::          DO_SURFACE_EMISSION

!  number of Stokes components

      INTEGER ::          NSTOKES, NSTOKESSQ

!   Number and index-list of bidirectional functions

      INTEGER ::          N_BRDF_KERNELS
      INTEGER ::          WHICH_BRDF ( MAX_BRDF_KERNELS )

!  Parameters required for Kernel families

      INTEGER ::          N_BRDF_PARAMETERS ( MAX_BRDF_KERNELS )
      DOUBLE PRECISION :: BRDF_PARAMETERS   ( MAX_BRDF_KERNELS, MAX_BRDF_PARAMETERS )

!  BRDF names

      CHARACTER (LEN=10) :: BRDF_NAMES ( MAX_BRDF_KERNELS )

!  Lambertian Surface control

      LOGICAL ::          LAMBERTIAN_KERNEL_FLAG ( MAX_BRDF_KERNELS )

!  Input kernel amplitude factors
!    2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

      LOGICAL          :: DO_SEPARATE_FACTORS ( MAX_BRDF_KERNELS ) 
      DOUBLE PRECISION :: BRDF_FACTORS        ( MAX_BRDF_KERNELS ) 
      DOUBLE PRECISION :: BRDF_VECTOR_FACTORS ( MAX_BRDF_KERNELS )

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!  Rob Fix 9/27/14. Output variable added

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING
      REAL(fpk) :: WSA_VALUE, BSA_VALUE

!  Shadowing effect flag (only for Cox-Munk type kernels)

      LOGICAL ::          DO_SHADOW_EFFECT

!  Solar sources + Observational Geometry flag 
!    -- 1/31/21. Version 2.8.3. Add doublet geometry flag

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!   Flag for the Direct-bounce term, replaces former "EXACT" variables 

      LOGICAL ::          DO_DBONLY

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents
!     -- default should be to use direct-bounce for this term

      LOGICAL ::          DO_FOURIER_USERSUN

!  Local angle control

      INTEGER ::          N_USER_OBSGEOMS
      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  7/28/21. 8/18/21. Half-range azimuth integration variable introduced. This is the default.
!   -- if in operation, azimuthal range is [0,pi], if not the range is [pi,pi]
!   -- first set as a hard-wired parameter, now determined according to type of kernel
!         (TRUE for all kernels with facet-isotropy, False for NewCM/NewGCM with facet anisotropy

!  BRDF azimuth quadrature streams

      LOGICAL ::          DO_HALF_RANGE
      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

!  Discrete ordinates (local, for Albedo scaling). Version 2.7.

      INTEGER          :: SCALING_NSTREAMS

!  White-sky and Black-sky albedos. Version 2.7.

      LOGICAL          :: DO_LOCAL_WSA, DO_LOCAL_BSA

!  New Cox-Munk Glint reflectance options (bypasses the usual Kernel system)
!  -------------------------------------------------------------------------

!  Overall flags for this option. 7/4/15 upgrade

      LOGICAL   :: DO_NewCMGLINT, DO_NewGCMGLINT

!  Input Salinity in [ppt]

      REAL(fpk) :: SALINITY

!  Input wavelength in [Microns]

      REAL(fpk) :: WAVELENGTH

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Local BRDF functions
!  ====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS, MAXSTREAMS_BRDF )

!  DB Kernel values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity
!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC      ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF)
      DOUBLE PRECISION :: USER_EBRDFUNC ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF)

!  Values for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0 ( MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  Local kernel Fourier components
!  ===============================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: LOCAL_BRDF_F   ( MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_BRDF_F_0 ( MAXSTOKES_SQ, MAXSTREAMS, MAXBEAMS   )

!  at user-defined stream directions

      DOUBLE PRECISION :: LOCAL_USER_BRDF_F   ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_USER_BRDF_F_0 ( MAXSTOKES_SQ, MAX_USER_STREAMS, MAXBEAMS   )

!  emissivities

      DOUBLE PRECISION :: LOCAL_EMISSIVITY      ( MAXSTOKES, MAXSTREAMS )
      DOUBLE PRECISION :: LOCAL_USER_EMISSIVITY ( MAXSTOKES, MAX_USER_STREAMS )

!  WSA/BSA scaling componnets, at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: SCALING_BRDF_F   ( MAXSTREAMS_SCALING, MAXSTREAMS_SCALING )
      DOUBLE PRECISION :: SCALING_BRDF_F_0 ( MAXSTREAMS_SCALING   )

!  Exception handling
!  ==================

!   New code, 02 April 2014. Version 2.7
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Other help variables
!  ====================

      DOUBLE PRECISION :: WSA_CALC (MAX_BRDF_KERNELS), TOTAL_WSA_CALC
      DOUBLE PRECISION :: BSA_CALC (MAX_BRDF_KERNELS), TOTAL_BSA_CALC

!  Local NewCM variables

      DOUBLE PRECISION :: Refrac_R, Refrac_I, WC_Reflectance, WC_Lambertian, Wc_Correction

!  other help

      INTEGER ::          K, B, I, I1, I2, J, IB, UI, UM, IA, M, O1, Q
      INTEGER ::          BRDF_NPARS, NMOMENTS
      DOUBLE PRECISION :: PARS ( MAX_BRDF_PARAMETERS ), BRDF_FACTOR_LOCAL(16)
      DOUBLE PRECISION :: HELP_A, SUM, SCALING, FF, REFRAC_INDEX
      DOUBLE PRECISION :: WAV, WAVMM, ALPHA, GAMMA, CHI, F1, F2, HELP
      LOGICAL ::          TRAWL
      INTEGER, PARAMETER :: LUM = 1   !@@
      INTEGER, PARAMETER :: LUA = 1   !@@

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!  Snow Model Ice Refractive index data. Formerly in the Kernel. Moved to here.
!   -- Initial data set: Table 1 from K&B, 9 values
!   -- mick mod 1/5/2021 - converted data and pars to "fpk"

      INTEGER, parameter :: nchidat = 9
      REAL(fpk) :: chilams(nchidat), Chivals(nchidat), Logchivals(nchidat)
      data chilams /  0.412_fpk,  0.443_fpk,  0.490_fpk,  0.565_fpk,  0.670_fpk,  0.865_fpk,    1.02_fpk,   1.24_fpk,   1.61_fpk /
      data chivals / 2.5e-9_fpk, 1.7e-9_fpk, 1.8e-9_fpk, 3.0e-9_fpk, 1.8e-8_fpk, 2.5e-7_fpk, 2.25e-6_fpk, 1.2e-5_fpk, 3.3e-4_fpk /

!  Initialize Exception handling
!  -----------------------------

      STATUS = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of VLIDORT BRDF Sup Master'

!  Copy from input type structure
!  ------------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = VBRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = VBRDF_Sup_In%BS_DO_BRDF_SURFACE
      DO_SURFACE_EMISSION = VBRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  Set number of stokes elements and streams

      NSTOKES   = VBRDF_Sup_In%BS_NSTOKES
      NSTOKESSQ = NSTOKES * NSTOKES
      NSTREAMS  = VBRDF_Sup_In%BS_NSTREAMS

!  Copy Geometry flags

      DO_SOLAR_SOURCES    = VBRDF_Sup_In%BS_DO_SOLAR_SOURCES
      DO_USER_OBSGEOMS    = VBRDF_Sup_In%BS_DO_USER_OBSGEOMS
      DO_DOUBLET_GEOMETRY = VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY

!   !@@ Observational/Doublet/Lattice Geometry + Solar sources Optionalities
!  -- Only need the numbers

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = VBRDF_Sup_In%BS_N_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
        ELSE
          NBEAMS         = 1
          N_USER_RELAZMS = 1
          N_USER_STREAMS = N_USER_OBSGEOMS
        ENDIF
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = VBRDF_Sup_In%BS_NBEAMS
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
          N_USER_RELAZMS    = N_USER_STREAMS
        ELSE
!  NOT ALLOWED
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = VBRDF_Sup_In%BS_NBEAMS
          N_USER_RELAZMS    = VBRDF_Sup_In%BS_N_USER_RELAZMS
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
        ELSE
          NBEAMS         = 1 
          N_USER_RELAZMS = 1
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
        ENDIF
      ENDIF

!  Copy BRDF inputs
!    -- 7/28/21. removed NSTREAMS_BRDF setting, now set later on after the half-range flag is set
!    2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

      N_BRDF_KERNELS         = VBRDF_Sup_In%BS_N_BRDF_KERNELS
      BRDF_NAMES             = VBRDF_Sup_In%BS_BRDF_NAMES
      WHICH_BRDF             = VBRDF_Sup_In%BS_WHICH_BRDF

!      N_BRDF_PARAMETERS      = VBRDF_Sup_In%BS_N_BRDF_PARAMETERS
      N_BRDF_PARAMETERS      = Book%N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = VBRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG

      DO_SEPARATE_FACTORS    = VBRDF_Sup_In%BS_DO_SEPARATE_FACTORS
      BRDF_FACTORS           = VBRDF_Sup_In%BS_BRDF_FACTORS
      BRDF_VECTOR_FACTORS    = VBRDF_Sup_In%BS_BRDF_VECTOR_FACTORS

      DO_SHADOW_EFFECT       = VBRDF_Sup_In%BS_DO_SHADOW_EFFECT

!  Scaling NSTREAMS

      SCALING_NSTREAMS       = MAXSTREAMS_SCALING

!  Rob Fix 9/25/14. Replaces DO_EXACT and DO_EXACTONLY

      DO_DBONLY              = VBRDF_Sup_In%BS_DO_DIRECTBOUNCE_ONLY

!  2/25/22. Version 2.8.5. Separate flag for the USER_BRDF_F_0 coomponents

      DO_FOURIER_USERSUN     = VBRDF_Sup_In%BS_DO_FOURIER_USERSUN

!  WSA and BSA scaling options.
!   Revised, 14 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!  Rob Fix 9/27/14. Output variable added

      DO_WSABSA_OUTPUT    = VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT
      DO_WSA_SCALING      = VBRDF_Sup_In%BS_DO_WSA_SCALING
      DO_BSA_SCALING      = VBRDF_Sup_In%BS_DO_BSA_SCALING
      WSA_VALUE           = VBRDF_Sup_In%BS_WSA_VALUE
      BSA_VALUE           = VBRDF_Sup_In%BS_BSA_VALUE

!  NewCM options. 7/4/15 upgrade

      DO_NewCMGLINT  = VBRDF_Sup_In%BS_DO_NewCMGLINT
      DO_NewGCMGLINT = VBRDF_Sup_In%BS_DO_NewGCMGLINT
      SALINITY     = VBRDF_Sup_In%BS_SALINITY
      WAVELENGTH   = VBRDF_Sup_In%BS_WAVELENGTH

      WINDSPEED = VBRDF_Sup_In%BS_WINDSPEED
      WINDDIR   = VBRDF_Sup_In%BS_WINDDIR

      DO_GlintShadow   = VBRDF_Sup_In%BS_DO_GlintShadow
      DO_FoamOption    = VBRDF_Sup_In%BS_DO_FoamOption
      DO_FacetIsotropy = VBRDF_Sup_In%BS_DO_FacetIsotropy 

!  Copy MSR inputs

!      DO_MSRCORR             = VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR
!      DO_MSRCORR_DBONLY      = VBRDF_Sup_In%BS_DO_GLITTER_MSRCORR_DBONLY 
!      MSRCORR_ORDER          = VBRDF_Sup_In%BS_GLITTER_MSRCORR_ORDER
!      N_MUQUAD               = VBRDF_Sup_In%BS_GLITTER_MSRCORR_NMUQUAD
!      N_PHIQUAD              = VBRDF_Sup_In%BS_GLITTER_MSRCORR_NPHIQUAD

!  Number of Fourier components to calculate

      IF ( DO_DEBUG_RESTORATION ) THEN
        NMOMENTS = NMOMENTS_INPUT
      ELSE
        NMOMENTS = 2 * NSTREAMS - 1
      ENDIF

!  Proxies from bookkeeping

       NBRDF_HALF = Book%NBRDF_HALF
       NSTREAMS_BRDF = Book%NSTREAMS_BRDF
       DO_HALF_RANGE = Book%DO_HALF_RANGE
       DO_LOCAL_WSA  = Book%DO_LOCAL_WSA
       DO_LOCAL_BSA  = Book%DO_LOCAL_BSA

!  Initialise ALL outputs
!  ----------------------

!  Zero Direct-Bounce BRDF

      VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC = ZERO

!  zero the BRDF Fourier components

      VBRDF_Sup_Out%BS_BRDF_F_0 = ZERO
      VBRDF_Sup_Out%BS_BRDF_F   = ZERO
      VBRDF_Sup_Out%BS_USER_BRDF_F_0 = ZERO
      VBRDF_Sup_Out%BS_USER_BRDF_F   = ZERO

!  Initialize surface emissivity
!    Set to zero if you are using Albedo Scaling

      if ( do_wsa_scaling .or. do_bsa_scaling ) then
         VBRDF_Sup_Out%BS_EMISSIVITY      = ZERO
         VBRDF_Sup_Out%BS_USER_EMISSIVITY = ZERO
      else
         VBRDF_Sup_Out%BS_EMISSIVITY      = ONE
         VBRDF_Sup_Out%BS_USER_EMISSIVITY = ONE
      endif

!  Initialize WSA/BSA albedos

      WSA_CALC = zero ; TOTAL_WSA_CALC = zero
      BSA_CALC = zero ; TOTAL_BSA_CALC = zero

!  Fill BRDF arrays
!  ----------------

!  1/31/21. Version 2.8.3. Add doublet geometry flag to all subroutine calls

!    -- 7/28/21. Introduce Half-range azimuth integration flag as an argument (all kernels except NewCM/NewGCM)
!    -- 2/25/22. Version 2.8.5. Introduce DO_FOURIER_USERSUN flag in all calls

!  Start kernel loop

      DO K = 1, N_BRDF_KERNELS

!  Local variables

        BRDF_NPARS = N_BRDF_PARAMETERS(K)
        DO B = 1, MAX_BRDF_PARAMETERS
          PARS(B) = BRDF_PARAMETERS(K,B)
        ENDDO

!  2/25/22. Version 2.8.5. BRDF_FACTOR_LOCAL array needs to be set here

        BRDF_FACTOR_LOCAL(1:NSTOKESSQ) = BRDF_FACTORS(K)
        IF ( DO_SEPARATE_FACTORS(K) .and. NSTOKES.gt.1) THEN
           BRDF_FACTOR_LOCAL(1) = ZERO
           BRDF_FACTOR_LOCAL(2:NSTOKESSQ) = BRDF_VECTOR_FACTORS(K)
        ENDIF

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    Set the refractive index parameataer according to those kernels that need scaling

         IF ( Book%DO_FRESNEL_SCALING(K) ) THEN
           IF ( BRDF_NAMES(K) .EQ. 'Cox-Munk  ' ) THEN
             REFRAC_INDEX = SQRT(BRDF_PARAMETERS(K,2))
           ELSE IF ( BRDF_NAMES(K) .EQ. 'GissCoxMnk' .or. BRDF_NAMES(K) .EQ. 'GCMcomplex' ) then
             REFRAC_INDEX = BRDF_PARAMETERS(K,2)
           ELSE IF ( BRDF_NAMES(K) .EQ. 'BPDF-Soil ' .or. BRDF_NAMES(K) .EQ. 'BPDF-Vegn ' .or. &
                     BRDF_NAMES(K) .EQ. 'BPDF-NDVI ' .or. BRDF_NAMES(K) .EQ. 'ModFresnel') then
             REFRAC_INDEX = BRDF_PARAMETERS(K,1)
           ENDIF        
         ENDIF

!  NewCM/NewGCM Kernels. New for Version 2.7. Upgrades 7/4/15
!  ----------------------------------------------------------

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!     -- Get the refractive index and foam_correction (these operations formerly inside the o lf kernel)
!     -- VBRDF_NewCM_MAKER call now performs just the necessary scaling 

!    Single kernel, Solar Sources only, No Surface Emission. 
!       No MSR (Multiple-surface reflections), No Scaling

!  Sequence is (1) Salinity, WhiteCap, Cox-Munk or  GCM

        IF ( ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX) .or. ( WHICH_BRDF(K) .EQ. NewGCMGLINT_IDX) ) THEN
  
!  Refractive index. Formerly INDWAT

          Call VBRDF_Water_RefracIndex  ( Wavelength, Salinity, Refrac_R, Refrac_I )

!  tempo (fake Cox-Munk)
!          Refrac_R = 1.34d0 ; Refrac_I = 0.0d0
!          Do_FoamOption = .false. ; DO_FacetIsotropy= .true.

!  Foam-reflectance correction.

          WC_Reflectance = zero ; WC_Lambertian = zero
          if ( Do_FoamOption ) then
             call VBRDF_WhiteCap_Reflectance &
               ( WindSpeed, Wavelength, WC_Reflectance, WC_Lambertian )
          endif

!  Whitecap correction to glint

          WC_correction = one - WC_Lambertian

! get the BRDFs
!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

          IF ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX ) then
            CALL VBRDF_NewCM_MAKER &
             ( Core, WC_correction, WC_Reflectance, Refrac_R, K,                           &
               DO_FOURIER_USERSUN, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, &
               DO_DBONLY, NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                           & ! Output
               BRDFUNC_0, USER_BRDFUNC_0  )                                         ! Output
          ELSE
            CALL VBRDF_NewGCM_MAKER &
             ( Core, WC_correction, WC_Reflectance, Refrac_R, K, NSTOKESSQ, &
               DO_FOURIER_USERSUN, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, &
               DO_DBONLY, NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS, &
               DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                            & ! Output
               BRDFUNC_0, USER_BRDFUNC_0  )                                          ! Output
          ENDIF

!  Kokhanovsky-Breon snow model
!  ----------------------------

        ELSE IF ( WHICH_BRDF(K) .eq. SNOWBRDF_Idx ) THEN 

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!     -- get ice refrac index (Log-linear interpolation from Data-set). Formerly inside the old kernel

          wav = PARS(3)
          Logchivals(1:nchidat) = Log(chivals(1:nchidat))
          if ( wav.le.chilams(1) ) then
             chi = chivals(1)
          else if ( wav.ge.chilams(nchidat) ) then
             chi = chivals(nchidat)
          else
            i1 = 0 ; trawl = .true.
            do while ( trawl .and. i1.lt.nchidat-1 )
              i1 = i1 + 1 ; i2 = i1 + 1
              if ( wav.gt.chilams(i1) .and. wav.le.chilams(i2) ) then
                trawl = .false. ; f1 = (chilams(i2) - wav ) / ( chilams(i2) - chilams(i1) ) ; f2 = One - f1
                chi = exp ( f1 * Logchivals(i1) + f2 * Logchivals(i2) )
              endif
            enddo
          endif

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!     -- Get the alpha coefficient (wavelenght dependent).  Formerly inside the old kernel

          WAVMM = wav / 1000.0_fpk            ! convert to [mm]
          Help  = PI4 / wavmm 
          gamma = Help * (chi + PARS(2) * 1.0e-8_fpk ) !gamma = Help * (chi + PARS(2) )
          ALPHA = SQRT ( PARS(1) * gamma)

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    -- Get kernel, Now just a matter of scaling the Core kernels

          CALL VBRDF_SnowKB_MAKER &
            ( Core, Alpha, K,                                                    & ! Main inputs
              DO_WSA_SCALING, DO_BSA_SCALING, DO_HALF_RANGE, DO_FOURIER_USERSUN, & ! New line, Version 2.7
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_DBONLY,    & ! Doublet-geometry flag added, Version 2.8.2
              DO_USER_STREAMS, DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF,       & ! Flags and Numbers
              SCALING_NSTREAMS, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,    & ! Numbers
              DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                           & ! output
              BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,                & ! output
              SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                                ! output, New line, Version 2.7

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    -- Other Kernels with Fresnel Scaling with refractive index
!    -- Bookkeeping inputs have already been set for control
          
        ELSE 

          CALL VBRDF_General_MAKER &
            ( Core, Book%DO_FRESNEL_SCALING(k), Book%DO_FRESNEL_SCALAR(k), REFRAC_INDEX, K, & ! Main inputs
              DO_WSA_SCALING, DO_BSA_SCALING, DO_HALF_RANGE, DO_FOURIER_USERSUN,  & ! New line, Version 2.7
              DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,  DO_DBONLY,   & ! Doublet-geometry flag added, Version 2.8.2
              DO_USER_STREAMS, DO_SURFACE_EMISSION, NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, & ! Flags and Numbers
              SCALING_NSTREAMS, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,    & ! Numbers
              DBKERNEL_BRDFUNC, BRDFUNC, USER_BRDFUNC,                            & ! output
              BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC,                 & ! output
              SCALING_BRDFUNC, SCALING_BRDFUNC_0  )                                 ! output, New line, Version 2.7

        ENDIF

!  Compute Direct Bounce BRDF
!  ==========================

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    -- Rest of calculation as before, EXCEPT now use precalculated bookkeeping quantities, e.g. Book%A_BRDF

!  RobFix 11/8/19. Use masking
!         DO Q = 1, NSTOKESSQ ; O1 = Book%QMASK(Q) instead of DO O1 = 1, NSTOKESSQ

!  Observational Geometry, Optionalities 12/31/12
!    -- 1/31/21. Version 2.8.3. Add doublet geometry option

!  2/25/22. Version 2.8.5. Apply BRDF_FACTOR_LOCAL array 

        IF ( DO_USER_OBSGEOMS ) THEN
          DO Q = 1, NSTOKESSQ ; O1 = Book%QMASK(Q)
            DO IB = 1, NBEAMS
              VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,LUM,LUA,IB) = &
              VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,LUM,LUA,IB) &
              + BRDF_FACTOR_LOCAL(Q) * DBKERNEL_BRDFUNC(O1,LUM,LUA,IB)
            ENDDO
          ENDDO
        ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
          DO Q = 1, NSTOKESSQ ; O1 = Book%QMASK(Q)
            DO IB = 1, NBEAMS
              DO UM = 1, N_USER_STREAMS
                VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,UM,LUA,IB) = &
                 VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,UM,LUA,IB) &
                + BRDF_FACTOR_LOCAL(Q) * DBKERNEL_BRDFUNC(O1,UM,LUA,IB)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO Q = 1, NSTOKESSQ ; O1 = Book%QMASK(Q)
            DO IA = 1, N_USER_RELAZMS
              DO IB = 1, NBEAMS
                DO UM = 1, N_USER_STREAMS
                  VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,UM,IA,IB) = &
                  VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(O1,UM,IA,IB) &
                  + BRDF_FACTOR_LOCAL(Q) * DBKERNEL_BRDFUNC(O1,UM,IA,IB)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF

!  Scaling Section. New code, 15 April 2014 for Version 2.7
!  ========================================================

!  Get the requisite Fourier 0 components
!   7/28/21. Have to introduce the Half-range flag here.

        IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN
          CALL SCALING_FOURIER_ZERO &
            ( DO_LOCAL_WSA, DO_LOCAL_BSA, LAMBERTIAN_KERNEL_FLAG(K),       &
              DO_HALF_RANGE, SCALING_NSTREAMS, NSTREAMS_BRDF, Book%A_BRDF, &
              SCALING_BRDFUNC, SCALING_BRDFUNC_0,                          &
              SCALING_BRDF_F, SCALING_BRDF_F_0 )
        ENDIF

!  White-sky Spherical albedo. Code Upgraded for Version 2.7
!  ---------------------------------------------------------
!mick fix 11/14/2019 - For these WSA computations, the BRDF kernel is assumed to be normalized by 1/Pi
!                      (which is not how the BRDF kernels are defined here in the supplement); thus,
!                      we divide the resulting integral by Pi to account for this) 

        IF ( DO_LOCAL_WSA ) THEN

!  Only for non-Lambertian kernels (trivially = 1 otherwise)

           WSA_CALC(K) = ONE
           IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              HELP_A = ZERO
              DO I = 1, SCALING_NSTREAMS
                 SUM = DOT_PRODUCT(SCALING_BRDF_F(I,1:SCALING_NSTREAMS),Book%SCALING_QUAD_STRMWTS(1:SCALING_NSTREAMS))
                 HELP_A = HELP_A + SUM * Book%SCALING_QUAD_STRMWTS(I)
              ENDDO
              !WSA_CALC(K) = HELP_A * FOUR
              WSA_CALC(K) = HELP_A * FOUR / PIE
           ENDIF
           TOTAL_WSA_CALC = TOTAL_WSA_CALC + BRDF_FACTORS(K) * WSA_CALC(K)

!  Perform consistency check on total white-sky spherical albedo
!    -- This is only done after the kernel summation is finished
!    -- If failed, go to 899 and the error output

           IF ( K.eq.N_BRDF_KERNELS ) then
              if ( TOTAL_WSA_CALC .le. zero ) then
                 STATUS = VLIDORT_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is Negative; examine BRDF Amplitudes'
              else if ( TOTAL_WSA_CALC .gt. one ) then
                 STATUS = VLIDORT_SERIOUS ; NMESSAGES = NMESSAGES + 1
                 MESSAGES(NMESSAGES) = 'Fatal error: Total White-sky albedo is > 1; examine BRDF Amplitudes'
              endif
              IF (STATUS.NE.vlidort_success) GO TO 899
           endif

        ENDIF

!  Black-sky Albedo, only for 1 solar beam. Code Upgraded for Version 2.7
!  ---------------------------------------

!  Compute it for non-Lambertian kernels
!     No check necessary, as the WSA is always checked (regardless of whether scaling is applied)
!mick fix 11/14/2019 - For this BSA computation, the BRDF kernel is assumed to be normalized by 1/Pi
!                      (which is not how the BRDF kernels are defined here in the supplement); thus,
!                      we divide the resulting integral by Pi to account for this) 

        IF (  DO_LOCAL_BSA ) THEN
           BSA_CALC(K) = ONE
           IF ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) ) THEN
              BSA_CALC(K) = TWO * DOT_PRODUCT(SCALING_BRDF_F_0(1:SCALING_NSTREAMS),Book%SCALING_QUAD_STRMWTS(1:SCALING_NSTREAMS))
           ENDIF
           BSA_CALC(K) = BSA_CALC(K) / PIE
           TOTAL_BSA_CALC = TOTAL_BSA_CALC + BRDF_FACTORS(K) * BSA_CALC(K)
        ENDIF

!  @@. Skip Fourier section, if Only the Direct bounce term is calculated.
!   Rob Fix 9/25/14. Direct-bounce flag DO_DBONLY, replaces DO_EXACTONLY

        IF ( DO_DBONLY ) go to 676

!  Fourier Components for the output BRDF arrays
!  =============================================

        DO M = 0, NMOMENTS

!  Call to Fourier routine
!     !@@ Rob Fix, 11/1/19. Bug: Cossin_Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)
!     !@@ Rob Fix, 11/1/19. Remove NSTOKESSQ argument, Cleaner presentation.

!    -- 7/28/21. Introduce Half-range azimuth integration flag as an argument
!    -- 7/28/21. THIS ROUTINE HAS HAD SOME VERY IMPORTANT CHANGES INSIDE

!  2/25/22. Version 2.8.5.  Include flag DO_FOURIER_USERSUN for calculation of LOCAL_USER_BRDF_F_0. 
!  3/28/22, Kernel FACTOR argument has been dropped. 
!       -- [Was only required for emission, factor now deployed outside this routine]

          CALL VBRDF_FOURIER &
            ( DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_USER_STREAMS, DO_SURFACE_EMISSION,  & ! Flags
              LAMBERTIAN_KERNEL_FLAG(K), DO_HALF_RANGE, DO_FOURIER_USERSUN,              & ! Flags
              M, NSTOKES, NBEAMS, NSTREAMS, N_USER_STREAMS, NSTREAMS_BRDF, NBRDF_HALF,   & ! Numbers
              Book%DELFAC(M), Book%BRDF_COSAZMFAC(M,:), Book%BRDF_SINAZMFAC(M,:), Book%A_BRDF, Book%BAX_BRDF,  & ! Surface/Azimuth factors
              BRDFUNC, USER_BRDFUNC, BRDFUNC_0, USER_BRDFUNC_0, EBRDFUNC, USER_EBRDFUNC, & ! Input BRDF Matrices
              LOCAL_BRDF_F, LOCAL_BRDF_F_0, LOCAL_USER_BRDF_F, LOCAL_USER_BRDF_F_0,      & ! Output Fourier components
              LOCAL_EMISSIVITY, LOCAL_USER_EMISSIVITY )                                    ! Output emissivities

!  Start Fourier addition

          IF ( Book%ADD_FOURIER(M,K) ) THEN

!  Rob Fix 11/8/19. Introduce masking
!     That is, replace DO Q = 1, NSTOKESSQ with DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)

!  2/25/22. Version 2.8.5. Use Local BRDF factors
!  Kernel combinations (for quadrature-quadrature reflectance)
!   !@@ Code separated

            DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)
              DO I = 1, NSTREAMS
                DO J = 1, NSTREAMS
                  VBRDF_Sup_Out%BS_BRDF_F(M,Q,I,J) = &
                  VBRDF_Sup_Out%BS_BRDF_F(M,Q,I,J) + BRDF_FACTOR_LOCAL(Q) * LOCAL_BRDF_F(Q,I,J)
if (isnan(LOCAL_BRDF_F(Q,I,J))) print*,' master ', M,Q,I,J, LOCAL_BRDF_F(Q,I,J)
                ENDDO
              ENDDO
            ENDDO

!  Kernel combinations (for Solar-quadrature reflectance)
!   !@@ Solar sources, Optionality 12/31/12

            IF ( DO_SOLAR_SOURCES ) THEN
              DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)
                DO I = 1, NSTREAMS
                  DO IB = 1, NBEAMS
                    VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,I,IB) = &
                    VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,I,IB) + BRDF_FACTOR_LOCAL(Q) * LOCAL_BRDF_F_0(Q,I,IB)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Kernel combinations (for Quadrature-to-Userstream reflectance)
!   !@@ Code separated 12/31/12

            IF ( DO_USER_STREAMS ) THEN
              DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)
                DO UM = 1, N_USER_STREAMS
                  DO J = 1, NSTREAMS
                    VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,UM,J) = &
                    VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,UM,J) + BRDF_FACTOR_LOCAL(Q) * LOCAL_USER_BRDF_F(Q,UM,J)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF

!  Kernel combinations (for Solar-to-Userstream reflectance)
!   !@@ Generally only required for a MS + SS Truncated calculation
!   !@@ Observational Goemetry and Solar sources, Optionalities 12/31/12
!   This is the Truncated Direct bounce calculation. 
!   -- 2/25/22. Version 2.8.5. Optional calculation now controlled by flag DO_FOURIER_USERSUN

            IF ( DO_USER_STREAMS.and.DO_SOLAR_SOURCES.and. DO_FOURIER_USERSUN ) THEN
              IF ( DO_USER_OBSGEOMS ) THEN
                DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)
                  DO IB = 1, NBEAMS
                    VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,LUM,IB) = &
                    VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,LUM,IB) &
                    + BRDF_FACTOR_LOCAL(Q) * LOCAL_USER_BRDF_F_0(Q,LUM,IB)
                  ENDDO
                ENDDO
              ELSE
                DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)
                  DO UM = 1, N_USER_STREAMS
                    DO IB = 1, NBEAMS
                      VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,UM,IB) = &
                      VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,UM,IB) &
                      + BRDF_FACTOR_LOCAL(Q) * LOCAL_USER_BRDF_F_0(Q,UM,IB)
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
            ENDIF

!  Total emissivities
!  Version 2.7, Surely need BRDF factors..........
!  mick fix 9/12/2014 - No, BRDF factors used in defining LOCAL_EMISSIVITY
!  3/28/22. Restore the kernel factors here 

            IF ( DO_SURFACE_EMISSION .AND. M .EQ. 0 ) THEN
              DO Q = 1, NSTOKES
                FF = BRDF_FACTOR_LOCAL(Q)
                DO I = 1, NSTREAMS
                  VBRDF_Sup_Out%BS_EMISSIVITY(Q,I) = &
                  VBRDF_Sup_Out%BS_EMISSIVITY(Q,I) - FF * LOCAL_EMISSIVITY(Q,I)
                ENDDO
                IF ( DO_USER_STREAMS ) THEN
                  DO UI = 1, N_USER_STREAMS
                    VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,UI) = &
                    VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,UI) - FF * LOCAL_USER_EMISSIVITY(Q,UI)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF

!  End Fourier addition

          ENDIF

!  End Fourier loop

        ENDDO

!  continuation point for skipping Fourier work. !@@

676     continue

!  End kernel loop

      ENDDO

!  Now perform normalizations and scaling with White-sky or Black-sky albedos. New section, 02-15 April 2014
!  =========================================================================================================

!  WSABSA OUTPUT. 9/27/14, Added output option.
!                10/16/14. Added Kernel output

      IF ( DO_WSABSA_OUTPUT ) THEN
         VBRDF_Sup_Out%BS_WSA_CALCULATED = TOTAL_WSA_CALC
         VBRDF_Sup_Out%BS_BSA_CALCULATED = TOTAL_BSA_CALC
         VBRDF_Sup_Out%BS_WSA_KERNELS    = WSA_CALC
         VBRDF_Sup_Out%BS_BSA_KERNELS    = BSA_CALC
      ENDIF

!  only if flagged.

      IF ( DO_WSA_SCALING .or. DO_BSA_SCALING ) THEN

!  set scaling factor

         if ( DO_WSA_SCALING ) then
            SCALING = WSA_VALUE / TOTAL_WSA_CALC
         else
            SCALING = BSA_VALUE / TOTAL_BSA_CALC
         endif

!  BRDF Scaling : Start loop over matrix entries
!  ---------------------------------------------

!  RobFix. 11/8/19. Introduce masking
!     That is, replace DO Q = 1, NSTOKESSQ with DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)

         DO O1 = 1, NSTOKESSQ ; Q = Book%QMASK(O1)

!  Exact Direct Beam BRDF  Observational Geometry, Optionalities 12/31/12
!    -- 1/31/21. Version 2.8.3. Add doublet geometry option

            IF ( DO_USER_OBSGEOMS ) THEN
               VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,LUM,LUA,1:NBEAMS) = &
                     SCALING * VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,LUM,LUA,1:NBEAMS)
            ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
               VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,1:N_USER_STREAMS,LUA,1:NBEAMS) = &
                     SCALING * VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,1:N_USER_STREAMS,LUA,1:NBEAMS)
            ELSE
               VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS) = &
                     SCALING * VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(Q,1:N_USER_STREAMS,1:N_USER_RELAZMS,1:NBEAMS)
            ENDIF

!  Fourier terms
!     Rob Fix 9/25/14. Change name to DBONLY

            IF ( .not. DO_DBONLY ) THEN
               DO M = 0, NMOMENTS
!        quadrature-quadrature    reflectance
                  VBRDF_Sup_Out%BS_BRDF_F(M,Q,1:nstreams,1:nstreams) = &
                   SCALING * VBRDF_Sup_Out%BS_BRDF_F(M,Q,1:nstreams,1:nstreams)
!        Solar-quadrature         reflectance
                  IF ( DO_SOLAR_SOURCES ) THEN
                     VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,1:nstreams,1:NBEAMS) = &
                        SCALING * VBRDF_Sup_Out%BS_BRDF_F_0(M,Q,1:nstreams,1:NBEAMS) 
                 ENDIF
!        Quadrature-to-Userstream reflectance
                  IF ( DO_USER_STREAMS ) THEN
                     VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,1:N_USER_STREAMS,1:nstreams) = &
                       SCALING * VBRDF_Sup_Out%BS_USER_BRDF_F(M,Q,1:N_USER_STREAMS,1:nstreams) 
                  ENDIF
!        Solar-to-Userstream      reflectance
                  IF ( DO_USER_STREAMS.and.DO_SOLAR_SOURCES ) THEN
                     IF ( DO_USER_OBSGEOMS ) THEN
                        VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,LUM,1:NBEAMS) = &
                         SCALING * VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,LUM,1:NBEAMS) 
                     ELSE
                       VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,1:N_USER_STREAMS,1:NBEAMS) = &
                         SCALING * VBRDF_Sup_Out%BS_USER_BRDF_F_0(M,Q,1:N_USER_STREAMS,1:NBEAMS) 
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF

!  End reflectance Matrix Loop

         ENDDO

!  Emissivity scaling. Unscaled Emissivity will be < 0
 
         IF ( DO_SURFACE_EMISSION ) THEN
            DO Q = 1, NSTOKES
               VBRDF_Sup_Out%BS_EMISSIVITY(Q,1:nstreams) = &
                    ONE + SCALING * VBRDF_Sup_Out%BS_EMISSIVITY(Q,1:nstreams)
               IF ( DO_USER_STREAMS ) THEN
                  VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,1:N_USER_STREAMS) = &
                       ONE + SCALING * VBRDF_Sup_Out%BS_USER_EMISSIVITY(Q,1:N_USER_STREAMS) 
               ENDIF
            ENDDO
         ENDIF

!  End WSA/BSA scaling option

      ENDIF

!  Continuation point for Error Finish from Consistency Check of WSA or BSA calculated values

899   continue

!  write Exception handling to output structure

      VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT   = STATUS
      VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES = NMESSAGES
      VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES  = MESSAGES

!  Finish

      RETURN
      END SUBROUTINE VBRDF_MAINMASTER

      END MODULE vbrdf_sup_masters_m
