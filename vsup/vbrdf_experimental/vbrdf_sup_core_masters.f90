
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
! # Subroutine in this Module                                   #
! #                                                             #
! #            VBRDF_CORE_MAINMASTER                            #
! #                                                             #
! ###############################################################

      MODULE vbrdf_sup_core_master_m

      PRIVATE
      PUBLIC :: VBRDF_CORE_MAINMASTER

      CONTAINS

!

      SUBROUTINE VBRDF_CORE_MAINMASTER ( VBRDF_Sup_In, Core, Book ) 

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
!     -- NEW. Develops kernel Cores and helper quantities, stored in "Core" type structure
!     -- Also develops and stores wavelengh-independent bookkeeping quantities ("Book")
!     -- This should be called outside of any hyperspectral loop.

! #####################################################################
! #####################################################################

      USE vlidort_pars_m

      USE vbrdf_sup_inputs_def_m
      USE vbrdf_sup_core_def_m

!  Revised Call, 3/17/17

      USE vbrdf_sup_aux_m, only : GETQUAD2,                 &
                                  BRDF_QUADRATURE_Gaussian, &
                                  BRDF_QUADRATURE_Trapezoid

      USE vbrdf_sup_kernels_experimental_m
      USE vbrdf_sup_core_routines_m

!  Implicit none

      IMPLICIT NONE

!  Inputs
!  ------


!  Input structure
!  ---------------

      TYPE(VBRDF_Sup_Inputs), INTENT(IN)   :: VBRDF_Sup_In

!  Output structures
!  -----------------

      TYPE(VBRDF_Sup_Core), INTENT(INOUT) :: Core
      TYPE(VBRDF_Sup_Book), INTENT(INOUT) :: Book

!  VLIDORT local variables
!  +++++++++++++++++++++++

!  Input arguments
!  ===============

!  User stream Control

      LOGICAL ::          DO_USER_STREAMS

!  Surface emission

      LOGICAL ::          DO_SURFACE_EMISSION

!  number of Stokes components

      INTEGER ::          NSTOKES

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

!  WSA and BSA scaling options.
!   Revised, 14-15 April 2014, first introduced 02 April 2014, Version 2.7
!      WSA = White-sky albedo. BSA = Black-sky albedo.
!  Rob Fix 9/27/14. Output variable added

      LOGICAL   :: DO_WSABSA_OUTPUT
      LOGICAL   :: DO_WSA_SCALING
      LOGICAL   :: DO_BSA_SCALING

!  Number of azimuth quadrature streams for BRDF

      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF

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

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: BEAM_SZAS   (MAXBEAMS)
      DOUBLE PRECISION :: USER_RELAZMS(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: USER_ANGLES (MAX_USER_STREAMS)

!  !@@ Local Observational Geometry control and angles

      INTEGER ::          N_USER_OBSGEOMS
      DOUBLE PRECISION :: USER_OBSGEOMS (MAX_USER_OBSGEOMS,3)

!  New Cox-Munk Glint reflectance options (bypasses the usual Kernel system)
!  -------------------------------------------------------------------------

!  Overall flags for this option. 7/4/15 upgrade

      LOGICAL   :: DO_NewCMGLINT, DO_NewGCMGLINT

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      REAL(fpk) :: WINDSPEED, WINDDIR ( MAXBEAMS )

!  Flags for glint shadowing, Foam Correction, facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FoamOption
      LOGICAL   :: DO_FacetIsotropy

!  Local angles, and cosine/sines/weights
!  ======================================

!  Azimuths

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

!  SZAs

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

!  Discrete ordinates (output)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_WEIGHTS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

!  Viewing zenith streams

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  BRDF azimuth quadrature streams

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations
!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: CXE_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Exception handling
!  ==================

!   New code, 02 April 2014. Version 2.7
!     Message Length should be at least 120 Characters

      INTEGER ::             STATUS
      INTEGER ::             NMESSAGES
      CHARACTER (LEN=120) :: MESSAGES ( 0:MAX_MESSAGES )

!  Other local variables
!  =====================

!  Discrete ordinates (local, for Albedo scaling). Version 2.7.

      INTEGER            :: SCALING_NSTREAMS
      DOUBLE PRECISION   :: SCALING_QUAD_STREAMS(MAXSTREAMS_SCALING)
      DOUBLE PRECISION   :: SCALING_QUAD_WEIGHTS(MAXSTREAMS_SCALING)
      DOUBLE PRECISION   :: SCALING_QUAD_SINES  (MAXSTREAMS_SCALING)

!  Local check of Albedo, for all regular kernel options

      LOGICAL :: MODIS_KERNEL_ALONE, DO_CHECK_ALBEDO

!  help. [RobFix 11/8/19, QMask introduced]
!    2/25/22. Version 2.8.5, Introduce Local factors BRDF_FACTOR_LOCAL

      INTEGER ::          K, B, I, IB, UM, IA, O1, M
      INTEGER ::          BRDF_NPARS, NSTOKESSQ
      DOUBLE PRECISION :: MUX, ARGUMENT, PARS ( MAX_BRDF_PARAMETERS )
      LOGICAL ::          DO_LOCAL_WSA, DO_LOCAL_BSA

      INTEGER, PARAMETER :: LUM = 1   !@@
      INTEGER, PARAMETER :: LUA = 1   !@@

!  Default, use Gaussian quadrature

      LOGICAL, PARAMETER :: DO_BRDFQUAD_GAUSSIAN = .true.

!  7/28/21. 8/18/21. Half-range azimuth integration variable introduced. This is the default.
!   -- if in operation, azimuthal range is [0,pi], if not the range is [pi,pi]
!   -- first set as a hard-wired parameter, now determined according to type of kernel
!         (TRUE for all kernels with facet-isotropy, False for NewCM/NewGCM with facet anisotropy

      LOGICAL :: DO_HALF_RANGE

!  Initialize Exception handling
!  -----------------------------

      STATUS = VLIDORT_SUCCESS
      MESSAGES(1:MAX_MESSAGES) = ' '
      NMESSAGES       = 0
      MESSAGES(0)     = 'Successful Execution of VLIDORT Core Master'

!  Copy from input type structure
!  ------------------------------

!  Copy Control inputs

      DO_USER_STREAMS     = VBRDF_Sup_In%BS_DO_USER_STREAMS
      !DO_BRDF_SURFACE     = VBRDF_Sup_In%BS_DO_BRDF_SURFACE
      DO_SURFACE_EMISSION = VBRDF_Sup_In%BS_DO_SURFACE_EMISSION

!  Set number of stokes elements and streams

      NSTOKES  = VBRDF_Sup_In%BS_NSTOKES
      NSTREAMS = VBRDF_Sup_In%BS_NSTREAMS

!  Copy Geometry results

!  !@@ New lines
!  1/31/21. Version 2.8.3. Add doublet geometry flag copy

      DO_SOLAR_SOURCES    = VBRDF_Sup_In%BS_DO_SOLAR_SOURCES
      DO_USER_OBSGEOMS    = VBRDF_Sup_In%BS_DO_USER_OBSGEOMS
      DO_DOUBLET_GEOMETRY = VBRDF_Sup_In%BS_DO_DOUBLET_GEOMETRY

!   !@@ Observational Geometry + Solar sources Optionalities
!   !@@ Either set from User Observational Geometry
!          Or Copy from Usual lattice input
!  -- 1/31/21. Version 2.8.3. Add doublet geometry option

      IF ( DO_USER_OBSGEOMS ) THEN
        N_USER_OBSGEOMS = VBRDF_Sup_In%BS_N_USER_OBSGEOMS
        USER_OBSGEOMS   = VBRDF_Sup_In%BS_USER_OBSGEOMS
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS          = N_USER_OBSGEOMS
          N_USER_STREAMS  = N_USER_OBSGEOMS
          N_USER_RELAZMS  = N_USER_OBSGEOMS
          BEAM_SZAS   (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,1)
          USER_ANGLES (1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
          USER_RELAZMS(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,3)
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS = N_USER_OBSGEOMS
          USER_ANGLES(1:N_USER_OBSGEOMS) = USER_OBSGEOMS(1:N_USER_OBSGEOMS,2)
        ENDIF
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = VBRDF_Sup_In%BS_NBEAMS
          BEAM_SZAS         = VBRDF_Sup_In%BS_BEAM_SZAS
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
          N_USER_RELAZMS    = N_USER_STREAMS
          USER_RELAZMS(1:N_USER_RELAZMS) = VBRDF_Sup_In%BS_USER_RELAZMS    (1:N_USER_RELAZMS)
          USER_ANGLES (1:N_USER_STREAMS) = VBRDF_Sup_In%BS_USER_ANGLES_INPUT(1:N_USER_STREAMS) 
        ELSE
!  NOT ALLOWED
!          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
!          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
!          N_USER_STREAMS = VBRDF_Sup_In%BS_N_USER_STREAMS
!          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ELSE
        IF ( DO_SOLAR_SOURCES ) THEN
          NBEAMS            = VBRDF_Sup_In%BS_NBEAMS
          BEAM_SZAS         = VBRDF_Sup_In%BS_BEAM_SZAS
          N_USER_RELAZMS    = VBRDF_Sup_In%BS_N_USER_RELAZMS
          USER_RELAZMS      = VBRDF_Sup_In%BS_USER_RELAZMS
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ELSE
          NBEAMS         = 1 ; BEAM_SZAS      = ZERO
          N_USER_RELAZMS = 1 ; USER_RELAZMS   = ZERO
          N_USER_STREAMS    = VBRDF_Sup_In%BS_N_USER_STREAMS
          USER_ANGLES = VBRDF_Sup_In%BS_USER_ANGLES_INPUT
        ENDIF
      ENDIF

!  Copy BRDF inputs
!    -- 7/28/21. removed NSTREAMS_BRDF setting, now set later on after the half-range flag is set
!    2/25/22. Version 2.8.5, Introduce Vector factors, and control for this.

      N_BRDF_KERNELS         = VBRDF_Sup_In%BS_N_BRDF_KERNELS
      BRDF_NAMES             = VBRDF_Sup_In%BS_BRDF_NAMES
      WHICH_BRDF             = VBRDF_Sup_In%BS_WHICH_BRDF
      N_BRDF_PARAMETERS      = VBRDF_Sup_In%BS_N_BRDF_PARAMETERS
      BRDF_PARAMETERS        = VBRDF_Sup_In%BS_BRDF_PARAMETERS
      LAMBERTIAN_KERNEL_FLAG = VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG

      DO_SHADOW_EFFECT       = VBRDF_Sup_In%BS_DO_SHADOW_EFFECT

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

!  NewCM options. 7/4/15 upgrade

      DO_NewCMGLINT  = VBRDF_Sup_In%BS_DO_NewCMGLINT
      DO_NewGCMGLINT = VBRDF_Sup_In%BS_DO_NewGCMGLINT

      WINDSPEED = VBRDF_Sup_In%BS_WINDSPEED
      WINDDIR   = VBRDF_Sup_In%BS_WINDDIR

      DO_GlintShadow   = VBRDF_Sup_In%BS_DO_GlintShadow
      DO_FoamOption    = VBRDF_Sup_In%BS_DO_FoamOption
      DO_FacetIsotropy = VBRDF_Sup_In%BS_DO_FacetIsotropy 

!  8/18/21. Set the Half range flag
!mick fix 8/18/2021 - adjusted IF condition to fix DO_HALF_RANGE bug
!                     (however, may not be most efficient fix)

      DO_HALF_RANGE = .true.
      !IF ( DO_NewCMGLINT .or. DO_NewGCMGLINT) THEN
      !  IF ( .not. DO_FacetIsotropy ) DO_HALF_RANGE = .false.
      !ENDIF
      DO K = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(K) .EQ. 'GissCoxMnk'  .OR. &  !10
              BRDF_NAMES(K) .EQ. 'GCMcomplex'  .OR. &  !11
              BRDF_NAMES(K) .EQ. 'BPDF-Soil '  .OR. &  !12
              BRDF_NAMES(K) .EQ. 'BPDF-Vegn '  .OR. &  !13
              BRDF_NAMES(K) .EQ. 'BPDF-NDVI '  .OR. &  !14
              BRDF_NAMES(K) .EQ. 'NewCMGlint'  .OR. &  !15
              BRDF_NAMES(K) .EQ. 'NewGCMGlit'  .OR. &  !16
              BRDF_NAMES(K) .EQ. 'ModFresnel' ) THEN   !18
           DO_HALF_RANGE = .false.
         ENDIF
      ENDDO

!  8/18/21. Move the NSTREAMS_BRDF local setting here, according to the Half-range flag

      IF ( DO_HALF_RANGE ) then
        NSTREAMS_BRDF = VBRDF_Sup_In%BS_NSTREAMS_BRDF / 2
      ELSE
        NSTREAMS_BRDF = VBRDF_Sup_In%BS_NSTREAMS_BRDF
      ENDIF


!  Main code
!  ---------
!mick mod 11/14/2019 - moved defining of these local flags from input section to here 

!  Define some local flags:

!  (1) Albedo check, 7/4/15 upgrade. Condition wrong, 1/4/16
!mick fix 11/14/2019 - added MODIS kernel condition to defining DO_CHECK_ALBEDO

      MODIS_KERNEL_ALONE = .FALSE.
      IF ( N_BRDF_KERNELS .EQ. 1 ) THEN
        IF ( ( WHICH_BRDF(1) .GE. 2 .AND. WHICH_BRDF(1) .LE. 5 ) .OR. WHICH_BRDF(1) .EQ. 7 ) &
          MODIS_KERNEL_ALONE = .TRUE.
      ENDIF

!     DO_CHECK_ALBEDO = ( .not.DO_NewCMGLINT .or.  .not.DO_NewGCMGLINT )
!     DO_CHECK_ALBEDO = ( .not.DO_NewCMGLINT .and. .not.DO_NewGCMGLINT )
!      DO_CHECK_ALBEDO = ( .not.DO_NewCMGLINT .and. .not.DO_NewGCMGLINT .and. .not.MODIS_KERNEL_ALONE )

      DO_CHECK_ALBEDO = .false.               ! Tempo  disable for debug
      DO_WSA_SCALING  = .false.               ! Tempo  disable for debug

!  (2) White & black sky albedo
!  Rob Fix 9/27/14. Output variable included

      DO_LOCAL_WSA = ( DO_WSA_SCALING .or. DO_WSABSA_OUTPUT ) .or.  DO_CHECK_ALBEDO
      DO_LOCAL_BSA = ( DO_BSA_SCALING .or. DO_WSABSA_OUTPUT ) .and. DO_SOLAR_SOURCES

!  Set up Quadrature streams for output
!    QUAD_STRMWTS dropped for Version 2.7 (now redefined for local WSA/BSA scaling)
!    Revised Call, 3/17/17

      CALL GETQUAD2 ( ZERO, ONE, NSTREAMS, QUAD_STREAMS, QUAD_WEIGHTS )
      DO I = 1, NSTREAMS
        QUAD_SINES(I) = SQRT(ONE-QUAD_STREAMS(I)*QUAD_STREAMS(I))
      enddo

!  Set up Quadrature streams for WSA/BSA Scaling. New code, Version 2.7
!    Revised Call, 3/17/17

      Book%SCALING_QUAD_STRMWTS = ZERO
      IF ( DO_LOCAL_WSA .or. DO_LOCAL_BSA ) THEN
         SCALING_NSTREAMS = MAXSTREAMS_SCALING
         CALL GETQUAD2 ( ZERO, ONE, SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_WEIGHTS )
         DO I = 1, SCALING_NSTREAMS
            SCALING_QUAD_SINES(I)   = SQRT(1.0d0-SCALING_QUAD_STREAMS(I)*SCALING_QUAD_STREAMS(I))
            Book%SCALING_QUAD_STRMWTS(I) = SCALING_QUAD_STREAMS(I) * SCALING_QUAD_WEIGHTS(I)
         enddo
      ENDIF

!  Number of Stokes components squared
!    ** Bookkeeping for surface kernel Cox-Munk types
!    ** Only the Giss CoxMunk kernel is vectorized (as of 19 January 2009)

!  Rob Extension 12/2/14. BPDF Kernels (replace BPDF2009)
!  Rob Fix, 14 March 2014. NSTOKESSQ > 1 for BPDF. Now, the BPDF2009 kernel is vectorized
!  Rob Extension 12/2/14. BPDF Kernels (replace BPDF2009)
!    ** Now, the BPDFVEGN, BPDFSOIL, BPDFNDVI kernels are vectorized

!   Additional code for complex RI Giss Cox-Munk, 15 march 2010.
!     3 parameters are PARS(1) = sigma_sq
!                      PARS(2) = Real (RI)
!                      PARS(3) = Imag (RI)

!   Version 2.8. Additional Modified-Fresnel Kernel has 4 parameters. 22 february 2016
!     4 parameters are PARS(1) = Real (RI)
!                      PARS(2) = sigma_sq
!                      PARS(3) = Scaling parameter 
!                      PARS(4) = Shadow Factor 

!  Fix Cox-Munk stuff

      DO K = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(K) .EQ. 'Cox-Munk  ' .OR. &
              BRDF_NAMES(K) .EQ. 'GissCoxMnk' ) THEN
            N_BRDF_PARAMETERS(K) = 3
            IF ( DO_SHADOW_EFFECT ) THEN
              BRDF_PARAMETERS(K,3) = ONE
            ELSE
              BRDF_PARAMETERS(K,3) = ZERO
            ENDIF
         ELSE IF ( BRDF_NAMES(K) .EQ. 'GCMcomplex' ) THEN
            N_BRDF_PARAMETERS(K) = 3
         ENDIF
      ENDDO

!  Half number of moments

      NBRDF_HALF = NSTREAMS_BRDF / 2

!  Usable solar beams. !@@ Optionality, added 12/31/12
!    Warning, this should be the BOA angle. OK for the non-refractive case.

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          MUX =  COS(BEAM_SZAS(IB)*DEG_TO_RAD)
          SZASURCOS(IB) = MUX
          SZASURSIN(IB) = SQRT(1.0D0-MUX*MUX)
        ENDDO
      ELSE
        SZASURCOS = 0.0D0 ; SZASURSIN = 0.0D0
      ENDIF

!  Viewing angles

      DO UM = 1, N_USER_STREAMS
        USER_STREAMS(UM) = COS(USER_ANGLES(UM)*DEG_TO_RAD)
        USER_SINES(UM)   = SQRT(ONE-USER_STREAMS(UM)*USER_STREAMS(UM))
      ENDDO

! Optionality, added 12/31/12
!   Rob Fix 9/25/14. Removed DO_EXACT; Contribution always required now....

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IA = 1, N_USER_RELAZMS
          PHIANG(IA) = USER_RELAZMS(IA)*DEG_TO_RAD
          COSPHI(IA) = COS(PHIANG(IA))
          SINPHI(IA) = SIN(PHIANG(IA))
        ENDDO
      ENDIF

!  Bookkeeping on Refractive index for Fresnel scaling

      Book%DO_FRESNEL_SCALING = .false.
      Book%DO_FRESNEL_SCALAR  = .false.
      DO K = 1, N_BRDF_KERNELS
         IF ( BRDF_NAMES(K) .EQ. 'Cox-Munk  ' ) THEN
            Book%DO_FRESNEL_SCALING(K) = .true.
            Book%DO_FRESNEL_SCALAR(K)  = .true.
         ELSE IF ( BRDF_NAMES(K) .EQ. 'GissCoxMnk' .or. BRDF_NAMES(K) .EQ. 'GCMcomplex' .or. &
                   BRDF_NAMES(K) .EQ. 'BPDF-Soil ' .or. BRDF_NAMES(K) .EQ. 'BPDF-Vegn ' .or. &
                   BRDF_NAMES(K) .EQ. 'BPDF-NDVI ' .or. BRDF_NAMES(K) .EQ. 'ModFresnel') then
            Book%DO_FRESNEL_SCALING(K) = .true.
         ENDIF
      ENDDO

!  BRDF quadrature
!  ---------------

!  Save these quantities for efficient coding
!    -- 7/28/21. Introduce Half-range azimuth integration flag as an argument

      IF ( DO_BRDFQUAD_GAUSSIAN ) then
        CALL BRDF_QUADRATURE_Gaussian &
           ( DO_SURFACE_EMISSION, DO_HALF_RANGE, NSTREAMS_BRDF, NBRDF_HALF, &
             X_BRDF, CX_BRDF, SX_BRDF, Book%A_BRDF, &
             Book%BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ELSE
        CALL BRDF_QUADRATURE_Trapezoid &
           ( DO_SURFACE_EMISSION, DO_HALF_RANGE, NSTREAMS_BRDF, NBRDF_HALF, &
             X_BRDF, CX_BRDF, SX_BRDF, Book%A_BRDF, &
             Book%BAX_BRDF, CXE_BRDF, SXE_BRDF )
      ENDIF

!  Save Bookkeeping stuff

      Book%N_BRDF_PARAMETERS = N_BRDF_PARAMETERS
      Book%NSTREAMS_BRDF = NSTREAMS_BRDF
      Book%NBRDF_HALF    = NBRDF_HALF
      Book%DO_HALF_RANGE = DO_HALF_RANGE
      Book%DO_LOCAL_WSA  = DO_LOCAL_WSA
      Book%DO_LOCAL_BSA  = DO_LOCAL_BSA

!  Reflectivity Mask. Introduced, 11/8/19.

      Book%QMask = 0
      IF ( nstokes.eq.1 ) then
         Book%QMask(1) = 1
      ELSE IF ( nstokes.eq.2 ) then
         Book%QMask(1) = 1 ; Book%QMask(2) = 2
         Book%QMask(3) = 5 ; Book%QMask(4) = 6
      ELSE IF ( nstokes.eq.3 ) then
         Book%QMask(1) = 1 ; Book%QMask(2) = 2  ; Book%QMask(3) = 3
         Book%QMask(4) = 5 ; Book%QMask(5) = 6  ; Book%QMask(6) = 7
         Book%QMask(7) = 9 ; Book%QMask(8) = 10 ; Book%QMask(9) = 11
      ELSE IF ( nstokes.eq.4 ) then
         do o1 = 1, 16
            Book%QMask(o1) = o1
         enddo
      ENDIF

!  Fourier stuff
!  -------------

!  initialize

      Book%ADD_FOURIER     = .false. ; Book%Delfac = ZERO
      Book%BRDF_COSAZMFAC = ZERO ; Book%BRDF_SINAZMFAC = ZERO

!  Only required if doing diffuse reflectance

      IF ( .not. DO_DBONLY ) then
        DO M = 0, 2*NSTREAMS - 1

!  Fourier addition flags

          DO K = 1, N_BRDF_KERNELS
            Book%ADD_FOURIER(M,K) = ( .NOT. LAMBERTIAN_KERNEL_FLAG(K) .OR. &
                                   (LAMBERTIAN_KERNEL_FLAG(K) .AND. M.EQ.0) )
          ENDDO
 
!  surface reflectance factors, Weighted Azimuth factors

          IF ( M .EQ. 0 ) THEN
            Book%DELFAC(M)  = ONE
            DO I = 1, NSTREAMS_BRDF
              Book%BRDF_COSAZMFAC(M,I) = Book%A_BRDF(I)
              Book%BRDF_SINAZMFAC(M,I) = ZERO
            ENDDO
          ELSE
            Book%DELFAC(M)   = TWO
            DO I = 1, NSTREAMS_BRDF
              ARGUMENT = DBLE(M) * X_BRDF(I)
              Book%BRDF_COSAZMFAC(M,I) = Book%A_BRDF(I) * COS ( ARGUMENT )
              Book%BRDF_SINAZMFAC(M,I) = Book%A_BRDF(I) * SIN ( ARGUMENT )
            ENDDO
          ENDIF

        ENDDO
      ENDIF

!  Fill BRDF CORE and HELP arrays
!  ==============================

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

!  Lambertian kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. LAMBERTIAN_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( LAMBERTIAN_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Ross thin kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHIN_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( ROSSTHIN_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Ross thick kernel, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROSSTHICK_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( ROSSTHICK_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  2/22/16. V2.8, add Ross-Thick Alternative with Hotspot, (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RTKHOTSPOT_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( ROSSTHICK_ALT_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Li Sparse kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LISPARSE_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( LISPARSE_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Li Dense kernel; 2 free parameters

        IF ( WHICH_BRDF(K) .EQ. LIDENSE_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( LIDENSE_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Hapke kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. HAPKE_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( HAPKE_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Rahman kernel (3 free parameters)

        IF ( WHICH_BRDF(K) .EQ. RAHMAN_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( RAHMAN_CORE, K, & 
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Roujean kernel (0 free parameters)

        IF ( WHICH_BRDF(K) .EQ. ROUJEAN_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( ROUJEAN_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Scalar-only original Cox-Munk kernel: (2 free parameters, Shadow = Third)
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. COXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
          CALL VBRDF_CORE_MAKER ( COXMUNK_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output

        ENDIF

!  GISS Vector Cox-Munk kernel: (2 free parameters, Shadow = Third). Real RI only
!    Distinguish between MS case.....

        IF ( WHICH_BRDF(K) .EQ. GISSCOXMUNK_IDX ) THEN
          IF ( DO_SHADOW_EFFECT ) PARS(3) = 1.0d0
          CALL VBRDF_CORE_MAKER ( GISSCOXMUNK_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Rob Extension 12/2/14. BPDF Kernels (replaces BPDF 2009 kernel)

        IF ( WHICH_BRDF(K) .EQ. BPDFVEGN_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( BPDFVEGN_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Rob Extension 12/2/14. BPDF Kernels (replaces BPDF 2009 kernel)

        IF ( WHICH_BRDF(K) .EQ. BPDFSOIL_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( BPDFSOIL_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  Rob Extension 12/2/14. BPDF Kernels (replaces BPDF 2009 kernel)

        IF ( WHICH_BRDF(K) .EQ. BPDFNDVI_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( BPDFNDVI_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  2/22/16. Add option for Modified-Fresnel Kernel. (4 free parameters)

        IF ( WHICH_BRDF(K) .EQ. MODFRESNEL_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( MODFRESNEL_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  1/31/21, Version 2.8.3. Kernel Call for Analytical Model for Snow BRDF.
!     -- First introduced to VLIDORT, 18 November 2020.
!     -- Input index is SNOWBRDF_IDX. Output SNOWMODELBRDF_CORE

        IF ( WHICH_BRDF(K) .EQ. SNOWBRDF_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( SNOWMODELBRDF_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter
!     -- Input index is SNOWRTS_IDX. Output SNOWRTSMODEL_CORE

        IF ( WHICH_BRDF(K) .EQ. SNOWRTS_IDX ) THEN
          CALL VBRDF_CORE_MAKER ( SNOWRTSMODEL_CORE, K, &
               DO_LOCAL_WSA, DO_LOCAL_BSA, DO_HALF_RANGE, DO_FOURIER_USERSUN,    & ! New line, Version 2.7
               DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,          &
               DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                  &
               NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                 &
               NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                 &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,               &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, PARS,               &
               SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,       & ! New line, Version 2.7
               X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                     &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE,                 & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP,                 & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE, Core%EBRDFUNC_CORE, Core%USER_EBRDFUNC_CORE, & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP, Core%EBRDFUNC_HELP, Core%USER_EBRDFUNC_HELP, & ! Help output
               Core%SCALING_BRDFUNC_CORE, Core%SCALING_BRDFUNC_0_CORE,                           & ! Core output
               Core%SCALING_BRDFUNC_HELP, Core%SCALING_BRDFUNC_0_HELP  )                           ! Help output
        ENDIF

!  NewCM Kernel. New for Version 2.7. Upgrades 7/4/15
!  --------------------------------------------------

!    Single kernel, Solar Sources only, No Surface Emission. 
!       No MSR (Multiple-surface reflections), No Scaling

!  Sequence is (1) Salinity, WhiteCap, Cox-Munk or  GCM

        IF ( ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX) .or. ( WHICH_BRDF(K) .EQ. NewGCMGLINT_IDX) ) THEN

!  Reverse angle effect !!!!!
!   NewCM Convention is opposite from VLIDORT, that is, Phi(6S) = 180 - Phi(VL)
!   Once this is realized, NewCM and Regular Cox-Munk will agree perfectly
!          ( Have to turn off Whitecaps in 6S and use Facet Isotropy, make sure RI same)
!   Rob Fix 9/25/14. Removed DO_EXACT

          IF ( DO_SOLAR_SOURCES ) THEN
             DO IA = 1, N_USER_RELAZMS
                PHIANG(IA) = PIE - PHIANG(IA)
                COSPHI(IA) = - COSPHI(IA)
             ENDDO
             DO I = 1, NSTREAMS_BRDF
                CX_BRDF(I) = - CX_BRDF(I)
             ENDDO
          ENDIF
  
!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

          IF ( WHICH_BRDF(K) .EQ. NewCMGLINT_IDX ) then
            CALL VBRDF_NewCM_CORE_MAKER ( K, &
               DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, DO_FOURIER_USERSUN, &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE, & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP, & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE,                         & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP )                          ! Help output
          ELSE IF ( WHICH_BRDF(K) .EQ. NewGCMGLINT_IDX ) then
            CALL VBRDF_NewGCM_CORE_MAKER ( K, &
               DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, DO_FOURIER_USERSUN, &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               Core%DBKERNEL_BRDFUNC_CORE, Core%BRDFUNC_CORE, Core%USER_BRDFUNC_CORE, & ! Core output
               Core%DBKERNEL_BRDFUNC_HELP, Core%BRDFUNC_HELP, Core%USER_BRDFUNC_HELP, & ! Help output
               Core%BRDFUNC_0_CORE, Core%USER_BRDFUNC_0_CORE,                         & ! Core output
               Core%BRDFUNC_0_HELP, Core%USER_BRDFUNC_0_HELP )                          ! Help output
          ENDIF

        ENDIF

!  End kernel loop

      ENDDO

!  End subroutine

      RETURN
      END SUBROUTINE VBRDF_CORE_MAINMASTER

!  End Module

END MODULE vbrdf_sup_core_master_m
