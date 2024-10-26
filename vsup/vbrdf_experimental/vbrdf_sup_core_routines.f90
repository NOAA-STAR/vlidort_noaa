
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
! # General Subroutines in this Module                          #
! #                                                             #
! #              VBRDF_MAKER                                    #
! #              VBRDF_GCMCRI_MAKER                             #
! #              SCALING_FOURIER_ZERO (New, Version 2.7)        #
! #              VBRDF_FOURIER                                  #
! #                                                             #
! # New Cox-Munk Subroutine in this Module  (Version 2.7)       #
! #                                                             #
! #              VBRDF_NewCM_MAKER                              #
! #              VBRDF_NewGCM_MAKER                             #
! #                                                             #
! ###############################################################

!  Changes for Version 2.8
!  -----------------------

!  Allow for presence of Two new kernels. Set-up calls for them,

!  Version 2p8p1. 2019 Overhaul.
!  -----------------------------

!  In the Fourier routines--->

!     !@@ Rob Fix, 3/20/19. Bug: Lambertian case wrongly used NSTOKESSQ instead of 1
!     !@@ Rob Fix, 11/1/19. Bug: Cossin_Mask wrongly used, code replaced. Thanks to X.Xu (UMBC)
!     !@@ Rob Fix, 11/1/19. Recoded the BRDF emissivity case
!     !@@ Rob Fix, 11/1/19. Used the Dot-Product command throughout. NSTOKESSQ Dropped.

!  Version 2p8p3.
!  --------------

!  1/31/21, Version 2.8.3. Analytical Model for Snow BRDF.
!     -- New VLIDORT BRDF Kernel. First introduced to VLIDORT, 18 November 2020.
!     -- Kokhanovsky and Breon, IEEE GeoScience and Remote Sensing Letters, Vol 9(5), 928-932 (2012)
!     -- The three parameters (L and M are free parameters) are
!        1. The L-value, related to the snow grain diameter size. Units [mm]
!        2. The M-value, "directly proportional to the mass concentration of pollutants"
!        3. The Wavelength in Microns --> Imaginary part of the refractive index.

!  2/28/21, Version 3.8.3. Doublet geometry option added

!  7/28/21. Version 2.8.3. Some Changes
!    -- Half-Range azimuth integration introduced, now the default option. Controlled by parameter setting
!    -- Adjustment of +/- signs in Fourier-component sine-series settings. Validates against TestBed
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/BAX/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Basic model is the ART (Asymptotic radiative transfer) model of Kokhanovsky and Breon
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!     -- NEW. Core Kernel Maker routines, including NewCM and NewGCM. Stand-Alone routines

      MODULE vbrdf_sup_core_routines_m

      USE vbrdf_sup_kernels_experimental_m, only : GenGlint_CORE, GenGlint_GCM_CORE

      PRIVATE
      PUBLIC :: VBRDF_CORE_MAKER, VBRDF_NewCM_CORE_MAKER, VBRDF_NewGCM_CORE_MAKER

      CONTAINS

      SUBROUTINE VBRDF_CORE_MAKER ( BRDF_VFUNCTION_CORE, NNN,                 &
           DO_WSA_SCALING, DO_BSA_SCALING, DO_HALF_RANGE, DO_FOURIER_USERSUN, & ! New line, Version 2.7
           DO_SOLAR_SOURCES, DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY,           & ! Doublet-geometry flag added, Version 2.8.2
           DO_DBONLY, DO_USER_STREAMS, DO_SURFACE_EMISSION,                   &
           NSTREAMS_BRDF, NBRDF_HALF, NSTOKESSQ, BRDF_NPARS,                  &
           NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,                  &
           QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
           SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI, BRDF_PARS,           &
           SCALING_NSTREAMS, SCALING_QUAD_STREAMS, SCALING_QUAD_SINES,        & ! New line, Version 2.7
           X_BRDF, CX_BRDF, SX_BRDF, CXE_BRDF, SXE_BRDF,                      &
           DBKERNEL_BRDFUNC_CORE, BRDFUNC_CORE, USER_BRDFUNC_CORE,                 & ! Core output
           DBKERNEL_BRDFUNC_HELP, BRDFUNC_HELP, USER_BRDFUNC_HELP,                 & ! Help output
           BRDFUNC_0_CORE, USER_BRDFUNC_0_CORE, EBRDFUNC_CORE, USER_EBRDFUNC_CORE, & ! Core output
           BRDFUNC_0_HELP, USER_BRDFUNC_0_HELP, EBRDFUNC_HELP, USER_EBRDFUNC_HELP, & ! Help output
           SCALING_BRDFUNC_CORE, SCALING_BRDFUNC_0_CORE,                           & ! Core output
           SCALING_BRDFUNC_HELP, SCALING_BRDFUNC_0_HELP  )                           ! Help output

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  7/28/21. Some changes
!    -- Half-Range azimuth integration control input introduced (DO_HALF_RANGE)
!    -- MAXSTHALF_BRDF Dimensioning for CXE/SXE/EBRDFUNC/USER_EBRDFUNC arrays, replaced by MAXSTREAMS_BRDF

!  include file of dimensions and numbers
!    -- 7/28/21.  MAXSTHALF_BRDF Dimensioning dropped

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      USE VLIDORT_PARS_m, only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                 MAXSTOKES_SQ, MAXSTREAMS, MAXSTREAMS_SCALING, & 
                                 MAX_BRDF_PARAMETERS, MAXSTREAMS_BRDF, ZERO

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Observational Geometry Inputs. Marked with !@@
!     Installed 31 december 2012.
!     Observation-Geometry input control.         (DO_USER_OBSGEOMS)
!     Added solar_sources flag for better control (DO_SOLAR_SOURCES)
!     Added Overall-exact flag for better control (DO_EXACT)

!  Input arguments
!  ===============

!  BRDF functions (external calls)

      EXTERNAL         BRDF_VFUNCTION_CORE

!  Kernel index

      INTEGER ::          NNN

!  White-sky and Black-sky albedo scaling flags. New Version 2.7

      LOGICAL ::          DO_WSA_SCALING
      LOGICAL ::          DO_BSA_SCALING

!  7/28/21. Half-range integration variable introduced
!   -- if set, azimuthal range is [0,pi], if not set the range is [pi,pi]

      LOGICAL ::          DO_HALF_RANGE

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      LOGICAL ::          DO_FOURIER_USERSUN

!  Solar sources + Observational Geometry flag
!    - 1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_SOLAR_SOURCES
      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Local flags

      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_SURFACE_EMISSION

!  Number of Azimuth waudrature streams

      INTEGER ::          NSTREAMS_BRDF
      INTEGER ::          NBRDF_HALF

!  Local number of Stokes component matrix entries
!    value = 1 for most kernels, except GISS Cox-Munk

      INTEGER ::          NSTOKESSQ

!  Local number of Kernel parameters

      INTEGER ::          BRDF_NPARS

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  Discrete ordinates (local, for Albedo scaling). New Version 2.7

      INTEGER          :: SCALING_NSTREAMS
      DOUBLE PRECISION :: SCALING_QUAD_STREAMS(MAXSTREAMS_SCALING)
      DOUBLE PRECISION :: SCALING_QUAD_SINES  (MAXSTREAMS_SCALING)

!  Local parameter array

      DOUBLE PRECISION :: BRDF_PARS ( MAX_BRDF_PARAMETERS )

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: CXE_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SXE_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF CORE functions
!  ==========================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC_CORE   ( 4,    MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_CORE ( 4,    MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_HELP   ( 4, 7, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_HELP ( 4, 7, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC_CORE   ( 4,    MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_CORE ( 4,    MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_HELP   ( 4, 7, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_HELP ( 4, 7, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_CORE ( 4,    MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_HELP ( 4, 7, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity
!   -- 7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC_CORE ( 4,   MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: EBRDFUNC_HELP ( 4, 7, MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC_CORE ( 4,   MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC_HELP ( 4, 7, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC_CORE   ( 4, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0_CORE ( 4, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_HELP   ( 4, 7, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0_HELP ( 4, 7, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

!  local variables
!  ---------------

      INTEGER ::            I, UI, J, K, KE, IB, NSQ, NKE_RANGE
      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!  Exact DB calculation
!  --------------------

!    !@@ Observational Geometry choice 12/31/12
!    !@@ Rob Fix 9/25/14. Always calculated for Solar Sources
!mick fix 9/19/2017 - initialize DBKERNEL_BRDFUNC

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      IF ( DO_SOLAR_SOURCES ) THEN
        IF ( DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
             CALL BRDF_VFUNCTION_CORE &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,   &
                 SZASURCOS(IB), SZASURSIN(IB),                 &
                 USER_STREAMS(IB), USER_SINES(IB), PHIANG(IB), &
                 COSPHI(IB), SINPHI(IB),                       &
                 DBKERNEL_BRDFUNC_CORE(NNN,LUM,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,LUM,LUA,IB) )
          ENDDO
        ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
          DO IB = 1, NBEAMS
             DO UI = 1, N_USER_STREAMS
                CALL BRDF_VFUNCTION_CORE &
                   ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,   &
                     SZASURCOS(IB), SZASURSIN(IB),                 &
                     USER_STREAMS(UI), USER_SINES(UI), PHIANG(UI), &
                     COSPHI(UI), SINPHI(UI),                       &
                     DBKERNEL_BRDFUNC_CORE(NNN,UI,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,LUA,IB) )
              ENDDO
           ENDDO
        ELSE
          DO K = 1, N_USER_RELAZMS
             DO IB = 1, NBEAMS
                DO UI = 1, N_USER_STREAMS
                   CALL BRDF_VFUNCTION_CORE &
                    ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                      SZASURCOS(IB), SZASURSIN(IB),                 &
                      USER_STREAMS(UI), USER_SINES(UI), PHIANG(K), &
                      COSPHI(K), SINPHI(K),                        &
                      DBKERNEL_BRDFUNC_CORE(NNN,UI,K,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,K,IB)  )
                ENDDO
             ENDDO
          ENDDO
        END IF
      ENDIF

!  SCALING OPTIONS (New Section, Version 2.7)
!  ------------------------------------------

!  White-sky albedo, scaling. Only requires the (1,1) component
!     Use Local "Scaling_streams", both incident and outgoing

      IF ( DO_WSA_SCALING ) THEN
         NSQ = 1
         DO I = 1, SCALING_NSTREAMS
            DO J = 1, SCALING_NSTREAMS
               DO K = 1, NSTREAMS_BRDF
                  CALL BRDF_VFUNCTION_CORE &
                     ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,      &
                       SCALING_QUAD_STREAMS(J), SCALING_QUAD_SINES(J),  &
                       SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I),  &         
                       X_BRDF(K), CX_BRDF(K), SX_BRDF(K),               &
                       SCALING_BRDFUNC_CORE(NNN,I,J,K), SCALING_BRDFUNC_HELP(NNN,1,I,J,K) )
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!  Black-sky albedo, scaling. Only requires the (1,1) component
!     Use Local "Scaling_streams" for outgoing, solar beam for incoming (IB = 1)

      IF ( DO_BSA_SCALING .and. DO_SOLAR_SOURCES ) THEN
         IB = 1 ; NSQ = 1
         DO I = 1, SCALING_NSTREAMS
            DO K = 1, NSTREAMS_BRDF
               CALL BRDF_VFUNCTION_CORE &
                   ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,     &
                     SZASURCOS(IB), SZASURSIN(IB),                   &
                     SCALING_QUAD_STREAMS(I), SCALING_QUAD_SINES(I), &
                     X_BRDF(K), CX_BRDF(K), SX_BRDF(K),              &
                     SCALING_BRDFUNC_0_CORE(NNN,I,K), SCALING_BRDFUNC_0_HELP(NNN,1,I,K) )
            ENDDO
         ENDDO
      ENDIF

!  Return if the Direct-bounce BRDF is all that is required (scaled or not!)

      IF ( DO_DBONLY ) RETURN

!  Quadrature outgoing directions
!  ------------------------------

!  Incident Solar beam
!    !@@  Solar Optionality. 12/31/12

      IF ( DO_SOLAR_SOURCES ) THEN
        DO IB = 1, NBEAMS
          DO I = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_CORE &
                ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,            &
                  SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),         &
                  QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                  BRDFUNC_0_CORE(NNN,I,IB,K), BRDFUNC_0_HELP(NNN,1,I,IB,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  incident quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL BRDF_VFUNCTION_CORE &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,            &
                 QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),       &
                 QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),      &
                 BRDFUNC_CORE(NNN,I,J,K), BRDFUNC_HELP(NNN,1,I,J,K) )
          ENDDO
        ENDDO
      ENDDO

!  Emissivity (optional) - BRDF quadrature input directions
!    - 7/28/21. Use the HALF_RANGE flag to set the number of azimuth points NKE_RANGE.

      IF ( DO_SURFACE_EMISSION ) THEN
        NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
        DO I = 1, NSTREAMS
          DO KE = 1, NKE_RANGE
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_CORE &
               ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,                 &
                 CXE_BRDF(KE), SXE_BRDF(KE), QUAD_STREAMS(I), QUAD_SINES(I), &
                 X_BRDF(K),  CX_BRDF(K), SX_BRDF(K),                         &
                 EBRDFUNC_CORE(NNN,I,KE,K), EBRDFUNC_HELP(NNN,1,I,KE,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  User-streams outgoing directions
!  --------------------------------

      IF ( DO_USER_STREAMS ) THEN

!  Incident Solar beam, Outgoing User-stream
!    !@@ Observational Geometry choice + Solar Optionality. 12/31/12

!  2/25/22. Version 2.8.5. This is now optional (add do_FOURIER_USERSUN flag)
 
        IF (DO_SOLAR_SOURCES .and. do_FOURIER_USERSUN) THEN
          IF ( DO_USER_OBSGEOMS ) THEN
            DO IB = 1, NBEAMS
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_VFUNCTION_CORE &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,            &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),        &
                   USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                   USER_BRDFUNC_0_CORE(NNN,LUM,IB,K), USER_BRDFUNC_0_HELP(NNN,1,LUM,IB,K) )
              ENDDO
            ENDDO
          ELSE
            DO IB = 1, NBEAMS
             DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_VFUNCTION_CORE &
                   ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,            &
                     SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),        &
                     USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K),     &
                     USER_BRDFUNC_0_CORE(NNN,UI,IB,K), USER_BRDFUNC_0_HELP(NNN,1,UI,IB,K) )
              ENDDO
             ENDDO
            ENDDO
          ENDIF
        ENDIF

!  incident quadrature directions, Outgoing User-stream

        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL BRDF_VFUNCTION_CORE &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                   QUAD_STREAMS(J), QUAD_SINES(J),   &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), &
                   CX_BRDF(K), SX_BRDF(K),                      &
                   USER_BRDFUNC_CORE(NNN,UI,J,K), USER_BRDFUNC_HELP(NNN,1,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

!  Emissivity (optional) - BRDF quadrature input directions
!  7/28/21. Use the HALF_RANGE flag to set the number of azimth points NKE_RANGE.

        IF ( DO_SURFACE_EMISSION ) THEN
          NKE_RANGE = NSTREAMS_BRDF ; if (.not. DO_HALF_RANGE ) NKE_RANGE = NBRDF_HALF
          DO UI = 1, N_USER_STREAMS
            DO KE = 1, NKE_RANGE
              DO K = 1, NSTREAMS_BRDF
                CALL BRDF_VFUNCTION_CORE &
                 ( MAX_BRDF_PARAMETERS, BRDF_NPARS, BRDF_PARS,  &
                   CXE_BRDF(KE), SXE_BRDF(KE),       &
                   USER_STREAMS(UI), USER_SINES(UI), X_BRDF(K), &
                   CX_BRDF(K), SX_BRDF(K),                      &
                   USER_EBRDFUNC_CORE(NNN,UI,KE,K), USER_EBRDFUNC_HELP(NNN,1,UI,KE,K) )
              ENDDO
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_CORE_MAKER

! 

      SUBROUTINE VBRDF_NewCM_CORE_MAKER ( NNN, &
               DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, DO_FOURIER_USERSUN, &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS,   &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               DBKERNEL_BRDFUNC_CORE, BRDFUNC_CORE, USER_BRDFUNC_CORE,            & ! Core output
               DBKERNEL_BRDFUNC_HELP, BRDFUNC_HELP, USER_BRDFUNC_HELP,            & ! Help output
               BRDFUNC_0_CORE, USER_BRDFUNC_0_CORE,                               & ! Core output
               BRDFUNC_0_HELP, USER_BRDFUNC_0_HELP )                                ! Help output

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m,      only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                      MAXSTREAMS, MAXSTREAMS_BRDF, ZERO, ONE, DEG_TO_RAD

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  Kernel index

      INTEGER   :: NNN

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      DOUBLE PRECISION::  WINDSPEED, WINDDIR ( MAXBEAMS )

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      LOGICAL ::          DO_FOURIER_USERSUN

!  Local flags
!    1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Number of Azimuth quadrature streams

      INTEGER ::          NSTREAMS_BRDF

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC_CORE   ( 4,    MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_CORE ( 4,    MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_HELP   ( 4, 7, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_HELP ( 4, 7, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC_CORE   ( 4,    MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_CORE ( 4,    MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_HELP   ( 4, 7, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_HELP ( 4, 7, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_CORE ( 4,    MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_HELP ( 4, 7, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  local variables
!  ---------------

      LOGICAL          :: DO_COEFFS, Local_Isotropy
      INTEGER          :: I, UI, J, K, IB
      DOUBLE PRECISION :: PHI_W(MAXBEAMS), CPHI_W(MAXBEAMS), SPHI_W(MAXBEAMS)
      DOUBLE PRECISION :: SUNGLINT_COEFFS(7)

      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
         DO IB = 1, nbeams
            PHI_W(IB)  = WINDDIR(IB)
            CPHI_W(IB) = cos(WINDDIR(IB) * deg_to_rad) 
            SPHI_W(IB) = sin(WINDDIR(IB) * deg_to_rad)
         ENDDO
      endif

!  Direct Bounce calculation
!  -------------------------
!mick fix 9/19/2017  - initialize DBKERNEL_BRDFUNC
!mick mod 11/14/2019 - trimmed dimensions initialized

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      IF ( .NOT. DO_USER_OBSGEOMS ) THEN
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO  UI = 1, N_USER_STREAMS
                  CALL GenGlint_CORE &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
                      WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),    &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  &
                      USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
                      SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,UI,K,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,K,IB) )
               ENDDO
            ENDDO
         ENDDO
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
         DO IB = 1, NBEAMS
            DO  UI = 1, N_USER_STREAMS
               CALL GenGlint_CORE &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
                   WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),       &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     &
                   USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), &
                   SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,UI,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,LUA,IB) )
            ENDDO
         ENDDO
      ELSE
         DO IB = 1, NBEAMS
            CALL GenGlint_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
               WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),       &
               SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     &
               USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), &
               SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,LUM,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,LUM,LUA,IB) )
         ENDDO
      ENDIF

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================
!mick fix 11/14/2019  - initialize BRDFUNC_0, USER_BRDFUNC_0

!  Quadrature outgoing directions

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GenGlint_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),     &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),    &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, BRDFUNC_0_CORE(NNN,I,IB,K), BRDFUNC_0_HELP(NNN,1,I,IB,K) )
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation
!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      IF ( DO_USER_STREAMS .and. DO_FOURIER_USERSUN ) THEN
        IF (.NOT. DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL GenGlint_CORE &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                   WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),      &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),    &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   SUNGLINT_COEFFS, USER_BRDFUNC_0_CORE(NNN,UI,IB,K), USER_BRDFUNC_0_HELP(NNN,1,UI,IB,K))
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GenGlint_CORE &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),      &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),    &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, USER_BRDFUNC_0_CORE(NNN,UI,LUM,K), USER_BRDFUNC_0_HELP(NNN,1,UI,LUM,K))
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!  Incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================
!mick fix 11/14/2019  - initialize BRDFUNC, USER_BRDFUNC

!   Can only be treated with 1 Wind direction.....
!     if ( NBEAMS > 1) MUST assume local Facet Isotropy
!            --> set up Local Wind-direction and re-set coefficients flag.
!     if ( NBAMS = 1 ) use the first wind direction, no need to re-calculate coefficients

      if ( NBEAMS .gt. 1 ) then
         local_Isotropy = .true.
!         local_Isotropy = .false.  ! Bug 9/27/14, fixed
         PHI_W      = zero 
         CPHI_W     = one
         SPHI_W     = zero
         DO_COEFFS  = .true.
      endif
 
!  Outgoing quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GenGlint_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               WINDSPEED, PHI_W(1), CPHI_W(1), SPHI_W(1),        &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, BRDFUNC_CORE(NNN,I,J,K), BRDFUNC_HELP(NNN,1,I,J,K)  )
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GenGlint_CORE &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 WINDSPEED, PHI_W(1), CPHI_W(1), SPHI_W(1),         &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, USER_BRDFUNC_CORE(NNN,UI,J,K), USER_BRDFUNC_HELP(NNN,1,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO

      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_NewCM_CORE_MAKER

!

      SUBROUTINE VBRDF_NewGCM_CORE_MAKER ( NNN, &
               DO_GlintShadow, DO_FacetIsotropy, WINDSPEED, WINDDIR, DO_FOURIER_USERSUN, &
               DO_USER_OBSGEOMS, DO_DOUBLET_GEOMETRY, DO_USER_STREAMS, DO_DBONLY, &
               NSTREAMS_BRDF, NSTREAMS, NBEAMS, N_USER_STREAMS, N_USER_RELAZMS  , &
               QUAD_STREAMS, QUAD_SINES, USER_STREAMS, USER_SINES,                &
               SZASURCOS, SZASURSIN, PHIANG, COSPHI, SINPHI,                      &
               X_BRDF, CX_BRDF, SX_BRDF,                                          &
               DBKERNEL_BRDFUNC_CORE, BRDFUNC_CORE, USER_BRDFUNC_CORE,            & ! Core output
               DBKERNEL_BRDFUNC_HELP, BRDFUNC_HELP, USER_BRDFUNC_HELP,            & ! Help output
               BRDFUNC_0_CORE, USER_BRDFUNC_0_CORE,                               & ! Core output
               BRDFUNC_0_HELP, USER_BRDFUNC_0_HELP )                                ! Help output

!  1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag added to argument list

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

!  include file of dimensions and numbers

      USE VLIDORT_PARS_m,      only : MAXBEAMS, MAX_USER_RELAZMS, MAX_USER_STREAMS, &
                                      MAXSTREAMS, MAXSTREAMS_BRDF, zero, one, DEG_TO_RAD

      IMPLICIT NONE

!  Prepares the bidirectional reflectance scatter matrices

!  Input arguments
!  ===============

!  Kernel index

      INTEGER   :: NNN

!  NewCM Glitter options (bypasses the usual Kernel system)
!  -------------------------------------------------------

!  Flags for glint shadowing, Facet Isotropy

      LOGICAL   :: DO_GlintShadow
      LOGICAL   :: DO_FacetIsotropy

!  Input Wind speed in m/s, and azimuth directions relative to Sun positions

      DOUBLE PRECISION:: WINDSPEED, WINDDIR ( MAXBEAMS )

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      LOGICAL ::          DO_FOURIER_USERSUN

!  Local flags
!    1/31/21. Version 2.8.3. DO_DOUBLET_GEOMETRY flag 

      LOGICAL ::          DO_USER_OBSGEOMS
      LOGICAL ::          DO_USER_STREAMS
      LOGICAL ::          DO_DOUBLET_GEOMETRY

!  Rob Fix 9/25/14. Two variables replaced
!   Flag for the Direct-bounce term, replaces former "EXACT" variables 
!   Exact flag (!@@) and Exact only flag --> no Fourier term calculations
!      LOGICAL ::          DO_EXACT
!      LOGICAL ::          DO_EXACTONLY
      LOGICAL ::          DO_DBONLY

!  Number of Azimuth quadrature streams

      INTEGER ::          NSTREAMS_BRDF

!  Local angle control

      INTEGER ::          NSTREAMS
      INTEGER ::          NBEAMS
      INTEGER ::          N_USER_STREAMS
      INTEGER ::          N_USER_RELAZMS

!  Local angles

      DOUBLE PRECISION :: PHIANG(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: COSPHI(MAX_USER_RELAZMS)
      DOUBLE PRECISION :: SINPHI(MAX_USER_RELAZMS)

      DOUBLE PRECISION :: SZASURCOS(MAXBEAMS)
      DOUBLE PRECISION :: SZASURSIN(MAXBEAMS)

      DOUBLE PRECISION :: QUAD_STREAMS(MAXSTREAMS)
      DOUBLE PRECISION :: QUAD_SINES  (MAXSTREAMS)

      DOUBLE PRECISION :: USER_STREAMS(MAX_USER_STREAMS)
      DOUBLE PRECISION :: USER_SINES  (MAX_USER_STREAMS)

!  azimuth quadrature streams for BRDF

      DOUBLE PRECISION :: X_BRDF  ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: CX_BRDF ( MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SX_BRDF ( MAXSTREAMS_BRDF )

!  Output BRDF functions
!  =====================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC_CORE   ( 4,    MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_CORE ( 4,    MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_HELP   ( 4, 7, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_HELP ( 4, 7, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC_CORE   ( 4,    MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_CORE ( 4,    MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_HELP   ( 4, 7, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_HELP ( 4, 7, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_CORE ( 4,    MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_HELP ( 4, 7, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  local variables
!  ---------------

      LOGICAL          :: DO_COEFFS, Local_Isotropy
      INTEGER          :: I, UI, J, K, IB
      DOUBLE PRECISION :: PHI_W(MAXBEAMS), CPHI_W(MAXBEAMS), SPHI_W(MAXBEAMS)
      DOUBLE PRECISION :: SUNGLINT_COEFFS(7)

      INTEGER, PARAMETER :: LUM = 1
      INTEGER, PARAMETER :: LUA = 1

!   Wind-direction and coefficient set-up

      DO_COEFFS = .true.
      PHI_W = zero ; CPHI_W = one ; SPHI_W = zero
      Local_Isotropy = DO_FacetIsotropy 
      if ( .not.Local_Isotropy ) then
         DO IB = 1, nbeams
            PHI_W(IB)  = WINDDIR(IB)
            CPHI_W(IB) = cos(WINDDIR(IB) * deg_to_rad) 
            SPHI_W(IB) = sin(WINDDIR(IB) * deg_to_rad)
         ENDDO
      endif

!  Direct Bounce calculation
!  -------------------------
!mick fix 11/14/2019  - initialize DBKERNEL_BRDFUNC
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!  1/31/21. Version 2.8.3. Add doublet geometry option, in which the Exact-DB calculation
!   has the Azimuth angle coupled with the user zenith angle in a doublet.

      IF ( .NOT. DO_USER_OBSGEOMS ) THEN
         DO K = 1, N_USER_RELAZMS
            DO IB = 1, NBEAMS
               DO  UI = 1, N_USER_STREAMS
                  CALL GenGlint_GCM_CORE &
                    ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,       &
                      WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),    &
                      SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),  &
                      USER_SINES(UI), PHIANG(K), COSPHI(K), SINPHI(K), &
                      SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,UI,K,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,K,IB) )
               ENDDO
            ENDDO
         ENDDO
      ELSE IF ( DO_DOUBLET_GEOMETRY ) THEN
         DO IB = 1, NBEAMS
            DO  UI = 1, N_USER_STREAMS
               CALL GenGlint_GCM_CORE &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
                   WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),       &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),     &
                   USER_SINES(UI), PHIANG(UI), COSPHI(UI), SINPHI(UI), &
                   SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,UI,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,UI,LUA,IB) )
            ENDDO
         ENDDO
      ELSE
         DO IB = 1, NBEAMS
            CALL GenGlint_GCM_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,          &
               WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),       &
               SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),     &
               USER_SINES(IB), PHIANG(IB), COSPHI(IB), SINPHI(IB), &
               SUNGLINT_COEFFS, DBKERNEL_BRDFUNC_CORE(NNN,LUM,LUA,IB), DBKERNEL_BRDFUNC_HELP(NNN,1,LUM,LUA,IB) )
         ENDDO
      ENDIF

!      pause'after direct bounce'

!  Return if this is all you require

      IF ( DO_DBONLY ) RETURN

!  Incident Solar beam
!  ===================
!mick fix 11/14/2019  - initialize BRDFUNC_0, USER_BRDFUNC_0
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!  Quadrature outgoing directions

      DO IB = 1, NBEAMS
        DO I = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GenGlint_GCM_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),     &
               SZASURCOS(IB), SZASURSIN(IB), QUAD_STREAMS(I),    &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, BRDFUNC_0_CORE(NNN,I,IB,K), BRDFUNC_0_HELP(NNN,1,I,IB,K)  )
          ENDDO
        ENDDO
      ENDDO

!  User-streams outgoing directions
!   This is the "Truncated" Direct Bounce calculation
!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Introduce flag DO_FOURIER_USERSUN for controlling calculation of USER_BRDF_F_0.

      IF ( DO_USER_STREAMS .and. DO_FOURIER_USERSUN ) THEN

        IF (.NOT. DO_USER_OBSGEOMS ) THEN
          DO IB = 1, NBEAMS
            DO UI = 1, N_USER_STREAMS
              DO K = 1, NSTREAMS_BRDF
                CALL GenGlint_GCM_CORE &
                 ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                   WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),      &
                   SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(UI),    &
                   USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                   SUNGLINT_COEFFS, USER_BRDFUNC_0_CORE(NNN,UI,IB,K), USER_BRDFUNC_0_HELP(NNN,1,UI,IB,K) )
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO IB = 1, NBEAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GenGlint_GCM_CORE &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 WINDSPEED, PHI_W(IB), CPHI_W(IB), SPHI_W(IB),      &
                 SZASURCOS(IB), SZASURSIN(IB), USER_STREAMS(IB),    &
                 USER_SINES(IB), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, USER_BRDFUNC_0_CORE(NNN,LUM,IB,K), USER_BRDFUNC_0_HELP(NNN,1,LUM,IB,K) )
            ENDDO
          ENDDO
        ENDIF

      ENDIF

!  Incident quadrature directions (MULTIPLE SCATTERING)
!  ==============================
!mick fix 11/14/2019  - initialize BRDFUNC, USER_BRDFUNC
!                     - changed dimensional extent of elements being passed from subroutine
!                       output arrays to calling routine output arrays from NSSQ to NSM

!   Can only be treated with 1 Wind direction.....
!     if ( NBEAMS > 1) MUST assume local Facet Isotropy
!            --> set up Local Wind-direction and re-set coefficients flag.
!     if ( NBAMS = 1 ) use the first wind direction, no need to re-calculate coefficients

      if ( NBEAMS .gt. 1 ) then
         local_Isotropy = .true.
!         local_Isotropy = .false.  ! Bug 9/27/14, fixed
         PHI_W      = zero 
         CPHI_W     = one
         SPHI_W     = zero
         DO_COEFFS  = .true.
      endif
 
!  Outgoing quadrature directions

      DO I = 1, NSTREAMS
        DO J = 1, NSTREAMS
          DO K = 1, NSTREAMS_BRDF
            CALL GenGlint_GCM_CORE &
             ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,        &
               WINDSPEED, PHI_W(1), CPHI_W(1), SPHI_W(1),        &
               QUAD_STREAMS(J), QUAD_SINES(J), QUAD_STREAMS(I),  &
               QUAD_SINES(I), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
               SUNGLINT_COEFFS, BRDFUNC_CORE(NNN,I,J,K), BRDFUNC_HELP(NNN,1,I,J,K)  )
          ENDDO
        ENDDO
      ENDDO

!  User stream outgoing directions

      IF ( DO_USER_STREAMS ) THEN
        DO UI = 1, N_USER_STREAMS
          DO J = 1, NSTREAMS
            DO K = 1, NSTREAMS_BRDF
              CALL GenGlint_GCM_CORE &
               ( Local_Isotropy, DO_GlintShadow, DO_Coeffs,         &
                 WINDSPEED, PHI_W(1), CPHI_W(1), SPHI_W(1),        &
                 QUAD_STREAMS(J), QUAD_SINES(J), USER_STREAMS(UI),  &
                 USER_SINES(UI), X_BRDF(K), CX_BRDF(K), SX_BRDF(K), &
                 SUNGLINT_COEFFS, USER_BRDFUNC_CORE(NNN,UI,J,K), USER_BRDFUNC_HELP(NNN,1,UI,J,K) )
            ENDDO
          ENDDO
        ENDDO
      ENDIF

!  Finish

      RETURN
      END SUBROUTINE VBRDF_NewGCM_CORE_MAKER

!  End module

      END MODULE vbrdf_sup_core_routines_m

