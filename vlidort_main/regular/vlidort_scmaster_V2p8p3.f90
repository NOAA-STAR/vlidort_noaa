
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
! #            SPHERICITY_CORRECTION_MASTER (master, public)    #
! #            MS2pt_SPHERICITY_CORRECTION_V3                   #
! #            MS3pt_SPHERICITY_CORRECTION_V1                   #
! #            MSMpt_SPHERICITY_CORRECTION_V1                   #
! #                                                             #
! #            BOATOA_2point_conversion                         #
! #            BOATOA_3point_conversion                         #
! #            BOATOA_Mpoint_conversion                         #
! #                                                             #
! ###############################################################

MODULE VLIDORT_SPHERCORR_ROUTINES_m

!  HERE IS THE HISTORY
!  ===================

!  2-point correction history --
!    V1: 11/20/19. Originally Coded 20-29 November 2019 for VLIDORT Version 2.8.1
!    V2: 12/18/19. Extension to include BOA downwelling situation as well as TOA upwelling.
!    V3: 03/01/20. Renamed, along with addition of 3-point and multi-point Corrections

!  3-point correction history --
!    V1: 03/01/20. New 3pt Correction

!  Multi-point correction history --
!    V1: 03/01/20. New  Correction

!  1/31/21. Developed for Official VLIDORT Version 2.8.3 Package
!     R. Spurr. RT Solutions Inc.

!  11/28/22. Upgrade to Version 2.8.3.
!    -- Gather together sphericity correction code all in one module
!    -- Use single master routine to call required correction (new)

!  HERE IS THE VLIDORT IMPLEMENTATION
!  ==================================

!  Module files for VLIDORT.

      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m
      USE VLIDORT_MASTERS_m
 
!  9/8/22. CODECTOOL Upgrade to Version 2.8.5.
!    -- NEW master routine is public
!    -- Geometry conversion subroutines must be public (needed elsewhere)

public  ::  SPHERICITY_CORRECTION_MASTER, &
            BOATOA_2point_conversion, &
            BOATOA_3point_conversion, &
            BOATOA_Mpoint_conversion

!  All three MSType routines are private

private  :: MS2pt_SPHERICITY_CORRECTION_V3, &
            MS3pt_SPHERICITY_CORRECTION_V1, &
            MSMpt_SPHERICITY_CORRECTION_V1

contains

SUBROUTINE SPHERICITY_CORRECTION_MASTER ( do_debug_input,                   &
             VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup, MSType,             & ! VLIDORT type structure Inputs
             FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR,     & ! Component outputs (FO/MS)
             VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,                 & ! VLIDORT and Complete outputs
             Twilight_Flag, Fail, Main_Message, Local_Message, Local_Action ) ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  Sphericity type ('2' = 2-point, '3' = 3-point, 'M' = Multi-point)

      CHARACTER*1, intent(in) :: MSTYPE

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Main output with explanations

      DOUBLE PRECISION, INTENT(OUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(OUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(OUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector with 2-point corr
      DOUBLE PRECISION, INTENT(OUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector for BOA geo.     4 = 1 + 2
      DOUBLE PRECISION, INTENT(OUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector with MS2pt corr. 5 = 1 + 3

!  Twilight flag and exception handling

      LOGICAL      , intent(out) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(out) :: Main_message, Local_Message, Local_Action

!  Start Code
!  ==========

!  Initialize

      FO_STOKES_BOAGEOM   = zero
      MS_STOKES_BOAGEOM   = zero
      MS_STOKES_SPHERCORR = zero
      STOKES_BOAGEOM      = zero
      STOKES_SPHERCORR    = zero

!  Set exceptions

      FAIL          = .false.
      Twilight_Flag = .false.
      Main_message  = ' '
      Local_Message = ' '
      Local_Action  = ' '

!  Call to dedicated sphericity correction routine (either 2-point, 3-point or multi-point)

      if ( MSType .eq. '2' ) then
        CALL MS2pt_SPHERICITY_CORRECTION_V3 ( do_debug_input, &
             VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                     & ! VLIDORT type structure Inputs
             FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR,     & ! Component outputs (FO/MS)
             VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,                 & ! VLIDORT and Complete outputs
             Twilight_Flag, Fail, Local_Message, Local_Action )               ! Exceptions and errors

      else if ( MSType .eq. '3' ) then
        CALL MS3pt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
             VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                     & ! VLIDORT type structure Inputs
             FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR,     & ! Component outputs (FO/MS)
             VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,                 & ! VLIDORT and Complete outputs
             Twilight_Flag, Fail, Local_Message, Local_Action )               ! Exceptions and errors

      else if ( MSType .eq. 'M' ) then
        CALL MSMpt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
             VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                     & ! VLIDORT type structure Inputs
             FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR,     & ! Component outputs (FO/MS)
             VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,                 & ! VLIDORT and Complete outputs
             Twilight_Flag, Fail, Local_Message, Local_Action )               ! Exceptions and errors

      endif

!  Exception handling

      IF ( Fail ) THEN
        Main_message = 'Failure from Call to MS'//MStype//'pt_SPHERICITY_CORRECTION'
      endif

!  return

      return
END SUBROUTINE SPHERICITY_CORRECTION_MASTER

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    P R I V A T E    R O U T I N E S
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine MS2pt_SPHERICITY_CORRECTION_V3 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations (pre-initialized)

      DOUBLE PRECISION, INTENT(INOUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(INOUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector for BOA geo
      DOUBLE PRECISION, INTENT(INOUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector with 2-point corr
      DOUBLE PRECISION, INTENT(INOUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector for BOA geo.     4 = 1 + 2
      DOUBLE PRECISION, INTENT(INOUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector with MS2pt corr. 5 = 1 + 3

!  Twilight flag and exception handling (pre-initialized)

      LOGICAL      , intent(INOUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(INOUT) :: Local_Message, Local_Action

!  Local
!  -----

!  Variables for doing the sphericity calculation

      DOUBLE PRECISION :: MU, MU0, MU1, D01, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      DOUBLE PRECISION :: EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA(3), OBSGEOMS_TOA(3), COSSCAT
      DOUBLE PRECISION :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NLAYERS, NS, O1, DIRIDX

!  Start Code
!  ==========

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 2-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 2
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]

!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set second geometry through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Conversion. 12/18/19. Add Direction index.

      Call BOATOA_2point_conversion ( DIRIDX, &
         EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3), & ! input TOA/BOA heights, BOA Geometry
         cosscat, OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3) )                   ! output TOA geometry, scattering angle

!  Debug
!write(*,*)DIRIDX,OBSGEOMS_BOA(1), OBSGEOMS_BOA(2), OBSGEOMS_BOA(3)
!write(*,*)DIRIDX,OBSGEOMS_TOA(1), OBSGEOMS_TOA(2), OBSGEOMS_TOA(3)
!write(*,*)DIRIDX,ACOS(COSSCAT)/DEG_TO_RAD
!stop

!  Reset the VLIDORT second geometrical input (TOA)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_TOA(3)

!  re-set geometry numbers. 2 Geometries

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 2
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = 2
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 2
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = 2

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_TOA(1).gt.90.0_fpk ) then
         Twilight_flag = .true. ; return
      endif

!  debug - Check this against the FO geometry value.
!WRITE(*,*)VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3) 
!WRITE(*,*)VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) 

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling (simply done here)

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS2pt_SPHERICITY_CORRECTION_V3'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store Stokes-vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS) = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS) = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  Stokes-vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,2,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,2,1:NS) 
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,2,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,2,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,2,1:NS,DNIDX)
      endif

!write(*,*)FO_STOKES_BOAGEOM(1),FO_STOKES_TOAGEOM(1), VLIDORT_Out%Main%TS_STOKES(1,1:2,1,DNIDX) ; stop

!  perform 2-POINT SPHERICITY CORRECTION, Radiative Transfer
!  =========================================================

!  set the surface MS source term

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  get the Layer MS source terms by mid-point interpolation with Cos(SZA)
!    -- Works, except for the twilight and azimuth-flip conditions
!      Mu0 = COS(OBSGEOMS_BOA(1)*DEG_TO_RAD)
!      Mu1 = COS(OBSGEOMS_TOA(1)*DEG_TO_RAD)
!      D01 = ONE / ( Mu1 - Mu0 )
!      do N = 1, nlayers
!         Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N-1)) )
!         F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
!         DO O1 = 1, NS
!            AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,O1,N)
!write(33,*)N,AMSST(N,O1),VLIDORT_Out%Main%TS_PATHGEOMS(1,N)/DEG_TO_RAD
!         ENDDO
!      ENDDO

!  Alternative: get the Layer MS source terms by end-point interpolation with Cos(VZA)
!     -- Always works, however may not be accurate around azimuth-flip conditions

      Mu0 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
      Mu1 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
      D01 = ONE / ( Mu1 - Mu0 )
      do N = 1, nlayers
         Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
         F0 = ( Mu1 - Mu ) * D01 ; F1 = ONE - F0
         DO O1 = 1, NS
           AMSST(N,O1) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,O1,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,O1,N)
         ENDDO
!write(*,*)N,F0,VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1,N),VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1,N)
      ENDDO

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final Radiance answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MS2pt_SPHERICITY_CORRECTION_V3

!

subroutine MS3pt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations

      REAL(fpk), INTENT(INOUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(INOUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(INOUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector vector with 3-point corr
      REAL(fpk), INTENT(INOUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector vector for BOA geo.     4 = 1 + 2.
      REAL(fpk), INTENT(INOUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector vector with MS3pt corr. 5 = 1 + 3.

!  Twilight flag and exception handling

      LOGICAL      , intent(INOUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(INOUT) :: Local_Message, Local_Action

!  Local
!  -----

!  Variables for doing the sphericity calculation

      REAL(fpk) :: MU, MU0, MU1, MU2, F0, F1, TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      REAL(fpk) :: EARTH_RADIUS, TOA_HEIGHT, MID_HEIGHT, OBSGEOMS_BOA(3), OBSGEOMS_MID(3), OBSGEOMS_TOA(3), COSSCAT
      REAL(fpk) :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NLAYERS, NS, DIRIDX

!  Start Code
!  ==========

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 3-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set other geometries through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Geometrical Conversion. 3/1/20. Add MID_HEIGHT value

      Call BOATOA_3point_conversion ( DIRIDX, &
         EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_BOA,         & ! input TOA height, BOA Geometry
         MID_HEIGHT, cosscat, OBSGEOMS_TOA, OBSGEOMS_MID ) ! output TOA/MID geometries, scattering angle, mid-height

!  Reset the VLIDORT second (mid) and third (TOA) geometrical input

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(2,1:3) = OBSGEOMS_MID(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(2)                = OBSGEOMS_MID(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(2)     = OBSGEOMS_MID(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(2)            = OBSGEOMS_MID(3)

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(3,1:3) = OBSGEOMS_TOA(1:3)
      VLIDORT_ModIn%MSunrays%TS_SZANGLES(3)                = OBSGEOMS_TOA(1)
      VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(3)     = OBSGEOMS_TOA(2)
      VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(3)            = OBSGEOMS_TOA(3)

!  re-set geometry numbers. 3 Geometries

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = 3
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = 3
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = 3
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = 3

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_TOA(1).gt.90.0_fpk ) then
         Twilight_flag = .true. ; return
      endif

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MS3pt_SPHERICITY_CORRECTION_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store STOKES Vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS)  = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  STOKES Vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,3,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,3,1:NS) 
         STOKES_TOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,3,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS)  = VLIDORT_Sup%SS%TS_STOKES_SS(1,3,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)     = VLIDORT_Out%Main%TS_STOKES(1,3,1:NS,DNIDX)
      endif

!  perform 3-point SPHERICITY CORRECTION, Recurrence Relation
!  ==========================================================

!  set the surface MS source term, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS) = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  Alternative: get the Layer MS source terms by 3-point interpolation with Cos(VZA) using exponential relaxation
!     -- Always works, however may not be accurate around azimuth-flip conditions
!     -- THIS DOES NOT WORK GOOD AND WAS ABANDONED................................
!      Mu2 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
!      Mu1 = COS(OBSGEOMS_MID(2)*DEG_TO_RAD) ! Average value
!      Mu0 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
!     write(*,*)'mu2,mu1,mu0',mu2, mu1, mu0
!     write(397,'(I2,f10.6)')0,COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,0))
!      D21 = ONE / (Mu2 - Mu1)
!      do N = 1, nlayers
!         Mu = COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N))
!            ZTOP = ( VLIDORT_Out%Main%TS_LAYER_MSSTS(3,1:NS,N) - VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) )
!            ZBOT = ( VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) - VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N) )
!            Z = ZTOP/ZBOT ; KAY = Log(Z) * D21
!            BCON = ZBOT * EXP ( KAY * MU1 )  / ( Z - ONE )
!            ACON = VLIDORT_Out%Main%TS_LAYER_MSSTS(2,N) - BCON * EXP ( - KAY * MU1 )
!            AMSST(N,1:NS) = ACON + BCON * EXP ( - KAY * MU )
!        write(397,'(i2,f10.6,1p3e16.6,2x,1pe16.6)')&
!N,COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)),VLIDORT_Out%Main%TS_LAYER_MSSTS(1:3,1:NS,N), AMSST(N,1:NS)
!        ENDDO

!  Alternative: get the Layer MS source terms by 3-point Linear interpolation with Cos(VZA) 
!     -- Always works, however may not be accurate around azimuth-flip conditions
!     -- in any layer, Mu is the layer average cosine, and we interpolate to this value.

      Mu2 = COS(OBSGEOMS_BOA(2)*DEG_TO_RAD)
      Mu1 = COS(OBSGEOMS_MID(2)*DEG_TO_RAD) ! Average value of MU0 and MU2
      Mu0 = COS(OBSGEOMS_TOA(2)*DEG_TO_RAD)
      do N = 1, nlayers
        Mu = 0.5_fpk * ( COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)) + COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N-1)) )
        IF ( Mu.gt.mu1 ) then
          F0 = ( Mu1 - Mu ) / ( Mu1 - Mu0 ) ; F1 = ONE - F0
          AMSST(N,1:NS) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(3,1:NS,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N)
        ELSE
          F0 = ( Mu2 - Mu ) / ( Mu2 - Mu1 ) ; F1 = ONE - F0
          AMSST(N,1:NS) = F0 * VLIDORT_Out%Main%TS_LAYER_MSSTS(2,1:NS,N) + F1 * VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N)
        ENDIF
      enddo

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final Radiance answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MS3pt_SPHERICITY_CORRECTION_V1

!

subroutine MSMpt_SPHERICITY_CORRECTION_V1 ( do_debug_input, &
           VLIDORT_FixIn, VLIDORT_ModIn, VLIDORT_Sup,                 & ! VLIDORT Inputs
           FO_STOKES_BOAGEOM, MS_STOKES_BOAGEOM, MS_STOKES_SPHERCORR, & ! Component outputs (FO/MS)
           VLIDORT_Out, STOKES_BOAGEOM, STOKES_SPHERCORR,             & ! VLIDORT and Complete  outputs
           Twilight_Flag, Fail, Local_Message, Local_Action )           ! Exceptions and errors

!  Implicit none

      IMPLICIT NONE

!  Debug input flag

      LOGICAL, intent(in) :: do_debug_input

!  VLIDORT structures (including supplements i/o)

      TYPE(VLIDORT_Fixed_Inputs)   , Intent(in)    :: VLIDORT_FixIn
      TYPE(VLIDORT_Modified_Inputs), Intent(inout) :: VLIDORT_ModIn
      TYPE(VLIDORT_Sup_InOut)      , Intent(inout) :: VLIDORT_Sup
      TYPE(VLIDORT_Outputs)        , Intent(inout) :: VLIDORT_Out

!  Output
!  ------

!  Main output with explanations

      REAL(fpk), INTENT(INOUT) :: FO_STOKES_BOAGEOM(4)   ! 1. First-order (SS/DB) STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(INOUT) :: MS_STOKES_BOAGEOM(4)   ! 2. MS STOKES Vector vector for BOA geo
      REAL(fpk), INTENT(INOUT) :: MS_STOKES_SPHERCORR(4) ! 3. MS STOKES Vector vector with 3-point corr
      REAL(fpk), INTENT(INOUT) :: STOKES_BOAGEOM(4)      ! 4. Full STOKES Vector vector for BOA geo.          4 = 1 + 2.
      REAL(fpk), INTENT(INOUT) :: STOKES_SPHERCORR(4)    ! 5. Full STOKES Vector vector with Multipoint corr. 5 = 1 + 3.

!  Twilight flag and exception handling

      LOGICAL      , intent(INOUT) :: Twilight_Flag, FAIL
      CHARACTER*(*), intent(INOUT) :: Local_Message, Local_Action

!  Local
!  -----

      INTEGER, parameter :: maxmults = MAXLAYERS + 1

!  Variables for doing the sphericity calculation

      REAL(fpk) :: TRANS, SMSST(4), AMSST(MAXLAYERS,4)
      REAL(fpk) :: EARTH_RADIUS, TOA_HEIGHT, OBSGEOMS_ALL(maxmults,3),  OBSGEOMS_BOA(3), COSSCAT
      REAL(fpk) :: FO_STOKES_TOAGEOM(4), STOKES_TOAGEOM(4)

!  Miscellaneous variables
!   -- DirIdx chooses TOA Upwelling (1) or BOA Downwelling (2)

      LOGICAL :: Verbose
      INTEGER :: N, NMULT, NG, NG1, NS, NLAYERS, DIRIDX

!  Start Code
!  ==========

!  Set Proxies (saved values)

      ns      = VLIDORT_FixIn%Cont%TS_nstokes
      nlayers = VLIDORT_FixIn%Cont%TS_nlayers

!  Direction index

      DirIdx = 1 ; if ( VLIDORT_FixIn%Bool%TS_DO_DNWELLING ) DirIdx = 2

!  Set debug output flag

      Verbose   = .false.
!      Verbose   = .true.

!  Set the DO_MSSTS flag for the 3-point sphericity correction
!  There are very strict conditions for this specialist option, as follows :==>
!     1. Either DO_UPWELLING or DO_DNWELLING must be set, Not Both !!!!
!     2. DO_FULLRAD_MODE must be set
!     3. DO_OBSERVATION_GEOMETRY must be set, with N_USER_OBSGEOMS = 3
!     4. DO_FOCORR and DO_FOCORR_OUTGOING must both be set
!     5a. Upwelling  : N_USER_LEVELS = 1, and USER_LEVELS(1) = 0.0            [ TOA output only ]
!     5b. Downwelling: N_USER_LEVELS = 1, and USER_LEVELS(1) = Real(nlayers)  [ BOA output only ]
!  These checks have been implemented inside VLIDORT.
!      VLIDORT_FixIn%Bool%TS_DO_MSSTS = .true.

!  Set other geometries through BOA-to-TOA conversion. [Overwrites config-file input]
!  --------------------------------------------------

!  inputs to conversion routine (BOA geometry, height grid, earth radius)

      EARTH_RADIUS      = VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS   
      TOA_HEIGHT        = VLIDORT_FixIn%Chapman%TS_height_grid(0)
      OBSGEOMS_BOA(1:3) = VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1,1:3)

!  Geometrical Conversion. 3/1/20. Multipoint calculation

      NMULT = NLAYERS + 1
      CALL BOATOA_Mpoint_conversion ( &
          maxmults, maxlayers, nlayers, nmult, DirIdx,                     & ! Input numbers
          EARTH_RADIUS, VLIDORT_FixIn%Chapman%TS_height_grid, OBSGEOMS_BOA, & ! input Heights, BOA Geometry
          cosscat, OBSGEOMS_ALL )                                            ! output all geometries, Cosine scattering angle

!  Multi-points: Reset the VLIDORT geometrical input at all layer boundaries from BOA to TOA
!    NMULT is TOA geometry

      DO N = 2, NMULT
        VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(N,1:3) = OBSGEOMS_ALL(N,1:3)
        VLIDORT_ModIn%MSunrays%TS_SZANGLES(N)                = OBSGEOMS_ALL(N,1)
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(N)     = OBSGEOMS_ALL(N,2)
        VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(N)            = OBSGEOMS_ALL(N,3)
      ENDDO

!  re-set geometry numbers. NMULT Geometries in all

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = NMULT
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES      = NMULT
      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = NMULT
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = NMULT

!  Twilight condition at TOA (SZA > 90); set dummy output flag and skip calculation. 

      if ( OBSGEOMS_ALL(NMULT,1).gt.90.0_fpk ) then
         Twilight_flag = .true. ; return
      endif

!  Call to VLIDORT

      CALL VLIDORT_MASTER ( do_debug_input, &
          VLIDORT_FixIn, &
          VLIDORT_ModIn, &
          VLIDORT_Sup,   &
          VLIDORT_Out )

!  Exception handling

      if ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK  .eq. VLIDORT_SERIOUS .or. &
           VLIDORT_Out%Status%TS_STATUS_CALCULATION .eq. VLIDORT_SERIOUS ) then
         Local_message = 'Some errors arising from VLIDORT CALL in MSMpt_SPHERICITY_CORRECTION_V1'
         Local_Action  = 'In Driver, use call to VLIDORT_WRITE_STATUS to examine errors'
         FAIL = .true. ; return
      Endif

!  Store STOKES Vector results for BOA-Geometry. THIS IS PART OF THE OUTPUT from this subroutine
!    ==> Choose between upwelling (DirIDx = 1) and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_BOAGEOM(1:NS)    = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,1,1:NS) 
         MS_STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)       = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,UPIDX)
      else
         FO_STOKES_BOAGEOM(1:NS)    = VLIDORT_Sup%SS%TS_STOKES_SS(1,1,1:NS,DNIDX)
         MS_STOKES_BOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)  - FO_STOKES_BOAGEOM(1:NS)
         STOKES_BOAGEOM(1:NS)       = VLIDORT_Out%Main%TS_STOKES(1,1,1:NS,DNIDX)
      endif

!  DEBUG: Check geometry
!      do n = 0, nlayers
!         NG = NMULT - n
!  write(*,*)VLIDORT_Out%Main%TS_PATHGEOMS(1,N)/DEG_TO_RAD,VLIDORT_Out%Main%TS_PATHGEOMS(2,N)/DEG_TO_RAD,OBSGEOMS_ALL(NG,1:2)
!      enddo
!stop

!  STOKES Vector results for TOA-geometry (only needed for debug purposes)
!    ==> Choose between upwelling and Downwelling scenarios.

      if ( DirIdx .eq. 1 ) then
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,NMULT,1:NS,UPIDX) + VLIDORT_Sup%SS%TS_STOKES_DB(1,NMULT,1:NS) 
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,NMULT,1:NS,UPIDX)
      else
         FO_STOKES_TOAGEOM(1:NS) = VLIDORT_Sup%SS%TS_STOKES_SS(1,NMULT,1:NS,DNIDX)
         STOKES_TOAGEOM(1:NS)    = VLIDORT_Out%Main%TS_STOKES(1,NMULT,1:NS,DNIDX)
      endif

!  perform multi-point SPHERICITY CORRECTION, Recurrence Relation
!  ==============================================================

!  set the surface MS source term, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
         SMSST(1:NS)  = VLIDORT_Out%Main%TS_SURF_MSSTS(1,1:NS)
      endif

!  Alternative: get the Layer MS source terms by layer averaging the multi-point values
!     -- Always works, however may not be accurate around azimuth-flip conditions

      do N = 1, nlayers
         NG = NMULT - N + 1 ; NG1 = NG - 1
         AMSST(N,1:NS) = HALF * ( VLIDORT_Out%Main%TS_LAYER_MSSTS(NG,1:NS,N) + VLIDORT_Out%Main%TS_LAYER_MSSTS(NG1,1:NS,N))
!            write(398,'(i2,2f10.6,2x,1p2e16.6,2x,1pe16.6)')N,&
!               COS(VLIDORT_Out%Main%TS_PATHGEOMS(1,N)),COS(VLIDORT_Out%Main%TS_PATHGEOMS(2,N)),&
!               VLIDORT_Out%Main%TS_LAYER_MSSTS(NMULT,1:NS,N), VLIDORT_Out%Main%TS_LAYER_MSSTS(1,1:NS,N), AMSST(N,1:NS)
      ENDDO

!  Verbose debug output

      if ( verbose ) then
        do NG = 1, NMULT
          write(399,'(i2,2f8.5,1x,1p23e11.4)')NG,&
               COS(OBSGEOMS_ALL(NG,1)*DEG_TO_RAD),COS(OBSGEOMS_ALL(NG,2)*DEG_TO_RAD),&
               VLIDORT_Out%Main%TS_LAYER_MSSTS(NG,1:NS,1:NLAYERS)
        ENDDO
        stop
      endif

!  Radiative Transfer recursion
!  ============================

!  Perform the radiative transfer recursion for MS_STOKES_SPHERCORR, the spherically-corrected MS field
!   ==> Start with the surface term SMSST, Upwelling scenario

      if ( DirIdx .eq. 1 ) then
!mick fix 1/5/2021 - trim dim passed from SMSST
         MS_STOKES_SPHERCORR(1:NS) = SMSST(1:NS)
         DO N = NLAYERS, 1, -1
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      else
         MS_STOKES_SPHERCORR(1:NS) = ZERO
         DO N = 1, NLAYERS
            TRANS = VLIDORT_Out%Main%TS_LOSTRANS(1,N)
            MS_STOKES_SPHERCORR(1:NS) = TRANS * MS_STOKES_SPHERCORR(1:NS) + AMSST(N,1:NS) 
         ENDDO
      endif

!  Final STOKES Vector answer: Add FO Outgoing result to MS field.

      STOKES_SPHERCORR(1:NS) = MS_STOKES_SPHERCORR(1:NS) + FO_STOKES_BOAGEOM(1:NS)

!  Verbose debug

      if ( verbose ) then
         write(242,*) OBSGEOMS_BOA(1),OBSGEOMS_BOA(2)
         write(242,*) 'STOKES SS/DB BOA, MS BOA, IBOA : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_BOAGEOM(1:NS),  STOKES_BOAGEOM(1:NS)
         write(242,*) 'STOKES SS/DB BOA, MS Sph, ISph : ',&
              FO_STOKES_BOAGEOM(1:NS),MS_STOKES_SPHERCORR(1:NS),STOKES_SPHERCORR(1:NS)
         write(242,*) 'STOKES SS/DB TOA, MS TOA, ITOA : ',&
              FO_STOKES_TOAGEOM(1:NS),STOKES_TOAGEOM(1:NS)-FO_STOKES_TOAGEOM(1:NS), STOKES_TOAGEOM(1:NS)
      endif

!  normal Finish

      return
end subroutine MSMpt_SPHERICITY_CORRECTION_V1

subroutine BOATOA_2point_conversion ( DirIdx, &
     Rearth, Hatmos, theta_boa, alpha_boa, phi_boa, & ! input TOA/BOA heights, BOA Geometry
     cosscat, theta_toa, alpha_toa, phi_toa )         ! output TOA geometry, scattering angle

!  Completely stand-alone 2-point conversion routine.
!   12/18/19. Add direction index for downwelling case

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx
   real(fpk), intent(in)  :: Rearth, Hatmos, theta_boa, alpha_boa, phi_boa
   real(fpk), intent(out) :: theta_toa, alpha_toa, phi_toa, cosscat

!  Local

   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3)
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha_toa, calpha_toa, stheta_toa, ctheta_toa, cphi_toa
   real(fpk) :: term1, term2, rtoa, cumangle
   real(fpk) :: alpha_toa_R, theta_toa_R, phi_toa_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  TOA calculation
!  ===============

!  Vza at TOA (sine rule)

   rtoa = rearth+hatmos
   salpha_toa = rearth * salpha_boa / rtoa
   calpha_toa = sqrt(one-salpha_toa*salpha_toa)
   alpha_toa_R = asin(salpha_toa)
   cumangle = alpha_boa_R - alpha_toa_R

!  Calculate SZA at TOA, using vector geometry

   px(1) = - rtoa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rtoa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_toa  = - b/rtoa
   stheta_toa  = sqrt(one-ctheta_toa*ctheta_toa)
   theta_toa_R = acos(ctheta_toa)

!  Calculate PHI at TOA using constancy of scattering angle

   cphi_toa = (cosscat+vsign*calpha_toa*ctheta_toa)/stheta_toa/salpha_toa
   if ( cphi_toa.gt.one)  cphi_toa = one
   if ( cphi_toa.lt.-one) cphi_toa = -one
   phi_toa_R = acos(cphi_toa)
   if ( phi_boa.gt.180.0_fpk) phi_toa_R = 2.0_fpk * Pie - phi_toa_R

!  Final answers for TOA geometry (in degrees)

   alpha_toa = alpha_toa_R / dtr
   theta_toa = theta_toa_R / dtr
   phi_toa   = phi_toa_R   / dtr

!  finish

   return
end subroutine BOATOA_2point_conversion

!

subroutine BOATOA_3point_conversion ( DirIdx, &
     Rearth, Hatmos, BOA_geometry,                & ! input TOA height, BOA Geometry
     Hmid, cosscat, TOA_geometry, MID_geometry )    ! output TOA/MID geometries, scattering angle, Midheight

!  3/1/20. Completely new, stand-alone 3-point geometry conversion routine.

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx
   real(fpk), intent(in)  :: Rearth, Hatmos, BOA_geometry(3)
   real(fpk), intent(out) :: TOA_geometry(3), MID_geometry(3), cosscat, Hmid

!  Local

   real(fpk) :: alpha_boa, theta_boa, phi_boa
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3), term1, term2, cumangle
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha_toa, calpha_toa, stheta_toa, ctheta_toa, cphi_toa
   real(fpk) :: rtoa, alpha_toa_R, theta_toa_R, phi_toa_R
   real(fpk) :: salpha_mid, calpha_mid, stheta_mid, ctheta_mid, cphi_mid
   real(fpk) :: rmid, alpha_mid_R, theta_mid_R, phi_mid_R

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   theta_boa = BOA_GEOMETRY(1)
   alpha_boa = BOA_GEOMETRY(2)
   phi_boa   = BOA_GEOMETRY(3)

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  TOA calculation
!  ===============

!  Vza at TOA (sine rule)

   rtoa = rearth+hatmos
   salpha_toa = rearth * salpha_boa / rtoa
   calpha_toa = sqrt(one-salpha_toa*salpha_toa)
   alpha_toa_R = asin(salpha_toa)
   cumangle = alpha_boa_R - alpha_toa_R

!  Calculate SZA at TOA, using vector geometry

   px(1) = - rtoa * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rtoa * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_toa  = - b/rtoa
   stheta_toa  = sqrt(one-ctheta_toa*ctheta_toa)
   theta_toa_R = acos(ctheta_toa)

!  Calculate PHI at TOA using constancy of scattering angle

   cphi_toa = (cosscat+vsign*calpha_toa*ctheta_toa)/stheta_toa/salpha_toa
   if ( cphi_toa.gt.one)  cphi_toa = one
   if ( cphi_toa.lt.-one) cphi_toa = -one
   phi_toa_R = acos(cphi_toa)
   if ( phi_boa.gt.180.0_fpk) phi_toa_R = 2.0_fpk * Pie - phi_toa_R

!  Final answers for TOA geometry (in degrees)

   TOA_GEOMETRY(1) = theta_toa_R / dtr
   TOA_GEOMETRY(2) = alpha_toa_R / dtr
   TOA_GEOMETRY(3) = phi_toa_R   / dtr

!  MID_HEIGHT calculation
!  ======================

!  determine mid height by averaging cos(VZA) of the TOA and BOA values.

   calpha_mid = ( calpha_toa + calpha_boa ) * 0.5_fpk
   salpha_mid = sqrt ( 1.0_fpk - calpha_mid * calpha_mid )

!  Use sine rule to get the mid-height value

   rmid  = rearth * salpha_boa / salpha_mid
   hmid  = rmid - rearth

!  VZA and earth-centered angles at Mid-height

   alpha_mid_R = asin(salpha_mid)
   cumangle = alpha_boa_R - alpha_mid_R

!  Calculate SZA at mid-height, using via vector geometry

   px(1) = - rmid * sin(CumAngle)
   px(2) = 0.0_fpk
   px(3) =   rmid * cos(CumAngle)
   b = DOT_PRODUCT(px,SolarDirection)
   ctheta_mid  = - b/rmid
   stheta_mid  = sqrt(one-ctheta_mid*ctheta_mid)
   theta_mid_R = acos(ctheta_mid)

!  Calculate PHI at mid-height using constancy of scattering angle

   cphi_mid = (cosscat+vsign*calpha_mid*ctheta_mid)/stheta_mid/salpha_mid
   if ( cphi_mid.gt.one)  cphi_mid = one
   if ( cphi_mid.lt.-one) cphi_mid = -one
   phi_mid_R = acos(cphi_mid)
   if ( phi_boa.gt.180.0_fpk) phi_mid_R = 2.0_fpk * Pie - phi_mid_R

!  Final answers for mid-height geometry (in degrees)

   MID_GEOMETRY(1) = theta_mid_R / dtr
   MID_GEOMETRY(2) = alpha_mid_R / dtr
   MID_GEOMETRY(3) = phi_mid_R   / dtr

!  finish

   return
end subroutine BOATOA_3point_conversion

!

subroutine BOATOA_Mpoint_conversion ( &
     maxgeoms, maxlayers, nlayers, ngeoms, DirIdx, & ! Input numbers
     Rearth, Heights, BOA_geometry,                & ! input Heights, BOA Geometry
     cosscat, ALL_geometry )                         ! output all geometries, Cosine scattering angle

!  3/2/20. Completely new, stand-alone Multi-point geometry conversion routine.

   implicit none

   integer, parameter :: fpk = selected_real_kind(15)

!  Inputs/Outputs

   integer  , intent(in)  :: DirIdx, maxgeoms, maxlayers, nlayers, ngeoms
   real(fpk), intent(in)  :: Rearth, Heights(0:maxlayers), BOA_geometry(3)
   real(fpk), intent(out) :: ALL_geometry(maxgeoms,3), cosscat

!  Local

   integer   :: n, ng
   real(fpk) :: alpha_boa, theta_boa, phi_boa
   real(fpk) :: alpha_boa_R, theta_boa_R, phi_boa_R, dtr, pie, zero, one
   real(fpk) :: SolarDirection(3), vsign, vsign_solar, b, px(3), term1, term2, cumangle
   real(fpk) :: salpha_boa, calpha_boa, stheta_boa, ctheta_boa, cphi_boa, sphi_boa
   real(fpk) :: salpha, calpha, stheta, ctheta, cphi
   real(fpk) :: rad, alpha_R, theta_R, phi_R

!   real(fpk) :: stheta_mid, ctheta_mid, cphi_mid
!   real(fpk) :: theta_mid_R, phi_mid_R
!   real(fpk) :: alpha_mid_R, calpha_mid, salpha_mid, rmid

!  constants

   zero = 0.0_fpk ; one = 1.0_fpk
   Pie = acos(-1.0_fpk)
   dtr = Pie / 180.0_fpk

!  BOA angles in radians (_R) and cos/sin

   All_Geometry(1,1:3) = BOA_Geometry(1:3)
   theta_boa = BOA_GEOMETRY(1)
   alpha_boa = BOA_GEOMETRY(2)
   phi_boa   = BOA_GEOMETRY(3)

   alpha_boa_R = alpha_boa * dtr
   salpha_boa  = sin(alpha_boa_R)
   calpha_boa  = cos(alpha_boa_R)

   theta_boa_R    = theta_boa * dtr
   stheta_boa     = sin(theta_boa_R)
   ctheta_boa     = cos(theta_boa_R)
      
   phi_boa_R   = phi_boa * dtr
   cphi_boa    = cos(phi_boa_R)
   sphi_boa    = sin(phi_boa_R)

!  Solar direction (unit vector)

   vsign       = 1.0_fpk ; if ( DirIdx .eq. 2 ) vsign = -1.0_fpk
   vsign_solar = vsign ! vsign_solar = 1.0_fpk
   SolarDirection(1) = - stheta_boa * cphi_boa * vsign_solar
   SolarDirection(2) = - stheta_boa * sphi_boa
   SolarDirection(3) = - ctheta_boa

!  scattering angle

   term1 = salpha_boa * stheta_boa * cphi_boa
   term2 = calpha_boa * ctheta_boa
   cosscat = -vsign * term2 + term1 

!  calculation
!  ===========

   do ng = 2, ngeoms
     n = ngeoms - ng

!  Vza at boundary-level height (sine rule)

     rad = rearth+heights(n)
     salpha = rearth * salpha_boa / rad
     calpha = sqrt(one-salpha*salpha)
     alpha_R = asin(salpha)
     cumangle = alpha_boa_R - alpha_R

!  Calculate SZA at this boundary, using vector geometry

     px(1) = - rad * sin(CumAngle)
     px(2) = 0.0_fpk
     px(3) =   rad * cos(CumAngle)
     b = DOT_PRODUCT(px,SolarDirection)
     ctheta  = - b/rad
     stheta  = sqrt(one-ctheta*ctheta)
     theta_R = acos(ctheta)

!  Calculate PHI using constancy of scattering angle

     cphi = (cosscat+vsign*calpha*ctheta)/stheta/salpha
     if ( cphi.gt.one)  cphi = one
     if ( cphi.lt.-one) cphi = -one
     phi_R = acos(cphi)
     if ( phi_boa.gt.180.0_fpk) phi_R = 2.0_fpk * Pie - phi_R

!  Final answers for TOA geometry (in degrees)

     ALL_GEOMETRY(ng,1) = theta_R / dtr
     ALL_GEOMETRY(ng,2) = alpha_R / dtr
     ALL_GEOMETRY(ng,3) = phi_R   / dtr

   enddo

!  finish

   return
end subroutine BOATOA_Mpoint_conversion

!  end module

END MODULE VLIDORT_SPHERCORR_ROUTINES_m
