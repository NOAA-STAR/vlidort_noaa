! ###############################################################
! #                                                             #
! #                    THE VLIDORT MODEL                        #
! #                                                             #
! #  Vectorized LInearized Discrete Ordinate Radiative Transfer #
! #  -          --         -        -        -         -        #
! #                                                             #
! ###############################################################

! ###############################################################
! #                                                             #
! #  Author :      Robert. J. D. Spurr                          #
! #                                                             #
! #  Address :     RT Solutions, inc.                           #
! #                9 Channing Street                            #
! #                Cambridge, MA 02138, USA                     #
! #                Tel: (617) 492 1183                          #
! #                                                             #
! #  Email :       rtsolutions@verizon.net                      #
! #                                                             #
! #  Versions     :   2.0, 2.2, 2.3, 2.4, 2.4R, 2.4RT, 2.4RTC,  #
! #                   2.5, 2.6, 2.7, 2.8                        #
! #  Release Date :   December 2005  (2.0)                      #
! #  Release Date :   March 2007     (2.2)                      #
! #  Release Date :   October 2007   (2.3)                      #
! #  Release Date :   December 2008  (2.4)                      #
! #  Release Date :   April 2009     (2.4R)                     #
! #  Release Date :   July 2009      (2.4RT)                    #
! #  Release Date :   October 2010   (2.4RTC)                   #
! #  Release Date :   March 2011     (2.5)                      #
! #  Release Date :   May 2012       (2.6)                      #
! #  Release Date :   August 2014    (2.7)                      #
! #  Release Date :   May 2017       (2.8)                      #
! #                                                             #
! #       NEW: TOTAL COLUMN JACOBIANS         (2.4)             #
! #       NEW: BPDF Land-surface KERNELS      (2.4R)            #
! #       NEW: Thermal Emission Treatment     (2.4RT)           #
! #       Consolidated BRDF treatment         (2.4RTC)          #
! #       f77/f90 Release                     (2.5)             #
! #       External SS / New I/O Structures    (2.6)             #
! #                                                             #
! #       SURFACE-LEAVING / BRDF-SCALING      (2.7)             #
! #       TAYLOR Series / OMP THREADSAFE      (2.7)             #
! #       New Water-Leaving Treatment         (2.8)             #
! #       LBBF & BRDF-Telescoping, enabled    (2.8)             #
! #       Several Performance Enhancements    (2.8)             #
! #                                                             #
! ###############################################################

!    ###########################################################
!    #                                                         #
!    # This is Version 2.8 of the VLIDORT software library.    #
!    # This library comes with the GNU General Public License, #
!    # Version 3.0. Please read this license carefully.        #
!    #                                                         #
!    #      Copyright (c) 2003-2017.                           #
!    #          Robert Spurr, RT Solutions Inc.                #
!    #                                                         #
!    # This file is part of VLIDORT Version 2.8.               #
!    #                                                         #
!    # VLIDORT is free software: you can redistribute it       #
!    # and/or modify it under the terms of the GNU General     #
!    # Public License as published by the Free Software        #
!    # Foundation, either version 3 of the License, or any     #
!    # later version.                                          #
!    #                                                         #
!    # VLIDORT is distributed in the hope that it will be      #
!    # useful, but WITHOUT ANY WARRANTY; without even the      #
!    # implied warranty of MERCHANTABILITY or FITNESS FOR A    #
!    # PARTICULAR PURPOSE.  See the GNU General Public License #
!    # for more details.                                       #
!    #                                                         #
!    # You should have received a copy of the GNU General      #
!    # Public License along with VLIDORT Version 2.8.          #
!    # If not, see <http://www.gnu.org/licenses/>.             #
!    #                                                         #
!    ###########################################################

! ###############################################################
! #                                                             #
! # Subroutines in this File                                    #
! #                                                             #
! #            Testpoint output subroutines starting            #
! #            with the prefixes TP*                            #
! #                                                             #
! ###############################################################

!Note that this file is NOT a module (and should never be one)

!******************************************************************************
!******************************************************************************
!                 "Test point prefix / VLIDORT module" Reference
!                 ---------------------------------------------

!Test point naming convention:  below, "*" is recommended to be either nothing or
!a capital letter (e.g. TP1, TP1A, TP1B, ... , TP11, TP11A, TP11B, etc...)

!Test point prefix  VLIDORT module
!-----------------  --------------
!TP1*               vlidort_aux.f90
!TP2*               vlidort_bvproblem.f90
!TP3*               vlidort_corrections.f90        (removed 2.8)
!TP4*               vlidort_geometry.f90
!TP5*               vlidort_getplanck.f90
!TP6*               vlidort_inputs.f90
!TP7*               vlidort_intensity.f90
!TP8*               vlidort_la_miscsetups.f90
!TP9*               vlidort_lbbf_jacobians_vector.f90
!TP10*              vlidort_lc_bvproblem.f90
!TP11*              vlidort_lc_corrections.f90     (removed 2.8)
!TP12*              vlidort_lc_miscsetups.f90
!TP13*              vlidort_lc_PostProcessing.f90
!TP14*              vlidort_lcs_masters_V2p8p3.f90
!TP15*              vlidort_lc_wfatmos.f90
!TP16*              vlidort_l_inputs.f90
!TP17*              vlidort_lp_bvproblem.f90
!TP18*              vlidort_lp_corrections.f90     (removed 2.8)
!TP19*              vlidort_lpc_solutions.f90
!TP20*              vlidort_lp_miscsetups.f90
!TP21*              vlidort_lp_PostProcessing.f90
!TP22*              vlidort_lps_masters_V2p8p3.f90
!TP23*              vlidort_lp_wfatmos.f90
!TP24*              vlidort_ls_corrections.f90     (removed 2.8)
!TP25*              vlidort_ls_wfsleave.f90
!TP26*              vlidort_ls_wfsurface.f90
!TP27*              vlidort_l_thermalsup.f90
!TP28*                     not used
!TP29*              vlidort_masters_V2p8p3.f90
!TP30*              vlidort_miscsetups.f90
!TP31*              vlidort_PostProcessing.f90
!TP32*              vlidort_solutions.f90
!TP33*              vlidort_Taylor.f90
!TP34*              vlidort_thermalsup.f90

!TP35*              vlidort_vfo_interface.f90      (added 2.8)
!TP36*              vlidort_vfo_lcs_interface.f90  (added 2.8)
!TP37*              vlidort_vfo_lps_interface.f90  (added 2.8)
!TP38*              vlidort_transflux_MK3.f90      (added 2.8 - removed 2.??)

! (TP39*-TP69* reserved for possible later modules common to LIDORT & VLIDORT) 

!TP71*              vlidort_multipliers.f90        (vlidort only - removed 2.??)
!TP72*              vlidort_la_corrections.f90     (vlidort only - removed 2.8)
!TP73*              vlidort_lc_solutions.f90       (vlidort only)
!TP74*              vlidort_lp_solutions.f90       (vlidort only)
!TP75*              vlidort_lpc_bvproblem.f90      (vlidort only)

!TP76*              vlidort_converge.f90           (vlidort only)
!TP77*              vlidort_mediaprops.f90         (vlidort only)
!TP78*              vlidort_lc_mediaprops.f90      (vlidort only)
!TP79*              vlidort_lcs_converge.f90       (vlidort only)
!TP80*              vlidort_lp_mediaprops.f90      (vlidort only)
!TP81*              vlidort_lps_converge.f90       (vlidort only)

!TP82*              vlidort_scmaster_V2p8p3.f90     (vlidort only - added 2.8.3)
!TP83*              vlidort_lcs_scmaster_V2p8p3.f90 (vlidort only - added 2.8.3)
!TP84*              vlidort_lps_scmaster_V2p8p3.f90 (vlidort only - added 2.8.3)

!None               vlidort_writemodules.f90
!None               vlidort_l_writemodules.f90

!None               vlidort_pack.f90
!None               vlidort_unpack.f90
!None               vlidort_l_pack.f90
!None               vlidort_l_unpack.f90
!None               vlidort_lc_pack.f90
!None               vlidort_lc_unpack.f90
!None               vlidort_lp_pack.f90
!None               vlidort_lp_unpack.f90

!******************************************************************************
!******************************************************************************
!Test point subroutine template
!------------------------------

!It is recommended that:
!(1) the template be copied to the appropriate VLIDORT module section below and
!    then the "#" be replaced with EITHER the assigned module number OR the
!    assigned module number & a letter for the current testpoint subroutine to be built
! or
!(2) copy a testpoint that has ALREADY been built, modify its subroutine label, and
!    grab items (looping, etc ...) from this template as needed


!SUBROUTINE TP# ()
!
!USE VLIDORT_pars_m
!IMPLICIT NONE
!
!INTEGER, intent(in)   :: 
!DOUBLE PRECISION, intent(in) :: 
!
!INTEGER :: I,J,K
!
!write(*,*)
!write(*,*) 'TP#: IN SUBROUTINE ??'
!
!write(*,*)
!write(*,*) ' = ',,'  = ',,'  = ',,'  = ',
!write(*,*) ' = ',
!write(*,*) ' = ',
!write(*,*) ' = ',
!write(*,*) ' = ',
!write(*,*) ' = ',
!write(*,*) ' = ',
!
!write(*,*)
!DO I=1,
!  write(*,*) 'I = ',I,' (I) = ',(I)
!ENDDO
!
!write(*,*)
!DO J=1,
!  DO I=1,
!    write(*,*) 'I = ',I,' J = ',J,' (I,J) = ',(I,J)
!  ENDDO
!ENDDO
!
!write(*,*)
!DO K=1,
!  DO J=1,
!    DO I=1,
!      write(*,*) 'I = ',I,' J = ',J,' K = ',K,' (I,J,K) = ',(I,J,K)
!    ENDDO
!  ENDDO
!ENDDO
!
!END SUBROUTINE TP#


!******************************************************************************
!******************************************************************************

INTEGER FUNCTION ATPI()

! Define the atmospheric test point index (ATPI) - the index of the atmospheric
!   wvl or wvn whose RT solution you want the test point to check

  !ATPI = -1
  !ATPI = 0
  !ATPI = 1
  ATPI = 7

END FUNCTION ATPI

!------------------------------------------------------------------------------
!                   Begin test point subroutine list
!------------------------------------------------------------------------------

!******************************************************************************
!******************************************************************************
!                   Section: TP1*  --> vlidort_aux.f90         
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP2*  --> vlidort_bvproblem.f90   
!******************************************************************************
!******************************************************************************

SUBROUTINE TP2A (AWI,NSTREAMS,NSTOKES,R2_BEAM)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,NSTREAMS,NSTOKES
DOUBLE PRECISION, intent(in) :: R2_BEAM ( MAXSTREAMS, MAXSTOKES )

INTEGER :: I,O1,ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP2A: IN SUBROUTINE BVP_SOLUTION_MASTER'

write(*,*)
DO O1=1,NSTOKES
  DO I=1,NSTREAMS
    write(*,*) 'O1 = ',O1,' I = ',I,' R2_BEAM(I,O1) = ',R2_BEAM(I,O1)
  ENDDO
ENDDO

END SUBROUTINE TP2A

!

SUBROUTINE TP2B (AWI,NSTKS_NSTRMS_2,IBEAM,COL2)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,NSTKS_NSTRMS_2,IBEAM
DOUBLE PRECISION, intent(in) :: COL2 (MAXSTRMSTKS_2,MAXBEAMS)

INTEGER :: I,ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP2B: IN SUBROUTINE BVP_SOLUTION_MASTER'

write(*,*)
write(*,*) 'IBEAM = ',IBEAM
DO I=1,NSTKS_NSTRMS_2
  write(*,*) 'I = ',I,' COL2(I,IBEAM) = ',COL2(I,IBEAM)
ENDDO

END SUBROUTINE TP2B

!

SUBROUTINE TP2C (AWI,NSTKS_NSTRMS,NLAYERS,LCON,MCON)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,NSTKS_NSTRMS,NLAYERS
DOUBLE PRECISION, intent(in) :: LCON(MAXSTRMSTKS,MAXLAYERS),&
                                MCON(MAXSTRMSTKS,MAXLAYERS)

INTEGER :: I,J,ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP2C: IN SUBROUTINE BVP_SOLUTION_MASTER'

write(*,*)
DO J=1,NLAYERS
  DO I=1,NSTKS_NSTRMS
    write(*,*) 'I = ',I,' J = ',J,' LCON(I,J) = ',LCON(I,J)
  ENDDO
ENDDO
write(*,*)
DO J=1,NLAYERS
  DO I=1,NSTKS_NSTRMS
    write(*,*) 'I = ',I,' J = ',J,' MCON(I,J) = ',MCON(I,J)
  ENDDO
ENDDO

END SUBROUTINE TP2C

!******************************************************************************
!******************************************************************************
!                   Section: TP3*  --> vlidort_corrections.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP4*  --> vlidort_geometry.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP5*  --> vlidort_getplanck.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP6*  --> vlidort_inputs.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP7*  --> vlidort_intensity.f90
!******************************************************************************
!******************************************************************************

SUBROUTINE TP7A (AWI,N,NC,UM,O1,LAYER_SOURCE,T_DELT_USERM,CUMSOURCE_UP)

! VLIDORT_UPUSER_INTENSITY cumulative source

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,N,NC,UM,O1
DOUBLE PRECISION, intent(in) :: LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES ),&
                                T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS ),&
                                CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7A: IN SUBROUTINE UPUSER_INTENSITY'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' O1 = ',O1
write(*,*) 'LAYER_SOURCE(UM,O1)      = ',LAYER_SOURCE(UM,O1)
write(*,*) 'T_DELT_USERM(N,UM)       = ',T_DELT_USERM(N,UM)
write(*,*) 'CUMSOURCE_UP(UM,O1,NC-1) = ',CUMSOURCE_UP(UM,O1,NC-1)

END SUBROUTINE TP7A

!

SUBROUTINE TP7B1 (AWI,UTA,UM,IBEAM,UT,NC,O1,FLUX_MULTIPLIER,&
                  LAYER_SOURCE,T_UTUP_USERM,CUMSOURCE_UP,STOKES_F)

! VLIDORT_UPUSER_INTENSITY Offgrid output

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,UTA,UM,IBEAM,UT,NC,O1
DOUBLE PRECISION, intent(in) :: FLUX_MULTIPLIER,&
                                LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES ),&
                                T_UTUP_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS ),&
                                CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS ),&
                                STOKES_F     ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7B1: IN SUBROUTINE UPUSER_INTENSITY'

write(*,*)
write(*,*) 'UTA = ',UTA,' UM = ',UM,' IBEAM = ',IBEAM,' UT = ',UT,' O1 = ',O1
write(*,*) 'FLUX_MULTIPLIER        = ',FLUX_MULTIPLIER
write(*,*) 'LAYER_SOURCE(UM,O1)    = ',LAYER_SOURCE(UM,O1)
write(*,*) 'T_UTUP_USERM(UT,UM)    = ',T_UTUP_USERM(UT,UM)
write(*,*) 'CUMSOURCE_UP(UM,O1,NC) = ',CUMSOURCE_UP(UM,O1,NC)

write(*,*) 'STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = ',STOKES_F(UTA,UM,IBEAM,O1,UPIDX)

END SUBROUTINE TP7B1

!

SUBROUTINE TP7B2 (AWI,UTA,UM,IBEAM,NC,O1,FLUX_MULTIPLIER,CUMSOURCE_UP,STOKES_F)

! VLIDORT_UPUSER_INTENSITY Ongrid output

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,UTA,UM,IBEAM,NC,O1
DOUBLE PRECISION, intent(in) :: FLUX_MULTIPLIER,&
                                CUMSOURCE_UP ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS ),&
                                STOKES_F     ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7B2: IN SUBROUTINE UPUSER_INTENSITY'

write(*,*)
write(*,*) 'UTA = ',UTA,' UM = ',UM,' IBEAM = ',IBEAM,' O1 = ',O1
write(*,*) 'FLUX_MULTIPLIER        = ',FLUX_MULTIPLIER
write(*,*) 'CUMSOURCE_UP(UM,O1,NC) = ',CUMSOURCE_UP(UM,O1,NC)

write(*,*) 'STOKES_F(UTA,UM,IBEAM,O1,UPIDX) = ',STOKES_F(UTA,UM,IBEAM,O1,UPIDX)

END SUBROUTINE TP7B2

!

!Note: not active
SUBROUTINE TP7C (N,NC,UM,LAYER_SOURCE,T_DELT_USERM,CUMSOURCE_DN)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,NC,UM
DOUBLE PRECISION, intent(in) :: LAYER_SOURCE ( MAX_USER_STREAMS ),&
                                T_DELT_USERM ( MAXLAYERS, MAX_USER_STREAMS ),&
                                CUMSOURCE_DN ( MAX_USER_STREAMS, 0:MAXLAYERS )

write(*,*)
write(*,*) 'TP7C: IN SUBROUTINE DNUSER_INTENSITY'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM
write(*,*) 'LAYER_SOURCE(UM)      = ',LAYER_SOURCE(UM)
write(*,*) 'T_DELT_USERM(N,UM)    = ',T_DELT_USERM(N,UM)
write(*,*) 'CUMSOURCE_DN(UM,NC-1) = ',CUMSOURCE_DN(UM,NC-1)

END SUBROUTINE TP7C

!

SUBROUTINE TP7D1 (AWI,UTA,UM,IBEAM,UT,NC,O1,FLUX_MULTIPLIER,&
                  LAYER_SOURCE,T_UTDN_USERM,CUMSOURCE_DN,STOKES_F)

! VLIDORT_DNUSER_INTENSITY Partial layer source

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,UTA,UM,IBEAM,UT,NC,O1
DOUBLE PRECISION, intent(in) :: FLUX_MULTIPLIER,&
                                LAYER_SOURCE ( MAX_USER_STREAMS, MAXSTOKES ),&
                                T_UTDN_USERM ( MAX_PARTLAYERS, MAX_USER_STREAMS ),&
                                CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS ),&
                                STOKES_F  ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7D1: IN SUBROUTINE DNUSER_INTENSITY'

write(*,*)
write(*,*) 'UTA = ',UTA,' UM = ',UM,' IBEAM = ',IBEAM,' UT = ',UT,' O1 = ',O1
write(*,*) 'FLUX_MULTIPLIER        = ',FLUX_MULTIPLIER
write(*,*) 'LAYER_SOURCE(UM,O1)    = ',LAYER_SOURCE(UM,O1)
write(*,*) 'T_UTDN_USERM(UT,UM)    = ',T_UTDN_USERM(UT,UM)
write(*,*) 'CUMSOURCE_DN(UM,O1,NC) = ',CUMSOURCE_DN(UM,O1,NC)

write(*,*) 'STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = ',STOKES_F(UTA,UM,IBEAM,O1,DNIDX)

END SUBROUTINE TP7D1

!

!Note: not active
SUBROUTINE TP7D2 (UTA,UM,IBEAM,NC,O1,FLUX_MULTIPLIER,CUMSOURCE_DN,STOKES_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: UTA,UM,IBEAM,NC,O1
DOUBLE PRECISION, intent(in) :: FLUX_MULTIPLIER,&
                                CUMSOURCE_DN ( MAX_USER_STREAMS, MAXSTOKES, 0:MAXLAYERS ),&
                                STOKES_F     ( MAX_USER_LEVELS, MAX_USER_STREAMS, MAXBEAMS, MAXSTOKES, MAX_DIRECTIONS )

write(*,*)
write(*,*) 'TP7D2: IN SUBROUTINE DNUSER_INTENSITY'

write(*,*)
write(*,*) 'UTA = ',UTA,' UM = ',UM,' IBEAM = ',IBEAM,' O1 = ',O1
write(*,*) 'FLUX_MULTIPLIER        = ',FLUX_MULTIPLIER
write(*,*) 'CUMSOURCE_DN(UM,O1,NC) = ',CUMSOURCE_DN(UM,O1,NC)

write(*,*) 'STOKES_F(UTA,UM,IBEAM,O1,DNIDX) = ',STOKES_F(UTA,UM,IBEAM,O1,DNIDX)

END SUBROUTINE TP7D2

!

SUBROUTINE TP7E (AWI,FOURIER,UT,V,O,W,STOKES,STOKES_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, INTENT(IN)   :: AWI,FOURIER,UT,V,O,W
DOUBLE PRECISION, INTENT (IN) :: &
                         STOKES    (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS ), &
                         STOKES_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

if ( ( ut==1) .and. ( v==1) .and. ( o==2) .and. ( w==1) ) then
write(*,*)
write(*,*) 'TP7E: IN SUBROUTINE VLIDORT_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' V = ',V,' O = ',O,' W = ',W
write(*,*) 'STOKES(UT,V,O,W) before = ',STOKES(UT,V,O,W)
write(*,*) 'STOKES_SS(UT,V,O,W)     = ',STOKES_SS(UT,V,O,W)
write(*,*) 'STOKES(UT,V,O,W) after  = ',STOKES(UT,V,O,W) + STOKES_SS(UT,V,O,W)
endif

END SUBROUTINE TP7E

!

SUBROUTINE TP7F (AWI,FOURIER,UT,V,O,STOKES,STOKES_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,FOURIER,UT,V,O
DOUBLE PRECISION, intent(in) :: &
                         STOKES    (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
                         STOKES_DB (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

if ( ( ut==1) .and. ( v==1) .and. ( o==2) ) then
write(*,*)
write(*,*) 'TP7F: IN SUBROUTINE VLIDORT_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' V = ',V,' O = ',O
write(*,*) 'STOKES(UT,V,O,UPIDX) before = ',STOKES(UT,V,O,UPIDX)
write(*,*) 'STOKES_DB(UT,V,O)           = ',STOKES_DB(UT,V,O)
write(*,*) 'STOKES(UT,V,O,UPIDX) after  = ',STOKES(UT,V,O,UPIDX) + STOKES_DB(UT,V,O)
endif

END SUBROUTINE TP7F

!

SUBROUTINE TP7G (AWI,FOURIER,UT,I,IBEAM,V,O,W,STOKES,STOKES_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,FOURIER,UT,I,IBEAM,V,O,W
DOUBLE PRECISION, intent(in) :: STOKES    (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
                                STOKES_F  (MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7G: IN SUBROUTINE VLIDORT_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' I = ',I,' IBEAM = ',IBEAM,' V = ',V,' O = ',O,' W = ',W
write(*,*) 'STOKES(UT,V,O,W) before  = ',STOKES(UT,V,O,W)
write(*,*) 'STOKES_F(UT,I,IBEAM,O,W) = ',STOKES_F(UT,I,IBEAM,O,W)

END SUBROUTINE TP7G

!

SUBROUTINE TP7E2 (AWI,FOURIER,UT,IBEAM,O,W,STOKES,STOKES_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,FOURIER,UT,IBEAM,O,W
DOUBLE PRECISION, intent(in) :: STOKES    (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
                                STOKES_SS (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS)

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7E2: IN SUBROUTINE LIDORT_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' IBEAM = ',IBEAM,' W = ',W
write(*,*) 'STOKES(UT,IBEAM,O,W) before = ',STOKES(UT,IBEAM,O,W)
write(*,*) 'STOKES_SS(UT,IBEAM,O,W)     = ',STOKES_SS(UT,IBEAM,O,W)
write(*,*) 'STOKES(UT,IBEAM,O,W) after  = ',STOKES(UT,IBEAM,O,W) + STOKES_SS(UT,IBEAM,O,W)

END SUBROUTINE TP7E2

!

SUBROUTINE TP7F2 (AWI,FOURIER,UT,IBEAM,O,STOKES,STOKES_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,FOURIER,UT,IBEAM,O
DOUBLE PRECISION, intent(in) :: STOKES    (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
                                STOKES_DB (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7F2: IN SUBROUTINE VLIDORT_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' IBEAM = ',IBEAM,' O = ',O
write(*,*) 'STOKES(UT,IBEAM,O,UPIDX) before = ',STOKES(UT,IBEAM,O,UPIDX)
write(*,*) 'STOKES_DB(UT,IBEAM,O)           = ',STOKES_DB(UT,IBEAM,O)
write(*,*) 'STOKES(UT,IBEAM,O,UPIDX) after  = ',STOKES(UT,IBEAM,O,UPIDX) + STOKES_DB(UT,IBEAM,O)

END SUBROUTINE TP7F2

!

SUBROUTINE TP7G2 (AWI,FOURIER,UT,UM,IBEAM,O,W,STOKES,STOKES_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: AWI,FOURIER,UT,UM,IBEAM,O,W
DOUBLE PRECISION, intent(in) :: STOKES   (MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
                                STOKES_F (MAX_USER_LEVELS,MAX_USER_VZANGLES,MAX_SZANGLES,MAXSTOKES,MAX_DIRECTIONS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP7G2: IN SUBROUTINE VLIDORT_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' UT = ',UT,' UM = ',UM,' IBEAM = ',IBEAM,' O = ',O,' W = ',W
write(*,*) 'STOKES(UT,IBEAM,O,W) before = ',STOKES(UT,IBEAM,O,W)
write(*,*) 'STOKES_F(UT,UM,IBEAM,O,W)  = ',STOKES_F(UT,UM,IBEAM,O,W)

END SUBROUTINE TP7G2

!******************************************************************************
!******************************************************************************
!                   Section: TP8*  --> vlidort_la_miscsetups.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP9*  --> vlidort_lbbf_jacobians_vector.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP10* --> vlidort_lc_bvproblem.f90
!******************************************************************************
!******************************************************************************

!Note: These 3 subroutines from "lidort_testpts.f90" not active

SUBROUTINE TP10 (N,I,AA,Q,NCON,PCON,XPOS,XNEG)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in) :: N,I,AA,Q

DOUBLE PRECISION, intent(in)  :: NCON ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
DOUBLE PRECISION, intent(in)  :: PCON ( MAXSTREAMS, MAXLAYERS, MAX_ATMOSWFS )
DOUBLE PRECISION, intent(in)  :: XPOS ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: XNEG ( MAXSTREAMS_2, MAXSTREAMS, MAXLAYERS )

if ( (N == 1) .and. (I == 1) .and. (AA == 1) .and. (Q == 1) ) then

write(*,*)
write(*,*) 'TP10: IN SUBROUTINE LC_BVP_FULLSOLUTION_MASTER'

write(*,*) 'N = ',N,' I = ',I,' AA = ',AA,' Q = ',Q
write(*,*)
write(*,*) 'NCON(AA,N,Q) = ',NCON(AA,N,Q)
write(*,*) 'PCON(AA,N,Q) = ',PCON(AA,N,Q)
write(*,*)
write(*,*) 'XPOS(I,AA,N) = ',XPOS(I,AA,N)
write(*,*) 'XNEG(I,AA,N) = ',XNEG(I,AA,N)
!stop 'TP10'
end if

END SUBROUTINE TP10

!

SUBROUTINE TP10A (NLAYERS,NSTREAMS_2,NTOTAL,N_WEIGHTFUNCS,COL2_WF)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in) :: NLAYERS,NSTREAMS_2,NTOTAL,N_WEIGHTFUNCS
DOUBLE PRECISION, intent(in) :: COL2_WF (MAXTOTAL,MAX_ATMOSWFS)

INTEGER :: I,Q

write(*,*)
write(*,*) 'TP10A: IN SUBROUTINE LC_BVP_COLUMN_SETUP'

write(*,*) 'NLAYERS       = ',NLAYERS
write(*,*) 'NSTREAMS_2    = ',NSTREAMS_2
write(*,*) 'NTOTAL        = ',NTOTAL
!write(*,*) 'N_WEIGHTFUNCS = ',N_WEIGHTFUNCS

DO I = 1, NTOTAL
  DO Q = 1, 1!N_WEIGHTFUNCS
    IF (MOD(I-1,NSTREAMS_2) == 0) write(*,*) 'LAYER = ',((I-1)/NSTREAMS_2) + 1
    !write(*,*) 'I = ',I,' Q = ',Q,' COL2_WF(I,Q) = ',COL2_WF(I,Q)
    IF (MOD(I-1,NSTREAMS_2) == 0) write(*,*) 'I = ',I,' Q = ',Q,' COL2_WF(I,Q) = ',COL2_WF(I,Q)
  ENDDO
ENDDO
!stop 'TP10A'

END SUBROUTINE TP10A

!

SUBROUTINE TP10B (N,I,AA,Q,NSTREAMS,L_GFUNC_UP,L_GFUNC_DN,L_XPOS,S_P_U,S_M_U,S_P_L,S_M_L)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,I,AA,Q,NSTREAMS

DOUBLE PRECISION, intent(in) :: L_GFUNC_UP(MAXSTREAMS,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: L_GFUNC_DN(MAXSTREAMS,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: L_XPOS(MAXSTREAMS_2,MAXSTREAMS,MAXLAYERS,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: S_P_U, S_P_L, S_M_U, S_M_L

INTEGER :: I1

if ( (N == 1) .and. (I == 1).and. (AA == 1) .and. (Q == 1) ) then
I1 = I + NSTREAMS

write(*,*)
write(*,*) 'TP10B: IN SUBROUTINE LC_BEAMSOLUTION_NEQK'

write(*,*) 'N = ',N,' AA = ',AA,' Q = ',Q
write(*,*) 'L_GFUNC_UP(AA,Q)  = ',L_GFUNC_UP(AA,Q)
write(*,*) 'L_GFUNC_DN(AA,Q)  = ',L_GFUNC_DN(AA,Q)
write(*,*) 'L_XPOS(I,AA,N,Q)  = ',L_XPOS(I,AA,N,Q)
write(*,*) 'L_XPOS(I1,AA,N,Q) = ',L_XPOS(I1,AA,N,Q)
write(*,*) 'S_P_U = ',S_P_U
write(*,*) 'S_P_L = ',S_P_L
write(*,*) 'S_M_U = ',S_M_U
write(*,*) 'S_M_L = ',S_M_L
!stop 'TP10B'
end if

END SUBROUTINE TP10B

!******************************************************************************
!******************************************************************************
!                   Section: TP11* --> vlidort_lc_corrections.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP12* --> vlidort_lc_miscsetups.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP13* --> vlidort_lc_PostProcessing.f90
!******************************************************************************
!******************************************************************************

!Note: These 4 subroutines from "lidort_testpts.f90" not active

SUBROUTINE TP13 (N,Q,UTA,I,IB,HOM1,HOM2,HOM3,HOM4,HOM5,FTERM,SPAR,SHOM,QUADCOLUMNWF )

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,Q,UTA,I,IB
DOUBLE PRECISION, intent(in) :: HOM1,HOM2,HOM3,HOM4,HOM5,FTERM,SPAR,SHOM
DOUBLE PRECISION, intent(in) :: QUADCOLUMNWF ( MAX_ATMOSWFS, MAX_USER_LEVELS, &
                                               MAXSTREAMS,   MAXBEAMS, MAX_DIRECTIONS )

if ( (Q == 1) .and. (UTA == 1) .and. (I == 1) .and. (IB == 1) ) then

write(*,*)
write(*,*) 'TP13: IN SUBROUTINE QUADCOLUMNWF_LEVEL_UP'

write(*,*) 'upwelling:'
write(*,*) 'N = ',N
write(*,*) 'Q = ',Q,' UTA = ',UTA,' I = ',I,' IB = ',IB
write(*,*)
write(*,*) 'HOM1 = ',HOM1
write(*,*) 'HOM2 = ',HOM2
write(*,*) 'HOM3 = ',HOM3
write(*,*) 'HOM4 = ',HOM4
write(*,*) 'HOM5 = ',HOM5
write(*,*)
write(*,*) 'FM   = ',FTERM
write(*,*) 'SPAR = ',SPAR
write(*,*) 'SHOM = ',SHOM
write(*,*) 'QUADCOLUMNWF(Q,UTA,I,IB,UPIDX) = ',QUADCOLUMNWF(Q,UTA,I,IB,UPIDX)
!stop 'TP13'
end if

END SUBROUTINE TP13

!

SUBROUTINE TP13A (N,Q,UM,IB,H1,H2,H3,H4,H5,H6,SHOM,L_LAYER_SOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,Q,UM,IB
DOUBLE PRECISION, intent(in) :: H1,H2,H3,H4,H5,H6,SHOM
DOUBLE PRECISION, intent(in) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

if ( (N == 1) .and. (Q == 1) .and. (UM == 1) .and. (IB == 1) ) then

write(*,*)
write(*,*) 'TP13A: IN SUBROUTINE LC_WHOLELAYER_STERM_DN'

write(*,*) 'dnwelling:'
write(*,*) 'N = ',N
write(*,*) 'Q = ',Q,' UM = ',UM,' IB = ',IB
write(*,*)
write(*,*) 'last H1 = ',H1
write(*,*) 'last H2 = ',H2
write(*,*) 'last H3 = ',H3
write(*,*) 'last H4 = ',H4
write(*,*) 'last H5 = ',H5
write(*,*) 'last H6 = ',H6
write(*,*)
write(*,*) 'SHOM = ',SHOM
write(*,*) 'L_LAYER_SOURCE(UM,Q) = ',L_LAYER_SOURCE(UM,Q)
!stop 'TP13A'
end if

END SUBROUTINE TP13A

!

SUBROUTINE TP13B (N,Q,AA,UM,IB,H1,H2,SPAR)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,Q,AA,UM,IB
DOUBLE PRECISION, intent(in) :: H1,H2,SPAR

if ( (N == 1) .and. (AA == 1) .and. (Q == 1) .and. (UM == 1) .and. (IB == 1) ) then

write(*,*)
write(*,*) 'TP13B: IN SUBROUTINE LC_WHOLELAYER_STERM_DN'

write(*,*) 'dnwelling:'
write(*,*) 'N = ',N
write(*,*) 'Q = ',Q,' AA = ',AA,' UM = ',UM,' IB = ',IB
write(*,*)
write(*,*) 'last H1 = ',H1
write(*,*) 'last H2 = ',H2
write(*,*)
write(*,*) 'SPAR = ',SPAR
!stop 'TP13B'
end if

END SUBROUTINE TP13B

!

SUBROUTINE TP13C (N,Q,LUM,UM,IB,IC,U_WNEG,LC_U_WNEG,EMULT_DN,LC_EMULT_DN,SFOR,L_LAYER_SOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: N,Q,LUM,UM,IB,IC
DOUBLE PRECISION, intent(in) :: U_WNEG(MAX_USER_STREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in) :: LC_U_WNEG(MAX_USER_STREAMS,MAXLAYERS,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,2)
DOUBLE PRECISION, intent(in) :: LC_EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS,2,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: SFOR
DOUBLE PRECISION, intent(in) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

!if ( (N == 1) .and. (Q == 1) .and. (UM == 1) .and. (IB == 1) ) then
if ( (N == 1) .and. (Q == 1) .and. (UM == 1) ) then

write(*,*)
write(*,*) 'TP13C: IN SUBROUTINE LC_WHOLELAYER_STERM_DN'

write(*,*) 'dnwelling:'
write(*,*) 'N = ',N
write(*,*) 'Q = ',Q,' UM = ',UM,' IB = ',IB,' IC = ',IC
write(*,*)
write(*,*) 'U_WNEG(UM,N)               = ',U_WNEG(UM,N)
write(*,*) 'LC_U_WNEG(UM,N,Q)          = ',LC_U_WNEG(UM,N,Q)
write(*,*) 'EMULT_DN(LUM,N,IB,IC)      = ',EMULT_DN(LUM,N,IB,IC)
write(*,*) 'LC_EMULT_DN(LUM,N,IB,IC,Q) = ',LC_EMULT_DN(LUM,N,IB,IC,Q)
write(*,*) 'SFOR                       = ',SFOR
write(*,*)
write(*,*) 'L_LAYER_SOURCE(UM,Q)       = ',L_LAYER_SOURCE(UM,Q)
!stop 'TP13C'
end if

END SUBROUTINE TP13C

!******************************************************************************
!******************************************************************************
!                   Section: TP14* --> vlidort_lcs_masters_V2p8.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP15* --> vlidort_lc_wfatmos.f90
!******************************************************************************
!******************************************************************************

!Note: TP15, TP15A, TP15B subroutines from "lidort_testpts.f90" not active

SUBROUTINE TP15 (IBEAM,N,NC,Q,UM,IC,T_DELT_USERM,L_T_DELT_USERM,L_LAYER_SOURCE,&
                 CUMSOURCE_DN,L_CUMUL_SOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: IBEAM,N,NC,Q,UM,IC

DOUBLE PRECISION, intent(in) :: T_DELT_USERM   (MAXLAYERS,MAX_USER_STREAMS,2)
DOUBLE PRECISION, intent(in) :: L_T_DELT_USERM (MAXLAYERS,MAX_USER_STREAMS,2,MAX_ATMOSWFS)

DOUBLE PRECISION, intent(in) :: L_LAYER_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)
DOUBLE PRECISION, intent(in) :: CUMSOURCE_DN   (MAX_USER_STREAMS,0:MAXLAYERS)
DOUBLE PRECISION, intent(in) :: L_CUMUL_SOURCE (MAX_USER_STREAMS,MAX_ATMOSWFS)

if ( (IBEAM == 1) .and. (N == 1) .and. (Q == 1) .and. (UM == 1) ) then

write(*,*)
write(*,*) 'TP15: IN SUBROUTINE DNUSER_COLUMNWF'

write(*,*) 'dnwelling:'
write(*,*) 'IBEAM = ',IBEAM
write(*,*) 'N = ',N,' Q = ',Q,' UM = ',UM,' IC = ',IC
write(*,*)
write(*,*) 'L_LAYER_SOURCE(UM,Q)      = ',L_LAYER_SOURCE(UM,Q)
write(*,*) 'T_DELT_USERM(N,UM,IC)     = ',T_DELT_USERM(N,UM,IC)
write(*,*) 'L_CUMUL_SOURCE(UM,Q)      = ',L_CUMUL_SOURCE(UM,Q)
write(*,*) 'L_T_DELT_USERM(N,UM,IC,Q) = ',L_T_DELT_USERM(N,UM,IC,Q)
write(*,*) 'CUMSOURCE_DN(UM,NC-1)     = ',CUMSOURCE_DN(UM,NC-1)
!stop 'TP15'
end if

END SUBROUTINE TP15

!

SUBROUTINE TP15A (IB,N,Q,ISCENE,L_BEAM,FAC,L_BOA_DBSOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: IB,N,Q,ISCENE

DOUBLE PRECISION, intent(in) :: L_BEAM
DOUBLE PRECISION, intent(in) :: FAC
DOUBLE PRECISION, intent(in) :: L_BOA_DBSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

!if ( (IB == 1) .and. (N == 23) .and. (Q == 1) ) then
if ( (N == 23) .and. (Q == 1) ) then

write(*,*)
write(*,*) 'TP15A: IN SUBROUTINE GET_LC_BOASOURCE'

write(*,*) 'upwelling:'
write(*,*) 'IB = ',IB
write(*,*) 'N  = ',N,' Q = ',Q,' ISCENE = ',ISCENE
write(*,*)
write(*,*) 'L_BEAM = ',L_BEAM
write(*,*) 'FAC    = ',FAC
write(*,*) 'L_BOA_DBSOURCE(IB,Q) = ',L_BOA_DBSOURCE(IB,Q)
!stop 'TP15A'
end if

END SUBROUTINE TP15A

!

SUBROUTINE TP15B (IB,Q,REFLEC,L_BOA_MSSOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: IB,Q

DOUBLE PRECISION, intent(in) :: REFLEC
DOUBLE PRECISION, intent(in) :: L_BOA_MSSOURCE(MAX_USER_STREAMS,MAX_ATMOSWFS)

if (Q == 1) then

write(*,*)
write(*,*) 'TP15B: IN SUBROUTINE GET_LC_BOASOURCE'

write(*,*) 'upwelling:'
write(*,*) 'IB = ',IB,' Q = ',Q
write(*,*)
write(*,*) 'REFLEC = ',REFLEC
write(*,*) 'L_BOA_MSSOURCE(IB,Q) = ',L_BOA_MSSOURCE(IB,Q)
!stop 'TP15B'
end if

END SUBROUTINE TP15B

!

SUBROUTINE TP15E (FOURIER,Q,UT,V,O1,W,COLUMNWF,COLUMNWF_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,V,O1,W
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF    (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_SS (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS)

if (  (Q == 2) .and. (UT == 4) .and. (V == 1) .and. (O1 == 1) .and. (W == 1) ) then

write(*,*)
write(*,*) 'TP15E: IN SUBROUTINE VLIDORT_LCS_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' V = ',V,' O1 = ',O1,' W = ',W
write(*,*) 'COLUMNWF(Q,UT,V,O1,W) before = ',COLUMNWF(Q,UT,V,O1,W)
write(*,*) 'COLUMNWF_SS(Q,UT,V,O1,W)     = ',COLUMNWF_SS(Q,UT,V,O1,W)
write(*,*) 'COLUMNWF(Q,UT,V,O1,W) after  = ',COLUMNWF(Q,UT,V,O1,W) + COLUMNWF_SS(Q,UT,V,O1,W)
end if

END SUBROUTINE TP15E

!

SUBROUTINE TP15F (FOURIER,Q,UT,V,O1,COLUMNWF,COLUMNWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,V,O1
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF    (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_DB (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

if (  (Q == 2) .and. (UT == 4) .and. (V == 1) .and. (O1 == 1) ) then

write(*,*)
write(*,*) 'TP15F: IN SUBROUTINE VLIDORT_LCS_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' V = ',V,' O1 = ',O1
write(*,*) 'COLUMNWF(Q,UT,V,O1,UPIDX) before = ',COLUMNWF(Q,UT,V,O1,UPIDX)
write(*,*) 'COLUMNWF_DB(Q,UT,V,O1)           = ',COLUMNWF_DB(Q,UT,V,O1)
write(*,*) 'COLUMNWF(Q,UT,V,O1,UPIDX) after  = ',COLUMNWF(Q,UT,V,O1,UPIDX) + COLUMNWF_DB(Q,UT,V,O1)
endif

END SUBROUTINE TP15F

!

SUBROUTINE TP15G (FOURIER,Q,UT,I,IBEAM,V,O1,W,COLUMNWF,COLUMNWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,I,IBEAM,V,O1,W
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF   (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_F (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

if (  (Q == 1) .and. (UT == 1) .and. (I == 1) .and. &
      (IBEAM == 1) .and. (V == 1) .and. (O1 == 1) .and. (W == 1) ) then

write(*,*)
write(*,*) 'TP15G: IN SUBROUTINE VLIDORT_LCS_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' I = ',I,' IBEAM = ',IBEAM,' V = ',V,' O1 = ',O1,' W = ',W
write(*,*) 'COLUMNWF(Q,UT,V,O1,W) before  = ',COLUMNWF(Q,UT,V,O1,W)
write(*,*) 'COLUMNWF_F(Q,UT,I,IBEAM,O1,W) = ',COLUMNWF_F(Q,UT,I,IBEAM,O1,W)
endif

END SUBROUTINE TP15G

!

SUBROUTINE TP15E2 (FOURIER,Q,UT,IB,O1,W,COLUMNWF,COLUMNWF_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,IB,O1,W
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF    (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_SS (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP15E2: IN SUBROUTINE VLIDORT_LC_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' IB = ',IB,' O1 = ',O1,' W = ',W
write(*,*) 'COLUMNWF(Q,UT,IB,O1,W) before = ',COLUMNWF(Q,UT,IB,O1,W)
write(*,*) 'COLUMNWF_SS(Q,UT,IB,O1,W)     = ',COLUMNWF_SS(Q,UT,IB,O1,W)
write(*,*) 'COLUMNWF(Q,UT,IB,O1,W) after  = ',COLUMNWF(Q,UT,IB,O1,W) + COLUMNWF_SS(Q,UT,IB,O1,W)

END SUBROUTINE TP15E2

!

SUBROUTINE TP15F2 (FOURIER,Q,UT,IB,O1,COLUMNWF,COLUMNWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,IB,O1
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF    (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_DB (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

write(*,*)
write(*,*) 'TP15F2: IN SUBROUTINE VLIDORT_LC_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' IB = ',IB,' O1 = ',O1
write(*,*) 'COLUMNWF(Q,UT,IB,O1,UPIDX) before = ',COLUMNWF(Q,UT,IB,O1,UPIDX)
write(*,*) 'COLUMNWF_DB(Q,UT,IB,O1)           = ',COLUMNWF_DB(Q,UT,IB,O1)
write(*,*) 'COLUMNWF(Q,UT,IB,O1,UPIDX) after  = ',COLUMNWF(Q,UT,IB,O1,UPIDX) + COLUMNWF_DB(Q,UT,IB,O1)

END SUBROUTINE TP15F2

!

SUBROUTINE TP15G2 (FOURIER,Q,UT,LUM,IB,O1,W,COLUMNWF,COLUMNWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,UT,LUM,IB,O1,W
DOUBLE PRECISION, intent(in) :: &
  COLUMNWF   (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  COLUMNWF_F (MAX_ATMOSWFS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP15G2: IN SUBROUTINE VLIDORT_LC_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' UT = ',UT,' LUM = ',LUM,' IB = ',IB,' O1 = ',O1,' W = ',W
write(*,*) 'COLUMNWF(Q,UT,IB,O1,W) before = ',COLUMNWF(Q,UT,IB,O1,W)
write(*,*) 'COLUMNWF_F(Q,UT,LUM,IB,O1,W)  = ',COLUMNWF_F(Q,UT,LUM,IB,O1,W)

END SUBROUTINE TP15G2

!******************************************************************************
!******************************************************************************
!                   Section: TP16* --> vlidort_l_inputs.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP17* --> vlidort_lp_bvproblem.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP18* --> vlidort_lp_corrections.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP19* --> vlidort_lpc_solutions.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP20* --> vlidort_lp_miscsetups.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP21* --> vlidort_lp_PostProcessing.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP22* --> vlidort_lps_masters_V2p8.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP23* --> vlidort_lp_wfatmos.f90
!******************************************************************************
!******************************************************************************

SUBROUTINE TP23E (FOURIER,Q,N,UT,V,O1,W,PROFILEWF,PROFILEWF_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,V,O1,W
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF    (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_SS (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS)

if ( (q == 2) .and. (n == 23) .and. (ut == 4) .and. (v == 1) .and. (o1 == 1) .and. (w == 1) ) then

write(*,*)
write(*,*) 'TP23E: IN SUBROUTINE VLIDORT_LP_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' V = ',V,' O1 = ',O1,' W = ',W
write(*,*) 'PROFILEWF(Q,N,UT,V,O1,W) before = ',PROFILEWF(Q,N,UT,V,O1,W)
write(*,*) 'PROFILEWF_SS(Q,N,UT,V,O1,W)     = ',PROFILEWF_SS(Q,N,UT,V,O1,W)
write(*,*) 'PROFILEWF(Q,N,UT,V,O1,W) after  = ',PROFILEWF(Q,N,UT,V,O1,W) + PROFILEWF_SS(Q,N,UT,V,O1,W)

endif

END SUBROUTINE TP23E

!

SUBROUTINE TP23F (FOURIER,Q,N,UT,V,O1,PROFILEWF,PROFILEWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,V,O1
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF    (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_DB (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

if ( (q == 2) .and. (n == 23) .and. (ut == 4) .and. (v == 1) .and. (o1 == 1) ) then

write(*,*)
write(*,*) 'TP23F: IN SUBROUTINE VLIDORT_LP_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' V = ',V,' O1 = ',O1
write(*,*) 'PROFILEWF(Q,N,UT,V,O1,UPIDX) before = ',PROFILEWF(Q,N,UT,V,O1,UPIDX)
write(*,*) 'PROFILEWF_DB(Q,N,UT,V,O1)           = ',PROFILEWF_DB(Q,N,UT,V,O1)
write(*,*) 'PROFILEWF(Q,N,UT,V,O1,UPIDX) after  = ',PROFILEWF(Q,N,UT,V,O1,UPIDX) + PROFILEWF_DB(Q,N,UT,V,O1)

endif

END SUBROUTINE TP23F

!

SUBROUTINE TP23G (FOURIER,Q,N,UT,I,IBEAM,V,O1,W,PROFILEWF,PROFILEWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,I,IBEAM,V,O1,W
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF   (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_F (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP23G: IN SUBROUTINE VLIDORT_LP_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' I = ',I,' IBEAM = ',IBEAM,' V = ',V,' O1 = ',O1,' W = ',W
write(*,*) 'PROFILEWF(Q,N,UT,V,O1,W) before  = ',PROFILEWF(Q,N,UT,V,O1,W)
write(*,*) 'PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W) = ',PROFILEWF_F(Q,N,UT,I,IBEAM,O1,W)

END SUBROUTINE TP23G

!

SUBROUTINE TP23E2 (FOURIER,Q,N,UT,IB,O1,W,PROFILEWF,PROFILEWF_SS)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,IB,O1,W
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF    (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_SS (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP23E2: IN SUBROUTINE VLIDORT_LP_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' IB = ',IB,' O1 = ',O1,' W = ',W
write(*,*) 'PROFILEWF(Q,N,UT,IB,O1,W) before = ',PROFILEWF(Q,N,UT,IB,O1,W)
write(*,*) 'PROFILEWF_SS(Q,N,UT,IB,O1,W)     = ',PROFILEWF_SS(Q,N,UT,IB,O1,W)
write(*,*) 'PROFILEWF(Q,N,UT,IB,O1,W) after  = ',PROFILEWF(Q,N,UT,IB,O1,W) + PROFILEWF_SS(Q,N,UT,IB,O1,W)

END SUBROUTINE TP23E2

!

SUBROUTINE TP23F2 (FOURIER,Q,N,UT,IB,O1,PROFILEWF,PROFILEWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,IB,O1
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF    (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_DB (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

write(*,*)
write(*,*) 'TP23F2: IN SUBROUTINE VLIDORT_LP_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' IB = ',IB,' O1 = ',O1
write(*,*) 'PROFILEWF(Q,N,UT,IB,O1,UPIDX) before = ',PROFILEWF(Q,N,UT,IB,O1,UPIDX)
write(*,*) 'PROFILEWF_DB(Q,N,UT,IB,O1)           = ',PROFILEWF_DB(Q,N,UT,IB,O1)
write(*,*) 'PROFILEWF(Q,N,UT,IB,O1,UPIDX) after  = ',PROFILEWF(Q,N,UT,IB,O1,UPIDX) + PROFILEWF_DB(Q,N,UT,IB,O1)

END SUBROUTINE TP23F2

!

SUBROUTINE TP23G2 (FOURIER,Q,N,UT,LUM,IB,O1,W,PROFILEWF,PROFILEWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Q,N,UT,LUM,IB,O1,W
DOUBLE PRECISION, intent(in) :: &
  PROFILEWF   (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  PROFILEWF_F (MAX_ATMOSWFS,MAXLAYERS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP23G2: IN SUBROUTINE VLIDORT_LP_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Q = ',Q,' N = ',N,' UT = ',UT,' LUM = ',LUM,' IB = ',IB,' O1 = ',O1,' W = ',W
write(*,*) 'PROFILEWF(Q,N,UT,IB,O1,W) before = ',PROFILEWF(Q,N,UT,IB,O1,W)
write(*,*) 'PROFILEWF_F(Q,N,UT,LUM,IB,O1,W)  = ',PROFILEWF_F(Q,N,UT,LUM,IB,O1,W)

END SUBROUTINE TP23G2

!******************************************************************************
!******************************************************************************
!                   Section: TP24* --> vlidort_ls_corrections.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP25* --> vlidort_ls_wfsleave.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP26* --> vlidort_ls_wfsurface.f90
!******************************************************************************
!******************************************************************************

!Note: The call to this subroutine actually found inside "vlidort_lc_wfatmos.f90"
SUBROUTINE TP26F (FOURIER,Z,UT,V,O1,SURFACEWF,SURFACEWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Z,UT,V,O1
DOUBLE PRECISION, intent(in) :: &
  SURFACEWF    (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  SURFACEWF_DB (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

write(*,*)
write(*,*) 'TP26F: IN SUBROUTINE VLIDORT_LCS_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Z = ',Z,' UT = ',UT,' V = ',V,' O1 = ',O1
write(*,*) 'SURFACEWF(Z,UT,V,O1,UPIDX) before = ',SURFACEWF(Z,UT,V,O1,UPIDX)
write(*,*) 'SURFACEWF_DB(Z,UT,V,O1)           = ',SURFACEWF_DB(Z,UT,V,O1)
write(*,*) 'SURFACEWF(Z,UT,V,O1,UPIDX) after  = ',SURFACEWF(Z,UT,V,O1,UPIDX) + SURFACEWF_DB(Z,UT,V,O1)

END SUBROUTINE TP26F

!

!Note: The call to this subroutine actually found inside "vlidort_lc_wfatmos.f90"
SUBROUTINE TP26G (FOURIER,Z,UT,I,IBEAM,V,O1,W,SURFACEWF,SURFACEWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Z,UT,I,IBEAM,V,O1,W
DOUBLE PRECISION, intent(in) :: &
  SURFACEWF   (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  SURFACEWF_F (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP26G: IN SUBROUTINE VLIDORT_LCS_CONVERGE'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Z = ',Z,' UT = ',UT,' I = ',I,' IBEAM = ',IBEAM,' V = ',V,' O1 = ',O1,' W = ',W
write(*,*) 'SURFACEWF(Z,UT,V,O1,W) before  = ',SURFACEWF(Z,UT,V,O1,W)
write(*,*) 'SURFACEWF_F(Z,UT,I,IBEAM,O1,W) = ',SURFACEWF_F(Z,UT,I,IBEAM,O1,W)

END SUBROUTINE TP26G

!

!Note: The call to this subroutine actually found inside "vlidort_lc_wfatmos.f90"
SUBROUTINE TP26F2 (FOURIER,Z,UT,IB,O1,SURFACEWF,SURFACEWF_DB)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Z,UT,IB,O1
DOUBLE PRECISION, intent(in) :: &
  SURFACEWF    (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  SURFACEWF_DB (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES)

write(*,*)
write(*,*) 'TP26F2: IN SUBROUTINE VLIDORT_LCS_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Z = ',Z,' UT = ',UT,' IB = ',IB,' O1 = ',O1
write(*,*) 'SURFACEWF(Z,UT,IB,O1,UPIDX) before = ',SURFACEWF(Z,UT,IB,O1,UPIDX)
write(*,*) 'SURFACEWF_DB(Z,UT,IB,O1)           = ',SURFACEWF_DB(Z,UT,IB,O1)
write(*,*) 'SURFACEWF(Z,UT,IB,O1,UPIDX) after  = ',SURFACEWF(Z,UT,IB,O1,UPIDX) + SURFACEWF_DB(Z,UT,IB,O1)

END SUBROUTINE TP26F2

!

!Note: The call to this subroutine actually found inside "vlidort_lc_wfatmos.f90"
SUBROUTINE TP26G2 (FOURIER,Z,UT,LUM,IB,O1,W,SURFACEWF,SURFACEWF_F)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)   :: FOURIER,Z,UT,LUM,IB,O1,W
DOUBLE PRECISION, intent(in) :: &
  SURFACEWF   (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_GEOMETRIES,MAXSTOKES,MAX_DIRECTIONS),&
  SURFACEWF_F (MAX_SURFACEWFS,MAX_USER_LEVELS,MAX_USER_STREAMS,MAXBEAMS,MAXSTOKES,MAX_DIRECTIONS)

write(*,*)
write(*,*) 'TP26G2: IN SUBROUTINE VLIDORT_LCS_CONVERGE_OBSGEO'

write(*,*)
write(*,*) 'FOURIER = ',FOURIER,' Z = ',Z,' UT = ',UT,' LUM = ',LUM,' IB = ',IB,' O1 = ',O1,' W = ',W
write(*,*) 'SURFACEWF(Z,UT,IB,O1,W) before  = ',SURFACEWF(Z,UT,IB,O1,W)
write(*,*) 'SURFACEWF_F(Z,UT,LUM,IB,O1,W)   = ',SURFACEWF_F(Z,UT,LUM,IB,O1,W)

END SUBROUTINE TP26G2

!******************************************************************************
!******************************************************************************
!                   Section: TP27* --> vlidort_l_thermalsup.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP28* --> not used
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP29* --> vlidort_masters_V2p8.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP30* --> vlidort_miscsetups.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP31* --> vlidort_PostProcessing.f90
!******************************************************************************
!******************************************************************************

SUBROUTINE TP31A1 (AWI,N,UM,LUM,K,O1,LCON,MCON,UHOM_UPDN,UHOM_UPUP,&
                   LUXR,MUXR,HMULT_1,HMULT_2,H_R)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)           :: AWI,N,UM,LUM,K,O1
DOUBLE PRECISION, intent(in)  :: LCON ( MAXSTRMSTKS, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: MCON ( MAXSTRMSTKS, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: UHOM_UPDN ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: UHOM_UPUP ( MAX_USER_STREAMS, MAXSTOKES, MAXEVALUES, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: LUXR,MUXR
DOUBLE PRECISION, intent(in)  :: HMULT_1 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: HMULT_2 ( MAXEVALUES, MAX_USER_STREAMS, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: H_R

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP31A1: IN SUBROUTINE WHOLELAYER_STERM_UP'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' LUM = ',LUM,' K = ',K,' O1 = ',O1
write(*,*) 'LCON(K,N)            = ',LCON(K,N)
write(*,*) 'MCON(K,N)            = ',MCON(K,N)
write(*,*) 'UHOM_UPDN(UM,O1,K,N) = ',UHOM_UPDN(UM,O1,K,N)
write(*,*) 'UHOM_UPUP(UM,O1,K,N) = ',UHOM_UPUP(UM,O1,K,N)
!write(*,*) 'LUXR                 = ',LUXR
!write(*,*) 'MUXR                 = ',MUXR
write(*,*) 'HMULT_1(K,UM,N)      = ',HMULT_1(K,UM,N)
write(*,*) 'HMULT_2(K,UM,N)      = ',HMULT_2(K,UM,N)
write(*,*) 'H_R                  = ',H_R

END SUBROUTINE TP31A1

!

SUBROUTINE TP31A2 (AWI,N,UM,O1,H_R,H_CR,LAYERSOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: AWI,N,UM,O1
DOUBLE PRECISION, intent(in)  :: H_R, H_CR
DOUBLE PRECISION, intent(in)  :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP31A2: IN SUBROUTINE WHOLELAYER_STERM_UP'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' O1 = ',O1
write(*,*) 'H_R                = ',H_R
write(*,*) 'H_CR               = ',H_CR
write(*,*) 'LAYERSOURCE(UM,O1) = ',LAYERSOURCE(UM,O1)

END SUBROUTINE TP31A2

!

SUBROUTINE TP31A3 (AWI,N,UM,LUM,IB,O1,UPAR_UP_2,EMULT_UP,LAYERSOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: AWI,N,UM,LUM,IB,O1
DOUBLE PRECISION, intent(in)  :: UPAR_UP_2 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: EMULT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
DOUBLE PRECISION, intent(in)  :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP31A3: IN SUBROUTINE WHOLELAYER_STERM_UP'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' LUM = ',LUM,' IB = ',IB,' O1 = ',O1
write(*,*) 'UPAR_UP_2(UM,O1,N) = ',UPAR_UP_2(UM,O1,N)
write(*,*) 'EMULT_UP(LUM,N,IB) = ',EMULT_UP(LUM,N,IB)
write(*,*) 'LAYERSOURCE(UM,O1) = ',LAYERSOURCE(UM,O1)

END SUBROUTINE TP31A3

!

SUBROUTINE TP31A4 (AWI,N,UM,K,SD,SU,ATERM_SAVE,BTERM_SAVE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: AWI,N,UM,K
DOUBLE PRECISION, intent(in)  :: SD,SU
DOUBLE PRECISION, intent(in)  :: ATERM_SAVE ( MAXEVALUES, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: BTERM_SAVE ( MAXEVALUES, MAXLAYERS )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP31A4: IN SUBROUTINE WHOLELAYER_STERM_UP'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' K = ',K
write(*,*) 'SD = ',SD
write(*,*) 'SU = ',SU
write(*,*) 'ATERM_SAVE(K,N) = ',ATERM_SAVE(K,N)
write(*,*) 'BTERM_SAVE(K,N) = ',BTERM_SAVE(K,N)

END SUBROUTINE TP31A4

!

SUBROUTINE TP31A5 (AWI,N,UM,LUM,IB,O1,UPAR_UP_1,EMULT_UP,LAYERSOURCE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: AWI,N,UM,LUM,IB,O1
DOUBLE PRECISION, intent(in)  :: UPAR_UP_1 ( MAX_USER_STREAMS, MAXSTOKES, MAXLAYERS )
DOUBLE PRECISION, intent(in)  :: EMULT_UP  ( MAX_USER_STREAMS, MAXLAYERS, MAXBEAMS )
DOUBLE PRECISION, intent(in)  :: LAYERSOURCE ( MAX_USER_STREAMS, MAXSTOKES )

INTEGER :: ATPI

IF (AWI /= ATPI()) RETURN

write(*,*)
write(*,*) 'TP31A5: IN SUBROUTINE WHOLELAYER_STERM_UP'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' LUM = ',LUM,' IB = ',IB,' O1 = ',O1
write(*,*) 'UPAR_UP_1(UM,O1,N) = ',UPAR_UP_1(UM,O1,N)
write(*,*) 'EMULT_UP(LUM,N,IB) = ',EMULT_UP(LUM,N,IB)
write(*,*) 'LAYERSOURCE(UM,O1) = ',LAYERSOURCE(UM,O1)

END SUBROUTINE TP31A5

!

SUBROUTINE TP31B1 (N,UM,AA,LCON,MCON,U_XNEG,U_XPOS,HMULT_1,HMULT_2)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: N,UM,AA
DOUBLE PRECISION, intent(in)  :: LCON(MAXSTREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: MCON(MAXSTREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: U_XPOS(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: U_XNEG(MAX_USER_STREAMS,MAXSTREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: HMULT_1(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: HMULT_2(MAXSTREAMS,MAX_USER_STREAMS,MAXLAYERS)

write(*,*)
write(*,*) 'TP31B1: IN SUBROUTINE WHOLELAYER_STERM_DN'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' AA = ',AA
write(*,*) 'LCON(AA,N)       = ',LCON(AA,N)
write(*,*) 'MCON(AA,N)       = ',MCON(AA,N)
write(*,*) 'U_XPOS(UM,AA,N)  = ',U_XPOS(UM,AA,N)
write(*,*) 'U_XNEG(UM,AA,N)  = ',U_XNEG(UM,AA,N)
write(*,*) 'HMULT_1(AA,UM,N) = ',HMULT_1(AA,UM,N)
write(*,*) 'HMULT_2(AA,UM,N) = ',HMULT_2(AA,UM,N)

END SUBROUTINE TP31B1

!

SUBROUTINE TP31B2 (N,UM,AA,SD,SU,ATERM_SAVE,BTERM_SAVE)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: N,UM,AA
DOUBLE PRECISION, intent(in)  :: SD,SU
DOUBLE PRECISION, intent(in)  :: ATERM_SAVE(MAXSTREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: BTERM_SAVE(MAXSTREAMS,MAXLAYERS)

write(*,*)
write(*,*) 'TP31B2: IN SUBROUTINE WHOLELAYER_STERM_DN'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' AA = ',AA
write(*,*) 'SD = ',SD
write(*,*) 'SU = ',SU
write(*,*) 'ATERM_SAVE(AA,N) = ',ATERM_SAVE(AA,N)
write(*,*) 'BTERM_SAVE(AA,N) = ',BTERM_SAVE(AA,N)

END SUBROUTINE TP31B2

!

SUBROUTINE TP31B3 (N,UM,LUM,IB,U_WNEG1,EMULT_DN)

USE VLIDORT_pars_m
IMPLICIT NONE

INTEGER, intent(in)    :: N,UM,LUM,IB
DOUBLE PRECISION, intent(in)  :: U_WNEG1(MAX_USER_STREAMS,MAXLAYERS)
DOUBLE PRECISION, intent(in)  :: EMULT_DN (MAX_USER_STREAMS,MAXLAYERS,MAXBEAMS)

write(*,*)
write(*,*) 'TP31B3: IN SUBROUTINE WHOLELAYER_STERM_DN'

write(*,*)
write(*,*) 'N = ',N,' UM = ',UM,' LUM = ',LUM,' IB = ',IB
write(*,*) 'U_WNEG1(UM,N)      = ',U_WNEG1(UM,N)
write(*,*) 'EMULT_DN(LUM,N,IB) = ',EMULT_DN(LUM,N,IB)

END SUBROUTINE TP31B3

!******************************************************************************
!******************************************************************************
!                   Section: TP32* --> vlidort_solutions.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP33* --> vlidort_Taylor.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP34* --> vlidort_thermalsup.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP35* --> vlidort_vfo_interface.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP36* --> vlidort_vfo_lcs_interface.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP37* --> vlidort_vfo_lps_interface.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP38* --> vlidort_transflux_MK3.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP71* --> vlidort_multipliers.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP72* --> vlidort_la_corrections.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP73* --> vlidort_lc_solutions.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP74* --> vlidort_lp_solutions.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP75* --> vlidort_lpc_bvproblem.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                   Section: TP76*  --> vlidort_converge.f90 
!******************************************************************************
!******************************************************************************


!******************************************************************************
!******************************************************************************
!                   Section: TP77*  --> vlidort_mediaprops.f90
!******************************************************************************
!******************************************************************************


!******************************************************************************
!******************************************************************************
!                   Section: TP78*  --> vlidort_lc_mediaprops.f90
!******************************************************************************
!******************************************************************************


!******************************************************************************
!******************************************************************************
!                   Section: TP79*  --> vlidort_lcs_converge.f90
!******************************************************************************
!******************************************************************************


!******************************************************************************
!******************************************************************************
!                   Section: TP80*  --> vlidort_lp_mediaprops.f90
!******************************************************************************
!******************************************************************************


!******************************************************************************
!******************************************************************************
!                   Section: TP81*  --> vlidort_lps_converge.f90
!******************************************************************************
!******************************************************************************



!******************************************************************************
!******************************************************************************
!                                END OF FILE
!******************************************************************************
!******************************************************************************
