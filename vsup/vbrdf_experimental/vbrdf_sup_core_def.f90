
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

!  1/31/21, Version 2.8.3. Restrict parameters, otherwise no changes

!  2/25/22. Version 2.8.5. Increased flexibility.
!    -- Scalar and vector amplitudes given separately, controlled with flag
!    -- remove calculation of USER_BRDF_F_0. Flag for this.
!    -- reduced kernel for SNOW BRDF

!  2/28/22. Version 2.8.5. Semi-analytical Snow BRDF kernel
!     -- RTS model (Ross-Thick-Snow) of A. Ding et al., Remote Sensing, 11, 1611 (2019).
!     -- Basic model is the ART (Asymptotic radiative transfer) model of Kokhanovsky and Breon
!     -- Use as replacement of the Roujean kernel in the usual MODIS 3-kernel configuration
!     -- Must be used in conjunction with Lambertian and Ross-Thick kernels
!     -- One free parameter

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!     -- NEW. Two new type structures ("Core" and "Book")
!     -- "Core" contains kernel cores and helper routines. Wavelength independent
!     -- "Book" contains pre-calculated useul bookkeeping quantities (wavelength independent)

module VBRDF_Sup_Core_def_m

!  This module contains the following structures:

!  VBRDF_Sup_Core - Intent(InOut) for VBRDFSup. These are Core VARIABLES

      use vlidort_pars_m, Only : fpk, MAXSTREAMS, MAXMOMENTS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF, &
                                 MAXBEAMS, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAX_BRDF_KERNELS

      implicit none

! #####################################################################
! #####################################################################

      type VBRDF_Sup_Core

!  Output BRDF CORE functions
!  ==========================

!  at quadrature (discrete ordinate) angles

      DOUBLE PRECISION :: BRDFUNC_CORE   ( MAX_BRDF_KERNELS,    MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_CORE ( MAX_BRDF_KERNELS,    MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_HELP   ( MAX_BRDF_KERNELS, 7, MAXSTREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: BRDFUNC_0_HELP ( MAX_BRDF_KERNELS, 7, MAXSTREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  at user-defined stream directions

      DOUBLE PRECISION :: USER_BRDFUNC_CORE   ( MAX_BRDF_KERNELS,    MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_CORE ( MAX_BRDF_KERNELS,    MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_HELP   ( MAX_BRDF_KERNELS, 7, MAX_USER_STREAMS, MAXSTREAMS, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_BRDFUNC_0_HELP ( MAX_BRDF_KERNELS, 7, MAX_USER_STREAMS, MAXBEAMS,   MAXSTREAMS_BRDF )

!  Exact DB values

      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_CORE ( MAX_BRDF_KERNELS,    MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )
      DOUBLE PRECISION :: DBKERNEL_BRDFUNC_HELP ( MAX_BRDF_KERNELS, 7, MAX_USER_STREAMS, MAX_USER_RELAZMS, MAXBEAMS )

!  Values for Emissivity
!   -- 7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: EBRDFUNC_CORE ( MAX_BRDF_KERNELS,   MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: EBRDFUNC_HELP ( MAX_BRDF_KERNELS, 7, MAXSTREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC_CORE ( MAX_BRDF_KERNELS,   MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: USER_EBRDFUNC_HELP ( MAX_BRDF_KERNELS, 7, MAX_USER_STREAMS, MAXSTREAMS_BRDF, MAXSTREAMS_BRDF )

!  Output for WSA/BSA scaling options. New, Version 2.7

      DOUBLE PRECISION :: SCALING_BRDFUNC_CORE   ( MAX_BRDF_KERNELS, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0_CORE ( MAX_BRDF_KERNELS, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_HELP   ( MAX_BRDF_KERNELS, 7, MAXSTREAMS_SCALING, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )
      DOUBLE PRECISION :: SCALING_BRDFUNC_0_HELP ( MAX_BRDF_KERNELS, 7, MAXSTREAMS_SCALING, MAXSTREAMS_BRDF )

      end type VBRDF_Sup_Core

! #####################################################################
! #####################################################################

      type VBRDF_Sup_Book

!  Output BRDF Bookkeeping
!  =======================

!  Masking

      INTEGER            :: QMASK(16)

!  Bookkeeping on the number of parameters

      INTEGER            :: N_BRDF_PARAMETERS(MAX_BRDF_KERNELS)

!  White-sky and Black-sky albedo scaling

      LOGICAL            :: DO_LOCAL_WSA, DO_LOCAL_BSA
      DOUBLE PRECISION   :: SCALING_QUAD_STRMWTS(MAXSTREAMS_SCALING)

!  Fresnel scaling information

      LOGICAL          :: DO_FRESNEL_SCALING(MAX_BRDF_KERNELS)
      LOGICAL          :: DO_FRESNEL_SCALAR (MAX_BRDF_KERNELS)

!  BRDF azimuth quadrature streams
!  7/28/21. 8/18/21. Half-range azimuth integration variable introduced. This is the default.
!   -- if in operation, azimuthal range is [0,pi], if not the range is [pi,pi]
!   -- first set as a hard-wired parameter, now determined according to type of kernel
!         (TRUE for all kernels with facet-isotropy, False for NewCM/NewGCM with facet anisotropy

      LOGICAL ::          DO_HALF_RANGE
      INTEGER ::          NSTREAMS_BRDF, NBRDF_HALF
      DOUBLE PRECISION :: A_BRDF  ( MAXSTREAMS_BRDF )

!  BRDF azimuth quadrature streams For emission calculations
!  7/28/21. Necessary now to use MAXSTREAMS_BRDF here instead of MAXSTHALF_BRDF

      DOUBLE PRECISION :: BAX_BRDF ( MAXSTREAMS_BRDF )

!  Azimuth factors

      LOGICAL          :: ADD_FOURIER(0:MAXMOMENTS,MAX_BRDF_KERNELS)
      DOUBLE PRECISION :: DELFAC     (0:MAXMOMENTS)
      DOUBLE PRECISION :: BRDF_COSAZMFAC(0:MAXMOMENTS,MAXSTREAMS_BRDF)
      DOUBLE PRECISION :: BRDF_SINAZMFAC(0:MAXMOMENTS,MAXSTREAMS_BRDF)

      end type VBRDF_Sup_Book

! #####################################################################
! #####################################################################

!  EVERYTHING PUBLIC HERE

   PRIVATE
   PUBLIC :: VBRDF_Sup_Core, VBRDF_Sup_Book

end module VBRDF_Sup_Core_def_m
