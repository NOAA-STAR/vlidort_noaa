
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

module vfzmat_Pre_Master_m

!mick mod 9/19/2017 - use new numerical subroutines

   Use vfzmat_Numerical_m, only : GETQUAD2
   Use vfzmat_Rotation_m
   Use vfzmat_DevelopCoeffs_New_m, Only : vfzmat_DevelopGSFs

public

contains

subroutine vfzmat_Pre_Master &
   ( max_geoms, max_szas, max_vzas, max_azms, dtr,                            & ! Input Dimensions (VLIDORT)
     do_upwelling, do_dnwelling, do_ObsGeoms, ncoeffs, nstokes, n_QuadAngles, & ! Input Flags and Control
     n_geoms, n_szas, n_vzas, n_azms, offsets, szas, vzas, azms, obsgeoms,    & ! Input Geometries
     C1_up, S1_up, C2_up, S2_up, C1_dn, S1_dn, C2_dn, S2_dn,                  & ! Output rotation angles
     COSSCAT_up, COSSCAT_dn, QuadAngles, QuadCosines, QuadWeights,            & ! Output Scatcosines and Quadrature
     GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved )             ! Output Saved GSFs

   implicit none

!  precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions (geometry only)

   INTEGER, INTENT(IN)   :: max_geoms, max_szas, max_vzas, max_azms

!  Directional flags

   LOGICAL, INTENT(IN)   :: do_upwelling, do_dnwelling

!  Flags for use of observational geometry
 
   LOGICAL, INTENT(IN)   :: do_obsgeoms

!  Input Quadrature number

   INTEGER  , INTENT(IN) :: N_QuadAngles

!  numbers (general)

   INTEGER, INTENT(IN)   :: nstokes, ncoeffs

!  Geometry numbers

   INTEGER, INTENT(IN)   :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER, INTENT(IN)   :: offsets(max_szas,max_vzas)

!  Angles. Convention as for  VLIDORT

   REAL(dpk), INTENT(IN) :: dtr
   REAL(dpk), INTENT(IN) :: szas(max_szas)
   REAL(dpk), INTENT(IN) :: vzas(max_vzas)
   REAL(dpk), INTENT(IN) :: azms(max_azms)
   REAL(dpk), INTENT(IN) :: obsgeoms(max_geoms,3)

!  Output
!  ------

!  Output Quadrature

   REAL(dpk), INTENT(inout) :: QuadAngles  ( N_QuadAngles )
   REAL(dpk), INTENT(inout) :: QuadCosines ( N_QuadAngles )
   REAL(dpk), INTENT(inout) :: QuadWeights ( N_QuadAngles )

!  Spherical functions for quadrature angles

   REAL(dpk), INTENT (inout) :: GSF_P00_Saved  (N_QuadAngles,0:ncoeffs)
   REAL(dpk), INTENT (inout) :: GSF_P02_Saved  (N_QuadAngles,0:ncoeffs)
   REAL(dpk), INTENT (inout) :: GSF_P2p2_Saved (N_QuadAngles,0:ncoeffs)
   REAL(dpk), INTENT (inout) :: GSF_P2m2_Saved (N_QuadAngles,0:ncoeffs)

!  Scattering angle cosines

   REAL(dpk), INTENT(inout) :: COSSCAT_up(max_geoms)
   REAL(dpk), INTENT(inout) :: COSSCAT_dn(max_geoms)

!  rotational angles

   REAL(dpk), INTENT(inout) :: C1_up(max_geoms), S1_up(max_geoms)
   REAL(dpk), INTENT(inout) :: C2_up(max_geoms), S2_up(max_geoms)
   REAL(dpk), INTENT(inout) :: C1_dn(max_geoms), S1_dn(max_geoms)
   REAL(dpk), INTENT(inout) :: C2_dn(max_geoms), S2_dn(max_geoms)

!  local variables
!  ---------------

!  Parameters

   REAL(dpk), PARAMETER :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   REAL(dpk), PARAMETER :: d_half  = 0.5_dpk, d_two  = 2.0_dpk

!  Directional sign

   REAL(dpk) :: vsign

!  Help variables

   INTEGER   :: ia, ib, k, k1, k2, um, v
   REAL(dpk) :: ctheta, stheta, calpha, salpha, cphi, x1, x2, w1, w2

!  STAGE 1. Develop Geometry inputs for  F/Z matrix computation
!  ============================================================

!  1. Get the scattering angle cosines
!  -----------------------------------

   IF ( .not. Do_Obsgeoms ) THEN
     DO IB = 1, n_szas
       ctheta = cos ( szas(ib) * dtr )
       stheta = sin ( szas(ib) * dtr )
       DO UM = 1, n_vzas
         calpha = cos ( vzas(um) * dtr )
         salpha = sin ( vzas(um) * dtr )
         DO IA = 1, n_azms
           cphi = cos ( azms(ia) * dtr )
           V = OFFSETS(IB,UM) + IA
           COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
           COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
         ENDDO
       ENDDO
     ENDDO
   ELSE
     DO V = 1, n_geoms
       CTHETA = cos(obsgeoms(v,1)*dtr)
       STHETA = sin(obsgeoms(v,1)*dtr)
       CALPHA = cos(obsgeoms(v,2)*dtr)
       SALPHA = sin(obsgeoms(v,2)*dtr)
       CPHI   = cos(obsgeoms(v,3)*dtr)
       COSSCAT_up(v) = - CTHETA * CALPHA + STHETA * SALPHA * CPHI
       COSSCAT_dn(v) = + CTHETA * CALPHA + STHETA * SALPHA * CPHI
     ENDDO
   ENDIF

!  3. Get the Rotataions
!  ---------------------

!  upwelling.

   if ( do_upwelling ) then
      vsign = - d_one
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1_up, S1_up, C2_up, S2_up )
   endif

!  Downwelling

   if ( do_dnwelling ) then
      vsign = + d_one
      Call vfzmat_Rotation &
       ( max_geoms, max_szas, max_vzas, max_azms, vsign, dtr,   & ! Inputs
         do_ObsGeoms, nstokes, n_geoms, n_szas, n_vzas, n_azms, & ! inputs
         offsets, szas, vzas, azms, obsgeoms,                   & ! Inputs
         C1_dn, S1_dn, C2_dn, S2_dn )
   endif

!  STAGE 2. Develop quadrature and GSFs for Greekmat
!  =================================================

!  1. Quadrature
!  -------------

   Call GETQUAD2 ( -d_one, d_one, N_QuadAngles, QuadCosines, QuadWeights )
   do k1 = 1, N_QuadAngles / 2
      k2 = N_QuadAngles + 1 - k1 
      x1 = QuadCosines(k1) ; x2 = QuadCosines(k2)
      w1 = QuadWeights(k1) ; w2 = QuadWeights(k2)
      QuadCosines(k1) = x2 ; QuadCosines(k2) = x1
      QuadWeights(k1) = w2 ; QuadWeights(k2) = w1
   enddo
   do k = 1, N_QuadAngles
      QuadAngles(k) = acos ( QuadCosines(k) ) / dtr
   enddo

!  2. Develop generalized spherical functions for this quadrature
!  --------------------------------------------------------------

   Call  vfzmat_DevelopGSFs &
   ( ncoeffs, n_QuadAngles, Quadcosines, Quadweights,&
     GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved )

!  Done

   return
end subroutine vfzmat_Pre_Master

!  done module

end module vfzmat_Pre_Master_m

