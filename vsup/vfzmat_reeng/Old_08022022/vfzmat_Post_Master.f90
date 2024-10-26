
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

module vfzmat_Post_Master_m

!mick mod 9/19/2017 - use new numerical subroutines

   Use vfzmat_Numerical_m, only : BSpline, Seval
   Use vfzmat_PhasMat_m,   only : vfzmat_PhasMat1
   Use vfzmat_DevelopCoeffs_New_m, only : vfzmat_DevelopCoeffs_Fast

public

contains

subroutine vfzmat_Post_Master &
   ( max_geoms, N_InAngles, InCosines, InFmatrices,                & ! input  Fmatrices
     do_upwelling, do_dnwelling, Sunlight,                         & ! input  Flags
     ncoeffs, nstokes, n_geoms, n_QuadAngles, QuadCosines,         & ! Input numbers and quad
     COSSCAT_up, COSSCAT_dn, C1_up, S1_up, C2_up, S2_up, C1_dn, S1_dn, C2_dn, S2_dn, & ! Input angles
     GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved,                   & ! Input Saved GSFs     
     OutFmatrices_up, OutFmatrices_dn, Zmatrices_up, Zmatrices_dn, FMatCoeffs )        ! Output

!  Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16

!  Purpose
!  -------

!  Stand-alone routine to develop Z-matrices and F-Matrix expansion coefficients,
!    given only a set of F-Matrix inputs on a regular scattering-angle grid

!  Allows for various VLIDORT-based geometry options, upwelling and/or downwelling.

!  Stage 1. (a) Interpolate F-matrices to values implied by geometrical input.
!           (b) Transformation from scattering plane (Fmat) to meridional plane 
!               (Zmat) follows the rotations given in the VLIDORT code.

!  Stage 2. (a) Interpolate F-matrix input to quadrature grid for integration
!           (b) develop Coefficients from integrations using spherical-functions

!  Fmatrix input convention (J is the angle grid on input, N the layer index)
!  ------------------------

!       InFmatrices(N,J,1) = F11(J)
!       InFmatrices(N,J,2) = F22(J)
!       InFmatrices(N,J,3) = F33(J)
!       InFmatrices(N,J,4) = F44(J)
!       InFmatrices(N,J,5) = F12(J)
!       InFmatrices(N,J,6) = F34(J)

!   For Mie scattering, F11 = F22, F33 = F44

!  Convention for using FMatCoeffs output to get VLIDORT Greekmat input
!  --------------------------------------------------------------------

!       GREEKMAT(L,N,1)  --> use + FMatCoeffs(L,N,1)
!       GREEKMAT(L,N,2)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,5)  --> use - FMatCoeffs(L,N,5)
!       GREEKMAT(L,N,6)  --> use + FMatCoeffs(L,N,2)
!       GREEKMAT(L,N,11) --> use + FMatCoeffs(L,N,3)
!       GREEKMAT(L,N,12) --> use - FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,15) --> use + FMatCoeffs(L,N,6)
!       GREEKMAT(L,N,16) --> use + FMatCoeffs(L,N,4)

!    all other entries zero.

!  PATCH Upgrade, 08 November 2019
!  ===============================

!   BSPLINE/Seval routines must be done together inside m-loop
!   formerly, BSPLINE output was only for the m = 6 value.

   implicit none

!  precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER, INTENT(IN)   :: max_geoms

!  Directional flags and sunlight

   LOGICAL, INTENT(IN)   :: do_upwelling, do_dnwelling, Sunlight

!  Input F-matrix stuff ( angles and 6 scattering matrix entries )

   INTEGER  , INTENT(IN) :: N_InAngles
   REAL(dpk), INTENT(IN) :: InCosines   ( N_InAngles )
   REAL(dpk), INTENT(IN) :: InFmatrices ( N_InAngles, 6 )

!  numbers (general)

   INTEGER, INTENT(IN)   :: nstokes, ncoeffs, n_QuadAngles, n_geoms

!  Quadrature

   REAL(dpk), INTENT(in) :: QuadCosines ( N_QuadAngles )

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

!  Output
!  ------

!  Output Fmatrices (Interpolated)

   REAL(dpk), INTENT(INOUT) :: OutFmatrices_up   ( Max_Geoms, 6 )
   REAL(dpk), INTENT(INOUT) :: OutFmatrices_dn   ( Max_Geoms, 6 )

!  Zmatrices

   REAL(dpk), INTENT(INOUT) :: Zmatrices_up(max_geoms,4,4)
   REAL(dpk), INTENT(INOUT) :: Zmatrices_dn(max_geoms,4,4)

!  Fmatrix coefficients

   REAL(dpk), INTENT(INOUT) :: FMatCoeffs(0:ncoeffs,6)

!  local variables
!  ---------------

!  Parameters

   REAL(dpk), PARAMETER :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   REAL(dpk), PARAMETER :: d_half  = 0.5_dpk, d_two  = 2.0_dpk

!  Local Output Quad Fmatrices (Interpolated)

   REAL(dpk) :: QuadFmatrices ( N_QuadAngles, 6 )

!  Local Input Fmatrices

   REAL(dpk) :: Local_InFmatrices ( N_InAngles, 6 )

!  Help variables

   LOGICAL   :: mask(6)
   INTEGER   :: k, k1, m, v
   REAL(dpk) :: x, y
   REAL(dpk) :: bbas( N_InAngles ),cbas( N_InAngles ),dbas( N_InAngles )

!  STAGE 1. Develop Z-matrices
!  ===========================

!  . Spline-Interpolate F-matrices to VLIDORT Geometrical grid
!  ------------------------------------------------------------

!  Local F-matrix, reverse order

   do k = 1, N_InAngles
      k1 = N_InAngles + 1 - k 
      Local_InFmatrices(k1,1:6) = InFmatrices(k,1:6)
   enddo

!  Mask

   mask(1) = .true.
   if ( nstokes.eq.1) then
      mask(2:6) = .false.
   else
      mask(2:3) = .true. ;  mask(5) = .true.
      if ( nstokes.eq.4 ) then
         mask(4) = .true. ; mask(6) = .true.
      else
         mask(4) = .true. ; mask(6) = .true.
      endif
   endif

!  Start coefficient loop

   do m = 1, 6
      if (mask(m) ) then 

!    - Set the End-point gradient (YPN) = input gradient at forward-peak
!    - This make a HUGE difference to the accuracy
         !yp1 = d_zero
         !ypn = ( Local_InFmatrices(N_InAngles,m) - Local_InFmatrices(N_InAngles-1,m) ) / &
         !        ( InCosines(N_InAngles) -  InCosines(N_InAngles-1) )
         !Call vfzmat_SPLINE(N_InAngles,InCosines,Local_InFmatrices(:,m),N_InAngles,yp1,ypn,y2(:,m))
         !Call vfzmat_SPLINT(N_InAngles,InCosines,Local_InFmatrices(:,m),y2(:,m),N_InAngles,x,y)

         Call BSpline (N_InAngles,N_InAngles,InCosines,Local_InFmatrices(1:n_InAngles,m),bbas,cbas,dbas)

!  Interpolate upwelling

         if ( do_upwelling ) then
           do v = 1, n_geoms
             x = cosscat_up(v)
             Call Seval (N_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(1:n_InAngles,m),bbas,cbas,dbas,y)
             OutFmatrices_up(v,m) = y
           enddo
         endif

!  Interpolate downwelling

         if ( do_dnwelling ) then
           do v = 1, n_geoms
             x = cosscat_dn(v)
             Call Seval (N_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(1:n_InAngles,m),bbas,cbas,dbas,y)
             OutFmatrices_dn(v,m) = y
           enddo
         endif

!  End coefficient type loop

      endif
   enddo

!  3. Get the Z-matrices
!  ---------------------

!  upwelling. [ C1/S1/C2/S2 are local, will be overwritten ]

   if ( do_upwelling ) then
      Call vfzmat_PhasMat1 &
       ( max_geoms, nstokes, n_geoms, Sunlight,       & ! inputs
         C1_up, S1_up, C2_up, S2_up, OutFmatrices_up, & ! Inputs
         Zmatrices_up )
   endif

!  Downwelling

   if ( do_dnwelling ) then
      Call vfzmat_PhasMat1 &
       ( max_geoms, nstokes, n_geoms, Sunlight,        & ! inputs
         C1_dn, S1_dn, C2_dn, S2_dn, OutFmatrices_dn,  & ! Inputs
         Zmatrices_dn )
   endif

!  STAGE 2. Develop Coefficients
!  =============================

!  1. Spline-Interpolate F-matrices to Quadrature grid for coefficients
!  --------------------------------------------------------------------

!  Start coefficient loop

   do m = 1, 6
      if (mask(m) ) then 
         Call BSpline (N_InAngles,N_InAngles,InCosines,Local_InFmatrices(1:n_InAngles,m),bbas,cbas,dbas)
         do k = 1, N_QuadAngles
            k1 = N_QuadAngles + 1 - k 
            x = QuadCosines(k1)
            Call Seval (N_InAngles,N_InAngles,x,InCosines,Local_InFmatrices(1:n_InAngles,m),bbas,cbas,dbas,y)
            QuadFmatrices(k1,m) = y
         enddo
      endif
   enddo

!  Get the coefficients from the DevelopCoeffs routine.

   Call vfzmat_DevelopCoeffs_Fast &
     ( ncoeffs, N_QuadAngles, nstokes, QuadFmatrices,               &
      GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved, &
      Fmatcoeffs )

!  Done

   return
end subroutine vfzmat_Post_Master

!  done module

end module vfzmat_Post_Master_m

