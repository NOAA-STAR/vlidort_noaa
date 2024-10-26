
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

module vfzmat_DevelopCoeffs_New_m

public

contains

subroutine vfzmat_DevelopGSFs &
   ( ncoeffs, nangles, cosines, weights,&
     GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved )

!  Stand-alone routine to develop GSFs on Quadrature grid
!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!  First Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16
!   This routine adapted, 02/11/22

!  ************************************************************************
!  *  Calculate the generalized spherical functions 
!  ************************************************************************

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER          , INTENT (IN) :: ncoeffs, nangles
 
   REAL    (KIND=dpk), INTENT (IN) :: cosines(nangles)
   REAL    (KIND=dpk), INTENT (IN) :: weights(nangles)

!  local variables

   REAL    (KIND=dpk), INTENT (INOUT) :: GSF_P00_Saved  (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (INOUT) :: GSF_P02_Saved  (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (INOUT) :: GSF_P2p2_Saved (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (INOUT) :: GSF_P2m2_Saved (nangles,0:ncoeffs)

!  Local

   REAL    (KIND=dpk) :: P00(nangles,2)
   REAL    (KIND=dpk) :: P02(nangles,2)
   REAL    (KIND=dpk) :: P2p2(nangles,2)
   REAL    (KIND=dpk) :: P2m2(nangles,2)

   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   real(dpk), parameter :: d_half  = 0.5_dpk, d_two  = 2.0_dpk
   real(dpk), parameter :: d_three = 3.0_dpk, d_four = 4.0_dpk

   INTEGER            :: i, l, lnew, lold, itmp
   REAL    (KIND=dpk) :: dl, dl2, qroot6, fac1, fac2, fac3, fl, &
                         sql4, sql41, twol1, tmp1, tmp2, denom, help

!  Initialization

  GSF_P00_Saved  = d_zero
  GSF_P02_Saved  = d_zero
  GSF_P2p2_Saved = d_zero
  GSF_P2m2_Saved = d_zero

  qroot6 = -0.25_dpk*SQRT(6.0_dpk)

!  first update generalized spherical functions, then calculate coefs. *
!  lold and lnew are pointer-like indices used in recurrence           *

  lnew = 1
  lold = 2

  DO l = 0, ncoeffs

    IF (l == 0) THEN

      dl   = d_zero
      DO  i=1, nangles
        P00(i,lold) = d_one
        P00(i,lnew) = d_zero
        P02(i,lold) = d_zero
        P2p2(i,lold) = d_zero
        P2m2(i,lold)= d_zero
        P02(i,lnew) = d_zero
        P2p2(i,lnew) = d_zero
        P2m2(i,lnew)= d_zero
      END DO

    ELSE

      dl   = DBLE(l)
      dl2  = dl * dl
      fac1 = (d_two*dl-d_one)/dl
      fac2 = (dl-d_one)/dl
      DO  i=1, nangles
        P00(i,lold) = fac1*cosines(i)*P00(i,lnew) - fac2*P00(i,lold)
      END DO

    ENDIF

    IF (l == 2) THEN

      DO  i=1, nangles
        P02(i,lold) = qroot6*(d_one-cosines(i)*cosines(i))
        P2p2(i,lold) = 0.25_dpk*(d_one+cosines(i))*(d_one+cosines(i))
        P2m2(i,lold)= 0.25_dpk*(d_one-cosines(i))*(d_one-cosines(i))
        P02(i,lnew) = d_zero
        P2p2(i,lnew) = d_zero
        P2m2(i,lnew) = d_zero
      END DO
      sql41 = d_zero

    ELSE IF (l > 2) THEN

      sql4  = sql41
      sql41 = dsqrt(dl2-d_four)
      twol1 = d_two*dl - d_one
      tmp1  = twol1/sql41
      tmp2  = sql4/sql41
      denom = (dl-d_one)*(dl2-d_four)
      fac1  = twol1*(dl-d_one)*dl/denom
      fac2  = d_four*twol1/denom
      fac3  = dl*((dl-d_one)*(dl-d_one)-d_four)/denom
      DO i=1, nangles
        P02 (i,lold) = tmp1*cosines(i)*P02(i,lnew)         - tmp2*P02(i,lold)
        P2p2(i,lold) = (fac1*cosines(i)-fac2)*P2p2(i,lnew) - fac3*P2p2(i,lold)
        P2m2(i,lold) = (fac1*cosines(i)+fac2)*P2m2(i,lnew) - fac3*P2m2(i,lold)
      END DO

    END IF

    itmp = lnew
    lnew = lold
    lold = itmp

    fl = dl+d_half
    do i=1, nangles
      help =  weights(i) * fl
      GSF_P00_saved (i,L) = help * P00 (i,lnew)
      GSF_P02_saved (i,L) = help * P02 (i,lnew)
      GSF_P2p2_saved(i,L) = help * P2p2(i,lnew) * d_half
      GSF_P2m2_saved(i,L) = help * P2m2(i,lnew) * d_half
    END DO

!  End coefficient loop

  END DO

!  finish

   return
end subroutine vfzmat_DevelopGSFs


subroutine vfzmat_DevelopCoeffs_Fast &
   ( ncoeffs, nangles, nstokes, Fmat,      &
     GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved, &
     Expcoeffs )

!  Stand-alone routine to develop coefficients from Scattering Matrix on Quadrature grid
!  Based on the Meerhoff Mie code (as found in RTSMie package), and adapted

!  First Programmed 03 february 2016 by R. Spurr, RT Solutions Inc.
!   "vfzmat" Supplement for VLIDORT, arranged 9/19/16
!   This routine adapted, 02/11/22. Uses pre-calculated GSFs.

!  ************************************************************************
!  *  Calculate the expansion coefficients of the scattering matrix in    *
!  *  generalized spherical functions by numerical integration over the   *
!  *  scattering angle.                                                   *
!  ************************************************************************

   implicit none

!  precision

   integer, parameter :: dpk = SELECTED_REAL_KIND(15)

!  input

   INTEGER          , INTENT (IN) :: ncoeffs, nangles, nstokes
 
!  Fmatrices

   REAL    (KIND=dpk), INTENT (IN) :: FMAT(nangles,6)

!  Modified GSFs

   REAL    (KIND=dpk), INTENT (IN) :: GSF_P00_Saved  (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (IN) :: GSF_P02_Saved  (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (IN) :: GSF_P2p2_Saved (nangles,0:ncoeffs)
   REAL    (KIND=dpk), INTENT (IN) :: GSF_P2m2_Saved (nangles,0:ncoeffs)

!  output (Initialized)

   REAL    (KIND=dpk), INTENT (OUT) :: expcoeffs(0:ncoeffs,6)

!  local variables


   real(dpk), parameter :: d_zero  = 0.0_dpk, d_one  = 1.0_dpk
   INTEGER              :: l, index_11, index_12, index_22, index_33, index_34, index_44 
   REAL    (KIND=dpk)   :: alfap, alfam, fmat_plus(nangles), fmat_minus(nangles)

!  Initialization

   expcoeffs = d_zero

!  New indexing consistent with Tmatrix output

  index_11 = 1
  index_22 = 2
  index_33 = 3
  index_44 = 4
  index_12 = 5
  index_34 = 6

  DO l = 0, ncoeffs
     expcoeffs(L,index_11) = dot_Product ( GSF_P00_Saved(1:nangles,L),fmat(1:nangles,1) )
  enddo

  if ( nstokes .gt. 1 ) then
    fmat_plus (1:nangles) = fmat(1:nangles,2)+fmat(1:nangles,3)
    fmat_minus(1:nangles) = fmat(1:nangles,2)-fmat(1:nangles,3)
    DO l = 0, ncoeffs
      alfap = dot_Product ( GSF_P2p2_Saved(1:nangles,L),fmat_plus (1:nangles) )
      alfam = dot_Product ( GSF_P2m2_Saved(1:nangles,L),fmat_minus(1:nangles) )
      expcoeffs(L,index_22) =  alfap + alfam
      expcoeffs(L,index_33) =  alfap - alfam
      expcoeffs(L,index_12) = dot_Product ( GSF_P02_Saved(1:nangles,L),fmat(1:nangles,5) )
    ENDDO
  endif

  if ( nstokes.eq.4 ) then
    DO l = 0, ncoeffs
      expcoeffs(L,index_44) =   dot_Product ( GSF_P00_Saved(1:nangles,L),fmat(1:nangles,4) )
      expcoeffs(L,index_34) = - dot_Product ( GSF_P02_Saved(1:nangles,L),fmat(1:nangles,6) )
    ENDDO
  endif

!  Phase function normalization

  expcoeffs(0,index_11)        = d_one

!  finish

  return
end subroutine vfzmat_DevelopCoeffs_Fast

!  End module

End Module vfzmat_DevelopCoeffs_New_m

