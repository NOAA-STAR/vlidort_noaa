Module GEMSTOOL_aerload_routines_m

!  Modified slightly for the PCA Software. 9/5/17.

!  Modules in GEMSTOOL_sourcecode/structures

   USE GEMSTOOL_pars_m, Only : GTPK, GTZERO, GTONE

private :: exp1_profile_only, exp1_profile_plus, &
           gdf1_profile_only, gdf1_profile_plus

public :: profiles_uniform,&
          profiles_expone, &
          profiles_gdfone

contains

subroutine profiles_uniform                                 &
          ( maxlayers, nlayers, heightgrid, do_derivatives, & ! Inputs
            z_upperlimit, z_lowerlimit, total_column,       & ! Inputs
            profile, d_profile_dcol,                        & ! output
            fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) :: maxlayers
      integer, intent(in) :: nlayers
      real(GTPK), intent(in) ::  heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of Uniform slab
!     z_lowerlimit = lowest    height value at bottom of Uniform slab

!  Notes:
!    (1) Must have  z_u >/= z_l. This is checked internally
!    (2) Must have  z_l >/= heightgrid(nlayers). Also checked internally

      real(GTPK),    intent(in) ::  z_upperlimit, z_lowerlimit

!  Total column between upper and lower limits

      real(GTPK),    intent(in) ::  total_column

!  Control for calculating derivatives

      logical, intent(in) :: do_derivatives

!  Output
!  ======

!  Umkehr profile

      real(GTPK),    intent(out) ::  profile ( maxlayers )

!  Umkehr profile derivatives  (w.r.t. total column )

      real(GTPK),    intent(out) :: d_profile_dcol ( maxlayers )

!  Exception handling

      logical      , intent(out) ::  fail
      character*(*), intent(out) ::  message, action

!  Local variables
!  ===============

!  help variables

      integer          :: n, nupper, nlower
      real(GTPK)       :: z1,z2,density, scaling, za,zb

!  Initial checks
!  ==============

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        message = 'Bad input: Uniform slab Lower limit is below ground'
        action  = ' increase z_lowerlimit > BOA, heightgrid(nlayers)'
        fail = .true. ; return
      endif

      if ( z_upperlimit .gt. heightgrid(0)  ) then
        message = 'Bad input: Uniform slab upperlimit is above TOA'
        action  = ' decrease z_upperlimit < TOA, heightgrid(0)'
        fail = .true. ; return
      endif

      if ( z_lowerlimit .gt. z_upperlimit ) then
        message = 'Bad input:  Uniform slab lower > upper limit'
        action  = ' increase z_upperlimit > z_Lowerlimit'
        fail = .true. ; return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output 

      do n = 1, nlayers
        profile(n) = GTZERO
      enddo
      if ( do_derivatives ) then
        do n = 1, nlayers
          d_profile_dcol(n) = GTZERO
        enddo
      endif

!  Basic settings

      z1 = z_upperlimit
      z2 = z_lowerlimit
      density = total_column / ( z1 - z2 )
 
!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit

      n = nupper
      za = z1
      if ( nlower .eq. nupper ) then
        scaling = GTONE
      else
        zb = heightgrid(nupper)
        scaling = ( z1 - zb ) / ( z1 - z2 )
      endif
      profile(n) = total_column * scaling
      if ( do_derivatives) d_profile_dcol(n) = scaling
 
!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        scaling = ( za - zb ) / ( z1 - z2 )
        profile(n)        = total_column * scaling
        if ( do_derivatives) d_profile_dcol(n) = scaling
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        scaling = ( za - zb ) / ( z1 - z2 )
        profile(n)        = total_column * scaling
        if ( do_derivatives) d_profile_dcol(n) = scaling
      endif

!  Debug Check integrated profile
!  ------------------------------

!      sum = GTZERO
!      sumd = GTZERO
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo
!      write(*,*)sum, total_column, sumd

!  Finish

      return
end subroutine profiles_uniform

!

subroutine profiles_expone                                         &
          ( maxlayers, nlayers, heightgrid, do_derivatives,        & ! Inputs
            z_upperlimit, z_lowerlimit, relaxation, total_column,  & ! Inputs
            profile, d_profile_dcol, d_profile_drxn,               & ! Output
            fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) :: maxlayers
      integer, intent(in) :: nlayers
      real(GTPK),    intent(in) ::  heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of EXP profile
!     z_lowerlimit = lowest    height value at bottom of EXP profile

!  Notes:
!    (1) Must have  z_u >/= z_l. This is checked internally
!    (2) Must have  z_l >/= heightgrid(nlayers). Also checked internally

      real(GTPK),    intent(in) ::  z_upperlimit, z_lowerlimit

!  Relaxation constant for the Exponential function

      real(GTPK),    intent(in) ::  Relaxation

!  Total column between upper and lower limits

      real(GTPK),    intent(in) ::  total_column

!  Control for calculating derivatives

      logical, intent(in) :: do_derivatives

!  Output
!  ======

!  Umkehr profile 

      real(GTPK),    intent(out) ::  profile ( maxlayers )

!  Umkehr profile derivatives
!    (w.r.t. total column and relaxation parameter)

      real(GTPK),    intent(out) :: d_profile_dcol ( maxlayers )
      real(GTPK),    intent(out) :: d_profile_drxn ( maxlayers )

!  Exception handling

      logical      , intent(out) ::  fail
      character*(*), intent(out) ::  message, action

!  Local variables
!  ===============

!  help variables

      integer          :: n, nupper, nlower
      real(GTPK)    :: k,b,ed,zd,kd,za,zb,term,omega,z1,z2,db_dk

!  Debug
!      logical          do_debug_profile
!      integer          m,msav,mdebug
!      parameter        ( mdebug = 1001 )
!      real(GTPK)    pp(mdebug),zp(mdebug), sum, sumd, z, q

!  Initial checks
!  ==============

!  set debug output

!      do_debug_profile = .false.

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is below ground'
        action  = ' increase z_lowerlimit > heightgrid(nlayers)'
        return
      endif

      if ( z_lowerlimit .gt. z_upperlimit ) then
        fail = .true.
        message = 'Bad input: GDF transition height above upper limit'
        action  = ' increase z_upperlimit > z_lowerlimit'
        return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output

      do n = 1, nlayers
        profile(n) = GTZERO
        if ( do_derivatives ) then
          d_profile_dcol(n) = GTZERO
          d_profile_drxn(n) = GTZERO
        endif
      enddo

!  set profile parameters

      k  = relaxation
      z1 = z_upperlimit
      z2 = z_lowerlimit
      omega = total_column
      zd = z1 - z2
      kd = k * zd
      ed = dexp ( - kd)
      term = kd - GTONE + ed
      b  = omega * k / term

!  Set profile derivative parameters

      db_dk = ( omega - b * zd * (1.0 - ed ) ) / term
 
!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit. Partly EXP

      n = nupper
      za = z1
      zb = heightgrid(nupper)
      if ( do_derivatives ) then
        call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
              profile(n), d_profile_dcol(n), d_profile_drxn(n))
      else
        call exp1_profile_only(za,zb,z1,b,k,profile(n))
      endif

!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        if ( do_derivatives ) then
          call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
              profile(n), d_profile_dcol(n), d_profile_drxn(n))
        else
          call exp1_profile_only(za,zb,z1,b,k,profile(n))
        endif
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        if ( do_derivatives ) then
          call exp1_profile_plus (za,zb,z1,b,k,omega,db_dk, &
             profile(n), d_profile_dcol(n), d_profile_drxn(n))
        else
          call exp1_profile_only(za,zb,z1,b,k,profile(n))
        endif
      endif

!  debug Check integrated profile
!  ------------------------------

!      sum = GTZERO
!      sumd = GTZERO
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo

!  Debug profile itself
!  --------------------

!      if ( do_debug_profile ) then
!        m = 0
!        z = z1
!        do while (z.ge.z2)
!          z = z - 0.01
!          m = m + 1
!          if(m.gt.mdebug)stop'debug dimension needs to be increased'
!          zp(m) = z
!          q = dexp ( - k * (z1 - z ) )
!          pp(m) = b * (z1-z) + (b*q/k)
!        enddo
!        msav = m
!        do m = 1, msav
!          write(81,'(13f10.4)')zp(m),pp(m)
!        enddo
!      endif

!  Finish

      return
end subroutine profiles_expone

!

subroutine exp1_profile_only(za,zb,z1,b,k,prof)

      implicit none

!  Profile assignation

!  Arguments

      real(GTPK),    intent(in)  :: za,zb,z1,b,k
      real(GTPK),    intent(out) :: prof

!  Local variables

      real(GTPK)    :: da,db,bok,ea,eb

!  Code

      prof = GTZERO
      da = z1 - za
      db = z1 - zb
      bok = b / k
      ea = dexp(-k*da)
      eb = dexp(-k*db)
      prof = b * ( za - zb ) - bok * ( ea - eb )

!  Finish

      return
end subroutine exp1_profile_only

!

subroutine exp1_profile_plus ( za,zb,z1,b,k,omega,db_dk, &
                               prof,deriv1,deriv2)

      implicit none

!  Profile and derivative assignations

!  Arguments

      real(GTPK),    intent(in)  :: za,zb,z1,b,k,omega,db_dk
      real(GTPK),    intent(out) :: prof,deriv1,deriv2

!  Local variables

      real(GTPK)    :: da,db,ea,eb,bok,term1,term2,term3

!  initialise

      prof   = GTZERO
      deriv1 = GTZERO
      deriv2 = GTZERO

!  Code for profile

      da = z1 - za
      db = z1 - zb
      bok = b / k
      ea = dexp(-k*da)
      eb = dexp(-k*db)
      prof = b * ( za - zb ) - bok * ( ea - eb )

!  Derivatives

      deriv1 = prof / omega
      term1  = prof * db_dk / b
      term2  = (ea - eb ) * bok / k
      term3  = bok * ( da*ea - db*eb )
      deriv2 = term1 + term2 + term3

!  Finish

      return
end subroutine exp1_profile_plus

!

subroutine profiles_gdfone                                                     &
         ( maxlayers, nlayers, heightgrid, do_derivatives,                     &
           z_upperlimit, z_peakheight, z_lowerlimit, half_width, total_column, &
           profile, d_profile_dcol, d_profile_dpkh, d_profile_dhfw,            &
           fail, message, action )

      implicit none

!  Inputs
!  ======

!  Height grid in [km]

      integer, intent(in) ::  maxlayers
      integer, intent(in) ::  nlayers
      real(GTPK),    intent(in) :: heightgrid(0:maxlayers)

!  height levels controlling profile shape. All in [km].

!     z_upperlimit = uppermost height value at start  of GDF profile
!     z_lowerlimit = lowest    height value at bottom of GDF profile
!     z_peakheight = height at which peak value of GDF profile occurs

!  Notes:
!    (1) Must have  z_u > z_p >/= z_l. This is checked internally
!    (2) Must have  z_l > heightgrid(nlayers). Also checked internally

      real(GTPK),    intent(in) :: z_upperlimit, z_lowerlimit,  z_peakheight

!  Half width

      real(GTPK),    intent(in) :: half_width

!  Total column between upper and lower limits

      real(GTPK),    intent(in) :: total_column

!  Control for calculating derivatives

      logical, intent(in)   :: do_derivatives

!  Output
!  ======

!  Umkehr profile 

      real(GTPK),    intent(out) :: profile ( maxlayers )

!  Umkehr profile derivatives
!    (w.r.t. total column and peak height and half-width)

      real(GTPK),    intent(out)  :: d_profile_dcol ( maxlayers )
      real(GTPK),    intent(out)  :: d_profile_dpkh ( maxlayers )
      real(GTPK),    intent(out)  :: d_profile_dhfw ( maxlayers )

!  Exception handling

      logical      , intent(out)     :: fail
      character*(*), intent(out)     :: message, action

!  Local variables
!  ===============

!  help variables
     
      real(GTPK)    :: omega,z1,z0,z2,h
      integer      :: n, nupper, nlower
      real(GTPK)    :: a,zt1,zt2,za,zb
      real(GTPK)    :: a_d,a_p,a_w,term,q1,q2,r1,r2,rsq1,rsq2

!  Debug

!      logical          do_debug_profile
!      integer          m,msav,mdebug
!      parameter        ( mdebug = 1001 )
!      real(GTPK)    pp(mdebug),zp(mdebug), sum, sumd, z, q, r
     
!  Initial checks
!  ==============

!  set debug output

!      do_debug_profile = .true.

!  Initialise exception handling

      fail = .false.
      message = ' '
      action  = ' '

!  Check input

      if ( z_lowerlimit .lt. heightgrid(nlayers) ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is below ground'
        print*,message
        action  = ' increase z_lowerlimit > heightgrid(nlayers)'
        return
      endif

      if ( z_lowerlimit .gt. z_peakheight ) then
        fail = .true.
        message = 'Bad input: GDF Lower limit is above peak height'
        print*,message 
        action  = ' increase z_peakheight > z_lowerlimit'
        return
      endif

      if ( z_peakheight .gt. z_upperlimit ) then
        fail = .true.
        message = 'Bad input: GDF transition height above upper limit'
        print*,message
        action  = ' increase z_upperlimit > z_peakheight'
        return
      endif

!  Find intercept layers
!  ---------------------

      n = 0
      do while ( heightgrid(n).ge.z_upperlimit  )
        n = n + 1
      enddo
      nupper = n

      do while (heightgrid(n).gt.z_lowerlimit)
        n = n + 1
      enddo
      nlower = n

!  Zero output 

      profile = GTZERO
      if ( do_derivatives ) then
        d_profile_dcol = GTZERO
        d_profile_dpkh = GTZERO
        d_profile_dhfw = GTZERO
      endif

!  Basic settings

      h  = half_width
      z1 = z_upperlimit
      z0 = z_peakheight
      z2 = z_lowerlimit
      omega = total_column

!  set profile parameters

      zt1 = z1 - z0
      zt2 = z0 - z2
      q1 = dexp ( - h * zt1)
      q2 = dexp ( - h * zt2)
      r1 = GTONE / ( GTONE + q1 )
      r2 = GTONE / ( GTONE + q2 )
      term = ( r1 + r2 - GTONE )
      a  = omega * h / term

!  Set profile derivative parameters

      if ( do_derivatives ) then
        rsq1 = r1 * r1 * q1
        rsq2 = r2 * r2 * q2
        a_d = a / omega
        a_p  = - a * h * ( rsq2 - rsq1 ) / term
        a_w = ( omega - a * ( rsq2*zt2 + rsq1*zt1 ) ) / term
       endif

!  construct profile Umkehrs
!  -------------------------

!  Layer containing z_upperlimit. Partly EXP

      n = nupper
      za = z1
      zb = heightgrid(nupper)
      if ( do_derivatives ) then
        call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                           & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
      else
        call gdf1_profile_only(za,zb,z0,h,a,profile(n))
      endif

!  Intermediate Layers

      do n = nupper+1,nlower-1
        za = heightgrid(n-1)
        zb = heightgrid(n)
        if ( do_derivatives ) then
          call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                         & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
        else
          call gdf1_profile_only(za,zb,z0,h,a,profile(n))
        endif
      enddo

!  Layer containing z_lowerlimit
!   ( Will not be done if nupper = nlower)

      if ( nupper .lt. nlower ) then
        n = nlower
        za = heightgrid(n-1)
        zb = z2
        if ( do_derivatives ) then
          call gdf1_profile_plus (za,zb,z0,h,a,a_d,a_p,a_w,                         & ! Inputs
              profile(n), d_profile_dcol(n),  d_profile_dpkh(n),d_profile_dhfw(n))    ! Output
        else
          call gdf1_profile_only(za,zb,z0,h,a,profile(n))
        endif
      endif

!  Debug Check integrated profile
!  ------------------------------

!      sum = GTZERO
!      sumd = GTZERO
!      do n = nupper, nlower
!        sum = sum + profile(n)
!        if(do_derivatives)sumd = sumd + d_profile_dcol(n)
!      enddo
!      write(*,*)sum, omega, sumd

!  Debug profile itself
!  --------------------

!      if ( do_debug_profile ) then
!        m = 0
!        z = z1
!        do while (z.ge.z2)
!         z = z - 0.01
!         m = m + 1
!         if(m.gt.mdebug)stop'debug dimension needs to be increased'
!          zp(m) = z
!          q = dexp ( - h * dabs(( z - z0 )) )
!          r = GTONE + q
!          pp(m) = a * q / r / r
!        enddo
!        msav = m
!        do m = 1, msav
!          write(81,'(13f10.4)')zp(m),pp(m)
!        enddo
!      endif

!  Finish

      return
end subroutine profiles_gdfone

!

subroutine gdf1_profile_only(za,zb,z0,h,a,prof)

      implicit none

!  Profile assignation

!  Arguments

      real(GTPK),    intent(in)   :: za,zb,z0,h,a
      real(GTPK),    intent(out)  :: prof

!  Local variables

      real(GTPK)     :: fac,ra,rb

!  initialise

      prof = GTZERO

      fac = a / h
      ra = GTONE / ( GTONE + dexp ( - h * dabs(za-z0)) ) 
      rb = GTONE / ( GTONE + dexp ( - h * dabs(zb-z0)) )
      if ( za.gt.z0 ) then
        if ( zb .ge. z0 ) then
          prof = fac * ( ra - rb )
        else 
          prof = fac * ( ra + rb - GTONE )
        endif
      else
        prof = fac * ( rb - ra )
      endif

!  Finish

      return
end subroutine gdf1_profile_only

!

subroutine gdf1_profile_plus ( za,zb,z0,h,a,a_d,a_p,a_w,   & ! Input
                              prof,deriv1,deriv2,deriv3)     ! Output

      implicit none

!  Profile and derivative assignation

!  Arguments

      real(GTPK),    intent(in)   :: za,zb,z0,h,a,a_d,a_p,a_w
      real(GTPK),    intent(out)  :: prof,deriv1,deriv2,deriv3

!  Local variables

      real(GTPK)    :: fac,qa,qb,ra,rb,mza,mzb,rsqa,rsqb,rab,fac_d
      real(GTPK)    :: fac_p,ra_p,rb_p,rab_p
      real(GTPK)    :: fac_w,ra_w,rb_w,rab_w

!  initialise

      prof   = GTZERO
      deriv1 = GTZERO
      deriv2 = GTZERO
      deriv3 = GTZERO

!  Type 1

      fac = a / h
      mza =  dabs(za-z0)
      mzb =  dabs(zb-z0)
      qa = dexp ( - h * mza)
      qb = dexp ( - h * mzb)
      ra = GTONE / ( GTONE + qa ) 
      rb = GTONE / ( GTONE + qb )
      fac_d = a_d/h
      fac_p = a_p/h
      fac_w = (a_w/h) - (fac/h)
      rsqa =  qa * ra * ra
      rsqb =  qb * rb * rb
      ra_p = rsqa * h
      rb_p = rsqb * h
      ra_w = rsqa * mza
      rb_w = rsqb * mzb
      if ( za.gt.z0 ) then
        if ( zb .ge. z0 ) then
          rab   =   ra   - rb
          rab_p = - ra_p + rb_p
          rab_w =   ra_w - rb_w
        else 
          rab   =   ra   + rb   - GTONE
          rab_p = - ra_p + rb_p
          rab_w =   ra_w + rb_w
        endif
      else
        rab    = rb   - ra
        rab_p  = rb_p - ra_p
        rab_w  = rb_w - ra_w
      endif
      prof   = fac   * rab
      deriv1 = fac_d * rab
      deriv2 = fac_p * rab + fac * rab_p
      deriv3 = fac_w * rab + fac * rab_w

!  Finish

      return
end subroutine gdf1_profile_plus

!  End module
 
End Module GEMSTOOL_aerload_routines_m

