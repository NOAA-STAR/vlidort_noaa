program test1_original

   use vfzmat_Master_m

   implicit none

!  precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER, parameter   :: Max_InAngles = 500, max_moments = 100, maxlayers = 1
   INTEGER, parameter   :: max_geoms = 125, max_szas = 5, max_vzas  = 5, max_azms = 5

!  Directional flags

   LOGICAL, parameter  :: do_upwelling = .true.,  do_dnwelling = .true.

!  Flags for use of observational geometry, Sunlight
 
   LOGICAL, parameter  :: do_obsgeoms = .false., Sunlight = .true.

!  numbers (general)

   INTEGER  :: nstokes, ncoeffs, nlayers, nstreams, n_Quadangles

!  Angles. Convention as for  VLIDORT

   REAL(dpk) :: dtr
   REAL(dpk) :: szas(max_szas)
   REAL(dpk) :: vzas(max_vzas)
   REAL(dpk) :: azms(max_azms)
   REAL(dpk) :: obsgeoms(max_geoms,3)
   INTEGER   :: n_geoms, n_szas, n_vzas, n_azms
   INTEGER   :: offsets(max_szas,max_vzas)

!  Input F-matrices

   INTEGER   :: N_InAngles
   REAL(dpk) :: InAngles          ( Max_InAngles )
   REAL(dpk) :: Local_InFmatrices ( 6 )
   REAL(dpk) :: InFmatrices ( Maxlayers, Max_InAngles, 6 )
   LOGICAL   :: Exist_InFmatrices ( Maxlayers )

!  Output Fmatrices (Interpolated)

   REAL(dpk)  :: OutFmatrices_up   ( Max_Geoms, MaxLayers, 6 )
   REAL(dpk)  :: OutFmatrices_dn   ( Max_Geoms, MaxLayers, 6 )

!  Zmatrices

   REAL(dpk)  :: Zmatrices_up(Max_Geoms, MaxLayers,4,4)
   REAL(dpk)  :: Zmatrices_dn(Max_Geoms, MaxLayers,4,4)

!  Fmatrix coefficients

   REAL(dpk)  :: FMatCoeffs(maxlayers,0:max_moments,6)


!  Other variables
!  ---------------

   INTEGER      :: mask(6)
   data mask / 1, 3, 4, 6, 2, 5 /
   logical       :: do_ice_cloud
   character*100 :: Phase_function_inputfile
   integer :: i, j, k, L, irun, nrun
   real    :: e1, e2

!   top level

   open(1,file='TOPLEVEL.inp', status = 'old' )
   read(1,*)nstokes
   read(1,*)nstreams ; ncoeffs = 2*nstreams
   read(1,*)n_Quadangles
   close(1)
   nrun = 50
   nlayers = 1 ; Exist_InFmatrices(1) = .true.

!  setups

   n_geoms = 125 ;  n_szas = 5 ; n_vzas = 5 ; n_azms = 5
   szas(1:5) = 30 
   dtr = acos(-1.0_dpk)/180.0_dpk
   Obsgeoms = 0.0_dpk
   do i = 1, n_szas
      szas(i) = 30.0_dpk + 5.0_dpk * real(i-1,dpk)
      do j = 1, n_vzas
         vzas(j) = 17.0_dpk + 4.0_dpk * real(j-1,dpk)
         offsets(i,j) = n_vzas * (j-1) + n_szas * n_vzas * (i-1)
         azms(1:5) = 55.0_dpk       
   enddo ; enddo

!  Cosines of input angles.
!    -- Reverse directions for the Splining (Monotonically increasing)

!  Real ice clouds

   do_ice_cloud = .true.
   Phase_function_inputfile = 'Ice_Fmatrices.txt'

   IF (do_ice_cloud) THEN
      N_Inangles = 498
      OPEN (1,FILE = Trim(Phase_function_inputfile),STATUS='old')
      DO k=1,15
         read(1,*)
      ENDDO
      DO k = 1 ,N_Inangles
         read(1,*)Inangles(k),Local_InFmatrices(1:6)
         do L = 1, 6
            InFmatrices(1,k,L) = Local_InFmatrices(mask(L))
         enddo
         do L = 2, 6
            InFmatrices(1,k,L) = InFmatrices(1,k,L) * InFmatrices(1,k,1)
         enddo
      ENDDO
      CLOSE(1)
   else  
     N_Inangles = Max_InAngles
     do k = 1, N_InAngles
       Inangles(k)        = 0.0d0 + 0.25d0*real(k-1,dpk)
       InFmatrices(1,k,1:6) = 1.0_dpk
     enddo
   endif

   call cpu_time(e1)
   do irun = 1, nrun
      Call vfzmat_Master &
        ( max_moments, max_geoms, max_szas, max_vzas, max_azms, maxlayers,    & ! input  Dimensions (VLIDORT)
          max_InAngles, N_InAngles, InAngles, InFmatrices, Exist_InFmatrices, & ! input  Fmatrices
          do_upwelling, do_dnwelling, do_ObsGeoms, Sunlight,                  & ! input  Flags
          ncoeffs, nlayers, nstokes, n_geoms, n_szas, n_vzas, n_azms,         & ! input  Numbers
          offsets, dtr, szas, vzas, azms, obsgeoms,                           & ! Input  Geometries
          OutFmatrices_up, OutFmatrices_dn, Zmatrices_up, Zmatrices_dn, FMatCoeffs )

   enddo
   call cpu_time(e2)
   write(*,77)'Nstokes/nstreams/nquads = ',nstokes,nstreams,2000,'; timing for 50 points, Orig Master = ',e2-e1
77 format(A,2i3,i5,A,f10.6)

   do k = 1, n_geoms
     write(101,78)k,OutFmatrices_up(k,1,1:6)
     write(102,78)k,OutFmatrices_dn(k,1,1:6)
   enddo
   do L = 0, ncoeffs
      write(103,78)L,FMatCoeffs(1,L,1:6)
   enddo
78 format(i5, 1p6e18.8)

!  done

   stop
end program test1_original

