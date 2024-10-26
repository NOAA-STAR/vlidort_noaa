program test1

   use vfzmat_Pre_Master_m
   use vfzmat_Post_Master_m

   implicit none

!  precision

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)

!  Inputs
!  ------

!  dimensions

   INTEGER, parameter   :: Max_InAngles = 500
   INTEGER, parameter   :: max_geoms = 125, max_szas = 5, max_vzas  = 5, max_azms = 5

!  Directional flags

   LOGICAL, parameter  :: do_upwelling = .true.,  do_dnwelling = .true.

!  Flags for use of observational geometry, Sunlight
 
   LOGICAL, parameter  :: do_obsgeoms = .false., Sunlight = .true.

!  numbers (general)

   INTEGER  :: nstokes, ncoeffs, nstreams, n_Quadangles

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
   REAL(dpk) :: InAngles    ( Max_InAngles )
   REAL(dpk) :: InCosines   ( Max_InAngles )
   REAL(dpk) :: InFmatrices ( Max_InAngles, 6 )
   REAL(dpk) :: Local_InFmatrices ( 6 )

!  Output from Pre-Master
!  ----------------------

!  Quadrature

   REAL(dpk), allocatable   :: QuadAngles  (:)
   REAL(dpk), allocatable   :: QuadCosines (:)
   REAL(dpk), allocatable   :: QuadWeights (:)

!  Spherical functions for quadrature angles

   REAL(dpk), allocatable :: GSF_P00_Saved  (:,:)
   REAL(dpk), allocatable :: GSF_P02_Saved  (:,:)
   REAL(dpk), allocatable :: GSF_P2p2_Saved  (:,:)
   REAL(dpk), allocatable :: GSF_P2m2_Saved  (:,:)

!  Scattering angle cosines, rotational angles 

   REAL(dpk)  :: COSSCAT_up(max_geoms)
   REAL(dpk)  :: COSSCAT_dn(max_geoms)
   REAL(dpk)  :: C1_up(max_geoms), S1_up(max_geoms)
   REAL(dpk)  :: C2_up(max_geoms), S2_up(max_geoms)
   REAL(dpk)  :: C1_dn(max_geoms), S1_dn(max_geoms)
   REAL(dpk)  :: C2_dn(max_geoms), S2_dn(max_geoms)

!  Output from Post-Master
!  -----------------------

!  Output Fmatrices (Interpolated)

   REAL(dpk)  :: OutFmatrices_up   ( Max_Geoms, 6 )
   REAL(dpk)  :: OutFmatrices_dn   ( Max_Geoms, 6 )

!  Zmatrices

   REAL(dpk)  :: Zmatrices_up(max_geoms,4,4)
   REAL(dpk)  :: Zmatrices_dn(max_geoms,4,4)

!  Fmatrix coefficients

   REAL(dpk), allocatable :: FMatCoeffs(:,:)

!  Other variables
!  ---------------

   INTEGER      :: mask(6)
   data mask / 1, 3, 4, 6, 2, 5 /
   logical       :: do_ice_cloud
   character*100 :: Phase_function_inputfile

   integer :: i, j, k, k1, L, irun, nrun, nrun0
   real    :: drun0, e1, e2, e3, e4


!   top level

   open(1,file='TOPLEVEL.inp', status = 'old' )
   read(1,*)nstokes
   read(1,*)nstreams ; ncoeffs = 2*nstreams
   read(1,*)n_Quadangles
   close(1)
   nrun = 50 ; nrun0 = 10 ; drun0 = real(nrun0)

!   Allocate

   Allocate (QuadAngles(n_Quadangles),QuadCosines(n_Quadangles),QuadWeights(n_Quadangles))
   Allocate (GSF_P00_Saved (n_Quadangles,0:ncoeffs),GSF_P02_Saved (n_Quadangles,0:ncoeffs))
   Allocate (GSF_P2p2_Saved(n_Quadangles,0:ncoeffs),GSF_P2m2_Saved(n_Quadangles,0:ncoeffs))
   Allocate (FMatCoeffs(0:ncoeffs,6))

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
         k1 = N_InAngles + 1 - k 
         read(1,*)Inangles(k),Local_InFmatrices(1:6)
         InCosines(k1) = cos ( InAngles(k) * dtr )
         do L = 1, 6
            InFmatrices(k,L) = Local_InFmatrices(mask(L))
         enddo
         do L = 2, 6
            InFmatrices(k,L) = InFmatrices(k,L) * InFmatrices(k,1)
         enddo
      ENDDO
      CLOSE(1)
   else  
     N_Inangles = Max_InAngles
     do k = 1, N_InAngles
       k1 = N_InAngles + 1 - k 
       Inangles(k)   = 0.0d0 + 0.25d0*real(k-1,dpk)
       InCosines(k1) = cos ( InAngles(k) * dtr )
       InFmatrices(k,1:6) = 1.0_dpk
     enddo
   endif

!  Zero output

   OutFmatrices_up = 0.0_dpk
   OutFmatrices_dn = 0.0_dpk
   Zmatrices_up    = 0.0_dpk
   Zmatrices_dn    = 0.0_dpk
   FMatCoeffs      = 0.0_dpk

!  premaster. perform this a few times

   call cpu_time(e1)
   do irun = 1, nrun0
      call vfzmat_Pre_Master &
      ( max_geoms, max_szas, max_vzas, max_azms, dtr,                            & ! Input Dimensions (VLIDORT)
        do_upwelling, do_dnwelling, do_ObsGeoms, ncoeffs, nstokes, n_QuadAngles, & ! Input Flags and Control
        n_geoms, n_szas, n_vzas, n_azms, offsets, szas, vzas, azms, obsgeoms,    & ! Input Geometries
        C1_up, S1_up, C2_up, S2_up, C1_dn, S1_dn, C2_dn, S2_dn,                  & ! Output rotation angles
        COSSCAT_up, COSSCAT_dn, QuadAngles, QuadCosines, QuadWeights,            & ! Output Scatcosines and Quadrature
        GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved )             ! Output Saved GSFs
   enddo
   call cpu_time(e2)
   write(*,'(A,f10.6)')'timing for single call, Premaster = ',(e2-e1)/drun0


!  Postmaster

   call cpu_time(e3)
   do irun = 1, nrun
      Call vfzmat_Post_Master &
      ( max_geoms, N_InAngles, InCosines(1:n_InAngles), InFmatrices(1:n_InAngles,:), & ! input  Fmatrices
        do_upwelling, do_dnwelling, Sunlight,                         & ! input  Flags
        ncoeffs, nstokes, n_geoms, n_QuadAngles, QuadCosines,         & ! Input numbers and quad
        COSSCAT_up, COSSCAT_dn, C1_up, S1_up, C2_up, S2_up, C1_dn, S1_dn, C2_dn, S2_dn, & ! Input angles
        GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved,                   & ! Input Saved GSFs     
        OutFmatrices_up, OutFmatrices_dn, Zmatrices_up, Zmatrices_dn, FMatCoeffs )        ! Output
   enddo
   call cpu_time(e4)
   write(*,'(A,f10.6)')'timing for 50 points,  Postmaster = ',e4-e3

   write(*,77)'Nstokes/nstreams/nquads = ',nstokes,nstreams,n_QuadAngles,'; timing for Both = ', ((e2-e1)/drun0) + e4 - e3
77 format(A,2i3,i5,A,f10.6)

   do k = 1, n_geoms
     write(201,78)k,OutFmatrices_up(k,1:6)
     write(202,78)k,OutFmatrices_dn(k,1:6)
   enddo
   do L = 0, ncoeffs
      write(203,78)L,FMatCoeffs(L,1:6)
   enddo
78 format(i5, 1p6e18.8)


!  done

   stop
end program test1

