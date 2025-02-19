
MODULE crtm_lbl_simulator

  ! Upgrade for Version 2.8.3, August 2023
  ! ---------------------------------------
  USE CRTM_Atmosphere_Define,     ONLY: CRTM_Atmosphere_type, &
                                        CRTM_Atmosphere_Destroy, &
                                        CRTM_Atmosphere_IsValid, &
                                        CRTM_Get_PressureLevelIdx
  USE CRTM_Surface_Define,        ONLY: CRTM_Surface_type, &
                                        CRTM_Surface_IsValid
  USE CRTM_Geometry_Define,       ONLY: CRTM_Geometry_type, &
                                        CRTM_Geometry_IsValid
  USE CRTM_Options_Define,        ONLY: CRTM_Options_type, &
                                        CRTM_Options_IsValid

  USE CRTM_RTSolution_Define  ,   ONLY: CRTM_RTSolution_type

  ! module files of vfzmat
  USE vfzmat_Rayleigh_m
  USE vfzmat_Pre_Master_m
  USE vfzmat_Post_Master_m


  USE VLIDORT_PARS_m
  USE VLIDORT_IO_DEFS_m

  USE VLIDORT_VBRDF_SUP_ACCESSORIES_m

  USE VLIDORT_AUX_m,    Only : VLIDORT_READ_ERROR, VLIDORT_WRITE_STATUS
  USE VLIDORT_INPUTS_m, Only : VLIDORT_INPUT_MASTER, VLIDORT_Sup_Init
  USE VLIDORT_MASTERS_m

  USE PCRTM_Sw_Aux_m
  USE GAS_OPT_m

  USE pcrtm_Interp_Utility
  USE pcrtm_File_Utility

  IMPLICIT NONE
  PRIVATE

  PUBLIC  ::  crtm_omps_simulator
  INTEGER, PUBLIC :: N_USER_Channels
 
  CONTAINS

    FUNCTION crtm_omps_simulator(  &
    Atmosphere , &  ! Input, M
    Surface    , &  ! Input, M
    Geometry   , &  ! Input, M
    RTSolution , &  ! Input, L,M
    Options    ) &  ! Optional input, M
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),        INTENT(IN OUT) :: Atmosphere
    TYPE(CRTM_Surface_type),           INTENT(IN)     :: Surface
    TYPE(CRTM_Geometry_type),          INTENT(IN)     :: Geometry
    TYPE(CRTM_RTSolution_type),        INTENT(IN OUT) :: RTSolution(:)
    TYPE(CRTM_Options_type), OPTIONAL, INTENT(IN)     :: Options
    ! Function result
    INTEGER :: Error_Status



    ! VLIDORT file inputs status structure
    TYPE(VLIDORT_Input_Exception_Handling) :: VLIDORT_InputStatus

    ! VLIDORT debug input control
    LOGICAL :: DO_DEBUG_INPUT

    ! VLIDORT input structures
    TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
    TYPE(VLIDORT_ModIFied_Inputs)          :: VLIDORT_ModIn

    ! VLIDORT supplements i/o structure
    TYPE(VLIDORT_Sup_InOut)                :: VLIDORT_Sup

    ! VLIDORT output structure
    TYPE(VLIDORT_Outputs)                  :: VLIDORT_Out

    ! Input F-matrix stuff ( angles and 6 scattering matrix entries )

    ! !   INTEGER, parameter :: Max_InAngles = 18001
    INTEGER   :: Ncount
    INTEGER   :: N_InAngles_Tamu,N_InAngles_Opac

    ! Output from Pre-Master
    ! ----------------------

    ! Quadrature
    DOUBLE PRECISION, ALLOCATABLE   :: QuadAngles  (:)
    DOUBLE PRECISION, ALLOCATABLE   :: QuadCosines (:)
    DOUBLE PRECISION, ALLOCATABLE   :: QuadWeights (:)

    ! Spherical functions for quadrature angles
    DOUBLE PRECISION, ALLOCATABLE :: GSF_P00_Saved  (:,:)
    DOUBLE PRECISION, ALLOCATABLE :: GSF_P02_Saved  (:,:)
    DOUBLE PRECISION, ALLOCATABLE :: GSF_P2p2_Saved  (:,:)
    DOUBLE PRECISION, ALLOCATABLE :: GSF_P2m2_Saved  (:,:)

    ! Scattering angle cosines, rotational angles
    DOUBLE PRECISION  :: COSSCAT_up(max_geometries)
    DOUBLE PRECISION  :: COSSCAT_dn(max_geometries)
    DOUBLE PRECISION  :: C1_up(max_geometries), S1_up(max_geometries)
    DOUBLE PRECISION  :: C2_up(max_geometries), S2_up(max_geometries)
    DOUBLE PRECISION  :: C1_dn(max_geometries), S1_dn(max_geometries)
    DOUBLE PRECISION  :: C2_dn(max_geometries), S2_dn(max_geometries)

    LOGICAL ::  quadrature_gen_flag

    ! Proxies for the VFZMAT supplement
    ! ---------------------------------

    ! Flags (observational geometry, Sunlight), Geometry numbers

    LOGICAL ::  DO_OBSGEOMS, DO_SUNLIGHT, DO_UPWELLING, DO_DNWELLING
    LOGICAL ::  DO_DOublets, DO_Planetary

    INTEGER ::          N_GEOMS, N_SZAS, N_VZAS, N_AZMS
    INTEGER   :: Lattice_offsets(MAX_SZANGLES, MAX_USER_VZANGLES)
    INTEGER   :: DOublet_offsets(MAX_SZANGLES)

    INTEGER, PARAMETER :: n_Quadangles = 5000
    INTEGER :: nstokes, ncoeffs, nlayers, nstreams

    ! Angles. Convention as for  VLIDORT

    DOUBLE PRECISION :: SZAS (MAX_SZANGLES)
    DOUBLE PRECISION :: VZAS (MAX_USER_VZANGLES)
    DOUBLE PRECISION :: AZMS (MAX_USER_RELAZMS)
    DOUBLE PRECISION :: OBSGEOMS (MAX_GEOMETRIES,3)


    ! The VFZMAT supplemental variables for Rayleigh
    ! Output Fmatrices (Calculated from Coefficients), Zmatrices, rayleigh coefficients
    ! --- The Rayleigh coefficients are "PROBLEM_RAY"
    DOUBLE PRECISION :: RayFmatrices_up   ( MAX_GEOMETRIES, 6 )
    DOUBLE PRECISION :: RayFmatrices_dn   ( MAX_GEOMETRIES, 6 )
    DOUBLE PRECISION :: RayZmatrices_up   ( MAX_GEOMETRIES, 4, 4 )
    DOUBLE PRECISION :: RayZmatrices_dn   ( MAX_GEOMETRIES, 4, 4 )
    DOUBLE PRECISION :: RayCoeffs(0:2,6)


    ! The VFZMAT supplemental variables for Aerosols, same for dIFferent layers
    ! Output Fmatrices (Interpolated), Zmatrices, Fmatrix coefficients
    DOUBLE PRECISION :: OutFmatrices_up   ( MAX_GEOMETRIES,  6 )
    DOUBLE PRECISION :: OutFmatrices_dn   ( MAX_GEOMETRIES,  6 )
    DOUBLE PRECISION :: Zmatrices_up      ( MAX_GEOMETRIES,  4, 4 )
    DOUBLE PRECISION :: Zmatrices_dn      ( MAX_GEOMETRIES,  4, 4 )

    DOUBLE PRECISION, ALLOCATABLE :: FmatCoeffs(:,:)

    ! Local Variables
    ! ===============

    ! Flag for opening error output file

    LOGICAL ::          OPENFILEFLAG
    INTEGER ::          K
    INTEGER ::          N,V
    DOUBLE PRECISION :: DEPOL,CO2_PPMV_MIXRATIO

    INTEGER::           k1,k2,ib,um

    ! VLIDORT standard input preparation
    LOGICAL ::          DO_FOCORR, DO_FOCORR_NADIR, DO_FOCORR_OUTGOING
    LOGICAL ::          DO_DELTAM_SCALING, DO_SOLUTION_SAVING, DO_BVP_TELESCOPING
    INTEGER ::         NFINELAYERS


    ! Optical proxies. Fmatrix proxies are new for Version 2.8
    DOUBLE PRECISION :: OMEGA_TOTAL_INPUT ( MAXLAYERS )
    DOUBLE PRECISION :: DELTAU_VERT_INPUT ( MAXLAYERS )
    DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
    DOUBLE PRECISION :: LAMBERTIAN_ALBEDO
    DOUBLE PRECISION :: FMATRIX_UP ( MAXLAYERS, MAX_GEOMETRIES, 6 )
    DOUBLE PRECISION :: FMATRIX_DN ( MAXLAYERS, MAX_GEOMETRIES, 6 )

    ! Proxies (Map no longer required for Version 2.8)
    INTEGER ::          N_USER_LEVELS
    DOUBLE PRECISION :: USER_LEVELS ( MAX_USER_LEVELS )

    ! VLIDORT standard output preparation
    INTEGER ::          N_GEOMETRIES
    INTEGER ::          N_SZANGLES

    INTEGER  :: indx_lambertian
    INTEGER   :: nbeams, N_USER_STREAMS, N_USER_RELAZMS
    REAL*8      ::  wavelen
    INTEGER :: i
    REAL    :: e1,e2

    ! !  OD LUT a
    REAL,DIMENSION(:,:), ALLOCATABLE ::odIO_tot

    REAL,DIMENSION(:,:,:),ALLOCATABLE::Iup_output,Idn_output,Qup_output,Qdn_output,Uup_output,Udn_output
    REAL,DIMENSION(:,:),ALLOCATABLE::Iup_output2,Qup_output2,Uup_output2
    REAL,DIMENSION(:,:),ALLOCATABLE::BTRANS_I,BTRANS_Q,BTRANS_U
    REAL,DIMENSION(:), ALLOCATABLE :: BSPHER
    REAL,DIMENSION(:),  ALLOCATABLE::SolTran_output
    REAL,DIMENSION(:,:),ALLOCATABLE::Iup_output2_usr,Qup_output2_usr,Uup_output2_usr
    REAL,DIMENSION(:,:),ALLOCATABLE::BTRANS_I_usr,BTRANS_Q_usr,BTRANS_U_usr
    REAL,DIMENSION(:), ALLOCATABLE :: BSPHER_usr
    REAL(8),DIMENSION(:), ALLOCATABLE :: wnum_usr

    INTEGER                          :: ngases
    INTEGER                          :: ndat,w,checkrun

    REAL(8), DIMENSION(:),ALLOCATABLE :: wavenums, lambdas
    INTEGER, DIMENSION(:),ALLOCATABLE :: waveindex

    REAL(8), DIMENSION(:), ALLOCATABLE :: rayleigh_xsec,rayleigh_depol

    REAL(8)                          :: wresol
    REAL(8), DIMENSION(1:MAX_SZANGLES)          :: szangles,user_relazms,user_vzangles

    REAL                             :: t_start,t_end !cpu time


    ! ! ATM profile
    INTEGER                          :: nProf, nprof_st,nprof_end

    ! ! Filename related
    CHARACTER(LEN=300)                      :: profile_data_filename,nameProf


    ! SfcRef variable
    INTEGER,PARAMETER :: nwav_sfcref =2152
    REAL, DIMENSION(1:nwav_sfcref, 1:1946) :: SfcRef          !sfcRef has 1946 cases
    REAL, DIMENSION(1:nwav_sfcref)         :: SfcRef1D
    REAL, DIMENSION(1:nwav_sfcref)         :: lambSfcRef
    REAL, DIMENSION(:), ALLOCATABLE        :: SfcRef_usr,lambdas_usr

    CHARACTER(LEN=4)                      :: nprof_s



    ! RanDOm ocean surface windspeed and salinity
    REAL, DIMENSION(nMerra)              :: rand_windspeed, rand_salinity   ! 8/16/2023

    INTEGER  :: st_aergrid,end_aergrid            !index of LUT/opac aer grid
    LOGICAL  :: DO_Aerosol, DO_Lambertian

    ! ! aerosol profile adapted from PCA-vlidort
    INTEGER :: nmode_aer_in,nmode_opac_in,nmode_dust_in
    INTEGER :: imod
    INTEGER :: flag_aer(maxtypeaer,MAXLAYERS)

    DOUBLE PRECISION, ALLOCATABLE :: FmatCoeffs_mode(:,:,:,:)
    ! rev: ming this dimension could be large, need to revise accordingly
    DOUBLE PRECISION, ALLOCATABLE ::FMatCoeffs_tot(:,:,:,:)!(MAXLAYERS,max_aerwave,MAX_ncoeff,6)

    DOUBLE PRECISION :: OutFmat_up_intp(MAXLAYERS,MAX_GEOMETRIES,6)
    DOUBLE PRECISION :: OutFmat_dn_intp(MAXLAYERS,MAX_GEOMETRIES,6)
    DOUBLE PRECISION, ALLOCATABLE ::FMatCoeffs_intp(:,:,:)

    ! Rev: Ming extinction coefficient [1/km] at reference wavelength for each aerosol type

    DOUBLE PRECISION ::aerext_intp(MAXLAYERS),aersca_intp(MAXLAYERS),waer_intp(MAXLAYERS)


    ! AFGL profile variable
    INTEGER  :: ATM_TYPE, j
    INTEGER  :: iband           !index of LUT band


    INTEGER :: RanDOm_Seed_times, nstep

    CHARACTER(LEN=80)  :: path0, path1
    REAL(8)                          :: wnStart, wnEnd
    INTEGER,DIMENSION(1:maxgases)          :: gasIdBand !gasIDall! !hitran ID of all abs. gases in the band

    INTEGER :: MerraProfID
    REAL    :: Psfc_in
    REAL    :: psfc_out
    REAL    :: temp_at_psfc
    INTEGER :: nlevel_usr
    REAL    :: pres_usr(1:maxlayers)
    REAL    :: vmr_usr(1:MAXLAYERS, 1:ngas_lut)
    REAL    :: temp_usr(1:MAXLAYERS)

    REAL    :: vmr_atm(1:maxAtmlev, 1:ngas_lut) 
    REAL    :: pres_atm  ( 1:maxAtmlev ) 
    REAL    :: temp_atm  ( 1:maxAtmlev ) 
    INTEGER :: nlevel_atm, ngas_atm 

    REAL    :: vmr_ref(1:maxAtmlev, 1:ngas_lut)
    REAL    :: pres_ref  ( 1:maxAtmlev )
    REAL    :: temp_ref  ( 1:maxAtmlev )
    INTEGER :: nlevel_ref, ngas_ref
    REAL    :: psfc_out_ref


    REAL, ALLOCATABLE :: Wgas_usr (:,:)
    REAL, ALLOCATABLE :: gph_usr (:)
    REAL, ALLOCATABLE :: airAmt_usr (:)
    REAL*8,ALLOCATABLE:: airClm_vld(:)
    REAL*8,ALLOCATABLE:: htgrid_vld(:)
    REAL*8,ALLOCATABLE:: odIOTot_vld(:,:)


    INTEGER             ::  st_band, end_band, N_DOublets, N_USER_OBSGEOMS
    REAL, DIMENSION(1:MAX_GEOMETRIES) :: Out_szas,Out_vzas, Out_azms   ! three angles for output
    CHARACTER(256) :: Buffer
    LOGICAL :: Test_hotspot, DO_Scalar_only
    INTEGER ::  st_wvnum_indx, end_wvnum_indx

    CHARACTER(8)         :: date
    CHARACTER(10)        :: time
    CHARACTER(5)         :: zone
    INTEGER,DIMENSION(8) :: values

    ! add user set stream, default =6 , 11/20/2022
    LOGICAL :: DO_User_Stream
    INTEGER :: userset_stream = 6

    CHARACTER(LEN=80)  :: PresGridFile
    LOGICAL :: DO_mapgrid,DO_usrgrid,DO_stdgrid
    LOGICAL :: DO_opac, DO_tamu

    INTEGER :: ndat_lut,nwv_usr

    REAL*8, DIMENSION(:),ALLOCATABLE :: wv_usr
    REAL*8, DIMENSION(:,:),ALLOCATABLE :: wv_band
    INTEGER, DIMENSION(:,:),ALLOCATABLE :: wvidx_band

    INTEGER:: ndat_band (nband_lut)

    CHARACTER(LEN=300)   :: file_usrwv

    INTEGER   :: nwav_st,nwav_end
    INTEGER   :: Ncount_usr
    INTEGER :: IFile

    REAL    :: H2o_user(1:maxAtmlev),O3_user(1:maxAtmlev),temp_user(maxAtmlev)  ! note this dIFference with temp_usr !check to make consistent later
    REAL*8 ::  WINDSPEED, Salinity, RefIdx_R, RefIdx_I

    ! icetest height at temp. lt 273k
    REAL    :: height_T0
    LOGICAL :: DO_icecld
    INTEGER :: NMESSAGES

    LOGICAL :: VLIDORT_CRTM = .true.

    CALL cpu_time(e1)

    ! working directory
    path0='./'
    path1='./PCRTM_VLIDORT_Config/'

    CALL date_and_time(date,time,zone,values)

    ! --SECTION 1: setting for VLIDORT run
    ! 1.1 User control flags
    CALL read_user_control(DO_mapgrid,DO_usrgrid,    & !pres grid
         ATM_TYPE, Psfc_in,                 & !wave grid
         DO_Lambertian,LAMBERTIAN_ALBEDO, &
         DO_Planetary,                 & !surface and planetary
         DO_Scalar_only,userset_stream,               & !stokes and stream
         N_USER_OBSGEOMS, &
         szangles,user_relazms,user_vzangles)

    Psfc_in = Atmosphere%pressure(size(Atmosphere%Pressure))
    IF (VLIDORT_CRTM) THEN
         N_USER_OBSGEOMS = 1
         szangles(1) = Geometry%Source_Zenith_Angle
         user_vzangles(1) = Geometry%Sensor_Zenith_Angle
         user_relazms(1) = Geometry%Source_Azimuth_Angle -  Geometry%Sensor_Azimuth_Angle 
    END IF

    ! 1.2 read in user spectral grid:
    file_usrwv='./PCRTM_VLIDORT_Config/user_wav_inputs.dat'
    CALL read_usr_wv(path0,file_usrwv,nwv_usr,wv_usr)  !xiong: add the filename
    ! setting band based on wv_usr
    CALL Set_Usr_wvgrid(path0,nwv_usr,wv_usr,&                                !input
         st_band,end_band,wv_band,wvidx_band,ndat_band)         !output
    IF (end_band .LT. st_band) THEN
      PRINT*, "Error in Band setting !"
      STOP
    ENDIF
    PRINT*,'DO LUT band from Band', st_band,' to band', end_band

    ! 1.3 read pressure grid
    IF (DO_mapgrid) THEN
      PresGridFile = 'map_pres_grid.dat'
      CALL read_Pres_Grid(path0,PresGridFile, nlevel_usr,pres_usr)
    ELSE IF (DO_usrgrid) THEN !default the 50 level mapping grid
      PresGridFile = 'user_pres_grid.dat_sample'
      CALL read_Pres_Grid(path0,PresGridFile,nlevel_usr,pres_usr)
    ENDIF

   
    ! 1.4 fixed setting
    nprof_st = 1
    nProf    = nprof_st
    WRITE(nprof_s, '(i4.4)') nProf
   
    ! 1.5 read atm profile (1~6AFGL atm profile; 0 user defined profile)
    IF (atm_type >=1 .AND. atm_type<=6) THEN
      CALL Read_AFGL_Prof(path0,    ATM_type,Psfc_in,&                              !input
           nlevel_ref,ngas_ref,pres_ref,temp_ref,vmr_ref,psfc_out_ref,&   !output
           nameProf,profile_data_filename)                            !output
    ELSE IF (atm_type .EQ. 0) THEN
      profile_data_filename = './data_and_control/profile_201702_290.500.dat.mid.land'
      psfc_in = 1013.
      nameProf='profile_201702_290.500.dat'
      nProf_s='000'
      ! first read defult profile
      CALL Read_Merra_Prof(profile_data_filename,psfc_in,&                            !input
           nlevel_ref,ngas_ref,pres_ref,temp_ref,vmr_ref,psfc_out_ref)
      
      IF(.NOT. VLIDORT_CRTM) THEN
        !THEN read user defined temperature, H2O and O3 profile
        profile_data_filename = './data_and_control/Inputs_User_profile_Geo.dat_NP_J2_midlat'
        CALL READ_UserProf_Geo_INPUTS(profile_data_filename, &  !latitude, longitude, &
                   temp_user, H2o_user,O3_user)
     
       temp_ref(:)  = temp_user(:)
       vmr_ref(:,1) = H2o_user(:)
       vmr_ref(:,3) = O3_user(:)
      ENDIF

    ELSE
      PRINT*,'warning: atm type must between 1 and 6!'
      STOP
    ENDIF

    IF (VLIDORT_CRTM) THEN
      ! 1.5.1 interpolating reference gas profile to CRTM  pressure grid
      nlevel_atm = size(Atmosphere%Pressure)
      ngas_atm = ngas_ref
      pres_atm(1:nlevel_atm) = Atmosphere%Pressure(:)
      CALL Interpolate_User_grid(pres_atm,nlevel_atm,nlevel_ref, psfc_out_ref,&             !input
         pres_ref, temp_ref,vmr_ref,&                          !input
         temp_atm,vmr_atm,nlayers,temp_at_psfc)                !output
      !ngas_atm = 2
      !vmr_atm = 0.0
      temp_atm(1:nlevel_atm) = Atmosphere%Temperature(:)
      vmr_atm(1:nlevel_atm, 1) = Atmosphere%Absorber(:,1) * (28.97/18.015)/1000.0 ! g/kg -> vmr
      vmr_atm(1:nlevel_atm, 3) = Atmosphere%Absorber(:,2) *1.0e-6 ! ppmv to vmr
      psfc_out = Atmosphere%Level_Pressure(nlevel_atm)
      nameProf = 'CRTM_Profile'
      !print*, size(Atmosphere%Pressure),size(Atmosphere%Level_Pressure),psfc_out

    ELSE
      nlevel_atm = nlevel_ref
      ngas_atm = ngas_ref
      pres_atm = pres_ref  
      temp_atm = temp_ref
      psfc_out = psfc_out_ref
    ENDIF


    ! 1.6 interpolating atm profile to user slected pressure grid
    CALL Interpolate_User_grid(pres_usr,nlevel_usr,nlevel_atm,psfc_out,&             !input
         pres_atm, temp_atm,vmr_atm,&                          !input
         temp_usr,vmr_usr,nlayers,temp_at_psfc)                !output

    ! 1.7 define layer quantities
    IF (allocated(Wgas_usr) ) deallocate(Wgas_usr)
    ALLOCATE ( Wgas_usr(1:nlayers, 1:ngas_lut) )

    IF (allocated(gph_usr) ) deallocate(gph_usr)
    ALLOCATE ( gph_usr(1:nlayers+1) )

    IF (allocated(airAmt_usr) ) deallocate(airAmt_usr)
    ALLOCATE ( airAmt_usr(1:nlayers) )

    IF (allocated(htgrid_vld) ) deallocate(htgrid_vld)     !height grid for vlidort
    ALLOCATE ( htgrid_vld(0:nlayers) )

    ! 1.8 calculate height grid based on geo-potential height
    CALL Calc_GPH(nlayers,psfc_out,temp_at_psfc,pres_usr,temp_usr,vmr_usr,DO_icecld,&    !input
         height_T0,gph_usr)      !output
    htgrid_vld(0:nlayers) =  dble(gph_usr(1:nlayers+1)/1000.)  ![km]


    ! 1.9 VLIDORT default control Read input, abort IF failed
    CALL VLIDORT_INPUT_MASTER ( &
         trim(path1)//'Vlidort_cfg/V2p8p3_VLIDORT_ReadInput.cfg', & ! Input
         VLIDORT_FixIn,      & ! Outputs
         VLIDORT_ModIn,      & ! Outputs
         VLIDORT_InputStatus ) ! Outputs

    VLIDORT_FixIn%Cont%TS_NSTREAMS = userset_stream   ! 11/20/2022    ! xiong
    IF (DO_Scalar_only) VLIDORT_FixIn%Cont%TS_NSTOKES = 1   ! 11/29/2022, Xiong

    IF (DO_PLANETARY) THEN
      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .TRUE.
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .TRUE.
    ELSE IF (DO_Lambertian ) THEN
      VLIDORT_FixIn%Bool%TS_DO_PLANETARY_PROBLEM = .FALSE.
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .TRUE.
    ELSE
      VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .FALSE.
    ENDIF

    IF (DO_Planetary) THEN
      PRINT*, "#### DOing PLANETARY Calculation ##### "
    ELSE IF (DO_Lambertian) THEN
      PRINT*, "#### DOing LAMBERTIAN SURFACE ##### "
    ELSE
      PRINT*, "#### DOing VBRDF SURFACE ##### "
    ENDIF

    ! set FO flags
    DO_FOCORR          = .TRUE.
    DO_FOCORR_NADIR    = .TRUE.   ! should be false for more accurate for higher vza > 35 ?? 5/19/2023
    DO_FOCORR_OUTGOING = .FALSE.

    DO_DELTAM_SCALING  = .TRUE.
    DO_SOLUTION_SAVING = .FALSE.  ! Xiong .TRUE.
    DO_BVP_TELESCOPING = .FALSE.
    NFINELAYERS        = 0         ! Non zero here
    IF (DO_FOCORR_OUTGOING)  THEN
      DO_FOCORR_NADIR = .FALSE.     ! 8/24/2022
      NFINELAYERS        = 4         ! Non zero here
    ENDIF

    VLIDORT_ModIn%MBool%TS_DO_FOCORR          = DO_FOCORR
    VLIDORT_ModIn%MBool%TS_DO_FOCORR_NADIR    = DO_FOCORR_NADIR
    VLIDORT_ModIn%MBool%TS_DO_FOCORR_OUTGOING = DO_FOCORR_OUTGOING
    VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING  = DO_DELTAM_SCALING
    VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING = DO_SOLUTION_SAVING
    VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING = DO_BVP_TELESCOPING
    VLIDORT_FixIn%Cont%TS_NFINELAYERS         = NFINELAYERS
    VLIDORT_ModIn%MBool%TS_DO_SSCORR_USEFMAT = .TRUE.

    IF ( VLIDORT_InputStatus%TS_STATUS_INPUTREAD .NE. VLIDORT_SUCCESS ) &
         CALL VLIDORT_READ_ERROR ( 'V2p8p3_VLIDORT_ReadInput.log', VLIDORT_InputStatus )
    VLIDORT_FixIn%Cont%TS_ASYMTX_TOLERANCE = 1.0d-10   ! xiong
    VLIDORT_FixIn%Bool%TS_DO_MSSTS = .FALSE.   !!*** should set it to .true., 5/19/2023
    VLIDORT_FixIn%Bool%TS_DO_FOURIER0_NSTOKES2 = .FALSE.  ! xiong

    CALL VLIDORT_Sup_Init ( VLIDORT_Sup )

    ! Set Output-level Proxies (saved values)
    nstokes       = VLIDORT_FixIn%Cont%TS_nstokes
    n_szangles    = VLIDORT_ModIn%MSunrays%TS_n_szangles

    ! 1.10 read surface reflectance data
    IF (VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) &
         CALL Read_LambSfc_Alb (path0,lambSfcRef,SfcRef)

    ! initial setting bf band loop
    nwav_st = 0
    nwav_end = 0
    quadrature_gen_flag = .FALSE.
    ncount_usr = 0

    ! 1.10 allocate outputs
    IF (allocated(Iup_output2_usr ) ) deallocate( Iup_output2_usr  )
    ALLOCATE (Iup_output2_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

    IF (allocated(Qup_output2_usr ) ) deallocate( Qup_output2_usr  )
    ALLOCATE (Qup_output2_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

    IF (allocated(Uup_output2_usr ) ) deallocate( Uup_output2_usr  )
    ALLOCATE (Uup_output2_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

    IF (allocated(wnum_usr) ) deallocate( wnum_usr)
    ALLOCATE (wnum_usr(1:nwv_usr) )

    IF ( DO_Planetary) THEN   ! 11-28-2023
      IF (allocated(BTRANS_I_usr) ) deallocate( BTRANS_I_usr  )
      ALLOCATE (BTRANS_I_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

      IF (allocated(BTRANS_Q_usr) ) deallocate( BTRANS_Q_usr  )
      ALLOCATE (BTRANS_Q_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

      IF (allocated(BTRANS_U_usr) ) deallocate( BTRANS_U_usr  )
      ALLOCATE (BTRANS_U_usr(1:N_USER_OBSGEOMS, 1:nwv_usr) )

      IF (allocated(BSPHER_usr) ) deallocate( BSPHER_usr)
      ALLOCATE (BSPHER_usr(1:nwv_usr) )
    ENDIF


    ! --SECTION 2:  B A N D    L O O P

    DO iband = st_band, end_band
      PRINT*,iband, 'band loop from', st_band, end_band

      ! 2.1 readin the gasTable to identIFy which gases are selected in each band
      CALL Set_Gas_Band(path0,iband,wnStart, wnEnd, wresol,ndat_lut, ngases,gasIdBand)

      ndat = ndat_band(iband)
      IF (allocated(wavenums) ) deallocate(wavenums)
      ALLOCATE( wavenums(1:ndat) )
      IF (allocated(lambdas) ) deallocate(lambdas)
      ALLOCATE( lambdas(1:ndat) )
      IF (allocated(waveindex) ) deallocate(waveindex)
      ALLOCATE( waveindex(1:ndat) )

      IF (ndat .EQ. 0)  goto 2017            !skip this band IF no user grid in the band

      DO n = 1, ndat
        wavenums(n) = wv_band(iband, n)
        lambdas(n) = 1.0D+07/wavenums(n)
        waveindex(n) = wvidx_band(iband, n)
        PRINT*,'wavenums(n)',wavenums(n),waveindex(n)
      ENDDO


      ! 2.2: gas related setting and calculating
      IF (allocated(RAYLEIGH_XSEC) ) deallocate(RAYLEIGH_XSEC)
      ALLOCATE ( RAYLEIGH_XSEC(ndat) )
      IF (allocated(RAYLEIGH_DEPOL) ) deallocate(RAYLEIGH_DEPOL)
      ALLOCATE ( RAYLEIGH_DEPOL(ndat) )
      IF (allocated(odIO_tot) ) deallocate(odIO_tot)
      ALLOCATE ( odIO_tot(1:nlayers, 1:ndat) )
      IF (allocated(airClm_vld) ) deallocate(airClm_vld)
      ALLOCATE ( airClm_vld(1:nlayers) )
      IF (allocated(odIOTot_vld) ) deallocate(odIOTot_vld)
      ALLOCATE ( odIOTot_vld(1:nlayers,1:ndat) )

      ! 2.2.1 rayleigh
      CO2_PPMV_MIXRATIO = 390.   ! old value
      CALL Rayleigh_function &
           ( ndat, CO2_PPMV_MIXRATIO, &
           ndat,   LAMBDAS, &
           RAYLEIGH_XSEC, RAYLEIGH_DEPOL )


      ! 2.2.2 get gas OD
      CALL gas_OD_LayerAmt (path0,nlayers,ngases,wnStart,wnEnd,ndat,ndat_lut,wavenums, waveindex,gasIdBand, & !input
           pres_usr,nlevel_atm,psfc_out,pres_atm,temp_atm,vmr_atm,&                          !input
           odIO_tot,airAmt_usr,Wgas_usr)                                                     !output

      airClm_vld(1:nlayers) = dble(airAmt_usr(1:nlayers))
      odIOTot_vld(1:nlayers,:) = dble(odIO_tot(1:nlayers,:) )



      ! 2.3: Assign vlidort variables
      VLIDORT_FixIn%Chapman%TS_PRESSURE_GRID(1:Nlayers+1) = pres_usr(1:Nlayers+1)
      VLIDORT_FixIn%Chapman%TS_TEMPERATURE_GRID(1:Nlayers+1) = temp_usr(1:Nlayers+1)
      VLIDORT_FixIn%Cont%TS_nlayers = NLAYERS
      ! set vlidrot output levels
      IF (NLAYERS > 5) THEN
        VLIDORT_FixIn%UserVal%TS_n_user_levels = 6
        VLIDORT_ModIn%MUserVal%TS_user_levels(1) = 0.0
        VLIDORT_ModIn%MUserVal%TS_user_levels(2) = 1.0
        VLIDORT_ModIn%MUserVal%TS_user_levels(3) = 2.0
        VLIDORT_ModIn%MUserVal%TS_user_levels(4) = NLAYERS-2
        VLIDORT_ModIn%MUserVal%TS_user_levels(5) = NLAYERS-1
        VLIDORT_ModIn%MUserVal%TS_user_levels(6) = NLAYERS
        ! Set Output-level Proxies (saved values)
        n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
        USEr_levels(1:n_user_levels) =VLIDORT_ModIn%MUserVal%TS_user_levels(1:n_user_levels)
      ELSE
        VLIDORT_FixIn%UserVal%TS_n_user_levels = 1
        VLIDORT_ModIn%MUserVal%TS_user_levels(1) = 0.0
        ! Set Output-level Proxies (saved values)
        n_user_levels = VLIDORT_FixIn%UserVal%TS_n_user_levels
        USEr_levels(1:n_user_levels) =VLIDORT_ModIn%MUserVal%TS_user_levels(1:n_user_levels)
      ENDIF

      ! 2.4: setting geometry
      ! ! update the VLIDORT read-input by user configure
      DO_Obsgeoms = .TRUE.
      DO_DOublets = .FALSE.
      VLIDORT_ModIn%MBool%TS_DO_OBSERVATION_GEOMETRY =DO_Obsgeoms
      VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY =DO_DOublets

      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,1) = szangles(1:N_USER_OBSGEOMS)
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,2) = user_vzangles(1:N_USER_OBSGEOMS)
      VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(1:N_USER_OBSGEOMS,3) = user_relazms(1:N_USER_OBSGEOMS)

      VLIDORT_ModIn%MSunrays%TS_szangles(1:N_USER_OBSGEOMS) = szangles(1:N_USER_OBSGEOMS)
      VLIDORT_ModIn%MUserVal%TS_user_relazms(1:N_USER_OBSGEOMS) = user_relazms(1:N_USER_OBSGEOMS)
      VLIDORT_ModIn%MUserVal%TS_user_vzangles_input(1:N_USER_OBSGEOMS) = user_vzangles(1:N_USER_OBSGEOMS)

      VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS = N_USER_OBSGEOMS

      VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES = N_USER_OBSGEOMS
      VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS = N_USER_OBSGEOMS
      VLIDORT_ModIn%MSunrays%TS_N_SZANGLES   = N_USER_OBSGEOMS

      ! !update driver variable input
      n_szas = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
      n_vzas = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
      n_azms = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
      szas = VLIDORT_ModIn%MSunrays%TS_szangles
      azms = VLIDORT_ModIn%MUserVal%TS_user_relazms
      vzas = VLIDORT_ModIn%MUserVal%TS_user_vzangles_input
      n_geometries = VLIDORT_ModIn%MUserVal%TS_N_USER_OBSGEOMS
      Out_szas(1) = VLIDORT_ModIn%MSunrays%TS_szangles(1)
      Out_azms(1) = VLIDORT_ModIn%MUserVal%TS_user_relazms(1)
      Out_vzas(1) = VLIDORT_ModIn%MUserVal%TS_user_vzangles_input(1)


      nstreams = VLIDORT_FixIn%Cont%TS_NSTREAMS
      ncoeffs = 2*nstreams

      DO_upwelling = VLIDORT_FixIn%Bool%TS_DO_UPWELLING
      DO_dnwelling = VLIDORT_FixIn%Bool%TS_DO_DNWELLING
      DO_Sunlight    = VLIDORT_ModIn%MBool%TS_DO_SOLAR_SOURCES  ! default for natural light


      ! initiate
      DOublet_Offsets = 0
      Lattice_Offsets = 0

      ! Angles for VFZMAT
      IF ( DO_ObsGeoms ) THEN
        n_geoms = 1       ! 3 -- I think it should be 1 not 3, Xiong, 7/15/2022
        DO v = 1, n_geoms
          obsgeoms(v,1:3) =VLIDORT_ModIn%MUserVal%TS_USER_OBSGEOMS_INPUT(v,1:3)
        ENDDO

      ELSE IF (DO_DOublets) THEN
        DO ib = 1, n_szas
          DOublet_offsets(ib) = n_vzas*(ib-1)
        ENDDO
        n_geoms = N_szas * n_vzas     ! added 8/19/2022
      ELSE
        DO ib = 1, n_szas
          DO um = 1, n_vzas
            Lattice_offsets(ib,um) = n_azms*n_vzas*(ib-1) + n_azms*(um-1)
          ENDDO
        ENDDO
        n_geoms = N_szas * n_vzas *n_azms    ! added 8/19/2022
      ENDIF


      IF (VLIDORT_ModIn%MBool%TS_DO_DOUBLET_GEOMETRY) THEN
        n_geoms = VLIDORT_ModIn%MUserVal%TS_N_USER_DOUBLETS
        PRINT*, "DO DOUBLET GEOMETRY", n_geoms
      ENDIF



      ! 2.5: surface reflectivity (spectral depended)
      IF ( VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE .AND. LAMBERTIAN_ALBEDO .GE. 1 .AND. &
           LAMBERTIAN_ALBEDO .LE. 10) THEN
        IF (allocated(SfcRef_Usr) ) deallocate(SfcRef_Usr)
        ALLOCATE ( SfcRef_Usr(1:ndat) )
        IF (allocated(lambdas_usr) ) deallocate(lambdas_usr)
        ALLOCATE ( lambdas_usr(1:ndat) )
        ! Interpolate SfcRef at User grid
        indx_lambertian = LAMBERTIAN_ALBEDO
        SfcRef1D(:) = SfcRef(:,indx_lambertian)
        lambdas_usr(:) = lambdas(1:ndat) * 0.001        ! make units as um, same as surface reflectance data
        CALL interp_linear(nwav_sfcref, lambSfcRef, SfcRef1D, ndat, lambdas_usr,SfcRef_Usr )
      ENDIF



      ! 2.6: setting and calculating Fmat using VFZMAT
      IF (Allocated(FMatCoeffs)) deAllocate(FMatCoeffs)
      ALLOCATE (FMatCoeffs(0:ncoeffs,6))

      IF (Allocated(FMatCoeffs_intp)) deAllocate(FMatCoeffs_intp)
      ALLOCATE (FMatCoeffs_intp(MAXLAYERS,0:ncoeffs,6))

      ! CALL pre-vfzmat under sample loop
      IF (.not. quadrature_gen_flag ) THEN

        ALLOCATE(QuadAngles(n_Quadangles),QuadCosines(n_Quadangles),QuadWeights(n_Quadangles))
        ALLOCATE (GSF_P00_Saved (n_Quadangles,0:ncoeffs),GSF_P02_Saved(n_Quadangles,0:ncoeffs))
        ALLOCATE(GSF_P2p2_Saved(n_Quadangles,0:ncoeffs),GSF_P2m2_Saved(n_Quadangles,0:ncoeffs))

        CALL vfzmat_Pre_Master &
             ( max_geometries, max_szangles, max_user_vzangles, max_user_relazms,deg_to_rad, & ! Input Dimensions (VLIDORT)
             DO_upwelling, DO_dnwelling, DO_ObsGeoms, DO_DOublets,ncoeffs, nstokes,n_QuadAngles, & ! Input Flags and Control
             n_geoms, n_szas, n_vzas, n_azms, Lattice_offsets, DOublet_offsets, szas, vzas, azms, obsgeoms,& ! Input Geometries
             C1_up, S1_up, C2_up, S2_up, C1_dn, S1_dn, C2_dn, S2_dn,& ! Output rotation angles
             COSSCAT_up, COSSCAT_dn, QuadAngles, QuadCosines, QuadWeights,& ! Output Scatcosines and Quadrature
             GSF_P00_Saved, GSF_P02_Saved, GSF_P2p2_Saved, GSF_P2m2_Saved )! Output Saved GSFs

        quadrature_gen_flag = .TRUE.
      ENDIF

      ! initilization
      OutFmatrices_up =  zero
      OutFmatrices_dn =   zero
      FMatCoeffs =    zero
      FMatCoeffs_mode =  zero
      FMatCoeffs_tot =   zero


      CALL cpu_time(t_start)

      ! 2.7: wave grid loop for vlidort CALL
      ncount = 0
      nstep =  1
      st_wvnum_indx = 1
      end_wvnum_indx = ndat
      PRINT*,'ndat',ndat,st_wvnum_indx,end_wvnum_indx
      DO w = st_wvnum_indx, end_wvnum_indx, nstep
        PRINT*,'wavenumber',w,wavenums(w)
        ! initilization
        deltau_vert_input    = zero
        omega_total_input    = zero
        greekmat_total_input = zero
        Fmatrix_up           = zero
        Fmatrix_dn           = zero
        RayFmatrices_up = zero
        RayFmatrices_dn  = zero
        RayCoeffs  = zero
        OutFmat_up_intp =  zero
        OutFmat_dn_intp =   zero
        FMatCoeffs_intp =  zero


        ! reyleigh VFZMAT
        depol = rayleigh_depol(w)

        CALL vfzmat_Rayleigh &
             ( max_geometries, max_szangles, max_user_vzangles, max_user_relazms,deg_to_rad, & !Input  Dimensions (VLIDORT)
             DO_upwelling, DO_dnwelling, DO_ObsGeoms, DO_DOublets,DO_Sunlight,              & !Input  Flags
             nstokes, n_geoms, n_szas, n_vzas, n_azms,                          & !Input  Numbers
             Lattice_offsets, DOublet_offsets, szas, vzas, azms, obsgeoms, DEPOL,            & !Input  Geometries + Depol
             RayFmatrices_up, RayFmatrices_dn, RayZmatrices_up, RayZmatrices_dn,RayCoeffs )

        ! intepolate Fmat at current wavelen
        wavelen = lambdas(w)*1.e-3        ! 11/10/2022


        IF (LAMBERTIAN_ALBEDO .GT. 1 ) THEN
          LAMBERTIAN_ALBEDO =  SfcRef_Usr(w)
        END IF

        IF (DO_Planetary) LAMBERTIAN_ALBEDO = 0.0

        ! assign all necessary vlidrot input for vlidort run
        IF (.not. DO_AEROSOL) flag_aer = 0
        CALL Assign_LayOpt_to_Vlidort(w,nlayers,ndat,nstokes,n_geoms,MAXLAYERS, MAX_GEOMETRIES,    & ! INPUTS
             ncoeffs,MAXMOMENTS_INPUT,MAXSTOKES_SQ, airClm_vld,  RAYLEIGH_XSEC,&        ! INPUTS
             RayCoeffs, RayFmatrices_up,RayFmatrices_dn, &                          ! INPUTS
             flag_aer,aerext_intp,aersca_intp, &                          ! INPUTS
             odIOTot_vld,                      &                          ! INPUTS
             OutFmat_up_intp,OutFmat_dn_intp,FMatCoeffs_intp,&              ! INPUTS
             deltau_vert_input,omega_total_input,&          ! OUTPUTS
             GREEKMAT_TOTAL_INPUT,fmatrix_up,fmatrix_dn)  ! OUTPUTS


        ! Copy to vlidort type-structure input
        VLIDORT_ModIn%MCont%TS_ngreek_moments_input   = ncoeffs  !NGREEK_MOMENTS_INPUT
        VLIDORT_FixIn%Chapman%TS_height_grid(0:nlayers) = htgrid_vld(0:nlayers)
        VLIDORT_FixIn%Optical%TS_lambertian_albeDO    = LAMBERTIAN_ALBEDO
        VLIDORT_FixIn%Optical%TS_deltau_vert_input    = deltau_vert_input
        VLIDORT_FixIn%Optical%TS_greekmat_total_input = greekmat_total_input
        VLIDORT_ModIn%MOptical%TS_omega_total_input   = omega_total_input
        VLIDORT_FixIn%Optical%TS_FMATRIX_UP = Fmatrix_up
        VLIDORT_FixIn%Optical%TS_FMATRIX_DN = Fmatrix_dn
        PRINT*,'lambertian', LAMBERTIAN_ALBEDO
        ! vlidort master CALL
        DO_debug_input = .FALSE.
        CALL VLIDORT_MASTER (  DO_debug_input,&
             VLIDORT_FixIn, &
             VLIDORT_ModIn, &
             VLIDORT_Sup,   &
             VLIDORT_Out )

        OPENFILEFLAG = .FALSE.
        CALL VLIDORT_WRITE_STATUS ( &
             'V2p8p3_VLIDORT_Execution.log', &
             VLIDORT_ERRUNIT, OPENFILEFLAG,VLIDORT_Out%Status )
        checkrun      = VLIDORT_Out%Status%TS_STATUS_INPUTCHECK

        ncount_usr = ncount_usr + 1
        DO v = 1, n_geometries
          Iup_output2_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_stokes(1,v,1,1)
          Qup_output2_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_stokes(1,v,2,1)
          Uup_output2_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_stokes(1,v,3,1)
          IF (DO_Planetary) THEN
            BTRANS_I_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(1,v)
            BTRANS_Q_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(2,v)
            BTRANS_U_usr(v,ncount_usr) = VLIDORT_Out%Main%TS_PLANETARY_TRANSTERM(3,v)
            IF (v .EQ. 1) BSPHER_usr(ncount_usr) = VLIDORT_Out%Main%TS_PLANETARY_SBTERM
          ENDIF
        ENDDO
        PRINT*,'I usr',Iup_output2_usr(1,ncount_usr)
        wnum_usr(ncount_usr) = wavenums(w)
        PRINT*,'ncount_usr',ncount_usr
      ENDDO
      ! end wave grid loop

      CALL cpu_time(t_end)
      PRINT*, 'cal time vlidort w loop', t_end-t_start


      999  CONTINUE
      2017 print*, " ###### end in band # ", iband


    ENDDO
    ! END  B A N D - L O O P ************************
    CALL cpu_time(e2)

    ! --SECTION 3: output at user wave grid
    IF (VLIDORT_CRTM) THEN
      RTSolution(1:ncount_usr)%Stokes(1)  = Iup_output2_usr(1,1:ncount_usr)
      RTSolution(1:ncount_usr)%Stokes(2)  = Qup_output2_usr(1,1:ncount_usr)
      RTSolution(1:ncount_usr)%Stokes(3)  = Uup_output2_usr(1,1:ncount_usr)
      RTSolution(1:ncount_usr)%Stokes(4)  = 0
      RTSolution(1:ncount_usr)%radiance  = Iup_output2_usr(1,1:ncount_usr)
      N_USER_Channels = ncount_usr
      Error_Status = ZERO
      PRINT*,'tot time',e2-e1
      RETURN 
    ENDIF

    IF ( DO_Planetary) THEN
      wnStart = wnum_usr(1)
      PRINT*,'wnStart',wnStart
      CALL Write_3components (date, Nstokes, ncount_usr, nstep,n_geometries,nProf_s,nameProf,LAMBERTIAN_ALBEDO,&
           Out_szas(1:n_geometries), Out_vzas(1:n_geometries), Out_azms(1:n_geometries),maxtypeaer,   &
           wnStart,wnEnd,wresol,nstreams,nprof,psfc_out,VLIDORT_ModIn,st_wvnum_indx, end_wvnum_indx,wnum_usr(1:ncount_usr),&
           BTRANS_I_usr(1:n_geometries,1:ncount_usr), BTRANS_Q_usr(1:n_geometries,1:ncount_usr),&
           BTRANS_U_usr(1:n_geometries,1:ncount_usr), BSPHER_usr(1:ncount_usr), &
           Iup_output2_usr(1:n_geometries,1:ncount_usr),Qup_output2_usr(1:n_geometries,1:ncount_usr),&
           Uup_output2_usr(1:n_geometries,1:ncount_usr))

    ELSE
      CALL Write_IQU (date,Nstokes, ncount_usr, nstep,n_geometries,nProf_s,nameProf,LAMBERTIAN_ALBEDO,&
           Out_szas(1:n_geometries), Out_vzas(1:n_geometries), Out_azms(1:n_geometries),maxtypeaer,   &
           wnStart,wnEnd,wresol,nstreams,nprof,psfc_out,VLIDORT_ModIn,st_wvnum_indx, end_wvnum_indx,&
           Iup_output2_usr(1:n_geometries,1:ncount_usr),Qup_output2_usr(1:n_geometries,1:ncount_usr),&
           Uup_output2_usr(1:n_geometries,1:ncount_usr))
    ENDIF
    PRINT*,'tot time',e2-e1

    ! #############################################
    ! MODIFIED LUT COUPLING DRIVER ENDED
    ! #############################################
    RETURN 

  END FUNCTION crtm_omps_simulator



  END MODULE crtm_lbl_simulator
