

program test_gas

 
      USE PCRTM_Sw_Aux_m
      USE GAS_OPT_m

      IMPLICIT NONE


!!  OD LUT a
    REAL,dimension(:,:), allocatable :: odIO_tot

    REAL,dimension(:,:,:),allocatable::Iup_output,Idn_output,Qup_output,Qdn_output,Uup_output,Udn_output
    REAL,dimension(:,:),allocatable::Iup_output2,Qup_output2,Uup_output2   
    REAL,dimension(:,:),allocatable::BTRANS_I,BTRANS_Q,BTRANS_U
    REAL,dimension(:), allocatable :: BSPHER  
    REAL,dimension(:),  allocatable::SolTran_output
    REAL,dimension(:,:),allocatable::Iup_output2_usr,Qup_output2_usr,Uup_output2_usr 
    REAL,dimension(:,:),allocatable::BTRANS_I_usr,BTRANS_Q_usr,BTRANS_U_usr
    REAL,dimension(:), allocatable :: BSPHER_usr
    REAL(8),dimension(:), allocatable :: wnum_usr

    INTEGER                          :: ngases
    INTEGER                          :: ndat,w,checkrun

    REAL(8), dimension(:),allocatable :: wavenums, lambdas
    INTEGER, dimension(:),allocatable :: waveindex

    REAL(8), dimension(:), allocatable :: rayleigh_xsec,rayleigh_depol

    REAL(8)                          :: wresol
    REAL(8), dimension(1:MAX_SZANGLES)          :: szangles,user_relazms,user_vzangles

    REAL                             :: t_start,t_end !cpu time


!! ATM profile      
    INTEGER                          :: nProf, nprof_st,nprof_end

!! Filename related
    CHARACTER(LEN=300)               :: profile_data_filename,nameProf


!

! AFGL profile variable
    INTEGER  :: ATM_TYPE, j
    INTEGER  :: iband           !index of LUT band


    INTEGER :: RanDOm_Seed_times, nstep

    CHARACTER(LEN=80) :: path0, path1
    REAL(8)                          :: wnStart, wnEnd
    INTEGER,dimension(1:maxgases)          :: gasIdBand !gasIDall! !hitran ID of all abs. gases in the band

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
    INTEGER :: nlevel_atm,ngas_atm


    REAL, allocatable :: Wgas_usr (:,:)
    REAL, allocatable :: gph_usr (:)
    REAL, allocatable :: airAmt_usr (:)
    REAL*8,allocatable:: airClm_vld(:)        
    REAL*8,allocatable:: htgrid_vld(:)
    REAL*8,allocatable:: odIOTot_vld(:,:)


    INTEGER             ::  st_band, end_band, N_DOublets, N_USER_OBSGEOMS
    CHARACTER(256) :: Buffer
    LOGICAL :: Test_hotspot, DO_Scalar_only
    INTEGER ::  st_wvnum_indx, end_wvnum_indx      

    CHARACTER(8)         :: date
    CHARACTER(10)        :: time
    CHARACTER(5)         :: zone
    INTEGER,dimension(8) :: values

! add user set stream, default =6 , 11/20/2022
    LOGICAL :: DO_User_Stream
    INTEGER :: userset_stream = 6

    CHARACTER(LEN=80) ::PresGridFile
    LOGICAL :: DO_mapgrid,DO_usrgrid,DO_stdgrid
    LOGICAL :: DO_opac, DO_tamu

    INTEGER :: ndat_lut,nwv_usr

    REAL*8, dimension(:),allocatable :: wv_usr
    REAL*8, dimension(:,:),allocatable :: wv_band
    INTEGER, dimension(:,:),allocatable :: wvidx_band

    INTEGER:: ndat_band (nband_lut)

    CHARACTER(LEN=300)  :: file_usrwv

    INTEGER   :: nwav_st,nwav_end
    INTEGER   :: Ncount_usr


   !working directory
   path0='/data/jcsda3/mchen/VLIDORT/vlidort_NOAA_V1'
   path1=trim(path0)//'/PCRTM_VLIDORT_Main/'



   !1.2 read in user spectral grid:  
   file_usrwv=trim(path0)//'/PCRTM_VLIDORT_Config/user_wav_inputs.dat'
   CALL read_usr_wv(path0,file_usrwv,nwv_usr,wv_usr) 
   !setting band based on wv_usr
   print*, wv_usr
   CALL Set_Usr_wvgrid(path0,nwv_usr,wv_usr,&                                !input
                       st_band,end_band,wv_band,wvidx_band,ndat_band)         !output
   IF(end_band .lt. st_band) THEN
    Print*, "Error in Band setting !"
    stop
   endIF
   print*,'DO LUT band from Band', st_band,' to band', end_band


end program test_gas


