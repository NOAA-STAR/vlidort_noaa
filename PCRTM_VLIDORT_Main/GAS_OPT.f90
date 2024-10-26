MODULE GAS_OPT_m

! By Ming Zhao, 6/8/2022

! ###############################################################
! #                                                             
! # SUBROUTINEs in this Module:
!                                  
! # Set_Gas_Band
! # gen_Usr_Pres_grid                                                            
! # read_Prefix_Pres_Grid
! # Read_US_Stand_Prof
! # Read_AFGL_Prof
! # Read_Merra_Prof
! # Interpolate_User_grid
! # Calc_Layer_Amount
! # Calc_Gas_OD
! # Get_MerraProf_list
! ###############################################################

USE VLIDORT_PARS_m, only: MAXLAYERS

  INTEGER, parameter :: ngas_afgl = 28                 !currently has 28 afglgas(hirtranordercopy from lihui)
  INTEGER, parameter :: nlev_afgl = 50,nlay_afgl=49
  INTEGER, parameter :: nlev_lut =101,nlay_lut = 100
  INTEGER, parameter :: nband_lut = 36                 !bands of LUT
  INTEGER, parameter :: ngas_lut = 45                  !max gases including some blanks
  INTEGER, parameter :: nMerra = 9690
  ! SETUP MAX DIMENSIONS
  INTEGER, parameter :: maxlambdas = 250001            !250001 !the max num ofdatapoint from all  lut block
  INTEGER, parameter :: maxgases   = 30                !max abs gases in one lut block
  INTEGER, parameter :: maxUsrWv   = 200000
  INTEGER, parameter :: maxAtmlev = 101                !max level of atm profile
  !c -- Parameters
  REAL, parameter :: r_equator = 6.378388e+06          ! Earth radius at equator (m)
  REAL, parameter :: r_polar   = 6.356911e+06          ! Earth radius at pole (m)
  REAL, parameter :: r_avg     = 0.5*(r_equator+r_polar)
  REAL, parameter :: g_sfc = 9.80665                   ! Gravity at surface (m/s/s)
  REAL, parameter :: rho_ref = 1.2027e-12              ! Reference air "density"
  REAL, parameter :: mw_dryair = 28.97                 ! Molec. wgt. of dry air(g/mol)
  REAL, parameter :: R_gas = 8.3143                    ! Ideal gas constant(J/mole/K)
  REAL, parameter :: R_air = 0.9975*R_gas              ! Gas constant for air (worst case)
  REAL, parameter :: AVOGAD = 6.02214199E+23
  REAL, parameter :: rair=2.12157985e+22



CONTAINS


SUBROUTINE Get_MerraProf_list(path0,logprof)

   implicit none

   CHARACTER(LEN=80), intent(in) :: path0
   CHARACTER(LEN=300),dimension(1:nMerra) ,intent(out):: logprof 
   INTEGER :: IFILE, i   
  
    CALL GETLUN(IFILE)

    open(unit = IFILE, file=trim(path0)//'/Profile_merra/Profile_Merra_list.txt',&
         action='read', status ='old')
    rewind (IFILE)
    DO i = 1, nMerra
       read (IFILE,*) logprof(i)
    ENDDO
    CLOSE(IFILE)

END SUBROUTINE Get_MerraProf_list



SUBROUTINE Set_Gas_Band(path0,iband,&                                !input
                        vUsr1, vUsr2, wresol,ndat, ngases,gasIdBand) !output

   implicit none

  Character(LEN=80),intent(in) :: path0
  Character(len=100)           :: fnameGasTable
  REAL*8, dimension(nband_lut) :: v1_lut,v2_lut,res_lut
  INTEGER  :: iband, i,j,IFile
  INTEGER  :: gasFlg(ngas_lut, nband_lut)
  INTEGER,intent(out)  :: ngases,ndat
  REAL*8, intent(out)  :: vUsr1, vUsr2,wresol
  INTEGER,intent(out)  :: gasIdBand(1:maxgases)
  


  fnameGasTable=trim(path0)//"/PCRTM_VLIDORT_Config/gasTable"
  CALL GETLUN(IFILE)
  open ( unit = IFile, file = trim(fnameGasTable), action = 'read', status ='old')
  rewind ( IFile )
  read (IFile, *)
  DO i = 1, nband_lut
     read (IFile, *) v1_lut(i), v2_lut(i), res_lut(i), (gasFlg(j,i), j=1,ngas_lut)
  ENDDO
  close ( IFile )

  !setting start/END point of the spectral
  vUsr1 = v1_lut(iband)
  vUsr2 = v2_lut(iband)
  wresol = res_lut(iband)!4.0D-3
  ndat = ceiling((vUsr2-vUsr1)/wresol)+1
  ngases = sum (gasFlg(:,iband))
  
  !assign gasID (HITRAN ID) for abs gases
  j=1
  DO i = 1, ngas_lut
     IF (gasFlg(i,iband) == 1) then
        gasIdBand(j) = i !1:H2O 2:CO2 ...44:O4 45:BRO
        j = j +1
     ENDIF
  ENDDO



END SUBROUTINE Set_Gas_Band


SUBROUTINE read_usr_wv(path0,filename,nwv_usr, wn_usr)

   implicit none

 CHARACTER(LEN=80) , intent(in)   :: path0


 REAL*8  :: st_wv,END_wv,reso
 INTEGER  :: IFILE, i
 Character(LEN=300) ,intent(inout)  :: filename
 REAL*8,             intent(out), dimension(:),allocatable  :: wn_usr
 INTEGER,            intent(out)  :: nwv_usr

! local INTEGER
 INTEGER :: nwv_usr0
 REAL*8  :: wl_usr

    CALL GETLUN(IFILE)
    open ( unit = IFILE, file = filename,action='read',status='old')
    rewind ( IFILE )
    DO i = 1, 5
      read(IFILE,*)
    ENDDO
    read(IFILE,*) nwv_usr0

    DO i = 1, 2
      read(IFILE,*)
    ENDDO
    read(IFILE,*) st_wv
    read(IFILE,*) END_wv

    DO i = 1, 2
      read(IFILE,*)
    ENDDO
    read(IFILE,*) reso
    IF (nwv_usr0 < 0) then
         nwv_usr = (END_wv-st_wv) / reso + 1
    ELSE
        nwv_usr = nwv_usr0
    ENDIF

    IF (allocated(wn_usr) ) deallocate(wn_usr)
    allocate(wn_usr(1:nwv_usr) )

    IF (nwv_usr0 > 0) then
       DO i = 1, 2
         read(IFILE,*)
       ENDDO
       DO i = 1, nwv_usr
          read(IFILE,*) wl_usr
          wn_usr(nwv_usr-i+1) = 1.D7/wl_usr
       END DO
    ELSE
       DO i =1, nwv_usr
         wl_usr  = st_wv + (i-1)*reso
         wn_usr (nwv_usr-i+1) = 1.D7/wl_usr
       ENDDO
    ENDIF

    close ( IFILE )
End SUBROUTINE read_usr_wv



!SUBROUTINE read_usr_wv(path0,filename,nwv_usr, wn_usr)
!
!   implicit none
!
! CHARACTER(LEN=80) , intent(in)   :: path0
!
!! REAL*8  :: wl_usr(maxUsrWv)
!
! REAL*8  :: st_wv,END_wv,reso 
! INTEGER  :: IFILE, i
! Character(LEN=300) ,intent(inout)  :: filename
!! REAL*8,             intent(out)  :: wn_usr(maxUsrWv)
! REAL*8,             intent(out), dimension(:),allocatable  :: wn_usr
! INTEGER,            intent(out)  :: nwv_usr
!
!! local INTEGER
! INTEGER :: nwv_usr0
! REAL*8, dimension(:),allocatable  :: wl_usr
! logical :: DO_wn
!
!   !    filename='PCRTM_VLIDORT_Config/user_wav_inputs.dat'
!   print*,filename
!   CALL GETLUN(IFILE)
!   open ( unit = IFILE, file = filename,action='read',status='old')
!   rewind ( IFILE )
!   read(IFILE,*)
!   read(IFILE,*) DO_wn
!   read(IFILE,*)
!   read(IFILE,*) nwv_usr0
!   read(IFILE,*)
!   read(IFILE,*) st_wv
!   read(IFILE,*)
!   read(IFILE,*) END_wv
!   read(IFILE,*)
!   read(IFILE,*) reso
!   close ( IFILE )
!
!   nwv_usr = nwv_usr0 
!
!   IF (allocated(wn_usr) ) deallocate(wn_usr)
!   allocate(wn_usr(1:nwv_usr)) 
!
!   IF (allocated(wl_usr) ) deallocate(wl_usr)
!   allocate(wl_usr(1:nwv_usr) ) 
!
!   IF (DO_wn) then
!      DO i = 1, nwv_usr
!         wn_usr(i)  = st_wv + (i-1)*reso
!         wl_usr (nwv_usr-i+1) = 1.D7/wn_usr(i)
!      ENDDO
!   ELSE
!      DO i = 1, nwv_usr
!         wl_usr(i)  = st_wv + (i-1)*reso
!         wn_usr (nwv_usr-i+1) = 1.D7/wn_usr(i)
!      ENDDO
!   ENDIF
!
!    close ( IFILE )
!
!End SUBROUTINE read_usr_wv



SUBROUTINE Set_Usr_wvgrid(path0,nwv_usr,wv_usr,&                                !input
                         st_band,END_band, wv_band,wvidx_band,ndat_band)        !output
   implicit none

  Character(LEN=80),intent(in) :: path0
  INTEGER,intent(in)  :: nwv_usr
  REAL*8 ,intent(in)  :: wv_usr(nwv_usr)

  Character(len=100)           :: fnameGasTable
  REAL*8, dimension(nband_lut) :: v1_lut,v2_lut,res_lut
  INTEGER  ::  i,j,k,n,IFILE,st_band,END_band
  INTEGER  :: gasFlg(ngas_lut, nband_lut)
  INTEGER  :: n1
  REAL*8   :: vUsr1, vUsr2,wresol

! local variable
  REAL(8), dimension(:),allocatable :: wavenums_1
  INTEGER:: ndat_band (nband_lut) 

  REAL*8, dimension(:,:),allocatable :: wv_band
  INTEGER, dimension(:,:),allocatable:: wvidx_band

  !fnameGasTable=trim(path0)//"/GasTable/gasTable"
  fnameGasTable=trim(path0)//"/PCRTM_VLIDORT_Config/gasTable"
  CALL GETLUN(IFILE)
  open ( unit = IFILE, file = trim(fnameGasTable), action = 'read', status ='old')
  rewind ( IFILE )
  read (IFILE, *)
  DO i = 1, nband_lut
     read (IFILE, *) v1_lut(i), v2_lut(i), res_lut(i), (gasFlg(j,i), j=1,ngas_lut)
  ENDDO
  close ( IFILE )

  DO i = 1, nband_lut-1
     IF ( wv_usr(1) .ge. v1_lut(i) .and. wv_usr(1) .lt. v1_lut(i+1)  ) then
         st_band = i
     ENDIF
  ENDDO


  DO i = 1, nband_lut-1
     IF ( wv_usr(nwv_usr) .ge. v1_lut(i) .and. wv_usr(nwv_usr) .lt. v1_lut(i+1)  ) then
         END_band = i
     ENDIF
  ENDDO

      IF (allocated(wv_band) ) deallocate(wv_band)
       allocate(wv_band(1:nband_lut,1:nwv_usr) )

      IF (allocated(wvidx_band) ) deallocate(wvidx_band)
       allocate(wvidx_band(1:nband_lut,1:nwv_usr) )


  wv_band(:,:) = 0.0D0
  wvidx_band(:,:) = 0.0D0
  ndat_band(:) = 0 
  DO i = st_band, END_band
     vUsr1 = v1_lut(i)
     vUsr2 = v2_lut(i)
     wresol = res_lut(i)
     n1 = ceiling((vUsr2-vUsr1)/wresol)+1

     IF (allocated(wavenums_1) ) deallocate(wavenums_1)
     allocate( wavenums_1(n1) )

     DO n = 1, n1
        wavenums_1(n) = vUsr1 + dble(n-1) * wresol
     ENDDO

     k = 1
     DO j = 1, nwv_usr
      IF (wv_usr(j) .ge. vUsr1 .and. wv_usr(j) .lt. vUsr2 ) then
             wv_band(i,k) = wv_usr(j)
             k = k+1
        ENDIF 
     ENDDO     
     ndat_band(i) = k-1
  
     DO j = 1, ndat_band(i)
         DO n = 1, n1-1
            IF (wv_band(i,j ) .ge. wavenums_1(n) .and. wv_band(i,j ) .lt. wavenums_1(n+1)) then
                wvidx_band(i,j) = n
            ENDIF 
         ENDDO
     ENDDO
       
  ENDDO

END SUBROUTINE Set_Usr_wvgrid





SUBROUTINE read_gasLUT_P(path0,pres_lut)

   implicit none

   CHARACTER(LEN=80) , intent(in)   :: path0
   Character(LEN=300)               :: filename

   INTEGER :: i,IFile
   INTEGER           :: level_index(1:nlev_LUT)
   REAL, intent(out) :: pres_lut(1:nlev_lut)

    filename =trim(path0)//'/PCRTM_VLIDORT_Config/GasLUT_PGrid.dat'
    CALL GETLUN(IFILE) 
    open ( unit = IFile, file = filename,action='read',status='old')
    rewind ( IFile )
    DO i = 1, nlev_LUT !reverse the grid: 1 is TOA
       read(IFile,*) level_index(nlev_LUT-i+1), pres_lut(nlev_LUT-i+1)
    END DO
    close ( IFile )

END SUBROUTINE read_gasLUT_P



SUBROUTINE read_Pres_Grid(path0,PresGridFile,& !input
                         nlevel_usr,pres_usr)  !output
 INTEGER :: i,IFile
 INTEGER :: level_index(1:maxlayers)
 CHARACTER(LEN=300) :: filename
 CHARACTER(LEN=80), intent(in)  :: path0,PresGridFile
 INTEGER, intent(out) :: nlevel_usr
 REAL   , intent(out) :: pres_usr(1:maxlayers)

    filename='/PCRTM_VLIDORT_Config/'//trim(PresGridFile)     
    CALL GETLUN(IFILE) 
    open ( unit = IFile, file = trim(path0)//filename,action='read',status='old')
    rewind ( IFile )
    read(IFile,*) nlevel_usr
    DO i = 1, nlevel_usr !reverse the grid: 1 is TOA
       read(IFile,*) level_index(i), pres_usr(i)
    END DO
    close ( IFile )

End SUBROUTINE read_Pres_Grid



SUBROUTINE Read_AFGL_Prof(path0,ATM_type, Psfc_in,&                                     !input 
                         nlevel_atm,ngas_atm,pres_atm,temp_atm,vmr_atm,psfc_out,&
                         nameProf,profile_data_filename)        !output 

   implicit none

 Character(LEN=80),intent(in) :: path0
 INTEGER          ,intent(in) :: atm_type
 REAL,             intent(inout) :: pSfc_in
 INTEGER ::  iafgl, jafgl,i,IFILE
 INTEGER, intent(out) :: nlevel_atm,ngas_atm 
 
 !local
 Character(Len=1)        :: AfglMod_s
 REAL,dimension(1:nlev_afgl)    :: level_afgl,alt_afgl
 REAL*8  :: afgl_grid   ( 0:maxAtmLev )

 REAL,                       intent(out)  :: vmr_atm(1:maxAtmLev, 1:ngas_lut) !vmr dirctily from atm profile 
 REAL, dimension(maxAtmLev), intent(out)  :: pres_atm,temp_atm
 REAL,                       intent(out)  :: psfc_out
 CHARACTER(LEN=300),         intent(out)  :: profile_data_filename,nameProf

     ngas_atm = ngas_afgl
     nlevel_atm =50
     ! AFGL ALT (All afgl model use the same height grid)
      alt_afgl=(/ 0.0,      1.0,       2.0,       3.0,       4.0,        &
                 5.0,       6.0,       7.0,       8.0,       9.0,        &
                 10.0,      11.0,      12.0,      13.0,      14.0,       &
                 15.0,      16.0,      17.0,      18.0,      19.0,       &
                 20.0,      21.0,      22.0,      23.0,      24.0,       &
                 25.0,      27.5,      30.0,      32.5,      35.0,       &
                 37.5,      40.0,      42.5,      45.0,      47.5,       &
                 50.0,      55.0,      60.0,      65.0,      70.0,       &
                 75.0,      80.0,      85.0,      90.0,      95.0,       &
                100.0,     105.0,     110.0,     115.0,     120.0/)

    write(AfglMod_s, '(i0)') ATM_TYPE
    nameProf= 'afgl_model'//trim(AfglMod_s)//'.dat'  
    profile_data_filename='data_and_control/Profile_afgl/afgl_model'//trim(AfglMod_s)//'.dat'

    CALL GETLUN(IFILE)
    open ( unit = IFILE, file = trim(path0)//profile_data_filename, action='read',status='old')
    rewind ( IFILE )
    DO iafgl = 1, nlev_afgl !reverse the grid: 1 is TOA
       read(IFILE,*) level_afgl(iafgl),pres_atm(nlev_afgl-iafgl+1),temp_atm(nlev_afgl-iafgl+1),(vmr_atm(nlev_afgl-iafgl+1,jafgl),jafgl=1,ngas_afgl)
    END DO
    close ( IFILE )

    print*,'use AFGL profile: ', profile_data_filename

    !set TOA height and 'reverse' alt_afgl to afgl_grid (to match the vliDOrt
    !grid)
    DO  iafgl =  nlev_afgl, 1, -1
       afgl_grid(nlev_afgl-iafgl+1)= alt_afgl(iafgl)
    END DO

    !avoid user input psfc_in > lower atm press
    IF (Psfc_in .gt. pres_atm(nlev_afgl)) Psfc_in=pres_atm(nlev_afgl)
    
    psfc_out = Psfc_in


END SUBROUTINE Read_AFGL_Prof


SUBROUTINE  Read_Merra_Prof(profile_data_filename,psfc_in,&                            !input
                            nlevel_atm,ngas_atm,pres_atm,temp_atm,vmr_atm,psfc_out)
   implicit none

    INTEGER, parameter :: ngasMerra = 45
    INTEGER, parameter :: nlevMerra = 101
    INTEGER, parameter :: nlayMerra = 100
 
    !input
    !CHARACTER(LEN=80), intent(in)    :: path0
    !INTEGER,           intent(in)    :: MerraProfID                     
    REAL,              intent(inout) :: Psfc_in

    !local
    CHARACTER(LEN=300),dimension(1:nMerra) :: logprof !MERRA PROFILE
    REAL                 :: lat, lon, Psfc, ODCld, ExAer, U10, V10
    INTEGER, dimension(1:nlevMerra)        :: level_prof
    REAL, dimension(1:nlevMerra, 1:50)     :: data_prof !45vmr+pres temp cld ql qi
    INTEGER :: ii, jj, i, IFile
    CHARACTER(LEN=100)  :: dummy
    REAL        :: vmr_atm_CO2(1:maxAtmLev)

    !output
    INTEGER,                   intent(out) :: nlevel_atm, ngas_atm 
    REAL,dimension(maxAtmLev), intent(out)  :: pres_atm,temp_atm
    REAL,                      intent(out)  :: vmr_atm(1:maxAtmLev, 1:ngas_lut)
    REAL,                      intent(out)  :: psfc_out   !ming 01-03-2023, combined with merra's psfc, the 'REAL' surface press used in vliDOrt
    CHARACTER(LEN=300)                      :: profile_data_filename,nameProf

     ngas_atm = ngasMerra
     nlevel_atm =101


    !! READ PRESSURE TEMPERATURE AND VMR FROM PROFILE AND INVERSE COORDINATE TO FIT  VLIDORT COORDINATE (1 IS TOP)
    !!setting some constant values for atm profile
    CALL GETLUN(IFILE)
    open ( unit = IFILE, file = profile_data_filename, action = 'read', status ='old')
    rewind ( IFILE )

    read (IFILE, *) dummy,lat
    read (IFILE, *) dummy,lon
    read (IFILE, *) dummy,Psfc
    read (IFILE, *) dummy,ODCld
    read (IFILE, *) dummy,EXAer
    read (IFILE, *) dummy,U10
    read (IFILE, *) dummy,V10
    DO ii = 1, 50                 ! total 50 vairables (with 101 level)
       read (IFILE, *)
       DO jj = 1, nlevMerra       !reverse coordinate: 1 is TOA
          read (IFILE, *) level_prof(nlevMerra-jj+1), data_prof(nlevMerra-jj+1, ii)
       END DO
    END DO
    close ( IFILE )

    !!changed it to the following way, it works. DO not know WHY ??? 1/20/2022, Xiong
    DO jj = 1, nlevMerra
       pres_atm(jj) = data_prof(jj,1)
       temp_atm(jj)  = data_prof(jj,2)
       vmr_atm(jj,1:ngas_lut) = data_prof(jj,6:ngas_lut+5)
    END DO

    !Change merra CO2 coordinate upside DOwn (CO2 profile is opposite to othe gases)            
    DO jj = 1, nlevMerra
       vmr_atm_CO2(jj) =  vmr_atm(nlevMerra-jj+1, 2)
    ENDDO
    vmr_atm(1:nlevMerra,2) = vmr_atm_CO2(1:nlevMerra)


    print*,'pfc in merra',Psfc

    IF (Psfc_in .gt. 0 .and. Psfc_in .le. Psfc) Psfc = Psfc_in       !IF user input a surface pressure, use this value

    print*,'psfc re-selected',psfc

    !ming 01-03-2023 setting the surface pressure as psfc_out 
    psfc_out = Psfc

END SUBROUTINE Read_Merra_Prof


SUBROUTINE  Interpolate_LUT_grid(pres_usr,nlevel_usr,nlevel_atm,psfc_out,&                   !input
                                  pres_atm, temp_atm, vmr_atm,&                               !input
                                  temp_usr,vmr_usr,nlayers,temp_at_psfc)                      !output 
   implicit none

 !input
 REAL, dimension(nlev_lut), intent(in)  :: pres_usr                   !user provided pres grid
 INTEGER,                    intent(in)  :: nlevel_usr,nlevel_atm
 REAL,                       intent(in)  :: psfc_out                   !the REAL psfc used in vliDOrt
 REAL, dimension(maxAtmLev), intent(in)  :: temp_atm, pres_atm
 REAL,                    intent(inout)  :: vmr_atm(1:maxAtmLev, 1:ngas_lut)

 !local
 INTEGER :: l,ilev, i, j
 REAL  ::  A
 REAL    :: vmr_usr_dry(1:nlev_lut, 1:ngas_lut)

 !output
 REAL, intent(out)    :: vmr_usr(1:nlev_lut, 1:ngas_lut)
 REAL, intent(out)    :: temp_usr(1:nlev_lut)
 REAL, intent(out)    :: temp_at_psfc               !temp at psfc  
 INTEGER,intent(out)  :: nlayers


   !ming 01-17-2023 find number of layers for  pres grid
   !note: in driver CALL, pres_usr is user selected pres grid
   !      in the gas_OD_LayerAmt CALL, pres_usr is the lut pres grid 
   DO i = nlevel_usr,1, -1
        IF (pres_usr(i) < psfc_out) exit
   ENDDO
   nlayers=i

   IF (nlayers > nlevel_usr -1) then
     print*,'stop!! check the surface press, must < lowest usr pres grid'
     print*,'reset the psfc_in <= pres_usr(nlevel_usr)'
     print*,'nlayers=',nlayers,'> nlenlevel_usr-1=',nlevel_usr-1
     stop
   ENDIF


   !ming0802  
   DO i = 1,maxAtmLev       !MAXLAYERS
      DO j = 1, ngas_lut
         IF (vmr_atm(i,j) .lt. 1.e-36) vmr_atm(i,j) = 1.E-36
      ENDDO
   ENDDO

   !ming 01-03-2023 interpolating all user level, instead of nlayers+1 
   DO l=nlevel_usr,1,-1
   ilev = nlevel_atm -1
18        IF (pres_atm(ilev) .gt. pres_usr(l) .and. ilev .lt. (nlevel_atm)) then
            ilev=ilev-1
            goto 18
          ENDIF
          IF (ilev .eq. 0 ) ilev =1                !mingmod when user pres grid out of atm grid range 
          A=log(Pres_usr(l)/pres_atm(ilev)) / log(pres_atm(ilev+1)/pres_atm(ilev))
          temp_usr(l)=temp_atm(ilev)*(temp_atm(ilev+1)/temp_atm(ilev)) **A
          vmr_usr_dry(l,:)= vmr_atm(ilev,:)*(vmr_atm(ilev+1,:)/vmr_atm(ilev,:)) **A

   ENDDO

   DO i =1, nlevel_usr
     vmr_usr(i,1)= vmr_Usr_dry (i,1)
     DO j=2, ngas_lut
       vmr_usr(i,j)= vmr_Usr_dry(i,j)*(1 - vmr_usr(i,1)) !convt. dry vmr to wet vmr
     ENDDO
   ENDDO

!ming 01-03-2023 cal interpolated temp at psfc_out (use it for cal sfc gp-height)

   ilev = nlevel_atm -1
19 IF (pres_atm(ilev) .gt. psfc_out .and. ilev .lt. (nlevel_atm)) then
     ilev=ilev-1
     goto 19
   ENDIF
   A=log(psfc_out/pres_atm(ilev)) / log(pres_atm(ilev+1)/pres_atm(ilev))
   temp_at_psfc=temp_atm(ilev)*(temp_atm(ilev+1)/temp_atm(ilev)) **A


END SUBROUTINE Interpolate_LUT_grid



SUBROUTINE  Interpolate_User_grid(pres_usr,nlevel_usr,nlevel_atm,psfc_out,&                   !input
                                  pres_atm, temp_atm, vmr_atm,&                               !input
                                  temp_usr,vmr_usr,nlayers,temp_at_psfc)                      !output 
   implicit none

 !input
 REAL, dimension(MAXLAYERS), intent(in)  :: pres_usr                   !user provided pres grid
 INTEGER,                    intent(in)  :: nlevel_usr,nlevel_atm
 REAL,                       intent(in)  :: psfc_out                   !the REAL psfc used in vliDOrt
 REAL, dimension(maxAtmLev), intent(in)  :: temp_atm, pres_atm 
 REAL,                    intent(inout)  :: vmr_atm(1:maxAtmLev, 1:ngas_lut)

 !local
 INTEGER :: l,ilev, i, j
 REAL  ::  A
 REAL    :: vmr_usr_dry(1:MAXLAYERS, 1:ngas_lut)

 !output
 REAL, intent(out)    :: vmr_usr(1:MAXLAYERS, 1:ngas_lut)
 REAL, intent(out)    :: temp_usr(1:MAXLAYERS)
 REAL, intent(out)    :: temp_at_psfc               !temp at psfc  
 INTEGER,intent(out)  :: nlayers


   !01-17-2023 find number of layers for  pres grid
   !note: in driver CALL, pres_usr is user selected pres grid
   !      in the gas_OD_LayerAmt CALL, pres_usr is the lut pres grid 
   DO i = nlevel_usr,1, -1
        IF (pres_usr(i) < psfc_out) exit 
   ENDDO
   nlayers=i

   IF (nlayers > nlevel_usr -1) then
     print*,'stop!! check the surface press, must < lowest usr pres grid'
     print*,'reset the psfc_in <= pres_usr(nlevel_usr)'
     print*,'nlayers=',nlayers,'> nlenlevel_usr-1=',nlevel_usr-1
     stop
   ENDIF    


   !ming0802  
   DO i = 1,maxAtmLev       !MAXLAYERS
      DO j = 1, ngas_lut
         IF (vmr_atm(i,j) .lt. 1.e-36) vmr_atm(i,j) = 1.E-36
      ENDDO
   ENDDO

   !ming 01-03-2023 interpolating all user level, instead of nlayers+1 
   DO l=nlevel_usr,1,-1
   ilev = nlevel_atm -1
18        IF (pres_atm(ilev) .gt. pres_usr(l) .and. ilev .lt. (nlevel_atm)) then
            ilev=ilev-1
            goto 18
          ENDIF
          IF (ilev .eq. 0 ) ilev =1                !mingmod when user pres grid out of atm grid range 
          A=log(Pres_usr(l)/pres_atm(ilev)) / log(pres_atm(ilev+1)/pres_atm(ilev))
          temp_usr(l)=temp_atm(ilev)*(temp_atm(ilev+1)/temp_atm(ilev)) **A
          vmr_usr_dry(l,:)= vmr_atm(ilev,:)*(vmr_atm(ilev+1,:)/vmr_atm(ilev,:)) **A

   ENDDO

   DO i =1, nlevel_usr
     vmr_usr(i,1)= vmr_Usr_dry (i,1)
     DO j=2, ngas_lut
       vmr_usr(i,j)= vmr_Usr_dry(i,j)*(1 - vmr_usr(i,1)) !convt. dry vmr to wet vmr
     ENDDO
   ENDDO

!01-03-2023 cal interpolated temp at psfc_out (use it for cal sfc gp-height)

   ilev = nlevel_atm -1
19 IF (pres_atm(ilev) .gt. psfc_out .and. ilev .lt. (nlevel_atm)) then
     ilev=ilev-1
     goto 19
   ENDIF
   A=log(psfc_out/pres_atm(ilev)) / log(pres_atm(ilev+1)/pres_atm(ilev))
   temp_at_psfc=temp_atm(ilev)*(temp_atm(ilev+1)/temp_atm(ilev)) **A


END SUBROUTINE Interpolate_User_grid


SUBROUTINE Calc_GPH(nlayers,psfc_out,temp_at_psfc,P,T,vmr_usr,DO_icecld,& !input
                    height_T0, z)                                         !output   


     implicit none

      !input
      INTEGER,                   intent(in) ::  nlayers
      REAL,                      intent(in) ::  psfc_out,temp_at_psfc
      REAL,dimension(nlayers+1), intent(in)  ::  P,T
      REAL,                      intent(in)  :: vmr_usr(1:MAXLAYERS, 1:ngas_lut)

      !local variable
      REAL                                      ::  zsfc,  r_hgt,rho1,rho2,c_avg,p1,p2,dp,w1,w2
      REAL,   dimension(1:nlayers+1)           ::  W        !w is H2o g/kg
      REAL,   dimension(1:nlayers)             ::  Pavg
      REAL,   dimension(1:nlayers+1)           ::  mw_air
      REAL,   dimension(1:nlayers+1)           ::  g,c,rho_air
      REAL                                     ::  g_avg
      INTEGER                                  ::  i_dir, Lsfc,lObs, l,k
      REAL                                     :: w_avg

      REAL ::  Wgas(1:nlayers, 1:ngas_lut)   !!Wgas(1) is TOA
      REAL ::  airAmt(1:nlayers)             ![#/cm^2]
      REAL,   dimension(1:nlayers) ::  Tavg,H2Oavg



      !output
      REAL,   intent(out) ::  z(1:nlayers+1)
      !icetest
      REAL,   intent(out) :: height_T0
      logical             :: DO_icecld

      !calculate surface height [m]
      zsfc = log(1013.25/psfc_out) * 8.3144598 * temp_at_psfc / (0.0289644*9.80665 )![m]

      !cal h2o mixing ratio in g/kg from vmr
      w(1:nlayers+1) = vmr_usr(1:nlayers+1,1) * 1000.0 * 18.015/mw_dryair

      !cal geopotential height for user grid, [m]
      i_dir=1 ! levle(1) = p(top)
      CALL gphite( p, t, w, zsfc, nlayers+1, i_dir, z)

      Lsfc = nlayers+1
      lObs = 1
      DO l = lSfc,lObs,-1

        mw_air(l) = ( ( 1.0 - vmr_usr(l,1) ) * mw_dryair ) + (vmr_usr(l,1) * 18.015 )
        c(l) = 0.001 * mw_air(l) / R_gas  ! 0.001 factor for g->kg (kg*K/J)
        rho_air(l) = c(l) * p(l) / t(l) !kg/m3
        r_hgt = r_avg + z(l)                            ! (m)
        g(l) = g_sfc*(r_avg*r_avg)/(r_hgt*r_hgt)        ! (m/s^2)
      END DO

!c -- LAYER quantities --
      DO  l = lSfc-1,lObs,-1
        rho1 = rho_air(l)
        rho2 = rho_air(l+1)
        c_avg = ( (rho1*c(l)) + (rho2*c(l+1)) ) / ( rho1 + rho2 )
        Tavg(l) =( (rho1*t(l)) + (rho2*t(l+1)) ) / ( rho1 + rho2 )
        !ming 01-23-2023 cal H2O avg vmr
        H2Oavg(l) = ( (rho1*vmr_usr(l,1)) + (rho2*vmr_usr(l+1,1)) ) / ( rho1 + rho2 )
        p1 = p(l)
        p2 = p(l+1)
        dp = p2 - p1
        Pavg(l) = 0.5*(p1+p2)   ! density weigting is used (LBLRTM)
        g_avg = 0.5*(g(l)+g(l+1))

!c----Calculate air amount in a layer (similar to LBLRTM approach):
        airAmt(l)=dp*avogad*10./g_avg*2./(mw_air(l)+mw_air(l+1))

        !Calculate LAYER amounts for atmospheric gases:
        !mingmod amax
        DO k=1,ngas_lut
           w1 = amax1(vmr_usr(l,k),1.0e-20)
           w2 = amax1(vmr_usr(l+1,k),1.0e-20)

!c---  Kcarta approach:
           w_avg =  ( (rho1*w1) + (rho2*w2) ) / ( rho1+rho2 )
           !Wgas1(l,k) = w_avg * airAmt(l)
!c---  LBLRTM approach:
           CALL integAmt(w1,w2,p1,p2,g(l),g(l+1),mw_air(l),mw_air(l+1),Wgas(l,k))
        ENDDO
      ENDDO


!icetest ming: find the level that temp below 273k
    height_T0 = 0. 
    IF (DO_icecld) then
      IF (temp_at_psfc .lt. 273.0) then
         height_T0 = z(nlayers+1) 
      ELSE  
        DO  l = nlayers, 1, -1
          IF (T(l) .lt. 273.0) exit
        ENDDO  
          height_T0 = z(l)
      ENDIF
      height_T0 = height_T0 /1000. ![km]
    ENDIF
END SUBROUTINE Calc_GPH



SUBROUTINE Calc_Layer_Amount(nlayers,psfc_out,temp_at_psfc,P,T,vmr_usr,& !input
                             Wgas,z,airAmt,Tavg,H2Oavg)                  !output   


     implicit none

      !input
      INTEGER,                   intent(in) ::  nlayers
      REAL,                      intent(in) ::  psfc_out,temp_at_psfc
      REAL,dimension(nlayers+1), intent(in)  ::  P,T 
      REAL,                      intent(in)  :: vmr_usr(1:nlev_lut, 1:ngas_lut)

      !local variable
      REAL                                      ::  zsfc,  r_hgt,rho1,rho2,c_avg,p1,p2,dp,w1,w2
      REAL,   dimension(1:nlayers+1)           ::  W        !w is H2o g/kg
      REAL,   dimension(1:nlayers)             ::  Pavg       
      REAL,   dimension(1:nlayers+1)           ::  mw_air
      REAL,   dimension(1:nlayers+1)           ::  g,c,rho_air
      REAL                                     ::  g_avg
      INTEGER                                  ::  i_dir, Lsfc,lObs, l,k 
      REAL                                     :: w_avg
   
      !output
      REAL,   intent(out) ::  z(1:nlayers+1)
      REAL,   intent(out) ::  Wgas(1:nlayers, 1:ngas_lut)   !!Wgas(1) is TOA
      REAL,   intent(out) ::  airAmt(1:nlayers)             ![#/cm^2]
      REAL,   dimension(1:nlayers),intent(out) ::  Tavg,H2Oavg



      !calculate surface height [m]
      zsfc = log(1013.25/psfc_out) * 8.3144598 * temp_at_psfc / (0.0289644*9.80665 )![m]

      !cal h2o mixing ratio in g/kg from vmr
      w(1:nlayers+1) = vmr_usr(1:nlayers+1,1) * 1000.0 * 18.015/mw_dryair

      !cal geopotential height for user grid, [m]
      i_dir=1 ! levle(1) = p(top)
      CALL gphite( p, t, w, zsfc, nlayers+1, i_dir, z) 

      !mingmod setting lsfc=nlev lobs=1, it will cal all layer's column density
      !        IF the psfc_in is <  the pres_atm(nlev), it should use the
      !        'nlayer' to select the layers above the psfc_in
      Lsfc = nlayers+1
      lObs = 1
      DO l = lSfc,lObs,-1

        mw_air(l) = ( ( 1.0 - vmr_usr(l,1) ) * mw_dryair ) + (vmr_usr(l,1) * 18.015 )
        c(l) = 0.001 * mw_air(l) / R_gas  ! 0.001 factor for g->kg (kg*K/J)
        rho_air(l) = c(l) * p(l) / t(l) !kg/m3
        r_hgt = r_avg + z(l)                            ! (m)
        g(l) = g_sfc*(r_avg*r_avg)/(r_hgt*r_hgt)        ! (m/s^2)
      END DO

!c -- LAYER quantities --
      DO  l = lSfc-1,lObs,-1
        rho1 = rho_air(l)
        rho2 = rho_air(l+1)
        c_avg = ( (rho1*c(l)) + (rho2*c(l+1)) ) / ( rho1 + rho2 )
        Tavg(l) =( (rho1*t(l)) + (rho2*t(l+1)) ) / ( rho1 + rho2 )
        !ming 01-23-2023 cal H2O avg vmr
        H2Oavg(l) = ( (rho1*vmr_usr(l,1)) + (rho2*vmr_usr(l+1,1)) ) / ( rho1 + rho2 )
        p1 = p(l)
        p2 = p(l+1)
        dp = p2 - p1
        Pavg(l) = 0.5*(p1+p2)   ! density weigting is used (LBLRTM)
        g_avg = 0.5*(g(l)+g(l+1))

!c----Calculate air amount in a layer (similar to LBLRTM approach):
        airAmt(l)=dp*avogad*10./g_avg*2./(mw_air(l)+mw_air(l+1))

        !Calculate LAYER amounts for atmospheric gases:
        !mingmod amax
        DO k=1,ngas_lut
           w1 = amax1(vmr_usr(l,k),1.0e-20)
           w2 = amax1(vmr_usr(l+1,k),1.0e-20)

!c---  Kcarta approach:
           w_avg =  ( (rho1*w1) + (rho2*w2) ) / ( rho1+rho2 )
           !Wgas1(l,k) = w_avg * airAmt(l)
!c---  LBLRTM approach:
           CALL integAmt(w1,w2,p1,p2,g(l),g(l+1),mw_air(l),mw_air(l+1),Wgas(l,k))
        ENDDO
      ENDDO


END SUBROUTINE Calc_Layer_Amount




SUBROUTINE integAmt(w1,w2,p1,p2,g1,g2,mwAir1,mwAir2,gasLayAmt)

   implicit none


      REAL w1,w2,p1,p2,g1,g2,mwAir1,mwAir2,gasLayAmt, zeta, y,alpha,x

      zeta=log(p1/p2)

      alpha=log(w1/w2*g2/g1*mwAir2/mwAir1)/(zeta+1.e-20)

      x=(1.0+alpha)*zeta

      IF(abs(x) > 0.0002) then
         y=(1.0-exp(x))/x
      ELSE
         y=-1.00
      ENDIF
      gasLayAmt=w2*avogad*10.0/g2/mwAir2*p2*zeta*y

END SUBROUTINE integAmt

SUBROUTINE gphite(p, t, w, z_sfc,n_levels, i_dir, z)
      implicit none
      REAL, dimension(1:n_levels),      intent(in)   :: p, t, w
      INTEGER,                          intent(in)   :: n_levels,i_dir
      REAL,                             intent(in)   :: z_sfc
      REAL, dimension(1:n_levels),      intent(out)  :: z
!-----------------------------------------------------------------------
!                         -- Local variables --
!-----------------------------------------------------------------------
! -- Parameters
! -- rog=R/g/mw_dryair=8.3143/9.80665/28.97=0.029265538(km/K)=29.2655(m/K)
!      REAL            ::    rog, fac
      REAL, parameter ::  rog = 29.2898, fac = 0.5 * rog
! -- Scalars
      INTEGER         :: i_start, i_END, l
      REAL            :: v_lower, v_upper, algp_lower, algp_upper, hgt


!***********************************************************************
!                         ** Executable code **
!***********************************************************************
!-----------------------------------------------------------------------
!  -- Calculate virtual temperature adjustment and exponential       --
!  -- pressure height for level above surface.  Also set integration --
!  -- loop bounds                                                    --
!-----------------------------------------------------------------------
      IF( i_dir > 0 ) then    ! Data stored top DOwn, i_dir=1
        v_lower = t(n_levels) * ( 1.0 + ( 0.00061 * w(n_levels) ) )
        algp_lower = log( p(n_levels) )
        i_start = n_levels-1
        i_END   = 1
      ELSE                      ! Data stored bottom up

        v_lower = t(1) * ( 1.0 + ( 0.00061 * w(1) ) )
        algp_lower = log( p(1) )
        i_start = 2
        i_END   = n_levels
      END IF
      hgt = z_sfc   ! Assign surface height

!-----------------------------------------------------------------------
!             -- Loop over layers always from sfc -> top --
!-----------------------------------------------------------------------
      DO l = i_start, i_END, -1*i_dir
!       ----------------------------------------------------
!       Apply virtual temperature adjustment for upper level
!       ----------------------------------------------------
        v_upper = t(l)
        IF( p(l) >= 300.0 )  v_upper = v_upper * ( 1.0 + ( 0.00061 * w(l) ) )
!       ----------------------------------------------------- 
!       Calculate exponential pressure height for upper layer
!       ----------------------------------------------------- 
        algp_upper = log( p(l) )
!       ----------------
!       Calculate height use Eq. d[ln(p)]=-[(mw_dryair*g_sfc/R_gas)/T]*d(hgt)
!       ----------------
        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )
!       -------------------------------
!       Overwrite values for next layer
!       -------------------------------
        v_lower = v_upper
        algp_lower = algp_upper
!       ---------------------------------------------
!       Store heights in same direction as other data
!       ---------------------------------------------
        z(l) = hgt
      END DO
      z(n_levels) = z_sfc
!      return
END SUBROUTINE gphite


SUBROUTINE gas_OD_LayerAmt (path0, nlayers,ngases,wnStart,wnEnd,ndat,ndat_band, wavenums,waveindex,gasIdBand, pres_usr,& !input
                       nlevel_atm,psfc_out,pres_atm,temp_atm,vmr_atm,&                           !input
                       odIO_tot,airAmt_usr,Wgas_usr)                                             !output

   implicit none

   CHARACTER(LEN=80), intent(in) :: path0
   INTEGER,           intent(in) :: ngases,nlayers, ndat,ndat_band
   REAL*8,            intent(in) :: wavenums(1:maxlambdas)
   INTEGER,           intent(in) :: waveindex(1:maxlambdas)
   INTEGER,           intent(in) :: gasIdBand(1:maxgases)
   REAL,              intent(in) :: pres_usr(1:maxlayers)
   REAL*8,            intent(in) :: wnStart, wnEnd
   INTEGER,           intent(in) :: nlevel_atm
   !ming 01-03-2023  add the sfc pressure and temp
   REAL,              intent(in) :: psfc_out
   REAL,              intent(in) :: temp_atm  ( 1:maxAtmLev)    !MAXLAYERS )
   REAL,              intent(in) :: pres_atm  ( 1:maxAtmLev)    !MAXLAYERS )
   REAL,              intent(inout) :: vmr_atm(1:maxAtmLev, 1:ngas_lut)





   !local
   INTEGER :: gasID!,nlevel_usr
   INTEGER               :: k,k1,k2, n, ii, jj,l,ilev
   REAL    :: pres_usr_local(1:maxlayers)
   REAL    :: temp_usr(1:maxlayers)
   REAL    :: vmr_usr(1:MAXLAYERS, 1:ngas_lut)
   REAL    :: gasOD(1:nlay_lut,1:ndat)
   REAL    :: gasOD_all(1:ngases,1:nlayers, 1:ndat)
   
   !ming0801 add
   REAL    :: gasOD_usr(1:nlayers,1:ndat)
   REAL, dimension(1:ndat) :: od_low,od_up,od_mid
   REAL    :: airAmt_low,airAmt_up,airAmt_mid
   REAL    :: wgas_low,wgas_up,wgas_mid
   
   
   !ming0801 local variable for find grid 
   REAL    :: pres_lutgrid(1:nlev_lut)              !MAXLAYERS)
   INTEGER :: nlevel_lutgrid,nlayers_lutgrid
   REAL    :: vmr_lutgrid(1:nlev_lut, 1:ngas_lut)
   REAL    :: temp_lutgrid(1:nlev_lut)
   REAL, allocatable :: Wgas_lutgrid (:,:)
   REAL, allocatable :: gph_lutgrid (:)
   REAL, allocatable :: airAmt_lutgrid (:)
   !REAL, allocatable :: odIO_tot_lutgrid(:,:)
   REAL, allocatable :: Tavg_lutgrid(:),H2Oavg_lutgrid(:)
   
   REAL  :: temp_at_psfc


   REAL, intent(out)    :: Wgas_usr(1:nlayers, 1:ngas_lut)
   REAL, intent(out)    :: airAmt_usr(1:nlayers)
   REAL, intent(out)    :: odIO_tot(1:nlayers, 1:ndat)


  !ming 01-23-2023 CALL in SUBROUTINE (not in driver) 
  CALL read_gasLUT_P(path0,pres_lutgrid)         !pres_lutgrid is from TOA to sfc

  !assign local variable
  nlevel_lutgrid = nlev_lut        !nlev_lut =101 


  !interpolate atm profile to LUT's press grid
  !nlayers_lutgrid is based on actual pres surface, psfc_out
  CALL Interpolate_LUT_grid(pres_lutgrid,nlevel_lutgrid,nlevel_atm,&               !input
                             psfc_out,pres_atm,temp_atm,vmr_atm, &                  !input
                             temp_lutgrid,vmr_lutgrid,nlayers_lutgrid,temp_at_psfc) !output


  IF (allocated(Wgas_lutgrid) ) deallocate(Wgas_lutgrid)
  allocate ( Wgas_lutgrid(1:nlayers_lutgrid, 1:ngas_lut) )

  IF (allocated(gph_lutgrid) ) deallocate(gph_lutgrid)
  allocate ( gph_lutgrid(1:nlayers_lutgrid+1) )

  IF (allocated(airAmt_lutgrid) ) deallocate(airAmt_lutgrid)
  allocate ( airAmt_lutgrid(1:nlayers_lutgrid) )

 IF (allocated(Tavg_lutgrid) ) deallocate(Tavg_lutgrid)
  allocate ( Tavg_lutgrid(1:nlayers_lutgrid) )

 IF (allocated(H2Oavg_lutgrid) ) deallocate(H2Oavg_lutgrid)
  allocate ( H2Oavg_lutgrid(1:nlayers_lutgrid) )

!cal gas amount at each lut's grid
  CALL Calc_Layer_Amount(nlayers_lutgrid,psfc_out,temp_at_psfc,pres_lutgrid,temp_lutgrid,vmr_lutgrid,&    !input
                         Wgas_lutgrid,gph_lutgrid,airAmt_lutgrid,Tavg_lutgrid,H2Oavg_lutgrid)      !output
 
  !note: Tavg_lutgrid,H2Oavg_lutgrid is layer averaged values (air density weighted) 

!must initialize
 gasOD_all(:,:,:) = 0.0

 DO n =1, ngases
    gasID = gasIdBand(n)          !gas ID in hitran order
    print*,'gasid',gasid

!cal OD at lut's pres grid (100 layers)
    CALL Calc_Gas_OD (path0,nlayers_lutgrid,wnStart,wnEnd,ndat,ndat_band,wavenums,waveindex,gasID,&              !input
                      nlevel_lutgrid,pres_lutgrid,temp_lutgrid,vmr_lutgrid,Wgas_lutgrid,& !input
                      Tavg_lutgrid,H2Oavg_lutgrid,&                                       !input
                      gasOD)                                                              !output


!assign local variable to avoid changing pres_usr 
!add patial layer IF the psfc between the LUT's press grid
    pres_usr_local = pres_usr
    pres_usr_local(nlayers+1) = psfc_out


!assign fine grid (100 layer) OD to coase user grid
    DO ilev = nlayers+1,2,-1
      DO l = nlayers_lutgrid+1,1,-1
          IF (pres_lutgrid(l) .lt. pres_usr_local(ilev) ) exit
      ENDDO
!revise above .le. to .lt.; below k1=l to k1=l+1
      k1=l+1

      DO l = nlayers_lutgrid+1,1,-1
          IF (pres_lutgrid(l) .le. pres_usr_local(ilev-1) ) exit
      ENDDO
      k2 = l


      IF (k1-1 .eq. k2) then
           gasOD_Usr(ilev-1,:) = gasOD(k1-1,:) * &
                                 log( pres_usr_local(ilev) /pres_usr_local(ilev-1) ) / log( pres_lutgrid(k1) / pres_lutgrid(k2)  )

           airAmt_usr(ilev-1) = airAmt_lutgrid(k1-1) * &
                                 log( pres_usr_local(ilev) /pres_usr_local(ilev-1) ) / log(pres_lutgrid(k1) / pres_lutgrid(k2)  )

           Wgas_usr(ilev-1,gasID) = Wgas_lutgrid(k1-1,gasID) * &
                                 log( pres_usr_local(ilev) /pres_usr_local(ilev-1) ) /log(pres_lutgrid(k1) / pres_lutgrid(k2)  )
      ELSE
           od_mid = 0.0
           airAmt_mid =0.0
           wgas_mid = 0.0 
           DO k = k1-2,k2+1,-1
             od_mid(:) = od_mid(:) + gasOD(k,:)  
             airAmt_mid = airAmt_mid + airAmt_lutgrid(k)  
             wgas_mid = wgas_mid + Wgas_lutgrid(k,gasID)
           ENDDO
           od_low(:) = gasOD(k1-1,:) * &
                       log( pres_usr_local(ilev)/pres_lutgrid(k1-1) ) / log( pres_lutgrid(k1)/pres_lutgrid(k1-1)  )
           airAmt_low = airAmt_lutgrid(k1-1) * &
                       log( pres_usr_local(ilev)/pres_lutgrid(k1-1) ) / log(pres_lutgrid(k1)/pres_lutgrid(k1-1)  )
           wgas_low = wgas_lutgrid(k1-1,gasID) * &
                       log( pres_usr_local(ilev)/pres_lutgrid(k1-1) ) /log(pres_lutgrid(k1)/pres_lutgrid(k1-1)  )


           od_up(:) = gasOD(k2, :) * &
                      log( pres_lutgrid(k2+1)/pres_usr_local(ilev-1) ) / log( pres_lutgrid(k2+1)/pres_lutgrid(k2)  )
           airAmt_up = airAmt_lutgrid(k2) * &
                      log( pres_lutgrid(k2+1)/pres_usr_local(ilev-1) ) / log(pres_lutgrid(k2+1)/pres_lutgrid(k2)  )
           wgas_up = wgas_lutgrid(k2,gasID) * &
                      log( pres_lutgrid(k2+1)/pres_usr_local(ilev-1) ) /log(pres_lutgrid(k2+1)/pres_lutgrid(k2)  )
            gasOD_Usr(ilev-1,:) = od_mid(:) + od_low(:) + od_up(:)
            airAmt_Usr(ilev-1) = airAmt_mid + airAmt_low + airAmt_up
            wgas_Usr(ilev-1,gasID) = wgas_mid + wgas_low + wgas_up        
      ENDIF

    ENDDO


    DO ii = 1, nlayers!nlevel_usr-1
       DO jj = 1, ndat
          gasOD_all(n,ii,jj) = gasOD_usr (ii,jj) !
       ENDDO
    ENDDO
 ENDDO

  !sum all gas OD, assign to odIO_tot for vliodrt input
  ! rev: ming this step take 12% of total gasOD time
  DO ii = 1, nlayers!nlevel_usr-1
    DO jj = 1, ndat
       odIO_tot (ii,jj) = sum( gasOD_all(1:ngases,ii,jj) )
    ENDDO
  ENDDO


END SUBROUTINE gas_OD_LayerAmt




SUBROUTINE Calc_Gas_OD (path0,nlayers, wnStart,wnEnd,ndat,ndat_band,wavenums,waveindex,gasID,&       !input
                       nlevel_usr,pres_usr,temp_usr,vmr_usr,Wgas_usr,&        !input
                       Tavg_usr, H2Oavg_usr,&                                 !input
                       odIO)                                                  !output

   implicit none

   !input
   CHARACTER(LEN=80), intent(in) :: path0
   INTEGER, intent(in)  :: nlayers,nlevel_usr
   REAL*8,  intent(in)  :: wnStart, wnEnd
   INTEGER, intent(in)  :: gasID,ndat,ndat_band
   REAL, dimension(nlev_lut), intent(in) :: temp_usr,pres_usr
   REAL,                       intent(in) :: vmr_usr(1:nlev_lut, 1:ngas_lut)
   REAL,                       intent(in) :: Wgas_usr(1:nlayers, 1:ngas_lut)
   REAL, dimension(nlayers),   intent(in) :: Tavg_usr,H2Oavg_usr
   INTEGER,                    intent(in) :: waveindex(1:maxlambdas) 
   REAL*8,                     intent(in) :: wavenums(1:maxlambdas)
   !local
   CHARACTER(len=10) :: wnStartStr, wnEndStr
   REAL*8  :: wresol
   INTEGER :: l, i, sigGas, npp,IFile
   REAL    :: fct, xint, yint, yyint
   REAL    :: dT(1:29)
   REAL,   dimension(1:3)         :: dH2o,xtab,ytab,yytab
   REAL,dimension(1:nlevel_usr-1) :: dTId, dh2oId
   CHARACTER(len=6)               :: gas(1:ngas_lut)
   INTEGER :: find,pp,pp_lut
   INTEGER :: tind1, tind2, tind
   REAL,dimension(1:ndat_band) :: odUsr_1lay,odUsr_1lay1
   REAL,dimension(1:ndat_band) :: odUsr_1lay2
   REAL,dimension(1:ndat_band) :: odUsr_1lay11,odUsr_1lay21
   REAL,dimension(1:ndat_band) :: odUsr_1lay12,odUsr_1lay22
   REAL,dimension(1:ndat_band) :: odUsr_1lay31,odUsr_1lay32
   REAL*4, dimension(1:nlay_lut,1:ndat_band) :: lut2D11, lut2D12, lut2D21,lut2D22
   REAL*4, dimension(1:nlay_lut,1:ndat_band) :: lut2D31, lut2D32
   CHARACTER(len=2)                :: tid1,tid2, tid
   CHARACTER(len=2)                :: wid1,wid2,wid3
   CHARACTER(len=100)              :: fname11,fname21,fname22,fname31
   REAL*4, dimension(1:ndat_band)  :: lut2D_1layer
   REAL,dimension(1:nlevel_usr-1)  :: Pavg_lut, Tavg_lut, H2Oavg_lut  !layer avg for lut 'std' 
   CHARACTER(LEN=80)               :: filename
   INTEGER :: idx 
   REAL*8, dimension(1:3) :: x_wvtab, y_odtab
   REAL*8 :: y_odint 
   !output
   REAL, intent(out)  :: odIO(1:nlay_lut,1:ndat)

!!  ## pseuDO-code: 
!!  ##   
!!  ##   read in temp, vmr_h2o, find the index of temp and H2O 
!!  ##   
!!  ##   IF gas .not. H2O then
!!  ##       DO_pres_interp (interpolate XS, then *column density, get OD) at temp index1 and index2
!!  ##       temp. linear interpolation
!!  ##   ELSE !gas =H2O
!!  ##       DO_pres_interp at water index 1 (0.1 ratio of standard), index2 (ratio=1),index3 (ration = 10) at temp index1
!!  ##       DO quad interpolation for index1,2,3, get the interpolated OD at temp index1
!!  ##   
!!  ##       DO pres_interp at water index 1 (0.1 ratio of standard), index2 (ratio=1), index3 (ration = 10) at temp index2
!!  ##       DO quad interpolation for index1,2,3, get the interpolated OD at temp index2
!!  ##       temp. liner interpolation
!!  ##   ENDIF
!

     write(wnStartStr, '(i5.5)') int(wnStart)
     write(wnEndStr, '(i5.5)') int(wnEnd)
     wresol = (wnEnd-wnStart)/(ndat_band-1)

! INTERPOLATING US STANDARD TEMP TO USR DEFINED PRESS GRID

   !ming 01-23-2023 read in layer averge T VMR_h2o for 'std' lut profile
   filename='PCRTM_VLIDORT_Config/lut_PTH2O_layerAvg.dat'     !user_pres_grid.dat'
   CALL GETLUN(IFILE)
   open ( unit = IFile, file = trim(path0)//filename,action='read',status='old')
   rewind ( IFile )
   DO i = 1, nlevel_usr-1                      !nlevel_usr =100
      read(IFile,*) Pavg_lut(i), Tavg_lut(i), H2Oavg_lut(i)
   END DO
   close ( IFile )

   !ming 01-23-2023 using the Tavg_lut and Tavg_usr to cal dTid
   !                using the h2oavg_lut and h2oavg_usr to cal dh2oid
   DO i = 1, nlayers
      dTId(i) = Tavg_usr(i) - Tavg_lut(i)
      dh2oId(i) = H2Oavg_usr(i)/H2Oavg_lut(i)
   ENDDO


! SET TEMP DIFF FACTOR
   DO i = 1, 29
     dT(i) = (i-15) * 5             !temp dt range:-70~70 K  
   ENDDO

! SET H2O VMR FACTOR
  dH2o = (/0.1,1.0,10.0/)


!! SETUP 45 GAS NAMES, USED TO IDENTIFY LUT'S FILENAME
!  1~37 is same Hitran order gases
!  42~45 is for ISOP....
! G20(H2CO), G21(HOCL),G24(CH3CL),G25(H2O2),G26(C2H2)
   gas = (/'H2O' ,'CO2', 'O33', 'N2O', 'COO', 'CH4' ,'OO2', 'NO1', 'SO2', 'NO2',&
           'NH3', 'G12', 'AOH', 'AHF', 'HCL', 'HBR', 'AHI', 'CLO', 'OCS', 'H2CO',    &
           'G21', 'AN2', 'HCN', 'G24', 'G25', 'G26', 'C2H6','PH3', 'COF2','SF6',&
           'H2S', 'H2CO2','HO2','O  ','CLONO2','NOPLUS','HOBR','BLK','BLK','BLK',  &
           'BLK', 'ISOP','CFC113','OO4 ','BRO   '/)

print*,'gas name:',gas(gasID)


!!!*************************************************************************************************************
!! M A I N   I N T E P O L A T I O N   S T A R T   H E R E
!!!*************************************************************************************************************

! START LUT INTERPOLATION AND LOOPING EACH USR DEFINED PRESS LAYER        
! Pres loop from TOA to Surface, nlayers is # of layers for LUT press grid
       DO pp = 1, nlayers   
         ! PREPARE TEMPERATURE  INTERPOLATION (FIND tid1 & tid2):
         ! DECIDE IF NEED TO DO TEMPERATURE INTERPOLATION
            find = 0
 find_loop: DO i = 1, 29
               IF ( abs(dT(i) - dTId(pp)) < 1.e-5 ) then
               !IF ( abs(dT(i) - dTId(pp)) < 1.e-8 ) then
                find = 1
                tind = i
                exit find_loop
               END IF
            END DO find_loop


         IF (find == 1) then                      !when the dTid is same as the dT
             write(tid,'(I2.2)') tind+1           !tid for filename
             ! SETTING TEMP ID FOR LUT FILE NAME
             write(tid1,'(I2.2)') tind             !tid for filename
             write(tid2,'(I2.2)') tind

         ! START SETTING TEMPERATURE INTERPOLATION
         ELSE  ! IF (find == 1) then 
             tind1 = -1
             tind2 = -1
   find_interval: DO i = 1, 29
                    IF ( dT(i) > dTId(pp) )then
                      tind2 = i
                      exit find_interval
                    END IF
                  END DO find_interval
             tind1 = tind2 - 1 ! index of dT
             fct =  (dTId(pp)-dT(tind1))/(dT(tind2)-dT(tind1))
             IF (tind2 == -1 .or. tind1 == -1) then
                print*,'Warning: temp was out of LUT range at layer:',pp
                tind1=29
                tind2=29
                print*,'tind1=tind2=29'!'STOPPING'
                fct = 1.0
                !stop
             ELSE IF (tind1 == 0) then
                print*,'Warning: temp was out of LUT range at layer:',pp
                tind1=1
                tind2=1
                print*,'tind1=tind2=1'
                 fct = 1.0
             END IF

             ! SETTING TEMP ID FOR LUT FILE NAME
             write(tid1,'(I2.2)') tind1 !tid for filename
             write(tid2,'(I2.2)') tind2
         ENDIF !IF (find == 1) then


          ! ming 01-23-2022 find the layer index for LUT (LUT 1 is suface)
          pp_lut = nlay_lut- pp +1 

          ! DO TEMP ONLY INTERP FOR GASES EXCEPT H2O
          IF ( gasID /= 1 ) then                              !gasID=0 --> H2O   
                fname21=trim(path0)//'data_and_control/GAS_XS/XS_'//trim(gas(gasID))//   &
                         '.'//trim(wnStartStr)//'.'//trim(wnEndStr)//'.W02_direct'
                CALL GETLUN(IFILE)
                open(unit=IFile,file=trim(fname21),form='unformatted',access='direct', recl=ndat_band*4)

                  read(IFile, rec=pp_lut+(tind1-1)*100), lut2D_1layer
                  lut2D21(pp,:) = lut2D_1layer(:)   !ming 09/20/2022 assign value cost long time?
                close (IFile)


                fname22=trim(path0)//'data_and_control/GAS_XS/XS_'//trim(gas(gasID))//&
                         '.'//trim(wnStartStr)//'.'//trim(wnEndStr)//'.W02_direct'
                CALL GETLUN(IFILE)
                open(unit=IFile,file=trim(fname22),form='unformatted',access='direct',recl=ndat_band*4)

                  read(IFile, rec=pp_lut+(tind2-1)*100), lut2D_1layer
                  lut2D22(pp,:) = lut2D_1layer(:)
                close (IFile)

               odUsr_1lay21(:) = lut2D21(pp,:) * Wgas_usr(pp, gasID) !rhoUsrLut(ilev,gasID) 
               odUsr_1lay22(:) = lut2D22(pp,:) * Wgas_usr(pp, gasID)
               odUsr_1lay = odUsr_1lay21 + ( odUsr_1lay22 - odUsr_1lay21  ) *fct


          ! END 'regular gas' temp INTERP, START FOR H2O
          ELSE

                ! READ SIX  H2O  TAPE5 FILE (WITH 3 WID AND 2 TID) 
                sigGas = 1 ! 1 for TAPE5_H2O; 0 for TAPE5_ALL
                wid1 = '01'
                wid2 = '02'
                wid3 = '03'!'03' for wfac=10; '06' for wfact=4

                !ming 09/22/2022 read in the direct access gas LUT (29 temp wrap) for layer k1 to k2-1
                !                one W01 has both tind1 and tind2
                fname11=trim(path0)//'data_and_control/GAS_XS/XS_'//'H2O'//   &
                         '.'//trim(wnStartStr)//'.'//trim(wnEndStr)//'.W01_direct'
                CALL GETLUN(IFILE)
                open(unit=IFile,file=trim(fname11),form='unformatted',access='direct', recl=ndat_band*4)

                  read(IFile, rec=pp_lut+(tind1-1)*100), lut2D_1layer
                  lut2D11(pp,:) = lut2D_1layer(:)   
                  read(IFile, rec=pp_lut+(tind2-1)*100), lut2D_1layer
                  lut2D12(pp,:) = lut2D_1layer(:)
                close (IFile)

                fname21=trim(path0)//'data_and_control/GAS_XS/XS_'//'H2O'//   &
                         '.'//trim(wnStartStr)//'.'//trim(wnEndStr)//'.W02_direct'
                CALL GETLUN(IFILE)
                open(unit=IFile,file=trim(fname21),form='unformatted',access='direct', recl=ndat_band*4)

                  read(IFile, rec=pp_lut+(tind1-1)*100), lut2D_1layer
                  lut2D21(pp,:) = lut2D_1layer(:)  
                  read(IFile, rec=pp_lut+(tind2-1)*100), lut2D_1layer
                  lut2D22(pp,:) = lut2D_1layer(:)
                close (IFile)


                fname31=trim(path0)//'data_and_control/GAS_XS/XS_'//'H2O'//   &
                         '.'//trim(wnStartStr)//'.'//trim(wnEndStr)//'.W03_direct'
                CALL GETLUN(IFILE) 
                open(unit=IFile,file=trim(fname31),form='unformatted',access='direct', recl=ndat_band*4)

                  read(IFile, rec=pp_lut+(tind1-1)*100), lut2D_1layer
                  lut2D31(pp,:) = lut2D_1layer(:)  
                 read(IFile, rec=pp_lut+(tind2-1)*100), lut2D_1layer
                 lut2D32(pp,:) = lut2D_1layer(:)
                close (IFile)


                !at tind2
                odUsr_1lay11(:) = lut2D11(pp,:) * Wgas_usr(pp, gasID)
                odUsr_1lay21(:) = lut2D21(pp,:) * Wgas_usr(pp, gasID)
                odUsr_1lay31(:) = lut2D31(pp,:) * Wgas_usr(pp, gasID)
                !at tind2
                odUsr_1lay12(:) = lut2D12(pp,:) * Wgas_usr(pp, gasID)
                odUsr_1lay22(:) = lut2D22(pp,:) * Wgas_usr(pp, gasID)
                odUsr_1lay32(:) = lut2D32(pp,:) * Wgas_usr(pp, gasID)

                ! 3-POINT LAGRANGE INTERPOLATION FOR H20
                xtab = (/0.1, 1.0, 10.0/)
                xint =  dh2oId(pp)
                DO i = 1, ndat_band
                  ytab(1) = odUSR_1lay11(i)
                  ytab(2) = odUSR_1lay21(i)
                  ytab(3) = odUSR_1lay31(i)
                  CALL quadterp(xtab,ytab,xint,yint)
                  odUsr_1lay1(i) = yint

                  yytab(1) = odUSR_1lay12(i)
                  yytab(2) = odUSR_1lay22(i)
                  yytab(3) = odUSR_1lay32(i)

                  CALL quadterp(xtab,yytab,xint,yyint)
                  odUsr_1lay2(i) = yyint
                END DO

                ! LINEAR TEMPERATURE INTERPOLATION:
                odUsr_1lay = odUsr_1lay1 + ( odUsr_1lay2 - odUsr_1lay1  ) * fct

          END IF ! IF ( gasID /= 0 )
          ! *********** FINISHED INTERPOLATION *******************


          ! STORE EACH LAYER'S OD INTO ODIO ARRAY
          IF (ndat .eq. ndat_band) then
              IF ( gasID == 44) then
                DO npp=1,ndat!0,pts-1
                   odUsr_1lay(npp)=odUsr_1lay(npp)/1.D46!O4 XS need to divied by1E46
                   odIO(pp,npp)  = odUsr_1lay(npp)
                ENDDO
              ELSE
                DO npp=1,ndat_band!0,pts-1
                   odIO(pp,npp)  = odUsr_1lay(npp)!/1.D46
                ENDDO
              ENDIF
          ELSE
              IF ( gasID == 44) then
                DO npp=1,ndat
                   idx = waveindex(npp)
                   IF (idx > 1) then
                      x_wvtab(1) = wnStart+wresol*(idx-2);x_wvtab(2) = wnStart+wresol*(idx-1);x_wvtab(3) = wnStart+wresol*(idx)
                      y_odtab(1) = odUsr_1lay(idx-1); y_odtab(2) = odUsr_1lay(idx);y_odtab(3) = odUsr_1lay(idx+1)
                   ELSE
                      x_wvtab(1) = wnStart+wresol*(idx-1);x_wvtab(2) = wnStart+wresol*(idx);x_wvtab(3) = wnStart+wresol*(idx+1)
                      y_odtab(1) = odUsr_1lay(idx); y_odtab(2) = odUsr_1lay(idx+1);y_odtab(3) = odUsr_1lay(idx+2)
                   ENDIF
                   CALL quadterp2(x_wvtab,y_odtab,wavenums(npp),y_odint)
                   odIO(pp,npp)  = y_odint/1.D46 
                ENDDO
              ELSE
                DO npp=1,ndat
                   idx = waveindex(npp)
                  
                   IF (idx > 1) then
                      x_wvtab(1) = wnStart+wresol*(idx-2);x_wvtab(2) = wnStart+wresol*(idx-1);x_wvtab(3) = wnStart+wresol*(idx)
                      y_odtab(1) = odUsr_1lay(idx-1); y_odtab(2) = odUsr_1lay(idx);y_odtab(3) = odUsr_1lay(idx+1)
                   ELSE
                      x_wvtab(1) = wnStart+wresol*(idx-1);x_wvtab(2) = wnStart+wresol*(idx);x_wvtab(3) = wnStart+wresol*(idx+1)
                      y_odtab(1) = odUsr_1lay(idx); y_odtab(2) = odUsr_1lay(idx+1);y_odtab(3) = odUsr_1lay(idx+2)
                   ENDIF
                   CALL quadterp2(x_wvtab,y_odtab,wavenums(npp),y_odint)
                   odIO(pp,npp)  = y_odint 
                ENDDO

              ENDIF
          ENDIF

       ENDDO !DO pp=0, nlayer
       !*********************   FINISHING LAYER LOOP 
END SUBROUTINE  Calc_Gas_OD





END MODULE GAS_OPT_m
