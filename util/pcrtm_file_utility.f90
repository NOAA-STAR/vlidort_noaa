MODULE PCRTM_File_Utility

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GETLUN
  PUBLIC :: Write_Inputs
  PUBLIC :: Write_IQU, Write_3components, write_stokes

  CONTAINS

  SUBROUTINE GETLUN(LUN)

    INTEGER, INTENT(OUT) :: LUN
    INTEGER              :: I
    LOGICAL              :: OPENUNIT

    DO I = 10,99
      INQUIRE(I,OPENED = OPENUNIT)
      IF ( .NOT. OPENUNIT) EXIT
    ENDDO

    LUN = I

  END SUBROUTINE GETLUN


    SUBROUTINE Write_IQU (date,Nstokes, ncount, nstep,n_geometries, nProf_s,nameProf,Lambertian_Albedo,  &
                        Out_szas, Out_vzas, Out_azms,maxtypeaer,                             &
                        wnStart,wnEnd,wresol,nstreams,nprof,Psfc_vld,VLIDORT_ModIn,st_wvnum_indx, END_wvnum_indx,&
			Iup_output2,Qup_output2,Uup_output2)


      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

    IMPLICIT NONE

!  VLIDORT output structure

    TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
    TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

    CHARACTER(len=10)                   :: wnStartStr, wnEndStr
    CHARACTER(len=4)                   :: nProf_s
    CHARACTER(len=300)                   ::nameProf
    CHARACTER(LEN=100)                   :: output_filename
    CHARACTER(8)  :: date
    CHARACTER(len=:),allocatable  :: sub_name

    REAL,dimension(1:n_geometries,1:ncount) :: Iup_output2,Qup_output2,Uup_output2  
    REAL :: Psfc_vld
    INTEGER :: Indx_lambertian, indx_brdf, indx_bpdf,k,v,maxtypeaer, nstep,st_wvnum_indx, END_wvnum_indx

     REAL(8)                          :: wnStart, wnEnd, wresol
     REAL(8),allocatable :: Out_Tmp(:)
     REAL,dimension(1:n_geometries) ::   Out_szas, Out_vzas, Out_azms

     INTEGER      :: Nrec_len, npt, ncount,nstokes, ifile, n_geometries,nstreams,nprof,nmode_aer_in

    REAL(8),dimension(maxtypeaer) :: AOD_aer_in
    REAL(8),dimension(maxtypeaer) :: Aerosol_upperboundary_in,Aerosol_lowerboundary_in

!    INTEGER,dimension(:) :: type_aer_in,Aerosol_prof_in
!    REAL(8),dimension(:)  ::expAerosol_relaxation_in,gdfAerosol_peakheight_in,gdfAerosol_halfwidth_in

     REAL(8) :: Lambertian_Albedo
     write(wnStartStr, '(i5.5)') int(wnStart)
     write(wnEndStr, '(i5.5)') int(wnEnd)


      print*,ncount,nstokes

      npt = ncount
      Nrec_len = npt*nstokes+100

      print*,npt,Nrec_len

      if (allocated( out_tmp ) ) deallocate( out_tmp  )
      allocate ( out_tmp(Nrec_len))

      k = len_trim(nameProf)
      sub_name = trim(nameProf)

     if(Nstokes .eq. 1) then
       output_filename = 'Output_IQU/I.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.' //  &
           sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 2) then
       output_filename = 'Output_IQU/IQ.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//  &
           sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 3) then
       output_filename = 'Output_IQU/IQU.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'// &
           sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else
       output_filename = 'Output_IQU/IQUV.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'// &
           sub_name(1:k-4)//'_'//trim(date)//'.dat'
     ENDif

     print*,output_filename

      CALL GETLUN(IFILE)

      open(UNIT=IFILE, file= output_filename,form='unformatted',&
               access='direct', recl=Nrec_len*8)

      K = 0    
      DO v = 1, n_geometries

        Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,&
                 nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
                 DBLE(out_szas(v)), DBLE(out_vzas(v)), DBLE(out_azms(v)),DBLE(st_wvnum_indx), DBLE(END_wvnum_indx)/)
 
        Out_tmp(21:21+12)=(/dble(nmode_aer_in),AOD_aer_in,Aerosol_lowerboundary_in,Aerosol_upperboundary_in/)

        Out_tmp(51:54)=(/Indx_lambertian*1.D0, indx_brdf*1.D0, indx_bpdf*1.D0, Lambertian_Albedo/)

        !  lambertian (index) (-999), brdf_index, bpdf_index, aerosol_loading (bot, top),aersol_vertial_proftype (mixing?),aerosol_opt,

        Out_Tmp(101:npt+100) = Iup_output2(v,1:ncount)

        IF( Nstokes .eq. 3) THEN

          Out_Tmp(npt+101:npt*2+100) = Qup_output2(v,1:ncount)
          Out_Tmp(npt*2+101:npt*3+100) = Uup_output2(v,1:ncount)

        ENDIF
        k = k +1
        write(IFILE, rec=k) Out_tmp(1:Nrec_len)

        Out_Tmp(1:Nrec_len) = 0.0D0

      ENDDO

      Close (IFILE)

      print*, " writint output to file", output_filename
    END SUBROUTINE Write_IQU


    SUBROUTINE Write_3components (date,Nstokes, ncount, nstep,n_geometries, nProf_s,nameProf, Lambertian_Albedo, &
                        Out_szas, Out_vzas, Out_azms,maxtypeaer,                            &
                        wnStart,wnEnd,wresol,nstreams,nprof,Psfc_vld,VLIDORT_ModIn,st_wvnum_indx, END_wvnum_indx,wavenums,&
                        BTRANS_I, BTRANS_Q,BTRANS_U, BSPHER, &
			Iup_output2,Qup_output2,Uup_output2)


      USE VLIDORT_PARS_m
      USE VLIDORT_IO_DEFS_m

    IMPLICIT NONE

!  VLIDORT output structure

    TYPE(VLIDORT_Fixed_Inputs)             :: VLIDORT_FixIn
    TYPE(VLIDORT_Modified_Inputs)          :: VLIDORT_ModIn

    CHARACTER(len=10)                   :: wnStartStr, wnEndStr
    CHARACTER(len=4)                   :: nProf_s
    CHARACTER(len=300)                   ::nameProf
    CHARACTER(LEN=80)                   :: output_filename
    CHARACTER(8)  :: date
    CHARACTER(len=:),allocatable  :: sub_name

    REAL,dimension(1:n_geometries,1:ncount) :: Iup_output2,Qup_output2,Uup_output2  
    REAL,dimension(1:n_geometries,1:ncount) :: BTRANS_I, BTRANS_Q, BTRANS_U
    REAL,dimension(1:ncount) :: BSPHER

    REAL :: Psfc_vld
    INTEGER :: k,v,maxtypeaer, nstep

     REAL(8)                          :: wnStart, wnEnd, wresol
     REAL(8)                          ::  wavenums(1:ncount)                  !12-01-2023
     REAL(8),allocatable :: Out_Tmp(:)
     REAL,dimension(1:n_geometries) ::   Out_szas, Out_vzas, Out_azms

     INTEGER      :: Nrec_len, npt, ncount,nstokes, ifile, n_geometries,nstreams,nprof,nmode_aer_in,st_wvnum_indx, END_wvnum_indx

     REAL(8),dimension(maxtypeaer) :: AOD_aer_in
     REAL(8),dimension(maxtypeaer) :: Aerosol_upperboundary_in,Aerosol_lowerboundary_in

     REAL(8) ::   Lambertian_Albedo

     write(wnStartStr, '(i5.5)') int(wnStart)
     write(wnEndStr, '(i5.5)') int(wnEnd)


      print*,ncount,nstokes

       npt = ncount
       Nrec_len = npt*nstokes+100

       print*,npt,Nrec_len

       if (allocated( out_tmp ) ) deallocate( out_tmp  )
       allocate ( out_tmp(Nrec_len))

      k = len_trim(nameProf)
      sub_name = trim(nameProf)

     if(Nstokes .eq. 1) then
       output_filename = 'Output_IQU/Planetary_I.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//   &
       sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 2) then
       output_filename = 'Output_IQU/Planetary_IQ.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//  &
       sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 3) then
       output_filename = 'Output_IQU/Planetary_IQU.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'// &
       sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else
       output_filename = 'Output_IQU/Planetary_IQUV.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'// &
       sub_name(1:k-4)//'_'//trim(date)//'.dat'
     ENDif

       print*,output_filename

        CALL GETLUN(IFILE)

        open(UNIT=IFILE, file= output_filename,form='unformatted',&
               access='direct', recl=Nrec_len*8)

      Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,&
                  nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
                  DBLE(out_szas(1)), DBLE(out_vzas(1)), DBLE(out_azms(1)),DBLE(st_wvnum_indx), DBLE(END_wvnum_indx)/)

      Out_Tmp(101:npt+100) = wavenums(1:ncount)
       write(IFILE, rec=1) Out_tmp(1:Nrec_len)

      k = 1 
! added 12-01-2023, add output of wavenumbers
      DO v = 1, n_geometries

        Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,&
               nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
               DBLE(out_szas(v)), DBLE(out_vzas(v)), DBLE(out_azms(v)),DBLE(st_wvnum_indx), DBLE(END_wvnum_indx)/)
 
        Out_tmp(21:21+12)=(/dble(nmode_aer_in),AOD_aer_in,Aerosol_lowerboundary_in,Aerosol_upperboundary_in/)

        Out_tmp(54)=Lambertian_Albedo

! --  I1 ---
        Out_Tmp(101:npt+100) = Iup_output2(v,1:ncount)

      IF( Nstokes .eq. 3) THEN

        Out_Tmp(npt+101:npt*2+100) = Qup_output2(v,1:ncount)
        Out_Tmp(npt*2+101:npt*3+100) = Uup_output2(v,1:ncount)

      ENDIF
         k = k +1
         write(IFILE, rec=k) Out_tmp(1:Nrec_len)

! -- Transmittance ---

        Out_Tmp(1:Nrec_len) = 0.0D0

        Out_Tmp(101:npt+100) = BTRANS_I(v,1:ncount)

      IF( Nstokes .eq. 3) THEN

        Out_Tmp(npt+101:npt*2+100) = BTRANS_Q(v,1:ncount)
        Out_Tmp(npt*2+101:npt*3+100) = BTRANS_U(v,1:ncount)

      ENDIF

         k = k +1
         write(IFILE, rec=k) Out_tmp(1:Nrec_len)

        ENDDO

! -- Spherical albedo ---
        Out_Tmp(1:Nrec_len) = 0.0D0

        Out_Tmp(101:npt+100) = BSPHER(1:ncount)

         k = k +1
         write(IFILE, rec=k) Out_tmp(1:Nrec_len)

	Close (IFILE)

       print*, " writint output to file", output_filename
    END SUBROUTINE Write_3components


!**************************************************************************
    SUBROUTINE Write_Inputs (Nstokes, n_geometries,nProf_s,nameProf,LAMBERTIAN_ALBEDO,  &
                        Out_szas, Out_vzas, Out_azms,&
                        nstreams,nprof,Psfc_vld,htgrid_vld,&
                        maxtypeaer,nmode_aer_in,AOD_aer_in,Aerosol_lowerboundary_in,Aerosol_upperboundary_in,&
                        Aerosol_prof_in,expAerosol_relaxation_in,gdfAerosol_peakheight_in,gdfAerosol_halfwidth_in, &
                        Loading_Aerosol_Aod,  type_aer_in, MAXLAYERS,&
!                        Indx_lambertian, indx_brdf, indx_bpdf, &
                        Indx_lambertian, indx_brdf, indx_bpdf,  WINDSPEED, Salinity,Nkernels,BRDF_kernels,&
                        ATM_type,nlayers,ngas_lut,pres_usr,temp_usr,vmr_usr)


    IMPLICIT NONE

!  VLIDORT output structure

    Integer :: atm_type, nLayers, ngas_lut, nprof,maxtypeaer, maxlayers
    REAL    :: pres_usr  ( 1:nlayers+1 )
    REAL  ::  temp_usr  ( 1:nlayers+1 )
    REAL*8 :: htgrid_vld(0:nlayers)                   ! 11/15/2022, changed from REAL to REAL*8, Xiong
    REAL :: vmr_usr(1:nlayers+1, 1:ngas_lut)
    REAL :: Psfc_vld

    CHARACTER(len=4)                   :: nProf_s
    CHARACTER(len=300)                   ::nameProf
    CHARACTER(LEN=80)                   :: output_filename

     INTEGER :: Indx_lambertian, indx_brdf, indx_bpdf,k,v

     REAL(8),allocatable :: Out_Tmp(:)
     REAL,dimension(1:n_geometries) ::   Out_szas, Out_vzas, Out_azms

     INTEGER      :: Nrec_len, nstokes, ifile,n_geometries,nstreams,nmode_aer_in,i

     REAL(8) :: LAMBERTIAN_ALBEDO
     REAL(8),dimension(maxtypeaer) :: AOD_aer_in
     REAL(8),dimension(maxtypeaer) ::Aerosol_upperboundary_in,Aerosol_lowerboundary_in
     REAL(8),dimension(maxtypeaer)::expAerosol_relaxation_in,gdfAerosol_peakheight_in,gdfAerosol_halfwidth_in
     INTEGER,dimension(maxtypeaer) :: type_aer_in,Aerosol_prof_in

     REAL*8, dimension(maxtypeaer,MAXLAYERS)    :: Loading_Aerosol_Aod
     REAL*8 ::  WINDSPEED, Salinity     ! 8/16/2023
 
     INTEGER :: min_val, max_val, nsza,nvza,nazm, val(1:n_geometries), unique(1:n_geometries)    ! Xiong, 11/18/2022
     INTEGER ::  Nkernels,BRDF_kernels(1:Nkernels)

       Nrec_len = 100

       print*,Nrec_len

! ----------------------

       if (allocated( out_tmp ) ) deallocate( out_tmp  )
       allocate ( out_tmp(Nrec_len))

       output_filename ='Output_IQU/VLIDORT_Inputs.'//trim(nProf_s)//'.'//trim(nameProf)

       print*,output_filename

        CALL GETLUN(IFILE)

        open(UNIT=IFILE, file= output_filename,form='unformatted',&
               access='direct', recl=Nrec_len*8)
!updated on 5/10/2023
        Out_tmp(1:11)=(/DBLE(ngas_lut),DBLE(nlayers),DBLE(nstokes),DBLE(n_geometries),&
	       nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
               Indx_lambertian*1.D0, indx_brdf*1.D0, indx_bpdf*1.D0, Lambertian_Albedo/)      
        Out_tmp(12:12+12)=(/dble(nmode_aer_in),AOD_aer_in,Aerosol_lowerboundary_in,Aerosol_upperboundary_in/)
        Out_tmp(25:25+15)=&
            (/dble(Aerosol_prof_in),expAerosol_relaxation_in,gdfAerosol_peakheight_in,gdfAerosol_halfwidth_in/)
        Out_tmp(41:44)= DBLE(type_aer_in)
        Out_tmp(45:45+n_geometries*3-1)= &
	       (/DBLE(out_szas(1:n_geometries)),DBLE(out_vzas(1:n_geometries)),DBLE(out_azms(1:n_geometries))/)
!note if n_geometries gt 18, it will make error !
        write(IFILE, rec=1) Out_tmp

       Out_tmp(1:3+nkernels) = (/WINDSPEED, Salinity, dble(Nkernels),dble(BRDF_kernels(1:Nkernels))/)
!  using the 1st twenty for ocean, dust etc. 8/16/2023  
        Out_tmp(21:nlayers+21) = DBLE(pres_usr(1:nlayers+1))
        write(IFILE, rec=2) Out_tmp
       Out_tmp(1:nlayers+1) = DBLE(temp_usr(1:nlayers+1))
        write(IFILE, rec=3) Out_tmp
       Out_tmp(1:nlayers+1) = DBLE(htgrid_vld(0:nlayers))
        write(IFILE, rec=4) Out_tmp
        k=4
        DO i=1,ngas_lut
          k = k+1
          Out_tmp(1:nlayers+1) = DBLE(vmr_usr(1:nlayers+1, i))
          write(IFILE, rec=k) Out_tmp
        ENDDO
        DO i=1,maxtypeaer
          k = k+1
          Out_tmp(1:nlayers) = Loading_Aerosol_Aod(i,1:nlayers)
          write(IFILE, rec=k) Out_tmp
        ENDDO
        Close (IFILE)

      print*, " writing truth to file", output_filename
    END SUBROUTINE Write_Inputs



!******************************************************************************
      SUBROUTINE write_stokes &
        ( un, o1, a1, max_tasks, n_geometries, n_szangles,  n_user_levels, nt, &
          user_levels, STOKES_PT, MEANST_PT, FLUX_PT )

      USE vlidort_pars_m, Only : max_user_levels, max_geometries, max_szangles, &
                                 maxstokes, max_directions, upidx, dnidx

      implicit none

!  input variables

      INTEGER, intent(in) ::           max_tasks, n_geometries, n_szangles, n_user_levels
      INTEGER, intent(in) ::           un, nt, o1
      CHARACTER (len=1), intent(in) :: a1
      DOUBLE PRECISION, intent(in)  :: user_levels ( max_user_levels )

      DOUBLE PRECISION, INTENT(IN)  :: STOKES_PT &
          ( MAX_USER_LEVELS, MAX_GEOMETRIES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN)  :: MEANST_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES, MAX_DIRECTIONS, MAX_TASKS )
      DOUBLE PRECISION, INTENT(IN)  :: FLUX_PT &
          ( MAX_USER_LEVELS, MAX_SZANGLES, MAXSTOKES ,MAX_DIRECTIONS, MAX_TASKS )

!  local variables

      INTEGER :: n, v, t

!  Write stokes component

      write(un,'(/T32,a/T32,a/)') &
        '  STOKES-'//a1//' , tasks 1-6', &
        '  ===================='
      write(un,'(a,T32,6(a16,2x)/)')'Geometry    Level/Output', &
                                    'No FOCorr No DM','No FOCorr + DM ', &
                                    'FO_ReglrPS + DM','FO_EnhncPs + DM', &
                                    '#4 + Solsaving ','#5 + BVPTelscpe'

      do v = 1, n_geometries
        do n = 1, n_user_levels
          write(un,366)v,'Upwelling @',user_levels(n), (stokes_PT(n,v,o1,upidx,t),t=1,nt)
        ENDdo
	do n = 1, n_user_levels
          write(un,366)v,'Dnwelling @',user_levels(n), (stokes_PT(n,v,o1,dnidx,t),t=1,nt)
        ENDdo
        write(un,*)' '
      ENDdo

!  Write fluxes

      write(un,'(/T32,a/T32,a/)') &
        'Stokes-'//a1//' ACTINIC + REGULAR FLUXES, tasks 1-2', &
        '============================================'
      write(un,'(a,T31,4(a16,2x)/)')' Sun SZA    Level/Output', &
                                    ' Actinic No DM ',' Actinic + DM  ', &
                                    ' Regular No DM ',' Regular + DM  '

      do v = 1, n_szangles
        do n = 1, n_user_levels
          write(un,367)v,'Upwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,upidx,t),t=1,2), &
                (FLUX_PT(n,v,o1,upidx,t),t=1,2)
        ENDdo
	do n = 1, n_user_levels
          write(un,367)v,'Dnwelling @',user_levels(n), &
                (MEANST_PT(n,v,o1,dnidx,t),t=1,2), &
                (FLUX_PT(n,v,o1,dnidx,t),t=1,2)
        ENDdo
        write(un,*)' '
      ENDdo

366   format(i5,T11,a,f6.2,2x,8(1x,1pe16.7,1x))
367   format(i5,T11,a,f6.2,2x,4(1x,1pe16.7,1x))

      END SUBROUTINE write_stokes


END MODULE PCRTM_File_Utility
