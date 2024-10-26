MODULE PCRTM_Sw_Aux_m

! By Xiaozhen Xiong, 6/10/2022

! ###############################################################
! #                                                             #
! # SUBROUTINEs in this Module                                  #
! #                                                             #
! #            Generate_Random_Inputs                           #
! #            Get_aircraft_outlevel                            #
! #            RAYLEIGH_FUNCTION  				#
! #            Write_IQU    					#
! #            Write_IQU_aircraft				#
! #            Write_3components	    	       		#
! #            Write_Inputs		    	       		#
! #            Write_Stokes    					#
! #            Add output
! #                                                             #
! ###############################################################

USE VLIDORT_PARS_m, Only :  MAX_SZANGLES


! SETUP SOME CONSTANT VARIABLES (from old driver of ming, useful ?)
    REAL(8), parameter ::  RHO_STAND = 2.68675D+19!#/cm3
    REAL(8), parameter ::  PZERO     = 1013.25D0
    REAL(8), parameter ::  TZERO     = 273.15D0
    REAL(8), parameter ::  RHO_ZERO  = RHO_STAND * TZERO / PZERO
     ! 10^5 is for unit change, to calculate column density
     ! heights_jpl is [km], = 10^5 cm
    REAL(8), parameter ::  CONST     = 1.0D+05 * RHO_ZERO
    REAL(8), parameter ::  O2RATIO   = 0.2095D0
     !CO2 PPMV mixing ratio (for the Rayleigh stuff)
!    DOUBLE PRECISION  CO2_PPMV_MIXRATIO
!    PARAMETER        ( CO2_PPMV_MIXRATIO = 390.0d0 )
INTEGER, parameter :: maxtypeaer=4

CONTAINS


SUBROUTINE   read_user_control(DO_mapgrid,DO_usrgrid,         & !pres grid
                     ATM_TYPE,Psfc_in,              & !wave grid
                     DO_Lambertian,LAMBERTIAN_ALBEDO,                  & !surface
                     Do_Planetary,                                     & !planetary
                     DO_Scalar_only,userset_stream,           & !stokes and stream 
                     N_USER_OBSGEOMS,&
                     szangles,user_relazms,user_vzangles)
                     


logical :: DO_mapgrid,DO_usrgrid,DO_stdgrid
logical :: Do_LUT_Allspectral,do_usr_wv
INTEGER :: ATM_TYPE,userset_stream,N_USER_OBSGEOMS
REAL(8), dimension(1:MAX_SZANGLES) :: szangles,user_relazms,user_vzangles
logical::  DO_Aerosol, DO_opac, DO_tamu,do_icecld !aerosol, cloud
logical::  DO_Lambertian, Do_Planetary            !surface and planetary
logical::  DO_Obsgeoms, Do_Doublets 
logical::  DO_Scalar_only,DO_User_Stream          !stokes and stream 
CHARACTER(LEN=300) :: filename
REAL   ::   Psfc_in
REAL(8)::   LAMBERTIAN_ALBEDO
   filename ='./PCRTM_VLIDORT_Config/user_control.cfg'
   print*,filename
   CALL GETLUN(IFILE)
   open ( unit = IFILE, file = filename,action='read',status='old')
   rewind ( IFILE )
   read(IFILE,*)
   read(IFILE,*)
   read(IFILE,*) DO_mapgrid
   read(IFILE,*)
   read(IFILE,*) DO_usrgrid
   read(IFILE,*)
   read(IFILE,*) ATM_TYPE
   read(IFILE,*)
   read(IFILE,*) do_icecld
   read(IFILE,*)
   read(IFILE,*) DO_Lambertian
   read(IFILE,*)
   read(IFILE,*) Do_Planetary
   read(IFILE,*)
   read(IFILE,*) DO_Scalar_only
   read(IFILE,*)
   read(IFILE,*) userset_stream
   read(IFILE,*)
   read(IFILE,*) N_USER_OBSGEOMS
   read(IFILE,*)
   do n =1, N_USER_OBSGEOMS
     read(IFILE,*) szangles(n),user_vzangles(n),user_relazms(n)
   ENDdo
   read(IFILE,*)
   read(IFILE,*) Psfc_in
   read(IFILE,*)
   read(IFILE,*) lambertian_albedo
   close ( IFILE )



END SUBROUTINE read_user_control


SUBROUTINE  Read_LambSfc_Alb (path0,lambSfcRef, SfcRef)

   implicit none


     CHARACTER(LEN=80) :: path0

     REAL, dimension(1:2152, 1:10) :: SfcRef
     REAL, dimension(1:2152)         :: lambSfcRef         !,SfcRef1D
     INTEGER :: IFile, n

! Read SfcRef Qiguang data
      CALL GETLUN(IFILE)

!     open ( unit = IFILE, file = trim(path0)//'data_and_control/surface_reflectance_sub.bin', action='read',&
!            status ='old',form='unformatted',access='stream')
!     read (IFILE) SfcRef
!
!      CLOSE(IFILE)

! SET  WAVENUM/LAMDAS FOR SfcRef
      DO n = 1, 2151
        lambSfcRef(n+1) = 0.350 + REAL(n-1) * 0.001
      ENDDO
      lambSfcRef(1) =0.3

END SUBROUTINE Read_lambSfc_Alb


  SUBROUTINE READ_UserProf_Geo_INPUTS(FILNAM,  &
              temp_user, H2o_user,O3_user)

      USE VLIDORT_PARS_m, Only :  MAXLAYERS, MAX_SZANGLES,   &
                                 MAX_MESSAGES, MAX_USER_RELAZMS, MAX_USER_VZANGLES, MAX_USER_LEVELS, &
                                 VLIDORT_SUCCESS, VLIDORT_SERIOUS, VLIDORT_INUNIT, ONE

      USE VLIDORT_Inputs_def_m

      USE VLIDORT_AUX_m , Only : GFINDPAR, FINDPAR_ERROR

      IMPLICIT NONE


   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15), Nlev =101


    REAL*8, dimension(1:5200)       :: user_sza
    REAL*8, dimension(1:5200)  :: user_vza
    REAL*8, dimension(1:5200)   :: user_azm

    REAL, dimension(1:Nlev)   :: pres_user,temp_user, h2o_user,O3_user           ! change to variable layer later

    REAL                 :: latitude, longitude,psfc
    REAL*8               :: alb

     CHARACTER (LEN=12), PARAMETER :: PREFIX = 'USER INPUT -'
      LOGICAL            :: ERROR
      CHARACTER (LEN=80) :: PAR_STR
      INTEGER            :: FILUNIT, NM, I,n_szas,n_vzas,n_azms, N_geos
      INTEGER, dimension(1:nlev)        :: level_prof

      CHARACTER (LEN=*) :: FILNAM
      INTEGER  ::                STATUS
      INTEGER  ::              NMESSAGES
      CHARACTER (LEN=120) ::    MESSAGES ( 0:MAX_MESSAGES )
      CHARACTER (LEN=120) ::    ACTIONS  ( 0:MAX_MESSAGES )

      CHARACTER(len=200),dimension(3)   :: message_Loading

     CALL GETLUN(FILUNIT)
      OPEN(FILUNIT,FILE=FILNAM,ERR=300,STATUS='OLD')

       nm = 0

      PAR_STR = 'Latitude'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) latitude
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Longitude'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) longitude
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Number of Solar Zenith Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_szas
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Solar Zenith Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,N_szas
         READ (FILUNIT,*,ERR=998) user_sza(i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Number of Satellite Viewing Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_vzas
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
      PAR_STR = 'Satellite Viewing Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,N_vzas
         READ (FILUNIT,*,ERR=998) user_vza(i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Number of Relative Azimuthal Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_azms
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Relative Azimuthal Angle'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,N_azms
         READ (FILUNIT,*,ERR=998) user_azm(i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Number of Observation Geometry inputs'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) N_Geos
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Observation Geometry inputs'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,N_Geos
         READ (FILUNIT,*,ERR=998) user_sza(i),user_vza(i),user_azm(i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Surface Albedo'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) alb
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

      PAR_STR = 'Surface Pressure'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) READ (FILUNIT,*,ERR=998) psfc
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )
     PAR_STR = 'Pressure'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,Nlev
         READ (FILUNIT,*,ERR=998) level_prof(i), pres_user(i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

     PAR_STR = 'Temperature Profile'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,Nlev
         READ (FILUNIT,*,ERR=998) level_prof(Nlev+1-i),temp_user(Nlev+1-i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

     PAR_STR = 'Water Vapor Profile (Kg/Kg)'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,Nlev
         READ (FILUNIT,*,ERR=998) level_prof(Nlev+1-i),H2o_user(Nlev+1-i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )


     PAR_STR = 'Ozone Profile'
      IF (GFINDPAR ( FILUNIT, PREFIX, ERROR, PAR_STR)) then
       DO I = 1,Nlev
         READ (FILUNIT,*,ERR=998) level_prof(Nlev+1-i),O3_user(Nlev+1-i)
       ENDDO
      ENDIF
      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )


300   CONTINUE
      STATUS = VLIDORT_SERIOUS
      NMESSAGES = NMESSAGES + 1
      MESSAGES(NMESSAGES) = 'openfile failure for '//trim(adjustl(FILNAM))
      ACTIONS(NMESSAGES)  = 'Find the Right input file!!'
      CLOSE(FILUNIT)
998   CONTINUE
      NM = NM + 1
      STATUS       = VLIDORT_SERIOUS
      MESSAGES(NM) = 'Read failure for entry below String: ' //Trim(Adjustl(PAR_STR))
      ACTIONS(NM)  = 'Re-set value: Entry wrongly formatted in Input file'
      NMESSAGES    = NM

!  Finish

      RETURN

      CALL FINDPAR_ERROR ( ERROR, PAR_STR, STATUS, NM, MESSAGES, ACTIONS )

    END SUBROUTINE READ_UserProf_Geo_INPUTS


 SUBROUTINE Assign_LayOpt_to_Vlidort(w,nlayers,ndat,nstokes,n_geoms, MAXLAYERS, MAX_GEOMETRIES,    & ! INPUTS
                            ncoeffs, MAXMOMENTS_INPUT,MAXSTOKES_SQ, airClm_vld,  RAYLEIGH_XSEC,&        ! INPUTS
                            RayCoeffs, RayFmatrices_up,RayFmatrices_dn, &                          ! INPUTS
                            flag_aer,aerext_intp,aersca_intp, &                          ! INPUTS
                            odIOTot_vld,                      &                          ! INPUTS
                            OutFmat_up_intp,OutFmat_dn_intp,FMatCoeffs_intp,&              ! INPUTS
                            deltau_vert_input,omega_total_input,&          ! OUTPUTS
                            GREEKMAT_TOTAL_INPUT,fmatrix_up,fmatrix_dn)  ! OUTPUTS

        implicit none

! INPUTS
       INTEGER :: nlayers,nstokes,n_geoms,MAXLAYERS, MAX_GEOMETRIES, MAXMOMENTS_INPUT,MAXSTOKES_SQ,ngreekmat_entries
!ming0801
       INTEGER :: flag_aer(maxtypeaer,MAXLAYERS)
       INTEGER :: ncoeffs,ndat
       DOUBLE PRECISION  :: airClm_vld(1:nlayers)
!ming0801 raycoeff define 
       DOUBLE PRECISION  :: RAYLEIGH_XSEC(1:ndat), RayCoeffs(0:2,6)!airClm_vld(:), RAYLEIGH_XSEC(:), RayCoeffs(:,:)
       !REAL              :: odIOTot_vld(:,:)
       DOUBLE PRECISION  :: odIOTot_vld(1:nlayers,1:ndat)
       DOUBLE PRECISION :: RayFmatrices_up   ( MAX_GEOMETRIES, 6 )
       DOUBLE PRECISION :: RayFmatrices_dn   ( MAX_GEOMETRIES, 6 )
!ming0801
       DOUBLE PRECISION :: aerext_intp(MAXLAYERS), aersca_intp(MAXLAYERS)
!ming0801
       DOUBLE PRECISION :: OutFmat_up_intp(MAXLAYERS,MAX_GEOMETRIES,6)
       DOUBLE PRECISION :: OutFmat_dn_intp(MAXLAYERS,MAX_GEOMETRIES,6)
       DOUBLE PRECISION :: FMatCoeffs_intp(MAXLAYERS,0:ncoeffs,6)!FMatCoeffs_intp(:,0:,:)

       DOUBLE PRECISION :: total_molsca,total_moltau
       DOUBLE PRECISION ::  aerext,  aersca,  molsca, totsca, raywt,aerwt, totext

!OUTPUTS
       DOUBLE PRECISION :: deltau_vert_input(nlayers),omega_total_input(nlayers)
       DOUBLE PRECISION :: GREEKMAT_TOTAL_INPUT ( 0:MAXMOMENTS_INPUT, MAXLAYERS, MAXSTOKES_SQ )
       DOUBLE PRECISION :: fmatrix_up( MAXLAYERS, MAX_GEOMETRIES, 6 ),fmatrix_dn( MAXLAYERS, MAX_GEOMETRIES, 6 )

! LOCAL : 

      INTEGER ::          K, CK, GK, RK, FK,CMASK(8), GMASK(8), SMASK(8),RMASK(8),FMASK(8), n,w, kk,ll, sk

      IF ( nstokes .eq. 1 ) ngreekmat_entries = 1
      IF ( nstokes .eq. 3 ) ngreekmat_entries = 5
      IF ( nstokes .eq. 4 ) ngreekmat_entries = 8


!------ below is from the V2p8p3_vfzmat_tester.f90, 5/1/2023
      GMASK = (/  1, 2, 5, 6, 11, 12, 15, 16 /)
      SMASK = (/  1, -1, -1, 1, 1, -1, 1, 1 /)
      FMASK = (/  1, 2, 2, 3, 4, 5, 5, 6 /)    !  Fmatrix filling

!  here, we are using output from the vfzmat supplement, so the same quantity.

      RMASK = (/  1, 5, 5, 2, 3, 6, 6, 4 /)    !  This   Rayleigh
      CMASK = RMASK   !  OUTPUT OF RTS MIE
!--------------

   DO n=1, nlayers
       !gas layer: assign opt prop
       IF (sum(flag_aer(:,n)) .eq. 0) THEN

            total_molsca = airClm_vld(n) * RAYLEIGH_XSEC(w)

            !ODIO need to set same vertical direction 1 is toa

            total_moltau = total_molsca + dble(odIOTot_vld(n,w))

            deltau_vert_input(n) = total_moltau

            IF(total_moltau .gt. 10) THEN
              !print*,'rayleigh',n, airClm_vld(n),RAYLEIGH_XSEC(w),total_moltau,odIOTot_vld(n,w)
            ENDIF

            omega_total_input(n) = total_molsca /total_moltau
            IF ( omega_total_input(n)  .gt. 0.999999d0) omega_total_input(n)  =0.999999d0
            IF ( omega_total_input(n)  .lt. 1.d-6)    omega_total_input(n)  = 1.d-6


            DO K = 1, ngreekmat_entries
              GK = gmask(k) ; rk = rmask(k) !; fk = fmask(k)
              GREEKMAT_TOTAL_INPUT(0:2,n,gk) = RayCoeffs(0:2,rk)
              GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
            ENDDO
            !  Fmatrices (efficient). Only need the first two entries with natural sunlight.
            DO k = 1, 2
              rk = rmask(k) ; fk = fmask(k)
              fmatrix_up(n,1:n_geoms,fk) = RayFmatrices_up(1:n_geoms,rk)
              fmatrix_dn(n,1:n_geoms,fk) = RayFmatrices_dn(1:n_geoms,rk)
            ENDDO

       !aerosol layer: assign opt prop
       ELSE
            aersca = aersca_intp(n)
            aerext = aerext_intp(n)
            molsca = airClm_vld(n) * RAYLEIGH_XSEC(w)

            totsca = molsca  + aersca
            raywt  = molsca / totsca
            aerwt  = aersca / totsca

            totext = molsca+ aerext+odIOTot_vld(n,w)!(49-n+1,w)
            deltau_vert_input(n) = totext
            omega_total_input(n) = totsca / totext


            IF ( omega_total_input(n)  .gt. 0.999999d0) omega_total_input(n)=0.999999d0
            IF ( omega_total_input(n)  .lt. 1.d-6)    omega_total_input(n)  =1.d-6

            ! below is the new one, 5/3/2023
            DO K = 1, ngreekmat_entries
              GK = gmask(k); sk = dble(smask(k)); ck = cmask(k) ; rk = rmask(k) ; fk = fmask(k)
              greekmat_total_input(0:2,n,gk) = raywt * RayCoeffs(0:2,rk) + sk * aerwt * FmatCoeffs_intp(n,0:2,ck)
              greekmat_total_input(3:ncoeffs,n,gk) = sk * aerwt * FmatCoeffs_intp(n,3:ncoeffs,ck)
            ENDDO
            GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0

            !  Fmatrices (efficient). Only need the first two entries with natural sunlight.
            DO k = 1, 2
              sk = dble(smask(k)); ck = cmask(k) ;rk = rmask(k) ; fk = fmask(k)
              fmatrix_up(n,1:n_geoms,fk) = raywt * RayFmatrices_up(1:n_geoms,rk) + sk * aerwt * OutFmat_up_intp(n,1:n_geoms,ck)
              fmatrix_dn(n,1:n_geoms,fk) = raywt * RayFmatrices_dn(1:n_geoms,rk) + sk * aerwt * OutFmat_dn_intp(n,1:n_geoms,ck)
            ENDDO
       ENDIF
       !rev: ming Must setting here to avoid input checking error
       GREEKMAT_TOTAL_INPUT(0,n,1) = 1.0d0
    ENDDO !end of layer loop

 END SUBROUTINE Assign_LayOpt_to_Vlidort



 

SUBROUTINE Generate_Random_Aer_Inputs(seed_times, &       !rand_az, rand_sz, rand_vz, rand_prof, rand_psfc,rand_Ref, &
           rand_aerlowbound,rand_aerthickness,rand_aerAOD, rand_aertype, rand_brdftype,rand_scaleH, &
           rand_dustlowbound,rand_dustthickness,rand_dustAOD,rand_dusttype,rand_dustSph)


   USE vlidort_pars_m, Only : MAX_GEOMETRIES, MAX_USER_VZANGLES, MAX_USER_RELAZMS, MAX_SZANGLES

   implicit none

   INTEGER, PARAMETER :: dpk = SELECTED_REAL_KIND(15)


   INTEGER   ::  seed

    REAL, dimension(9690)                :: rProf
    INTEGER                              :: n, seed_times

!aeresol random number
    REAL,dimension(9690) ::rand_aerlowbound,rand_aerthickness,rand_aerAOD
    REAL,dimension(9690) ::rand_dustlowbound,rand_dustthickness,rand_dustAOD
    INTEGER,dimension(9690) ::rand_aertype, rand_brdftype,rand_dusttype
!ming09/12/2022 add random scale height for exp. distribution
    INTEGER,dimension(9690) ::rand_scaleH
    REAL,dimension(9690) :: rand_dustSph
    INTEGER :: Max_prof
    INTEGER :: icemod_ind
!  setups

    !icetest, mannually setting (change to auto setting in driver?)
    icemod_ind = 4000   !4000 thm; 5000 mc5; 6000 mc6

    Max_prof = 9690

    seed = 760013 ;   seed = seed + seed_times



   !!ming rev genereate random aerosol profile propertites
     seed =  961028 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
   !ming 09/02/2022 setting aer low boundary 0.~1km
     rand_aerlowbound=(rProf*1.0)


     seed =  961029 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
  !ming 09/02/2022 setting aer low boundary 3.0~6km
     rand_aerthickness=(rProf*3.0+3.0)


     seed =  961030 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
     !rand_aerAOD=(rProf*5)
     !rand_aerAOD=(rProf*0.4+0.1)   !Ming 09/01/2022
     rand_aerAOD=(rProf*1.95+0.05)   !Ming 09/16/2022

     seed =  961031 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
     !rand_aertype=mod(ceiling(rProf*90),10)*10
     rand_aertype=ceiling(rProf*10)+1000 ! Ming 09/01/2022 select Modtran Aertype 1001~1010

     seed =  961032 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
!     rand_brdftype=mod(ceiling(rProf*19),10)*10

     rand_brdftype=ceiling(rProf*Max_prof)   ! Xiong, 12/8/2022


    seed =  961033 ;   seed = seed + seed_times
     do n = 1, Max_prof
       rProf(n) = ran ( seed )
     ENDdo
     rand_scaleH=ceiling(rProf*5+1)

!ming 11-01-2022 add dust random
     seed =  961034 ;   seed = seed + seed_times
     do n = 1, 9690
       rProf(n) = ran ( seed )
     ENDdo
     
     if (icemod_ind .eq. 4000) rand_dusttype=ceiling(rProf*67)+4001  !THM ice model
     if (icemod_ind .eq. 5000) rand_dusttype=ceiling(rProf*47)+5001  !MC5 ice model
     if (icemod_ind .eq. 6000) rand_dusttype=ceiling(rProf*83)+6001  !MC6 ice model

     seed =  961035 ;   seed = seed + seed_times
     do n = 1, 9690
       rProf(n) = ran ( seed )
     ENDdo
     !rand_aerAOD=(rProf*5)
     !rand_dustAOD=(rProf*10.+2.0)   !dust 2~12
     rand_dustAOD=(rProf*60.+0.05)   !icetest increase ice AOD to 60 

     seed =  961036 ;   seed = seed + seed_times
     do n = 1, 9690
       rProf(n) = ran ( seed )
     ENDdo
   !ming 09/02/2022 setting aer low boundary 0.~0.1km
     !rand_dustlowbound=(rProf*0.1)
   !ming icetest change to 5~10 (before is 10-12km)
     !rand_dustlowbound=(rProf*5.+5.)
   !ming chnage lower boundary 12-25km
     rand_dustlowbound=(rProf*12.+13.)

     seed =  961037 ;   seed = seed + seed_times
     do n = 1, 9690
       rProf(n) = ran ( seed )
     ENDdo
  !ming 09/02/2022 setting aer low boundary 3.0~6km
     !rand_dustthickness=(rProf*3.0+6.0)
  !ming ice test
     rand_dustthickness=(rProf*2.0+0.3)

     seed =  961038 ;   seed = seed + seed_times
     do n = 1, 9690
       rProf(n) = ran ( seed )
     ENDdo
     rand_dustSph=ceiling(rProf*10)*0.01+0.68
END SUBROUTINE Generate_Random_Aer_Inputs




SUBROUTINE RAYLEIGH_FUNCTION &
          ( FORWARD_MAXLAMBDAS, CO2_PPMV_MIXRATIO, &
            FORWARD_NLAMBDAS,   FORWARD_LAMBDAS, &
            RAYLEIGH_XSEC, RAYLEIGH_DEPOL )

!  Rayleigh cross sections and depolarization ratios
!     Bodhaine et. al. (1999) formulae
!     Module is stand-alone.

      IMPLICIT NONE

!  Input arguments
!  ---------------

!  wavelength

      INTEGER          FORWARD_MAXLAMBDAS, FORWARD_NLAMBDAS
      DOUBLE PRECISION FORWARD_LAMBDAS ( FORWARD_MAXLAMBDAS )

!  CO2 mixing ratio

      DOUBLE PRECISION CO2_PPMV_MIXRATIO

!  Output arguments
!  ----------------

!  cross-sections and depolarization output

      DOUBLE PRECISION RAYLEIGH_XSEC  ( FORWARD_MAXLAMBDAS )
      DOUBLE PRECISION RAYLEIGH_DEPOL ( FORWARD_MAXLAMBDAS )

!  Local variables
!  ---------------

      INTEGER          W
      DOUBLE PRECISION MASS_DRYAIR
      DOUBLE PRECISION NMOL, PI, CONS
      DOUBLE PRECISION MO2,MN2,MARG,MCO2,MAIR
      DOUBLE PRECISION FO2,FN2,FARG,FCO2,FAIR
      DOUBLE PRECISION LAMBDA_C,LAMBDA_M,LPM2,LP2
      DOUBLE PRECISION N300M1,NCO2M1,NCO2
      DOUBLE PRECISION NCO2SQ, NSQM1,NSQP2,TERM
      DOUBLE PRECISION S0_A, S0_B
      DOUBLE PRECISION S1_A, S1_B, S1_C, S1_D, S1_E
      DOUBLE PRECISION S2_A
      DOUBLE PRECISION S3_A, S3_B, S3_C, S3_D, S3_E

! data statements and parameters
!  ------------------------------

      DATA MO2  / 20.946D0 /
      DATA MN2  / 78.084D0 /
      DATA MARG / 0.934D0 /

      PARAMETER        ( S0_A = 15.0556D0 )
      PARAMETER        ( S0_B = 28.9595D0 )

      PARAMETER        ( S1_A = 8060.51D0 )
      PARAMETER        ( S1_B = 2.48099D+06 )
      PARAMETER        ( S1_C = 132.274D0 )
      PARAMETER        ( S1_D = 1.74557D+04 )
      PARAMETER        ( S1_E = 39.32957D0 )

      PARAMETER        ( S2_A = 0.54D0 )

      PARAMETER        ( S3_A = 1.034D0 )
      PARAMETER        ( S3_B = 3.17D-04 )
      PARAMETER        ( S3_C = 1.096D0 )
      PARAMETER        ( S3_D = 1.385D-03 )
      PARAMETER        ( S3_E = 1.448D-04 )

!  Start of code
!  -------------

!  constants

      NMOL = 2.546899D19
      PI   = DATAN(1.0D0)*4.0D0
      CONS = 24.0D0 * PI * PI * PI

!  convert co2

      MCO2 = 1.0D-06 * CO2_PPMV_MIXRATIO

!  mass of dry air: Eq.(17) of BWDS

      MASS_DRYAIR = S0_A * MCO2 + S0_B

!  start loop


      DO W = 1, FORWARD_NLAMBDAS

!  wavelength in micrometers

      LAMBDA_M = 1.0D-03 * FORWARD_LAMBDAS(W)
      LAMBDA_C = 1.0D-07 * FORWARD_LAMBDAS(W)
      LPM2     = 1.0D0 / LAMBDA_M / LAMBDA_M

!  step 1: Eq.(18) of BWDS

      N300M1 = S1_A + ( S1_B / ( S1_C - LPM2 ) ) + &
                      ( S1_D / ( S1_E - LPM2 ) )
      N300M1 = N300M1 * 1.0D-08

!  step 2: Eq.(19) of BWDS

      NCO2M1 = N300M1 * ( 1.0D0 + S2_A * ( MCO2  - 0.0003D0 ) )
      NCO2   = NCO2M1 + 1
      NCO2SQ = NCO2 * NCO2

!  step 3: Eqs. (5&6) of BWDS (Bates' results)

      FN2  = S3_A + S3_B * LPM2
      FO2  = S3_C + S3_D * LPM2 + S3_E * LPM2 * LPM2

!  step 4: Eq.(23) of BWDS
!     ---> King factor and depolarization ratio

      FARG = 1.0D0
      FCO2 = 1.15D0
      MAIR = MN2 + MO2 + MARG + MCO2
      FAIR = MN2*FN2 + MO2*FO2 + MARG*FARG + MCO2*FCO2
      FAIR = FAIR / MAIR
      RAYLEIGH_DEPOL(W) = 6.0D0*(FAIR-1.0D0)/(3.0D0+7.0D0*FAIR)

! step 5: Eq.(22) of BWDS
!     ---> Cross section

      LP2  = LAMBDA_C * LAMBDA_C
      NSQM1 = NCO2SQ - 1.0D0
      NSQP2 = NCO2SQ + 2.0D0
      TERM = NSQM1 / LP2 / NMOL / NSQP2
      RAYLEIGH_XSEC(W) =  CONS * TERM * TERM * FAIR

!  END loop

      ENDDO

!  finish
!  ------

      RETURN
END SUBROUTINE RAYLEIGH_FUNCTION

END MODULE PCRTM_Sw_Aux_m



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
       output_filename = 'Output_IQU/I.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 2) then
       output_filename = 'Output_IQU/IQ.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 3) then
       output_filename = 'Output_IQU/IQU.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else
       output_filename = 'Output_IQU/IQUV.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     ENDif

     print*,output_filename

      CALL GETLUN(IFILE)

      open(UNIT=IFILE, file= output_filename,form='unformatted',&
               access='direct', recl=Nrec_len*8)

      K = 0    
      DO v = 1, n_geometries

        Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
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
       output_filename = 'Output_IQU/Planetary_I.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 2) then
       output_filename = 'Output_IQU/Planetary_IQ.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else if (Nstokes .eq. 3) then
       output_filename = 'Output_IQU/Planetary_IQU.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     else
       output_filename = 'Output_IQU/Planetary_IQUV.'//trim(wnStartStr)//'_'//trim(wnEndStr)//'.'//trim(nProf_s)//'.'//sub_name(1:k-4)//'_'//trim(date)//'.dat'
     ENDif

       print*,output_filename

        CALL GETLUN(IFILE)

        open(UNIT=IFILE, file= output_filename,form='unformatted',&
               access='direct', recl=Nrec_len*8)

      Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
                    DBLE(out_szas(1)), DBLE(out_vzas(1)), DBLE(out_azms(1)),DBLE(st_wvnum_indx), DBLE(END_wvnum_indx)/)

      Out_Tmp(101:npt+100) = wavenums(1:ncount)
       write(IFILE, rec=1) Out_tmp(1:Nrec_len)

      k = 1 
! added 12-01-2023, add output of wavenumbers
      DO v = 1, n_geometries

        Out_tmp(1:16)=(/DBLE(npt), DBLE(Nrec_len),DBLE(nstokes),DBLE(n_geometries),wnStart,wnEnd,wresol,nstep*1.D0,nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
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
        Out_tmp(1:11)=(/DBLE(ngas_lut),DBLE(nlayers),DBLE(nstokes),DBLE(n_geometries),nstreams*1.D0,nprof*1.D0,Psfc_vld*1.D0, &
                     Indx_lambertian*1.D0, indx_brdf*1.D0, indx_bpdf*1.D0, Lambertian_Albedo/)      
        Out_tmp(12:12+12)=(/dble(nmode_aer_in),AOD_aer_in,Aerosol_lowerboundary_in,Aerosol_upperboundary_in/)
        Out_tmp(25:25+15)=(/dble(Aerosol_prof_in),expAerosol_relaxation_in,gdfAerosol_peakheight_in,gdfAerosol_halfwidth_in/)
        Out_tmp(41:44)= DBLE(type_aer_in)
        Out_tmp(45:45+n_geometries*3-1)= (/DBLE(out_szas(1:n_geometries)), DBLE(out_vzas(1:n_geometries)),DBLE(out_azms(1:n_geometries))/)
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





     SUBROUTINE quadterp (xtab,ytab,xint,yint)
          implicit none
          REAL, dimension(0:2)   :: xtab, ytab
          REAL                   :: xint,yint
          REAL                   :: x0,x1,x2
          REAL                   :: y0,y1,y2

          yint = ( (xint-xtab(1)) * ( xint-xtab(2)) )/ ( (xtab(0)-xtab(1)) * (xtab(0)-xtab(2)) ) * ytab(0) + &
                 ( (xint-xtab(0)) * ( xint-xtab(2)) )/ ( (xtab(1)-xtab(0)) * (xtab(1)-xtab(2)) ) * ytab(1) + &
                 ( (xint-xtab(0)) * ( xint-xtab(1)) )/ ( (xtab(2)-xtab(0)) * (xtab(2)-xtab(1)) ) * ytab(2)


    END SUBROUTINE quadterp


    SUBROUTINE quadterp2 (xtab,ytab,xint,yint)
          implicit none
          REAL*8, dimension(0:2)   :: xtab, ytab
          REAL*8                   :: xint,yint
          REAL*8                   :: x0,x1,x2
          REAL*8                   :: y0,y1,y2

          yint = ( (xint-xtab(1)) * ( xint-xtab(2)) )/ ( (xtab(0)-xtab(1)) * (xtab(0)-xtab(2)) ) * ytab(0) + &
                 ( (xint-xtab(0)) * ( xint-xtab(2)) )/ ( (xtab(1)-xtab(0)) * (xtab(1)-xtab(2)) ) * ytab(1) + &
                 ( (xint-xtab(0)) * ( xint-xtab(1)) )/ ( (xtab(2)-xtab(0)) * (xtab(2)-xtab(1)) ) * ytab(2)


    END SUBROUTINE quadterp2






