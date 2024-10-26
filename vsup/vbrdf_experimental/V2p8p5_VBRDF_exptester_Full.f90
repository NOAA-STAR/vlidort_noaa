program V2p8p5_VBRDF_exptester_Full

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    -- Full BRDF calculations with 36 geometries. Kernels 4/9/10/14/19.
!    -- Tests and outputs timing. Validate against similar tester in Old Code.

!  Pars file

   USE VLIDORT_PARS_m

!  Module files for VBRDF

   USE VBRDF_SUP_AUX_m, Only : VBRDF_READ_ERROR
   USE VBRDF_SUP_MOD_m

   implicit none

!  VBRDF supplement file inputs status structure

   TYPE(VBRDF_Input_Exception_Handling)   :: VBRDF_Sup_InputStatus

!  VBRDF supplement input structures

   TYPE(VBRDF_Sup_Inputs)                 :: VBRDF_Sup_In

!  VBRDF supplement output structure

   TYPE(VBRDF_Sup_Outputs)                :: VBRDF_Sup_Out
   TYPE(VBRDF_Sup_Core)                   :: Core
   TYPE(VBRDF_Sup_Book)                   :: Book
   TYPE(VBRDF_Output_Exception_Handling)  :: VBRDF_Sup_OutputStatus

!  Local Variables
!  ===============

   CHARACTER(LEN=100) :: VBRDF_INPUT_FILE
   INTEGER            :: BS_NMOMENTS_INPUT
   LOGICAL            :: DO_DEBUG_RESTORATION

!  help

   INTEGER            :: nruns, irun, I, m, mruns, choose_kernel, nstokes, nssq, o1, ounit
   DOUBLE PRECISION   :: WINDSPEED, REFIDX, PARS(4)
   REAL               :: e1, e2, e_core(10), e_scale(10), druns,dmruns
   LOGICAL            :: Screen_output = .false.
!   LOGICAL            :: Screen_output = .true.

!  Read Input file

!   VBRDF_INPUT_FILE = 'VBRDF_ReadInput_experimental.cfg'
   VBRDF_INPUT_FILE = 'VBRDF_ReadInput_xiong.cfg'
   CALL VBRDF_INPUTMASTER ( VBRDF_INPUT_FILE, &
          VBRDF_Sup_In,         & ! Outputs
          VBRDF_Sup_InputStatus ) ! Outputs

!  Exception handling

   IF ( VBRDF_Sup_InputStatus%BS_STATUS_INPUTREAD .ne. VLIDORT_SUCCESS ) THEN
      CALL VBRDF_READ_ERROR ( 'VBRDF_ReadInput_experimental.log', VBRDF_Sup_InputStatus )
      STOP 'Stop program after input read failure'
   ENDIF

!  Open file

   OPEN(1,file='toplevel.inp',status='old')
   read(1,*)choose_kernel
   read(1,*)nstokes ; NSSQ = NSTOKES * NSTOKES
   close(1)

!  Fresnel

   WINDSPEED = 5.0d0
   REFIDX    = 1.334d0

!  Local

   BS_NMOMENTS_INPUT = 2*VBRDF_Sup_In%BS_NSTREAMS - 1
   DO_DEBUG_RESTORATION = .false.

!  timing

   nruns = 1000  ; druns = real(nruns)
   mruns = 5    ; dmruns = real(mruns)

!  Initialize

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 0
   VBRDF_Sup_In%BS_BRDF_NAMES     = '          '
   VBRDF_Sup_In%BS_WHICH_BRDF     = 0
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS = 0
   VBRDF_Sup_In%BS_BRDF_PARAMETERS   = ZERO
   VBRDF_Sup_In%BS_LAMBERTIAN_KERNEL_FLAG = .false.
   VBRDF_Sup_In%BS_DO_SEPARATE_FACTORS    = .false.
   VBRDF_Sup_In%BS_BRDF_FACTORS        = ZERO
   VBRDF_Sup_In%BS_BRDF_VECTOR_FACTORS = ZERO

!  Settings
!  ========

if ( choose_kernel .eq. 9 ) then

   PARS(1) = 0.00512d0*WINDSPEED + 0.003d0
   PARS(2) = REFIDX * REFIDX
   PARS(3) = ONE

!  setting

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
   VBRDF_Sup_In%BS_BRDF_NAMES(1)          = 'Cox-Munk  '
   VBRDF_Sup_In%BS_WHICH_BRDF(1)          = 9
   VBRDF_Sup_In%BS_BRDF_FACTORS(1)        = ONE
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)   = 3
   VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1:3) = PARS(1:3)

else if ( choose_kernel .eq. 10 ) then

   PARS(1) = 0.5d0 * ( 0.00512d0*WINDSPEED + 0.003d0 )
   PARS(2) = REFIDX 
   PARS(3) = ONE

!  setting

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
   VBRDF_Sup_In%BS_BRDF_NAMES(1)          = 'GissCoxMnk'
   VBRDF_Sup_In%BS_WHICH_BRDF(1)          = 10
   VBRDF_Sup_In%BS_BRDF_FACTORS(1)        = ONE
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)   = 3
   VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1:3) = PARS(1:3)

else if ( choose_kernel .eq. 14 ) then

   PARS(1) = REFIDX
   PARS(2) = 0.8d0 ! NDVI
   PARS(3) = ONE   ! SCALING

!  setting

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
   VBRDF_Sup_In%BS_BRDF_NAMES(1)          = 'BPDF-NDVI '
   VBRDF_Sup_In%BS_WHICH_BRDF(1)          = 14
   VBRDF_Sup_In%BS_BRDF_FACTORS(1)        = ONE
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)   = 3
   VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1:3) = PARS(1:3)

else if ( choose_kernel .eq. 19 ) then

   PARS(1) = 3.6d0   ! L-parameter (snow grain size)
   PARS(2) = 5.5d0   ! M-parameter
   PARS(3) = 1.02d0  ! Wavelength in Microns

!  setting

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
   VBRDF_Sup_In%BS_BRDF_NAMES(1)          = 'SnowModel '
   VBRDF_Sup_In%BS_WHICH_BRDF(1)          = 19
   VBRDF_Sup_In%BS_BRDF_FACTORS(1)        = ONE
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)   = 3
   VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1:3) = PARS(1:3)

else if ( choose_kernel .eq. 4 ) then

   PARS(1) = 2.0d0   
   PARS(2) = 1.0d0  

!  setting

   VBRDF_Sup_In%BS_N_BRDF_KERNELS = 1
   VBRDF_Sup_In%BS_BRDF_NAMES(1)          = 'Li-sparse '
   VBRDF_Sup_In%BS_WHICH_BRDF(1)          = 4
   VBRDF_Sup_In%BS_BRDF_FACTORS(1)        = ONE
   VBRDF_Sup_In%BS_N_BRDF_PARAMETERS(1)   = 2
   VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,1:2) = PARS(1:2)

else if ( choose_kernel .eq. 0 ) then  ! below is added by Xiong

   CALL VBRDF_INPUTMASTER ( VBRDF_INPUT_FILE, &
          VBRDF_Sup_In,         & ! Outputs
          VBRDF_Sup_InputStatus ) ! Outputs

   print*, " --------------"
   print*,VBRDF_Sup_In%BS_N_BRDF_KERNELS
   print*,VBRDF_Sup_In%BS_BRDF_NAMES(1:VBRDF_Sup_In%BS_N_BRDF_KERNELS)


      DO_DEBUG_RESTORATION = .false.

! A normal calculation will require

      BS_NMOMENTS_INPUT = 2 * VBRDF_Sup_In%BS_NSTREAMS - 1

endif

!  RUNS
!  ====

   do m = 1, mruns

!  Core calculation

   call CPU_TIME(e1)
   do irun = 1, nruns
      CALL VBRDF_CORE_MAINMASTER ( VBRDF_Sup_In, Core, Book )
   enddo
   call CPU_TIME(e2) ; e_Core(m) = (e2 - e1) !/druns

!   write(*,*)'Core Result and CPU Time  = ',Core%DBKERNEL_BRDFUNC_CORE(1,1,1,1), e_core(m)
!   do o1 = 1, nstokes
!      write(*,*)Core%DBKERNEL_BRDFUNC_HELP(1,o1,1,1,1)
!   enddo
!  First run to see if everything is OK
!   CALL VBRDF_MAINMASTER ( &
!        DO_DEBUG_RESTORATION,     & ! Inputs
!        BS_NMOMENTS_INPUT,        & ! Inputs
!        VBRDF_Sup_In, Core, Book, & ! Inputs
!        VBRDF_Sup_Out,            & ! Outputs
!        VBRDF_Sup_OutputStatus )    ! Output Status
!   IF ( VBRDF_Sup_OutputStatus%BS_STATUS_OUTPUT .ne. vlidort_success ) then
!     write(*,'(/A/)') '      **************** FIRST RUN # '
!     write(*,*)'Number of messages from VBRDF MAINMASTER = ',VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
!     do i = 1, VBRDF_Sup_OutputStatus%BS_NOUTPUTMESSAGES
!       write(*,'(a,i2,a,a)')' Message # ', i, ': ',(trim(VBRDF_Sup_OutputStatus%BS_OUTPUTMESSAGES(i))
!     enddo
!     write(*,*)' ' ; Stop ' BRDF calculation failed, stop program'
!   endif

!  Scaling and final calculations

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL VBRDF_MAINMASTER ( &
        DO_DEBUG_RESTORATION,     & ! Inputs
        BS_NMOMENTS_INPUT,        & ! Inputs
        VBRDF_Sup_In, Core, Book, & ! Inputs
        VBRDF_Sup_Out,            & ! Outputs
        VBRDF_Sup_OutputStatus )    ! Output Status
   enddo
   call CPU_TIME(e2) ; e_scale(m) = (e2-e1)
!   call CPU_TIME(e2) ; e_scale(m) = (e2-e1)/druns  ! xiong
   write(*,*)'Final DB-1-1-1 Result and CPU times = ',&
                VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(1,1,1,1), e_core(m), e_scale(m),e_core(m) + e_scale(m)

!  mruns

   enddo

!  Stop if debugging

   if ( mruns.eq.1.and.nruns.eq.1 ) stop 'debugging........'

!  final results ( Only if flagged )

   ounit = 10*choose_kernel + nstokes ; if ( Screen_output ) ounit = 6
   write(ounit,'(/A/)')'DBONLY results (36 geometries) ===> '
   do o1 = 1, nssq
      write(ounit,'(A,i3)')'  -- Reflectance Matrix element # ', o1
      do m = 1, 4
         write(ounit,'(10x,i2,3(3x,3e12.4))')m,VBRDF_Sup_Out%BS_DBOUNCE_BRDFUNC(o1,1:3,1:3,m)
      enddo
   enddo
   write(ounit,'(/A/)')' Fourier results (I/J) ===> '
   do o1 = 1, nssq
      write(ounit,'(A,i3)')'  -- Reflectance Matrix element # ', o1
      do m = 0, BS_NMOMENTS_INPUT
         write(ounit,'(10x,i2,4(3x,4e12.4))')m,VBRDF_Sup_Out%BS_BRDF_F(m,o1,1:4,1:4)
      enddo
   enddo
   write(ounit,'(/A/)')' Fourier results (SUN/J) ===> '
   do o1 = 1, nssq
      write(ounit,'(A,i3)')'  -- Reflectance Matrix element # ', o1
      do m = 0, BS_NMOMENTS_INPUT
         write(ounit,'(10x,i2,4(3x,4e12.4))')m,VBRDF_Sup_Out%BS_BRDF_F_0(m,o1,1:4,1:4)
      enddo
   enddo

!  Final average timing

   if ( mruns.gt.1 ) then
     write(*,'(/A,3F12.6/)')'   @@@ Average 5-run times = ',sum(e_core(1:mruns))/dmruns,sum(e_scale(1:mruns))/dmruns,&
                  (sum(e_core(1:mruns))/dmruns)+(sum(e_scale(1:mruns))/dmruns)
   endif

!  done

   stop
end program V2p8p5_VBRDF_exptester_Full


