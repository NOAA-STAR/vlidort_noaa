program V2p8p5_VBRDF_exptester_Konly

!  5/5/22. Version 2.8.5. Experimental Version for NASA-LARC PCRTM team
!    -- Single kernel calculations with 1 geometry. Coxmunk and GissCoxMunk only.
!    -- Compares performance of Old kernel (complete) vs New kernel (Core/scaling)

   use vlidort_pars_m
   use vbrdf_sup_kernels_experimental_m, Only : COXMUNK_CORE, GISSCOXMUNK_CORE, ADD_FRESNEL_SCALING
   use vbrdf_sup_kernels_m             , Only : COXMUNK_VFUNCTION, GISSCOXMUNK_VFUNCTION
   implicit none

   INTEGER, parameter :: MAXPARS = 4
   INTEGER            :: I, J, NPARS, NSSQ, nruns, irun, choose_kernel, NSTOKES
   DOUBLE PRECISION   :: PARS ( MAXPARS )
   DOUBLE PRECISION   :: XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI
   DOUBLE PRECISION   :: KERNELCORE, KERNELHELP(7), KERNEL_Old(16), KERNEL_New(16)
   DOUBLE PRECISION   :: WINDSPEED, REFIDX, DRUNS
   REAL	              :: e1, e2, e3, e4, e_old, e_new

!  Open file

   OPEN(1,file='toplevel.inp',status='old')
   read(1,*)choose_kernel
   read(1,*)nstokes ; NSSQ = NSTOKES * NSTOKES
   close(1)

!  Geometry

   PHI   = 10.0d0
   CPHI  = cos(Phi*DEG_TO_RAD)
   SKPHI = sin(Phi*DEG_TO_RAD)
   XI = COS(35.0d0 * DEG_TO_RAD) ; SXI = sqrt(ONE*ONE-XI*XI)
   XJ = COS(40.0d0 * DEG_TO_RAD) ; SXJ = sqrt(ONE*ONE-XJ*XJ)

!  Fresnel

   WINDSPEED = 5.0d0
   REFIDX    = 1.334d0

!  others

   NPARS = 4 ; PARS = ZERO

!  timing

   nruns = 10000 ; druns = dble(nruns)

!  CoxMunk # 9 
!  ===========

if ( choose_kernel .eq. 9 ) then

   PARS(1) = 0.00512d0*WINDSPEED + 0.003d0
   PARS(2) = REFIDX * REFIDX
   PARS(3) = ONE

!  Old

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL COXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, nssq, XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI, &
           KERNEL_Old )
   enddo
   call CPU_TIME(e2) ; e_old = e2 - e1

   write(*,*)'Old Result and CPU Time  = ',KERNEL_OLD(1), e_old
   write(*,*) nruns
!  New, Core run

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL COXMUNK_CORE &
         ( MAXPARS, NPARS, PARS, XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI, &
           KERNELCORE, KERNELHELP )
   enddo
!   call CPU_TIME(e2) ; e3 = (e2-e1)/druns
   call CPU_TIME(e2) ; e3 = (e2-e1)  ! xiong

!  New,  scaling

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL ADD_FRESNEL_SCALING ( REFIDX, .true., 1, KERNELCORE, KERNELHELP, KERNEL_New )
   enddo
   call CPU_TIME(e2) ; e4 = (e2-e1) ; e_new = e3 + e4
   write(*,*)'New Result, Time/Speedup = ',KERNELCORE, KERNEL_New(1), e_new, e_old/e_new, e4/e3

endif

!  Giss Cox Munk
!  =============

if ( choose_kernel .eq. 10 ) then

   PARS(1) = 0.5d0 * ( 0.00512d0*WINDSPEED + 0.003d0 )
   PARS(2) = REFIDX 
   PARS(3) = ONE

!  Old

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL GISSCOXMUNK_VFUNCTION &
         ( MAXPARS, NPARS, PARS, NSSQ, XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI, &
           KERNEL_Old )
   enddo
   call CPU_TIME(e2) ; e_old = e2 - e1
   write(*,'(/A,1pe17.7,2x,F12.7)')'Old Result and CPU Time  = ', KERNEL_Old(1), e_old
   do I = 1, NSTOKES
      write(*,'(i2,2x,1p4e17.7)')I,(KERNEL_OLD(4*(I-1)+J),j=1,nstokes)
   enddo

!  New, Core run

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL GISSCOXMUNK_CORE &
         ( MAXPARS, NPARS, PARS, XI, SXI, XJ, SXJ, PHI, CPHI, SKPHI, &
           KERNELCORE, KERNELHELP )
   enddo
   call CPU_TIME(e2) ; e3 = (e2-e1)/druns

!  New,  scaling

   call CPU_TIME(e1)
   do irun = 1, nruns
     CALL ADD_FRESNEL_SCALING ( REFIDX, .false., NSSQ, KERNELCORE, KERNELHELP, KERNEL_New )
   enddo
   call CPU_TIME(e2) ; e4 = (e2-e1) ; e_new = e3 + e4
   write(*,'(/A,1p2e17.7,2x,2F12.7)')'New Result, Time/Speedup = ',KERNELCORE, KERNEL_New(1), e_new, e_old/e_new
   do I = 1, NSTOKES
      write(*,'(i2,2x,1p4e17.7)')I,(KERNEL_NEW(4*(I-1)+J),j=1,nstokes)
   enddo

endif

!  done

   stop
end program V2p8p5_VBRDF_exptester_Konly


