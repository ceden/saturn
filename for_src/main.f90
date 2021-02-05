
program main
 use main_module
 use timing_module   
 implicit none
 integer :: otaum1,enditt,ierr,iargc
 character (len=80) :: arg

 call mpi_init(ierr)
 call my_mpi_init()

 n_pes_i = 1; n_pes_j = 1
 if (n_pes>1) then
   if (iargc() < 2) then
      call halt_stop(' not enough command line input')
   endif
   call getarg(1,arg); read(arg,*) n_pes_i
   call getarg(2,arg); read(arg,*) n_pes_j
 endif
 if (my_pe==0) print'(/a,i4,a,i4,a)','using ',n_pes_i,' x ',n_pes_j ,' PEs'


 call tic('setup')
 call set_parameter()
 call pe_decomposition()
 call allocate_main_module()
 call set_initial_conditions
 if (my_pe==0) then
     print*,' '
     print'(a,i4,a,i4)',' nx x ny = ',nx,' x ',ny
     print'(a,e12.6,a)',' Delta t = ',dt,' s'
     print'(a,e12.6,a)',' Delta x = ',dx,' m'
     print'(a,e12.6,a)',' Delta y = ',dy,' m'
     print'(a,e12.6,a)',' c       = ',sqrt(csqr),' m/s'
     print'(a,e12.6,a)',' Ro      = ',Ro
     print*,' '
 endif
 call toc('setup')

 !--------------------------------------------------------------
 ! forward integration: inititalize output and time integration
 !--------------------------------------------------------------
 call init_snap_cdf
 enditt = itt+int(runlen/dt)

 if (my_pe==0) then
     print*,' '
     print'(a,e8.2,a)',' setup integration for ',runlen,' s'
     print'(a,i10,a,i10)',' from time step ',itt,' to ',enditt
     print'(a,e8.2,a)',' snapshot intervall is ',snapint,' s'
     print'(a,i8,a)',' this is any ',int(snapint/dt),' time steps'
     print*,' '
 endif


 !--------------------------------------------------------------
 ! start time stepping loop
 !--------------------------------------------------------------
 do while (itt < enditt) 
 
       call tic('integrate')
       call set_forcing
       call integrate
       call toc('integrate')
       
       call tic('diagnose')     
       if ( mod(itt,max(1,int(snapint/dt)))  == 0 .or. itt == 0)  call diagnose
       call toc('diagnose')
       
       otaum1=taum1; taum1= tau; tau  = taup1; taup1= otaum1
       itt=itt+1
 enddo
 !--------------------------------------------------------------
 !     show timing results here
 !--------------------------------------------------------------
 !do n = 0,n_pes
     call fortran_barrier
     if (my_pe == 0) then
        print'(/,a,i4)','Timing summary for PE #',my_pe 
        print'(a,f15.2,a)',' costs for measuring      = ',timing_secs('tictoc'),' s'
        print'(a,f15.2,a)',' setup time summary       = ',timing_secs('setup'),' s'
        print'(a,f15.2,a)',' integrate                = ',timing_secs('integrate'),' s'
        print'(a,f15.2,a)',' diagnostics              = ',timing_secs('diagnose') ,' s'
       endif
 !enddo


 if (my_pe==0) print'(/a/)','cancelling MPI service'
 call mpi_finalize(ierr)

end program main



