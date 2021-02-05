


subroutine pe_decomposition
  use main_module   
  implicit none
  integer :: n,tag=0,iloc(20),ierr
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status

! ----------------------------------
!      domain decomposition for each PE
! ----------------------------------
   if (n_pes>1) then
      if (n_pes_i*n_pes_j /= n_pes ) call halt_stop(' n_pes_i times n_pes_j not equal number of PEs')
      i_blk = (nx-1)/n_pes_i + 1    ! i-extent of each block
      j_blk = (ny-1)/n_pes_j + 1    ! j-extent of each block
      my_blk_i = mod(my_pe,n_pes_i)+1! number of PE in i-dir.
      my_blk_j = (my_pe)/n_pes_i + 1 ! number of PE in j-dir.
      is_pe = (my_blk_i-1)*i_blk + 1 ! start index in i-dir of this PE
      ie_pe = min(my_blk_i*i_blk,nx)
      js_pe = (my_blk_j-1)*j_blk + 1
      je_pe = min(my_blk_j*j_blk,ny)
! ----------------------------------
!     last block might have been truncated
! ----------------------------------
      j_blk = je_pe-js_pe+1      
      i_blk = ie_pe-is_pe+1      
! ----------------------------------
! ----------------------------------
!     check for incorrect domain decomposition
! ----------------------------------
      if (my_blk_j==n_pes_j .and. js_pe>je_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in j-direction'
       print*,' choose other number of PEs in j-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_i==n_pes_i .and. is_pe>ie_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in i-direction'
       print*,' choose other number of PEs in i-direction'
       call halt_stop(' in pe_decomposition')
      endif
   else
       n_pes_j = 1; n_pes_i = 1
       i_blk = nx; j_blk = ny
       my_blk_j = 1 ; my_blk_i = 1 
       js_pe = 1; je_pe = ny
       is_pe = 1; ie_pe = nx
   endif
! ----------------------------------
!      print out the PE decomposition, let PE 0 talk
! ----------------------------------
   do n=0,n_pes-1
     if (n==0) then
       iloc(1:6) = (/my_blk_i,my_blk_j,is_pe,ie_pe,js_pe,je_pe/)
     else
       if (my_pe==n) then
          iloc(1:6) = (/my_blk_i,my_blk_j,is_pe,ie_pe,js_pe,je_pe/)
          call mpi_send(iloc,6,mpi_integer,0,tag,my_comm,ierr)
       endif
       if (my_pe==0) call mpi_recv(iloc,6,mpi_integer,n,tag,my_comm,status,ierr)
     endif
     if (my_pe==0) print'(a,i4,a,i4,a,i4,a,i4,a,i4)','domain of PE #',n,' i=',iloc(3),':',iloc(4),'   j=',iloc(5),':',iloc(6)
   enddo
   if (my_pe==0) print*,' '
   call fortran_barrier
end subroutine pe_decomposition




 subroutine my_mpi_init()
!--------------------------------------------------------------
!     intitialize mpi system for model
!--------------------------------------------------------------
  use main_module   
  implicit none
  integer :: nlen,ierr
  include "mpif.h"
      character (len=MPI_MAX_PROCESSOR_NAME) :: pname
      !if (comm_ == MPI_COMM_NULL) then
      !  print *, 'You passed MPI_COMM_NULL !!!'
      !  return
      ! end if
      ! my_comm=comm_
       my_comm = MPI_COMM_WORLD
       call MPI_Comm_rank(my_comm, my_pe, ierr)
       call MPI_Comm_size(my_comm, n_pes, ierr)
       call MPI_Get_processor_name(pname, nlen, ierr)
       call my_mpi_test(my_comm)
 end subroutine my_mpi_init


subroutine halt_stop(string)
!--------------------------------------------------------------
!     controlled stop, should not be called from python
!--------------------------------------------------------------
      implicit none
      character*(*) :: string
      integer :: ierr,code,my_pe
      include "mpif.h"
      call mpi_comm_rank(MPI_COMM_WORLD,my_pe,ierr)
      print*,' global pe #',my_pe,' : ',string
      print*,' global pe #',my_pe,' aborting '
      code=99
      call MPI_ABORT(mpi_comm_world, code, IERR)
end subroutine halt_stop



subroutine fortran_barrier
!--------------------------------------------------------------
!     A barrier for the local sub domain
!     for use in fortran part only
!--------------------------------------------------------------
  use main_module   
  implicit none
  integer :: ierr
  call mpi_barrier(my_comm, ierr)
end subroutine fortran_barrier



subroutine my_mpi_test(my_comm)
!--------------------------------------------------------------
!     test some basic mpi routines
!--------------------------------------------------------------
      implicit none
      integer :: my_comm
      integer :: my_pe=-1,all_pes,xint,xint2,ierr
      real*8    :: xreal,xreal2
      include "mpif.h"
!   get some mpi infos first
      call mpi_comm_rank(my_comm       ,my_pe,ierr)
      if (my_pe==0) print*,' testing mpi routines'
      call mpi_comm_size(my_comm       ,all_pes,ierr)
!   try first global barrier
      call mpi_barrier(my_comm       , ierr)
!   try broadcasting
      xreal = 1.0
      call mpi_bcast(xreal,1,mpi_real8,0,my_comm       ,ierr)
      xint = 1
      call mpi_bcast(xint,1,mpi_integer,0,my_comm       ,ierr)
!   check results of broadcasting
      if (xreal /= 1.0 ) then
       print*,'fatal: MPI test failed on broadcasting reals for PE #',my_pe
       stop
      endif
      if (xint /= 1 ) then
       print*,'fatal: MPI test failed on broadcasting integer for PE #',my_pe
       stop
      endif
      call mpi_barrier(my_comm       , ierr)
!   try global sum
      xreal = 2.0
      call mpi_allreduce(xreal,xreal2,1,mpi_real8,MPI_SUM,my_comm       ,ierr)
      xint = 2
      call mpi_allreduce(xint,xint2,1,mpi_integer,MPI_SUM,my_comm       ,ierr)
!   check results 
      xreal = xreal2/all_pes
      if (xreal /= 2.0 ) then
       print*,'fatal: MPI test failed on global sum (real) for PE #',my_pe
       stop
      endif
      xint = xint2/all_pes
      if (xint /= 2.0 ) then
       print*,'fatal: MPI test failed on global sum (int) for PE #',my_pe
       stop
      endif
      call mpi_barrier(my_comm       , ierr)
end subroutine my_mpi_test


subroutine pe0_bcast_int(a,len)
!--------------------------------------------------------------
!     Broadcast an integer vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      integer, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_integer,0,my_comm,ierr)
end subroutine pe0_bcast_int


subroutine pe0_bcast(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      real*8, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_real8,0,my_comm,ierr)
end subroutine pe0_bcast



subroutine bcast_real(x,len,pe)
!--------------------------------------------------------------
!     Broadcast a real vector from PE pe to others
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer :: len,ierr,pe
      real*8 :: x(len)
      include "mpif.h"
      call mpi_barrier(my_comm, ierr)
      call mpi_bcast(x,len,mpi_real8,pe,my_comm,ierr)
end subroutine bcast_real

subroutine bcast_integer(x,len,pe)
!--------------------------------------------------------------
!     Broadcast an integer vector from PE pe to others
!--------------------------------------------------------------
      use main_module
      implicit none
      integer :: len,ierr,pe
      integer :: x(len)
      include "mpif.h"
      call mpi_barrier(my_comm, ierr)
      call mpi_bcast(x,len,mpi_integer,pe,my_comm,ierr)
end subroutine bcast_integer



subroutine global_max(x)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MAX,my_comm       ,ierr)
      x = x_sym2
 end subroutine global_max
 



subroutine global_min(x)
!--------------------------------------------------------------
!     Get the min of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MIN,my_comm       ,ierr)
      x = x_sym2
end subroutine global_min


subroutine global_sum(x)
!--------------------------------------------------------------
!     Do a sum of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_SUM,my_comm       ,ierr)
      x = x_sym2
end subroutine global_sum






subroutine global_max_int(x)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MAX,my_comm       ,ierr)
      x = x_sym2
 end subroutine global_max_int
 

subroutine global_min_int(x)
!--------------------------------------------------------------
!     Get the min of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MIN,my_comm       ,ierr)
      x = x_sym2
end subroutine global_min_int


subroutine global_sum_int(x)
!--------------------------------------------------------------
!     Do a sum of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_SUM,my_comm       ,ierr)
      x = x_sym2
end subroutine global_sum_int


 
 subroutine global_max_int2(x,len)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x(len)
      integer,intent(in) :: len
      integer    :: x_sym(len),x_sym2(len),ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,len,mpi_integer,MPI_MAX,my_comm,ierr)
      x = x_sym2
 end subroutine global_max_int2




 subroutine global_max2(x,len)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x(len)
      integer,intent(in) :: len
      real*8    :: x_sym(len),x_sym2(len)
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,len,mpi_real8,MPI_MAX,my_comm ,ierr)
      x = x_sym2
 end subroutine global_max2




subroutine border_exchg(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
! and also for periodic boundaries
!--------------------------------------------------------------
  use main_module   
  implicit none
  real*8, intent(inout)  :: a(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  integer  ::  tag=0, ierr,i,j,len,north,south,east,west
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status

  south = my_pe-n_pes_i 
  north = my_pe+n_pes_i
  if (my_blk_j==1)       south = my_pe+ n_pes_i*(n_pes_j-1)
  if (my_blk_j==n_pes_j) north = my_pe- n_pes_i*(n_pes_j-1)

  if ( n_pes_j > 1) then
     len=(ie_pe-is_pe+1+2*onx)
     do j=1,onx
       if (my_blk_j /=1 )       call mpi_send(a(:,js_pe+j-1),len,mpi_real8,south,tag,my_comm,ierr)
       if (my_blk_j /= n_pes_j) call mpi_recv(a(:,je_pe+j  ),len,mpi_real8,north,tag,my_comm,status,ierr)
       if (enable_periodic_y) then
        if (my_blk_j ==1 )       call mpi_send(a(:,js_pe+j-1),len,mpi_real8,south,tag,my_comm,ierr)
        if (my_blk_j == n_pes_j) call mpi_recv(a(:,je_pe+j  ),len,mpi_real8,north,tag,my_comm,status,ierr)
       endif 
     enddo
     do j=1,onx
       if (my_blk_j /= n_pes_j) call mpi_send(a(:,je_pe-j+1),len,mpi_real8,north,tag,my_comm,ierr)
       if (my_blk_j /=1 )       call mpi_recv(a(:,js_pe-j  ),len,mpi_real8,south,tag,my_comm,status,ierr)
       if (enable_periodic_y) then
        if (my_blk_j == n_pes_j) call mpi_send(a(:,je_pe-j+1),len,mpi_real8,north,tag,my_comm,ierr)
        if (my_blk_j ==1 )       call mpi_recv(a(:,js_pe-j  ),len,mpi_real8,south,tag,my_comm,status,ierr)
       endif 
     enddo
  else
   if (enable_periodic_y) then
     do j=1,onx
      a(:,ny+j) = a(:,j  )
      a(:,1-j ) = a(:,ny-j+1) 
     enddo
   endif  
  endif

  west = my_pe-1
  if (my_blk_i==1      ) west = my_pe+ n_pes_i-1
  east = my_pe+1
  if (my_blk_i==n_pes_i) east = my_pe- (n_pes_i-1)

  if ( n_pes_i > 1) then
     len=(je_pe-js_pe+1+2*onx)
     do i=1,onx
       if (my_blk_i /=1 )       call mpi_send(a(is_pe+i-1,:),len,mpi_real8,west,tag,my_comm,ierr)
       if (my_blk_i /= n_pes_i) call mpi_recv(a(ie_pe+i  ,:),len,mpi_real8,east,tag,my_comm,status,ierr)
       if (enable_periodic_x) then
        if (my_blk_i ==1 )       call mpi_send(a(is_pe+i-1,:),len,mpi_real8,west,tag,my_comm,ierr)
        if (my_blk_i == n_pes_i) call mpi_recv(a(ie_pe+i  ,:),len,mpi_real8,east,tag,my_comm,status,ierr)
       endif 
     enddo
     do i=1,onx
       if (my_blk_i /= n_pes_i) call mpi_send(a(ie_pe-i+1,:),len,mpi_real8,east,tag,my_comm,ierr)
       if (my_blk_i /=1 )       call mpi_recv(a(is_pe-i  ,:),len,mpi_real8,west,tag,my_comm,status,ierr)
       if (enable_periodic_x) then
        if (my_blk_i == n_pes_i) call mpi_send(a(ie_pe-i+1,:),len,mpi_real8,east,tag,my_comm,ierr)
        if (my_blk_i ==1 )       call mpi_recv(a(is_pe-i  ,:),len,mpi_real8,west,tag,my_comm,status,ierr)
       endif 
     enddo
  else
   if (enable_periodic_x) then
    do i=1,onx
     a(nx+i,:) = a(i     ,:)
     a(1-i ,:) = a(nx-i+1,:) 
    enddo
   endif 
  endif

end subroutine border_exchg









subroutine pe0_recv_2D(a)
!--------------------------------------------------------------
!     all PEs send their data of a 2D array to PE0
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8, intent(inout) :: a(nx,ny)
      integer                     :: js,je,iproc,is,ie
      integer                     :: tag=0, ierr,len
      include "mpif.h"
      integer, dimension(MPI_STATUS_SIZE) :: Status
      js=js_pe; je=je_pe; is=is_pe; ie=ie_pe
      do iproc=1,n_pes-1
       call mpi_barrier(my_comm       ,ierr)
       if ( my_pe == iproc ) then
        call mpi_send(js,1,mpi_integer,0,tag,my_comm,ierr)
        call mpi_send(je,1,mpi_integer,0,tag,my_comm,ierr)
        call mpi_send(is,1,mpi_integer,0,tag,my_comm,ierr)
        call mpi_send(ie,1,mpi_integer,0,tag,my_comm,ierr)
        len=(ie-is+1)*(je-js+1)
        call mpi_send(a(is:ie,js:je),len,mpi_real8,0,tag,my_comm,ierr)
       endif
       if ( my_pe == 0 ) then
        call mpi_recv(js,1,mpi_integer,iproc,tag,my_comm,Status,ierr)
        call mpi_recv(je,1,mpi_integer,iproc,tag,my_comm,Status,ierr)
        call mpi_recv(is,1,mpi_integer,iproc,tag,my_comm,Status,ierr)
        call mpi_recv(ie,1,mpi_integer,iproc,tag,my_comm,Status,ierr)
        len=(ie-is+1)*(je-js+1)
        call mpi_recv(a(is:ie,js:je),len,mpi_real8,iproc,tag,my_comm,Status,ierr)
       endif
       call mpi_barrier(my_comm       ,ierr)
      enddo
end subroutine pe0_recv_2D


