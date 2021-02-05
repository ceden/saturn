

subroutine pe_decomposition
  use main_module   
  implicit none
! ----------------------------------
!      domain decomposition for each PE
! ----------------------------------
     n_pes = 1
     n_pes_j = 1
     n_pes_i = 1
     j_blk = ny
     i_blk = nx
     my_blk_j = 1 
     my_blk_i = 1 
     js_pe = 1
     je_pe = ny
     is_pe = 1
     ie_pe = nx
     my_comm = 0
     my_pe = 0
end subroutine pe_decomposition


subroutine mpi_finalize(ierr)
  use main_module   
  implicit none
  integer, intent(out) :: ierr
  ierr=0
end subroutine mpi_finalize

subroutine mpi_init(ierr)
  use main_module   
  implicit none
  integer, intent(out) :: ierr
  ierr=0
end subroutine mpi_init


 subroutine my_mpi_init(comm_)
!--------------------------------------------------------------
!     intitialize mpi system for model
!--------------------------------------------------------------
  use main_module   
  implicit none
  integer :: comm_
  comm_=0
 end subroutine my_mpi_init



subroutine halt_stop(string)
!--------------------------------------------------------------
!     controlled stop, should not be called from python
!--------------------------------------------------------------
      implicit none
      character*(*) :: string
      print*,string
      stop
end subroutine halt_stop


subroutine fortran_barrier
end subroutine fortran_barrier


subroutine pe0_bcast_int(a,len)
!--------------------------------------------------------------
!     Broadcast an integer vector from pe0 to all other pe
!--------------------------------------------------------------
      implicit none
      integer, intent(in) :: len
      integer, intent(inout) :: a(len)
end subroutine pe0_bcast_int


subroutine pe0_bcast(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      implicit none
      integer, intent(in) :: len
      real*8, intent(inout) :: a(len)
end subroutine pe0_bcast



subroutine bcast_real(x,len,pe)
!--------------------------------------------------------------
!     Broadcast a real vector from PE pe to others
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer :: len,ierr,pe
      real*8 :: x(len)
end subroutine bcast_real

subroutine bcast_integer(x,len,pe)
!--------------------------------------------------------------
!     Broadcast an integer vector from PE pe to others
!--------------------------------------------------------------
      use main_module
      implicit none
      integer :: len,ierr,pe
      integer :: x(len)
end subroutine bcast_integer


subroutine global_max(x)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      real*8, intent(inout)    :: x
end subroutine global_max


subroutine global_min(x)
!--------------------------------------------------------------
!     Get the min of real x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      real*8, intent(inout)    :: x
end subroutine global_min


subroutine global_sum(x)
!--------------------------------------------------------------
!     Do a sum of real x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      real*8, intent(inout)    :: x
end subroutine global_sum


subroutine global_max_int(x)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      integer,intent(inout)    :: x
 end subroutine global_max_int


subroutine global_min_int(x)
!--------------------------------------------------------------
!     Get the min of integer x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      integer,intent(inout)    :: x
end subroutine global_min_int


subroutine global_sum_int(x)
!--------------------------------------------------------------
!     Do a sum of integer x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      integer,intent(inout)    :: x
end subroutine global_sum_int


 
 subroutine global_max_int2(x,len)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!-------------------------------------------------------------
      implicit none
      integer,intent(inout)    :: x(len)
      integer,intent(in) :: len
 end subroutine global_max_int2




 subroutine global_max2(x,len)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      implicit none
      real*8,intent(inout)    :: x(len)
      integer,intent(in) :: len
 end subroutine global_max2






subroutine border_exchg(a)
!--------------------------------------------------------------
!     Exchange overlapping areas of 3D array a in all PEs of sub 
!     domain. Number of overlapping indicees are given by jx.
!--------------------------------------------------------------
  use main_module   
  implicit none
  real*8, intent(inout)  :: a
end subroutine border_exchg



subroutine setcyclic(p1)
!--------------------------------------------------------------
!       set cyclic boundary conditions for 2D array
!--------------------------------------------------------------
  use main_module   
  implicit none
  real*8, intent(inout) :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
  integer :: i,j

  do i=1,onx
     p1(nx+i,:) = p1(i  ,:)
     p1(1-i,:)  = p1(nx-i+1,:) 
  enddo

  do j=1,onx
     p1(:,ny+j) = p1(:,j  )
     p1(:,1-j ) = p1(:,ny-j+1) 
  enddo
end subroutine setcyclic



subroutine pe0_recv_2D(a)
!--------------------------------------------------------------
!     all PEs send their data of a 2D array to PE0
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8, intent(inout)  :: a(nx,ny)
end subroutine pe0_recv_2D

