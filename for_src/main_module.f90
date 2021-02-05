

module main_module
 implicit none

!---------------------------------------------------------------------------------
! fixed parameters
!---------------------------------------------------------------------------------
  real*8, parameter :: version = 0.1
  real*8, parameter :: pi      = 3.14159265358979323846264338327950588
!---------------------------------------------------------------------------------
! model parameter
!---------------------------------------------------------------------------------
  integer :: nx         ! number of grid points in zonal direction       
  integer :: ny         ! number of grid points in meridional direction       
  real*8 :: Ahbi = 0.d0 ! biharmonic viscosity in m^4/s
  real*8 :: Khbi = 0.d0 ! biharmonic diffusivity in m^4/s
  real*8 :: csqr = 0.d0 ! 
  real*8 :: Ro   = 1.d0 ! Rossby number
  real*8 :: dx          ! zonal extent of grid cell in m
  real*8 :: dy          ! meridional extent of grid cell in m
  real*8 :: dt          ! time step
!---------------------------------------------------------------------------------
! switches for configuration
!---------------------------------------------------------------------------------
  logical :: enable_periodic_x = .false.                ! set periodic in x direction 
  logical :: enable_periodic_y = .false.                ! set periodic in y direction 
  logical :: enable_check_energy_conservation = .false. ! only for testing 
!---------------------------------------------------------------------------------
! model arrays
!---------------------------------------------------------------------------------
  real*8, allocatable :: h(:,:), u(:,:), v(:,:)                 ! thickness, zonal and meridiona velocity
  real*8, allocatable :: dh(:,:,:), du(:,:,:), dv(:,:,:)        !   time dendencies
  real*8, allocatable :: dh_mix(:,:), du_mix(:,:), dv_mix(:,:)  ! time tendencies due to mixing terms
  real*8, allocatable :: f(:,:)          ! Coriolis parameter
  real*8, allocatable :: fe(:,:),fn(:,:) ! east/northward fluxes
  real*8, allocatable :: q(:,:),K(:,:)   ! potential vorticity, kinetic energy
  real*8, allocatable :: hc(:,:)         ! total thickness
!---------------------------------------------------------------------------------
! time integration stuff
!---------------------------------------------------------------------------------
  integer :: taum1 = 1, tau = 2, taup1 = 3 ! last, current and furture time step index
  integer :: itt = 0                       ! current time step number
  real*8 :: runlen = 10*360.               ! length of integration in s
  real*8 :: snapint = 0.                   ! snapshot interval in s
  real*8 :: AB_eps = 0.01
  character*80 :: snap_file = 'snap.cdf'
!---------------------------------------------------------------------------------
! Parallel domain setup
!---------------------------------------------------------------------------------
  integer :: n_pes     ! total number of processors
  integer :: my_pe     ! index of this processor from 0 to n_pes-1
  integer :: n_pes_i   ! total number of processors in x direction
  integer :: n_pes_j   ! total number of processors in y direction
  integer :: my_blk_i  ! index of this processor in x direction from 1 to n_pes_i
  integer :: my_blk_j  ! index of this processor in y direction from 1 to n_pes_j
  integer :: i_blk     ! grid points of domain decompostion in x direction 
  integer :: j_blk     ! grid points of domain decompostion in y direction
  integer :: is_pe     ! start index of grid points in x direction of this processor
  integer :: ie_pe     ! end index of grid points in x direction of this processor
  integer :: js_pe     ! start index of grid points in y direction of this processor
  integer :: je_pe     ! end index of grid points in y direction of this processor
  integer :: onx=2     ! number of overlapping points in x and y direction
  integer :: my_comm=0 ! communicator for MPI library
end module main_module



subroutine  allocate_main_module
 use main_module
 implicit none
 allocate(   h(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate(   u(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate(   v(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 h=0; u=0; v=0
 allocate(   dh(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) )
 allocate(   du(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) )
 allocate(   dv(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx,3) )
 dh=0; du=0; dv=0
 allocate(   dh_mix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate(   du_mix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate(   dv_mix(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 dh_mix=0; du_mix=0; dv_mix=0
 allocate(   f(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 f=0.
 allocate( fe(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate( fn(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate( q(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate( K(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 allocate( hc(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx) )
 
end subroutine  allocate_main_module
