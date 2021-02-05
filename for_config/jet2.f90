

module config_module
 implicit none
 real*8, allocatable :: h_forc(:)
end module config_module



subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=2
 
 nx=int(50*fac)
 ny=int(50*fac)

 dx = 10d3/fac
 dy = 10d3/fac
 dt = 100./fac
 csqr = 9.81*5./1000.*1000.
 
 Ro = 1.0 
 Ahbi = 0.1e10/fac**4
 Khbi = 0e9
 !AB_eps = 
 
 snapint =  3600.*24
 runlen =  snapint*2000
 snap_file =  'snap.cdf'
 enable_periodic_x = .true.
 
 !enable_check_energy_conservation = .true.
end subroutine set_parameter




subroutine set_initial_conditions
 use main_module 
 use config_module
 implicit none
 integer :: i,j
 real*8 :: x,y,Ly,Lx
 real*8 :: h_ini(ny), u_ini(ny), f_ini(ny) 
 
 allocate(  h_forc(ny) )
 Ly = dy*ny
 Lx = dx*nx
 
 ! set Coriolis parameter
 do j=1,ny
  y = (j-1)*dy
  f_ini(j) = 1e-4 + 6e-11*(y-Ly/2)
 enddo   

 ! set initial velocity
 do j=1,ny
  y = (j-1)*dy
  u_ini(j) = 1*exp( -( (y-Ly/2) )**2/(Lx*0.02)**2 ) 
 enddo 
 
 ! h from approximate geostrophic balance h_y = - f u 
 h_ini(1) = 0.
 do j=2,ny
  y = (j-1)*dy
  h_ini(j) = h_ini(j-1) - dy*u_ini(j-1)*f_ini(j)
 enddo
 
 ! transfer to model arrays
 do j=js_pe,je_pe
  y = (j-1)*dy
  do i=is_pe,ie_pe
   x = (i-1)*dx
   f(i,j) = f_ini(j)
   u(i,j) = u_ini(j)
   h(i,j) = h_ini(j) + 0.2*sin(X/Lx*10*pi)*cos(Y/Ly*8*pi)
  enddo
 enddo  
 call border_exchg(f)
 call border_exchg(u) 
 call border_exchg(v) 
 call border_exchg(h) 
 
 ! contruct the forcing field
 do j=1,ny
  y = (j-1)*dy
  !h_ini(j) = -10*tanh( (y/Ly-0.5)*20. )
  h_forc(j) = -10*(y/Ly-0.5) 
 enddo 
 
end subroutine set_initial_conditions



subroutine set_forcing
 use main_module 
 use config_module
 implicit none
 integer :: i,j
 
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   h(i,j) = h(i,j) + dt/(3600.*24*10)*(h_forc(j)-h(i,j))
  enddo
 enddo
end subroutine set_forcing