


subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=2
 
 nx = int(100*fac)
 ny = int(50*fac)

 dx = 2d3/fac
 dy = 2d3/fac
 dt = 50./fac
 
 csqr = 9.81*5./1000.*1000.  ! g delta rho/rho_0 *H
 Ro = 1.0 
 
 Ahbi = 4e6/fac**4
 Khbi = 4e6/fac**4
 AB_eps = 0.1
 
 snapint =  3600.
 runlen =  snapint*2000
 snap_file =  'snap.cdf'
 
 enable_periodic_x = .true.
 
end subroutine set_parameter




subroutine set_initial_conditions
 use main_module 
 implicit none
 integer :: i,j
 real*8 :: x,y,Ly,Lx,u_ini(ny),h_ini(ny),f_ini(ny)
 
 Ly = dy*ny
 Lx = dx*nx
 
 ! set Coriolis parameter
 do j=1,ny
  f_ini(j) = 1e-4
 enddo
 
 ! set initial velocity
 do j=1,ny
  y = (j-1)*dy
  u_ini(j) = 1.*exp( -( (y-Ly/2) )**2/(Lx*0.02)**2 ) 
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
   h(i,j) = h_ini(j) + 0.1*sin(X/Lx*10*pi)*cos(Y/Ly*8*pi)
  enddo
 enddo  
 ! exchange over domain borders and PEs
 call border_exchg(f)
 call border_exchg(u) 
 call border_exchg(v) 
 call border_exchg(h) 

end subroutine set_initial_conditions



subroutine set_forcing
end subroutine set_forcing


