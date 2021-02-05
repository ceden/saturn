


 
subroutine integrate
 !---------------------------------------------------------------------------------
 ! Sadourny (1975) energy conserving scheme 
 !---------------------------------------------------------------------------------
 use main_module
 use timing_module   
 implicit none
 integer :: i,j
 real*8 :: ke,ke_abs,fxa

 !---------------------------------------------------------------------------------
 ! total thickness 
 !--------------------------------------------------------------------------------- 
 hc = Ro*h + csqr
 
 !---------------------------------------------------------------------------------
 ! thickness fluxes
 !--------------------------------------------------------------------------------- 
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx-1
   fe(i,j) =  u(i,j)*0.5*(hc(i,j)+hc(i+1,j))
   fn(i,j) =  v(i,j)*0.5*(hc(i,j)+hc(i,j+1))
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! thickness tendency
 !--------------------------------------------------------------------------------- 
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
    dh(i,j,tau)= - ( (fe(i,j)-fe(i-1,j))/dx +(fn(i,j)-fn(i,j-1))/dy )
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! planetary and relative vorticity
 !--------------------------------------------------------------------------------- 
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx-1
    q(i,j) = f(i,j) + Ro*(  (v(i+1,j)-v(i,j))/dx - (u(i,j+1)-u(i,j))/dx )
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! potential vorticity
 !--------------------------------------------------------------------------------- 
 do j=js_pe-onx,je_pe+onx-1
  do i=is_pe-onx,ie_pe+onx-1
    q(i,j) = q(i,j)/( 0.25*(hc(i,j)+hc(i+1,j)+hc(i,j+1)+hc(i+1,j+1)) )
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! momentum tendency by -qV, qU 
 !--------------------------------------------------------------------------------- 
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,tau) = +0.5*(   q(i,j)*0.5*(fn(i,j)+fn(i+1,j)) + q(i,j-1)*0.5*(fn(i,j-1)+fn(i+1,j-1)) )
   dv(i,j,tau) = -0.5*(   q(i,j)*0.5*(fe(i,j)+fe(i,j+1)) + q(i-1,j)*0.5*(fe(i-1,j)+fe(i-1,j+1)) )
  enddo
 enddo
 
 !---------------------------------------------------------------------------------
 ! kinetic energy
 !--------------------------------------------------------------------------------- 
 do j=js_pe-onx+1,je_pe+onx
   do i=is_pe-onx+1,ie_pe+onx
    K(i,j) = 0.5*(  0.5*(u(i,j)**2 + u(i-1,j)**2) + 0.5*(v(i,j)**2 + v(i,j-1)**2) )
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! momentum tendency by kinetic energy gradient
 !--------------------------------------------------------------------------------- 
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
    du(i,j,tau) = du(i,j,tau) - Ro*(K(i+1,j)-K(i,j))/dx
    dv(i,j,tau) = dv(i,j,tau) - Ro*(K(i,j+1)-K(i,j))/dy
   enddo
 enddo

 !---------------------------------------------------------------------------------
 ! momentum tendency by thickness gradient
 !---------------------------------------------------------------------------------
 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   du(i,j,tau) = du(i,j,tau) - (h(i+1,j)-h(i,j))/dx
   dv(i,j,tau) = dv(i,j,tau) - (h(i,j+1)-h(i,j))/dy
  enddo
 enddo

 !---------------------------------------------------------------------------------
 ! zero out momentum time tendencies in case of closed boundary conditions
 !---------------------------------------------------------------------------------
 if (my_blk_i == n_pes_i .and. .not.enable_periodic_x) du(nx,:,tau) = 0.
 if (my_blk_j == n_pes_j .and. .not.enable_periodic_y) dv(:,ny,tau) = 0.

 if (Ahbi>0d0) call biharmonic(u,Ahbi,du_mix,'u')
 if (Ahbi>0d0) call biharmonic(v,Ahbi,dv_mix,'v')
 if (Khbi>0d0) call biharmonic(h,Khbi,dh_mix,'h')
 
 !---------------------------------------------------------------------------------
 ! check total energy conservation, only for testing, only without friction
 !---------------------------------------------------------------------------------
 if (enable_check_energy_conservation) then
  ke=0;ke_abs=0
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     fxa = Ro*( fe(i,j)*du(i,j,tau) + fn(i,j)*dv(i,j,tau) )
     ke     = ke     +      dx*dy*fxa
     ke_abs = ke_abs + abs( dx*dy*fxa )
   enddo
  enddo
  do j=js_pe,je_pe
   do i=is_pe,ie_pe
     fxa = (Ro*h(i,j)+csqr+Ro**2*K(i,j))*dh(i,j,tau)
     ke     = ke     +      dx*dy*fxa 
     ke_abs = ke_abs + abs( dx*dy*fxa )
   enddo
  enddo
  call global_sum(ke)
  call global_sum(ke_abs)
  if (my_pe==0) print*,'tot.En.change = ',ke,' relative error = ',ke/(1e-22+ke_abs)
 endif

 !---------------------------------------------------------------------------------
 ! time stepping with Adams Bashforth 
 !---------------------------------------------------------------------------------
 if (itt==0) then
  h = h + dt*dh(:,:,tau) + dt*dh_mix
  u = u + dt*du(:,:,tau) + dt*du_mix
  v = v + dt*dv(:,:,tau) + dt*dv_mix
 else 
  h = h + dt*( (1.5+AB_eps)*dh(:,:,tau) - ( 0.5+AB_eps)*dh(:,:,taum1)) + dt*dh_mix
  u = u + dt*( (1.5+AB_eps)*du(:,:,tau) - ( 0.5+AB_eps)*du(:,:,taum1)) + dt*du_mix 
  v = v + dt*( (1.5+AB_eps)*dv(:,:,tau) - ( 0.5+AB_eps)*dv(:,:,taum1)) + dt*dv_mix
 endif
 
 call border_exchg(u) 
 call border_exchg(v) 
 call border_exchg(h) 

end subroutine integrate













subroutine biharmonic(p1,Hmix,dp,typ)
!---------------------------------------------------------------------------------
! biharmonic mixing/friction
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,j
 real*8 :: p1(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx),Hmix
 real*8 :: dp(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 real*8 :: del2(is_pe-onx:ie_pe+onx,js_pe-onx:je_pe+onx)
 character*(*) :: typ
 
 do j=js_pe-1,je_pe+1
  do i=is_pe-2,ie_pe+1
   fe(i,j) = -Hmix*(p1(i+1,j)-p1(i,j))/dx
  enddo
 enddo

 if (my_blk_i == 1      .and. .not.enable_periodic_x)               fe(1-onx:0,:)   = 0d0
 if (my_blk_i == n_pes_i.and. .not.enable_periodic_x)               fe(nx:nx+onx,:) = 0d0
 if (my_blk_i == n_pes_i.and. .not.enable_periodic_x.and. typ=='u') fe(nx-1,:) = 0d0
 
 
 do j=js_pe-2,je_pe+1
  do i=is_pe-1,ie_pe+1
   fn(i,j) = -Hmix*(p1(i,j+1)-p1(i,j))/dy
  enddo
 enddo

 if (my_blk_j == 1       .and. .not. enable_periodic_y)             fn(:,1-onx:0)   = 0d0
 if (my_blk_j == n_pes_j .and. .not. enable_periodic_y)             fn(:,ny:ny+onx) = 0d0
 if (my_blk_i == n_pes_i.and. .not.enable_periodic_x.and. typ=='v') fn(:,ny-1) = 0d0
 
 do j=js_pe-1,je_pe+1
  do i=is_pe-1,ie_pe+1
   del2(i,j)= (fe(i,j)-fe(i-1,j))/dx +(fn(i,j)-fn(i,j-1))/dy 
  enddo
 enddo

 do j=js_pe,je_pe
  do i=is_pe-1,ie_pe
   fe(i,j) = (del2(i+1,j)-del2(i,j))/dx
  enddo
 enddo

 if (my_blk_i == 1       .and. .not. enable_periodic_x)             fe(1-onx:0,:)   = 0d0
 if (my_blk_i == n_pes_i .and. .not. enable_periodic_x)             fe(nx:nx+onx,:) = 0d0
 if (my_blk_i == n_pes_i.and. .not.enable_periodic_x.and. typ=='u') fe(nx-1,:) = 0d0
 
 do j=js_pe-1,je_pe
  do i=is_pe,ie_pe
   fn(i,j) = (del2(i,j+1)-del2(i,j))/dy
  enddo
 enddo

 if (my_blk_j == 1       .and. .not. enable_periodic_y)             fn(:,1-onx:0)   = 0d0
 if (my_blk_j == n_pes_j .and. .not. enable_periodic_y)             fn(:,ny:ny+onx) = 0d0
 if (my_blk_i == n_pes_i.and. .not.enable_periodic_x.and. typ=='v') fn(:,ny-1) = 0d0

 do j=js_pe,je_pe
  do i=is_pe,ie_pe
   dp(i,j)= (fe(i,j)-fe(i-1,j))/dx +(fn(i,j)-fn(i,j-1))/dy 
  enddo
 enddo
end subroutine biharmonic
