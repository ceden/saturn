





subroutine init_snap_cdf
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,i,j
 integer :: lon_tdim,lat_tdim,itimedim
 integer :: lon_tid,lat_tid,itimeid,id
 real*8 :: x(nx),y(ny)

 if (my_pe==0) then
  !ncid = nccre ('snap.cdf', NCCLOB, iret)
  iret = NF_CREATE (snap_file,OR(NF_SHARE,NF_64BIT_OFFSET),ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  lon_tdim  = ncddef(ncid, 'xt', nx, iret)
  Lat_tdim  = ncddef(ncid, 'yt', ny, iret)
  iTimedim  = ncddef(ncid, 'time', nf_unlimited, iret)
  Lon_tid  = ncvdef (ncid,'xt',nf_double,1,lon_tdim,iret)
  Lat_tid  = ncvdef (ncid,'yt',nf_double,1,lat_tdim,iret)
  itimeid  = ncvdef (ncid,'time', nf_double,1,itimedim,iret)
  
  iret = nf_put_att_text (ncid,lon_tid,'long_name',1,'x')
  iret = nf_put_att_text (ncid,lat_tid,'long_name',1,'y')
  iret = nf_put_att_text (ncid,itimeid,'long_name',4,'time')
  iret = nf_put_att_int (ncid,itimeid,'time_origin',nf_int,1,0)
  
  id  = ncvdef (ncid,'h', nf_double ,3,(/Lon_tdim, lat_tdim,iTimedim/),iret)
  iret = nf_put_att_text (ncid,id,'long_name',9,'thickness')
  id  = ncvdef (ncid,'u', nf_double ,3,(/Lon_tdim, lat_tdim,iTimedim/),iret)
  iret = nf_put_att_text (ncid,id,'long_name',8,'velocity')
  id  = ncvdef (ncid,'v', nf_double ,3,(/Lon_tdim, lat_tdim,iTimedim/),iret)
  iret = nf_put_att_text (ncid,id,'long_name',8,'velocity')

  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dx',nf_double,1,dx)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dy',nf_double,1,dy)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'dt',nf_double,1,dt)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Ahbi',nf_double,1,Ahbi)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Khbi',nf_double,1,Khbi)
  !iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'f0',nf_double,1,f0)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'csqr',nf_double,1,csqr)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'snapint',nf_double,1,snapint)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'runlen',nf_double,1,runlen)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'Ro',nf_double,1,Ro)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'AB_eps',nf_double,1,AB_eps)
 
  call ncendf(ncid, iret)
  do i=1,nx
   x(i)=(i-0.5)*dx
  enddo
  do j=1,ny
   y(j)=(j-0.5)*dy
  enddo
  iret= nf_put_vara_double(ncid,lon_Tid,1,nx,x)
  iret= nf_put_vara_double(ncid,lat_tid,1,ny,y)
  call ncclos (ncid, iret)
 endif
 call fortran_barrier()

end subroutine init_snap_cdf


subroutine diagnose
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret, itdimid,ilen,id,itimeid,i,j
 real*8 :: fxa,umax,hmax
 real*8 :: aloc(nx,ny)

 fxa = itt*dt
 if (my_pe==0) print'(a,i6,a,e12.6,a)','writing snapshot at itt=',itt,', ',fxa,' s'

 umax=0;hmax=0;
 do j=js_pe,je_pe
   do i=is_pe,ie_pe
     umax=max(umax,sqrt(u(i,j)**2+v(i,j)**2))
     hmax=max(hmax,abs(h(i,j)))
   enddo
 enddo
 call global_max(umax)
 call global_max(hmax)

 if (my_pe==0) then
   print*,'max CFL(u)  = ',umax*dt/dx
   print*,'max CFL(h)  = ',hmax*dt/dx
   print*,'max PECL(u) = ',umax*dx**3/(1e-12+Ahbi)
   print*,'max PECL(h) = ',hmax*dx**3/(1e-12+Khbi)
 endif



 if (my_pe==0) then
   iret=nf_open(snap_file,NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   iret=nf_inq_varid(ncid,'time',itimeid)
   ilen=ilen+1
   print*,'time steps in file : ',ilen
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,fxa)
 endif

 aloc(is_pe:ie_pe,js_pe:je_pe) = u(is_pe:ie_pe,js_pe:je_pe)
 call pe0_recv_2D(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'u',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),aloc)
 endif

 aloc(is_pe:ie_pe,js_pe:je_pe) = v(is_pe:ie_pe,js_pe:je_pe)
 call pe0_recv_2D(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'v',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),aloc)
 endif

 aloc(is_pe:ie_pe,js_pe:je_pe) = h(is_pe:ie_pe,js_pe:je_pe)
 call pe0_recv_2D(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'h',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,ilen/),(/nx,ny,1/),aloc)
 endif

 if (my_pe==0) then
   call ncclos (ncid, iret)
 endif


end subroutine diagnose




