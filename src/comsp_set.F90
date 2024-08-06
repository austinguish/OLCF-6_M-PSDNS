subroutine comsp_set
 
   use comp
 
   implicit none
 
   real(b8) beta_min

   beta_min=min(beta1,beta2,beta3)
   mxyz=max(nx*beta1,ny*beta2,nz*beta3)/beta_min
   if(taskid.eq.0) then 
      write (6,*) 'comsp_set: nxyz=',nx,ny,nz
      write (6,*) 'comsp_set: mxyz=',mxyz
      write (6,*) 'comsp_set: beta=',beta1,beta2,beta3,beta_min
      call custom_flush(6)
   end if

   allocate (kx2(nxh),ky2(ny),kz2(nz))
   allocate (tfact(nxh))
   tfact(:)=2.
   tfact(1)=1.
 
   allocate (kt(ncp))
 
   ! energy and dissipation sectrum
   allocate(eijk(mxyz,ncp))
   allocate(dijk(mxyz,ncp))
   allocate(sijk(mxyz,ncp))
   allocate(grijk(mxyz,ncp,3))
   allocate(ek(mxyz))
   allocate(dk(mxyz))
   allocate(sk(mxyz))
   allocate(corr(ncp))
   allocate(tmij(3+nc,3))
   allocate(lijk(ncp,3+nc))
   allocate(laak(3+nc,3))
   allocate(taylor(3+nc,3+nc,3))
   allocate(taak(3+nc,3))
   allocate(rms(3+nc))

   return

end
