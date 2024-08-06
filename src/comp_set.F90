subroutine comp_set

   use comp

   implicit none

   allocate (mask(xisz,ny,zjsz))
   allocate (xfac(nxh),yfac(ny),zfac(nz))
  
   if (taskid.eq.0) then
      write (6,*) 'nxhpad,nypad,nzpad,nu=',nxhpad,nypad,nzpad,nu
   end if

   allocate (un(nxpad,zisz,yjsz,3+nc),stat=ierr)
   if(taskid.eq.0) write (6,*) 'After allocate un'
   call custom_flush(6)

   if (ierr.ne.0) call abrt('error allocating un in comp_set')

   allocate (u(nxpad,zisz,yjsz,nu),stat=ierr)
   if(taskid.eq.0) write (6,*) 'After allocate u'
   call custom_flush(6)

   if (ierr.ne.0) call abrt('error allocating u in comp_set')

   un(:,:,:,:)=0.
   call CPUmeminfo ('after allocating un',19)
   u(:,:,:,:)=0.
   call CPUmeminfo ('after allocating u',18)

   if (i_press.eq.1) then
      allocate (up(nxpad,zisz,yjsz),stat=ierr)
      if(taskid.eq.0) write (6,*) 'After allocate up'
      call CPUmeminfo ('after allocating up',19)
   end if

   return
end
