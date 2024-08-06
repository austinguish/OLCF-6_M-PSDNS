   subroutine param_set
!
      use param
      use mpicom, only : taskid
      implicit none
 
      ! note that nx, ny, nz, nc, kfor are to be read in
      ! (set up a subroutine input.f)
 
      !  derived parameters                             
                                                  
      nxh=nx/2 ; nxhp=nxh+1 ; nxhp2=nxh+2
      nyh=ny/2 ; nyhp=nyh+1 ; nyhp2=nyh+2
      nzh=nz/2 ; nzhp=nzh+1 ; nzhp2=nzh+2


      !! NPM edit: set dealiasing in input file

      if(nxpad.lt.nx) then !user disable dealiasing for x
         nxpad = nx
      endif
        
      if(nypad.lt.ny) then !user disable dealiasing for y
         nypad = ny
      endif
        
      if(nzpad.lt.nz) then !user disable dealiasing for z
         nzpad = nz     
      endif        
      if (taskid.eq.1) write (6,*) 'param_set: pads=',nxpad,nypad,nzpad

      ! If we are not using in-place transforms we can set nxpad=nx
      ! Otherwise, for in-place transforms nxpad=nx+2

      nxhpad=nxpad/2 ; nxhppad=nxhpad+1 ; nxhp2pad=nxhpad+2
      nyhpad=nypad/2 ; nyhppad=nyhpad+1 
      nzhpad=nzpad/2 ; nzhppad=nzhpad+1 
                                                 
      nu = 5

      ! is1=4+nc ; is2=5+nc ; ipsi=4+nc
                                                 
      ncpp=nc*nc+7*nc+12
      ncp=ncpp/2

      return
   end
