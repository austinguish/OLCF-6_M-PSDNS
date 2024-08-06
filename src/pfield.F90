      subroutine pfield (unxr,icall)
      use com
        implicit none
!     
!
!     set initial velocity field in physical space
!     
!     sinusoidal velocity field (products of sines and cosines)
!     
!     u = au * sin(x)*cos(y)*cos(z)
!     v = av * cos(x)*sin(y)*cos(z)
!     w = aw * cos(x)*cos(y)*sin(z)
!     
!     continuity requires au + av + aw = 0.
!     
      real(b8) :: unxr(nx,zisz,yjsz,3)
!
      integer icall,m
      real diff1,diff2,diff3,difmax,globmax
      real au,av,aw,dx,dy,dz,termu,termv,termw,sinz,cosz,umax

     
      real(b8), allocatable :: sinx(:),siny(:),cosx(:),cosy(:)

      real cpu
      real*8 rtime1,rtime2

      integer ithr,OMP_GET_THREAD_NUM
!
      
      allocate (sinx(nx),stat=ierr)
      allocate (cosx(nx),stat=ierr)
      allocate (siny(ny),stat=ierr)
      allocate (cosy(ny),stat=ierr)

      data au,av/1.,1./

!     print to standard output the initial conditions in use
      if (taskid.eq.0) then
         write (6,*) 'enter pfield, taskid=',taskid, nx,zisz,yjsz
         write (6,*) "using hook pfield: sinusoidal"
      endif
      
        
      aw=-(au+av)

      pi=atan(1.)*4.
      dx=2.*pi/nx
      dy=2.*pi/ny
      dz=2.*pi/nz

      do x=1,nx
         sinx(x)=sin((x-1)*dx)
         cosx(x)=cos((x-1)*dx)
      enddo
      do y=1,ny
         siny(y)=sin((y-1)*dy)
         cosy(y)=cos((y-1)*dy)
      enddo

      umax=0.

      if (icall.eq.1) then
!
         rtime1=MPI_WTIME()

         !$OMP PARALLEL private (ithr,y,z,sinz,cosz,termu,termv,termw)
         ithr=OMP_GET_THREAD_NUM()

         !$OMP DO COLLAPSE (2)
         do yp=1,yjsz
            do zp=1,zisz

               y=yjst+yp-1
               z=zist+zp-1

               sinz=sin((z-1)*dz)
               cosz=cos((z-1)*dz)
     
               termu=cosz*cosy(y)
               termv=siny(y)*cosz
               termw=cosy(y)*sinz
            
               do x=1,nx
                  unxr(x,zp,yp,1)=au*sinx(x)*termu
                  unxr(x,zp,yp,2)=av*cosx(x)*termv
                  unxr(x,zp,yp,3)=aw*cosx(x)*termw
               enddo               !done with x

            enddo                  !done with zp
         end do
         !$OMP END DO

         !$OMP END PARALLEL

         rtime2=MPI_WTIME()
         cpu=rtime2-rtime1

      else if (icall.eq.2) then
 
         difmax=0.
 
         do yp=1,yjsz
            y=yjst+yp-1
            do zp=1,zisz
               z=zist+zp-1
               sinz=sin((z-1)*dz)
               cosz=cos((z-1)*dz)
     
               termu=cosz*cosy(y)
               termv=siny(y)*cosz
               termw=cosy(y)*sinz
            
               do x=1,nx
                  diff1=abs(unxr(x,zp,yp,1)-au*sinx(x)*termu)
                  diff2=abs(unxr(x,zp,yp,2)-av*cosx(x)*termv)
                  diff3=abs(unxr(x,zp,yp,3)-aw*cosx(x)*termw)
                  difmax=amax1(difmax,diff1,diff2,diff3)
               enddo               !done with x

            enddo                  !done with zp
         enddo                     !done with yp
 
 
         call MPI_REDUCE (difmax,globmax,1,mpireal,MPI_MAX, &
                          0,MPI_COMM_WORLD,mpierr)
         if (taskid.eq.0) write (6,*) 'max global diff.=',globmax
 
      end if
 
      deallocate (sinx,cosx,siny,cosy)
      
      return
      end
