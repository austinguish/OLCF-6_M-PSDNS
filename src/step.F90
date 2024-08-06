      subroutine step

      ! Created by PKY, 11/6/2021
 
      ! For flows with no mean shear only
 
 
      use comp
      implicit none
 
      real, allocatable :: logx(:),logy(:),logz(:),bk12(:)
      real, allocatable :: vmaxz(:)
 
      real tmpvmaxz
      real logxy
 
      integer ithr,omp_get_thread_num

      real*8 dtout_d, time_d, dt_d, tfout_d, strof_d, t1_d, t2_d


      real term1,term2,term3,strof,dt0,term,t1,t2,cn
      integer it1,it2,ist1,ist2,ita,itb,i,ierror


      data term1,term2,logxy/3*0./
      data strof/1.e-6/
 
      allocate (logx(nxhpad))
      allocate (logy(nypad))
      allocate (logz(nzpad))
      allocate (bk12(nxhpad))
      allocate (vmaxz(numtasks))
 
      if (dtout.gt.0) iostep=99999999

      if (jstep.eq.1.and.cfl.le.0.) dt0=dt
 
      call MPI_ALLREDUCE (velmax,tmpvmaxz,1,mpireal,MPI_MAX, &
                          MPI_COMM_WORLD,mpierr)
      velmax=tmpvmaxz

      if (velmax.lt.0.001) then
         write (6,*) 'Abort due to velmax too small'
         call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if
      if (velmax.gt.1.e8) then
         write (6,*) 'Abort due to velmax too large'
         call MPI_ABORT (MPI_COMM_WORLD,ierror)
      end if


      if (cfl.gt.0.) then
         dt=cfl/velmax
      end if
 
      ! advance mesh metric  --------------------------------------------
      ! (constant irrotational strain or pure shear)
 
      ! b11 = b(1,1)**2, etc. , b12 = b(1,2) / b(2,2)
 
      b11(1)=b11(2)
      b22(1)=b22(2)
      b33(1)=b33(2)
      b12(1)=b12(2)
 
      ioflag=0

      ! set ioflag at istep equal to multiples of iostep
      ioflag = 0
      if (iostep.gt.0) then
         if (mod(istep-istep0,iostep).eq.0) ioflag = 1
      end if
 
      ! set tfout and tfsave to time0 if they were initialized with
      ! negative values
 
      if (jstep.eq.1.and.tfout.lt.0.) tfout=time0
      if (jstep.eq.1.and.tfsave.lt.0.) tfsave=time0
 
      ! adjust for output at times equal to multiples of dtout since tfout
      dtout_d=dtout
      tfout_d=tfout
      time_d=time
      dt_d=dt
      strof_d=1.e-5
      if (dtout.gt.0) then
         t1_d=time_d
         t2_d=time_d+dt_d
         it1=(t1_d-tfout_d)/dtout_d
         it2=(t2_d-tfout_d)/dtout_d
         ita=it1
         itb=ita+1
         if (abs(t1_d-itb*dtout_d+tfout_d).le.strof_d) it1=itb
         ita=it2
         itb=ita+1
         if (abs(t2_d-itb*dtout_d+tfout_d).le.strof_d) it2=itb
         ioflag=-1
         if (it2.gt.it1) then
            dt_d=(it1+1)*dtout_d+tfout_d-time_d
            ioflag=1
            dt=dt_d
            if (taskid.eq.0) write (6,*) ' dt adjusted in dtout loop, dt=',dt,dt_d
         end if
      end if
 
      ! set the "save" flag at times equal to multiples of dtsave since tfsave
 
      isflag=0
 
      if (dtsave.gt.0.) then
         t1=time
         t2=time+dt
         if (abs(mod(t2-tfsave,dtsave)).le.strof) isflag=1
         it1=(t1-tfsave)/dtsave
         it2=(t2-tfsave)/dtsave
         ita=it1
         itb=ita+1
         if (abs(t1-itb*dtsave+tfsave).le.strof) it1=itb
         ita=it2
         itb=ita+1
         if (abs(t2-itb*dtsave+tfsave).le.strof) it2=itb
         if (it2.gt.it1) then
            ist2=it1+1
            dt=(it2*dtsave+tfsave)-time
            ioflag=1
            isflag=1
            if (taskid.eq.0)write (6,*)' dt adjusted in dtsave loop, dt=',dt
         end if
      end if
 
      b11(2)=b11(1)*exp(-2.*a11*dt)
      b22(2)=b22(1)*exp(-2.*a22*dt)
      b33(2)=b33(1)*exp(-2.*a33*dt)
      b12(2)=b12(1)
 
      ! increment the time
 
      t1=time
      time=time+dt
      if (dt.eq.0.) then
         write (6,*) ' abort: dt unexpectedly zero in sub. step, istep=', istep
      end if
      t2=time
      dt2=dt/2.
 
        
      ! print time step or courant number
 
      if ( cfl .gt. 0. ) then
 
         !  dt from cfl
         if (taskid.eq.0) write (6,601) istep,dt,time,cfl,velmax
 601     format (' istep,dt,time,cfl,velmax=',i6,1p,4e12.5)
      
         if (dt.gt.10..or.dt.lt.1.e-7) then
            write (6,*) 'dt gotten out of control in sub. step, dt=',dt
            call MPI_ABORT (MPI_COMM_WORLD,ierror)
         end if
      
      else if( cfl .lt. 0. ) then
 
         !  constant dt, monitor cfl
         cn=velmax*dt/2./pi
         if (taskid.eq.0) write (6,602) istep,dt,time,cn,velmax
 602     format (' istep,dt,time,courant no=',i6,1p,4e12.5)

         if (cn.gt.1) then 
            call MPI_FINALIZE (mpierr)
            stop 'courant no. limit exceeded'
         end if
      
      endif
      
      !$OMP PARALLEL DEFAULT(NONE) PRIVATE(x,y,z,i) SHARED(term,&
      !$OMP b11,b22,b33,dt,nxh,ny,nz,nc,sc,logx,logy,logz,xfac,&
      !$OMP yfac,zfac,xdif,ydif,zdif,viscos,kx2,ky2,kz2)

      term=b11(1)*dt
      !$OMP DO
      do x=1,nxh
         logx(x)=-viscos*kx2(x)*term
         xfac(x)=exp(logx(x))
      end do
      !$OMP END DO

      term=b22(1)*dt
      !$OMP DO
      do y=1,ny
         logy(y)=-viscos*ky2(y)*term
         yfac(y)=exp(logy(y))
      end do
      !$OMP END DO

      term=b33(1)*dt
      !$OMP DO
      do z=1,nz
         logz(z)=-viscos*kz2(z)*term
         zfac(z)=exp(logz(z))
      end do
      !$OMP END DO

      if(nc.gt.0) then
         !$OMP DO
         do x=1,nxh
            do i=1,nc
               xdif(x,i)=xfac(x)**(1./sc(i))
            end do
         end do
         !$OMP END DO

         !$OMP DO
         do y=1,ny
            do i=1,nc
               ydif(y,i)=yfac(y)**(1./sc(i))
            end do
         end do
         !$OMP END DO

         !$OMP DO
         do z=1,nz
            do i=1,nc
               zdif(z,i)=zfac(z)**(1./sc(i))
            end do
         end do
         !$OMP END DO
      end if

      !$OMP END PARALLEL
 
      deallocate (logx,logy,logz,bk12,vmaxz)

      return
      end
