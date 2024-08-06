subroutine sptr_omp (uny)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! By Kiran, 1/29/2022
! Compute the 3D energy spectrum and dissipation spectrum
! Compute other parameters printed in eulstat file
! Note:- Code is written for beta_i=1
!        Combine all communication calls into two MPI_REDUCEs,
!        one for 3D spectra, the other for 1D spectra
!        Use Reduce instead of Allreduce
!        Assume nxhp=nyhp=nzhp
! Currently uses CPU threads only, since GPU has issues with
! OMP REDUCTION on arrays. On Crusher, if we use 8 CPU threads
! per MPI task, time taken for each call to this routine
! is estimated to be approx twice the wall time of a single
! regular (non-output) time step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rohini helped generalize to case of nx> ny=nz, 2/4/2022
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   use comp
   implicit none
   
   complex(kind=b8) :: uny(xisz,ny,zjsz,3+nc)
   real(kind=b8), allocatable :: eijky(:,:),ekky(:)
   real(kind=b8), allocatable :: dijky(:,:)
   real(kind=b8), allocatable :: sijky(:,:)
   real(kind=b8), allocatable :: gijky(:,:),gijk(:,:)
   real(kind=b8), allocatable :: grijky(:,:,:)
   real(kind=b8) :: der2y(ncp),der2ij(ncp)
   real(kind=b8), allocatable :: xrijy(:,:)
   real(kind=b8), allocatable :: yrijy(:,:)
   real(kind=b8), allocatable :: zrijy(:,:)
   real(kind=b8) :: ssum,ssq(3+nc),s(3+nc),s12,s31,s23,s123,factor
   real(kind=b8) :: bkx2(nxh),bky2(ny),bkz2(nz),beta_min
   integer :: ik,xst,i,j,k,ij,nyp2,nzp2,iy,iz,ijc,a
   real(kind=b8) :: rk2,uiuj(6),term,dterm,fact
   ! for scalar statistics
   real(kind=b8) :: q2,tau
   integer :: is,nshell
   real(kind=8) :: rtime1,rtime2,rtime3,rtime4
   ! buffers for Reduce
   integer :: ncount
   real(kind=b8), allocatable :: sndbuf1(:,:,:),rcvbuf1(:,:,:)

   if(taskid.eq.0)write(6,*)"Enter Sptr"
   rtime1=MPI_WTIME()

   nyp2=ny+2
   nzp2=nz+2

   rtime2=MPI_WTIME()
   allocate(eijky(mxyz,ncp))
   allocate(ekky(mxyz))
   allocate(ekk(mxyz))
   allocate(dijky(mxyz,ncp))
   allocate(sijky(mxyz,ncp))
   allocate(grijky(mxyz,ncp,3))
   allocate(gijky(mxyz,ncp))
   allocate(gijk(mxyz,ncp))
   allocate(xrijy(nxhp,ncp))
   allocate(xrij(nxhp,ncp))
   allocate(yrijy(nyhp,ncp))
   allocate(yrij(nyhp,ncp))
   allocate(zrijy(nzhp,ncp))
   allocate(zrij(nzhp,ncp))
   rtime2=MPI_WTIME()-rtime2
   ! if (taskid.eq.0) write(6,"(' Sptr: Allocate (secs) :',1p,e12.4)") rtime2

   beta_min=min(beta1,beta2,beta3)
   ssq(1)=b11(2)
   ssq(2)=b22(2)
   ssq(3)=b33(2)

   if(nc.gt.0) then
      do i=4,3+nc
         ssq(i)=1.
      end do
   end if

   kx2(:)=kx(:)**2
   ky2(:)=ky(:)**2
   kz2(:)=kz(:)**2
   ky2(nyhp)=0.
   kz2(nzhp)=0.

   bkx2(:)=b11(2)*kx2(:)
   bky2(:)=b22(2)*ky2(:)
   bkz2(:)=b33(2)*kz2(:)

   eijky(:,:)=0.
   eijk(:,:)=0.
   ekky(:)=0.
   ekk(:)=0.
   dijky(:,:)=0.
   dijk(:,:)=0.
   sijky(:,:)=0.
   sijk(:,:)=0.
   gijky(:,:)=0.
   gijk(:,:)=0.
   grijky(:,:,:)=0.
   grijk(:,:,:)=0.
   xrijy(:,:)=0.
   xrij(:,:)=0.
   yrijy(:,:)=0.
   yrij(:,:)=0.
   zrijy(:,:)=0.
   zrij(:,:)=0.
   der2y(:)=0.
   der2ij(:)=0.
   uiuj=0.

   ij=0
   do i=1,3+nc
      do j=i,3+nc
         ij=ij+1
         if (i.eq.j) kt(i)=ij
      end do
   end do

   rtime2=MPI_WTIME()
   !$OMP PARALLEL DO DEFAULT(NONE) &
   !$OMP PRIVATE(zp,z,y,xst,xp,x,rk2,ik,ij,i,j,uiuj,term,dterm,iz,iy) &
   !$OMP SHARED(zjsz,zjst,xisz,xist,taskid,nzhp,ny,nyhp,nxh,bkx2,bky2, &
   !$OMP bkz2,nc,mask,uny,ssq,nyp2,nzp2,beta_min,tfact,istep) REDUCTION(+:eijky, &
   !$OMP dijky,der2y,sijky,gijky,grijky,ekky,xrijy,yrijy,zrijy) 

   do zp=1,zjsz

      z=zp+zjst-1
      do y=1,ny

         do xp=1,xisz

            if(mask(xp,y,zp).eq.0.) cycle

            x=xp+xist-1

            rk2=bkx2(x)+bky2(y)+bkz2(z)
            ik=sqrt(rk2)/beta_min+1.5

            uiuj(1)=real(uny(xp,y,zp,1)*conjg(uny(xp,y,zp,1)))
            uiuj(2)=real(uny(xp,y,zp,1)*conjg(uny(xp,y,zp,2)))
            uiuj(3)=real(uny(xp,y,zp,1)*conjg(uny(xp,y,zp,3)))
            uiuj(4)=real(uny(xp,y,zp,2)*conjg(uny(xp,y,zp,2)))
            uiuj(5)=real(uny(xp,y,zp,2)*conjg(uny(xp,y,zp,3)))
            uiuj(6)=real(uny(xp,y,zp,3)*conjg(uny(xp,y,zp,3)))

            ij=0
            do i=1,3
               do j=i,3
                  ij=ij+1
                  term=tfact(x)*uiuj(ij)
                  dterm=rk2*term

                  eijky(ik,ij)=eijky(ik,ij)+term
                  dijky(ik,ij)=dijky(ik,ij)+dterm
                  der2y(ij)=der2y(ij)+bky2(y)*term
                  sijky(ik,ij)=sijky(ik,ij)+dterm*rk2
                  gijky(ik,ij)=gijky(ik,ij)+dterm*rk2**2

                  if(i.eq.j.and.i.le.3.and.rk2.ne.0) then
                     ekky(ik)=ekky(ik)+term/sqrt(rk2)*ssq(i)
                  end if

                  ! spectra of individual gradients
                  grijky(ik,ij,1)=grijky(ik,ij,1)+term*bkx2(x)
                  grijky(ik,ij,2)=grijky(ik,ij,2)+term*bky2(y)
                  grijky(ik,ij,3)=grijky(ik,ij,3)+term*bkz2(z)

                  xrijy(x,ij)=xrijy(x,ij)+uiuj(ij)

               end do
            end do

            if (y.eq.1) then
               do ij=1,6
                  yrijy(y,ij)=yrijy(y,ij)+uiuj(ij)+(tfact(x)-1)*uiuj(ij)
               end do
            else if(y.lt.nyhp) then
               do ij=1,6
                  yrijy(y,ij)=yrijy(y,ij)+uiuj(ij)
               end do
            else
               do ij=1,6
                  iy=nyp2-y
                  yrijy(iy,ij)=yrijy(iy,ij)+(tfact(x)-1)*uiuj(ij)
               end do
            endif

            if(z.eq.1) then
               do ij=1,6
                  zrijy(z,ij)=zrijy(z,ij)+uiuj(ij)+(tfact(x)-1)*uiuj(ij)
               end do
            else if(z.lt.nzhp) then
               do ij=1,6
                  zrijy(z,ij)=zrijy(z,ij)+uiuj(ij)
               end do
            else
               do ij=1,6
                  iz=nzp2-z
                  zrijy(iz,ij)=zrijy(iz,ij)+(tfact(x)-1)*uiuj(ij)
               end do
            endif

         end do ! x
      end do ! y
   end do ! z
   !$OMP END PARALLEL DO

   rtime2=MPI_WTIME()-rtime2
   if (taskid.eq.0) write(6,"(' Sptr: Parallel region (secs) :',1p,e12.4)") rtime2
   rtime4=rtime2

   dijky(:,:)=dijky(:,:)*2.*viscos
   sijky(:,:)=sijky(:,:)*2.*viscos
   gijky(:,:)=gijky(:,:)*2.*viscos

   rtime2=MPI_WTIME()

   ! pack data
   rtime3=MPI_WTIME()
   allocate (sndbuf1(mxyz,ncp,8),rcvbuf1(mxyz,ncp,8))
   do i=1,mxyz
      do j=1,ncp
         sndbuf1(i,j,1)=eijky(i,j)
         sndbuf1(i,j,2)=dijky(i,j)
         sndbuf1(i,j,3)=sijky(i,j)
         sndbuf1(i,j,4)=gijky(i,j)
         sndbuf1(i,j,5)=grijky(i,j,1)
         sndbuf1(i,j,6)=grijky(i,j,2)
         sndbuf1(i,j,7)=grijky(i,j,3)
      end do
      sndbuf1(i,1,8)=ekky(i)
   end do
   do j=1,ncp
      sndbuf1(j,2,8)=der2y(j)
   end do
   rtime3=MPI_WTIME()-rtime3

   rtime3=MPI_WTIME()
   call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   rtime3=MPI_WTIME()-rtime3

   rtime3=MPI_WTIME()
   ncount=mxyz*ncp*7+mxyz+ncp
   call MPI_REDUCE (sndbuf1,rcvbuf1,ncount,mpireal,MPI_SUM,0,&
      MPI_COMM_WORLD,mpierr)
   rtime3=MPI_WTIME()-rtime3

   ! unpack data
   rtime3=MPI_WTIME()
   do i=1,mxyz
      do j=1,ncp
         eijk(i,j)=rcvbuf1(i,j,1)
         dijk(i,j)=rcvbuf1(i,j,2)
         sijk(i,j)=rcvbuf1(i,j,3)
         gijk(i,j)=rcvbuf1(i,j,4)
         grijk(i,j,1)=rcvbuf1(i,j,5)
         grijk(i,j,2)=rcvbuf1(i,j,6)
         grijk(i,j,3)=rcvbuf1(i,j,7)
      end do
      ekk(i)=rcvbuf1(i,1,8)
   end do
   do j=1,ncp
      der2ij(j)=rcvbuf1(j,2,8)
   end do
   deallocate (sndbuf1,rcvbuf1)
   rtime3=MPI_WTIME()-rtime3

   ! pack data, assume nxhp=nyhp=nzhp
   rtime3=MPI_WTIME()
   allocate (sndbuf1(nxhp,ncp,3),rcvbuf1(nxhp,ncp,3))
   sndbuf1=0.
   do i=1,ncp
      do x=1,nxhp
         sndbuf1(x,i,1)=xrijy(x,i)
      end do
      do y=1,nyhp
         sndbuf1(y,i,2)=yrijy(y,i)
      end do
      do z=1,nzhp
         sndbuf1(z,i,3)=zrijy(z,i)
      end do
   end do
   rtime3=MPI_WTIME()-rtime3

   rtime3=MPI_WTIME()
   call MPI_BARRIER(MPI_COMM_WORLD,mpierr)
   rtime3=MPI_WTIME()-rtime3

   rtime3=MPI_WTIME()
   ncount=nxhp*ncp*3
   call MPI_REDUCE (sndbuf1,rcvbuf1,ncount,mpireal,MPI_SUM,0,&
      MPI_COMM_WORLD,mpierr)
   rtime3=MPI_WTIME()-rtime3

   ! unpack data
   rtime3=MPI_WTIME()
   do i=1,ncp
      do x=1,nxhp
         xrij(x,i)=rcvbuf1(x,i,1)
      end do
      do y=1,nyhp
         yrij(y,i)=rcvbuf1(y,i,2)
      end do
      do z=1,nzhp
         zrij(z,i)=rcvbuf1(z,i,3)
      end do
   end do
   deallocate (sndbuf1,rcvbuf1)
   rtime3=MPI_WTIME()-rtime3

   rtime2=MPI_WTIME()-rtime2
   if (taskid.eq.0) write(6,"(' Sptr: MPI (secs) :',1p,e12.4)") rtime2

   s(1)=sqrt(b11(2))
   s(2)=sqrt(b22(2))
   s(3)=sqrt(b33(2))
   s23=s(2)*s(3)
   s31=s(3)*s(1)
   s12=s(1)*s(2)
   s123=s(1)*s(2)*s(3)
   if(nc.gt.0) then
      do i=4,3+nc
         s(i)=1.
      end do
   end if

   ! Account for non-unity beta
   ij=0
   do i=1,3+nc
      do j=i,3+nc
         ij=ij+1
         factor=s(i)*s(j)

         der2ij(ij)=der2ij(ij)*factor
         eijk(:,ij)=eijk(:,ij)*factor
         dijk(:,ij)=dijk(:,ij)*factor
         sijk(:,ij)=sijk(:,ij)*factor
         gijk(:,ij)=gijk(:,ij)*factor

         xrij(:,ij)=factor*xrij(:,ij)/s123
         yrij(:,ij)=factor*yrij(:,ij)/s123
         zrij(:,ij)=factor*zrij(:,ij)/s123
      end do
   end do
   if (taskid.eq.0) write(6,*) 'Sptr: After accounting for non-unity beta'
   ek(:)=(eijk(:,1)+eijk(:,4+nc)+eijk(:,6+2*nc))*0.5 ! energy spectrum
   ek(1)=0.
   dk(:)=(dijk(:,1)+dijk(:,4+nc)+dijk(:,6+2*nc))*0.5 ! dissipation spectrum
   dk(1)=0.
   sk(:)=(sijk(:,1)+sijk(:,4+nc)+sijk(:,6+2*nc))*0.5 ! dissipation skewness spectrum
   ekl=0.
   do x=1,nxh
      ekl=ekl+ekk(x)/2.
   end do
   if (taskid.eq.0) write(6,*) 'Sptr: After forming energy/diss/diss skew spec'
   tke=sum(ek) ! turbulent kinetic energy
   epslon=sum(dk) ! dissipation rate
   if(nc.gt.0) then
      do i=1,nc
         dijk(:,kt(i+3))=dijk(:,kt(i+3))/pr(i)
         scdiss(i)=sum(dijk(:,kt(i+3))) ! scalar dissipation rate
      end do
   end if
   skew=sum(sk)
   klen=(viscos**3/epslon)**0.25 ! kolmogorov length scale
   kvel=(viscos*epslon)**0.25 ! kolmogorov velocity scale
   ktime=(viscos/epslon)**0.5 ! kolmogorov time scale
   skew=skew*2./35.*(15.*viscos/epslon)**1.5 ! dissipation skewness
   ekl=1.5*pi/(2*tke)*ekl ! integral scale as integral of e(k)/k

   do ij=1,ncp

      xrij(:,ij)=s23*xrij(:,ij)
      yrij(:,ij)=s31*yrij(:,ij)
      zrij(:,ij)=s12*zrij(:,ij)

      ssum=real(xrij(1,ij))
      do x=2,nxh
         ssum=ssum+2.*real(xrij(x,ij))
      end do
      corr(ij)=s(1)*ssum
   end do
   
   !!!! Compute rms value of velocity and scalar fluctuations
   do i=1,3+nc
      rms(i)=sqrt(corr(kt(i)))
   end do

   !!!! Compute integral length scales
   lijk(:,:)=0.
   do ij=1,ncp
      if(corr(ij).ne.0.) then
         lijk(ij,1)=pi*real(xrij(1,ij))/corr(ij)
         lijk(ij,2)=pi*real(yrij(1,ij))/corr(ij)
         lijk(ij,3)=pi*real(zrij(1,ij))/corr(ij)
      end if
   end do

   do j=1,3
      do i=1,3+nc
         laak(i,j)=lijk(kt(i),j)
      end do
   end do

   !!!! Compute the eddy-turnover times
   do i=1,3
      if(rms(i).ne.0) then
         ett(i)=laak(i,i)/rms(i)
      end if
   end do

   !!!! Compute integral scale reynolds number
   if(viscos.gt.0) then
      do i=1,3
         raak(i)=laak(i,i)*rms(i)/viscos
      end do
   end if

   !!!! Compute the mean square of fluctuating strain rate (dui/dxj)
   do i=1,3+nc

      ssum=0.
      do x=2,nxh
         ssum=ssum+kx2(x)*real(xrij(x,kt(i)))
      end do
      tmij(i,1)=2.*ssum*s(1)**3

      tmij(i,2)=der2ij(kt(i))

      ssum=0.
      do z=2,nzh
         ssum=ssum+kz2(z)*real(zrij(z,kt(i)))
      end do
      tmij(i,3)=2.*ssum*s(3)**3

   end do

   !!!! Compute taylor microscales based corr & tmij
   do i=1,3+nc
      do j=1,3+nc
         fact=1.
         if(i.gt.3.or.j.gt.3) fact=sqrt(2.)
         do k=1,3
            taylor(i,j,k)=0.
            if(tmij(j,k).ne.0) taylor(i,j,k)=sqrt(corr(kt(i))/tmij(j,k))*fact
         end do
      end do
   end do

   do j=1,3
      do i=1,3+nc
         taak(i,j)=taylor(i,j,j)
      end do
   end do

   !!!! Compute the taylor microscale reynolds numbers
   if(viscos.gt.0) then
      do i=1,3
         tmre(i)=taak(i,i)*rms(i)/viscos
      end do
   end if

   !!!! Compute Scalar statistics
   if(nc.gt.0) then

      q2=corr(kt(1))+corr(kt(2))+corr(kt(3)) ! u^2+v^2+w^2
      tau=q2/epslon ! mechanical time scale
      ! line below commented out to match SHELL_DK compiler flag which is
      ! commonly used in PSDNS_HOMO_OMP code
      !nshell=nxh*min(beta1,beta2,beta3)
      nshell=nxh

      do is=1,nc
         i=3+is

         ! mean of scalar field
         if(taskid.eq.0) scmean(is)=real(uny(1,1,1,i))

         ! scalar time scale and mechanical-to-scalar time scale ratio
         sctime(is)=0.
         mtosratio(is)=0.
         if(scdiss(is).gt.0) then
            sctime(is)=corr(kt(i))/scdiss(is)
            mtosratio(is)=tau/sctime(is)
         end if

         ! Batchelor and Obukhov-Corrsin scales
         batchelor(is)=klen/sqrt(pr(is))
         obukc(is)=klen/pr(is)**.75

         ! calculate mean of fourth and sixth order derivative of the scalar
         spm4(is)=0.
         spm6(is)=0.
         do j=1,nshell
            spm4(is)=spm4(is)+sijk(j,kt(i))
            spm6(is)=spm6(is)+gijk(j,kt(i))
         end do
         spm4(is)=spm4(is)/2./viscos
         spm6(is)=spm6(is)/2./viscos
      end do

      ! compute scalar gradient variances in wavenumber space
      ij=kt(4)-1
      ijc=0
      do i=1,nc
         do j=i,nc
            ij=ij+1
            ijc=ijc+1
            do a=1,3
               scgvar(ijc,a)=0.
               do k=1,nxh
                  scgvar(ijc,a)=scgvar(ijc,a)+grijk(k,ij,a)
               end do
            end do
         end do
      end do
      write (30+taskid,*) "scgvar(:,1)=",scgvar(:,1)
      write (30+taskid,*) "scgvar(:,2)=",scgvar(:,2)
      write (30+taskid,*) "scgvar(:,3)=",scgvar(:,3)

   end if


   !!!! Write out calculated quantities to eulstat file
   if (taskid.eq.0) then
      call eulout()
   end if

   deallocate(eijky)
   deallocate(ekky)
   deallocate(ekk)
   deallocate(dijky)
   deallocate(sijky)
   deallocate(gijky)
   deallocate(gijk)
   deallocate(grijky)
   deallocate(xrijy)
   deallocate(xrij)
   deallocate(yrijy)
   deallocate(yrij)
   deallocate(zrijy)
   deallocate(zrij)

   rtime1=MPI_WTIME()-rtime1
   if (taskid.eq.0) write(6,"(' Exit Sptr (secs) :',1p,e12.4)") rtime1

return
end
