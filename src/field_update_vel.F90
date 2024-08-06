subroutine field_update_vel (uy,uny,nv)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute updated velocities and enforce continuity equation
! in wavenumber space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use comp
   implicit none

   integer nv
   complex(kind=b8) :: uy(xisz,ny,zjsz,nu)
   complex(kind=b8) :: uny(xisz,ny,zjsz,nv)
   complex(kind=b8) :: untemp,uytemp(3)
   real(kind=b8) :: factor,bk12,factxz
   integer :: i, yst

   ipsi=4+nc

   !$OMP TARGET DATA MAP(tofrom:uy,uny) MAP(to:mask,xfac,yfac,zfac, &
   !$OMP kx,ky,kz,b11,b22,bk1,bk3,bkk3) MAP(alloc:uytemp) 

   if (kstep.eq.1) then

      ! update the modified velocities (per Rogallo's scheme)
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(4) &
      !$OMP PRIVATE(i,xp,zp,y,x,z,factor,untemp,factxz) SHARED(uy, &
      !$OMP uny,mask,ny,zjsz,xisz,xist,zjst,xfac,yfac,zfac)
      do i=1,3
         do zp=1,zjsz
            do y=1,ny
               do xp=1,xisz
                  x=xist+xp-1
                  z=zjst+zp-1
                  factxz=xfac(x)*zfac(z)
                  factor=xfac(x)*yfac(y)*zfac(z)*mask(xp,y,zp)
                  untemp=uny(xp,y,zp,i)
                  uny(xp,y,zp,i)=factor*(untemp+uy(xp,y,zp,i))
                  uy(xp,y,zp,i)=factor*(untemp+uy(xp,y,zp,i)*2.)
               end do
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   else if(kstep.eq.2) then

      ! update the velocities
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(4) &
      !$OMP PRIVATE(i,xp,zp,y,x,z) SHARED(uy,uny,mask,ny,zjsz,xisz,xist, &
      !$OMP zjst,nc)
      do i=1,3+nc
         do zp=1,zjsz
            do y=1,ny
               do xp=1,xisz
                  uy(xp,y,zp,i) = &
                     (uny(xp,y,zp,i)+uy(xp,y,zp,i))*mask(xp,y,zp)
               end do
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   end if

   if (taskid.eq.0) then ! x=y=z=1
      !$OMP TARGET
      uytemp(1)=uy(1,1,1,1)
      uytemp(2)=uy(1,1,1,2)
      uytemp(3)=uy(1,1,1,3)
      !$OMP END TARGET
   end if

   ! form true velocities enforcing continuity
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
   !$OMP PRIVATE(xp,zp,y,x,z,bk12) SHARED(uy,uny,ny,zjsz,xisz, &
   !$OMP xist,zjst,kx,ky,kz,b22,bk1,bk3,bkk3,ipsi)
   do zp=1,zjsz
      do y=1,ny
         do xp=1,xisz
            x=xist+xp-1
            z=zjst+zp-1
            ! detrmine modified pressure from poisson equation
            ! Solve for i*pressure
            bk12=bk1(x,2)*kx(x)+b22(2)*ky(y)**2
            uy(xp,y,zp,ipsi)=-(bk1(x,2)*uy(xp,y,zp,1) + &
                        b22(2)*ky(y)*uy(xp,y,zp,2) + &
                        bk3(z,2)*uy(xp,y,zp,3))/(bk12+bkk3(z,2))
            uy(xp,y,zp,1)=uy(xp,y,zp,1)+kx(x)*uy(xp,y,zp,ipsi)
            uy(xp,y,zp,2)=uy(xp,y,zp,2)+ky(y)*uy(xp,y,zp,ipsi)
            uy(xp,y,zp,3)=uy(xp,y,zp,3)+kz(z)*uy(xp,y,zp,ipsi)
         end do
      end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   if (taskid.eq.0) then
      !$OMP TARGET
      uy(1,1,1,1)=uytemp(1)
      uy(1,1,1,2)=uytemp(2)
      uy(1,1,1,3)=uytemp(3)
      uy(1,1,1,ipsi)=cmplx(0.,0.)
      !$OMP END TARGET
   end if

   ! At last RK-substage only

   if (kstep.eq.2) then
      ! Set uny equal to uy
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(4) &
      !$OMP PRIVATE(i,xp,zp,y) SHARED(uny,uy,nv,xisz,zjsz,ny)
      do i=1,nv
         do zp=1,zjsz
            do y=1,ny
               do xp=1,xisz
                  uny(xp,y,zp,i)=uy(xp,y,zp,i)
               end do
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      ! enforce zero mean
      if(taskid.eq.0) then
         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) PRIVATE(i) &
         !$OMP SHARED(uy,uny,nv)
         do i=1,nv
            uy(1,1,1,i)=cmplx(0.,0.)
            uny(1,1,1,i)=cmplx(0.,0.)
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
      end if

   end if

   !$OMP END TARGET DATA

   ! if(taskid.eq.0) write(6,*)"Exit Field Update"

return
end
