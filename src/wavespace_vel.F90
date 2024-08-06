subroutine wavespace_vel (uy,m)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Take derivative in y direction.
! Form the three convective terms for wavenumbers less than the cut off
! wavenumber (defined by mask array)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use comp
   implicit none

   integer :: m,i
   real(kind=b8) :: norm
   complex(kind=b8) :: uy(xisz,ny,zjsz,nu)
   complex(kind=b8) :: utmp1,szx,sxyz
   real(kind=8) :: temp1,rtmap,rtconv

   is1=4+nc

   norm=dt2/ny/nz

   !$OMP TARGET DATA MAP(tofrom:uy) MAP(to:mask,sx,sz,sy,bk1,bk3,b22)

   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
   !$OMP PRIVATE(zp,z,y,szx,x,xp,i,sxyz,utmp1) SHARED(uy,mask,ny,nxh, &
   !$OMP nc,kstep,taskid,sx,sy,sz,norm,imagi,b22,bk1,bk3,is1,m,ky,xist, &
   !$OMP xisz,zjst,zjsz,istep)
   do zp=1,zjsz
      do y=1,ny
         do xp=1,xisz

            z=zjst+zp-1
            x=xist+xp-1

            szx=sz(z,kstep)*sx(x,kstep)
            sxyz=-norm*conjg(szx*sy(y,kstep))*mask(xp,y,zp) ! -dt/2
            utmp1=imagi*b22(m)*ky(y)*sxyz ! i*b22*ky
               
            ! c1 = -i*dt/2*(kx*(u^2-v^2)+ky*u*v+kz*u*w)
            uy(xp,y,zp,1)=sxyz*uy(xp,y,zp,1)+utmp1*uy(xp,y,zp,2)
            ! c2 = -i*dt/2*(kx*u*v+kz*v*w)
            uy(xp,y,zp,2)=sxyz*imagi*(bk1(x,m)*uy(xp,y,zp,2)+&
               bk3(z,m)*uy(xp,y,zp,is1))
            ! c3 = -i*dt/2*(kx*u*w+ky*v*w+kz*(w^2-v^2))
            uy(xp,y,zp,3)=sxyz*uy(xp,y,zp,3)+utmp1*uy(xp,y,zp,is1) 

         end do
      end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO   

   !$OMP END TARGET DATA

   return
end
