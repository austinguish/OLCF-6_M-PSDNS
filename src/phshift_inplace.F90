subroutine phshift_inplace (uy,nv,ev)
 
   use comp
   implicit none
 
   integer nv,ev
   complex(b8) uy(xisz,ny,zjsz,nv)

   ! apply phase-shifting only
 
   complex(b8) sxyz
 
   integer i
 

   !$OMP TARGET DATA MAP(tofrom:uy) MAP(to:mask,sx,sy,sz)

   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) DEFAULT(NONE) PRIVATE(xp,x,zp,z,y, &
   !$OMP i,sxyz) SHARED (xist,xisz,zjst,zjsz,ny,ev,mask,sx,sy,sz,nv,uy)
   do i=1,nv
      do zp=1,zjsz
         do y=1,ny
            do xp=1,xisz
               z=zjst+zp-1
               x=xist+xp-1
               sxyz=sz(z,ev)*sx(x,ev)*sy(y,ev)*mask(xp,y,zp)
               uy(xp,y,zp,i)=uy(xp,y,zp,i)*sxyz
            end do
         end do
      end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   !$OMP END TARGET DATA

   return
end
