subroutine waveno
   !
   !  routine to form differential operators
   !
   !       bk1 = b(1,1)**2 * kx
   !       bk2 = b(2,2)**2 * ky + b(1,2) * b(2,2) * kx
   !       bk3 = b(3,3)**2 * kz
   !
   !       bkk12 + bkk3 = laplacian   
   !
   !  loop over meshes
   !

   use comp
   implicit none

   integer m

   do m=1,2
      ! form wave numbers for differentiation
      do x=1,nxh
         bk1(x,m)=b11(m)*kx(x)
      end do

      do y=1,ny
         bk2(y,m)=b22(m)*ky(y)
      end do

      do z=1,nz
         bk3(z,m)=b33(m)*kz(z)
         bkk3(z,m)=bk3(z,m)*kz(z)
      end do

   end do

   do x=1,nxh
      kx2(x)=kx(x)**2
   end do
   do y=1,ny
      ky2(x)=ky(y)**2
   end do
   do z=1,nz
      kz2(x)=kz(z)**2
   end do

   return
end
