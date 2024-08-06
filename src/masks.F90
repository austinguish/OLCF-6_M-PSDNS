subroutine masks

   ! pencils version, wavenumber indices in the order 'y,z,x'.

   use comp
   implicit none
   real(b8) :: sqkx,sqkz,sqk,two9,rnxi,rnyi,rnzi
 
   rnxi = 1./nx
   rnyi = 1./ny
   rnzi = 1./nz
   two9 = 2./9.

 
   do xp=1,xisz
      x=xp+xist-1
      sqkx = ( kx(x)*rnxi )**2 
      do zp=1,zjsz
         z=zp+zjst-1
         sqkz = ( kz(z)*rnzi )**2 + sqkx
         do y=1,ny
            sqk = ( ky(y)*rnyi )**2 + sqkz               
            if( sqk .gt. two9 .or. y .eq. nyhp .or. z .eq. nzhp) then
               mask(xp,y,zp)=0.
            else
               mask(xp,y,zp)=1.
            end if
         end do
      end do
   end do
 
   return

end subroutine masks
