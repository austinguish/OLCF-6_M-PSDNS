subroutine blanks (char,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! routine to take out embedded blanks in a character string, and
! keep count of the number of resulting leading non-blank characters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   character(len=*) :: char
   integer :: n,j,i

   n=len(char)
   j=0
   do i=1,n
      if (char(i:i).ne.' ') then
         j=j+1
         char(j:j)=char(i:i)
      end if
   end do

return
end
