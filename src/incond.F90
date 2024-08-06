
subroutine incond 

   ! routine to manage setting of initial conditions 

   use comp
   implicit none

   if (kinit(1).eq.-2) then
      if (taskid.eq.0) write (6,*) 'incond: calling pfield'
      call pfield (un,1)
   end if

   return

end
