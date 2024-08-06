real function ranu (ii)
   ! routine to generate independent random numbers
   ! uniformly distributed in the exclusive interval (0,1).
 
   ! the argument ii is a dummy

   use com
   use ranseed

   save icount
   data icount/0/
   
   if (iseed.gt.0) then
      ! generate random number
      iseed=iseed*65539
      if (iseed.lt.0) iseed=iseed+2147483647+1
      ranu=iseed*0.46566127e-9
      return
   else
      write(6,*)' seed iseed=',iseed,' must be positive',taskid
      write(6,*)' stopped in ranu '
      call MPI_ABORT (MPI_COMM_WORLD,ierror)
   end if

end

real function cancel (delta)
   ! routine to generate the grid shift to cancel the grid shift from the
   ! previous step

   real :: delta

   cancel = mod( delta + 0.5 , 1. )
   return
end
