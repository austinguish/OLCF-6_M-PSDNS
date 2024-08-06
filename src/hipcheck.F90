subroutine hipCheck(hipError_t)
   use hipfort
   implicit none
   integer(kind(hipSuccess)) :: hipError_t
   if(hipError_t /= hipSuccess)then
      write(*,*) "HIP ERROR: Error code = ", hipError_t
      call exit(hipError_t)
   end if
end subroutine hipCheck
