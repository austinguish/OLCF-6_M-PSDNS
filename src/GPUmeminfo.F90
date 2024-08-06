subroutine GPUmeminfo(str, n)

   use mpicom
   use ISO_C_BINDING

   use hipfort
   use hipfort_check
   implicit none

   integer :: n
   character(len=n) :: str
   integer(kind=c_size_t) :: free, total

   free = 0
   total = 0

   call hipCheck(hipMemGetInfo_(free,total))

   if (taskid.eq.0) then
      open (8,file="GPUMemInfo")
      write (8,"(a25,'- Used= ',f8.4,' Free= ',f8.4,' Total= ',f8.4)") &
         str(1:n), (total-free)/1024./1024./1024., &
         free/1024./1024./1024., total/1024./1024./1024.
      call custom_flush(8)
      write (0,"(a25,'- Used= ',f8.4,' Free= ',f8.4,' Total= ',f8.4)") &
         str(1:n), (total-free)/1024./1024./1024., &
         free/1024./1024./1024., total/1024./1024./1024.
   end if

return
end
