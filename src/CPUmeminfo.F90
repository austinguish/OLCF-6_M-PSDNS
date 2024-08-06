subroutine CPUmeminfo(str, n)
   
   use mpicom
   use ISO_C_BINDING
   implicit none

   integer :: n
   character(len=n) :: str
   integer :: valueRSS, sum_valueRSS
   character(len=200):: filename=' '
   character(len=80) :: line
   character(len=8)  :: pid_char=' '
   integer :: pid
   logical :: ifxst
   integer(kind=4) :: getpid

   valueRSS=-1    ! return negative number if not found

   !--- get process ID

   pid=getpid()
   write(pid_char,'(I8)') pid
   filename='/proc/'//trim(adjustl(pid_char))//'/status'
   if(taskid.eq.0) write(6,*) "proc status file name= ",filename

   !--- read system file
   
   inquire (file=filename,exist=ifxst)
   if (.not.ifxst) then
      write (*,*) 'system file does not exist'
      return
   end if

   open(unit=1000, file=filename, action='read')
   do
      read (1000,'(a)',end=120) line
      if (line(1:6).eq.'VmRSS:') then
         read (line(7:),*) valueRSS
         exit
      end if
   end do

120 continue

   close(1000)

   call MPI_REDUCE(valueRSS, sum_valueRSS, 1, MPI_INTEGER, MPI_SUM, &
      0, MPI_COMM_WORLD, ierr)

   if (taskid.eq.0) then
      open (9,file="CPUMemInfo")
      write (9,"(a25,'- Used= ',f6.2)") str(1:n), valueRSS/1024./1024.
      call custom_flush(9)
      write (0,"(a25,'- Used= ',f6.2)") str(1:n), valueRSS/1024./1024.
   end if

   return
end subroutine CPUmeminfo
