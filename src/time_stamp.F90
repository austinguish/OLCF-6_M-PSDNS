subroutine time_stamp (string)
 
   use mpicom, only: taskid
 
   implicit none

   integer dmy(3),hms(3),date_time(8), nchar
 
   character*(*) string
 
   call date_and_time (values=date_time)
   dmy(1:3)=date_time(3:1:-1)
   hms(1:3)=date_time(5:7)
 
   nchar=min(len(string),25)
   if (taskid.eq.0) then
      write (6,10) string(1:nchar),taskid,dmy(2),dmy(1),dmy(3),hms
10    format (a25,2x,i5, ' date & time is  ',i2,'/',i2, &
              '/',i4,2x,i2,':',i2,':',i2)
   end if
 
   return
end
 
