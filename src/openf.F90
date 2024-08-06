subroutine openf

   use comp
   implicit none
   
   character(len=5) :: numer
   character(len=16) :: fn
   integer :: nchar,i,ii
   character(len=11) :: form,unform
   data form,unform/'formatted','unformatted'/
   character(len=20) :: caux

   if (taskid.eq.0) then

      caux='eulstat'
      io_unit_eulstat=16
      call fopen1 (io_unit_eulstat,caux,form,0)
      close(io_unit_eulstat)

      if(i_dissenst.eq.1) then
         caux='dns_max_epsenst'
         io_unit_maxepsenst=19
         call fopen1 (io_unit_maxepsenst,caux,form,0)
         !close(io_unit_maxepsenst)
 
         caux='dns_pdf_epsenst'
         io_unit_pdfepsenst=20
         call fopen1 (io_unit_pdfepsenst,caux,form,0)
         !close(io_unit_pdfepsenst)
       end if

   end if

return
end

subroutine fopen1 (lu,fn,fmt,flag)
   
   ! to open file named fn as logical unit lu,
   ! with fmt as formatted/unformatted option
   ! flag=1, logical unit value is returned by program
   ! other values of flag, logical unit value prvided by user is used
   
   use mpicom
   implicit none

   character(len=*) :: fn
   character(len=40) :: nam
   character(len=11) :: fmt
   logical :: opn
   integer :: lu,flag
   integer :: jj,jjnam,ios
   data opn/ .false. /

   ! flag=1 : assign available logical unit to file
   if(flag.eq.1) then

      call blanks (fn,jj)

      open (file=fn(1:jj), form=fmt, iostat=ios, newunit=lu)

   else

      ! check if the targetted logical unit is already opened   
      inquire (lu, opened=opn, name=nam)
      ! if the associated internal filename is the same as the argument fn,
      ! then just rewind the file
   
      ! take out any embedded blanks
      call blanks (fn,jj)
      call blanks (nam,jjnam)
   
      if (opn.and.nam(1:jjnam).eq.fn(1:jj)) then
         rewind (lu)
         return
      end if

      ! gives a warning if the associated internal filename is not
      ! the same as the argument fn
      if (opn.and.nam(1:jjnam).ne.fn(1:jj)) then
         write (6,*) 'warning: logical unit',lu,'is already connected'
         write (6,*) 'internal filename was   ',nam
      end if

      open (lu, file=fn(1:jj), form=fmt, iostat=ios,status='replace')

   end if

   if (taskid.le.1)  &
      write (6,*) 'taskid,lu,fopen1: filename=',taskid,lu,fn(1:jj)
   
   if (ios.gt.0) then
      write (6,*) 'error condition in opening logical unit',lu
      write (6,900) fn(1:jj),fmt
      write (6,*) 'stopped in fopen1'
      call MPI_ABORT (MPI_COMM_WORLD,ierr)
   end if

900 format (' fn, fmt=',a10,a14)

return
end
