program MPI_stats
   ! gfortran -O3 MPI_stats.F90 -o MPI_stats.x

   implicit none

   integer :: nsteps, numtasks, taskid, iproc, jproc, ngroups, groupsize
   integer :: i, j, k, nchar, ivar, istep, igroup
   character(len=30) :: mpi_fn_name
   real(kind=4), allocatable :: rt_a2a_all(:,:,:,:,:), bw(:,:,:,:,:)
   real(kind=4), allocatable :: time_stats(:,:,:,:,:), bw_stats(:,:,:,:,:)
   real(kind=4) :: temp(4), msg_size(4)

   write(6,*) "Enter number of MPI tasks"
   read(5,*) numtasks
   write(6,*) "Enter iproc, jproc"
   read(5,*) iproc, jproc

   ! get number of steps and message size
   write(mpi_fn_name,"('A2A.g0')") 
   call blanks(mpi_fn_name, nchar)
   open(20,file=mpi_fn_name(1:nchar))
   read(20,"(7x,i4,24x,i5)") nsteps, groupsize
   ngroups=numtasks/groupsize
   read(20,"(26x,f12.2,8x,f12.2,8x,f12.2,8x,f12.2)") msg_size(1:4)
   close(20)
   write(6,*) "nsteps=",nsteps
   write(6,*) "msgsize=",msg_size(:)

   allocate(rt_a2a_all(numtasks,4,5,2,nsteps))
   allocate(bw(numtasks,4,5,2,nsteps))

   ! read the A2A timings
   do i=0,ngroups-1
      write(mpi_fn_name,"('A2A.g',i5)") i
      call blanks(mpi_fn_name,nchar)
      open(20,file=mpi_fn_name(1:nchar))
      read(20,"(1x)") ! skip
      read(20,"(1x)") ! skip
      read(20,"(1x)") ! skip
      do istep=1,nsteps
         do j=1,2
            do igroup=1,groupsize
               do ivar=1,5
                  read(20,"(59x,4e12.4,i5)") temp(:), taskid
                  rt_a2a_all(taskid+1,:,ivar,j,istep) = temp(:)
                  ! calculate bandwidths
                  ! kx-col
                  k = 1
                  if (rt_a2a_all(taskid+1,k,ivar,j,istep).ne.0.) then
                     bw(taskid+1,k,ivar,j,istep) = &
                        2*msg_size(k)*jproc/rt_a2a_all(taskid+1,k,ivar,j,istep)/1024./1024.
                  end if
                  ! kx-row
                  k = 2
                  if (rt_a2a_all(taskid+1,k,ivar,j,istep).ne.0.) then
                     bw(taskid+1,k,ivar,j,istep) = &
                        2*msg_size(k)*iproc/rt_a2a_all(taskid+1,k,ivar,j,istep)/1024./1024.
                  end if
                  ! xk-row
                  k = 3
                  if (rt_a2a_all(taskid+1,k,ivar,j,istep).ne.0.) then
                     bw(taskid+1,k,ivar,j,istep) = &
                        2*msg_size(k)*iproc/rt_a2a_all(taskid+1,k,ivar,j,istep)/1024./1024.
                  end if
                  ! xk-col
                  k = 4
                  if (rt_a2a_all(taskid+1,k,ivar,j,istep).ne.0.) then
                     bw(taskid+1,k,ivar,j,istep) = &
                        2*msg_size(k)*jproc/rt_a2a_all(taskid+1,k,ivar,j,istep)/1024./1024.
                  end if
               end do
            end do
         end do
      end do
      close(20)
   end do

   allocate(time_stats(3,4,5,2,nsteps))
   allocate(bw_stats(3,4,5,2,nsteps))
   ! calculate statistics
   do istep=1,nsteps
      do j=1,2
         do ivar=1,5
            do k=1,4
               ! minimum
               time_stats(1,k,ivar,j,istep) = minval(rt_a2a_all(:,k,ivar,j,istep))
               bw_stats(1,k,ivar,j,istep) = minval(bw(:,k,ivar,j,istep))
               ! maximum
               time_stats(3,k,ivar,j,istep) = maxval(rt_a2a_all(:,k,ivar,j,istep))
               bw_stats(3,k,ivar,j,istep) = maxval(bw(:,k,ivar,j,istep))
               ! average
               time_stats(2,k,ivar,j,istep) = sum(rt_a2a_all(:,k,ivar,j,istep))/numtasks
               bw_stats(2,k,ivar,j,istep) = sum(bw(:,k,ivar,j,istep))/numtasks
            end do
         end do
      end do
   end do

   ! write out time statistics
   open(21,file='kxcol.Time.stats')
   open(22,file='kxrow.Time.stats')
   open(23,file='xkrow.Time.stats')
   open(24,file='xkcol.Time.stats')
   do istep=1,nsteps
      do j=1,2
         do ivar=1,5
            write(21,"('istep=',i4,' kstep=',i2,' ivar=',i2,' kx-col (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, time_stats(:,1,ivar,j,istep)
            write(22,"('istep=',i4,' kstep=',i2,' ivar=',i2,' kx-row (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, time_stats(:,2,ivar,j,istep)
            write(23,"('istep=',i4,' kstep=',i2,' ivar=',i2,' xk-row (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, time_stats(:,3,ivar,j,istep)
            write(24,"('istep=',i4,' kstep=',i2,' ivar=',i2,' xk-col (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, time_stats(:,4,ivar,j,istep)
         end do
      end do
      write(21,"('------------------------------------------------------------------------------------------------')")
      write(22,"('------------------------------------------------------------------------------------------------')")
      write(23,"('------------------------------------------------------------------------------------------------')")
      write(24,"('------------------------------------------------------------------------------------------------')")
   end do
   close(21)
   close(22)
   close(23)
   close(24)

   ! write out bandwidth statistics
   open(21,file='kxcol.BW.stats')
   write(21,"('Bandwidth in GB/s     BW = 2*msgsize*jproc/T')")
   open(22,file='kxrow.BW.stats')
   write(22,"('Bandwidth in GB/s     BW = 2*msgsize*iproc/T')")
   open(23,file='xkrow.BW.stats')
   write(23,"('Bandwidth in GB/s     BW = 2*msgsize*iproc/T')")
   open(24,file='xkcol.BW.stats')
   write(24,"('Bandwidth in GB/s     BW = 2*msgsize*jproc/T')")
   do istep=1,nsteps
      write(21,"('------------------------------------------------------------------------------------------------')")
      write(22,"('------------------------------------------------------------------------------------------------')")
      write(23,"('------------------------------------------------------------------------------------------------')")
      write(24,"('------------------------------------------------------------------------------------------------')")
      do j=1,2
         do ivar=1,5
            write(21,"('istep=',i4,' kstep=',i2,' ivar=',i2,' kx-col (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, bw_stats(:,1,ivar,j,istep)
            write(22,"('istep=',i4,' kstep=',i2,' ivar=',i2,' kx-row (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, bw_stats(:,2,ivar,j,istep)
            write(23,"('istep=',i4,' kstep=',i2,' ivar=',i2,' xk-row (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, bw_stats(:,3,ivar,j,istep)
            write(24,"('istep=',i4,' kstep=',i2,' ivar=',i2,' xk-col (min.avg,max)',1p,3e12.3)") &
               istep, j, ivar, bw_stats(:,4,ivar,j,istep)
         end do
      end do
   end do
   close(21)
   close(22)
   close(23)
   close(24)

end program

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
