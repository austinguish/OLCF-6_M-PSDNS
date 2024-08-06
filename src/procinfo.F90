subroutine procinfo()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Collect name of node on which each task is running
! Collect associated device name (GPU)
! Collect CPU-ID on which threads are running
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use omp_lib
   use hipfort
   use mpicom
   implicit none

   interface
      function getcpuid() bind(C,name='sched_getcpu')
         integer :: getcpuid
      end function
   end interface

   character(len=MPI_MAX_PROCESSOR_NAME) :: nodename
   integer :: resultlen
   integer :: itask
   integer(kind=c_int) :: iDev
   integer :: ithr,st,en,stride
   integer, allocatable :: cpuid(:),cpuid_all(:,:)
   character(len=256) :: gpuname
   integer(kind=c_int), allocatable :: iDev_all(:)
   character(len=:), allocatable :: nodename_all(:)
   real(kind=4), allocatable :: gpumem(:,:),gpumem_all(:,:,:)
   integer(kind(hipSuccess)) :: hip_ierr

   !!!! Get name of node on which each task is running 
   call MPI_GET_PROCESSOR_NAME(nodename,resultlen,ierr)

   allocate(character(resultlen) :: nodename_all(numtasks))

   call MPI_GATHER(nodename,resultlen,MPI_CHAR,nodename_all,resultlen,&
      MPI_CHAR,0,MPI_COMM_WORLD,ierr)

   !!!! Get cpuid of each thread
   num_thr=0
   !$OMP PARALLEL
   
   !$ num_thr=OMP_GET_NUM_THREADS()

   !$OMP END PARALLEL

   allocate(cpuid(num_thr))
   allocate(cpuid_all(num_thr,numtasks))

   ithr=0
   cpuid(:)=0
   !$OMP PARALLEL PRIVATE(ithr)

   !$ ithr=OMP_GET_THREAD_NUM()
   cpuid(ithr+1)=getcpuid()

   !$OMP END PARALLEL

   call MPI_GATHER(cpuid,num_thr,MPI_INTEGER,cpuid_all,num_thr,&
      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   !!!! Get gpuid for each task
   hip_ierr=hipGetDevice_(iDev)
   !resultlen=verify(prop%name,' ',.true.)
   !gpuname(1:resultlen) = prop%name(1:resultlen)
   !write(6,*) taskid, gpuname(1:resultlen)

   allocate(iDev_all(numtasks))

   call MPI_GATHER(iDev,1,MPI_INTEGER,iDev_all,1,&
      MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

   !!!! Print to file procinfo
   if(taskid.eq.0) then
      open(61,file='procinfo')
      write(61,"('TASKID',2x,' ipid  jpid ',1x,'START',2x,'END',1x,2x,'STRIDE', &
     &   6x,'GPU-ID',8x,'NODE-ID')")
      do itask=1,numtasks
         st=cpuid_all(1,itask)
         en=cpuid_all(num_thr,itask)
         stride=0
         if(num_thr.gt.1) stride=(en-st)/(num_thr-1)
         write(61,"(i6,2x,4i6,2x,i6,2x,i8,a20)")itask-1,ipid_all(itask-1), &
            jpid_all(itask-1),st,en,stride,iDev_all(itask),nodename_all(itask)

         !!if (mod(itask-1,8).eq.0) write(6,"(i6,a20)")itask-1,nodename_all(itask)

      end do
      close(61)
   end if

return
end
