! 'New' 3D DNS code on GPUs
! 
!  2D (pencils) decomposition, stride 1
!  mpi_alltoall, OpenMP implemented in pfield.f only
!
!  Intentionally simplified from prior versions, relative to past
!  directives:
!
!  --- No cylindrical truncation, no balance period, etc.
!  

!
program DNS_PENCIL

   use omp_lib
   use comp

   use timers_comm
   use timers_comp
   use timers_tran
   use timers_rkstep
   use timers_io
   use timers_fom
   use ranseed

   implicit none

   integer dmy(3),hms(3),date_time(8)
   integer i,n,j,m,iaux,nv
   real twopi
   integer idim(4),idimc(2)
   real sum1,raux
   real(8) rtime1,rtime2,difftime
   real acccpu,totcpu,tcpu,cpumin,cpumax,avcpu,acccpu2
   real(8) tt_comm,gtt_comm(2)

   real, allocatable :: tmp1(:),tmp2(:)
   
   real*8 totcpu_fft,totcpu_other,rtime0
   real*8, allocatable :: cpucomm(:),cpucomp(:,:),cputot(:)
   
   real*8 rtime3
   real(8) tmin,tav,tmax

   integer itask,ii
   integer prov, idev

   integer disp_int

   real(b8) msqvel
   real(kind=b8) :: msg1, msg2
   character(len=30) :: mpi_fn_name
   real(8) fom_timer
   integer fom_iters

   integer sizebs,sizebsg, maxsizebs,maxsizebsg

   !! local values that are accumulated each time step
   fom_iters = 0
   fom_timer = 0.0d0

   !! values in module timers_fom
   fom_avg_time = 0.0d0
   fom_problem_size = 0.0d0
   computed_fom = 0.0d0
   fom_fft_time = 0.0d0
   fom_pack_time = 0.0d0
   fom_comm_time = 0.0d0
   fom_other_time = 0.0d0
   ifom = 1



#ifdef OPENMP
   prov=0
   call MPI_INIT_THREAD (MPI_THREAD_FUNNELED,prov,ierr)
#else
   call MPI_INIT (ierr)
#endif
   call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
   call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)

   call date_and_time (values=date_time)
   dmy(1:3)=date_time(3:1:-1)
   hms(1:3)=date_time(5:7)
   if (taskid.eq.0.or.taskid.eq.numtasks-1) then
      write (6,601) taskid,dmy(2),dmy(1),dmy(3),hms
601   format (' starting taskid=',i6, ' date & time is  ',i2,'/',i2, &
              '/',i4,2x,i2,':',i2,':',i2)
   end if
!
   one = 1
   zero = 0
   twopi=2.*pi

   call input
   b11=1.; b22=1. ; b33=1. 
  
   ! Initializing
   call param_set

   if (taskid.eq.0) then
      write (6,"('nx,ny,nz,nc,nu,numtasks=',3i6,2i3,i6)") &
         nx,ny,nz,nc,nu,numtasks
   end if

   rtime1=MPI_WTIME()  ! time initialization tasks
   rwall0=rtime1


   call mpisetup
   call CPUmeminfo ('after mpisetup',14)
   call GPUmeminfo ('before com_set',14)
 
   rkstages=2

   call com_set
   call CPUmeminfo ('after com_set',13)

   call custom_flush(6)

   call comp_set
   call CPUmeminfo ('after comp_set',14)

   call custom_flush(6)

   call comsp_set
   call CPUmeminfo ('after comsp_set',15)
   if (taskid.eq.0) write(6,*) "After comsp_set"

   call custom_flush(6)

   call openf
   call procinfo
   
   call GPUmeminfo ('after procinfo',14)

   if(taskid .eq. 0) then
      print *,'Using processor grid ',iproc,' x ',jproc
   end if
 
!ccccccccc for test purposes
   istep0=0
   istop=0
!ccccccccccccccccccccccccccc	

   tcpu_fft=0.
   tcpu_other=0.
        
   call masks

   a11=0.
   a22=0.
   a33=0.

   call GPUmeminfo ('before epfftw',13)

   !! create the fft plans
   call epfftw (u,u,u,u)

   call CPUmeminfo ('after epfftw',12)


   call random_seed (size=seed_size)
   allocate (rseed(seed_size))

   call incond ()
   if (kinit(1).eq.-1) call sptvar_omp (un)

   call custom_flush(6)
!
   if (taskid.eq.0) then
      open (7,file='log')
      write (7,711) nx,ny,nz
711   format ('DNS-PEN-GPU turbulence code, nx,ny,nz=',3i6)
      write (7,712) numtasks,iproc,jproc
712   format ('No. of MPI tasks=',i7, '   iproc, proc=',i4,' x ',i4)
      write (7,"('pencils (no cyl trunc)')")
   end if


#ifdef DETAIL_MPI
   call mpistat_setup

   if (mpi_detail.eq.1) then

      ! Switch to 1 file per group, instead of 1 file per MPI process
      if (ipid_stat.eq.0) then
         write(mpi_fn_name,"('MPI_timings/A2A.g',i5)") jpid_stat
         call blanks(mpi_fn_name, i)
         !write(6,*) "taskid=",taskid,"opening file",mpi_fn_name(1:i)
         open(file=mpi_fn_name(1:i), newunit=io_mpi_grp)
         write(io_mpi_grp,"('nsteps=',i4,' taskid=',i5,' groupsize=',i5)") nsteps, taskid, mpistat_groupsize
         msg1=real(xisz,b8)*yjsz*zjsz*8./1024.
         msg2=real(xisz,b8)*yjsz*zisz*8./1024.
         write(io_mpi_grp,"('message size (KB): kx-col=',f12.2,' kx-row=',f12.2,' xk-row=',f12.2,' xk-col=',f12.2)") &
            msg1, msg2, msg2, msg1
         write(io_mpi_grp,"('------------------------------------------------------------------------------------------------')")
         write(io_mpi_grp,*) " " 
         write(io_mpi_grp,"('istep   kstep   ivar        kx-col        kx-row        xk-row        xk-col     rank')")
         write(io_mpi_grp,"('-------------------------------------------------------------------------------------')")
      end if

   end if
#endif

   if (taskid.eq.0) then
      write (6,"('shape(un)          ',5i6)") shape(un)
      write (6,"('shape(u)=          ',5i6)") shape(u)
      write (6,"('shape(sndbuf_single) =',5i6)") shape(sndbuf_single)
      write (6,"('shape(rcvbuf_single) =',5i6)") shape(rcvbuf_single)
      write (6,"('shape(buf)=        ',5i6)") shape(buf)
      write (6,"('shape(mask)=       ',5i6)") shape(mask)
   end if

   call CPUmeminfo ('before OMP in main.F90',22)
   call GPUmeminfo ('before OMP in main.F90',22)

   !$OMP TARGET DATA MAP (tofrom:u, un) MAP(alloc:sndbuf_single, rcvbuf_single, buf) &
   !$OMP MAP(to:b11, b22, b33, kx, ky, kz, bk1, bk3, bk1i, bk3i, &
   !$OMP sx, sy, sz, mask)

   call GPUmeminfo ('after OMP in main.F90',21)

   ! transform 'pfield' initial conditions to wavenumber space
   if (kinit(1).eq.-2) then

      do i=1,3+nc
         call xktran_gpu (un(1,1,1,i),un(1,1,1,i),un(1,1,1,i),1)
      end do

      !$OMP TARGET UPDATE FROM (un)
      call sptvar_omp (un)

   end if

   call sptr_omp (un)
   if (taskid.eq.0) write(6,*) "After sptr"

   call custom_flush(6)

   call MPI_BARRIER (MPI_COMM_WORLD, mpierr)

   ! update u
   u(:,:,:,1)=un(:,:,:,1)
   u(:,:,:,2)=un(:,:,:,2)
   u(:,:,:,3)=un(:,:,:,3)

   !$OMP TARGET UPDATE TO (u)

   nv = 1
   if (taskid.eq.0) write (6,*) 'main: MULTVAR directive is Not Available, nv=',nv

   rtime0=MPI_WTIME()
   rtime1=MPI_WTIME()

   ncpusteps=0
   t_rks=0.
   t_xkcomm1=0 ; t_xkcomm2=0.  ; t_trans=0.
   t_kxcomm1=0 ; t_kxcomm2=0.  ; t_itrans=0.

   tot_pack=0. ; tot_unpack=0.

   call time_stamp ('before time-stepping  ')

   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   
   istep0=0
   istop=0

   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
   
   !!! start the timesteps ...

11 istep=istep+1
   jstep=istep-istep0

   call date_and_time (values=date_time)
   dmy(1:3)=date_time(3:1:-1)
   hms(1:3)=date_time(5:7)
   if (taskid.eq.0) then
      write (6,15) istep,taskid,dmy(2),dmy(1),dmy(3),hms
      write (0,15) istep,taskid,dmy(2),dmy(1),dmy(3),hms
      write (7,15) istep,taskid,dmy(2),dmy(1),dmy(3),hms
   end if
15 format ('istep=',i6, ' taskid',i6, ' date & time is  ',& 
        i2,'/',i2, '/',i4,2x,i2,':',i2,':',i2)

   rtime1=MPI_WTIME()

   call shifts ()
   !$OMP TARGET UPDATE TO (sx,sy,sz)

   t_pack=0. ; t_unpack=0.

   ncpusteps=ncpusteps+1

   ifom = 1
     
   do kstep=1,2

      call rksubstep ()

   end do

   tot_pack(:)=tot_pack(:)+t_pack(:)
   tot_unpack(:)=tot_unpack(:)+t_unpack(:)
 

   rtime2=MPI_WTIME()
   tcpu=rtime2-rtime1
   call MPI_REDUCE (tcpu,totcpu,1,mpireal,MPI_SUM,0, &
         MPI_COMM_WORLD,ierr)
   call MPI_REDUCE (tcpu,cpumin,1,mpireal,MPI_MIN,0, &
         MPI_COMM_WORLD,ierr)
   call MPI_REDUCE (tcpu,cpumax,1,mpireal,MPI_MAX,0, &
         MPI_COMM_WORLD,ierr)

   if (taskid.eq.0) then
      avcpu=totcpu/real(numtasks)

      !! for FOM:
      !!     exclude timesteps that collect data, compute statistics, 
      !!     and/or report 
      if (ifom .eq. 1) then
         fom_iters = fom_iters + 1
         fom_timer = fom_timer + dble(avcpu)
      endif

      write (7,20) istep,cpumin,avcpu,cpumax
20    format ('cpu at step',i7,' :min/ave/max=',3f11.3, '  secs.')
   end if

   !! do this here to avoid corrupting the runtime timer with 
   !! the time to compute the statistics
   if (ioflag.eq.1) then
      !$OMP TARGET UPDATE FROM (un)
      call sptr_omp (un)
   endif


   if (istop.eq.0) go to 11

   !$OMP END TARGET DATA

   t_alltoall=t_alltoall/nsteps
   tcpu_fft=tcpu_fft/nsteps
   tcpu_other=tcpu_other/nsteps


   call write_timings

   rtime0=(MPI_WTIME()-rtime0)/nsteps

   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)



90 call date_and_time (values=date_time)
   dmy(1:3)=date_time(3:1:-1)
   hms(1:3)=date_time(5:7)
   if (taskid.eq.0) then
      write (6,25) taskid,dmy(2),dmy(1),dmy(3),hms
25    format (' DONE, taskid=',i5, ' date & time is  ',i2,'/',i2, &
              '/',i4,2x,i2,':',i2,':',i2)
   end if

   if (taskid.eq.0) then
      rwall0 = MPI_WTIME() - rwall0
      write (7,30) rwall0
30    format ('overall wall time for task 0 =',f8.1, ' secs.')
   end if

   if (taskid.eq.0) close(7)

   !! compute and report the FOM
   if (taskid.eq.0) then
      fom_problem_size = dble(nx)*dble(ny)*dble(nz)
      fom_avg_time = fom_timer/(dble(fom_iters))
      computed_fom = fom_problem_size/fom_avg_time

      fom_other_time = fom_avg_time - fom_fft_time - fom_pack_time - fom_comm_time

      open (20,file='FOM')

      write(20,40) fom_iters, numtasks
40    format ('Times averaged over ',i5,' compute-only timesteps and over ',i7,' mpi ranks')
      write(20,*) " " 

      write(20,"('    FFT (s)     Pack and Unpack (s)       MPI (s)       Other (s)       Total (s)')")
      write(20,"('----------------------------------------------------------------------------------')")

      write(20,43) fom_fft_time, fom_pack_time, fom_comm_time, fom_other_time, fom_avg_time
43    format (F10.4,8x,F10.4,10x,F10.4,5x,F10.4,6x,F10.4)

      write(20,"('----------------------------------------------------------------------------------')")
      write(20,*) " "

      write(20,45) fom_avg_time
45    format ('Average runtime per timestep = ',F10.4)
      write(20,50) fom_problem_size
50    format ('Problem size = nx*ny*nz = ',F23.2)
      write(20,55) computed_fom
55    format ('FOM = (problem size)/(average runtime) = ',ES10.4)

      close(20)
   end if

   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

   call MPI_FINALIZE (ierr)

end program
!=========================================================
