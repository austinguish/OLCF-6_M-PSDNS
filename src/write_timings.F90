        subroutine write_timings 

        use comp
        use timers_comm
        use timers_comp
        use timers_io
        use timers_rkstep
        use timers_tran
        use timers_fom

        implicit none

        integer itask,i,ii,jj
        real(8) sum,comm2
        real(8) fom_timer
        real(8), allocatable :: t_rks_all(:,:),t_itrans_all(:,:)
        integer, allocatable :: numal_all(:)
        real(8), allocatable :: t_xkcomm1_all(:,:)
        real(8), allocatable :: t_xkcomm2_all(:,:)
        real(8), allocatable :: t_kxcomm1_all(:,:)
        real(8), allocatable :: t_kxcomm2_all(:,:)
        real(8), allocatable :: t_trans_all(:,:)
        real(8), allocatable :: t_pack_all(:,:)
        real(8), allocatable :: t_unpack_all(:,:)
!
        integer itask1,itask2
        character*20 string
        character*6 filepos
!
! detailed instrumentation timings, added by PKY, 2/3/2012
!
	if (taskid.eq.0)   & 
        write (6,*) 'enter write_timings: ncpusteps=', ncpusteps

	if (ncpusteps.eq.0) return
!
	allocate (t_rks_all(10,0:numtasks-1))
	allocate (t_itrans_all(4,0:numtasks-1))
	allocate (t_trans_all(4,0:numtasks-1))
	allocate (t_kxcomm1_all(4,0:numtasks-1))
	allocate (t_kxcomm2_all(4,0:numtasks-1))
	allocate (t_xkcomm1_all(4,0:numtasks-1))
	allocate (t_xkcomm2_all(4,0:numtasks-1))


	allocate (t_pack_all(4,0:numtasks-1))
	allocate (t_unpack_all(4,0:numtasks-1))

	t_rks(:)=t_rks(:)/ncpusteps
        call  MPI_ALLGATHER (t_rks(1),10,MPI_DOUBLE_PRECISION, &
                             t_rks_all,10,MPI_DOUBLE_PRECISION,  &
                             MPI_COMM_WORLD,mpierr)
	t_itrans(:)=t_itrans(:)/ncpusteps
        call  MPI_ALLGATHER (t_itrans,4,MPI_DOUBLE_PRECISION,    &
                             t_itrans_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)

	t_trans(:)=t_trans(:)/ncpusteps
        call  MPI_ALLGATHER (t_trans,4,MPI_DOUBLE_PRECISION,       &
                             t_trans_all,4,MPI_DOUBLE_PRECISION,   &
                             MPI_COMM_WORLD,mpierr)

	t_kxcomm1(:)=t_kxcomm1(:)/ncpusteps
        call  MPI_ALLGATHER (t_kxcomm1,4,MPI_DOUBLE_PRECISION,  &
                             t_kxcomm1_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
!
 	t_kxcomm2(:)=t_kxcomm2(:)/ncpusteps
        call  MPI_ALLGATHER (t_kxcomm2,4,MPI_DOUBLE_PRECISION, &
                             t_kxcomm2_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
!
	t_xkcomm1(:)=t_xkcomm1(:)/ncpusteps
        call  MPI_ALLGATHER (t_xkcomm1,4,MPI_DOUBLE_PRECISION, &
                             t_xkcomm1_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
!
	t_xkcomm2(:)=t_xkcomm2(:)/ncpusteps
        call  MPI_ALLGATHER (t_xkcomm2,4,MPI_DOUBLE_PRECISION, &
                             t_xkcomm2_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
!
	t_pack(:)=tot_pack(:)/ncpusteps
        call  MPI_ALLGATHER (t_pack,4,MPI_DOUBLE_PRECISION, &
                             t_pack_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
	t_unpack(:)=tot_unpack(:)/ncpusteps
        call  MPI_ALLGATHER (t_unpack,4,MPI_DOUBLE_PRECISION, &
                             t_unpack_all,4,MPI_DOUBLE_PRECISION, &
                             MPI_COMM_WORLD,mpierr)
!
	t_rks(:)=t_rks(:)*ncpusteps
	t_itrans(:)=t_itrans(:)*ncpusteps
	t_trans(:)=t_trans(:)*ncpusteps
	t_kxcomm1(:)=t_kxcomm1(:)*ncpusteps
	t_kxcomm2(:)=t_kxcomm2(:)*ncpusteps
	t_xkcomm1(:)=t_xkcomm1(:)*ncpusteps
	t_xkcomm2(:)=t_xkcomm2(:)*ncpusteps
!
	if (taskid.gt.0) go to 90
!

	itask1=2*iproc-1
	itask2=numtasks-2*iproc
	write (6,*) 'write_timings, itask1,itask2=',itask1,itask2

	filepos='rewind'
	write (6,*) 'write_timings, ichkpt=',ichkpt
	if (ichkpt.gt.1) filepos='append'

 	open (77,file='cpu_rk_all',position=filepos)
      write (77,711) nx,ny,nz,nc,iproc, jproc
 711  format ('Pencils code, nx,ny,nz,nc=',3i6,i3, &
              'iproc, jproc=',i4,' x ',i6)
        write (77,"('num_thr=',i3,'  single prec,  rkmethod=',i3)") num_thr, &
          rkmethod

!	write (77,"(a20)") string
 	write (77,"('gpumpi=',i3)") gpumpi
	write (77,"('Detailed breakdown of CPU costs per time step')")
	write (77,"('averaged over', i4,' RK2 time steps (exclusive of I/O)')") ncpusteps
	write (77,"('taskid ipid jpid           itransform realspace', &
                  & ' transform           overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
	sum=0.
	do i=1,5
	sum=sum+t_rks_all(i,itask)
	end do
        write (77,630) itask,ipid_all(itask),jpid_all(itask), &
                   (t_rks_all(i,itask),i=1,6)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
 630    format (i6,2i5,1x,1p,6e10.3)
	close (77)
!
	open (77,file='aftrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, after sub. transform, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid    step     wavespace   advanc     barrier    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask), &
                   (t_rks_all(i,itask),i=7,10),t_rks_all(5,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
	end if
	end if
	end do
 631    format (i6,2i5,1p,5e11.3)
	close (77)
!
	open (77,file='t_itrans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. itransform_vel, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	if (nc.gt.0) then
	write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    scalars    total')")
	else
	write (77,"('taskid ipid jpid   kxcomm1    fft(z)     kxcomm2    velgrad    total')")
	end if
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),  &
                   (t_itrans_all(i,itask),i=1,4),              &
                   t_itrans_all(1,itask)+t_itrans_all(2,itask) &
                   +t_itrans_all(3,itask)+t_itrans_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
!
	open (77,file='t_trans_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. transform_vel, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid   xkcomm1    fft(z)     convec     xkcomm2    total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask), &
                   (t_trans_all(i,itask),i=1,4),              &
                   t_trans_all(1,itask)+t_trans_all(2,itask) &
                   +t_trans_all(3,itask)+t_trans_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
!
	open (77,file='t_kxcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then          
        write (77,631) itask,ipid_all(itask),jpid_all(itask), &
                   (t_kxcomm1_all(i,itask),i=1,3),            &
                   t_kxcomm1_all(1,itask)+t_kxcomm1_all(2,itask) &
                   +t_kxcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
!
	open (77,file='t_kxcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. kxcomm2 per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total      overall')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask),      &
                   (t_kxcomm2_all(i,itask),i=1,3),                 &
                   t_kxcomm2_all(1,itask)+t_kxcomm2_all(2,itask)   &
                   +t_kxcomm2_all(3,itask),t_kxcomm2_all(4,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do

	close (77)
!
	open (77,file='t_xkcomm1_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm1_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall  comp_fft   comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask), &
                   (t_xkcomm1_all(i,itask),i=1,3),            &
                   t_xkcomm1_all(1,itask)+t_xkcomm1_all(2,itask) &
                   +t_xkcomm1_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	end do
	close (77)
!
	open (77,file='t_xkcomm2_all',position=filepos)
	write (77,"('Detailed breakdown of CPU costs, sub. xkcomm2_clean per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid  alltoall   comp_fft  comp_other   total')")
        do itask=0,numtasks-1
	if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,631) itask,ipid_all(itask),jpid_all(itask), &
                   (t_xkcomm2_all(i,itask),i=1,3),            &
                   t_xkcomm2_all(1,itask)+t_xkcomm2_all(2,itask) &
                   +t_xkcomm2_all(3,itask)
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
	end if
	end do
	close (77)

        fom_timer = 0.0d0
!
	open (77,file='t_comm_all',position=filepos)
	write (77,"('Detailed breakdown of alltoall communication costs, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")


        do itask=0,numtasks-1

	comm2=t_kxcomm2_all(1,itask)
        sum=t_kxcomm1_all(1,itask)+comm2+    &
                   t_xkcomm1_all(1,itask)+t_xkcomm2_all(1,itask)

        if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,632) itask,ipid_all(itask),jpid_all(itask),     &
                   t_kxcomm1_all(1,itask),comm2,                  &
                   t_xkcomm1_all(1,itask),t_xkcomm2_all(1,itask), &
      		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
        fom_timer = fom_timer + sum
	end do
	close (77)
 632    format (i6,2i5,2x,1p,5e10.3,0p,f7.1,'%')

        fom_comm_time = fom_timer/dble(numtasks)
!
        fom_timer = 0.0d0

	open (77,file='fft_pct_all',position=filepos)
	write (77,"('Detailed breakdown of 1-D FFT costs, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid    fft-Y     fft-Z     fft-X     fft-X     fft-Z     fft-Y       Total      %Code')")
        do itask=0,numtasks-1
	sum=t_kxcomm1_all(2,itask)+t_kxcomm2_all(2,itask) & 
	     +t_xkcomm1_all(2,itask)+t_xkcomm2_all(2,itask) &
	     +t_itrans_all(2,itask)+t_trans_all(2,itask)
        if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,"(i6,1x,1p,6e10.3,2x,1p,e10.3,2x,0p,2x,f4.1,' %')") itask, &
        t_kxcomm1_all(2,itask),t_itrans_all(2,itask),t_kxcomm2_all(2,itask), &
        t_xkcomm1_all(2,itask),t_trans_all(2,itask),t_xkcomm2_all(2,itask),  &
      		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
        fom_timer = fom_timer + sum
	end do
	close (77)

        fom_fft_time = fom_timer/dble(numtasks)

        !! accumulate the xk_pack and kx_pack in a single timer for reporting
        fom_timer = 0.0d0

	open (77,file='kx_pack_pct_all',position=filepos)
	write (77,"('Detailed breakdown of pack+unpack (kxcomm1,2)  per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid  pack-kx1   unpack   pack-kx2   unpack       Total      %Code')")
        do itask=0,numtasks-1
	sum=t_pack_all(1,itask)+t_pack_all(2,itask) & 
	   +t_unpack_all(1,itask)+t_unpack_all(2,itask) 
        if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,"(i6,1x,1p,4e10.3,2x,1p,e10.3,2x,0p,2x,f4.1,' %')") itask, &
        t_pack_all(1,itask),t_unpack_all(1,itask),t_pack_all(2,itask),t_unpack_all(2,itask), &
      		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
        fom_timer = fom_timer + sum
	end do
	close (77)
!
	open (77,file='xk_pack_pct_all',position=filepos)
	write (77,"('Detailed breakdown of pack+unpack (xkcomm1,2)  per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid  pack-xk1   unpack   pack-xk2   unpack       Total      %Code')")
        do itask=0,numtasks-1
	sum=t_pack_all(3,itask)+t_pack_all(4,itask) & 
	   +t_unpack_all(3,itask)+t_unpack_all(4,itask) 
        if (itask.le.itask1.or.itask.ge.itask2) then
        write (77,"(i6,1x,1p,4e10.3,2x,1p,e10.3,2x,0p,2x,f4.1,' %')") itask, &
        t_pack_all(3,itask),t_unpack_all(3,itask),t_pack_all(4,itask),t_unpack_all(4,itask), &
      		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
        fom_timer = fom_timer + sum
	end do
	close (77)

        fom_pack_time = fom_timer/dble(numtasks)

!
!
!
#ifdef NOT_YET
	open (77,file='t_loctsp_all',position=filepos)
	write (77,"('Detailed breakdown of local-transpose costs, per time step')")
	write (77,"('averaged over', i4,'  non-output steps')") ncpusteps
	write (77,"('taskid ipid jpid    kxcomm1   kxcomm2  xkcomm1   xkcomm2    total      %Code')")
	if (itask.le.itask1.or.itask.ge.itask2) then
	comm2=t_kxcomm2_all(3,itask)
        sum=t_kxcomm1_all(3,itask)+comm2+   &
                   t_xkcomm1_all(3,itask)+t_xkcomm2_all(3,itask)
        write (77,632) itask,ipid_all(itask),jpid_all(itask),     &
                   t_kxcomm1_all(3,itask),comm2,                  &
                   t_xkcomm1_all(3,itask),t_xkcomm2_all(3,itask), &
      		sum,sum/t_rks_all(6,itask)*100.
        if (mod(itask+1,iproc).eq.0) then
        write (77,"('--------------------------------------------')")
        end if
        end if
	close (77)
!
	end if

#endif

 90     continue

	deallocate (t_rks_all)
	deallocate (t_itrans_all)
	deallocate (t_trans_all)
	deallocate (t_kxcomm1_all)
	deallocate (t_xkcomm1_all)
	deallocate (t_xkcomm2_all)
	deallocate (t_kxcomm2_all)
	deallocate (t_pack_all)
	deallocate (t_unpack_all)

	return
	end
