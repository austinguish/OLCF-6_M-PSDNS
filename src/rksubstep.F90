subroutine rksubstep

   use comp
   use timers_comm
   use timers_rkstep
   implicit none

   integer i,ii
   real(8) rtime0, rtime1, tt(10)

   ! A single R-K substep

   rt_a2a(:,:)=0.0
   rtime0=MPI_WTIME()

   if (kstep.eq.1) tt=0.

   ! Value of kstep is available through the modules

   ! apply phase shift on the first 3 variables (the velocity components)

   !$OMP TARGET DATA MAP (tofrom:u, un) MAP(alloc:sndbuf_single, rcvbuf_single, buf) &
   !$OMP MAP(to:b11, b22, b33, kx, ky, kz, bk1, bk3, bk1i, bk3i, &
   !$OMP sx, sy, sz, mask)
 
   rtime1=MPI_WTIME()

   call phshift_inplace (u,3,kstep)
 
   tt(1)=tt(1)+MPI_WTIME()-rtime1

   rtime1=MPI_WTIME()

   call itransform_vel (u,u,u)

   tt(2)=tt(2)+MPI_WTIME()-rtime1

   rtime1=MPI_WTIME()
        
   ! form nonlinear terms (and obtain velmax)
   call realspace_vel (u,2,kstep)
 
   tt(3)=tt(3)+MPI_WTIME()-rtime1

   rtime1=MPI_WTIME()

   ! transform back the partially formed convective terms to wavenumber space

   if (kstep.eq.1) then 
      call step
      if (entime.le.0.and.jstep.eq.nsteps) istop=1
      if (istop.eq.1) then
         ioflag=1
         if (isave.ge.0) isflag=1
      end if

      if (taskid.eq.0) write (6,"('istep,ioflag,isflag,istop=',i6,3i3)") istep,ioflag,isflag,istop
   end if

   rtime1=MPI_WTIME()

   call transform (u,u,u,2)

   tt(4)=tt(4)+MPI_WTIME()-rtime1

   ! complete the convective terms

   rtime1=MPI_WTIME()

   call wavespace_vel (u,2)
   call field_update_vel (u,un,3)

   tt(5)=tt(5)+MPI_WTIME()-rtime1

   !$OMP END TARGET DATA

   tt(6)=tt(6)+MPI_WTIME()-rtime0

#ifdef DETAIL_MPI
   if (mpi_detail.eq.1) then

      call MPI_GATHER (rt_a2a,20,MPI_REAL8,rt_a2a_grp,20,MPI_REAL8, &
                         0,mpi_comm_stat_within, ierr)

       if (ipid_stat.eq.0) then
          do ii=0,mpistat_groupsize-1
             do i=1,5
                write(io_mpi_grp,"(i5,3x,i5,3x,i4,2x,e12.4,2x,e12.4,2x,e12.4,2x,e12.4,2x,i7)") &
                      istep, kstep, i, rt_a2a_grp(:,i,ii),itask_grp(ii,jpid_stat)
             end do
          end do
       end if

   end if
#endif

   if (kstep.eq.2) then
      t_rks(:)=t_rks(:)+tt(:)
   end if

   return
end 
