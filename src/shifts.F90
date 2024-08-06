subroutine shifts
! (adapted from Kiran's CUDA code)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kshift=1 - random shifts, in which the aliasing error on the
!            predictor step approximately cancels that of the previous
!            corrector step (rogallo).  over a fixed time period,
!            the residual aliasing error in u is of order dt**2
!            but over a single step it is of order dt.  not suitable
!            for extraction of multi-time statistics.
! kshift=2 - random shifts, in which the aliasing error on the
!            corrector step cancels that of the previous
!            predictor step.  over a fixed time period,
!            the residual aliasing error in u is of order dt,
!            but over a single step it is of order dt**2.
! kshift=3 - fixed shifts.  error on any two adjacent steps cancels.
!            over a fixed time period, and over a single step, the
!            residual aliasing error in u is of order dt**2, but
!            but the residual error is not random.
! kshift=4 - zero shifts: order one aliasing errors. (for frozen fields)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   use com
   implicit none

   integer :: i,k
   complex(kind=b8) :: phase,argc
   real(kind=b8) :: ranu,cancel,dum

   integer, parameter :: nk=20
   common/seqcom/iseq(nk),klast
   integer :: iseq,klast
   integer :: iseqpdf(nk)

   !if(taskid.eq.0) write(6,*) "Enter shifts"

   if (taskid.eq.0) then

      if (istart.le.0.and.istep.eq.istep0+1.and.&
         (kshift.eq.1.or.kshift.eq.3)) then
         write (6,*) 'shifts: ksran=',ksran
         call ranseq (ksran,0)
         do i = 1,3
            gsh(i,2) = ranu(ksran)
         end do

         ! For RK4
         do i = 1,3
            gsh(i,4) = ranu(ksran)
         end do

      end if


      ! generate new random numbers for grid shifts if desired
      if (kshift.eq.1.or.kshift.eq.2) call ranseq (ksran,0)

      do i=1,3

         if(kshift.eq.1) then
            gsh(i,1) = cancel( gsh(i,2) )
            gsh(i,2) = ranu(ksran)
         else if(kshift.eq.2) then
            gsh(i,1) = ranu(ksran)
            gsh(i,2) = cancel( gsh(i,1) )
         else if(kshift.eq.3) then
            gsh(i,1) = cancel( gsh(i,2) )
            gsh(i,2) = cancel( gsh(i,1) )
         else if(kshift.eq.4) then
            gsh(i,1) = real(0.,b8)
            gsh(i,2) = real(0.,b8)
         end if

      end do


      if (rkstages.eq.4) then

         do i=1,3

            if(kshift.eq.1) then
               gsh(i,3) = cancel( gsh(i,4) )
               gsh(i,4) = ranu(ksran)
   
            else if(kshift.eq.2) then
               gsh(i,3) = ranu(ksran)
               gsh(i,4) = cancel( gsh(i,3) )
   
            else if(kshift.eq.3) then
               gsh(i,3) = cancel( gsh(i,4) )
               gsh(i,4) = cancel( gsh(i,3) )
   
            else if(kshift.eq.4) then
               gsh(i,3) = real(0.,b8)
               gsh(i,4) = real(0.,b8)
            end if
   
         end do

      end if

   end if

   ! Broadcast the grid shifts to all tasks
   call MPI_BCAST (gsh,3*4,mpireal,0,MPI_COMM_WORLD,mpierr)

   !if (taskid.eq.1) then
   !write (91,"('istep,gsh(1,1),gsh(1,2)=',i5,1p,2e12.4)")  &
   !            istep,gsh(1,1),gsh(1,2)
   !end if

   ! compute phase shift factors from grid shifts
   !if (.not.allocated(sx)) then
   !   allocate(sx(nxh,rkstages),sy(ny,rkstages),sz(nz,rkstages))
   !end if

   ! x shifts
   phase = imagi*2*pi/nx
   do k=1,rkstages
      do x=1,nxh
         argc = phase*kx(x)*gsh(1,k)
         sx(x,k) = exp( argc )
      end do
   end do

   ! y shifts
   phase = imagi*2*pi/ny
   do k=1,rkstages
      do y=1,ny
         argc = phase*ky(y)*gsh(2,k)
         sy(y,k) = exp( argc )
      end do
   end do

   ! z shifts
   phase = imagi*2*pi/nz
   do k=1,rkstages
      do z=1,nz
         argc = phase*kz(z)*gsh(3,k)
         sz(z,k) = exp( argc )
      end do
   end do

   !!!! Copy phase shifts to device
   !do iDev=0,nDevices-1
   !   ierr = cudaSetDevice (iDev)
   !   if (.not.allocated(iDevice(iDev)%d_sx)) then
   !      allocate (iDevice(iDev)%d_sx(nxh,rkstages), &
   !                iDevice(iDev)%d_sy(ny,rkstages), iDevice(iDev)%d_sz(mz,rkstages))
   !   end if
   !   iDevice(iDev)%d_sx(:,:) = sx(:,:)
   !   iDevice(iDev)%d_sy(:,:) = sy(:,:)
   !   iDevice(iDev)%d_sz(:,:) = sz(:,:)
   !end do
   !ierr = cudaSetDevice (0)
 
return
end
