!
! inverse transform of velocity from wavenumber space to physical space
! 3D FFT inverse transform with 2D domain decomposition

      subroutine itransform_vel (ux,uy,uz)

      use comp
      use timers_comm
      use timers_comp
      use timers_tran
      use timers_rkstep
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none

      real(b8) :: ux(nx,zisz,yjsz,nu)
      complex(b8), target :: uz(xisz,nz,yjsz,nu)
      complex(b8) :: uy(xisz,ny,zjsz,nu)

      real(b8) factor
      integer :: i,nv
      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip
 
      real*8 rtime1,rtime2,rtime3,rtime0

      rtime0=MPI_WTIME()
      rtime1=MPI_WTIME()

      if (istep.eq.1.and.kstep.eq.1) t_itrans=0.

      ! transform in y direction and switch into z-lines
 
      !$OMP TARGET DATA MAP (tofrom:ux,uz,uy)

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
      !$OMP PRIVATE(zp,xp,i) SHARED(uy,nyhp,nv,xisz,zjsz)
      do i=1,3
         do zp=1,zjsz
            do xp=1,xisz
               uy(xp,nyhp,zp,i)=cmplx(0.,0.)
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      do i=1,3
         call kxcomm1_gpu (uy(1,1,1,i),uz(1,1,1,i),1)
         rt_a2a(1,i)=t_a2a
      end do

      rtime1=MPI_WTIME()-rtime1
      t_itrans(1)=t_itrans(1)+rtime1
      t_kx(1)=t_kx(1)+rtime1
      t_kxcomm1(4)=t_kxcomm1(4)+rtime1

      call MPI_BARRIER (MPI_COMM_WORLD,mpierr)

      ! Transform in z dimension
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
      !$OMP DEFAULT(NONE) PRIVATE(i,x,y) SHARED(nv,xisz,yjsz,nzhp,uz)
      do i=1,3
         do y=1,yjsz
            do x=1,xisz
               uz(x,nzhp,y,i)=cmplx(0.,0.)
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
 

      rtime2=MPI_WTIME()

      !$OMP TARGET DATA USE_DEVICE_PTR(uz)
      do i=1,3
         do y=1,yjsz

#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecZ2Z(hip_plan_z, c_loc(uz(1,1,y,i)), &
                                     c_loc(uz(1,1,y,i)), HIPFFT_BACKWARD)
#else
            ierr_hip = hipfftExecC2C(hip_plan_z, c_loc(uz(1,1,y,i)), &
                                     c_loc(uz(1,1,y,i)), HIPFFT_BACKWARD)
#endif

         end do
      end do
      !$OMP END TARGET DATA

      ierr = hipDeviceSynchronize()

      call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
      rtime2=MPI_WTIME()-rtime2

      t_itrans(2)=t_itrans(2)+rtime2
      tcpu_fft=tcpu_fft+rtime2
      t_kx(2)=t_kx(2)+rtime2

 
      ! switch into x-lines and take complex-real transform in x	
 
      rtime3=MPI_WTIME()

      do i=1,3
         call kxcomm2_gpu (uz(1,1,1,i),ux(1,1,1,i),1)
         rt_a2a(2,i)=t_a2a
      end do

      rtime3=MPI_WTIME()-rtime3
      t_itrans(3)=t_itrans(3)+rtime3

      !$OMP END TARGET DATA

      rtime0=MPI_WTIME()-rtime0
 
      t_kx(3)=t_kx(3)+rtime3
      t_kxcomm2(4)=t_kxcomm2(4)+rtime3
      t_kx(4)=t_kx(4)+rtime0
 
      return
      end
