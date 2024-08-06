!  3D FFT inverse transform with 2D domain decomposition
!
!  This version uses MPI_Alltoallv to exchange data between processors
!  In the second step, y planes are sent separately
!  The order of array elements in memory is the same in all stages: (x,z,y) 

! Input: XZYg - comlpex array, with y dimension contained entirely,
!               while x and z are block-distributed among processors in 2D grid
! XZgY - complex array (auxiliary), z dimension contained entirely
!               while x and y are block-distributed among processors in 2D grid
! Output: XgZY - an array of real, x dimension is contained entirely within processors memory  
!               while z and y are block-distributed among processors in 2D grid

! !!! CAUTION: In this version: all arrays occupy the same memory space
!

      subroutine kxtran_gpu (XgZY,XZgY,XZYg,nv)
      use com
      use timers_comm
      use timers_comp
      use timers_tran
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none
!

 
      real(b8) :: XgZY(nxpad,zisz,yjsz,nv)
      complex(b8), target :: XZgY(xisz,nz,yjsz,nv)
      complex(b8) :: XZYg(xisz,ny,zjsz,nv)

      real(b8) factor
      integer :: i,nv
      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip
 
      real*8 rtime1,rtime2,rtime3,rtime0

      rtime0=MPI_WTIME()
      rtime1=MPI_WTIME()

      ! Transform in y dimension for all x and z
      !
      ! transform in y direction and switch into z-lines
 
      !$OMP TARGET DATA MAP (tofrom:XgZY,XZgY,XZYg)

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
      !$OMP PRIVATE(zp,xp,i) SHARED(XZYg,nyhp,nv,xisz,zjsz)
      do zp=1,zjsz
         do xp=1,xisz
            do i=1,nv
               XZYg(xp,nyhp,zp,i)=cmplx(0.,0.)
            end do
         end do
       end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      do i=1,nv
         call kxcomm1_gpu (XZYg(1,1,1,i),XZgY(1,1,1,i),1)
      end do
      rtime1=MPI_WTIME()-rtime1
      t_kx(1)=t_kx(1)+rtime1
      t_kxcomm1(4)=t_kxcomm1(4)+rtime1


      ! Transform in z dimension for all x, one y-plane at a time

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
      !$OMP DEFAULT(NONE) PRIVATE(i,x,y) SHARED(nv,xisz,yjsz,nzhp,XZgY)
      do i=1,nv
         do y=1,yjsz
            do x=1,xisz
               XZgY(x,nzhp,y,i)=cmplx(0.,0.)
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
 
      rtime2=MPI_WTIME()

      !$OMP TARGET DATA USE_DEVICE_PTR(XZgY)
      do i=1,nv
         do y=1,yjsz

#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecZ2Z(hip_plan_z, c_loc(XZgY(1,1,y,i)), &
                                     c_loc(XZgY(1,1,y,i)), HIPFFT_BACKWARD)
#else
            ierr_hip = hipfftExecC2C(hip_plan_z, c_loc(XZgY(1,1,y,i)), &
                                     c_loc(XZgY(1,1,y,i)), HIPFFT_BACKWARD)
#endif

         end do
      end do
      !$OMP END TARGET DATA

      rtime2=MPI_WTIME()-rtime2
      tcpu_fft=tcpu_fft+rtime2
      t_kx(2)=t_kx(2)+rtime2

 
      ! switch into x-lines and take complex-real transform in x	
      rtime3=MPI_WTIME()
       do i=1,nv
         call kxcomm2_gpu (XZgY(1,1,1,i),XgZY(1,1,1,i),1)
      end do

      !$OMP END TARGET DATA
      rtime3=MPI_WTIME()-rtime3
      rtime0=MPI_WTIME()-rtime0
 
      t_kx(3)=t_kx(3)+rtime3
      t_kxcomm2(4)=t_kxcomm2(4)+rtime3
      t_kx(4)=t_kx(4)+rtime0
 
 
      return
      end
