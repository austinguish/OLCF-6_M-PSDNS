      subroutine kxcomm1_gpu (source,dest,nv)

      ! Start with y-pencils (y,z,x)
      ! Perform forward inverse transform in y, then transpose to z pencils

      ! internally, one variable at a time


      use com
      use timers_comm
      use timers_comp
      use timers_tran
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none

      integer iv,nv
      ! STRIDED VERSION
      complex(b8), target :: source(xisz,nypad,zjsz,nv)
      complex(b8) :: dest(xisz,nzpad,yjsz,nv)

       real(b8) tcpu_pack,tcpu_unpack

      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

      integer i,j,n
      real*8 rtime1,rtime2,rtime3,rtime4

      t3_comm = t3_comm - MPI_Wtime() 

      ! This case is for using MPI_Alltoall 
      ! Assume nzpad=nz

      !$OMP TARGET DATA MAP(to:source) MAP(from:dest) &
      !$OMP MAP(alloc:sndbuf_col,rcvbuf_col)

      tcpu_pack=0.
      tcpu_unpack=0.

      do iv=1,nv

         ! zero the nyhp Fourier mode (for consistency with conjugate symmetry)

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
         !$OMP DEFAULT(NONE) PRIVATE(z,x) SHARED(zjsz,xisz,nyhp,iv,source)
         do z=1,zjsz
            do x=1,xisz
               source (x,nyhp,z,iv)=cmplx(0.,0.)
            end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         ! take inverse C-C transform in y

         rtime2=MPI_WTIME()

         !$OMP TARGET DATA USE_DEVICE_PTR(source)
         do zp=1,zjsz

#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecZ2Z(hip_plan_y, c_loc(source(1,1,zp,iv)), &
                                     c_loc(source(1,1,zp,iv)), HIPFFT_BACKWARD)
#else
            ierr_hip = hipfftExecC2C(hip_plan_y, c_loc(source(1,1,zp,iv)), &
                                     c_loc(source(1,1,zp,iv)), HIPFFT_BACKWARD)
#endif

         end do
         !$OMP END TARGET DATA

         !! sync device after the loop ... fewer sync calls
         ierr = hipDeviceSynchronize()

         rtime2=MPI_WTIME()-rtime2
         tcpu_fft=tcpu_fft+rtime2
         t_kxcomm1(2)=t_kxcomm1(2)+rtime2


         rtime1=MPI_WTIME()

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,xp,yp,zp,y) &
         !$OMP SHARED(jproc,yjsz,xisz,zjsz,iv, &
         !$OMP jjst,sndbuf_col,source)
         do i=1,jproc
            do zp=1,zjsz
               do yp=1,yjsz
                  do xp=1,xisz
                     y=jjst(i-1)+yp-1
                     sndbuf_col(xp,yp,zp,iv,i) = source(xp,y,zp,iv)
                  end do
               end do
            end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         if (gpumpi.eq.0) then
            !$OMP TARGET UPDATE FROM(sndbuf_col)
         end if

         rtime1=MPI_WTIME()-rtime1
         tcpu_pack=tcpu_pack+rtime1
         tcpu_other=tcpu_other+rtime1
         t_kxcomm1(3)=t_kxcomm1(3)+rtime1

      end do ! do nv

      ! Exchange the data (send the buffers)

#ifdef DETAIL_MPI
      if (mpi_detail.eq.1) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
#endif

      t_alltoall = t_alltoall - MPI_Wtime()
      rtime1=MPI_WTIME()
      if (gpumpi.eq.0) then
         call mpi_alltoall (sndbuf_col, xisz*yjsz*zjsz*nv, mpicomplex, &
                        rcvbuf_col, xisz*yjsz*zjsz*nv, mpicomplex, &
                     mpi_comm_col, ierr)
      else
         !$OMP TARGET DATA USE_DEVICE_PTR(sndbuf_col, rcvbuf_col)
         call mpi_alltoall (sndbuf_col, xisz*yjsz*zjsz*nv, mpicomplex, &
                        rcvbuf_col, xisz*yjsz*zjsz*nv, mpicomplex, &
                     mpi_comm_col, ierr)
         !$OMP END TARGET DATA
      end if

      rtime2=MPI_WTIME()
      t_a2a=rtime2-rtime1
      t_alltoall = t_alltoall + MPI_Wtime()
      t_kxcomm1(1)=t_kxcomm1(1)+rtime2-rtime1

!
      rtime1=MPI_WTIME()

      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE TO(rcvbuf_col)
      end if

      do iv=1,nv

         ! Unpack the receive buffer. Set z=nzhp mode to zero

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,yp,xp,zp,z) &
         !$OMP SHARED(jproc,yjsz,xisz,zjsz,iv, &
         !$OMP kjst,rcvbuf_col,dest)
         do i=1,jproc
            do yp=1,yjsz
               do zp=1,zjsz
                  do xp=1,xisz
                     z=kjst(i-1)+zp-1
                     dest(xp,z,yp,iv) = rcvbuf_col(xp,yp,zp,iv,i)
                  end do
               end do
            end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
         !$OMP DEFAULT(NONE) PRIVATE(yp,xp) SHARED(yjsz,xisz,nzhp, &
         !$OMP iv,dest)
         do yp=1,yjsz
            do xp=1,xisz
                  dest(xp,nzhp,yp,iv) = 0.0
            end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         rtime1=MPI_WTIME()-rtime1
         tcpu_unpack=tcpu_unpack+rtime1
         tcpu_other=tcpu_other+rtime1
         t_kxcomm1(3)=t_kxcomm1(3)+rtime1

      end do ! do nv
      

      t3_comm = t3_comm + MPI_Wtime() 
      i3_comm = i3_comm + 1

      !$OMP END TARGET DATA

      t_pack(1)=t_pack(1)+tcpu_pack
      t_unpack(1)=t_unpack(1)+tcpu_unpack

      return
      end
