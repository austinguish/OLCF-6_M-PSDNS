! Input: X-pencils (real)
! Output: Z-pencils (complex)
!
! This routine performs real-to-complex FFT,
! truncates the X dimension from nxpad to nx, and
! transposes the data to arrange Z-pencils while 
! interchanging the order of X and Z indices
!
! This version performs alltoall one variable at a time,
! thus allowing smaller send and receive buffers
! 2D domain decomposition: row_communicator only in this routine

      subroutine xkcomm1_gpu (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none

      ! assume nxpad=nx and nzpad=nz

      real(b8), target :: source(nxpad,zisz,yjsz,nv)
      complex(b8) dest(xisz,nzpad,yjsz,nv)

      integer i,j,nv
      real(b8) factor
!
      real*8 rtime1,rtime2,rtime0
      real(b8) tcpu_pack,tcpu_unpack

      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

      factor=1./real(nxpad,b8)

      tp1_comm = tp1_comm - MPI_Wtime() 

      !$OMP TARGET DATA MAP(tofrom:source, dest) &
      !$OMP MAP(alloc:buf, sndbuf_row, rcvbuf_row)
!
      tcpu_pack=0.
      tcpu_unpack=0.

      do j=1,nv

         do yp=1,yjsz

         rtime1=MPI_WTIME()
         ! Transform R2C

         !$OMP TARGET DATA USE_DEVICE_PTR(source,buf)

#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecD2Z(hip_plan_r2c, c_loc(source(1,1,yp,j)), &
                                     c_loc(buf(1,1)))
#else
            ierr_hip = hipfftExecR2C(hip_plan_r2c, c_loc(source(1,1,yp,j)), &
                                     c_loc(buf(1,1)))
#endif

         !$OMP END TARGET DATA

         !! sync device
         ierr = hipDeviceSynchronize()


         rtime2=MPI_WTIME()
         tcpu_fft=tcpu_fft+(rtime2-rtime1)
         t_xkcomm1(2)=t_xkcomm1(2)+(rtime2-rtime1)

         
! Pack the send buffer for exchanging z and x (within a given y plane ) into sndbuf_row
         rtime1=MPI_WTIME()

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
         !$OMP DEFAULT(NONE) PRIVATE(i,zp,xp,x) &
         !$OMP SHARED(iist,sndbuf_row,buf, &
         !$OMP j,iproc,zisz,xisz,yp,factor)
         do i=1,iproc
            do zp=1,zisz
               do xp=1,xisz
                  x=iist(i-1)+xp-1
                  sndbuf_row(xp,zp,yp,j,i)= buf(x,zp)*factor
               end do
            end do
         end do
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         rtime2=MPI_WTIME()
         tcpu_pack=tcpu_pack+(rtime2-rtime1)
         tcpu_other=tcpu_other+(rtime2-rtime1)
         t_xkcomm1(3)=t_xkcomm1(3)+(rtime2-rtime1)

         end do ! do yp

      end do ! do nv

      ! perform alltoall
      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE FROM(sndbuf_row)
      end if

#ifdef DETAIL_MPI
      if (mpi_detail.eq.1) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
#endif

      t_alltoall = t_alltoall - MPI_Wtime()
      t1_comm = t1_comm - MPI_Wtime() 
      rtime1=MPI_WTIME()

      if (gpumpi.eq.0) then
         call MPI_ALLTOALL(sndbuf_row, xisz*zisz*yjsz*nv, mpicomplex, &
            rcvbuf_row, xisz*zisz*yjsz*nv, mpicomplex, &
            mpi_comm_row, ierr)
      else
         !$OMP TARGET DATA USE_DEVICE_PTR(sndbuf_row,rcvbuf_row)
         call MPI_ALLTOALL(sndbuf_row, xisz*zisz*yjsz*nv, mpicomplex, &
            rcvbuf_row, xisz*zisz*yjsz*nv, mpicomplex, &
            mpi_comm_row, ierr)
         !$OMP END TARGET DATA
      end if

      rtime2=MPI_WTIME()
      t_a2a=rtime2-rtime1
      t1_comm = t1_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
      t_xkcomm1(1) = t_xkcomm1(1) + rtime2-rtime1

      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE TO(rcvbuf_row)
      end if

      do j=1,nv

         rtime1=MPI_WTIME()
         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,yp,zp,z,xp) &
         !$OMP SHARED(kist,rcvbuf_row,dest, &
         !$OMP iproc,yjsz,zisz,xisz,j)
         do i=1,iproc
            do yp=1,yjsz
               do zp=1,zisz
                  do xp=1,xisz
                     z=kist(i-1)+zp-1
                     dest(xp,z,yp,j) = rcvbuf_row(xp,zp,yp,j,i)
                  enddo
               enddo
            enddo
         enddo ! do i
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
         rtime2=MPI_WTIME()
         tcpu_unpack=tcpu_unpack+(rtime2-rtime1)
         tcpu_other=tcpu_other+(rtime2-rtime1)
         t_xkcomm1(3) = t_xkcomm1(3) + rtime2-rtime1

      end do ! do nv

      !$OMP END TARGET DATA

      t_pack(3)=t_pack(3)+tcpu_pack
      t_unpack(3)=t_unpack(3)+tcpu_unpack

      return
      end
