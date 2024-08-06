      subroutine xkcomm2_gpu (source,dest,nv)
!
      use com
      use timers_comm
      use timers_comp
      use timers_tran
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none

      ! assume nxpad=nx and nzpad=nz

      integer i,iv,nv
      complex(b8) :: source(xisz,nzpad,yjsz,nv)
      complex(b8), target :: dest(xisz,nypad,zjsz,nv)

      real(b8) tcpu_pack,tcpu_unpack,tcpu_0

      real*8 rtime1,rtime2

      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

      rtime1=MPI_WTIME()

!     if (taskid.eq.0) write(6,*) "enter xkcomm2_gpu, nv=",nv
 
      ! transpose from z pencils to y pencils


      t2_comm = t2_comm - MPI_Wtime() 

      !$OMP TARGET DATA MAP(tofrom:source,dest) &
      !$OMP MAP(alloc:sndbuf_col,rcvbuf_col)

      tcpu_pack=0.
      tcpu_unpack=0.
      tcpu_0=t_xkcomm2(3)

      do iv=1,nv

         ! Pack the sndbuf_col

         rtime1=MPI_WTIME()

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,yp,xp,zp,z) &
         !$OMP SHARED(jproc,yjsz,xisz,zjsz,iv, &
         !$OMP kjst,sndbuf_col,source)
         do i=1,jproc
            do yp=1,yjsz
               do zp=1,zjsz
                  do xp=1,xisz
                     z=kjst(i-1)+zp-1
                     sndbuf_col(xp,yp,zp,iv,i) = source(xp,z,yp,iv)
                  enddo
               enddo
            enddo
         enddo
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         rtime1=MPI_WTIME()-rtime1
         tcpu_pack=tcpu_pack+rtime1
         tcpu_other=tcpu_other+rtime1
         t_xkcomm2(3)=t_xkcomm2(3)+rtime1

      end do ! do nv
   
      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE FROM(sndbuf_col)
      end if

      ! Now send the buffers (exchange data)

#ifdef DETAIL_MPI
      if (mpi_detail.eq.1) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
#endif

      t_alltoall = t_alltoall - MPI_Wtime()
      rtime1=MPI_WTIME()
 
      if (gpumpi.eq.0) then
         call mpi_alltoall (sndbuf_col, xisz*zjsz*yjsz*nv, mpicomplex, &
                           rcvbuf_col, xisz*zjsz*yjsz*nv, mpicomplex, &
                            mpi_comm_col, ierr)
      else
         !$OMP TARGET DATA USE_DEVICE_PTR(sndbuf_col, rcvbuf_col)
         call mpi_alltoall (sndbuf_col, xisz*zjsz*yjsz*nv, mpicomplex, &
                           rcvbuf_col, xisz*zjsz*yjsz*nv, mpicomplex, &
                            mpi_comm_col, ierr)
         !$OMP END TARGET DATA
      end if

 
      rtime2=MPI_WTIME()
      t_a2a=rtime2-rtime1
      t_alltoall = t_alltoall + MPI_Wtime()
      t_xkcomm2(1)=t_xkcomm2(1)+rtime2-rtime1

      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE TO(rcvbuf_col)
      end if

      do iv=1,nv

         rtime1=MPI_WTIME()
         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,xp,yp,zp,y) &
         !$OMP SHARED(jproc,xisz,yjsz,zjsz,iv, &
         !$OMP jjst,rcvbuf_col,dest)
         do i=1,jproc 
            do yp=1,yjsz
               do zp=1,zjsz
                  do xp=1,xisz
                     y=jjst(i-1)+yp-1
                     dest(xp,y,zp,iv) = rcvbuf_col(xp,yp,zp,iv,i)
                  enddo
               enddo
            enddo
         enddo
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
         rtime1=MPI_WTIME()-rtime1
         tcpu_unpack=tcpu_unpack+rtime1
         tcpu_other=tcpu_other+rtime1
         t_xkcomm2(3)=t_xkcomm2(3)+rtime1

         ! perform FFT in y 

         rtime1=MPI_WTIME()
         !$OMP TARGET DATA USE_DEVICE_PTR(dest)
         do zp=1,zjsz

#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecZ2Z(hip_plan_y, c_loc(dest(1,1,zp,iv)), &
                                     c_loc(dest(1,1,zp,iv)), HIPFFT_FORWARD)
#else
            ierr_hip = hipfftExecC2C(hip_plan_y, c_loc(dest(1,1,zp,iv)), &
                                     c_loc(dest(1,1,zp,iv)), HIPFFT_FORWARD)
#endif

         end do
         !$OMP END TARGET DATA

         ierr = hipDeviceSynchronize()

         t_xkcomm2(2)=t_xkcomm2(2)+MPI_WTIME()-rtime1

      end do ! do nv

      !$OMP END TARGET DATA

      
      t2_comm = MPI_Wtime() + t2_comm
      i2_comm = i2_comm + 1
      
      t_pack(4)=t_pack(4)+tcpu_pack
      t_unpack(4)=t_unpack(4)+tcpu_unpack

      return
      end

