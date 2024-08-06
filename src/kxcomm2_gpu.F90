!
! This routine reorders X,Z indices while packing the send buffer
! (using loop blocking),
! exchanges data to arrange X pencils, expandsX dimension to nxhppad
! and performs inverse complex-to-real FFT
! Input: Z pencils, complex
! Output: X-pencils, real 
!
! Multivariable version
!

      subroutine kxcomm2_gpu (source,dest,nv)

      use com
      use timers_comm
      use timers_comp
      use timers_tran
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none

      real(b8), target :: dest(nxpad,zisz,yjsz,nv)
      complex(b8) :: source(xisz,nzpad,yjsz,nv)

      real(b8) tcpu_pack,tcpu_unpack

      integer i,j,k,n,nv
      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

      real*8 rtime1,rtime2,rtime3,rtime0

      rtime0=MPI_WTIME()
      rtime1=MPI_WTIME()

      rtime2=MPI_WTIME()

      !!$OMP TARGET DATA MAP(to:source) MAP(from:dest) &
      !!$OMP MAP(alloc:buf,sndbuf_row,rcvbuf_row)

      !$OMP TARGET DATA MAP(to:source) &
      !$OMP MAP(alloc:buf,sndbuf_row,rcvbuf_row)

      tcpu_pack=0.
      tcpu_unpack=0.

      ! Pack send buffer

      ! Use MPI_Alltoall
       rtime1=MPI_WTIME()

      do j=1,nv

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(4) &
         !$OMP DEFAULT(NONE) PRIVATE(i,z,x,xp,zp,yp) &
         !$OMP SHARED(j,nv,yjsz,iproc,xisz,zisz, &
         !$OMP kist,sndbuf_row,source)
         do i=1,iproc
            do yp=1,yjsz
               do zp=1,zisz
                  do xp=1,xisz
                     z=kist(i-1)+zp-1
                     sndbuf_row(xp,zp,yp,j,i) = source(xp,z,yp,j)
                  enddo
               enddo
            enddo
         enddo ! do yp
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

         rtime2=MPI_WTIME()
         tcpu_pack=tcpu_pack+(rtime2-rtime1)
         tcpu_other=tcpu_other+(rtime2-rtime1)
         t_kxcomm2(3)=t_kxcomm2(3)+(rtime2-rtime1)

      end do ! do nv

      ! Exchange x-z buffers

#ifdef DETAIL_MPI
      if (mpi_detail.eq.1) then
         call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      end if
#endif

      t_alltoall = t_alltoall - MPI_Wtime()
      t4_comm = t4_comm - MPI_Wtime() 
      rtime3=MPI_WTIME()

      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE FROM(sndbuf_row)
      end if

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


      rtime3=MPI_WTIME()-rtime3
      t_a2a=rtime3
      t4_comm = t4_comm + MPI_Wtime() 
      t_alltoall = t_alltoall + MPI_Wtime()
      t_kxcomm2(1)=t_kxcomm2(1)+rtime3

      if (gpumpi.eq.0) then
         !$OMP TARGET UPDATE TO(rcvbuf_row)
      end if

      do j=1,nv

         ! Unpack receive buffers 

         !! use unstructed enter/exit data movement for dest
         !$OMP TARGET ENTER DATA MAP(alloc:dest)

         do yp=1,yjsz

            rtime1=MPI_WTIME()
            !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
            !$OMP DEFAULT(NONE) PRIVATE(i,z,x,xp,zp) &
            !$OMP SHARED(j,iproc,xisz,zisz,yjsz,yp, &
            !$OMP iist,buf,rcvbuf_row)
            do i=1,iproc
               do zp=1,zisz
                  do xp=1,xisz
                     x=iist(i-1)+xp-1
                     buf(x,zp) = rcvbuf_row(xp,zp,yp,j,i)
                  enddo
               enddo
            enddo
            !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
         
         
            ! Add and zero extra elements in X
            !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) &
            !$OMP DEFAULT(NONE) PRIVATE(zp,x) &
            !$OMP SHARED(zisz,nxhp,nxhppad,buf)
            do zp=1,zisz
               do x=nxhp,nxhppad
                  buf(x,zp) = cmplx(0.0,0.0)
               enddo
            enddo
            !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
            rtime2=MPI_WTIME()
            tcpu_unpack=tcpu_unpack+(rtime2-rtime1)
            tcpu_other=tcpu_other+(rtime2-rtime1)
            t_kxcomm2(3)=t_kxcomm2(3)+(rtime2-rtime1)
            ! C2R Transform    
            rtime1=MPI_WTIME()

            !$OMP TARGET DATA USE_DEVICE_PTR(buf,dest)
#ifdef DOUBLE_PREC
            ierr_hip = hipfftExecZ2D(hip_plan_c2r, c_loc(buf(1,1)), &
                                     c_loc(dest(1,1,yp,j)))
#else
            ierr_hip = hipfftExecC2R(hip_plan_c2r, c_loc(buf(1,1)), &
                                     c_loc(dest(1,1,yp,j)))
#endif
            !$OMP END TARGET DATA

            ierr = hipDeviceSynchronize()


            rtime2=MPI_WTIME()
            tcpu_fft=tcpu_fft+(rtime2-rtime1)
            t_kxcomm2(2)=t_kxcomm2(2)+(rtime2-rtime1)

         enddo ! do yp

         !$OMP TARGET EXIT DATA MAP(from:dest)

      end do ! do nv

      !$OMP END TARGET DATA
         
      t_pack(2)=t_pack(2)+tcpu_pack
      t_unpack(2)=t_unpack(1)+tcpu_unpack

      return
      end
