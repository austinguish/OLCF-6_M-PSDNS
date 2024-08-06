! Operations from physical space to wavenumber space
subroutine transform (ux,uy,uz,m) 
   use comp
   use timers_comm
   use timers_tran
   use timers_rkstep
   use hipfort
   use hipfort_hipfft
   use hipfort_hipfft_enums
   implicit none

   real(b8) factor
   complex(b8), allocatable :: buf2(:,:),buf1(:,:,:)
   real(b8)    :: ux(nx,zisz,yjsz,nu)
   complex(b8), target :: uz(xisz,nz,yjsz,nu)
   complex(b8) :: uy(xisz,ny,zjsz,nu)
   integer :: m,i, novt,rkstep,i2f,x2,iz,ix,iz2,dnz,iy,l
   real(b8) s1,s2,s3

   integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

   integer ithr,omp_get_thread_num,yp1,yp2
   
   integer i1,ii
   real(8) rtime1

   if (taskid.eq.0.and.istep.eq.1) write (6,*) 'call transform_gpu, kstep=',kstep

   s1=sqrt(b11(m))
   s2=sqrt(b22(m))
   s3=sqrt(b33(m))
   is1=4+nc
   is2=5+nc
   factor=1./real(nx,b8)

   ! transform to x-wavenumber space, and take transpose 
   !
   novt=5+nc
   
   rtime1=MPI_WTIME()
   !
   !$OMP TARGET DATA MAP(tofrom:ux, uy, uz) MAP(to:bk1i, bk3i) &
   !$OMP MAP(alloc:buf, sndbuf_single, rcvbuf_single)
   !
   ! Transform R2C in X, transpose x <-> Z, truncate in X if needed

   do i=1,novt
      call xkcomm1_gpu (ux(1,1,1,i),uz(1,1,1,i),1)
      rt_a2a(3,i)=t_a2a
   end do

   t_trans(1)=t_trans(1)+MPI_WTIME()-rtime1

   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
    
   ! transform in z
    
   rtime1=MPI_WTIME()
 
   !$OMP TARGET DATA USE_DEVICE_PTR(uz)
   do i=1,novt
      do yp=1,yjsz

#ifdef DOUBLE_PREC
         ierr_hip = hipfftExecZ2Z(hip_plan_z, c_loc(uz(1,1,yp,i)), &
                                  c_loc(uz(1,1,yp,i)), HIPFFT_FORWARD)
#else
         ierr_hip = hipfftExecC2C(hip_plan_z, c_loc(uz(1,1,yp,i)), &
                                  c_loc(uz(1,1,yp,i)), HIPFFT_FORWARD)
#endif

      end do
   end do
   !$OMP END TARGET DATA

   ierr = hipDeviceSynchronize()
 
   call MPI_BARRIER (MPI_COMM_WORLD,mpierr)
   t_trans(2)=t_trans(2)+MPI_WTIME()-rtime1

   !#################################################
   ! operations for convective terms
   !#################################################

   rtime1=MPI_WTIME()
   !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
   !$OMP PRIVATE(yp,xp,z,y,x) SHARED(uz,bk1i,bk3i,yjst,xist,yjsz, &
   !$OMP xisz,nz,nzhp,is2)
   do yp=1,yjsz
      do z=1,nz
         do xp=1,xisz
            y=yjst+yp-1
            x=xist+xp-1
            uz(xp,z,yp,1)=bk1i(x,1)*uz(xp,z,yp,1)+bk3i(z,1)*uz(xp,z,yp,is2)
            uz(xp,z,yp,3)=bk1i(x,2)*uz(xp,z,yp,is2)+bk3i(z,2)*uz(xp,z,yp,3)
         end do
      end do
   end do
   !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   t_trans(3)=t_trans(3)+MPI_WTIME()-rtime1

   !###################################################
        
   i2f=4+nc   

   ! Truncate in Z if needed, Transpose ZXY -> YZX (Z- to Y-pencils), transform in Y 
   rtime1=MPI_WTIME()

   do i=1,i2f
      call xkcomm2_gpu (uz(1,1,1,i),uy(1,1,1,i),1)
      rt_a2a(4,i)=t_a2a
   end do

   t_trans(4)=t_trans(4)+MPI_WTIME()-rtime1

   !$OMP END TARGET DATA

   return
end subroutine transform
