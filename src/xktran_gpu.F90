      ! Strided xktran

      subroutine xktran_gpu (Source,buf1,Dest,nv)
      use com
      use hipfort
      use hipfort_hipfft
      use hipfort_hipfft_enums
      implicit none


      real(b8) Source(nxpad,zisz,yjsz,nv)
      complex(b8), target :: buf1(xisz,nzpad,yjsz,nv)
      complex(b8) Dest(xisz,nypad,zjsz,nv)
      
      real(b8) factor,factor_yz
      integer i,m,iz,ix,x2,iy,a
      integer :: nv
      real*8 rtime1,rtime2,rtime3,rtime0
      integer ii

      integer(kind(HIPFFT_SUCCESS)) :: ierr_hip

      factor_yz=1./real(nypad*nzpad,b8)
      factor = factor_yz / nxpad


      ! take real-to-complex transform in x, then switch to x-lines
      !$OMP TARGET DATA MAP(tofrom:Source, buf1, Dest)

      call xkcomm1_gpu (Source,buf1,nv)

      ! Perform FFT in z for all x for a given y plane and normalize

      do i=1,nv

         !$OMP TARGET DATA USE_DEVICE_PTR(buf1)
         do y=1,yjsz

#ifdef DOUBLE_PREC
             ierr_hip = hipfftExecZ2Z(hip_plan_z, c_loc(buf1(1,1,y,i)), &
                                      c_loc(buf1(1,1,y,i)), HIPFFT_FORWARD)
#else
             ierr_hip = hipfftExecC2C(hip_plan_z, c_loc(buf1(1,1,y,i)), &
                                      c_loc(buf1(1,1,y,i)), HIPFFT_FORWARD)
#endif

         end do
         !$OMP END TARGET DATA

         !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(3) &
         !$OMP DEFAULT(NONE) PRIVATE(y,x,z) SHARED(yjsz,xisz,nzpad, &
         !$OMP buf1,i,factor_yz)
         do y=1,yjsz
            do z=1,nzpad
               do x=1,xisz
                  buf1(x,z,y,i) = buf1(x,z,y,i) * factor_yz
               enddo
            enddo
         enddo
         !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      enddo
      

      ! switch to y-lines and take FFT in y

      call xkcomm2_gpu (buf1,Dest,nv)
      

      !$OMP END TARGET DATA


      return
      end

