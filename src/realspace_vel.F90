! Operations in physical space

subroutine realspace_vel (ux,m,rkstep)

   use comp
   implicit none
   real(b8)    :: ux(nx,zisz,yjsz,nu)

   integer :: m,i, rkstep,icall
   real(b8) s1,s2,s3,xnorm,ynorm,znorm,vel,bv2
   real(b8) tmpu,tmpv,tmpw,tmpvmax
   
   integer ithr,OMP_GET_THREAD_NUM,yp1,yp2
   real(b8), allocatable :: vmax(:)


   is1=4+nc
   is2=5+nc
   
   s1=sqrt(b11(m))
   s2=sqrt(b22(m))
   s3=sqrt(b33(m))

   xnorm=b11(m)*nx
   ynorm=b22(m)*ny
   znorm=b33(m)*nz

   !$OMP TARGET DATA MAP(tofrom:ux) MAP(to:b11,b22,b33,xnorm,ynorm,znorm)

   ! Courant number and convective terms in physical space
 
   if (rkstep.eq.1) then
      
      velmax=0.

      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
      !$OMP PRIVATE(yp,zp,x,bv2,vel) SHARED(ux,b11,b22,b33,is1,is2,xnorm, &
      !$OMP ynorm,znorm,yjsz,zisz,nx,m) &
#ifdef USE_MAP
      !$OMP MAP(tofrom:velmax) &
#endif
      !$OMP REDUCTION(max:velmax)
      do yp=1,yjsz
         do zp=1,zisz
            do x=1,nx
               bv2=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
               ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
               ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
               vel=xnorm*abs(ux(x,zp,yp,1)) + &
                   ynorm*abs(ux(x,zp,yp,2)) + &
                   znorm*abs(ux(x,zp,yp,3))
               velmax=max(vel,velmax)
               ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2
               ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
               ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

      tmpvmax=velmax
      call MPI_ALLREDUCE (tmpvmax,velmax,1,mpireal,MPI_MAX, &
         MPI_COMM_WORLD,mpierr)

   else
      
      !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO DEFAULT(NONE) COLLAPSE(3) &
      !$OMP PRIVATE(yp,zp,x,bv2) SHARED(ux,b11,b22,b33,is1,is2,xnorm, &
      !$OMP ynorm,znorm,yjsz,zisz,nx,m)
      do yp=1,yjsz
         do zp=1,zisz
            do x=1,nx
               bv2=b22(m)*ux(x,zp,yp,2)*ux(x,zp,yp,2)
               ux(x,zp,yp,is2)=ux(x,zp,yp,1)*ux(x,zp,yp,3)
               ux(x,zp,yp,is1)=ux(x,zp,yp,2)*ux(x,zp,yp,3)
               ux(x,zp,yp,3)=b33(m)*ux(x,zp,yp,3)*ux(x,zp,yp,3)-bv2
               ux(x,zp,yp,2)=ux(x,zp,yp,1)*ux(x,zp,yp,2)
               ux(x,zp,yp,1)=b11(m)*ux(x,zp,yp,1)*ux(x,zp,yp,1)-bv2
            end do
         end do
      end do
      !$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO

   end if

   !$OMP END TARGET DATA

   return
end
