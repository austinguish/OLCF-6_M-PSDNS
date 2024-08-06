subroutine com_set

   use com
   use hipfort

   implicit none

   integer :: m, nv
   integer(c_size_t) :: n_points, nr_points, nc_points

   allocate (kx(nxh),ky(ny),kz(nz))
   allocate (bk1(nxh,2))
   allocate (bk2(ny,2))
   allocate (bk3(nz,2))
   allocate (bkk3(nz,2))

   xiszST = int(xisz,kind=c_size_t)
   ziszST = int(zisz,kind=c_size_t)
   yjszST = int(yjsz,kind=c_size_t)
   zjszST = int(zjsz,kind=c_size_t)
   iprocST = int(iproc,kind=c_size_t)
   jprocST = int(jproc,kind=c_size_t)

   nv = 1
   nvST = int(nv,kind=c_size_t)
   
   allocate(buf(nxhppad,zisz))

   !! nr_points should be the same as nc_points, but
   !! to make sure, I perform a sanity test and use 
   !! the largest value when allocating the array
   !! NOTE: these values are c_size_t to make sure 
   !!       that they can hold very large integer values
   nr_points = xiszST*ziszST*yjszST*nvST*iprocST
   nc_points = zjszST*xiszST*yjszST*nvST*jprocST
   n_points = nr_points
   if (n_points < nc_points) n_points = nc_points
   allocate(sndbuf_single(n_points))
   allocate(rcvbuf_single(n_points))

   !! get creative with the single set of send and receive buffers
   !! that have been allocated since the row and column comms
   !! have a different shape
   !!    First, get the c address of the buffers
   sndbuf_buffer = c_loc(sndbuf_single(1))
   rcvbuf_buffer = c_loc(rcvbuf_single(1))
   !!    Second, set the shape of the fortran pointers
   call c_f_pointer(sndbuf_buffer,sndbuf_row,[xisz,zisz,yjsz,nv,iproc])
   call c_f_pointer(rcvbuf_buffer,rcvbuf_row,[xisz,zisz,yjsz,nv,iproc])
   call c_f_pointer(sndbuf_buffer,sndbuf_col,[xisz,yjsz,zjsz,nv,jproc])
   call c_f_pointer(rcvbuf_buffer,rcvbuf_col,[xisz,yjsz,zjsz,nv,jproc])
  
   ! define some used constants
   imagi=cmplx(0.,1.)
   pi=4.*atan(1.)

   ! initialize wavenumbers
   do x=1,nxhpad
      kx(x)=x-1
   end do
   do y=1,nypad
      ky(y)=y-1
      if(y.gt.nyhppad) ky(y)=ky(y)-ny
   end do
   do z=1,nzpad
      kz(z)=z-1
      if(z.gt.nzhppad) kz(z)=kz(z)-nz
   end do

   kmax=sqrt(2.)/3.*max(nxpad*beta1,nypad*beta2,nzpad*beta3)

   do m=1,2
      !form wave numbers for differentiation
      do x=1,nxh
         bk1(x,m)=b11(m)*kx(x)
      end do
      do y=1,ny
         bk2(y,m)=b22(m)*ky(y)
      end do
      do z=1,nz
         bk3(z,m)=b33(m)*kz(z)
         bkk3(z,m)=bk3(z,m)*kz(z)
      end do
   end do

   allocate(bk1i(nxh,2))
   allocate(bk3i(nz,2))
   allocate(bk2i(ny,2))
   m=2
   do x=1,nxh
      bk1i(x,1)=imagi*kx(x) ! i*kx
      bk1i(x,2)=imagi*bk1(x,m) ! i*b11*kx
   end do
   do z=1,nz
      bk3i(z,1)=imagi*bk3(z,m) ! i*b33*kz
      bk3i(z,2)=imagi*kz(z) ! i*kz
   end do
   do y=1,ny
      bk2i(y,1)=imagi*ky(y) ! i*ky
      bk2i(y,2)=imagi*bk2(y,m) ! i*b22*ky
   end do

   allocate (sx(nxh,rkstages),sy(ny,rkstages),sz(nz,rkstages))

   if (taskid.eq.0) write (6,*) ' exit com_set'

return
end
