subroutine epfftw(ux,uxc,uy,uz)
   use comp
   use ISO_C_BINDING
   use hipfort
   use hipfort_hipfft
   use hipfort_hipfft_enums
   implicit none

   real(b8) :: ux(nxpad,zisz,yjsz,nu)
   complex(b8) :: uxc(nxhpad,zisz,yjsz,nu)
   complex(b8) :: uy(xisz,nypad,zjsz,nu)
   complex(b8) :: uz(xisz,nzpad,yjsz,nu)

   complex(b8), allocatable :: buf1(:,:,:)
   complex(b8), allocatable :: buf2(:,:)

   integer(c_int) :: rank, batch, istride, ostride, idist, odist
   integer(c_int), target :: n, inembed, onembed
   integer(kind(HIPFFT_R2C)) :: transform_type
   integer(kind(HIPFFT_SUCCESS)) :: ierr_hip
   real*8 rtime1,rtime2

   rtime1=MPI_WTIME()

   !! use_work_buffer and hipfft_auto_alloc_mem are initialized in
   !! module.F90 and are used only in this subroutine

   !! checking use_work_buffer and hipfft_auto_alloc_mem
   if (use_work_buffer .lt. 0) use_work_buffer = 0
   if (use_work_buffer .gt. 1) use_work_buffer = 1
   if (hipfft_auto_alloc_mem .lt. 0) hipfft_auto_alloc_mem = 0
   if (hipfft_auto_alloc_mem .gt. 1) hipfft_auto_alloc_mem = 1

   !! only one of them can be set to 1, but not both
   if ((use_work_buffer + hipfft_auto_alloc_mem) .gt. 1) then
       hipfft_auto_alloc_mem = 1
       use_work_buffer = 0
   endif

   if (taskid .eq. 0) then
       write(6,*) "In epfftw: use_work_buffer = ", use_work_buffer
       write(6,*) "In epfftw: hipfft_auto_alloc_mem = ",hipfft_auto_alloc_mem
   endif

   !! x transform from real to complex (stride-1)
   ierr_hip = hipfftCreate(hip_plan_r2c)

   !! enable or disable auto allocation of work buffer
   ierr_hip = hipfftSetAutoAllocation(hip_plan_r2c, hipfft_auto_alloc_mem);

   rank=1
   n=nx
   inembed=nxpad
   istride=1
   idist=nxpad
   onembed=nxhp
   ostride=1
   odist=nxhp
   batch=zisz

#ifdef DOUBLE_PREC
   transform_type = HIPFFT_D2Z
#else
   transform_type = HIPFFT_R2C
#endif

   ierr_hip = hipfftMakePlanMany(hip_plan_r2c, rank, c_loc(n), &
                                 c_loc(inembed), istride, idist, &
                                 c_loc(onembed), ostride, odist, &
                                 transform_type, batch, c_loc(wbufsize_r2c))

   if (wbufsize_r2c .gt. wbufsize_single) wbufsize_single = wbufsize_r2c

   if (taskid.eq.0) write(6,*) "hipfft plan R2C: ierr = ",ierr_hip
   call custom_flush(6)

   !! z complex transform in xk and kx (strided)
   ierr_hip = hipfftCreate(hip_plan_z)

   !! enable or disable auto allocation of work buffer
   ierr_hip = hipfftSetAutoAllocation(hip_plan_z, hipfft_auto_alloc_mem);

   rank=1
   n=nz
   inembed=nz
   istride=xisz
   idist=1
   onembed=nz
   ostride=xisz
   odist=1
   batch=xisz

#ifdef DOUBLE_PREC
   transform_type = HIPFFT_Z2Z
#else
   transform_type = HIPFFT_C2C
#endif

   ierr_hip = hipfftMakePlanMany(hip_plan_z, rank, c_loc(n), &
                                 c_loc(inembed), istride, idist, &
                                 c_loc(onembed), ostride, odist, &
                                 transform_type, batch, c_loc(wbufsize_z))

   if (taskid.eq.0) write(6,*) "hipfft plan z: ierr = ", ierr_hip
   call custom_flush(6)

   if (wbufsize_z .gt. wbufsize_single) wbufsize_single = wbufsize_z

   !! y complex transform in xk and kx (strided)
   ierr_hip = hipfftCreate(hip_plan_y)

   !! enable or disable auto allocation of work buffer
   ierr_hip = hipfftSetAutoAllocation(hip_plan_y, hipfft_auto_alloc_mem);

   rank=1
   n=ny
   inembed=ny
   istride=xisz
   idist=1
   onembed=ny
   ostride=xisz
   odist=1
   batch=xisz

#ifdef DOUBLE_PREC
   transform_type = HIPFFT_Z2Z
#else
   transform_type = HIPFFT_C2C
#endif

   ierr_hip = hipfftMakePlanMany(hip_plan_y, rank, c_loc(n), &
                                 c_loc(inembed), istride, idist, &
                                 c_loc(onembed), ostride, odist, &
                                 transform_type, batch, c_loc(wbufsize_y))

   if (taskid.eq.0) write(6,*) "hipfft plan y: ierrn = ", ierr_hip
   call custom_flush(6)

   if (wbufsize_y .gt. wbufsize_single) wbufsize_single = wbufsize_y

   !! x transform from complex to real (stride-1)
   ierr_hip = hipfftCreate(hip_plan_c2r)

   !! enable or disable auto allocation of work buffer
   ierr_hip = hipfftSetAutoAllocation(hip_plan_c2r, hipfft_auto_alloc_mem);


   rank=1
   n=nx
   inembed=nxhppad
   istride=1
   idist=nxhppad
   onembed=nxpad
   ostride=1
   odist=nxpad
   batch=zisz

#ifdef DOUBLE_PREC
   transform_type = HIPFFT_Z2D
#else
   transform_type = HIPFFT_C2R
#endif

   ierr_hip = hipfftMakePlanMany(hip_plan_c2r, rank, c_loc(n), &
                                 c_loc(inembed), istride, idist, &
                                 c_loc(onembed), ostride, odist, &
                                 transform_type, batch, c_loc(wbufsize_c2r))

   if (taskid.eq.0) write(6,*) "hipfft plan C2R: ierr = ", ierr_hip
   call custom_flush(6)

   if (wbufsize_c2r .gt. wbufsize_single) wbufsize_single = wbufsize_c2r

   !! Now that I know the maximum size needed for the work buffer, allocate the 
   !! single work buffer and assign it to the appropriate plan if needed. 

   !! short-circuit the creation/use of the work buffer either because you don't 
   !! want to use it or because you're letting hip automatically manage it
   if (use_work_buffer .eq. 0) wbufsize_single = zeroST

   if (wbufsize_single .gt. zeroST) then

      ierr_hip = hipMalloc(wbuffer_single, wbufsize_single)
      if (taskid.eq.0) then
         write(6,*) "Creating work buffer with size (in GB) = ",(dble(wbufsize_single)/1024.0/1024.0/1024.0)
      endif

      if (wbufsize_r2c .gt. zeroST) then
         ierr_hip = hipfftSetWorkArea(hip_plan_r2c, wbuffer_single);
         if (taskid.eq.0) write(6,*) "r2c using fft work buffer"
      end if

      if (wbufsize_z .gt. zeroST) then
         ierr_hip = hipfftSetWorkArea(hip_plan_z, wbuffer_single);
         if (taskid.eq.0) write(6,*) "z using fft work buffer"
      end if

      if (wbufsize_y .gt. zeroST) then
         ierr_hip = hipfftSetWorkArea(hip_plan_y, wbuffer_single);
         if (taskid.eq.0) write(6,*) "y using fft work buffer"
      end if

      if (wbufsize_c2r .gt. zeroST) then
         ierr_hip = hipfftSetWorkArea(hip_plan_c2r, wbuffer_single);
         if (taskid.eq.0) write(6,*) "c2r using fft work buffer"
      end if

   endif

   return
end
