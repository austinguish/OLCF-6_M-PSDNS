subroutine mpisetup()
   use com
   use omp_lib
   implicit none
        
   ! re-developed by PKY, 8/31/2021
   ! strictly pencils only, with no cylindrical truncation.

   integer i,j,k,n1,ixp,lu,ii
   logical iex,exs
   real rnx
   integer ithr,itask

   ! mpi process info

   integer dims(2),  cartid(2)
   logical periodic(2),remain_dims(2)
   logical reorder

   integer, allocatable :: zist_all(:),zjst_all(:)
   integer, allocatable :: yjst_all(:),xist_all(:)

   dims(1) = 0
   dims(2) = 0

   ! numtasks is divided into a iproc x jproc stencil
 
   call MPI_Dims_create(numtasks,2,dims,ierr)

#ifdef ONED 
   dims(2) = numtasks
   dims(1) = 1
#endif	
   if (taskid.eq.0) write (6,*) 'dims from mpi=',dims(1),dims(2)

   if(dims(1) .gt. dims(2)) then
      dims(1) = dims(2)
      dims(2) = numtasks / dims(1)
   endif
   if (taskid.eq.0) write (6,*) 'dims from mpi=',dims(1),dims(2)

   inquire(file='dims',exist=iex)
   if (iex) then
      if (taskid.eq.0) print *, 'Reading grid from file dims'
      open (999,file='dims')
      read (999,*) dims(1), dims(2)
      close (999)
   else
      if (taskid.eq.0) print *, 'Creating grid with mpi_dims_create'
   endif  


   iproc = dims(1)
   jproc = dims(2)

   if (iproc*jproc.ne.numtasks) then
      if (taskid.eq.0) then
         write (6,*)  'ABORT: invalid user-specified choice of iproc x jproc!'
         write (6,*)  'Correct choices in the dims file'
         write (6,*) 'iproc,jproc,numtasks=',iproc,jproc,numtasks
         call MPI_ABORT (MPI_COMM_WORLD,ierr)
      end if
   end if

#ifdef DETAIL_MPI
   ! set default group size for reporting of MPI statistics
   ! but take user-specified value if available.

   mpistat_groupsize=iproc
   mpi_detail=1
   if (taskid.eq.0) then
      inquire (file='input.mpistat',exist=exs)
      if (exs) then
         open (file='input.mpistat',newunit=lu)
         read (lu,*) mpistat_groupsize
         read (lu,*,end=2) mpi_detail
 2       continue
         close (lu)
      end if
      write  (6,*) 'mpistat_groupsize=',mpistat_groupsize
      write  (6,*) 'mpi_detail=',mpi_detail
   end if

   call MPI_BCAST (mpistat_groupsize,1,MPI_INTEGER,0, &
                   MPI_COMM_WORLD,mpierr)
   call MPI_BCAST (mpi_detail,1,MPI_INTEGER,0, &
                   MPI_COMM_WORLD,mpierr)
#endif

   ! manipulation necessary for arguments into  MPI_CART_CREATE,
   ! which assumes row-major ordering

   i = dims(1)  
   dims(1) = dims(2)
   dims(2) = i

   periodic(1) = .false.
   periodic(2) = .false.
   reorder = .false.
   ! creating cartesian processor grid
   call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic, &
         reorder,mpi_comm_cart,ierr)
   ! Obtaining process ids in the cartesian grid
   call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)

   ! manipulation necessary for output from MPI_CART_COORDS,
   ! which needs to be flipped from row-major ordering to column-major

   ipid = cartid(2)
   jpid = cartid(1)

   allocate (ipid_all(0:numtasks-1))
   allocate (jpid_all(0:numtasks-1))

   call  MPI_ALLGATHER (ipid,1,MPI_INTEGER,ipid_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)
   call  MPI_ALLGATHER (jpid,1,MPI_INTEGER,jpid_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)

   ! Proceed to create the row and column communicators

   ! Again, because MPI_CART_CREATE and MPI_CART_COORDS use row-major
   ! ordering, basically flipping the first and second dimensions,
   ! here we contructs the column sub-communicator from the
   ! 'first dimension', and row communicator from the 'second dimension'

   remain_dims(1) = .true.
   remain_dims(2) = .false.
   call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)

   remain_dims(1) = .false.
   remain_dims(2) = .true.
   call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)


   ! mapping i onto iproc, i onto jproc, j onto iproc etc.
   allocate (iist(0:iproc-1))
   allocate (iisz(0:iproc-1))
   allocate (iien(0:iproc-1))
   allocate (jjst(0:jproc-1))
   allocate (jjsz(0:jproc-1))
   allocate (jjen(0:jproc-1))
   allocate (kist(0:iproc-1))
   allocate (kisz(0:iproc-1))
   allocate (kien(0:iproc-1))
   allocate (kjst(0:jproc-1))
   allocate (kjsz(0:jproc-1))
   allocate (kjen(0:jproc-1))
 
   !Mapping 3-D data arrays onto 2-D process grid
   ! (nx+2,ny+2,nz) => (iproc,jproc)      
  
   call MapDataToProc(nxh,iproc,iist,iien,iisz)
   call MapDataToProc(nypad,jproc,jjst,jjen,jjsz)
   call MapDataToProc(nzpad,iproc,kist,kien,kisz)
   call MapDataToProc(nz,jproc,kjst,kjen,kjsz)
 
   allocate(mymap(0:iproc-1))
   allocate(inverse_map(0:iproc-1))

   do ii=0,iproc-1
      mymap(ii)=ii
   end do

   !xist = iist(ipid)
   !yjst = jjst(jpid)
   !zist = kist(ipid)
   !zjst = kjst(jpid)
   !xisz = iisz(ipid)
   !yjsz = jjsz(jpid)
   !zisz = kisz(ipid)
   !zjsz = kjsz(jpid)
   !xien = iien(ipid)
   !yjen = jjen(jpid)
   !zien = kien(ipid)
   !zjen = kjen(jpid)
 
   allocate (zist_all(0:numtasks-1))
   allocate (zjst_all(0:numtasks-1))
   allocate (yjst_all(0:numtasks-1))
   allocate (xist_all(0:numtasks-1))

   xisz = nx/2/iproc
   xist = ipid*xisz+1
   xien = xist+xisz-1
   zisz = nz/iproc
   zist = ipid*zisz+1
   zien = zist+zisz-1
   zjsz = nz/jproc
   zjst = jpid*zjsz+1
   zjen = zjst+zjsz-1
   yjsz = ny/jproc
   yjst = jpid*yjsz+1
   yjen = yjst+yjsz-1

   call  MPI_ALLGATHER (zist,1,MPI_INTEGER,zist_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)
   call  MPI_ALLGATHER (zjst,1,MPI_INTEGER,zjst_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)
   call  MPI_ALLGATHER (yjst,1,MPI_INTEGER,yjst_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)
   call  MPI_ALLGATHER (xist,1,MPI_INTEGER,xist_all,1,MPI_INTEGER, &
                        MPI_COMM_WORLD,mpierr)
   if (taskid.eq.1) then
      open (file='ipid_xist_info',newunit=lu)
      write (lu,"('numtasks, iproc,jproc=',4i6)") numtasks,iproc,jproc
      do itask=0,numtasks-1
         write (lu,"('taskid,ipid,jpid,xist,yjst,zist,zjst=',7i5)") itask, &
                ipid_all(itask),jpid_all(itask), &
                xist_all(itask),yjst_all(itask), &
                zist_all(itask),zjst_all(itask)
      end do
      close (lu)
   end if

   do i=0,iproc-1
      iist(i)=i*xisz+1
      iisz(i)=xisz
      iien(i)=iist(i)+xisz-1
      kist(i)=i*zisz+1
      kisz(i)=zisz
      kien(i)=kist(i)+zisz-1
   end do

   do j=0,jproc-1
      jjst(j)=j*yjsz+1
      jjsz(j)=yjsz
      jjen(j)=jjst(j)+yjsz-1
      kjst(j)=j*zjsz+1
      kjsz(j)=zjsz
      kjen(j)=kjst(j)+zjsz-1
   end do


#ifdef OPENMP
   !$OMP PARALLEL private(ithr)
   num_thr = OMP_GET_NUM_THREADS()
   ithr = OMP_GET_THREAD_NUM()
   if (taskid.eq.0) then
      write (6,*) 'mpisetup: taskid,ithr,num_thr=', taskid,ithr,num_thr
   end if
   !$OMP END PARALLEL
#else
   ithr=0
   num_thr = 1
#endif


       
end subroutine

!==================================================================       

subroutine MapDataToProc (data,proc,st,en,szs)
     
   implicit none 
   integer data,proc,st(0:proc-1),en(0:proc-1),szs(0:proc-1)
   integer i,size,nadd,size2
   size=data/proc
   nadd=mod(data,proc)
   size2=size
   if(nadd.ne.0) size2= size2+1
   st(0) = 1
   szs(0) = size2
   en(0) = size2
   if (proc .gt. 1) then
      do i=1,proc-1
         size2=size
         if (i.lt.nadd) size2=size2+1
         st(i) = st(i-1) + szs(i-1)
         szs(i) = size2
         en(i) = en(i-1) + size2
      enddo
      en(proc-1)= data 
      szs(proc-1)= data-st(i-1)+1
   endif
 
end subroutine
!==================================================================       


       
