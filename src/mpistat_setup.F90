      subroutine mpistat_setup()
        use mpicom
        implicit none

        ! Added by PKY, 1/24/2023
#ifdef DETAIL_MPI
        ! set up sub-communicators for MPI tasks that will collect and
        ! write MPI statistics info.

        logical flag
        logical periodic(2),remain_dims(2)
        integer pcartid(2)
        integer, allocatable :: itask(:,:)


        if (mpi_detail.eq.0) return

        mpistat_ngroups=numtasks/mpistat_groupsize

        mpistat_dims(2)=mpistat_groupsize
        mpistat_dims(1)=mpistat_ngroups

        call MPI_DIMS_CREATE (numtasks,2,mpistat_dims,ierr)


        ! creating cartesian processor grid
        periodic(1) = .false.
        periodic(2) = .false.
        call MPI_CART_CREATE (MPI_COMM_WORLD, 2, mpistat_dims, periodic, &
            .false., mpi_comm_stat,ierr)

        ! Obtaining process ids with in the cartesian grid
        call MPI_CART_COORDS (mpi_comm_stat,taskid,2,pcartid,ierr)


        allocate (itask(0:mpistat_groupsize-1,0:mpistat_ngroups-1))
        allocate (itask_grp(0:mpistat_groupsize-1,0:mpistat_ngroups-1))

        ipid_stat = pcartid(2)
        jpid_stat = pcartid(1)
        remain_dims(1) = .true.
        remain_dims(2) = .false.
        call MPI_CART_SUB (mpi_comm_stat, remain_dims, &
                          mpi_comm_stat_across, ierr)
        remain_dims(1) = .false.
        remain_dims(2) = .true.
        call MPI_CART_SUB (mpi_comm_stat, remain_dims, & 
                           mpi_comm_stat_within, ierr)

        itask=0
        itask(ipid_stat,jpid_stat)=taskid

        call MPI_ALLREDUCE (itask,itask_grp,numtasks,MPI_INTEGER,MPI_SUM, &
                            MPI_COMM_WORLD,mpierr)
        
        deallocate (itask)

        allocate (rt_a2a_grp(4,5,0:mpistat_groupsize-1))
#endif

      end subroutine
!==================================================================       


       
