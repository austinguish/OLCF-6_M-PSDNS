subroutine input

    use com


    implicit none
    integer lu,i,j,nchar, ixxx
    character*120 caux
    real, allocatable :: rsendbuf(:)
    integer, allocatable :: isendbuf(:)

    integer no_err, no_end, io_status
    logical exs
    integer nchar2,nchar3
    character*2 tail
    real tl_hat,kmaxeta


    ! Perform alltoalls on the GPU


    if (taskid.eq.0) write (6,*) 'enter input'

    no_err=0
    no_end=0
 
    allocate (isendbuf(30),rsendbuf(30))
    isendbuf=0
    rsendbuf=0.

    !need to add variables to com module 

    !--------------------------------------------------------
    if (taskid.eq.0) then ! only task=0 opens/reads input

       open (file='input',newunit=lu)

       read (lu,fmt='(1x)')
       read (lu,*) nx,ny,nz,nc, gpumpi
       write (6,"('nx,ny,nz,nc,gpumpi=',5i6)") nx,ny,nz,nc,gpumpi
       isendbuf(1)=nx
       isendbuf(2)=ny
       isendbuf(3)=nz
       isendbuf(4)=nc

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) nsteps,iostep,istart,kstop,isave
       write (6,"('input : nsteps,iostep =',5i4)") nsteps,iostep,istart,&
          kstop,isave
 
       isendbuf(5)=nsteps; isendbuf(6)=iostep;
       isendbuf(7)=istart; isendbuf(8)=kstop;
       isendbuf(9)=isave;
 
       isendbuf(10)=gpumpi

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) entime,time0
       rsendbuf(2)=entime; rsendbuf(3)=time0;

       allocate (kinit(3+nc))
       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) (kinit(i),i=1,3)
       write (6,*) 'after read kinit...'
       call custom_flush(6)
       isendbuf(11)=kinit(1)
 
       allocate (fninit(3+nc))
       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) fninit(1)
       write (6,*) 'after read fninit(1)...'
       call custom_flush(6)
 
       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,fmt='(a120)') caux
       write (6,*) 'after read indir_fn...',caux
       call blanks(caux,nchar)
       indir_fn=caux(1:nchar)

       ! lines added 2/24/13 to detect identify last checkpoint
       ! (this is to enhance robustness of sequence of chained jobs,
       ! in case a prior run did not reach its own last intended checkpoint)
 
       caux=fninit(1)
       call blanks (caux,nchar3)
       if (nchar3.ne.7) go to 53

       if (fninit(1)(7:7).ne.'c')  then
 
          if (kinit(1).gt.0.and.istart.gt.0) then
             do ichkpt=20,1,-1
                write (tail,"(i2)") ichkpt
                caux=indir_fn(1:nchar-7)//'/chkptou.'//tail
                call blanks (caux,nchar2)
                inquire (file=caux(1:nchar2),exist=exs)
                write (6,*) ichkpt,caux(1:nchar2),exs
                if (exs) go to 51
             end do
 51          continue
             if (exs) then
                ! fninit(1)(7:7)=tail
                fninit(1)(7:8)=tail
                write (6,*) 'revised fninit=',fninit(1)
                write (6,*) caux(1:nchar2)
             end if
 52          continue
          end if
 
       end if
 
 53    continue
 
       ichkpt=0

       do i=2,3+nc
          fninit(i)=fninit(1)
       enddo


       read (lu,fmt='(1x)')
       read (lu,*,err=10,end=20) irz
       write (6,*) 'after read irz...'
       call custom_flush(6)
       isendbuf(13)=irz

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) viscos
       write (6,*) 'after read viscos...',viscos
       rsendbuf(4)=viscos
 

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) a11,a22,a33 !!! ,shear
       write (6,*) 'after read a11...'
       rsendbuf(5)=a11
       rsendbuf(6)=a22
       rsendbuf(7)=a33

       ! shear flows not considered in this version
       !rsendbuf(8)=shear
 
       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) beta1,beta2,beta3,bet12
       write (6,*) 'after read beta1...'
       rsendbuf(9)=beta1
       rsendbuf(10)=beta2
       rsendbuf(11)=beta3
       rsendbuf(12)=bet12

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) cfl,dt
       write (6,*) 'after read cfl...'
       rsendbuf(13)=cfl; rsendbuf(14)=dt

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) dtout,tfout
       write (6,*) 'after read dtout...'
       rsendbuf(15)=dtout; rsendbuf(16)=tfout

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) dtsave,tfsave
       write (6,*) 'after read dtsave...'
       rsendbuf(17)=dtsave; rsendbuf(18)=tfsave

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20) kshift
       write (6,*) 'after read kshift...'
       isendbuf(14)=kshift

       ! No forcing

       forc_type=0

       read (lu,fmt='(1x)',err=10,end=20)
       read (lu,*,err=10,end=20)  seed_input
       write (6,*) 'after read seed_input...'
       call custom_flush(6)

       isendbuf(17)=seed_input

       i_dissenst=0
       isendbuf(18)=i_dissenst

       call custom_flush(6)

 5     continue

       close (lu)

    end if

    ! designate different random number sequences for different purposes

    luran1=19 ; luran2=20
    kranf1=5  ;     kranf2=5
    ksran=1


    ! PK Yeung, 1/25/10: make all tasks stop if task 0 encounters
    ! error in reading 'input'
 
    go to 30
 10 no_err=1
    go to 30
 20 no_end=1
 30 continue
    call MPI_BCAST (no_err,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    call MPI_BCAST (no_end,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpierr)
    if (no_err.gt.0) then
       if (taskid.eq.0) write (6,*) &
          'Code stops because of error in reading ''input'' file'
       stop
    end if
    if (no_end.gt.0) then
       if (taskid.eq.0) write (6,*) &
          'Code stops because of end-of-file in reading ''input'' file'
       stop
    end if


    call MPI_BCAST(isendbuf, 30, MPI_INTEGER, 0, MPI_COMM_WORLD, mpierr)
    call MPI_BCAST(rsendbuf, 30, mpireal, 0, MPI_COMM_WORLD, mpierr)


    if (taskid.ne.0) then

       nx=isendbuf(1)
       ny=isendbuf(2)
       nz=isendbuf(3)
       nc=isendbuf(4)
       nsteps=isendbuf(5)
       iostep=isendbuf(6)
       istart=isendbuf(7)
       kstop=isendbuf(8)
       isave=isendbuf(9)
       gpumpi=isendbuf(10)

       allocate (kinit(3+nc))

       kinit(1)=isendbuf(11)
       kinit(2)=kinit(1)
       kinit(3)=kinit(1)

       irz=isendbuf(13)
       kshift=isendbuf(14)
       seed_input=isendbuf(17)
       i_dissenst=isendbuf(18)

       viscos=rsendbuf(4)
       a11=rsendbuf(5); a22=rsendbuf(6)
       a33=rsendbuf(7); !!! shear=rsendbuf(8)
       beta1=rsendbuf(9)
       beta2=rsendbuf(10)
       beta3=rsendbuf(11)
       bet12=rsendbuf(12)

       cfl=rsendbuf(13)
       dt=rsendbuf(14)
       dtout=rsendbuf(15)
       tfout=rsendbuf(16)
       dtsave=rsendbuf(17)
       tfsave=rsendbuf(18)

    end if

    return
end
