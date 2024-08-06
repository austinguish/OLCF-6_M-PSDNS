module precision
   use ISO_C_BINDING
        
   include 'mpif.h'

#ifndef DOUBLE_PREC	 
   integer, parameter:: b8 = 4
   integer, parameter:: mpireal = MPI_REAL
   integer, parameter:: mpicomplex = MPI_COMPLEX
#else
   integer, parameter:: b8 = 8
   integer, parameter:: mpireal = MPI_REAL8
   integer, parameter:: mpicomplex = MPI_COMPLEX16
#endif
 
end module precision
 
!--------------------------------------------------
! global grid resolution in 3 dimensions 
 
module param
       
   use precision
 
   integer  :: zero,one
   integer  :: nx,ny,nz,nxh,nxhp,nyh,nyhp,nzh,nzhp,nu,nxhp2pad,nyhpad,nyhppad
   integer  :: nxpad,nxhpad,nxhppad,nypad,nzpad,nyxhppad,nzhpad,nzhppad
   integer  :: nut,nc,ncpp,ncp
   integer  :: nxhp2,nyhp2,nzhp2

   integer :: gpumpi, ipack, mvar

end module param
 
!----------------------------------------------------
! global variables related to MPI

module mpicom
   use param

     integer, public :: num_thr

   ! mpi process info
   integer ierr !, dims(3),  cartid(3)
   ! logical periodic(3),remain_dims(3)
   integer numtasks,taskid,mpierr
   integer iproc, jproc, ipid, jpid

 
   integer, dimension(:), allocatable :: iist,iien,iisz
   integer, dimension(:), allocatable :: jjst,jjen,jjsz
   integer, dimension(:), allocatable :: kist,kien,kisz
   integer, dimension(:), allocatable :: kistpad,kienpad,kiszpad
   integer, dimension(:), allocatable :: kjst,kjen,kjsz
   integer, save :: mystart_x,mystart_z

   integer :: xist,xjst,xisz,xjsz,xien,xjen
   integer :: yist,yjst,yisz,yjsz,yien,yjen
   integer :: zist,zisz,zien
   integer :: zjst,zjsz,zjen
   integer :: num_al,max_al_x
   integer, allocatable :: max_al_z(:),cut_z(:),cut_z_i(:,:)
   integer, allocatable :: num_al_i(:)

   integer mpi_comm_cart
   integer mpi_comm_row, mpi_comm_col

   integer,dimension(:,:),allocatable:: status
   integer,dimension(:), allocatable:: mymap,inverse_map

   integer, allocatable :: ipid_all(:),jpid_all(:)

   complex(b8), allocatable, target :: buf(:,:)
   complex(b8), allocatable, target :: sndbuf_single(:)
   complex(b8), allocatable, target :: rcvbuf_single(:)

   ! time to get creative with the send/receive buffers
   type(c_ptr), target :: sndbuf_buffer = c_null_ptr
   type(c_ptr), target :: rcvbuf_buffer = c_null_ptr
   complex(b8), pointer, dimension(:,:,:,:,:) :: sndbuf_row
   complex(b8), pointer, dimension(:,:,:,:,:) :: rcvbuf_row
   complex(b8), pointer, dimension(:,:,:,:,:) :: sndbuf_col
   complex(b8), pointer, dimension(:,:,:,:,:) :: rcvbuf_col

   integer(c_size_t) :: xiszST,ziszST,yjszST,zjszST,nvST
   integer(c_size_t) :: iprocST,jprocST

   ! IO file units
   integer :: io_unit_eulstat, io_unit_escout
   integer :: io_unit_corr1d, io_unit_sptr1d, io_unit_msqvg
   integer :: io_unit_maxepsenst, io_unit_pdfepsenst, io_unit_pdfvelgrad
   integer :: io_unit_maxchi, io_unit_pdfchi
   integer :: io_unit_scgrad, io_unit_scgsptr
   integer :: io_unit_scflxsp, io_unit_scgpll1, io_unit_scgppp1
   integer :: io_unit_vgm, io_unit_vgskew, io_unit_scgm
   integer :: io_unit_forcinp

#ifdef DETAIL_MPI
   integer mpistat_groupsize,mpistat_ngroups,mpi_detail
   integer mpistat_dims(2),ipid_stat,jpid_stat,stat_cartid(2)
   integer mpi_comm_stat
   integer mpi_comm_stat_within
   integer mpi_comm_stat_across
   integer, allocatable :: itask_grp(:,:)
   real(8), allocatable :: rt_a2a_grp(:,:,:)
#endif

   contains

      subroutine custom_flush(n)
         implicit none
         integer :: n

#ifdef CALL_FLUSH
         call flush(n)
#else
         flush(n)
#endif

      end subroutine custom_flush

end module mpicom
 
!--------------------------------------------
 
! variables corresponding to old "com"
 
module com
 
   use mpicom

   ! rwall0 is initialized to clock value when code starts

   real(8) rwall0
 
   integer :: x,y,z,xp,yp,zp
   integer :: is1,is2,ipsi
 

   ! input/output logical units and filenames
   integer, allocatable :: luinit(:),kinit(:)
   character(8), allocatable ::  fninit(:)


   ! wavenumber components
   real(b8), allocatable :: kx(:),ky(:),kz(:)
   integer :: kmax


   complex(b8) imagi

   integer :: istep,nsteps,istep0, iostep,jstep,kstep
   integer :: istop,isave,ioflag,kstop
   integer :: istart,isflag
   integer :: rkstages
   real(b8) :: cfl,dt,entime,dtout,dtsave,tfout,tfsave,time,time0

   ! grid metric factors and grid shifts

   real(b8) :: beta1, beta2, beta3, bet12, pi
   real(b8) :: b11(2),b22(2),b33(2),b12(2),gsh(3,4)
   real(b8), allocatable :: bk1(:,:),bk2(:,:), bk3(:,:), bkk3(:,:)
   complex(kind=b8), allocatable :: bk1i(:,:), bk3i(:,:), bk2i(:,:)
   integer :: kshift
   complex(kind=b8), allocatable :: sx(:,:),sy(:,:),sz(:,:)

   ! viscosity, mean velocity gradients, Schmidt numbers
   real(b8) :: viscos
   real(b8), allocatable :: sc(:)
   real(b8) :: a11, a22, a33
   real(kind=b8) :: klen,kvel,ktime,ekl
 
   real(b8) velmax, dt2


   integer :: lustr,luran1,luran2

   ! random numbers stuff
   integer :: seed_input
   integer :: ksran,kranf1,kranf2

   integer irz
   character*256 indir_fn

   ! to determine endianness of the file to read in with
   integer endian

   type(c_ptr) :: hip_plan_r2c, hip_plan_c2r, &
                  hip_plan_z, hip_plan_y

   ! for work buffers
   integer :: use_work_buffer = 1
   integer :: hipfft_auto_alloc_mem = 0
   integer(kind=c_size_t),parameter :: zeroST = 0_c_size_t
   integer(kind=c_size_t), target :: wbufsize_r2c = zeroST
   integer(kind=c_size_t), target :: wbufsize_z = zeroST
   integer(kind=c_size_t), target :: wbufsize_y = zeroST
   integer(kind=c_size_t), target :: wbufsize_c2r = zeroST
   integer(kind=c_size_t) :: wbufsize_single = zeroST
   type(c_ptr) :: wbuffer_r2c = c_null_ptr
   type(c_ptr) :: wbuffer_z = c_null_ptr
   type(c_ptr) :: wbuffer_y = c_null_ptr
   type(c_ptr) :: wbuffer_c2r = c_null_ptr
   type(c_ptr) :: wbuffer_single = c_null_ptr

   ! integer to determine rkmethod
   integer(4) :: rkmethod

   ! index for numerical sequence of checkpoints written by the code
   integer ichkpt

   ! Forcing parameters
   integer forc_type
   ! FEK_FORC parameters
   integer kf_shell
   real, allocatable ::  ek_nf1(:),ek_nf2(:),ekinf(:)
   real, allocatable ::  ekinf_input(:)
   real, allocatable ::  fekz(:),feks(:)

   ! EPFOR parameters
   real(b8) :: tforce, epsfor, kforce, tfiuo
   real(b8), allocatable :: suo(:),svo(:)
   integer  nrproc

   integer :: kfor,k2fo
   complex(b8), allocatable :: for(:,:,:,:)
   complex(b8), allocatable :: velinc(:,:,:,:)

   real(b8), allocatable :: efki(:,:),efkz(:,:),efk(:)
   real(b8) :: eirate(3),erate



   ! i_dissent=1 if dissent routine is to be called (at output steps only)
   integer i_dissenst

   ! i_press=1 if pressure is to be calculated (at output steps only)
   integer i_press

end module com
 
!-----------------------------------------------------------
 
! variables corresponding to old "comp"
 
module comp

   use com
 
   ! the main field arrays
 
   integer mxyz
   real(b8), allocatable :: un(:,:,:,:)
   real(b8), allocatable :: u(:,:,:,:)

   ! next array used only if i_press=1
   real(b8), allocatable :: up(:,:,:)
 
   ! mask for truncation of aliased modes
 
   real(kind=b8), allocatable :: mask(:,:,:)


   real(b8) :: tke,epslon, tfact_x
   real(b8), allocatable :: ek(:),dk(:),sk(:)
   real(b8), allocatable :: kx2(:),ky2(:),kz2(:)
   real(b8), allocatable :: tfact(:)

   real(kind=b8), allocatable :: corr(:)
   real(kind=b8) :: skew,ett(3),raak(3),tmre(3)
   real(kind=b8), allocatable :: rms(:),tmij(:,:)
   real(kind=b8), allocatable :: lijk(:,:), laak(:,:), taylor(:,:,:),&
      taak(:,:)
   ! store indices for correlation terms needed in sptr
   integer, allocatable :: kt(:)
   real(kind=b8), allocatable :: ekk(:)
   real(kind=b8), allocatable :: eijk(:,:),dijk(:,:),sijk(:,:)
   real(kind=b8), allocatable :: xrij(:,:),yrij(:,:),zrij(:,:)

   ! for scalars
   real(kind=b8), allocatable :: pr(:),grad(:,:),cb(:,:)
   real(kind=b8), allocatable :: scmean(:)
   real(kind=b8), allocatable :: scdiss(:),sctime(:),mtosratio(:)
   real(kind=b8), allocatable :: batchelor(:),obukc(:)
   real(kind=b8), allocatable :: spm4(:),spm6(:)
   ! scalar gradients
   real(kind=b8), allocatable :: grijk(:,:,:)
   real(kind=b8), allocatable :: scgmsq(:,:),scgcor(:,:),scgcov(:,:),&
      scgvar(:,:)

   ! for 1d correlation
   complex(kind=b8), allocatable :: qijx(:,:),qijy(:,:),qijz(:,:)

   real(b8), allocatable :: xfac(:),yfac(:),zfac(:)
   real(b8), allocatable :: xdif(:,:),ydif(:,:),zdif(:,:)

end module comp


! random_number_seed ------------------------------------------
module ranseed
  integer :: iseed,seed_size
  integer, allocatable :: rseed(:)
end module ranseed

 
! -------------------------------------	

module timers_comm
   real(8) :: t1_comm,t2_comm,t3_comm,t4_comm,t4t_comm,tp1_comm
   real(8) :: gt1_comm,gt2_comm,gt3_comm,gt4_comm,gt4t_comm,gtp1_comm
   integer :: io_mpi, io_mpi_grp
   real(8) :: t_alltoall, t_a2a, rt_a2a(4,5)
   integer :: i1_comm, i2_comm, i3_comm, i4_comm, i4t_comm, ip_comm
end module timers_comm

module timers_comp
   real(8) :: tcpu_fft,tcpu_other
end module timers_comp
 
module timers_io
   integer :: iread_io,iwrite_io
   real(8) :: tread_io,twrite_io
end module timers_io

module timers_tran
   real(8) t_xk(4),t_kx(4)
   real(8) t_xkcomm1(4),t_xkcomm2(4)
   real(8) t_kxcomm1(4),t_kxcomm2(4)

   real(8) t_pack(4),t_unpack(4)
   real(8) tot_pack(4),tot_unpack(4)
end module timers_tran

module timers_rkstep
   integer ncpusteps,ncpusteps_io
   real(8) t_rks(10)
   real(8) t_itrans(4)
   real(8) t_trans(4)
end module timers_rkstep

module timers_fom
   integer ifom
   real(8) fom_problem_size, fom_avg_time, computed_fom
   real(8) fom_fft_time, fom_pack_time, fom_comm_time, fom_other_time
end module timers_fom

