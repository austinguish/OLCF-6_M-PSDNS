      subroutine sptvar_omp (uny)

      use comp
      implicit none

      complex(b8) uny(xisz,ny,zjsz,3)
 
      real, allocatable :: ekz(:),eks(:)
 
      real(b8) term,term1,term2,term3,sum,rk2
      integer ik,lu

      real sqkx,sqky,sqkz,sqk,two9
        

 
      if (taskid.eq.0) write (6,*) 'enter sptvar, OMP-CPU version'
 
      allocate (ekz(nxh),eks(nxh))
 
      kx2(:)=kx(:)**2
      ky2(:)=ky(:)**2
      kz2(:)=kz(:)**2
 
      ky2(nyhp)=0.
      kz2(nzhp)=0.

      two9=2./9.
 
 
      ekz=0.
 
      !$OMP PARALLEL DO DEFAULT(NONE) &
      !$OMP PRIVATE (zp,z,sqkz,y,sqky,xp,x,tfact_x,term1,term2,term3,rk2,ik, &
      !$OMP term,sqk) &
      !$OMP SHARED (zjsz,zjst,nz,kx,ky,kz,kx2,ky2,kz2,mask,nyhp,nzhp,uny, &
      !$OMP nx,ny,xisz,xist,two9,taskid,kinit) &
      !$OMP REDUCTION(+:ekz)


      do 100 zp=1,zjsz
         z=zp+zjst-1
         sqkz=(kz(z)/nz)**2
 
         do 10 y=1,ny

            sqky=(ky(y)/ny)**2+sqkz

            do xp=1,xisz
               x=xp+xist-1
               if(x .eq. 1) then
                  tfact_x=1
               else
                  tfact_x=2
               endif
               sqk=sqky+(kx(x)/nx)**2
 
               if (sqk.gt.two9.or.y.eq.nyhp.or.z.eq.nzhp) go to 10
 
               if (mask(xp,y,zp).eq.0.) go to 10
 
               term1=real(uny(xp,y,zp,1)*conjg(uny(xp,y,zp,1)))
               term2=real(uny(xp,y,zp,2)*conjg(uny(xp,y,zp,2)))
               term3=real(uny(xp,y,zp,3)*conjg(uny(xp,y,zp,3)))
 
               rk2=kx2(x)+ky2(y)+kz2(z)
               ik=sqrt(rk2)+1.5

               term=tfact_x*(term1+term2+term3)*0.5
               if (kinit(1).eq.-2.and.term.gt.1.e-4) then
                  write (6,"('sptvar_omp,x,y,z,ik,term,taskid,zp,xp=',4i8,1p,e12.4,3i8)") &
                          x,y,z,ik,term,taskid,zp,xp
                  write (6,"('sptvar_omp: term1,term2,term3=',1p,3e12.4)") term1,term2,term3
               end if
 
               ekz(ik)=ekz(ik)+term

            end do
 
 10	 continue
 
 100  continue
 
      !$OMP END PARALLEL DO

      !call MPI_REDUCE (zsum,sum,1,mpireal, &
      !      MPI_SUM,0,MPI_COMM_WORLD,mpierr)
      call MPI_REDUCE (ekz,eks,nxh,mpireal, &
             MPI_SUM,0,MPI_COMM_WORLD,mpierr)
 
      if (taskid.eq.0) then
         if (istep.eq.0) open (20,file='spectrum')
         write (20,*) 'istep,kstep=',istep,kstep
         write (20,201) (eks(ik),ik=1,nxh)
 201     format ((1p,10e13.5))
         write (20,*) 'sum of spectrum=',sum(eks)
         call custom_flush(20)

         if (istep.eq.0) then
            open (file='ek0',newunit=lu)
            do ik=1,nxh
               write (lu,"('ik,e(k)=',i6,1p,e12.4)") ik,eks(ik)
            end do
         end if

      end if
 
      deallocate (ekz,eks)
 
      return
      end
