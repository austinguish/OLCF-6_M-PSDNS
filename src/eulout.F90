subroutine eulout
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! output of Eulerian velocity statistiscs (calculated by sub. sptr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   use comp
   implicit none       
   integer lu,i,j
   integer :: k,ios
   integer :: nshell
   real(kind=b8) :: shell

   !!!!!!!!!!!!!!!!!!!!!!!! eulstat !!!!!!!!!!!!!!!!!!!!!!!!!!

   ! set logical unit to the already open eulstat file
   lu=io_unit_eulstat
   open(lu,file='eulstat',position='append')

   if (taskid.ne.0) return

   write (lu,210) istep,time,sqrt(b11(2)),sqrt(b22(2)),sqrt(b33(2))

   write (lu,211) 'u',(laak(1,j),j=1,3),ekl
   write (lu,211) 'v',(laak(2,j),j=1,3)
   write (lu,211) 'w',(laak(3,j),j=1,3)

   write (lu,212) (taak(i,i),i=1,3)
   write (lu,213) (raak(i),i=1,3)
   write (lu,214) (tmre(i),i=1,3),(tmre(1)+tmre(2)+tmre(3))/3.
   write (lu,215) (rms(i),i=1,3),corr(2)

   write (lu,216) tke
   write (lu,217) epslon, kmax*klen
   write (lu,218) skew
   write (lu,219) klen,kvel,ktime
   write (lu,231) (ett(i),i=1,3)

   shell=min(beta1,beta2,beta3)
   nshell=max(nxh,nyh,nzh)

   write (lu,2101) shell
   write (lu,220) (ek(k),k=1,nshell)
   write (lu,2102)
   write (lu,220) (dk(k),k=1,nshell)

   close(lu)

   !!!!!!!!!!!!!!!!!!!!!!!! msqvelgrad !!!!!!!!!!!!!!!!!!!!!!!!!!

   ! set logical unit to the already open eulstat file
   lu=io_unit_msqvg
   open(lu,file='msqvelgrad',position='append')

   if(istep.eq.0) then
      write (lu,244)
   end if

   write(lu,245) time,((tmij(i,j),j=1,3),i=1,3)


   close(lu)



210  format (/'istep=',i6,2x,'time=',1p,e13.5,'  betas=',1p,3e13.5)
211  format (a1,'integral length scales:',1p,4e12.4)
212  format (' taylor microscales    :',1p,3e12.4)
213  format (' int. scale reynolds no:',1p,3e12.4)
214  format (' taylor reynolds no.   :',1p,4e12.4)
215  format (' rms velocities & uv   :',1p,4e12.4)
216  format (' turb. kinetic energy  :',1p,e12.4)
217  format (' dissipation rate      :',1p,e12.4,2x,'kmax.eta=',1p,e12.4)
218  format (' dissipation skewness  :',1p,e12.4)
219  format (' kol. scales (l,v,t)   :',1p,3e12.4)
220  format ((3x,1p,10e13.5))
231  format (' eddy turnover times   :',1p,3e12.4)
2101 format ('3-d energy spectrum, shell thickness=',1p,e12.4)
2102 format ('3-d dissipation spectrum')

244 format ('mean-squares of velocity gradients'/14x,'   du/dx',&
   '      du/dy      du/dz      dv/dx      dv/dy   ',&
   '   dv/dz      dw/dx      dw/dy      dw/dz')
245 format ('t=',1p,e11.4,1x,1p,9e11.4)

   return 
end
