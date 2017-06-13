MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer(xmax,ymax,zmax,nphotons, numproc)

   use constants,only : nxg,nyg,nzg,fileplace
   use iarray,only : jmeanGLOBAL, rhokap

   implicit none

   integer           :: nphotons, i, u, numproc
   real              :: xmax, ymax, zmax !,jmeanf(nzg)
   ! character(len=20) :: fn
   
   

   
   !maybe not right
   jmeanGLOBAL = jmeanGLOBAL * (1./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))


   inquire(iolength=i)jmeanGLOBAL

   open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) jmeanGLOBAL
   close(u)

   inquire(iolength=i)rhokap

   open(newunit=u,file=trim(fileplace)//'jmean/rhokap.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) rhokap
   close(u)


   ! jmeanf=0.
   ! do i=1,nxg
   !    do j=1,nyg
   !       do k=1,nzg
   !          jmeanf(k) = jmeanf(k) + jmeanGLOBAL(i,j,k)
   !       end do
   !    end do
   ! end do
   ! jmeanf = jmeanf/(nxg*nyg)
   ! write(fn,"(i3.3,a,i3.3,a,i5.5,a)") nxg,'*',nyg,'*',nzg,'.dat'
   ! print*,trim(fn)
   ! open(newunit=u,file=trim(fn))
   ! do i=1,nzg        !for negative direction 1,nzg   +ive nzg,1,-1
   !    write(u,*)real(nzg-i)*(2./nzg),jmeanf(i)    !for -ive dir i-1*   for +ive nzg-i
   ! end do
   ! close(u)
   

   
   end subroutine writer
end MODULE writer_mod
