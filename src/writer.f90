MODULE writer_mod

implicit none
save

CONTAINS
   subroutine writer(xmax,ymax,zmax,nphotons, numproc)

   use constants, only : nxg,nyg,nzg,fileplace
   use iarray,    only : jmeanGLOBAL, rhokap

   implicit none

   integer :: nphotons, u, numproc
   real    :: xmax, ymax, zmax   

   !maybe not right
   jmeanGLOBAL = jmeanGLOBAL * (1./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))

   open(newunit=u,file=trim(fileplace)//'jmean/jmean.dat',access='stream',status='REPLACE',form='unformatted')
   write(u) jmeanGLOBAL
   close(u)


   open(newunit=u,file=trim(fileplace)//'jmean/rhokap.dat',access='stream',status='REPLACE',form='unformatted')
   write(u) rhokap
   close(u)
   end subroutine writer
end MODULE writer_mod
