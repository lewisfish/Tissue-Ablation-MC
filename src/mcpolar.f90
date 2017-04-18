program mcpolar

use mpi

!shared data
use constants
use photon_vars
use iarray
use opt_prop

!subroutines
use subs
use gridset_mod
use sourceph_mod
use inttau2
use ch_opt
use stokes_mod
use writer_mod
use threedfinite

implicit none

integer          :: nphotons,iseed,j,xcell,ycell,zcell
logical          :: tflag
DOUBLE PRECISION :: nscatt, nscattGLOBAL
real             :: xmax,ymax,zmax, ran,delta, start, finish, ran2
integer          :: id, error, numproc

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray


call MPI_init(error)

call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)

call MPI_Comm_rank(MPI_COMM_WORLD, id, error)

!**** Read in parameters from the file input.params
open(10,file=trim(resdir)//'input.params',status='old')
   read(10,*) nphotons
   read(10,*) xmax
   read(10,*) ymax
   read(10,*) zmax
   read(10,*) n1
   read(10,*) n2
   close(10)

! set seed for rnd generator. id to change seed for each process
iseed=-95648324+id

!****** setup up arrays and bin numbers/dimensions

!***** Set up constants, pi and 2*pi  ********************************

iseed=-abs(iseed)  ! Random number seed must be negative for ran2

call init_opt4

if(id == 0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(xmax,ymax,zmax,id)
!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta = 1.e-8*(2.*zmax/nzg)
nscatt = 0


call cpu_time(start)

!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',id
! do j = 1, nphotons

!    call init_opt4

!    tflag=.FALSE.

!    if(mod(j,10000) == 0)then
!       print *, j,' scattered photons completed on core: ',id
!       ! where(absorb > 5.)!do ablation
!       !    rhokap = 0.
!       ! elsewhere(absorb == 2.)
!       !    rhokap = rhokap/2.
!       ! end where
!    end if
    
! !***** Release photon from point source *******************************
!    call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed,j)

! !****** Find scattering location

!    call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
      
! !******** Photon scatters in grid until it exits (tflag=TRUE) 
!    do while(tflag.eqv..FALSE.)
!       ran = ran2(iseed)
      
!       if(ran < albedo)then!interacts with tissue

!          call stokes(iseed)
!          nscatt = nscatt + 1

!       else
!          absorb(xcell, ycell, zcell) = absorb(xcell, ycell, zcell) + 1.
!          tflag = .true.
!          exit
!       end if

! !************ Find next scattering location
!       call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)          
!    end do




! end do      ! end loop over nph photons


! call cpu_time(finish)
! if(finish-start.ge.60.)then
!  print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
! else
!       print*, 'time taken ~',floor(finish-start/60.),'s'
! end if

! call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
! call MPI_BARRIER(MPI_COMM_WORLD, error)


! call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
! call MPI_BARRIER(MPI_COMM_WORLD, error)

! if(id == 0)then
!    print*,'Average # of scatters per photon:',(nscattGLOBAL/(nphotons*numproc))
!    !write out files
!    call writer(xmax,ymax,zmax,nphotons, numproc)
!    print*,'write done'
! end if


call MPI_BARRIER(MPI_COMM_WORLD,error)
if(id==0)then
   call heat2d()
end if
call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_Finalize(error)
end program mcpolar
