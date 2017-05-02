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
use Heat

implicit none

integer          :: nphotons,iseed,j,xcell,ycell,zcell, N, counter
logical          :: tflag, flag
DOUBLE PRECISION :: nscatt, nscattGLOBAL
real             :: xmax,ymax,zmax, ran,delta, start, finish, ran2
real, allocatable :: tissue(:,:,:)
integer          :: id, error, numproc

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray

N = 50 ! points for heat sim
allocate(tissue(N,N,N))
tissue = 0.
counter = 0.
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
flag = .true.
!loop over photons 
call MPI_Barrier(MPI_COMM_WORLD, error)
print*,'Photons now running on core: ',id
do while(flag)
   do j = 1, nphotons

      call init_opt4

      tflag=.FALSE.

      if(mod(j,10000) == 0)then
         print *, j,' scattered photons completed on core: ',id
      end if
       
   !***** Release photon from point source *******************************
      call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed,j)

   !****** Find scattering location

      call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
         
   !******** Photon scatters in grid until it exits (tflag=TRUE) 
      do while(tflag.eqv..FALSE.)
         ran = ran2(iseed)
         
         if(ran < albedo)then!interacts with tissue

            call stokes(iseed)
            nscatt = nscatt + 1

         else
            tflag = .true.
            exit
         end if

   !************ Find next scattering location
         call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)          
      end do




   end do      ! end loop over nph photons
   call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
   call MPI_BARRIER(MPI_COMM_WORLD, error)

   if(id == 0)then
      jmeanGLOBAL = jmeanGLOBAL * (1./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
      call heat_sim_3d(jmeanGLOBAL, tissue, N, counter)

      counter = counter + 1
      if(counter == 10)flag = .false.
   end if
   call MPI_BARRIER(MPI_COMM_WORLD, error)
   call mpi_bcast(flag,1,MPI_LOGICAL, 0, MPI_COMM_WORLD, error)
   if(.not. flag)exit
   jmean = 0.
end do

call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)


call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,error)
call MPI_BARRIER(MPI_COMM_WORLD, error)

if(id == 0)then
   print*,'Average # of scatters per photon:',(nscattGLOBAL/(nphotons*numproc))
   !write out files
   call writer(xmax,ymax,zmax,nphotons, numproc)
   print*,'write done'
end if

call MPI_BARRIER(MPI_COMM_WORLD, error)

call MPI_Finalize(error)
end program mcpolar
