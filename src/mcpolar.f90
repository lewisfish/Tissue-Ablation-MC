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
use utils

implicit none

integer          :: nphotons,iseed,j,xcell,ycell,zcell, N, counter, u
logical          :: tflag, flag, end
DOUBLE PRECISION :: nscatt, nscattGLOBAL
real             :: xmax,ymax,zmax, ran,delta, start, finish, ran2
real, allocatable :: tissue(:,:,:), temp(:,:,:), tissueGLOBAL(:,:,:)

! mpi variables
integer :: new_comm, error, right, left, id, numproc, dims(2), ndims, tag, recv_status(mpi_status_size), comm
logical :: periods(1), reorder

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array
call zarray

N = 200 ! points for heat sim
allocate(tissue(nxg, nyg, nzg), tissueGLOBAL(nxg,nyg,nzg))
allocate(temp(0:N+1, 0:N+1, 0:N+1))

temp = 0.
tissue = 0.
counter = 0
comm = MPI_COMM_WORLD

call MPI_init(error)

call MPI_Comm_size(comm, numproc, error)

!setup topology variables
tag = 1
dims = 0
periods = .false.
reorder = .true.
ndims = 1

!create cartesian topology on cpu
call mpi_dims_create(numproc, ndims, dims, error)
call mpi_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
call mpi_comm_rank(new_comm, id, error)

!get neighbours
call mpi_cart_shift(new_comm, 0, 1, left, right, error)

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

call init_opt1

if(id == 0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(xmax, ymax, zmax, id)
!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta = 1.e-8*(2.*zmax/nzg)
nscatt = 0
call MPI_Barrier(MPI_COMM_WORLD, error)
call cpu_time(start)
flag = .true.
end = .true.
!loop over photons 
print*,'Photons now running on core: ',id


do while(end)
   do j = 1, nphotons

      call init_opt1

      tflag=.FALSE.

      if(mod(j,100000) == 0)then
         print *, str(j)//' scattered photons completed on core: '//str(id)
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
   call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm,error)

   if(id==0)jmeanGLOBAL = jmeanGLOBAL * (1./(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
call heat_sim_3d(jmeanGLOBAL, tissue, temp, N, flag, id, numproc, error, new_comm, tag, recv_status, right, left, counter)

   counter = counter + 1
   if(counter == 1)end = .false.

   
   if(id == 0)then
      where(tissue> 1000)
            rhokap = 0.
      end where
   end if

   call MPI_Bcast(rhokap, size(rhokap), MPI_DOUBLE_PRECISION ,0 , new_comm, error)

   if(id == 0)then
      open(newunit=u,file=trim(fileplace)//'jmean/rhokap-'//str(counter)//'.dat',access='stream',form='unformatted')
      write(u)rhokap
   end if
   if(.not. end)exit
   jmean = 0.
end do
call MPI_Barrier(new_comm, error)
call cpu_time(finish)
if(finish-start.ge.60.)then
 print*,floor((finish-start)/60.)+mod(finish-start,60.)/100.
else
      print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm,error)

call MPI_REDUCE(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,new_comm,error)

if(id == 0)then
   print*,'Average # of scatters per photon:',(nscattGLOBAL/(nphotons*numproc))
   !write out files
   call writer(xmax, ymax, zmax, nphotons, numproc)
   print*,'write done'
end if

call MPI_Finalize(error)
end program mcpolar
