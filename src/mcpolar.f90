program mcpolar

use mpi_f08

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

integer           :: nphotons ,iseed, j, xcell, ycell, zcell, N, counter, u, iters,i,k
logical           :: tflag, flag, end
double precision  :: nscatt, nscattGLOBAL
real              :: xmax, ymax, zmax, ran, delta, start, finish, ran2, power
real, allocatable :: tissue(:,:,:), temp(:,:,:), tissueGLOBAL(:,:,:), q(:,:,:)

! mpi variables
type(mpi_comm)   :: comm, new_comm
type(MPI_Status) :: recv_status
integer          :: right, left, id, numproc, dims(2), ndims, tag
logical          :: periods(1), reorder


call MPI_init()
comm    = MPI_COMM_WORLD
call MPI_Comm_size(comm, numproc)


!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array(numproc)
call zarray

N = 104 ! points for heat sim
allocate(tissue(nxg, nyg, nzg), tissueGLOBAL(nxg,nyg,nzg), q(nxg,nyg,nzg))
allocate(temp(0:N+1, 0:N+1, 0:N+1))

Q = 0.


temp    = 0.
tissue  = 0.
counter = 0
!setup topology variables
tag     = 1
dims    = 0
periods = .false.
reorder = .true.
ndims   = 1

!create cartesian topology on cpu
call mpi_dims_create(numproc, ndims, dims)
call mpi_cart_create(comm, ndims, dims, periods, reorder, new_comm)
call mpi_comm_rank(new_comm, id)

!get neighbours
call mpi_cart_shift(new_comm, 0, 1, left, right)

!**** Read in parameters from the file input.params
open(newunit=u,file=trim(resdir)//'input.params',status='old')
   read(u,*) nphotons
   read(u,*) xmax
   read(u,*) ymax
   read(u,*) zmax
   read(u,*) n1
   read(u,*) n2
   read(u,*) iters
   read(u,*) power
   close(u)

! set seed for rnd generator. id to change seed for each process
iseed = -95648324 + id

!****** setup up arrays and bin numbers/dimensions

!***** Set up constants, pi and 2*pi  ********************************

iseed = -abs(iseed)  ! Random number seed must be negative for ran2

call init_opt1

if(id == 0)then
   print*, ''      
   print*,'# of photons to run',nphotons*numproc
end if

!***** Set up density grid *******************************************
call gridset(xmax, ymax, zmax, id)
!***** Set small distance for use in optical depth integration routines 
!***** for roundoff effects when crossing cell walls
delta  = 1.e-8*(2.*zmax/nzg)
nscatt = 0
call MPI_Barrier(MPI_COMM_WORLD)
call cpu_time(start)
flag = .true.
end = .true.
!loop over photons 
print*,'Photons now running on core: ',id


temp = 37 + 273.
temp(:,:,size(temp,3)-1) = 25. + 273.

do while(end)
   do j = 1, nphotons

      call init_opt1

      tflag=.FALSE.

      if(mod(j,1000000) == 0)then
         print *, str(j)//' scattered photons completed on core: '//str(id)
      end if
       
   !***** Release photon from point source *******************************
      call sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

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
   call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm)

   jmeanGLOBAL = jmeanGLOBAL * (power/(nphotons*numproc*(2.*xmax/nxg)*(2.*ymax/nyg)*(2.*zmax/nzg)))
   if(id == 0)then
      print*,counter
   end if
call heat_sim_3d(jmeanGLOBAL, tissue, temp, N, flag, id, numproc, new_comm, right, left, counter, Q)

   counter = counter + 1
   if(counter == iters)end = .false.

   ! if(id == 0)then
      where(tissue >= 6.)
         rhokap = 0.  
      end where
   ! end if

   ! call MPI_Bcast(rhokap, size(rhokap), MPI_DOUBLE_PRECISION ,0 , new_comm)

   if(id == 0)then
      open(newunit=u,file=trim(fileplace)//'jmean/rhokap-'//str(counter-1)//"-"//str(power,3)//'.dat',access='stream',&
           form='unformatted',status='replace')
      write(u)rhokap
      close(u)
   end if
   if(.not. end)exit
   jmean = 0.
end do

call cpu_time(finish)
if(finish-start.ge.60.)then
    print*,str(floor((finish-start)/60.)+mod(finish-start,60.)/100.,5)//' mins'
else
    print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call mpi_reduce(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm)

call mpi_reduce(nscatt,nscattGLOBAL,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,new_comm)

if(id == 0)then
   print*,'Average # of scatters per photon:',(nscattGLOBAL/(nphotons*numproc))
   call writer(xmax, ymax, zmax, nphotons, numproc)
   print*,'write done'
end if

call MPI_Finalize()
end program mcpolar
