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
use writer_mod
use stokes_mod
use Heat, only : power, delt, energyPerPixel, laser_flag, laserOn, loops, pulseCount, pulsesToDo,pulseFlag,getPWr,&
                 repetitionCount, repetitionRate_1, time, total_time, initThermalCoeff, heat_sim_3d,realpulseLength
use utils
use memoryModule, only : checkallocate

implicit none

integer           :: nphotons ,iseed, j, xcell, ycell, zcell, N, counter, u, q, w, e
logical           :: tflag
double precision  :: nscatt
real              :: xmax, ymax, zmax, delta, start, finish, ablateTemp, ran, ran2

! mpi variables
type(mpi_comm)   :: comm, new_comm
integer          :: right, left, id, numproc, dims(2), ndims, tag
logical          :: periods(1), reorder

!start MPI
call MPI_init()
comm = MPI_COMM_WORLD
call MPI_Comm_size(comm, numproc)


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

!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array(numproc, id)
call zarray

N = nzg ! points for heat sim

tissue = 0.d0

time            = 0.
pulseCount      = 0.
repetitionCount = 0.
laserOn         = 1.
laser_flag      = .TRUE.
pulseFlag       = .FALSE.

counter = 0

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
   read(u,*) total_time
   read(u,*) loops
   read(u,*) repetitionRate_1
   read(u,*) power
   read(u,*) energyPerPixel
   read(u,*) ablateTemp
   read(u,*) pulsesToDo
   read(u, *) pulsetype
   close(u)

! set seed for rnd generator. id to change seed for each process
iseed = -95648324 + id
iseed = -abs(iseed)  ! Random number seed must be negative for ran2


! call init_opt4

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

!barrier to ensure correct timingd
call MPI_Barrier(MPI_COMM_WORLD)
call cpu_time(start)

!loop over photons 
print*,'Photons now running on core: ',id


temp          = 37.d0 + 273.d0
temp(N+1,:,:) = 37.+273.  ! side face
temp(0,:,:)   = 37.+273.    ! side face
temp(:,0,:)   = 37.+273. ! front face
temp(:,N+1,:) = 37.+273.  ! back face
temp(:,:,0)   = 37.+273.  ! bottom face
temp(:,:,N+1) = 37.+273.  ! top face 
call initThermalCoeff(delt, N, xmax, ymax, zmax, numproc)


if(id == 0)then
   print*,"Energy,      total loops,  realpulse length,    delt,            laser on,   total sim time"
   print'(F8.1,1x,I10.5,9x,F5.3,16x,E11.5,1x,I10.5,5x,F7.3)',energyPerPixel,int(total_time/delt),realpulselength,delt,&
                                                            int(realpulselength/delt),total_time
end if

do while(time <= total_time)
   if(laser_flag)then

      do j = 1, nphotons
         tflag=.FALSE.
          
      !***** Release photon *******************************
         call sourcephCO2(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

      !****** Find scattering/absorb location
         call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
            
      !******** Photon scatters in grid until it exits (tflag=TRUE) 
         do while(tflag.eqv..FALSE.)
            ran = ran2(iseed)
            
            if(ran < albedo(zcell))then!interacts with tissue
               call stokes(iseed)
               nscatt = nscatt + 1
            else
               tflag = .true.
               exit
            end if
      !************ Find next scattering location
            call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)          
         end do
      end do      ! end loop over j photons

      !reduce jmean from all processess
      call MPI_allREDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,new_comm)
      jmeanGLOBAL = jmeanGLOBAL * ((getPwr()/81.d0)/(nphotons*numproc*(2.*xmax*1.d-2/nxg)*(2.*ymax*1.d-2/nyg)*(2.*zmax*1.d-2/nzg)))
   end if

   !do heat simulation
   call heat_sim_3d(jmeanGLOBAL, temp, N, id, numproc, new_comm, right, left, counter)

   counter = counter + 1
   jmean = 0.
   if(id == 0)print*,counter
end do
   
   !write out results
   if(id == 0)then
      call writer(ablateTemp, temp, tissue, xmax, ymax, zmax)
   end if

call cpu_time(finish)
if(finish-start.ge.60.)then
    print*,str(floor((finish-start)/60.)+mod(finish-start,60.)/100.,5)//' mins'
else
    print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_Finalize()
end program mcpolar
