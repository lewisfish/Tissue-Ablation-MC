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
use Heat, only : power, delt, energyPerPixel, laser_flag, laserOn, loops, pulseCount, pulsesToDo,pulseFlag,getPWr,&
                 repetitionCount, repetitionRate_1, time, total_time, initThermalCoeff, heat_sim_3d,setupThermalCoeff, arrhenius&
                 ,realpulseLength
use utils

implicit none

integer           :: nphotons ,iseed, j, xcell, ycell, zcell, N, counter, u, q, w, e
logical           :: tflag
double precision  :: nscatt
real              :: xmax, ymax, zmax, delta, start, finish, ablateTemp
real, allocatable :: temp(:,:,:), tissue(:,:,:), ThresTime(:,:,:,:), ThresTimeGLOBAL(:,:,:,:)

! mpi variables
type(mpi_comm)   :: comm, new_comm
integer          :: right, left, id, numproc, dims(2), ndims, tag
logical          :: periods(1), reorder
character(len=3) :: gfd

!start MPI
call MPI_init()
comm = MPI_COMM_WORLD
call MPI_Comm_size(comm, numproc)


!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array(numproc)
call zarray

N = nzg ! points for heat sim
allocate(tissue(nxg, nyg, nzg))
allocate(temp(0:N+1, 0:N+1, 0:N+1))
allocate(ThresTime(nxg, nyg, nzg, 3))

tissue = 0.d0
ThresTime = 0.d0!

time            = 0.
pulseCount      = 0.
repetitionCount = 0.
laserOn         = 1.
laser_flag      = .TRUE.
pulseFlag       = .FALSE.

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
   read(u,*) total_time
   read(u,*) loops
   read(u,*) repetitionRate_1
   read(u,*) power
   read(u,*) energyPerPixel
   read(u,*) ablateTemp
   read(u,*) pulsesToDo
   close(u)

! set seed for rnd generator. id to change seed for each process
iseed = -95648324 + id
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

!barrier to ensure correct timingd
call MPI_Barrier(MPI_COMM_WORLD)
call cpu_time(start)

!loop over photons 
print*,'Photons now running on core: ',id


temp          = 5.d0 + 273.d0
temp(N+1,:,:) = 5.+273.  ! side face
temp(0,:,:)   = 5.+273.    ! side face
temp(:,0,:)   = 5.+273. ! front face
temp(:,N+1,:) = 5.+273.  ! back face
temp(:,:,0)   = 25.+273.  ! bottom face
temp(:,:,N+1) = 25.+273.  ! top face 
call initThermalCoeff(delt, N, xmax, ymax, zmax)

if(id == 0)print*,energyPerPixel,int(total_time/delt),realpulselength,delt,int(realpulselength/delt),total_time

!override total_time if set too low for laser to finish 1 pulse
if(int(total_time/delt) <= int(realpulseLength/delt))then
   total_time = delt * (realpulselength/delt + 2000.)
end if

do while(time <= total_time)
   if(laser_flag)then

      do j = 1, nphotons

         tflag=.FALSE.

         if(mod(j,1000000) == 0)then
            print *, str(j)//' scattered photons completed on core: '//str(id)
         end if
          
      !***** Release photon *******************************
         call sourcephCO2(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

      !****** Find scattering/absorb location
         call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
            
      !******** Photon scatters in grid until it exits (tflag=TRUE) 
         do while(tflag.eqv..FALSE.)
            tflag = .true.
            exit
         end do
      end do      ! end loop over j photons

      !reduce jmean from all processess
      jmeanGLOBAL = 0.d0
      call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm)

      jmeanGLOBAL = jmeanGLOBAL * ((getPwr()/81.d0)/(nphotons*numproc*(2.*xmax*1.d-2/nxg)*(2.*ymax*1.d-2/nyg)*(2.*zmax*1.d-2/nzg)))
   end if

   !do heat simulation
   call heat_sim_3d(jmeanGLOBAL, temp, N, id, numproc, new_comm, right, left, counter)

   call arrhenius(temp, delt, tissue, ThresTime, 1, N, N)
   !update thermal/optical properties
   call setupThermalCoeff(temp, N, ablateTemp)

   counter = counter + 1
   jmean = 0.
end do
   
   allocate(ThresTimeGLOBAL(nxg, nyg, nzg, 3))
   ThresTimeGLOBAL = 0.d0
   call MPI_REDUCE(ThresTime, ThresTimeGLOBAL, nxg*nyg*nzg*3, MPI_DOUBLE_PRECISION, mpi_min, 0, new_comm)

   !write out results
   if(id == 0)then
      !delete damage info about ablation crater
      !as no tissue left to damage!
      do q = 1, nxg
         do w = 1, nyg
            do e = 1, nzg
               if(rhokap(q,w,e) <= 0.1)then
                  tissue(q,w,e) = -1.d0
               end if
            end do
         end do
      end do
      call writer(ablateTemp, temp, tissue, xmax, ymax, zmax, ThresTimeGlobal)
   end if

call cpu_time(finish)
if(finish-start.ge.60.)then
    print*,str(floor((finish-start)/60.)+mod(finish-start,60.)/100.,5)//' mins'
else
    print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_Finalize()
end program mcpolar
