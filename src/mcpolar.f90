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
use Heat, only : power, delt, energyPerPixel, laser_flag, laserOn, loops, pulseCount, pulselength, pulsesToDo&
                 , repetitionCount, repetitionRate_1, time, total_time, initThermalCoeff, heat_sim_3d, watercontent,loops_left,q,&
                 setupThermalCoeff
use utils

implicit none

integer           :: nphotons ,iseed, j, xcell, ycell, zcell, N, counter, u
logical           :: tflag
double precision  :: nscatt
real              :: xmax, ymax, zmax, delta, start, finish, ablateTemp
real, allocatable :: temp(:,:,:), tissue(:,:,:), tissueGLOBAL(:,:,:)

! mpi variables
type(mpi_comm)   :: comm, new_comm
integer          :: right, left, id, numproc, dims(2), ndims, tag,r,w,e
logical          :: periods(1), reorder, ablateFlag=.FALSE.


call MPI_init()
comm    = MPI_COMM_WORLD
call MPI_Comm_size(comm, numproc)


!set directory paths
call directory

!allocate and set arrays to 0
call alloc_array(numproc)
call zarray

N = nzg ! points for heat sim
allocate(tissue(nxg, nyg, nzg), tissueGLOBAL(nxg,nyg,nzg))!, q(nxg,nyg,nzg))
allocate(temp(0:N+1, 0:N+1, 0:N+1))


time            = 0.
pulseCount      = 0.
repetitionCount = 0.
laserOn         = 1.
laser_flag      = .TRUE.

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
   read(u,*) energyPerPixel
   read(u,*) ablateTemp
   read(u,*) pulsesToDo
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
!loop over photons 
print*,'Photons now running on core: ',id

! pulsesDone=0
! open(newunit=u,file="/home/lewis/phdshizz/ablation/bin/temp-mid.dat",access='stream',form='unformatted')
! read(u)temp(0:n+1,0:n+1,0:n+1)
! close(u)

! temp = temp + 273.d0
! call initThermalCoeff(delt, N, xmax, ymax, zmax)


! open(newunit=u,file="/home/lewis/phdshizz/ablation/bin/rhokap-mid.dat",access='stream',form='unformatted')
! read(u)rhokap
! close(u)

! open(newunit=u,file="/home/lewis/phdshizz/ablation/bin/water-mid.dat",access='stream',form='unformatted')
! read(u)watercontent
! close(u)

! open(newunit=u,file="/home/lewis/phdshizz/ablation/bin/q-mid.dat",access='stream',form='unformatted')
! read(u)q
! close(u)

temp = 5.d0 + 273.d0
temp(N+1,:,:) = 5.+273.  ! side face
temp(0,:,:) = 5.+273.    ! side face
temp(:,0,:) = 5.+273. ! front face
temp(:,N+1,:) = 5.+273.  ! back face
temp(:,:,0) = 25.+273.  ! bottom face
temp(:,:,N+1) = 25.+273.  ! top face 
call initThermalCoeff(delt, N, xmax, ymax, zmax)


print*,energyPerPixel,int(total_time/delt),pulselength,delt,int(pulselength/delt)

do while(time <= total_time)
   if(laser_flag)then

      do j = 1, nphotons

         call init_opt1

         tflag=.FALSE.

         if(mod(j,1000000) == 0)then
            print *, str(j)//' scattered photons completed on core: '//str(id)
         end if
          
      !***** Release photon from point source *******************************
         ! call sourcephERYAG(xmax,ymax,zmax,xcell,ycell,zcell,iseed)
         call sourcephCO2(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

      !****** Find scattering location
         call tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
            
      !******** Photon scatters in grid until it exits (tflag=TRUE) 
         do while(tflag.eqv..FALSE.)
               tflag = .true.
               exit
         end do
      end do      ! end loop over nph photons
      jmeanGLOBAL = 0.d0
      call MPI_REDUCE(jmean, jmeanGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,0,new_comm)

      jmeanGLOBAL = jmeanGLOBAL * ((power/81.)/(nphotons*numproc*(2.*xmax*1.d-2/nxg)*(2.*ymax*1.d-2/nyg)*(2.*zmax*1.d-2/nzg)))
   end if
   call heat_sim_3d(jmeanGLOBAL, temp, tissue, N, id, numproc, new_comm, right, left, counter)


   ! tissueGLOBAL = 0.

   ! call MPI_allREDUCE(tissue, tissueGLOBAL, (nxg*nyg*nzg),MPI_DOUBLE_PRECISION, MPI_SUM,new_comm)

   do r = 1, N
      do w= 1, N
         do e = 1, N
            if(temp(r,w,e) >= ablateTemp + 273.d0)then
               ablateFlag = .TRUE.
               rhokap(r,w,e) = 0.d0
               temp(r,w,e) = 273.d0+25.d0
            end if
         end do
      end do
   end do

  call setupThermalCoeff(temp, N)

   counter = counter + 1
   jmean = 0.
end do
   if(id == 0)then

         open(newunit=u,&
      file=trim(fileplace)//"ErYAG/jmean-timestep-deltdiv100-"&
           //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
      write(u)jmeanGLOBAL
      close(u)

      open(newunit=u,&
      file=trim(fileplace)//"ErYAG/rhokap-timestep-deltdiv100-"&
           //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
      write(u)rhokap(1:nxg, 1:nyg, 1:nzg)
      close(u)

      open(newunit=u,&
      file=trim(fileplace)//"ErYAG/temp-timestep-deltdiv100-"&
         //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//"-"//str(pulsesToDo)//".dat" &
          ,access="stream",form="unformatted", status="replace")
      write(u)temp - 273.
      close(u)

      open(newunit=u,&
      file=trim(fileplace)//"ErYAG/water-timestep-deltdiv100-"&
         //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
      write(u)watercontent
      close(u)
   end if

call cpu_time(finish)
if(finish-start.ge.60.)then
    print*,str(floor((finish-start)/60.)+mod(finish-start,60.)/100.,5)//' mins'
else
    print*, 'time taken ~',floor(finish-start/60.),'s'
end if

call MPI_Finalize()
end program mcpolar
