
!*****************************************************************************************
! * HEAT2D Example - Parallelized FORTRAN Version
! * FILE: 2d_heat_equation.f
! * DESCRIPTIONS: This example is based on a simplified two-dimensional heat
!       * equation domain decomposition. The initial temperature is computed to be
!       * high in the middle of the domain and zero at the boundaries. The
!       * boundaries are held at zero throughout the simulation. During the
!       * time-stepping, an array containing two domains is used; these domains
!       * alternate between old data and new data.
!       *
!       * In this parallelized version, the grid is decomposed by the master
!       * process and then distributed by rows to the worker processes. At each
!       * time step, worker processes must exchange border data with neighbors,
!       * because a grid point's current temperature depends upon it's previous
!       * time step value plus the values of the neighboring grid points. Upon
!       * completion of all time steps, the worker processes return results
!       * to the master process.
!       *
!       * AUTHOR: Blaise Barney - adapted from D. Turner's serial version
!       * CONVERTED TO MPI: George L. Gusciora (1/25/95)
!       * MODIFIED BY: C. B. Connor (6/6/02)
!               * CONVERTED TO FORTRAN 77 BY: R. M. Mason (6/11/02)
!       ****************************************************************************

module common

    implicit none

    ! include 'param.in'
    integer, PARAMETER :: NXPROB=50
    integer, PARAMETER :: NYPROB=50 
    double precision :: U(0:NXPROB-1,0:NYPROB-1,0:1)
    integer :: iz

end module common



      program heat
      
      use mpi
      use common

      implicit none
      
      ! include 'param.in'

      integer :: status(MPI_STATUS_SIZE)

      integer :: TIME_STEPS,MAXWORKER,MINWORKER,BEGIN
      integer :: NGHBOR1,NGHBOR2,NONE,DONE,MASTER
      double precision :: CX,CY

      ! double precision U(0:NXPROB-1,0:NYPROB-1,0:1)
!         array for grid
      integer taskid
!         this task's unique id
      integer numworkers 
!         number of worker processes
      integer numtasks
!         number of tasks
      integer min_number_rows,number_rows,offset,extra_rows
!         for sending rows of data
      integer destination,source
!         to - from for message send-receive
      integer worker_number,neighbor1,neighbor2
!         neighbor tasks
      integer message_tag
!         for message types
      integer ierr,istart,iend
!         error, start and end for loops on each node
      integer i,ix,iy,it
!         loop variables


      TIME_STEPS=500
!          number of time steps
        MAXWORKER=8
!          maximum number of worker tasks 
       MINWORKER=3
!          minimum number of worker tasks
       BEGIN=1
!          message type
       NGHBOR1=2
!          message type
       NGHBOR2=3
!          message type
       NONE=0
!         indicates no neighbour
       DONE=4
!         message type
       MASTER=0
!         taskid of first process
       CX=0.1
       CY=0.1
!         diffusivity coefficients in x and y


!     First, find out my taskid and how many tasks are running 
      
      call mpi_init(ierr)

      call mpi_comm_size(MPI_COMM_WORLD, numtasks, ierr)

      call mpi_comm_rank(MPI_COMM_WORLD, taskid, ierr)

!     one of the processors is completely reserved for coordination

      numworkers = numtasks-1

!     the MASTER node subdivides the problem by grid decomposition,
!     in this case by columns. The MASTER also initializes the starting 
!     values, sets the boundary conditions on the problem , 
!     and receives results from the nodes after TIME_STEPS

!     NOTE that columns are used rather than rows, since FORTRAN stores arrays
!     by column.  THe opposite is true in C.

! ****************** master code ***************************

      iz=0
         if((numworkers.gt.maxworker).or.(numworkers.lt.minworker)) then
            if(taskid==master)write (6,'(a,i2,a,i2,a)') 'number of processors needs to be between ',minworker+1,&
            ' and ',maxworker+1,' for this exercise'
            call mpi_finalize(ierr)
            stop
         end if
      if(taskid.eq.MASTER) then

!     Check if numworkers is within range - quit if not


      
!     initialise data

         write(6,'(a,i3,a,i3,a,i4)') 'grid size X = ',NXPROB,' Y = ',NYPROB,' Time steps = ',TIME_STEPS

         print*, 'Initialising grid and writing initial.dat file'

         call inidat(NXPROB,NYPROB)
         call prtdat(NXPROB,NYPROB,0)
      
!     Distribute work to workers. Must first figure out how many rows to
!     send and what to do with extra rows.
         min_number_rows = NYPROB/numworkers
         extra_rows = NYPROB-min_number_rows*numworkers
         offset = 0

         do i=1,numworkers
         
!     The following is a particularly compact and efficient  way of distributing the 
!     grid. It assures that the number of rows received by each node is only up to one more
!     than any other node. 

            if(i.le.extra_rows) then
               number_rows=min_number_rows+1
            else
               number_rows=min_number_rows
            end if

!     Tell each worker who its neighbors are, since they must exchange
!     data with each other.

            if(i.eq.1) then
               neighbor1 = NONE
            else
               neighbor1 = i - 1
            end if

            if(i.eq.numworkers) then
               neighbor2 = NONE
            else
               neighbor2 = i + 1
            end if

!     Now send startup information to each worker Note that this 
!     information is "tagged" as the BEGIN message

            destination = i
            worker_number = i
            message_tag = BEGIN
            
!     Send the required information to each node
            
            call MPI_Send(worker_number, 1, MPI_INTEGER, destination, &
                message_tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(offset, 1, MPI_INTEGER, destination,message_tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(number_rows, 1, MPI_INTEGER, destination, message_tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(neighbor1, 1, MPI_INTEGER, destination,message_tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(neighbor2, 1, MPI_INTEGER, destination, message_tag, MPI_COMM_WORLD,ierr)
            call MPI_Send(U(0,offset,iz),NXPROB*number_rows, MPI_DOUBLE_PRECISION, destination, message_tag, MPI_COMM_WORLD, ierr)

!     let the world know how the problem has been divided up

            write(6,300) destination, offset, number_rows, neighbor1,neighbor2            
300  format('Sent to = ', I2, ' offset = ', I3, ' number_rows = ', I3,' neighbor 1 = ', I2, ' neighbor 2 = ', I2)
            
!     increment the offset by the number_rows so the next node will
!     know where its grid begins
            
            offset = offset+number_rows
            
!     continue doing the above for each node, i
            
         end do

!     Now wait for results from all worker tasks
         
         do i=1,numworkers
            source = i
            message_tag = DONE
            call MPI_Recv(offset, 1, MPI_INTEGER, source, message_tag, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(number_rows, 1, MPI_INTEGER, source, message_tag, MPI_COMM_WORLD, status, ierr)
            call MPI_Recv(U(0,offset,iz),NXPROB*number_rows, MPI_DOUBLE_PRECISION, source, &
                message_tag, MPI_COMM_WORLD, status, ierr)
         end do
         
         
         call prtdat(NXPROB, NYPROB,1)
         
      end if
      
!     End of master code
      
      if(taskid.ne.MASTER) then
         
!     ************************* worker code*********************************
         
!     Initialize everything - including the borders - to zero
         
         do ix=0,nxprob-1
            do iy=0,nyprob-1
               u(ix,iy,0)=0.d0
            end do
         end do

!     Now receive my offset, rows, neighbors and grid partition from master

         source = MASTER
         message_tag = BEGIN
         call MPI_Recv(worker_number, 1, MPI_INTEGER, source,message_tag, MPI_COMM_WORLD, status, ierr)
         call MPI_Recv(offset, 1, MPI_INTEGER, source,message_tag, MPI_COMM_WORLD, status, ierr)
         call MPI_Recv(number_rows, 1, MPI_INTEGER, source,message_tag, MPI_COMM_WORLD, status, ierr)
         call MPI_Recv(neighbor1, 1, MPI_INTEGER, source,message_tag, MPI_COMM_WORLD, status, ierr)
         call MPI_Recv(neighbor2, 1, MPI_INTEGER, source,message_tag, MPI_COMM_WORLD, status, ierr)
         call MPI_Recv(U(0,offset,iz),NXPROB*number_rows, MPI_DOUBLE_PRECISION, source, message_tag, MPI_COMM_WORLD, status, ierr)
         
!     Determine border elements. This takes into account that the
!     first and last rows have fixed temperature
         
         if(offset.eq.0) then
            istart=1
!     do not include row zero
         else
            istart=offset
         end if
         if ((offset+number_rows).eq.NYPROB) then
            iend=offset+number_rows-2  
!     do not include the last row
         else
            iend = offset + number_rows-1
         end if
         
         
!     take a look at how the work is partioned among processors
         
         write(6,400) worker_number, offset,number_rows,istart,iend
 400     format('Worker number = ', I2, ' offset = ', I3, ' number_rows = ', I3, ' start = ', I3, ' end = ', I3)
         
!     Begin doing TIME_STEPS iterations. Must communicate border rows with
!     neighbors. If I have the first or last grid row, then I only need to
!     communicate with one neighbor
         
         do it = 1, time_steps
            if(neighbor1 /= NONE) then 
           call MPI_Send(U(0,offset,iz), NXPROB,MPI_DOUBLE_PRECISION, neighbor1, NGHBOR2, MPI_COMM_WORLD, ierr)
               source = neighbor1
               message_tag = NGHBOR1
           call MPI_Recv(U(0,offset-1,iz), NXPROB, MPI_DOUBLE_PRECISION, source, message_tag,MPI_COMM_WORLD, status, ierr)
        end if
            if (neighbor2 /= NONE)then
           call MPI_Send(U(0,offset+number_rows-1,iz), NXPROB, MPI_DOUBLE_PRECISION, neighbor2, NGHBOR1, MPI_COMM_WORLD, ierr)
               source = neighbor2
               message_tag = NGHBOR2
           call MPI_Recv(U(0,offset+number_rows,iz),  NXPROB, MPI_DOUBLE_PRECISION, source, message_tag, &
                      MPI_COMM_WORLD, status, ierr)
            end if
            
!     Now call update to update the value of grid points 
            
            call update(istart,iend,NXPROB,cx,cy)
            
         end do
         
!     Finally, send my portion of final results back to master
         
         call MPI_Send(offset, 1, MPI_INTEGER, MASTER, DONE, MPI_COMM_WORLD, ierr)
         call MPI_Send(number_rows, 1, MPI_INTEGER, MASTER, DONE, MPI_COMM_WORLD, ierr)
         call MPI_Send(U(0,offset,0), NXPROB*number_rows, MPI_DOUBLE_PRECISION, MASTER, DONE, MPI_COMM_WORLD, ierr)
      end if
      
!     gracefully exit MPI
      
      call MPI_Finalize(ierr)
      
      contains
!**************************************************************************
! * subroutine update
!**************************************************************************

      subroutine update(istart,iend,nx,cx,cy)
      
      use common

      implicit none
      
      ! include 'param.in'
      
      integer istart,iend,ix,iy,nx
      double precision cx,cy
      ! double precision U(0:NXPROB-1,0:NYPROB-1,0:1)

      do iy=istart,iend
         do ix=1,nx-2
            U(ix,iy,1-iz)=U(ix,iy,iz)+cx*(U(ix+1,iy,iz)+U(ix-1,iy,iz) &
                -2.d0*U(ix,iy,iz))+cy*(U(ix,iy+1,iz)+U(ix,iy-1,iz)- &
                2.d0*U(ix,iy,iz))
         end do
      end do
 
      iz=1-iz

      end subroutine update

!*****************************************************************************
!     * subroutine inidat
!*****************************************************************************/
      
      subroutine inidat(nx,ny)
      
        use common

      implicit none
      
      ! include 'param.in'

      integer nx,ny,ix,iy
      ! double precision U(0:NXPROB-1,0:NYPROB-1,0:1)
            
      do ix=0,nx-1
         do iy=0,ny-1
            u(ix,iy,0)=dble(ix*(nx-ix-1)*iy*(ny-iy-1))/dble(4*nx*ny)
         end do
      end do
      end subroutine inidat
      
!**************************************************************************
!* subroutine prtdat
!**************************************************************************

      subroutine prtdat(nx,ny,initorfin)

        use common

      implicit none

      ! include 'param.in'
      
      integer nx,ny,ix,iy,initorfin
      ! double precision U(0:NXPROB-1,0:NYPROB-1,0:1)
            
      if(initorfin.eq.0) then
         open(10, file='initial.dat')
      else
         open(10, file='final.dat')
      end if

      do iy=ny-1,0,-1
         do ix=0,nx-1
            write(10,'(I3,1X,I3,1X,E12.5)') iy,ix,u(ix,iy,0)
         end do
            write(10,*)' '
      end do
      close(10)      
      end subroutine prtdat
    end program
