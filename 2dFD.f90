module threedfinite

    implicit none

    private
    public :: heat2D

    contains

    subroutine heat2D(jmean)

        ! use mpi
        use shrinkarray

        implicit none
        ! 
        real              :: u_xx, u_yy, diagx, diagy, weightx, weighty, tmp(50,50)
        real              :: delt, time, k0, hx, hy, hx2, hy2!, tmp, pwr, r
        real, allocatable :: T0(:,:), T(:,:), tissue(:,:), jmean(:,:,:), laser(:,:)
        integer           :: N, i, j, p, u, size_x, size_y, o

        ! integer :: error, numproc, id, comm, new_comm, ndims, dims(2), me, min_comm, col_id, col_coords(2)
        ! logical :: periods(2), reorder

        ! dims = 0
        ! periods = .false.
        ! reorder = .true.
        ! comm = MPI_Comm_world
        ! ndims = 2



    ! Create cartesian topology for processes
    !       dims(1) = nv
    !       dims(2) = mv
    !       call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, 
    !      &       period, reorder, grid_comm, ierr)
    !       call MPI_Comm_rank(grid_comm, me, ierr)
    !       call MPI_Cart_coords(grid_comm, me, ndim, coords, ierr)

    ! Create row subgrids ==>  3 subgrids of (1x2)
    !       remain(0) = .false.
    !       remain(1) = .true.
    !       call MPI_Cart_sub(grid_comm, remain, row_comm, ierr)

    !       call MPI_Comm_rank(row_comm, row_id, ierr)     
    !       call MPI_Cart_coords(row_comm, row_id, 1, row_coords, ierr)

    ! Create column subgrids ==>  2 subgrids of (3x1)
    !       remain(0) = .true.
    !       remain(1) = .false.
    !       call MPI_Cart_sub(grid_comm, remain, col_comm, ierr)

    !       call MPI_Comm_rank(col_comm, col_id, ierr)     
    !       call MPI_Cart_coords(col_comm, col_id, 1, col_coords, ierr)






        ! call MPI_init(error)

        ! call MPI_Comm_size(comm, numproc, error)

        ! call MPI_Comm_rank(comm, id, error)

        ! call MPI_dims_create(numproc, ndims, dims, error)
        ! call MPI_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
        ! ! call MPI_Comm_rank(new_comm, id, error)

        ! call MPI_Cart_sub(new_comm, [.true.,.false.], min_comm, error)
        ! call MPI_Comm_rank(min_comm, col_id, error)     
        ! call MPI_Cart_coords(min_comm, col_id, 1, col_coords, error)
        ! print*,col_id, col_coords

        ! call MPI_Finalize(error)
        ! call exit(0)

        ! if(mod(numproc, 2) /= 0 .and. numproc /= 1)then
        !     print*,'Need even number of process'
        !     ! call MPI_Finalize(error)
        !     stop
        ! else

        ! end if

        N = 50

        k0 = 1.

        size_x = N
        size_y = N

        hx = 1. / (size_x + 2)
        hy = 1. / (size_y + 2)
        print*,hx*50

        hx2 = 1. / hx**2
        hy2 = 1. / hy**2.

        delt = 0.125 * (min(hx,hy)**3)/k0
        
        !allocate mesh
        allocate(T(0:size_x+1, 0:size_y+1))
        allocate(jmean(200,200,200), laser(50,50))
        allocate(T0(0:size_x+1, 0:size_y+1))
        allocate(tissue(200, 200))

        inquire(iolength=i)jmean
        open(newunit=u,file='jmean.dat',access='direct',form='unformatted',recl=i)
        read(u,rec=1)jmean
        close(u)
        ! open(newunit=u,file='laser.dat')
        ! read(u,*)laser
        ! close(u)

        ! ! do i = 1, 200
            tissue(:,:) = jmean(100,:,:)
        ! ! end do
        call shrink(tissue, laser)
        deallocate(tissue)
        allocate(tissue(50,50))

        tissue = 0.

        t = 0.
        t0 = 37.
        t0(:,0) = 25. ! front face
        t0(:,N+1) = 37.  ! back face
        t0(N+1,:) = 37.  ! side face
        t0(0,:) = 37.    ! side face

        ! pwr = 2 !watts
        ! r = 0.1!in m

        ! do i = 0,51
        !     tmp = 2.*pwr/(4.*atan(1.)*r**2)*exp(-2.*((i-25)*hx)**2/r**2)
        !     if(tmp > 25)then
        !         t0(i,0) = tmp
        !     else
        !         t0(i,0) = 25.
        !     end if
        ! end do

        open(newunit=u,file='init.dat')
        do i = 0, n+1
            write(u,*)(t0(i,j),j=0,N+1)
        end do
       close(u)

        diagx = -2. + hx*hx/(3.*k0*delt)
        diagy = -2. + hy*hy/(3.*k0*delt)

        weightx = k0*delt/(hx*hx)
        weighty = k0*delt/(hy*hy)

        o = int(1./delt)

        do i = 1,50
            do j = 50, 1, -1
                tmp(i,j) = laser(i,j)
            end do
        end do
        laser = tmp

            open(newunit=u,file='laser.dat')
        do i = 1,50
            write(u,*)(laser(i,j),j=1,50)
        end do
        close(u)
        ! call exit(0)
        do p = 1, o
     !        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
     ! & " Percent Complete: ", (real(p)/real(o))*100.0, "%"
            if(mod(p,10000)==0)print*,real(p)/real(o)*100.,time
            do i = 1, size_x
                do j = 1,size_y
                    ! if(tissue(i,j) < 1.)then
                        u_xx = (t0(i+1,j)     - 2.*t0(i,j) + t0(i-1, j))  * hx2
                        u_yy = (t0(i,  j+1)   - 2.*t0(i,j) + t0(i,   j-1))  * hy2
                            t0(i,j) = t0(i,j) + k0*delt*(u_xx + u_yy)!+ k0*laser(i,j)
                    ! else
                    !     t0(i,j) = 37
                    ! end if
                end do
            end do
        time = time + delt
        ! call Arrhenius(t0, time, tissue)
        end do
        print*,
        print*,'time, loops, grid points'
        print*,time,p,n

        open(newunit=u,file='temperature.dat')
        do i = 1,n
            write(u,*)(t0(i,j),j=1,n)
        end do
        close(u)

        open(newunit=u,file='tissue.dat')
        do i = 1,n
            write(u,*)(tissue(i,j),j=1,n)
        end do
        close(u)
        ! call MPI_Finalize(error)

    end subroutine heat2D

    subroutine Arrhenius(Temp, time, tissue)

        implicit none

        real, intent(IN)    :: temp(:,:), time
        real, intent(INOUT) :: tissue(:,:)
        real :: a, g, r
        integer :: x, y

        a = 2.9e37
        g = 2.4e5
        r = 8.314

        do x = 1, 50
            do y = 1, 50
                if(temp(x, y) >= 42)then
                    tissue(x,y) =  time*A*exp(-G/(R*(temp(x,y)+273)))
                end if
            end do
        end do

    end subroutine  Arrhenius

end module threedfinite