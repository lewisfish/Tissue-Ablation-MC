program threedfinite

    use mpi

    implicit none
    ! 
    real              :: u_xx, u_yy, diagx, diagy, weightx, weighty
    real              :: delt, time, k0, hx, hy, hx2, hy2
    real, allocatable :: T0(:,:), T(:,:)
    integer           :: N, i, j, p, u, size_x, size_y, o

    integer :: error, numproc, id, comm, new_comm, ndims, dims(2), me, min_comm, col_id, col_coords(2)
    logical :: periods(2), reorder

    dims = 0
    periods = .false.
    reorder = .true.
    comm = MPI_Comm_world
    ndims = 2



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






    call MPI_init(error)

    call MPI_Comm_size(comm, numproc, error)

    call MPI_Comm_rank(comm, id, error)

    call MPI_dims_create(numproc, ndims, dims, error)
    call MPI_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
    ! call MPI_Comm_rank(new_comm, id, error)

    call MPI_Cart_sub(new_comm, [.true.,.false.], min_comm, error)
    call MPI_Comm_rank(min_comm, col_id, error)     
    call MPI_Cart_coords(min_comm, col_id, 1, col_coords, error)
    print*,col_id, col_coords

    call MPI_Finalize(error)
    call exit(0)

    if(mod(numproc, 2) /= 0 .and. numproc /= 1)then
        print*,'Need even number of process'
        call MPI_Finalize(error)
        stop
    else

    end if
    N = 50

    k0 = 1.

    size_x = N
    size_y = N

    hx = 1. / (size_x + 2)
    hy = 1. / (size_y + 2)

    hx2 = 1. / hx**2
    hy2 = 1. / hy**2.

    delt = 0.125 * (min(hx,hy)**3)/k0
    
    !allocate mesh
    allocate(T(0:size_x+1, 0:size_y+1))
    allocate(T0(0:size_x+1, 0:size_y+1))

    t = 0.
    t0 = 37.
    t0(:,0) = 0. ! front face
    t0(:,N+1) = 0.  ! back face
    t0(N+1,:) = 0.  ! side face
    t0(0,:) = 0.    ! side face

   inquire(iolength=i)t0(1:N,1:N)

   open(newunit=u,file='init.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) t0(1:N,1:N)
   close(u)


    diagx = -2. + hx*hx/(3.*k0*delt)
    diagy = -2. + hy*hy/(3.*k0*delt)

    weightx = k0*delt/(hx*hx)
    weighty = k0*delt/(hy*hy)

    o = int(0.1/delt)

    do p = 1, o
 !        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
 ! & " Percent Complete: ", (real(p)/real(o))*100.0, "%"
        do i = 1, size_x
            do j = 1,size_y
                u_xx = (t0(i+1,j)     - 2.*t0(i,j) + t0(i-1, j))  * hx2
                u_yy = (t0(i,  j+1)   - 2.*t0(i,j) + t0(i,   j-1))  * hy2
                t0(i,j) = t0(i,j) + k0*delt*(u_xx + u_yy)
            end do
        end do
    time = time + delt
    end do
    print*,
    print*,time,p,n

   open(newunit=u,file='temperature.dat')!,access='direct',status='REPLACE',form='unformatted',&
    do i =1,n
        write(u,*)(t0(i,j),j=1,n)
    end do
    close(u)
    call MPI_Finalize(error)

end program threedfinite