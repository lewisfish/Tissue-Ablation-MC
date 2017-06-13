program ring

    use mpi

    implicit none

    ! variables for topology
    integer :: error, numproc, id, comm, dims(2), ndims, new_comm
    integer :: left, right, i, tag, j
    integer :: recv_status(MPI_STATUS_SIZE)
    logical :: periods(1), reorder

    !variables for heat
    integer :: size_x, size_y, n, p, o, u, numpoints, lo, hi
    real ::  k0, hx, hx2, hy, hy2, delt, u_xx, u_yy
    real, allocatable :: t0(:,:), tmp(:,:)
    character(len=1) :: f

    !init all mpi shizz
    comm = mpi_comm_world
    call MPI_init(error)
    call MPI_Comm_size(comm, numproc, error)

    tag = 1

    dims = 0
    periods = .false.
    reorder = .true.
    ndims = 2

    ! create cart topology based on number of processors and create new comm group
    call mpi_dims_create(numproc, ndims, dims, error)
    call mpi_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
    call MPI_Comm_rank(new_comm, id, error)

                        !comm, direction, shift, source, dest, error
    call mpi_cart_shift(new_comm, 0, 1, left, right, error)


    numpoints = 8

    if(id == 0)then
        !do decomp
        if(mod(numpoints,numproc) /= 0)then
            print('(I2.1,a,I2.1)'),numpoints,' not divisable by : ', numproc
            call mpi_finalize(error)
            error stop
        else
            N = numpoints / numproc
        end if
    end if

    call mpi_barrier(new_comm, error)
    call MPI_Bcast(N, 1, MPI_integer ,0 , new_comm, error)


!init heat shizz
    k0 = 100.
    size_x = N
    size_y = N
#ifdef serial
    hx = 1. / (size_x + 2)
    hy = 1. / (size_y + 2)

#else
    hx = 1. / (numpoints + 2)
    hy = 1. / (numpoints + 2)
#endif
    hx2 = 1. / hx**2
    hy2 = 1. / hy**2

    delt = 0.125 * min(hx,hy)**2/k0

    allocate(t0(0:size_x+1, 0:size_y+1))
    allocate(tmp(0:n*numproc+1, 0:n*numproc+1))
    t0 = 0.
    !boundary conditions
    if(id == 0)then
        tmp = 100.
        tmp(0,:) = 0.
        tmp(:,0) = 0.
        tmp(numpoints+1,:) = 0.
        tmp(:,numpoints+1) = 0.
        do i = 1, numproc-1
            lo = i*size_y
            hi = lo + size_y+1
            call mpi_send(tmp(:, lo:hi), size(t0), mpi_real, i, tag, new_comm, error)
        end do
        t0 = tmp(0:size_x+1, 0:size_y+1)
    else
        call mpi_recv(t0, size(t0), mpi_real, 0, tag, new_comm, recv_status, error)
    end if

    o = int(1./delt)

    !do heat sim
    do p = 1, o
        do i = 1, n
            do j = 1, n
                u_xx = (t0(i+1,j)     - 2.*t0(i,j) + t0(i-1, j))  / hx2
                u_yy = (t0(i,  j+1)   - 2.*t0(i,j) + t0(i,   j-1))  / hy2
                t0(i,j) = t0(i,j) + k0 * delt *(u_xx + u_yy)
            end do
        end do

        !exchange data
        if(mod(id, 2)== 0)then
            call MPI_Sendrecv(t0(:,size_y), size(t0(:,size_y)), mpi_real, right, tag, &
                              t0(:,size_y+1), size(t0(:,size_y+1)), mpi_real, right, tag, &
                              new_comm, recv_status, error)
        else
            call MPI_Sendrecv(t0(:,1), size(t0(:,1)), mpi_real, left, tag, &
                              t0(:,0), size(t0(:,0)), mpi_real, left, tag, &
                              new_comm, recv_status, error)
        end if

    end do

    if(id /= 0 )then
        call mpi_send(t0, size(t0), mpi_real, 0, tag, new_comm, error)
    else
        tmp = 0.
        do i = 1, numproc-1
            lo = i*size_x
            hi = lo + size_x+1
            call mpi_recv(tmp(:, lo:hi), size(t0), mpi_real, i, tag, new_comm, recv_status, error)
        end do
        tmp(:, 0:size_y+1) = t0

        write(f,'(I1.1)') numproc
        open(newunit=u, file='2d_parrallel_'//f//'.dat')
        do i = 0, size(tmp,1)-1
            write(u,*)tmp(i,:)
        end do
        close(u)
    end if

    call mpi_finalize(error)
end program ring