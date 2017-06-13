program ring

    use mpi

    implicit none

    ! variables for topology
    integer :: error, numproc, id, comm, dims(1), ndims, new_comm
    integer :: left, right, i, tag
    integer :: recv_status(MPI_STATUS_SIZE)
    logical :: periods(1), reorder

    !variables for heat
    integer :: size_x, n, p, o, u, numpoints, lo, hi
    real ::  k0, hx, hx2, delt
    real, allocatable :: t0(:), tmp(:), tmp2(:)
    character(len=1) :: f

!init all mpi shizz
    comm = mpi_comm_world
    call MPI_init(error)
    call MPI_Comm_size(comm, numproc, error)

    tag = 1

    dims = 0
    periods = .false.
    reorder = .true.
    ndims = 1

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
#ifdef serial
    hx = 1. / (size_x + 2)
#else
    hx = 1. / (numpoints + 2)
#endif
    hx2 = 1. / hx**2
    delt = 0.125 * hx/k0

    allocate(t0(0:size_x+1))
    allocate(tmp(0:n*numproc+1), tmp2(0:n*numproc+1))

    !boundary conditions
    if(id == 0)then
        tmp = 0.
        tmp(0) = 100.
        do i = 0, numproc-1
            lo = i*size_x
            hi = lo + size_x+1
            call mpi_send(tmp(lo:hi), size(t0), mpi_real, i, tag, new_comm, error)
        end do
        t0 = tmp(0:size_x+1)
    else
        lo = i*size_x
        hi = lo + size_x+1
        call mpi_recv(t0, size(t0), mpi_real, 0, tag, new_comm, recv_status, error)
    end if

    o = int(1./delt)

    !do heat sim
    do p = 1, o
        do i = 1, n
            t0(i) = t0(i) + k0 * delt * (t0(i+1) - 2.*t0(i) + t0(i-1))/ hx2
        end do

        !exchange data
        if(mod(id, 2)== 0)then
            call MPI_Sendrecv(t0(size_x), 1, mpi_real, right, tag, &
                              t0(size_x+1), 1, mpi_real, right, tag, &
                              new_comm, recv_status, error)
        else
            call MPI_Sendrecv(t0(1), 1, mpi_real, left, tag, &
                              t0(0), 1, mpi_real, left, tag, &
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
            call mpi_recv(tmp(lo:hi), size(t0), mpi_real, i, tag, new_comm, recv_status, error)
        end do
        tmp(0:size_x+1) = t0

        write(f,'(I1.1)') numproc
        open(newunit=u, file='1d_parrallel_'//f//'.dat')
        do i = 0, size(tmp)-1
            write(u,*)tmp(i)
        end do
        close(u)
    end if

    call mpi_finalize(error)
end program ring