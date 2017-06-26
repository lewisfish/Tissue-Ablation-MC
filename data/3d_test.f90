program ring

    use mpi

    use utils, only : str

    implicit none

    ! variables for topology
    integer :: error, numproc, id, comm, dims(2), ndims, new_comm
    integer :: left, right, i, tag, j, k, xi, xf, yi, yf, zi, zf
    integer :: recv_status(MPI_STATUS_SIZE)
    logical :: periods(1), reorder

    !variables for heat
    integer :: size_x, size_y, size_z, n, p, o, u, numpoints, lo, hi
    real ::  k0, hx, hx2, hy, hy2, hz, hz2, delt, u_xx, u_yy, u_zz,start, finish
    real, allocatable :: tmp(:,:,:), out(:,:,:)

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

    numpoints = 48

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

    if(id==0)call cpu_time(start)

    !init heat shizz
    k0 = 1.
    size_x = numpoints
    size_y = numpoints
    size_z = N

    hx = 1. / (numpoints + 2)
    hy = 1. / (numpoints + 2)
    hz = 1. / (numpoints + 2)

    hx2 = 1. / hx**2
    hy2 = 1. / hy**2
    hz2 = 1. / hz**2

    delt = 0.125 * min(hx,hy,hz)**3/k0

    xi = 1
    xf = size_x

    yi = 1
    yf = size_y

    zi = 1
    zf = size_z

    !allocate strips to each processor in z dir for speed/memory
    allocate(tmp(0:numpoints+1, 0:numpoints+1, zi-1:zf+1))

    !boundary conditions
    tmp = 100.
    tmp(xi-1,:,:) = 0. !side face
    tmp(xf+1,:,:) = 0. !side face
    tmp(:,yi-1,:) = 0. !front face
    tmp(:,yf+1,:) = 0. !back face
    if(id == numproc-1)tmp(:,:,zf+1) = 0. !top face
    if(id == 0)tmp(:,:,zi-1) = 0. !bottom face

    !wrtie out init grid
    if(id == 0)then
        allocate(out(0:numpoints+1,0:numpoints+1, 0:numpoints+1))
        out = 0.

        do i = 1, numproc-1
            lo = i*size_z+1
            hi = lo + size_z-1
            call mpi_recv(out(:,:, lo:hi), size(out(:,:,lo:hi)), mpi_real, i, tag, new_comm, recv_status, error)
        end do

        out(:,:,zi:zf) = tmp(:,:,zi:zf)

        open(newunit=u,file='init_3d.dat',access='stream',form='unformatted',status='replace')
        write(u)out
        close(u)
    else
        call mpi_send(tmp(:, :, :), size(tmp(:, :,zi:zf)), mpi_real, 0, tag, new_comm, error)
    end if

    !work out # of steps for sim of tim x => $\frac{x}{\del t}$
    o = int(.1/delt)
    !do heat sim
    do p = 1, o
        if(mod(p, 1000) == 0 .and. id == 0) write(*,FMT="(A8,t21)",ADVANCE="NO") achar(13)//str(real(p)/real(o)*100., 5)//' %'
        do i = xi, xf
            do j = yi, yf
                do k = zi, zf
                u_xx = (tmp(i+1,j,   k)   - 2.*tmp(i,j,k) + tmp(i-1, j,   k))  * hx2
                u_yy = (tmp(i,  j+1, k)      - 2.*tmp(i,j,k) + tmp(i,   j-1, k))  * hy2
                u_zz = (tmp(i,   j,  k+1) - 2.*tmp(i,j,k) + tmp(i,   j,   k-1))* hz2

                tmp(i,j,k) = tmp(i,j,k) + k0 * delt *(u_xx + u_yy + u_zz)
                end do
            end do
        end do
        call mpi_barrier(new_comm, error)

        !send_recv data to right
        call MPI_Sendrecv(tmp(:,:,zf), size(tmp(:,:,zf)), mpi_real, right, tag, &
                          tmp(:,:,zf+1), size(tmp(:,:,zf+1)), mpi_real, right, tag, &
                          new_comm, recv_status, error)
        
        !send_recv data to left
        call MPI_Sendrecv(tmp(:,:,zi), size(tmp(:,:,zi)), mpi_real, left, tag, &
                          tmp(:,:,zi-1), size(tmp(:,:,zi-1)), mpi_real, left, tag, &
                          new_comm, recv_status, error)

        call mpi_barrier(new_comm, error)
    end do



    if(id==0)then
        print*,
        call cpu_time(finish)
        print*,finish-start
    end if

    !send data to master process and do I/O
    if(id == 0)then
        out = 0.
        do i = 1, numproc-1
            lo = i*size_z+1
            hi = lo + size_z-1
            print*,lo,hi
            call mpi_recv(out(:,:, lo:hi), size(out(:,:,lo:hi)), mpi_real, i, tag, new_comm, recv_status, error)
        end do

        out(:,:,zi:zf) = tmp(:,:,zi:zf)

        open(newunit=u, file='3d_parrallel_'//str(numproc)//'.dat',access='stream',form='unformatted',status='replace')
        write(u)out
        close(u)
        print*,'3d_parrallel_'//str(numproc)//'.dat'
    else
        call mpi_send(tmp(:, :, zi:zf), size(tmp(:, :,zi:zf)), mpi_real, 0, tag, new_comm, error)
    end if

    call mpi_finalize(error)
end program ring