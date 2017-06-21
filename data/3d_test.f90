program ring

    use mpi

    use utils, only : str

    implicit none

    ! variables for topology
    integer :: error, numproc, id, comm, dims(2), ndims, new_comm
    integer :: left, right, i, tag, j, k, xi, xf, yi, yf, zi, zf
    integer :: recv_status(MPI_STATUS_SIZE)
    logical :: periods(1), reorder

    !vars for subarray
    integer :: starts(3), oldsize(3), newsize(3), subarr

    !variables for heat
    integer :: size_x, size_y, size_z, n, p, o, u, numpoints, lo, hi
    real ::  k0, hx, hx2, hy, hy2, hz, hz2, delt, u_xx, u_yy, u_zz
    real, allocatable :: tmp(:,:,:)

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

    numpoints = 40

    if(id == 0)then
        !do decomp
        if(mod(numpoints, numproc) /= 0)then
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
    k0 = 1.
    size_x = numpoints
    size_y = N
    size_z = numpoints

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

    ! yi = id*size_y+1
    ! if(yi == 0)yi = 1
    ! yf = yi + size_y-1

    zi = 1
    zf = size_z

    ! allocate(tmp(0:n*numproc+1, 0:n*numproc+1, 0:n*numproc+1)) !old
    allocate(tmp(xi-1:xf+1,yi-1:yf+1,zi-1:zf+1))

    !create mpi sub array to avoid temp array creation
    oldsize = [size(tmp,1), size(tmp,2), size(tmp,3)]
    newsize = [[size(tmp,1)-2, 1, size(tmp,3)-2]]
    starts = [xi, yf, zi]

    print*,oldsize,newsize,starts

    call mpi_type_create_subarray(3, oldsize, newsize, starts, mpi_order_fortran, &
                                  mpi_real, subarr, error)
    call mpi_type_commit(subarr, error)


    !boundary conditions
    tmp = 100.
    tmp(0,:,:) = 0. !side face
    tmp(:,:,0) = 0. !bottom face
    tmp(numpoints+1,:,:) = 0. !side face
    tmp(:,:,numpoints+1) = 0. !top face

    if(id==numproc-1)tmp(:,yf+1,:) = 0. !back face
    if(id==0)tmp(:,0,:) = 0. !front face


    o = int(.1/delt)
    !do heat sim
    do p = 1, o
        if(mod(p, 1000) == 0 .and. id == 0)print*,str(real(p)/real(o)*100.),' %'
        do i = xi, xf
            do j = yi, yf!1, n
                do k = zi, zf
                u_xx = (tmp(i+1,j,   k)   - 2.*tmp(i,j,k) + tmp(i-1, j,   k))  * hx2
                u_yy = (tmp(i,  j+1, k)      - 2.*tmp(i,j,k) + tmp(i,   j-1, k))  * hy2
                u_zz = (tmp(i,   j,  k+1) - 2.*tmp(i,j,k) + tmp(i,   j,   k-1))* hz2

                tmp(i,j,k) = tmp(i,j,k) + k0 * delt *(u_xx + u_yy + u_zz)
                end do
            end do
        end do
        call mpi_barrier(new_comm, error)

        !need to add subarray send here. might have wrong subarray definition though...

        !send_recv data to right
        call MPI_Sendrecv(tmp(:,yf,:), size(tmp(:,yf,:)), mpi_real, right, tag, &
                          tmp(:,yf+1,:), size(tmp(:,yf+1,:)), mpi_real, right, tag, &
                          new_comm, recv_status, error)
        
        !send_recv data to right
        call MPI_Sendrecv(tmp(:,yi,:), size(tmp(:,yi,:)), mpi_real, left, tag, &
                          tmp(:,yi-1,:), size(tmp(:,yi-1,:)), mpi_real, left, tag, &
                          new_comm, recv_status, error)

        call mpi_barrier(new_comm, error)
    end do

    !send data to master process
    if(id == 0)then
        do i = 1, numproc-1
            lo = i*size_y+1
            hi = lo + size_y-1
            call mpi_recv(tmp(:, lo:hi,:), size(tmp(:,lo:hi,:)), mpi_real, i, tag, new_comm, recv_status, error)
        end do
    else
        call mpi_send(tmp(:, yi:yf,:), size(tmp(:, yi:yf,:)), mpi_real, 0, tag, new_comm, error)
    end if

    !write data out on master process
    if(id == 0)then
        open(newunit=u, file='3d_parrallel_'//str(numproc)//'.dat',access='stream',form='unformatted')
        ! do i = 0, size(tmp,1)
        !     write(u,*)tmp(i,:)
        ! end do
        write(u)tmp
        close(u)
        print*,'3d_parrallel_'//str(numproc)//'.dat'
    end if
    call mpi_finalize(error)
end program ring