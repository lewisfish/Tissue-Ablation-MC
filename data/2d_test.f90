program ring

    use mpi

    use utils, only : str

    implicit none

    ! variables for topology
    integer :: error, numproc, id, comm, dims(2), ndims, new_comm
    integer :: left, right, i, tag, j, xi, xf, yi, yf
    integer :: recv_status(MPI_STATUS_SIZE)
    logical :: periods(1), reorder

    !variables for heat
    integer :: size_x, size_y, n, p, o, u, numpoints, lo, hi
    real ::  k0, hx, hx2, hy, hy2, delt, u_xx, u_yy
    real, allocatable :: t0(:,:), tmp(:,:), out(:,:)

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


!     print*,id,right,left
! call exit(0)
    numpoints = 40

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
    k0 = 1.
    size_x = numpoints
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

    delt = 0.125 * min(hx,hy)**3/k0
    allocate(t0(0:n*numproc+1, 0:size_y+1))
    allocate(tmp(0:n*numproc+1, 0:n*numproc+1))
    allocate(out(0:n*numproc+1, 0:n*numproc+1))

    xi = 1
    xf = size_x

    yi = id*size_y+1
    if(yi == 0)yi = 1
    yf = yi + size_y-1

    print*,xi,xf,yi,yf,id
    out = 0.
    t0 = 0.
    !boundary conditions
    ! if(id == 0)then
        tmp = 0.
        tmp = 100.
        tmp(0,:) = 0.
        tmp(:,0) = 0.
        tmp(numpoints+1,:) = 0.
        tmp(:,numpoints+1) = 0.
    !     do i = 1, numproc-1
    !         lo = i*size_y
    !         hi = lo + size_y+1
    !         ! call mpi_send(tmp(:, lo:hi), size(t0), mpi_real, i, tag, new_comm, error)
    !         call mpi_send(tmp, size(tmp), mpi_real, i, tag, new_comm, error)

    !     end do
    !     ! t0 = tmp(:, 0:size_y+1)
    ! else
    !     call mpi_recv(tmp, size(tmp), mpi_real, 0, tag, new_comm, recv_status, error)
    !     ! call mpi_recv(t0, size(t0), mpi_real, 0, tag, new_comm, recv_status, error)
    ! end if

    o = int(.1/delt)
    !do heat sim
    do p = 1, o
        if(mod(p, 10000) == 0 .and. id == 0)print*,str(real(p)/real(o)*100.),' %'
        do i = 1, numpoints
            do j = yi, yf!1, n
                u_xx = (tmp(i+1,j)     - 2.*tmp(i,j) + tmp(i-1, j))  * hx2
                u_yy = (tmp(i,  j+1)   - 2.*tmp(i,j) + tmp(i,   j-1))  * hy2
                tmp(i,j) = tmp(i,j) + k0 * delt *(u_xx + u_yy)
            end do
        end do
        call mpi_barrier(new_comm, error)

        !exchange data
        if(mod(id, 2)== 0)then
            ! print*,yf
            call MPI_Sendrecv(tmp(:,20), size(tmp(:,20)), mpi_real, right, tag, &
                              tmp(:,21), size(tmp(:,21)), mpi_real, right, tag, &
                              new_comm, recv_status, error)
        else
            call MPI_Sendrecv(tmp(:,21), size(tmp(:,21)), mpi_real, left, tag, &
                              tmp(:,20), size(tmp(:,20)), mpi_real, left, tag, &
                              new_comm, recv_status, error)
        end if

        call mpi_barrier(new_comm, error)

        ! if(mod(id, 2) /= 0)then
        !     call MPI_Sendrecv(tmp(:,yf), size(tmp(:,yf)), mpi_real, right, tag, &
        !                       tmp(:,yf+1), size(tmp(:,yf+1)), mpi_real, right, tag, &
        !                       new_comm, recv_status, error)
        ! else
        !     call MPI_Sendrecv(tmp(:,yi+1), size(tmp(:,yi+1)), mpi_real, left, tag, &
        !                       tmp(:,yi), size(tmp(:,yi)), mpi_real, left, tag, &
        !                       new_comm, recv_status, error)
        ! end if

    end do

    call mpi_allreduce(tmp, out, size(tmp), mpi_real, mpi_sum, new_comm, error)
    if(id == 0)then
        out=out/real(numproc)
        open(newunit=u, file='2d_parrallel_'//str(numproc)//'.dat')
        do i = 0, size(out,1)
            write(u,*)out(i,:)
        end do
        close(u)


        open(newunit=u,file='2d_par_slice_'//str(numproc)//'.dat')
        do i = 1, size(out,1)-1
            write(u,*)out((size(out,1)-1)/2,i)
        end do
        close(u)
    end if




    ! if(id /= 0 )then
    !     call mpi_send(t0, size(t0), mpi_real, 0, tag, new_comm, error)
    ! else
    !     tmp = 0.
    !     do i = 1, numproc-1
    !         lo = i*size_y
    !         hi = lo + size_y+1
    !         call mpi_recv(tmp(:, lo:hi), size(t0), mpi_real, i, tag, new_comm, recv_status, error)
    !     end do
    !     tmp(:, 0:size_y+1) = t0

    !     open(newunit=u, file='2d_parrallel_'//str(numproc)//'.dat')
    !     do i = 0, size(tmp,1)-1
    !         write(u,*)tmp(i,:)
    !     end do
    !     close(u)


    !     open(newunit=u,file='2d_par_slice_'//str(numproc)//'.dat')
    !     do i = 0, size(tmp,1)-1
    !         write(u,*)tmp((size(tmp,1)-1)/2,i)
    !     end do
    !     close(u)
    ! end if

    call mpi_finalize(error)
end program ring