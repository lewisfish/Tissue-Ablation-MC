program tester

    use mpi

    implicit none
    
    logical :: periods(1), reorder
    integer :: error, numproc, id, comm, dims(1), ndims, new_comm, size_y, lo, hi
    integer :: left, right, i, j, u, xi, xf, yi, yf, coords(1), numpoints, recv_status(MPI_STATUS_SIZE)
    real, allocatable :: array(:,:), fin(:,:)

    comm = mpi_comm_world
    call mpi_init(error)
    call mpi_comm_size(comm, numproc, error)

    dims = 0
    periods = .false.
    reorder = .true.
    ndims = 1

    call mpi_dims_create(numproc, ndims, dims, error)
    call mpi_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
    call mpi_comm_rank(new_comm, id, error)

    call mpi_cart_coords(new_comm, id, ndims, coords, error)


    numpoints = 8


    size_y = numpoints / numproc

    xi = 1 
    xf = numpoints

    yi = 1          !id*size_y+1
    yf = size_y     !yi + size_y-1

    print*,xi,xf,yi,yf,id

    allocate(array(xi-1:xf+1, yi-1:yf+1))

    array = id

    if(id == 0)then
        allocate(fin(1:numpoints, 1:numpoints))
        do i = 1, numproc-1
            lo = i*size_y+1
            hi = lo + size_y-1
            print*,hi-lo,xf-xi,yf-yi!lo,hi,' ',xi,xf,yi,yf
            call mpi_recv(fin(1:numpoints,lo:hi), size(array(xi:xf,yi:yf)), mpi_real, i, 1, new_comm, recv_status, error)
        end do
        fin(1:numpoints,yi:yf)=array(xi:xf,yi:yf)
        open(newunit=u,file='test_decomp.dat')
        do i = 1,numpoints
            write(u,*)(fin(i,j),j=1,numpoints)
        end do
        close(u)
    else
        call mpi_send(array(1:xf,1:yf), size(array(xi:xf,yi:yf)), mpi_real, 0, 1, new_comm, error)
    end if




    call mpi_finalize(error)

end program tester