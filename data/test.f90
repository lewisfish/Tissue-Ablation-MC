program tester

    use mpi

    implicit none
    
    logical :: periods(1), reorder
    integer :: error, numproc, id, comm, dims(1), ndims, new_comm, size_y, lo, hi, zi, zf
    integer :: i, u, xi, xf, yi, yf, coords(1), numpoints, recv_status(MPI_STATUS_SIZE)
    real, allocatable :: array(:,:,:), fin(:,:,:)

    integer :: starts(3), oldsize(3), newsize(3), arr

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

    numpoints = 100

    if(mod(numpoints, numproc) /= 0)then
        call mpi_finalize(error)
        error stop 1
    end if

    size_y = numpoints / numproc

    xi = 1 
    xf = numpoints

    yi = 1
    yf = numpoints

    zi = 1
    zf = size_y

    print'(7(I3.1,1X))',xi,xf,yi,yf,zi,zf,id

    allocate(array(xi-1:xf+1, yi-1:yf+1,zi-1:zf+1))


    oldsize = [size(array,1), size(array,2), size(array,3)]
    newsize = [[size(array,1)-2, size(array,2)-2, size(array,3)-2]]
    starts = [xi, yi,zi]
                                  !ndims, array_of_sizes[], array_of_subsizes,array_of_starts,order,oldtype,newtype)
    call mpi_type_create_subarray(3, oldsize, newsize, starts, mpi_order_fortran, &
                                  mpi_real, arr, error)
    call mpi_type_commit(arr,error)


    array = (id+1)**2

    
    if(id == 0)then
        allocate(fin(1:numpoints, 1:numpoints, 1:numpoints))
        fin = 100.
        do i = 1, numproc-1
            lo = i*size_y+1
            hi = lo + size_y-1
            print*,lo,hi
            call mpi_recv(fin(:,:,lo:hi), size(array(xi:xf,yi:yf,zi:zf)), &
                          mpi_real, i, 1, new_comm, recv_status, error)
        end do
        fin(:,:,zi:zf)=array(xi:xf,yi:yf,zi:zf)
        open(newunit=u,file='test_decomp.dat',access='stream',form='unformatted',status='replace')
        write(u)fin
        close(u)

        open(newunit=u,file='test_slice.dat',status='replace')
        do i = 1, numpoints
            write(u,*)fin(numpoints/2,numpoints/2,i)
        end do
        close(u)
    else
        call mpi_send(array, 1, arr, 0, 1, new_comm, error)
    end if




    call mpi_finalize(error)

end program tester