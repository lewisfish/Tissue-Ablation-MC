program tester

    use mpi

    implicit none

    integer :: array(5,5)
    integer :: error, id, numproc, newtype, i, j, count, nbrbot, nbrtop
    
    integer :: new_comm, ndims, dims(2), coords(2)
    logical :: periods(2), reorder

    call MPI_init(error)
    call MPI_Comm_size(MPI_COMM_WORLD, numproc, error)
    call MPI_Comm_rank(MPI_COMM_WORLD, id, error)


    dims = 0
    periods = .false.
    reorder = .true.
    ndims = 2

    ! create cart topology based on number of processors and create new comm group
    call MPI_dims_create(numproc, ndims, dims, error)
    call MPI_cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, new_comm, error)

    !get coords of each processor
    !MPI_Cart_coords(comm, rank, maxdims, coords())

    call MPI_cart_get(new_comm, 1, dims, periods, coords, error)
    print*,coords(1),coords(2),id
    ! call MPI_cart_get(new_comm, 2, dims, periods, coords, error)
    ! print*,coords(1),coords(2),id
    nbrbot = coords(1)
    nbrtop = coords(2)
    if(id==1)print*,nbrtop, nbrbot

    ! call mpi_cart_shift(comm, direction, shift, src, dest, error)
    call mpi_cart_shift(new_comm, 0, 1, nbrbot, nbrtop, error)
    if(id==1)print*,nbrtop, nbrbot




    ! dims = 0
    ! periods = .false.
    ! reorder = .true.
    ! comm = MPI_Comm_world
    ! ndims = 2

    ! call MPI_dims_create(numproc, ndims, dims, error)
    ! call MPI_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
    ! ! call MPI_Comm_rank(new_comm, id, error)

    ! call MPI_Cart_sub(new_comm, [.true.,.false.], min_comm, error)
    ! call MPI_Comm_rank(min_comm, col_id, error)     
    ! call MPI_Cart_coords(min_comm, col_id, 1, col_coords, error)
    ! print*,col_id, col_coords

    call MPI_Finalize(error)
call exit(0)




    array = 0

    if(id ==0)then
        count = 1
        do i = 1,5
            do j = 1, 7
                array(i,j) = count
                count = count +1
            end do
        end do
    end if
    ! call MPI_TYPE_VECTOR(count, blocklength, stride, oldtype, newtype, error)
    call MPI_TYPE_VECTOR(1, 5, 5, MPI_INTEGER, newtype, error)
    call mpi_type_commit(newtype, error)

    call mpi_barrier(MPI_COMM_WORLD, error)

    call mpi_bcast(array(1,3), 1, newtype, 0, MPI_COMM_WORLD, error)
    if(id==1)then
        do i = 1,7
            print'(7(I2.1,1X))',array(:,i)
        end do
    end if
    call MPI_Finalize(error)
end program tester