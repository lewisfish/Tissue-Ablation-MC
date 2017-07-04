program run

    use mpi
    use heat

    implicit none
    

    ! double precision :: jmean(104,104,104)
    ! real:: j(104,104,104), t(104,104,104)
    ! integer :: n, counter, u

    real, allocatable :: jmean(:,:,:), tissue(:,:,:), temp(:,:,:)
    integer :: numpoints, counter, u, i
    logical :: flag 

    ! mpi variables
    integer :: new_comm, error, right, left, id, numproc, dims(2), ndims, tag, recv_status(mpi_status_size), comm
    logical :: periods(1), reorder

    !init mpi
    comm = mpi_comm_world
    call mpi_init(error)
    call mpi_comm_size(comm, numproc, error)

    !setup topology variables
    tag = 1
    dims = 0
    periods = .false.
    reorder = .true.
    ndims = 1

    !create cartesian topology on cpu
    call mpi_dims_create(numproc, ndims, dims, error)
    call mpi_cart_create(comm, ndims, dims, periods, reorder, new_comm, error)
    call mpi_comm_rank(new_comm, id, error)

    !get neighbours
    call mpi_cart_shift(new_comm, 0, 1, left, right, error)


    counter = 90
    numpoints = 200
    flag = .true.

    allocate(jmean(numpoints,numpoints,numpoints), tissue(numpoints,numpoints,numpoints), &
             temp(0:numpoints+1,0:numpoints+1,0:numpoints+1))

    jmean = 0.
    tissue = 0.
    temp = 0.

    !   to do:
    !         integrate with full mc code
    !         check running multiple times works i.e in a loop
    !         optimise memory and shizz in loop
    !         tidy up
    !         email matthew
    !           


    open(newunit=u,file='testjmean.dat',access='stream',form='unformatted')
    read(u)jmean
    close(u)
    jmean=jmean*100.
    do i = 1, 3

        call heat_sim_3d(jmean, tissue, temp, numpoints, flag, id, numproc, error, new_comm, tag, recv_status, right, left, i)
    end do
    call mpi_finalize(error)

end program run