program cart

    use mpi

    implicit none
    
    integer :: error, numproc, comm_old, comm_new, id, ndims, dims(2), coords(2), n, i, u
    integer :: start, finish, xcell, ycell, x_domains, y_domains, size_x, size_y, j, k
    integer, allocatable :: xe(:), xs(:), ye(:), ys(:)
    logical :: periods(2), reorder

    real, allocatable :: temp(:,:),fin(:,:),finGLOBAL(:,:)

    n = 20
    allocate(fin(n, n),finGLOBAL(n, n))
    ndims = 2
    dims = 0

    reorder = .true.
    periods = .false.
    comm_old = MPI_Comm_world

    call mpi_init(error)
    call MPI_Comm_size(comm_old, numproc, error)        !init mpi
    call MPI_Comm_rank(comm_old, id, error)

    allocate(xs(0:numproc-1))       !arrays for coords of sub grids
    allocate(xe(0:numproc-1))
    allocate(ys(0:numproc-1))
    allocate(ye(0:numproc-1))

    if(numproc /= 4)then
        print*,'numproc != 4'
        stop
    end if

    call MPI_dims_create(numproc, ndims, dims, error)
    call MPI_cart_create(comm_old, ndims, dims, periods, reorder, comm_new, error)  !init topo

    x_domains = dims(1)
    y_domains = dims(2)
    size_x = n
    size_y = n

    CALL MPI_Comm_rank(comm_new,id,error)
    CALL MPI_Cart_get(comm_new,ndims,dims,periods,coords,error)         !divide cart grid into sub grids

    call mpi_barrier(comm_new, error)

    allocate(temp(2*(n/numproc), 2*(n/numproc)))

    if(id == 0)then
        temp = 0.
    elseif(id == 1)then
        temp = 1.
    elseif(id == 2)then
        temp = 2.
    else
        temp = 3.
    end if

    call mpi_barrier(comm_new, error)

    xcell = size_x/x_domains
    ycell = size_y/y_domains

    xs(0:y_domains-1) = 1
    xe(0:y_domains-1) = xs(0:y_domains-1)+xcell-1
    
    do i=1,(x_domains-1)
        xs(i*y_domains:(i+1)*(y_domains)-1) = xs((i-1)*y_domains:i*(y_domains)-1)+xcell
        xe(i*y_domains:(i+1)*(y_domains)-1) = xs(i*y_domains:(i+1)*(y_domains)-1)+xcell-1
    end do
    ! print*,x_domains,y_domains
    ! if(id==0)print*,xe
    !     if(id==0)print*,xs


    do i=1,x_domains
        ys(i*y_domains-1) = 1
        ye(i*y_domains-1) = 2+ycell-2
    end do
    do j=1,x_domains
        do k=0,y_domains-2
            ys(j*y_domains-2-k) = ys(j*y_domains-2-k+1)+ycell
            ye(j*y_domains-2-k) = ys(j*y_domains-2-k)+ycell-1
        end do
    end do

    fin(xs(id):xe(id),ys(id):ye(id)) = temp

    call MPI_REDUCE(fin, finGLOBAL, n*n,mpi_real, MPI_SUM,0,comm_new,error)

    if(id == 0)then
        open(newunit=u,file='test.dat')
        do i = 1, n
            write(u,*) finGLOBAL(i,:)
        end do
    end if
    call mpi_finalize(error)

end program cart
