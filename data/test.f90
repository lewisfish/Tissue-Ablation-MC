program neighbors
!
! Run with 12 processes
!
    use mpi
    
    implicit none

    integer :: err, rank, size
    integer :: vu
    integer :: dim(2)
    logical :: period(2),reorder
    integer :: up,down,right,left

    CALL MPI_INIT(err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err)

    dim(1)=4
    dim(2)=3
    period(1)=.true.
    period(2)=.false.
    reorder=.true.

    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dim,period,reorder,vu,err)

    if(rank.eq.9) then
        call MPI_CART_SHIFT(vu,0,1,left,right,err)
        call MPI_CART_SHIFT(vu,1,1,up,down,err)
        print*,'P:',rank,' neighbors (rdlu)are',right,down,left,up
    end if

    CALL MPI_FINALIZE(err)
end program