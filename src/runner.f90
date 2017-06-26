program run

    ! use mpi
    use heat

    implicit none
    

    ! double precision :: jmean(104,104,104)
    ! real:: j(104,104,104), t(104,104,104)
    ! integer :: n, counter, u

    real, allocatable :: jmean(:,:,:), tissue(:,:,:), temp(:,:,:)
    integer :: numpoints, counter, u
    logical :: flag 


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
    ! do i = 1, 3
        call heat_sim_3d(jmean, tissue, temp, numpoints, counter, flag)
    ! end do
end program run