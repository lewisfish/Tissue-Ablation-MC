program p

    use shrinkarray

    implicit none
    

    real :: old(50,50,50), new(200,200,200)
    integer :: i,j,k,u

old = 0.
new =0

    do i = 1, 50
        do j =1, 50
            do k = 1, 50
                if(i > 20 .and. j < 40 .and. k == 20)old(i,j,k) = 3.14
                if(i < 20)old(i,j,k) = 10.
            end do
        end do
    end do

    call unshrink(old, new)

    inquire(iolength=i)old
    open(newunit=u,file='oldshrink.dat',access='direct',form='unformatted',recl=i)
    write(u,rec=1)old
    close(u)

    inquire(iolength=i)new
    open(newunit=u,file='newshrink.dat',access='direct',form='unformatted',recl=i)
    write(u,rec=1)new
    close(u)
end program p