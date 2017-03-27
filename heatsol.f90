program heatersoln

    implicit none
    

    real, allocatable :: array(:,:,:)
    real              :: pi, fac, x, y, z, delx, time, pi2, pi3
    integer           :: i, j, k, n, m, p, u, pts


    pts = 26
    delx = 1. / real(pts)
    allocate(array(pts+1,pts+1,pts+1))

    pi = 4.*atan(1.)
    pi2 = pi*pi
    pi3 = pi2 * pi
    array = 0.

    x = 0.
    y = 0.
    z = 0.      

    time = 0.1
        x = 0.
        do i = 1, pts+1
            print*,x
            y = 0.
            do j = 1, pts+1
                z = 0.
                do k = 1, pts+1
                    do n = 1, 11,2
                        do m = 1, 11,2
                            do p = 1, 11,2
                                fac = 37.*64./(pi3*real(n*m*p))
array(i,j,k) = array(i,j,k)+fac * sin(n*pi*x)*sin(m*pi*y)*sin(p*pi*z)*exp(-pi2*(m+n+p)*1.)
                            end do
                        end do
                    end do
                z = z + delx
                end do
            y = y + delx
            end do
        x = x + delx
        end do

    inquire(iolength=i)array
    open(newunit=u, file='soltemp.dat', access='direct',form='unformatted',status='replace',recl=i)
    write(u,rec=1)array
end program heatersoln