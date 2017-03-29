module shrinkarray

    implicit none
    
    contains

    subroutine shrink(old, new)

        implicit none

        double precision, intent(IN)    :: old(:,:,:)
        double precision, intent(OUT)   :: new(:,:,:)

        double precision    :: tmp
        integer :: i, j, k, p, l, o, div, n, m


        ! i,j,k iterate over new
        ! p,l,o iterate over old

        n = size(old,1)
        m = size(new,1)
        div = n/m

        new = 0.

        do i = 1, m
            do j = 1, m
                do k = 1, m
                    tmp = 0.
                            !black magic that I WILL forget...
                    do l = i*div-(div/2)-1, (i*div-(div/2)-1)+div-1
                        do o = j*div-(div/2)-1, (j*div-(div/2)-1)+div-1
                            do p = k*div-(div/2)-1, (k*div-(div/2)-1)+div-1
                                tmp = tmp + old(l,o,p)
                            end do
                        end do
                    end do
                    !divisor needs to be the same order as spatial dimensions i.e for 3D needs to be **3
                    new(i,j,k) = tmp/div**3.
                end do
            end do
        end do
    end subroutine shrink
end module shrinkarray