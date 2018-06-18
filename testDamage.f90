program testDamage

    implicit none
    

    real    :: temp(208,208,208),a,g,r,delt,tissue(208,208,208),time
    integer :: x,y,z, numpoints,u,t

    temp = 0.
    tissue = 0.

    numpoints = 208
    a = 3.1d98
    g = 6.27d5
    r = 8.314
    delt = 7.94d-4


    open(newunit=u,file="data/deposit/temp-27-0.222.dat", access="stream",form="unformatted")
    read(u)temp
    close(u)

    time = 0.
    temp = temp/maxval(temp)
    temp = temp * 65.
    temp = temp + 273
    do t = 1, 1260
        do z = 1, numpoints
            do y = 1, numpoints
                do x = 1, numpoints
                    if(temp(x, y, z) >= 44+273)then
                        tissue(x, y, z) = tissue(x, y, z) + delt*A*exp(-G/(R*(temp(x,y,z))))
                    end if
                end do
            end do
        end do
    time = time + delt
    print*,time
    end do



    open(newunit=u,file="damage.dat", access="stream",form="unformatted")
    write(u)tissue
    close(u)

end program testDamage