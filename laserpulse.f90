program timing

    implicit none
    
    real    :: time, laserOn, pulselength, repetitionRate_1, pulseCount, repetitionCount, delt
    logical :: laser_flag
    integer :: o, lo, p, numpoints
    real :: rho, c_heat, dx, h, kappa, alpha, betax


    rho = 1.07 ! g/cm^3
    c_heat = 3.4 !J/g K

    numpoints = 104
    rho = rho/1000.
    c_heat = c_heat*1000.
    dx = 1. / (numpoints + 2)
    kappa = 0.00209 !0.0056 ! W/cm K
    alpha = kappa / (rho * c_heat)
    h = 10.

    betax = 1. + (dx*h/kappa)
    delt = dx**2/(6.*alpha*betax)

    !pulseWidth = 1.2d-3
    !spotSize = 250d-6
    !repition rate: 0.5Hz to 5Hz

    time = 0.
    laserOn = 1.
    pulselength = 200d-3
    repetitionRate_1 = 1./5.

    if(pulselength < delt)then
        print*,"pulselength smaller than timestep, adjusting..."
        delt = pulselength / 2.
    end if

    pulseCount = 0.
    repetitionCount = 0.
    laser_flag = .true.
    o = int(1./delt)
    lo = o / 100

    do p = 1, o
        if(pulseCount >= pulseLength .and. laser_flag)then!turn laser off
            laser_flag = .false.
            laserOn = 0.
            pulseCount = 0.
            repetitionCount = 0.
        elseif((repetitionCount >= repetitionRate_1) .and. (.not.laser_flag))then
            laser_flag = .true.
            laserOn = 1.
            pulseCount = 0.
            repetitionCount = 0.
        end if
        print*,laserOn,laser_flag,time
        pulseCount = pulseCount + delt
        repetitionCount = repetitionCount + delt
        time = time + delt
    end do

end program timing