program ring

    implicit none

    integer :: i

    !variables for heat
    integer :: size_x, n, p, o
    real ::  alpha, dx, delt, kappa, rho, c_heat,  h, r, gamma, beta, T_air, eps, sigma, eta, a
    real, allocatable :: t0(:)


    !init heat shizz
    N = 10
    kappa = 0.00209!0.0056 ! W/cm C !.209 W m-1 K-1
    rho = 1.07 ! g/cm^3
    c_heat = 3.4 !J g-1 K-1
    T_air = (25.+273.)
    h = 10.
    sigma = 5.670373e-8 !Wm-2 K-4
    eps = 0.98
    a = 1.


    rho = rho/1000.
    c_heat = c_heat*1000.
    alpha = kappa / (rho * c_heat)

    size_x = N
    dx = 1. / (size_x + 1)
    allocate(t0(0:size_x+1))    

    beta = 1. + (dx * h / kappa)
    gamma = (dx * h * T_air) /kappa
    eta = sigma*dx*eps*beta
    delt = dx**2 / (2. * alpha * beta)
    r = (alpha * delt) / (dx**2)

    !boundary conditions
    t0        = 37.+273.
    t0(0)     = 0.+273.
    t0(n + 1) = 37.+273.

    o = int(10./delt)

    !do heat sim
    do p = 1, o
        do i = 1 , n
            if(i==1)then
                t0(i) = (1.-2.*r*beta) * t0(i) + (2. * r * t0(i+1)) +2.*r*gamma - 2.*r*eta*(t0(i)**4-T_air**4)
            else
                t0(i) = r*t0(i-1) + (1.-2.*r)*t0(i) + r*t0(i+1)
            end if
        end do
    end do
    do i = 1, n
        print*,t0(i)-273.
    end do

end program ring