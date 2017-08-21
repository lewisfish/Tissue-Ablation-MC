program threedfinite

    implicit none
    
    real              :: u_xx, u_yy, betax, gammax, t_air, h, betay, gammay, T_air4, eta
    real              :: delt, time, alpha, dx, dy, kappa, rho, c_heat, rx, ry, eps, sigma
    real, allocatable :: T0(:,:)
    integer           :: N, i, j, p, u, size_x, size_y, o


    N = 100

    h = 10.
    kappa = 0.00209! W/cm K !0.0056 ! W/cm C
    rho = 1.07 ! g/cm^3
    c_heat = 3.4 !J/g K
    t_air = 25.+273
    T_air4 = t_air**4
    eps = 0.98
    sigma = 5.670373e-8

    rho = rho / 1000.
    c_heat = c_heat*1000.
    alpha = kappa / (rho * c_heat)


    size_x = N
    size_y = N

    dx = 1. / (size_x + 1)
    dy = 1. / (size_y + 1)
    
    !allocate mesh
    allocate(T0(0:size_x+1, 0:size_y+1))

    !heat eqn constants
    betax  = 1. + (dx * h / kappa)
    gammax = (dx * h * T_air) /kappa
    delt = dx**2 / (4. * alpha * betax) !assume square

    betay  = 1. + (dy * h / kappa)
    gammay = (dy * h * T_air) /kappa
    eta = eps*sigma*dx*dy*betax

    rx = alpha * delt/dx**2
    ry = alpha * delt/dy**2

    !init conditions
    t0        = 37.+273
    t0(:,0)   = 37.+273    ! front face
    t0(:,N+1) = 37.+273    ! back face
    t0(N+1,:) = 37.+273    ! side face
    t0(0,:)   = 37.+273    ! side face

    o = nint(10./delt)
    time = 0.

    do p = 1, o
        do j = 1, size_y
            do i = 1,size_x
                if(i == 1)then
                    u_xx = (1.-2.*rx*betax) * t0(i,j) + (2. * rx * t0(i+1,j)) + (2. * rx * gammax) - 2.*rx*eta*(t0(i,j)**4-T_air4)
                    u_yy = ry*(t0(i,    j - 1) - 2. * t0(i, j) + t0(i    ,j + 1)) 
                    t0(i,j) = u_xx + u_yy
                else
                    u_xx     = rx*(t0(i - 1,j    ) - 2. * t0(i, j) + t0(i + 1,j    ))
                    u_yy     = ry*(t0(i,    j - 1) - 2. * t0(i, j) + t0(i    ,j + 1)) 
                    t0(i, j) = t0(i, j) + u_xx + u_yy
                end if
            end do
        end do
    time = time + delt
    end do
    print*,
    print*,time,p

   open(newunit=u,file='temperature.dat')
    do j =1,n
        write(u,*)(t0(i,j)-273.,i=1,n)
    end do
    close(u)

end program threedfinite