program threedfinite

    implicit none
    
    real              :: u_xx, u_yy, u_zz, h, kappa, rho, c_heat, t_air, alpha, rx, ry, rz, eta, eps, sigma, t_air4
    real              :: delt, time, dx, dy, dz, betax, betay, betaz, gammax, gammay, gammaz
    real, allocatable :: T0(:,:,:)
    integer           :: N, i, j, k, p, u, size_x, size_y, size_z, o


    N = 50

    h = 10.
    kappa = 0.00209 !0.0056 ! W/cm K
    rho = 1.07 ! g/cm^3
    c_heat = 3.4 !J/g K
    t_air = 25.+273
    t_air4 = t_air**4
    eps = 0.98
    sigma = 5.670373e-8

    rho = rho / 1000.
    c_heat = c_heat*1000.
    alpha = kappa / (rho * c_heat)

    size_x = N
    size_y = N
    size_z = N

    dx = 1. / (size_x + 1)
    dy = 1. / (size_y + 1)
    dz = 1. / (size_z + 1)

    !allocate mesh
    allocate(T0(0:size_x+1, 0:size_y+1, 0:size_z+1))

    !heat eqn constants
    betax  = 1. + (dx * h / kappa)
    gammax = (dx * h * T_air) /kappa
    delt = dx**2 / (6. * alpha * betax) !assume square

    betay  = 1. + (dy * h / kappa)
    betaz  = 1. + (dz * h / kappa)
    gammay = (dy * h * T_air) /kappa
    gammaz = (dz * h * T_air) /kappa

    eta = eps*sigma*dx*dy*betax

    rx = alpha * delt/dx**2
    ry = alpha * delt/dy**2
    rz = alpha * delt/dz**2
    
    t0 = 37.
    t0(:,:,0) = 37.  ! bottom face
    t0(:,0,:) = 37. ! front face
    t0(:,N+1,:) = 37.  ! back face
    t0(N+1,:,:) = 37.  ! side face
    t0(0,:,:) = 25.    ! side face
    t0(:,:,N+1) = 37.  ! top face 
    t0 = t0 + 273.

    o = nint(10./delt)

    do p = 1, o
        do k = 1, size_z
            do j = 1, size_y
                do i = 1,size_x
                    if(i == 1)then
                        u_xx = (1.-2.*rx*betax) * t0(i,j,k) + (2. *rx * t0(i+1,j,k)) + (2. * rx * gammax) &
                               - 2.*rx*eta*(t0(i,j,k)**4-T_air4)
                        u_yy = ry*(t0(i, j - 1, k    ) - 2. * t0(i, j, k) + t0(i, j + 1, k))
                        u_zz = rz*(t0(i, j,     k - 1) - 2. * t0(i, j, k) + t0(i, j, k + 1))
                        t0(i,j,k) = u_xx + u_yy + u_zz
                    else
                        u_xx = rx*(t0(i - 1, j, k    ) - 2. * t0(i, j, k) + t0(i + 1, j, k))
                        u_yy = ry*(t0(i, j - 1, k    ) - 2. * t0(i, j, k) + t0(i, j + 1, k))
                        u_zz = rz*(t0(i, j,     k - 1) - 2. * t0(i, j, k) + t0(i, j, k + 1))
                        t0(i,j,k) = t0(i,j,k) + (u_xx + u_yy + u_zz)
                    end if
                end do
            end do
        end do
    time = time + delt
    end do
    print*,time,p,n

    t0 = t0 - 273.

   inquire(iolength=i)t0(1:n,1:n,1:n)
   open(newunit=u,file='temperature.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) t0(1:n,1:n,1:n)
   close(u)

end program threedfinite