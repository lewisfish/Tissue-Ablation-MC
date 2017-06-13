program threedfinite

    use shrinkarray

    implicit none
    
    real              :: u_xx, u_yy, u_zz, diagx, diagy, diagz, weightx, weighty, weightz
    real              :: delt, time, k0, hx, hy, hz, hx2, hy2, hz2, jmean(50,50,50), old(200,200,200)
    real              :: rho, kappa, c_heat
    real*8:: popo(200,200,200)
    real, allocatable :: T0(:,:,:), T(:,:,:), tissue(:,:,:)
    integer           :: N, i, j, k, p, u, size_x, size_y, size_z, o,q,w
    character(len=3)  :: fn


    N = 50

    inquire(iolength=i)popo
    open(newunit=u,file='jmean.dat',access='direct',form='unformatted',recl=i)
    read(u,rec=1)popo
    close(u)

    old = popo

    call shrink(old, jmean)
    jmean = jmean*10.
    inquire(iolength=i)jmean
    open(newunit=u,file='old.dat',access='direct',form='unformatted',recl=i,status='REPLACE')
    write(u,rec=1)jmean 
    close(u)


    kappa = 0.0056 ! W/cm C
    rho = 1.07 ! g/cm^3
    c_heat = 3.4 !J/g C

    rho = rho/1000.
    c_heat = c_heat*1000.

    k0 = kappa / (rho * c_heat)
    !k0 = 1.07d-3/(0.0056*3400.)!1.

    size_x = N
    size_y = N
    size_z = N

    hx = 1. / (size_x + 2)
    hy = 1. / (size_y + 2)
    hz = 1. / (size_z + 2)  

    hx2 = 1./hx**2
    hy2 = 1./hy**2.
    hz2 = 1./hz**2.

    delt = 0.125 * (min(hx,hy,hz)**3)/k0
    
    !allocate mesh
    allocate(T(0:size_x+1, 0:size_y+1, 0:size_z+1))
    allocate(T0(0:size_x+1, 0:size_y+1, 0:size_z+1))
    allocate(tissue(N, N, N))

    tissue = 0.

    t = 0.
    t0 = 37.
    t0(:,:,0) = 37.  ! bottom face
    t0(:,0,:) = 37. ! front face
    t0(:,N+1,:) = 37.  ! back face
    t0(N+1,:,:) = 37.  ! side face
    t0(0,:,:) = 37.    ! side face
    t0(:,:,N+1) = 25.  ! top face 

   inquire(iolength=i)t0(1:N,1:N,1:N)

   open(newunit=u,file='init.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) t0(1:N,1:N,1:N)
   close(u)


    diagx = -2. + hx*hx/(3.*k0*delt)
    diagy = -2. + hy*hy/(3.*k0*delt)
    diagz = -2. + hz*hz/(3.*k0*delt)

    weightx = k0*delt/(hx*hx)
    weighty = k0*delt/(hy*hy)
    weightz = k0*delt/(hz*hz)


    o = int(100./delt)
    q = 0
   inquire(iolength=w)tissue

    do p = 1, o
        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
 & " Percent Complete: ", (real(p)/real(o))*100.0, "%"
        do k = 1, size_z
            do i = 1, size_x
                do j = 1,size_y
                    u_xx = (t0(i+1, j,   k)   - 2.*t0(i,j,k) + t0(i-1, j,   k))  * hx2
                    u_yy = (t0(i,   j+1, k)   - 2.*t0(i,j,k) + t0(i,   j-1, k))  * hy2
                    u_zz = (t0(i,   j,   k+1) - 2.*t0(i,j,k) + t0(i,   j,   k-1))* hz2
                    t0(i,j,k) = t0(i,j,k) + k0*delt*(u_xx + u_yy + u_zz) + k0*delt*jmean(i,j,k)/rho
                end do
            end do
        end do
    ! if(mod(p,100) == 0)then
    !     write(fn,'(I3)') q
    !    open(newunit=u,file='temp_'//trim(adjustl(fn))//'.dat',access='direct',status='REPLACE',form='unformatted',&
    !         recl=w)
    !    write(u,rec=1) t0(1:N,1:N,1:N)
    !    close(u)
    !    q = q + 1
    ! end if

    time = time + delt
    call Arrhenius(t0, delt, tissue)
    if(mod(p,1000) == 0)then
        write(fn,'(I3)') q
        open(newunit=u,file='tissue_'//trim(adjustl(fn))//'.dat',access='direct',status='REPLACE',form='unformatted',&
        recl=w)
        write(u,rec=1) tissue
        close(u)
        q = q + 1
    end if
    end do
    print*,
    print*,time,p,n

   inquire(iolength=i)t0

   open(newunit=u,file='temperature.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) t0
   close(u)


   inquire(iolength=i)tissue

   open(newunit=u,file='tissue3D.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) tissue
   close(u)


   contains
    subroutine Arrhenius(Temp, delt, tissue)

        implicit none

        real, intent(IN)    :: temp(:,:,:), delt
        real, intent(INOUT) :: tissue(:,:,:)
        real*8 :: a, g, r
        integer :: x, y, z

        a = 3.1d91!2.9e27
        g = 6.28e5!2.4e5
        r = 8.314

        do x = 1, 50
            do y = 1, 50
                do z = 1, 50
                    if(temp(x, y, z) >= 44)then
                        tissue(x, y, z) = tissue(x, y, z) + delt*A*exp(-G/(R*(temp(x,y,z)+273)))
                        ! print*,delt*A*exp(-G/(R*(temp(x,y,z)+273)))
                    end if
                end do
            end do
        end do

    end subroutine  Arrhenius
end program threedfinite