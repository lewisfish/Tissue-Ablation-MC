program threedfinite

    use mpi

    implicit none
    
    real              :: u_xx, u_yy, diagx, diagy, weightx, weighty
    real              :: delt, time, k0, hx, hy, hx2, hy2
    real, allocatable :: T0(:,:), T(:,:)
    integer           :: N, i, j, p, u, size_x, size_y, o


    N = 50

    k0 = 1.

    size_x = N
    size_y = N

    hx = 1. / (size_x + 2)
    hy = 1. / (size_y + 2)

    hx2 = 1./hx**2
    hy2 = 1./hy**2.

    delt = 0.125 * (min(hx,hy)**3)/k0
    
    !allocate mesh
    allocate(T(0:size_x+1, 0:size_y+1))
    allocate(T0(0:size_x+1, 0:size_y+1))

    t = 0.
    t0 = 37.
    t0(:,0) = 0. ! front face
    t0(:,N+1) = 0.  ! back face
    t0(N+1,:) = 0.  ! side face
    t0(0,:) = 0.    ! side face

   inquire(iolength=i)t0(1:N,1:N)

   open(newunit=u,file='init.dat',access='direct',status='REPLACE',form='unformatted',&
   recl=i)
   write(u,rec=1) t0(1:N,1:N)
   close(u)


    diagx = -2. + hx*hx/(3.*k0*delt)
    diagy = -2. + hy*hy/(3.*k0*delt)

    weightx = k0*delt/(hx*hx)
    weighty = k0*delt/(hy*hy)

    o = int(0.1/delt)

    do p = 1, o
        write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
 & " Percent Complete: ", (real(p)/real(o))*100.0, "%"
        do i = 1, size_x
            do j = 1,size_y
                u_xx = (t0(i+1,j)     - 2.*t0(i,j) + t0(i-1, j))  * hx2
                u_yy = (t0(i,  j+1)   - 2.*t0(i,j) + t0(i,   j-1))  * hy2
                t0(i,j) = t0(i,j) + k0*delt*(u_xx + u_yy)
            end do
        end do
    time = time + delt
    end do
    print*,
    print*,time,p,n

   open(newunit=u,file='temperature.dat')!,access='direct',status='REPLACE',form='unformatted',&
    do i =1,n
        write(u,*)(t0(i,j),j=1,n)
    end do
    close(u)

end program threedfinite