Module Heat

    implicit none

    private
    public :: heat_sim_3D

    contains

    subroutine heat_sim_3D(oldj, oldt, N, counter)

        use shrinkarray
        use constants, only : fileplace

        implicit none
        
        real,   intent(IN)    :: oldj(:,:,:)
        real,   intent(INOUT) :: oldt(:,:,:)
        integer, intent(IN)   :: N, counter

        real              :: u_xx, u_yy, u_zz, diagx, diagy, diagz, weightx, weighty, weightz
        real              :: delt, time, k0, hx, hy, hz, hx2, hy2, hz2
        real              :: rho, kappa, c_heat
        real, allocatable :: T0(:,:,:), T(:,:,:), jmean(:,:,:), tissue(:,:,:)
        integer           :: i, j, k, p, u, size_x, size_y, size_z, o,q,w
        character(len=3)  :: fn

        allocate(jmean(N,N,N))
        allocate(tissue(n,n,n))

        if(size(jmean,1) /= size(oldj,1))then
            call shrink(oldj, jmean)
        else
            jmean = oldj
        end if

        if(size(tissue,1) /= size(oldt,1))then
            call shrink(oldt, tissue)
        else
            tissue = oldt
        end if


        jmean = jmean*1.

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

        t = 0.
        t0 = 37.
        t0(:,:,0) = 37.  ! bottom face
        t0(:,0,:) = 37. ! front face
        t0(:,N+1,:) = 37.  ! back face
        t0(N+1,:,:) = 37.  ! side face
        t0(0,:,:) = 37.    ! side face
        t0(:,:,N+1) = 25.  ! top face 

        diagx = -2. + hx*hx/(3.*k0*delt)
        diagy = -2. + hy*hy/(3.*k0*delt)
        diagz = -2. + hz*hz/(3.*k0*delt)

        weightx = k0*delt/(hx*hx)
        weighty = k0*delt/(hy*hy)
        weightz = k0*delt/(hz*hz)

        o = int(100.e-3/delt)
        q = 0
        inquire(iolength=w)tissue
        print*,o

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
        ! if(mod(p,1000) == 0)then
        !     write(fn,'(I3)') q
        !     open(newunit=u,file='tissue_'//trim(adjustl(fn))//'.dat',access='direct',status='REPLACE',form='unformatted',&
        !     recl=w)
        !     write(u,rec=1) tissue
        !     close(u)
        !     q = q + 1
        ! end if
        end do
        print*,
        print*,time,p,n

        write(fn,'(I3.3)') counter
        inquire(iolength=w)tissue
        open(newunit=u,file=trim(fileplace)//'jmean/tissue-'//trim(fn)//'.dat',access='direct',status='REPLACE', &
            form='unformatted', recl=w)
        write(u,rec=1) tissue
        close(u)


        if(size(tissue,1) /= size(oldt,1))then
            call unshrink(tissue, oldt)
        end if


        inquire(iolength=w)t0
        open(newunit=u,file=trim(fileplace)//'jmean/temperature-'//trim(fn)//'.dat',access='direct',status='REPLACE', &
            form='unformatted', recl=w)
        write(u,rec=1) t0
        close(u)

        deallocate(T)
        deallocate(T0)
        deallocate(jmean)

   end subroutine heat_sim_3D

    subroutine Arrhenius(Temp, delt, tissue)

        implicit none

        real, intent(IN)    :: temp(:,:,:), delt
        real, intent(INOUT) :: tissue(:,:,:)
        double precision :: a, g, r
        integer          :: x, y, z

        a = 3.1d91!2.9e27
        g = 6.28e5!2.4e5
        r = 8.314

        do x = 1, size(temp,1)
            do y = 1, size(temp,1)
                do z = 1, size(temp,1)
                    if(temp(x, y, z) >= 44)then
                        tissue(x, y, z) = tissue(x, y, z) + delt*A*exp(-G/(R*(temp(x,y,z)+273)))
                        ! print*,delt*A*exp(-G/(R*(temp(x,y,z)+273)))
                    end if
                end do
            end do
        end do

    end subroutine  Arrhenius
end Module Heat