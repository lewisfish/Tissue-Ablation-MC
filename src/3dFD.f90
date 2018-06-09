Module Heat

    implicit none

    private
    public :: heat_sim_3D

    contains

subroutine heat_sim_3D(jmean, tissue, temp, numpoints, flag, id, numproc, new_comm, right, left,counter, Q)

        use mpi_f08
        use utils,     only : str
        use constants, only : fileplace, nxg, nyg, nzg

        implicit none
        
        real,           intent(IN)    :: jmean(:,:,:)
        integer,        intent(IN)    :: numpoints, right, left, id, numproc,counter
        type(mpi_comm), intent(IN)    :: new_comm
        real,           intent(INOUT) :: tissue(:,:,:), q(nxg,nyg,nzg), temp(0:numpoints+1,0:numpoints+1,0:numpoints+1)
        logical,        intent(INOUT) :: flag

        type(MPI_Status) :: recv_status

        !heat variables
        real              :: u_xx, u_yy, u_zz, delt, time, alpha, dx, dy, dz
        real              :: rho, kappa, c_heat, laserOn, pulselength, repetitionRate_1, pulseCount, repetitionCount
        real              :: betax, gammax, betay, gammay, betaz, gammaz, h, t_air, t_air4, rx, ry, rz, eps, sigma, eta
        real              :: tim_fin, tim_srt, coeff, Q_vapor, L_w, mass_voxel
        real, allocatable :: T0(:,:,:), jtmp(:,:,:), Ttmp(:,:,:), qtmp(:,:,:)
        integer           :: i, j, k, p, u, size_x, size_y, size_z, xi, yi, zi, xf, yf, zf, lo, N, o, tag=1
        logical :: laser_flag

        if(id == 0)then
            !do decomp
            if(mod(numpoints,numproc) /= 0)then
                print('(I2.1,a,I2.1)'),numpoints,' not divisable by : ', numproc
                call mpi_finalize()
                error stop
            else
                N = numpoints / numproc
            end if
        end if

        !send N to all processes
        call MPI_Bcast(N, 1, MPI_integer ,0 , new_comm)

        !init heat variables for medium
        kappa = 0.00209 !0.0056 ! W/cm K
        rho = 1.07 ! g/cm^3
        c_heat = 3.4 !J/g K
        eps = 0.98
        sigma = 5.670373e-8
        t_air = 25.+273.
        t_air4 = t_air**4
        h = 10.
        L_w = 2256.d3 !J/kg
        !assume 70% water
        mass_voxel = 997.*(2.*0.55d-2/nxg)**3.!kg cm-3
        Q_vapor = L_w * mass_voxel !joules
        print*,Q_vapor


        rho = rho/1000.
        c_heat = c_heat*1000.

        alpha = kappa / (rho * c_heat)
        !alpha = 1.07d-3/(0.0056*3400.)!1.

        !init grid
        size_x = numpoints
        size_y = numpoints
        size_z = N

        dx = 1. / (numpoints + 2)
        dy = 1. / (numpoints + 2)
        dz = 1. / (numpoints + 2)  

        !init heat variables for air
        ! kappa = 0.0260e-2 ! W/cm C
        ! rho = 0.0012 ! g/cm^3
        ! c_heat = 1.012 !J/g C

        ! rho = rho/1000.
        ! c_heat = c_heat*1000.

        xi = 1 
        xf = size_x

        yi = 1
        yf = size_y

        zi = 1
        zf = size_z

        !allocate mesh
        allocate(T0, source=temp)
        allocate(jtmp(numpoints, numpoints, zi:zf))
        allocate(Ttmp, source=tissue)
        allocate(qtmp, source=Q)
        Ttmp=0.

        betax = 1. + (dx*h/kappa)
        betay = 1. + (dy*h/kappa)
        betaz = 1. + (dz*h/kappa)

        delt = dx**2/(6.*alpha*betax)

        gammax = dx*h*t_air/kappa
        gammay = dy*h*t_air/kappa
        gammaz = dz*h*t_air/kappa

        eta = eps*sigma*dx*dy*betax
        coeff = alpha*delt/kappa

        rx = alpha * delt/(dx**2) 
        ry = alpha * delt/(dy**2) 
        rz = alpha * delt/(dz**2) 

        jtmp = 0.
        call mpi_scatter(jmean, size(jmean(:,:,zi:zf)), mpi_double_precision, jtmp(:,:,zi:zf), size(jtmp(:,:,zi:zf)), &
                         mpi_double_precision, 0, new_comm)

        time = 0.
        laserOn = 1.
        pulselength = 200d-3
        repetitionRate_1 = 1.d0/1.d0

        if(pulselength < delt)then
            if(id==0)print*,"pulselength smaller than timestep, adjusting..."
            delt = pulselength / 2.
        end if

        pulseCount = 0.
        repetitionCount = 0.
        laser_flag = .true.
        o = int(1./delt)
        lo = o / 100
        if(id == 0)print*,delt,o
        call cpu_time(tim_srt)
        do p = 1, o
            
            if(mod(p,lo) == 0 .and. id == 0) print*,p!write(*,FMT="(A8,t21)",advance='no') achar(13)//str(real(p)/real(o)*100., 5)//' %'
            if(p == 500 .and. id == 0)then
                call cpu_time(tim_fin)
                print*, 'est. time for loop ~ ',str((tim_fin-tim_srt)*o/100.,5)//'s'
            end if
            do k = zi, zf
                do j = yi, yf
                    do i = xi, xf
                        if(k == zf .and. id == numproc-1)then!B.Cs
                        u_zz = (1.-2.*rz*betax) * t0(i,j,k) + (2. *rz * t0(i,j,k+1)) + (2. * rz * gammax) &
                               - 2.*rz*eta*(t0(i,j,k)**4-T_air4)
                        else
                            u_zz = rz*t0(i, j,     k - 1) + (1.-2.*rz) * t0(i, j, k) + rz*t0(i, j, k + 1)
                        end if
                        u_yy = ry*t0(i, j - 1, k    ) + (1.-2.*ry) * t0(i, j, k) + ry*t0(i, j + 1, k)
                        u_xx = rx*t0(i - 1, j, k    ) + (1.-2.*rx) * t0(i, j, k) + rx*t0(i + 1, j, k)
                        ! if(t0(i,j,k) >= 100.+273 .and. t0(i,j,k) <= 101.+273 .and. qtmp(i,j,k) < Q_vapor)then
                        !     qtmp(i,j,k) = qtmp(i,j,k) + laserOn*delt*coeff*jtmp(i,j,k)
                        !     t0(i,j,k) = 100. + 273.
                        ! else
                            t0(i,j,k) = (u_xx + u_yy + u_zz)/3. + laserOn*coeff*jtmp(i,j,k)
                        ! end if
                    end do
                end do
            end do

            ! where(q >= Q_vapor)
            !     tissue = 7.  
            ! end where

            !send_recv data to right
            call MPI_Sendrecv(t0(:,:,zf), size(t0(:,:,zf)), mpi_double_precision, right, tag, &
                              t0(:,:,zf+1), size(t0(:,:,zf+1)), mpi_double_precision, right, tag, &
                              new_comm, recv_status)
            
            !send_recv data to left
            call MPI_Sendrecv(t0(:,:,zi), size(t0(:,:,zi)), mpi_double_precision, left, tag, &
                              t0(:,:,zi-1), size(t0(:,:,zi-1)), mpi_double_precision, left, tag, &
                              new_comm, recv_status)

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
            pulseCount = pulseCount + delt
            repetitionCount = repetitionCount + delt
            time = time + delt
           call Arrhenius(t0(xi:xf,yi:yf,zi:zf), time, tissue(xi:xf,yi:yf,zi:zf), zi, zf, numpoints)

        end do

        call mpi_gather(t0(:,:,zi-1:zf+1), size(t0(:,:,zi-1:zf+1)), mpi_double_precision, &
                        temp(:,:,zi-1:zf+1), size(temp(:,:,zi-1:zf+1)), mpi_double_precision, 0, new_comm)

        call mpi_gather(tissue(:,:,zi:zf), size(tissue(:,:,zi:zf)), mpi_double_precision, &
                        Ttmp(:,:,zi:zf), size(Ttmp(:,:,zi:zf)), mpi_double_precision, 0, new_comm)

        tissue = Ttmp
        call mpi_allreduce(qtmp, Q, size(Q), mpi_double_precision, mpi_sum, new_comm)

        !send data to master process and do I/O
        if(id == 0)then
            open(newunit=u,file=trim(fileplace)//'deposit/temp-'//str(counter)//'.dat',access='stream', &
                 form='unformatted',status='replace')
            write(u)temp - 273.
            close(u)

            open(newunit=u,file=trim(fileplace)//'deposit/tissue-'//str(counter)//'.dat',access='stream', &
                 form='unformatted',status='replace')
            write(u)tissue
            close(u)

            open(newunit=u,file=trim(fileplace)//'deposit/Q-'//str(counter)//'.dat',access='stream', &
                 form='unformatted',status='replace')
            write(u)Q
            close(u)
        end if
        ! call mpi_barrier(new_comm)
        deallocate(t0, jtmp, Ttmp, qtmp)
   end subroutine heat_sim_3D

    subroutine Arrhenius(temp, delt, tissue, zi, zf, numpoints)

        implicit none

        integer, intent(IN)    :: zi,zf, numpoints
        real,    intent(IN)    :: temp(numpoints,numpoints,zi:zf), delt
        real,    intent(INOUT) :: tissue(numpoints,numpoints,zi:zf)

        double precision :: a, g, r
        integer          :: x, y, z

        a = 3.1d91!2.9e27
        g = 6.28e5!2.4e5
        r = 8.314

        where(temp >= 44+273 .and. temp < 100+273)
            tissue = 0.!tissue + time*A*exp(-G/(R*(temp+273.)))
        elsewhere(temp >= 100+273 .and. temp < 200.+273)
            tissue = 1.
        elsewhere(temp >= 200.+273. .and. temp < 300.+273)
            tissue = 2.
        elsewhere(temp >= 300.+273. .and. temp < 400.+273)
            tissue = 3.
        elsewhere(temp >= 400.+273. .and. temp < 500.+273)
            tissue = 4.
        elsewhere(temp >= 500.+273. .and. temp < 600.+273)
            tissue = 5.
        elsewhere(temp >= 600.+273. .and. temp < 700.+273)
            tissue = 6.
        end where
        ! do z = zi, zf
        !     do y = 1, numpoints
        !         do x = 1, numpoints
        !             if(temp(x, y, z) >= 44)then
        !                 tissue(x, y, z) = tissue(x, y, z) + delt*A*exp(-G/(R*(temp(x,y,z)+273)))
        !             end if
        !         end do
        !     end do
        ! end do

    end subroutine  Arrhenius
end Module Heat