Module Heat

    implicit none
    ! pointer to power function for simulation: tophat, gaussian etc
    procedure(getPwrTriangular), pointer :: getPwr => null()

    real              :: pulseCount, repetitionCount, time, laserOn, total_time, repetitionRate_1, energyPerPixel
    real              :: Power, pulselength, delt, realPulseLength
    real              :: dx, dy, dz, massVoxel, volumeVoxel
    real, allocatable :: coeff(:,:,:), kappa(:,:,:), density(:,:,:), heatcap(:,:,:), WaterContent(:,:,:), alpha(:,:,:)
    logical           :: laser_flag, pulseFlag
    integer           :: loops, pulsesToDo, pulsesDone, loops_left

    private
    public :: power, delt, energyPerPixel, laser_flag, laserOn, loops, pulseCount, pulselength, pulsesToDo, repetitionCount
    public :: repetitionRate_1, time,total_time, initThermalCoeff, heat_sim_3d, watercontent
    public :: pulseFlag, realPulseLength, getPwr

    contains

subroutine heat_sim_3D(jmean, temp, numpoints, id, numproc, new_comm, right, left, counter)

        use mpi_f08
        use utils,     only : str
        use thermalConstants, only : getSkinHeatCap, getSkinDensity, getSkinThermalCond

        implicit none
        
        !heat input variables
        real,           intent(IN)    :: jmean(:,:,:)
        integer,        intent(IN)    :: numpoints, right, left, id, numproc,counter
        type(mpi_comm), intent(IN)    :: new_comm
        real,           intent(INOUT) :: temp(0:numpoints+1,0:numpoints+1,0:numpoints+1)

        type(MPI_Status) :: recv_status

        !heat local variables
        real              :: u_xx, u_yy, u_zz, tempIncrease,kappaMinHalf, kappaPlusHalf, energyIncrease
        real              :: heatcapMinHalf, densityPlusHalf, densityMinHalf, heatcapPlusHalf, a, b, d
        real, allocatable :: T0(:,:,:), jtmp(:,:,:), tn(:,:,:)
        integer           :: i, j, k, p, size_x, size_y, size_z, xi, yi, zi, xf, yf, zf, N, tag,zsta,zfin

        !calculate size of domain
        !discretizes along z axis
        if(id == 0)then
            !do decomp
            if(mod(numpoints, numproc) /= 0)then
                print('(I2.1,a,I2.1)'),numpoints,' not divisable by : ', numproc
                call mpi_abort(new_comm, 1)
                error stop
            else
                N = numpoints / numproc
            end if
        end if

        tag = 1
        !send N to all processes
        call MPI_Bcast(N, 1, MPI_integer ,0 , new_comm)

        !init grid
        size_x = numpoints
        size_y = numpoints
        size_z = N

        xi = 1 
        xf = size_x

        yi = 1
        yf = size_y

        zi = 1
        zf = size_z

        !allocate mesh
        allocate(T0(0:numpoints+1, 0:numpoints+1, zi-1:zf+1))
        allocate(Tn(0:numpoints+1, 0:numpoints+1, zi-1:zf+1))
        allocate(jtmp(numpoints, numpoints, zi:zf))
        t0 = 0.

        !split temp/jmean/Q up over all processes
        if(id == 0)then
            do i = 1, numproc-1
                zsta = (i)*zf
                zfin = zsta + zf + 1
                call mpi_send(temp(:,:,zsta:zfin), size(temp(:,:,zsta:zfin)), mpi_double_precision, i, tag, new_comm)
            end do
            t0(:,:,:) = temp(:,:,0:zf+1)
        else
            call mpi_recv(t0, size(t0), mpi_double_precision, 0, tag, new_comm, recv_status)
        end if
        tn = t0
        jtmp = 0.
        call mpi_scatter(jmean, size(jmean(:,:,zi:zf)), mpi_double_precision, jtmp(:,:,zi:zf), size(jtmp(:,:,zi:zf)), &
                         mpi_double_precision, 0, new_comm)


        if(pulselength < delt)then
            if(id == 0)print*,"pulselength smaller than timestep, adjusting..."
            delt = pulselength / 100.
        end if

        loops_left = int(total_time/(real(loops)*delt)) - counter
        if(id == 0 .and. mod(int(total_time/(real(loops)*delt)) - counter, 100) == 0)then
            print"(a,F9.5,1x,a,I11,1x,F9.5)","Elapsed Time: ",time, "Loops left: ",loops_left,getPwr()
        end if

        !time loop
        do p = 1, loops
            !heat sim loops
            do k = zi, zf
                do j = yi, yf
                    do i = xi, xf
                        !!adapted from Finite Differecne Methods in Heat Transfer Chapter 9. M. Ozisik
                        kappaPlusHalf = .5d0 * (kappa(i,j,k) + kappa(i,j,k+1))
                        kappaMinHalf = .5d0 * (kappa(i,j,k) + kappa(i,j,k-1))

                        densityPlusHalf = .5d0 * (density(i,j,k) + density(i,j,k+1))
                        densityMinHalf = .5d0 * (density(i,j,k) + density(i,j,k-1))

                        heatcapPlusHalf = .5d0 * (heatcap(i,j,k) + heatcap(i,j,k+1))
                        heatcapMinHalf = .5d0 * (heatcap(i,j,k) + heatcap(i,j,k-1))

                        a = 0.5d0 * (kappaMinHalf/(densityMinHalf * heatcapMinHalf)) * (1.d0/dz**2)
                        d = 0.5d0 * (kappaPlusHalf/(densityPlusHalf * heatcapPlusHalf)) * (1.d0/dz**2)
                        b = 0.5d0 * (a + d)

                        u_zz = a*t0(i,j,k-1) - 2.d0*b*t0(i,j,k) + d*t0(i,j,k+1)


                        kappaPlusHalf = 0.5d0 * (kappa(i,j,k) + kappa(i,j+1,k))
                        kappaMinHalf = 0.5d0 * (kappa(i,j,k) + kappa(i,j-1,k))

                        densityPlusHalf = 0.5d0 * (density(i,j,k) + density(i,j+1,k))
                        densityMinHalf = 0.5d0 * (density(i,j,k) + density(i,j-1,k))

                        heatcapPlusHalf = 0.5d0 * (heatcap(i,j,k) + heatcap(i,j+1,k))
                        heatcapMinHalf = 0.5d0 * (heatcap(i,j,k) + heatcap(i,j-1,k))

                        a = 0.5d0 * (kappaMinHalf/(densityMinHalf * heatcapMinHalf)) * (1.d0/dy**2)
                        d = 0.5d0 * (kappaPlusHalf/(densityPlusHalf * heatcapPlusHalf)) * (1.d0/dy**2)
                        b = 0.5d0 * (a + d)

                        u_yy = a*t0(i,j-1,k) - 2.d0*b*t0(i,j,k) + d*t0(i,j+1,k)


                        kappaPlusHalf = .5d0 * (kappa(i,j,k) + kappa(i+1,j,k))
                        kappaMinHalf = .5d0 * (kappa(i,j,k) + kappa(i-1,j,k))

                        densityPlusHalf = .5d0 * (density(i,j,k) + density(i+1,j,k))
                        densityMinHalf = .5d0 * (density(i,j,k) + density(i-1,j,k))

                        heatcapPlusHalf = .5d0 * (heatcap(i,j,k) + heatcap(i+1,j,k))
                        heatcapMinHalf = .5d0 * (heatcap(i,j,k) + heatcap(i-1,j,k))

                        a = 0.5d0 * (kappaMinHalf/(densityMinHalf * heatcapMinHalf)) * (1.d0/dx**2)
                        d = 0.5d0 * (kappaPlusHalf/(densityPlusHalf * heatcapPlusHalf)) * (1.d0/dx**2)
                        b = 0.5d0 * (a + d)

                        u_xx = a*t0(i-1,j,k) - 2.d0*b*t0(i,j,k) + d*t0(i+1,j,k)

                        tempIncrease = delt * (u_xx + u_yy + u_zz) 
                        energyIncrease = laserOn*jtmp(i,j,k)*delt*volumeVoxel + heatcap(i,j,k)*massVoxel*tempIncrease

                        tn(i,j,k) =  tn(i,j,k) + tempIncrease + laserOn*coeff(i,j,k)*jtmp(i,j,k)
                        !check result is physical
                        if(temp(i,j,k) < 0.d0)then
                            print*,id,i,j,k,t0(i,j,k),loops_left
                            call mpi_abort(new_comm, 1)
                        end if
                    end do
                end do
            end do
            t0 = tn
            !halo swap
            !send_recv data to right
                call MPI_Sendrecv(t0(:,:,zf), size(t0(:,:,zf)), mpi_double_precision, right, tag, &
                              t0(:,:,zf+1), size(t0(:,:,zf+1)), mpi_double_precision, right, tag, &
                              new_comm, recv_status)
            
            !send_recv data to left
                call MPI_Sendrecv(t0(:,:,zi), size(t0(:,:,zi)), mpi_double_precision, left, tag, &
                              t0(:,:,zi-1), size(t0(:,:,zi-1)), mpi_double_precision, left, tag, &
                              new_comm, recv_status)

            if(pulseCount >= realpulseLength .and. laser_flag)then!turn laser off
                laser_flag = .false.
                laserOn = 0.
                pulseCount = 0.
                pulsesDone = pulsesDone + 1
                repetitionCount = 0.
            elseif((repetitionCount >= repetitionRate_1) .and. (.not.laser_flag) .and. (pulsesDone < pulsesToDo))then
                laser_flag = .true.
                laserOn = 1.
                pulseCount = 0.
                repetitionCount = 0.
            end if

            pulseCount = pulseCount + delt
            repetitionCount = repetitionCount + delt
            time = time + delt
        end do

        !collate results
        call mpi_allgather(t0(:,:,zi:zf), size(t0(:,:,zi:zf)), mpi_double_precision, &
                        temp(:,:,zi:zf), size(temp(:,:,zi:zf)), mpi_double_precision, new_comm)

        deallocate(T0)
        deallocate(Tn)
        deallocate(jtmp)

    end subroutine heat_sim_3D


    subroutine initThermalCoeff(delt, numpoints, xmax, ymax, zmax, numproc)
    !init all thermal variables

        use thermalConstants
        use constants, only : nxg, nyg, nzg, spotsPerRow, spotsPerCol, pulsetype
        use memoryModule, only : checkallocate

        implicit none

        real,    intent(INOUT) :: delt
        real,    intent(IN)    :: xmax, ymax, zmax
        integer, intent(IN)    :: numpoints, numproc

        real :: densitytmp, alphatmp, kappatmp, heatCaptmp, constd

        dx = (2.d0 * xmax * 1.d-2) / (numpoints + 2.d0)
        dy = (2.d0 * ymax * 1.d-2) / (numpoints + 2.d0)
        dz = (2.d0 * zmax * 1.d-2) / (numpoints + 2.d0)  

        call checkallocate(coeff,   [nxg+1, nyg+1, nzg+1], "coeff",   numproc, [0,0,0])
        call checkallocate(alpha,   [nxg+1, nyg+1, nzg+1], "alpha",   numproc, [0,0,0])
        call checkallocate(kappa,   [nxg+1, nyg+1, nzg+1], "kappa",   numproc, [0,0,0])
        call checkallocate(density, [nxg+1, nyg+1, nzg+1], "density", numproc, [0,0,0])
        call checkallocate(heatCap, [nxg+1, nyg+1, nzg+1], "heatcap", numproc, [0,0,0])

        call checkallocate(watercontent, [nxg, nyg, nzg], "watercontent", numproc)

        skinDensityInit = getSkinDensity(watercontentInit)
        WaterContent = watercontentInit
        heatCaptmp =  getSkinHeatCap(watercontentInit)
        densitytmp  = getSkinDensity(watercontentInit)
        kappatmp = getSkinThermalCond(watercontentInit, densitytmp)
        alphatmp = kappatmp / (densitytmp * getSkinHeatCap(watercontentInit))

        alpha = alphatmp
        alpha(:,:,nzg+1) = airThermalCond(25.d0 + 273.d0, 0) / (airDensity(25.d0 + 273.d0) * airHeatCap) 

        kappa = airThermalCond(25.+273., 0)
        kappa(1:nxg,1:nyg,1:nzg) = getSkinThermalCond(watercontentInit, densitytmp)

        density = densitytmp
        heatcap = heatcaptmp

        constd = (1.d0/dx**2) + (1.d0/dy**2) + (1.d0/dz**2)
        delt   = 1d-5!1.d0 / (1.d0*alphatmp*constd)

        coeff = 0.d0
        coeff(1:nxg,1:nyg,1:nzg) = alphatmp * delt / kappatmp

        pulseLength = (energyPerPixel * 1.d-3 * real(spotsPerRow* spotsPerCol)) / Power!pulselength above avg pwr

        volumeVoxel = (2.d0*xmax*1.d-2/nxg) * (2.d0*ymax*1.d-2/nyg) * (2.d0*zmax*1.d-2/nzg)
        massVoxel = densitytmp*volumeVoxel
        QVapor = lw * massVoxel

        select case(trim(pulsetype))
        case ("tophat")
            getPwr => getPwrTopHat
            realPulseLength = pulseLength !total pulse length
        case ("gaussian")
            getPwr => getPwrGaussian
            realPulseLength = 20000.d0 * pulseLength !total pulse length
        case("triangular")
            getPwr => getPwrTriangular
            realPulseLength = 2.d0 * pulseLength !total pulse length
        case default
            call mpi_finalize()
            error stop "no pulse type"
        end select

    end subroutine initThermalCoeff

    !power functions for laser pulse
    real function getPwrGaussian() result (getPwr)

        implicit none

        real :: mu, sig, fact

        fact = (2. * sqrt(2. * log(2.)))
        mu = fact * pulseLength
        sig = pulseLength / fact ! convert FWHM to sigma 

        getPwr = power * exp(-(time - mu)**2 / (2.d0*sig**2))

    end function getPwrGaussian


    real function getPwrTopHat() result (getPwr)

        implicit none

        if(.not. laser_flag)then
            getPwr = 0.d0
            return
        else
            getPwr = power
            return
        end if

    end function getPwrTopHat


    real function getPwrTriangular() result (getPwr)

        implicit none

        real :: m, c

        m = power / pulseLength
        c = 2.d0 * power

        if(.not. laser_flag)then
            getPwr = 0.d0
        else
            if(pulseFlag)then
                getPwr = -m * time + c
                if(getPwr < 0.d0)getPwr=0.d0
            else
                if(time >= pulseLength)then
                    pulseFlag=.true.
                    getPwr = -m * time + c
                    if(getPwr < 0.d0)getPwr=0.d0
                else
                    getPwr = m * time
                end if
            end if
        end if

    end function getPwrTriangular

end Module Heat