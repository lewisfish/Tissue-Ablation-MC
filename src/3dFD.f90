Module Heat

    implicit none

    real              :: pulseCount, repetitionCount, time, laserOn, total_time, repetitionRate_1, energyPerPixel
    real              :: Power=60.d0, pulselength, delt
    real              :: dx, dy, dz, eta, massVoxel, volumeVoxel
    real, allocatable :: coeff(:,:,:), kappa(:,:,:), density(:,:,:), heatcap(:,:,:), WaterContent(:,:,:), Q(:,:,:), alpha(:,:,:)
    logical           :: laser_flag
    integer           :: loops, pulsesToDo, pulsesDone, loops_left

    private
    public :: heat_sim_3D, pulseCount, repetitionCount, laserOn, laser_flag, time, initThermalCoeff, repetitionRate_1
    public :: coeff, dx, dy, dz, eta, loops, total_time, energyPerPixel, pulselength,Power, kappa,alpha
    public :: delt, massVoxel, WaterContent, Q, pulsesToDo, pulsesDone, loops_left,setupThermalCoeff

    contains

subroutine heat_sim_3D(jmean, temp, tissue, numpoints, id, numproc, new_comm, right, left, counter)

        use mpi_f08
        use utils,     only : str
        use thermalConstants, only : QVapor, getSkinHeatCap, getSkinDensity, getSkinThermalCond

        implicit none
        
        real,           intent(IN)    :: jmean(:,:,:)
        integer,        intent(IN)    :: numpoints, right, left, id, numproc,counter
        type(mpi_comm), intent(IN)    :: new_comm
        real,           intent(INOUT) :: temp(0:numpoints+1,0:numpoints+1,0:numpoints+1), tissue(:,:,:)

        type(MPI_Status) :: recv_status

        !heat variables
        real              :: u_xx, u_yy, u_zz, tempIncrease,kappaMinHalf, kappaPlusHalf
        real :: heatcapMinHalf,densityPlusHalf,densityMinHalf,heatcapPlusHalf,a,b,d
        real, allocatable :: T0(:,:,:), jtmp(:,:,:), qtmp(:,:,:)!, Ttmp(:,:,:)
        integer           :: i, j, k, p, size_x, size_y, size_z, xi, yi, zi, xf, yf, zf, N, tag,zsta,zfin

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
        ! allocate(Ttmp(numpoints, numpoints, zi:zf))
        allocate(jtmp(numpoints, numpoints, zi:zf))
        allocate(Qtmp(numpoints, numpoints, zi:zf))
        t0 = 0.
        Qtmp = 0.

        !split temp up over all processes
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

        jtmp = 0.
        call mpi_scatter(jmean, size(jmean(:,:,zi:zf)), mpi_double_precision, jtmp(:,:,zi:zf), size(jtmp(:,:,zi:zf)), &
                         mpi_double_precision, 0, new_comm)

        ! call mpi_scatter(tissue, size(tissue(:,:,zi:zf)), mpi_double_precision, Ttmp, size(Ttmp(:,:,zi:zf)), &
        !                  mpi_double_precision, 0, new_comm)


        call mpi_scatter(Q, size(Q(:,:,zi:zf)), mpi_double_precision, Qtmp(:,:,zi:zf), size(Qtmp(:,:,zi:zf)), &
                 mpi_double_precision, 0, new_comm)




        if(pulselength < delt)then
            if(id == 0)print*,"pulselength smaller than timestep, adjusting..."
            delt = pulselength / 100.
        end if
        loops_left = int(total_time/(real(loops)*delt)) - counter
        if(id == 0 .and. mod(int(total_time/(real(loops)*delt)) - counter, 100) == 0)then
        print"(a,F9.5,1x,a,I11)","Elapsed Time: ",time, "Loops left: ",loops_left!,&
            ! laser_flag,pulsesDone
        end if

        do p = 1, loops
            
            do k = zi, zf
                do j = yi, yf
                    do i = xi, xf
                        ! if((k == zf .and. id == numproc-1))then!B.Cs
                        !     u_zz = (2.d0/dz**2) * (alpha(i,j,k+1)*t0(i,j,k+1) - alpha(i,j,k)*t0(i,j,k) &
                        !            + (dz*h*alpha(i,j,k)/kappa(i,j,k)) * (tempAir - t0(i,j,k)))
                        ! else
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

                            ! u_zz = (1.d0 / dz**2) * (alpha(i,j,k-1) * t0(i,j,k-1) - 2. * alpha(i,j,k) * t0(i,j,k) + &
                            !     alpha(i,j,k+1) * t0(i,j,k+1))
                        ! end if

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





                        ! u_xx = (1.d0 / dx**2) * (alpha(i-1,j,k) * t0(i-1,j,k) - 2. *alpha(i,j,k) * t0(i,j,k) + &
                        !         alpha(i+1,j,k) * t0(i+1,j,k))

                        ! u_yy = (1.d0 / dy**2) * (alpha(i,j-1,k) * t0(i,j-1,k) - 2. *alpha(i,j,k) * t0(i,j,k) + &
                        !         alpha(i,j+1,k) * t0(i,j+1,k))

                        tempIncrease = delt * (u_xx + u_yy + u_zz) 

                        if(t0(i,j,k) >= 100. + 273. .and. Qtmp(i,j,k) < QVapor)then
                            Qtmp(i,j,k) = min(Qtmp(i,j,k) + &
                            laserOn*jtmp(i,j,k)*delt*volumeVoxel, Qvapor) !+ &
                            !getSkinHeatCap(waterContent(i,j,k))*massVoxel*tempIncrease, Qvapor)
                            t0(i,j,k) = 100.d0 + 273.d0
                        else
                            t0(i,j,k) =  t0(i,j,k) + tempIncrease + laserOn*coeff(i,j,k)*jtmp(i,j,k)
                            if(t0(i,j,k) < 0.d0)then
                                print*,id,i,j,k,t0(i,j,k),loops_left
                                call mpi_abort(new_comm, 1)
                            end if
                        end if
                    end do
                end do
            end do
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
            ! call Arrhenius(t0, time, Ttmp, zi, zf, numpoints)
        end do

        call mpi_allgather(t0(:,:,zi:zf), size(t0(:,:,zi:zf)), mpi_double_precision, &
                        temp(:,:,zi:zf), size(temp(:,:,zi:zf)), mpi_double_precision, new_comm)

        ! call mpi_allgather(Ttmp, size(Ttmp), mpi_double_precision, &
        !                 tissue(:,:,zi:zf), size(tissue(:,:,zi:zf)), mpi_double_precision, new_comm)

        call mpi_allgather(Qtmp(:,:,zi:zf), size(Qtmp(:,:,zi:zf)), mpi_double_precision, &
                        Q(:,:,zi:zf), size(Q(:,:,zi:zf)), mpi_double_precision, new_comm)

   end subroutine heat_sim_3D


    subroutine initThermalCoeff(delt, numpoints, xmax, ymax, zmax)

        use thermalConstants
        use constants, only : nxg, nyg, nzg, spotsPerRow, spotsPerCol

        implicit none

        real,    intent(INOUT) :: delt
        real,    intent(IN)    :: xmax, ymax, zmax
        integer, intent(IN)    :: numpoints

        real :: densitytmp, alphatmp, kappatmp, heatCaptmp

        dx = (2.d0 * xmax * 1.d-2) / (numpoints + 2)
        dy = (2.d0 * ymax * 1.d-2) / (numpoints + 2)
        dz = (2.d0 * zmax * 1.d-2) / (numpoints + 2)  

        skinBeta = 1.d0 + (dz*h/skinThermalCond)

        eta = eps*S_B_Constant*dx*dy !$\varepsilon \sigma A$

        allocate(coeff(0:nxg+1, 0:nyg+1, 0:nzg+1))
        allocate(alpha(0:nxg+1, 0:nyg+1, 0:nzg+1))
        allocate(kappa(0:nxg+1, 0:nyg+1, 0:nzg+1))
        allocate(density(0:nxg+1, 0:nyg+1, 0:nzg+1))
        allocate(heatCap(0:nxg+1, 0:nyg+1, 0:nzg+1))

        allocate(Q(nxg, nyg, nzg))
        allocate(watercontent(nxg, nyg, nzg))

        Q = 0.d0
        WaterContent = watercontentInit
        heatCaptmp =  getSkinHeatCap(watercontentInit)
        densitytmp  = getSkinDensity(watercontentInit)
        kappatmp = getSkinThermalCond(watercontentInit, densitytmp)
        alphatmp = kappatmp / (densitytmp * getSkinHeatCap(watercontentInit))

        alpha = alphatmp
        alpha(:,:,nzg+1) = airThermalCond(25.+273, 0) / (airDensity(25.+273.) * airHeatCap) 

        kappa = airThermalCond(25.+273., 0)
        kappa(1:nxg,1:nyg,1:nzg) = getSkinThermalCond(watercontentInit, densitytmp)

        density = densitytmp
        heatcap = heatcaptmp

        delt = (dx**2 * dz**2 ) / (2.* alphatmp * (dx**2 + 2.*dz**2))
        delt = delt / 1.

        coeff = 0.
        coeff(1:nxg,1:nyg,1:nzg) = alphatmp * delt / kappatmp

        pulseLength = (energyPerPixel * 1d-3 * real(spotsPerRow* spotsPerCol)) / Power

        volumeVoxel = (2.*xmax*1.d-2/nxg) * (2.*ymax*1.d-2/nyg) * (2.*zmax*1.d-2/nzg)
        massVoxel = densitytmp*volumeVoxel
        QVapor = lw * massVoxel
    end subroutine initThermalCoeff


    subroutine setupThermalCoeff(temp, numpoints)

        use constants, only : nxg, nyg, nzg
        use iarray,    only : rhokap
        use opt_prop,  only : mua
        use thermalConstants

        implicit none

        integer, intent(IN) :: numpoints
        real,    intent(IN) :: temp(0:numpoints+1, 0:numpoints+1, 0:numpoints+1)

        integer :: i, j, k

        WaterContent = getWaterContent(Q)

        do k = 1, nzg
            do j = 1, nyg
                do i = 1, nxg
                    if(rhokap(i,j,k) > 0.)then
                        density(i,j,k) = getSkinDensity(WaterContent(i,j,k))
                        rhokap(i,j,k) = (watercontent(i,j,k) * mua)
                        heatCap(i,j,k) = getSkinHeatCap(waterContent(i,j,k))

                        kappa(i,j,k) = getSkinThermalCond(WaterContent(i,j,k), density(i,j,k))
                        coeff(i,j,k) = delt/ (density(i,j,k) * heatcap(i,j,k))
                    elseif(rhokap(i,j,k) <= 0.)then
                        density(i,j,k) = airDensity(temp(i,j,k))
                        heatcap(i,j,k) = airThermalCond(temp(i,j,k), loops_left)

                        kappa(i,j,k) = getSkinThermalCond(WaterContent(i,j,k), density(i,j,k))
                        alpha(i,j,k) = kappa(i,j,k) / (density(i,j,k) * heatcap(i,j,k))
                        coeff(i,j,k) = delt/ (airDensity(temp(i,j,k)) * airHeatCap)
                    end if
                end do
            end do
        end do
    end subroutine setupThermalCoeff


    subroutine Arrhenius(temp, delt, tissue, zi, zf, numpoints)

        implicit none

        integer, intent(IN)    :: zi,zf, numpoints
        real,    intent(IN)    :: temp(:,:,:), delt
        real,    intent(INOUT) :: tissue(:,:,:)

        double precision :: a, g, r
        integer          :: x,y,z

        a = 3.1d91!2.9e27
        g = 6.28e5!2.4e5
        r = 8.314
        do z = zi, zf
            do y = 1, numpoints
                do x = 1, numpoints
                    if(temp(x, y, z) >= 44+273.)then
                        tissue(x, y, z) = tissue(x, y, z) + delt*A*exp(-G/(R*(temp(x,y,z))))
                    end if
                end do
            end do
        end do

    end subroutine  Arrhenius
end Module Heat