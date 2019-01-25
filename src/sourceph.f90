module sourceph_mod

    implicit none
    save

    contains
    
        subroutine sourcephCO2(xmax,ymax,zmax,xcell,ycell,zcell,iseed)
        ! get phton entry location for C02 laser

            use constants, only : nxg,nyg,nzg,twopi,spotsPerRow, spotsPerCol
            use photon_vars

            implicit none


            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed
            real,    intent(IN)    :: xmax, ymax, zmax
            real                   :: ran2

            real :: spotSize, gridWidth, gridHeight, spotGapCol, spotGapRow, x, y, adjustFactorRow, adjustFactorCol, ranx, rany
            real ::  cenx, ceny, theta, r
            integer :: spotsTotal

            gridWidth   = 2. * xmax!cm
            gridHeight   = 2. * ymax!cm
            spotsTotal  = spotsPerRow * spotsPerCol
            spotSize    = 250d-4 !um
            spotGapRow  = (gridWidth - (spotsPerRow * spotSize/2.)) / (spotsPerRow)
            spotGapCol  = (gridHeight - (spotsPercol * spotSize/2.)) / (spotsPerCol) 

            !centre pixels
            adjustFactorRow = abs((nint(spotsPerRow/2.) * spotGapRow) - xmax)
            adjustFactorCol = abs((nint(spotsPerCol/2.) * spotGapcol) - ymax)

            !get random integer between 1 and spotsPer
            ranx = 4!int(ran2(iseed) * real(spotsPerRow)) + 1
            rany = 4!int(ran2(iseed) * real(spotsPerCol)) + 1

            !get centre of spot
            x = (ranx * spotGapRow) + (spotSize/2.)   
            y = (rany * spotGapCol) + (spotSize/2.)
            cenx = x - xmax + adjustFactorRow
            ceny = y - ymax + adjustFactorCol

            !http://mathworld.wolfram.com/DiskPointPicking.html
            !sample circle uniformly
            !sample radius between [0,r^2]
            ! r = ran2(iseed) * (spotSize/2.)**2
            ! theta = ran2(iseed) * twopi
            ! x = sqrt(r) * cos(theta)
            ! y = sqrt(r) * sin(theta)
            x = rang(0.d0, spotsize/4.d0, iseed)
            y = rang(0.d0, spotsize/4.d0, iseed)

            ! x = ranu(-spotSize/2.d0, spotSize/2.d0, iseed)
            ! y = ranu(-spotSize/2.d0, spotSize/2.d0, iseed)


            !get final point
            xp = x !- cenx
            yp = y !- ceny
            zp = zmax-(1.e-8*(2.*zmax/nzg))

            phi = twopi * ran2(iseed)
            cosp = cos(phi)
            sinp = sin(phi)          
            sint = 0.
            cost = -1.

            nxp = sint * cosp  
            nyp = sint * sinp
            nzp = cost

            !*************** Linear Grid *************************
            xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
            ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
            zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
            !*****************************************************
        end subroutine sourcephCO2


    subroutine sourcephERYAG(xmax, ymax, zmax, xcell, ycell, zcell, iseed)
    ! get phton entry location for ER:YAG laser

        use constants, only : nxg, nyg, nzg, twopi, spotsPerRow, spotsPerCol
        use photon_vars

        implicit none


        real,    intent(IN)    :: xmax, ymax, zmax
        integer, intent(OUT)   :: xcell, ycell, zcell
        integer, intent(INOUT) :: iseed

        real    :: width, height, ran2, spotSize, spotGapRow, spotGapCol, adjust, xPos, yPos
        real    :: r, theta, cenx, ceny, x, y
        integer :: i, j

        width  = 10d-3
        height = 10d-3
        iseed = -32489843

        spotSize = 250d-6
        spotGapRow = width / spotsPerRow
        spotGapCol = height / spotsPerCol
        adjust = spotGapRow / 2.


        do
            i = int(ran2(iseed) * real(spotsPerRow)) + 1
            j = int(ran2(iseed) * real(spotsPerCol)) + 1
            if((i == 1 .or. i == 5) .and. (j == 1 .or. j == 5))cycle 
            exit
        end do 
        xPos = (i - 1) * spotGapRow + (spotSize/2.)
        yPos = (j - 1) * spotGapcol + (spotSize/2.)
        cenx = xPos - xmax + adjust
        ceny = yPos - ymax + adjust

        r = ran2(iseed) * (spotSize/2.)**2
        theta = ran2(iseed) * twopi
        x = sqrt(r) * cos(theta)
        y = sqrt(r) * sin(theta)
        xp = x - cenx
        yp = y - ceny
        zp = zmax-(1.e-8*(2.*zmax/nzg))

        phi = twopi * ran2(iseed)
        cosp = cos(phi)
        sinp = sin(phi)          
        sint = 0.
        cost = -1.

        nxp = sint * cosp  
        nyp = sint * sinp
        nzp = cost

        !*************** Linear Grid *************************
        xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
        ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
        zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
        !*****************************************************
    end subroutine sourcephERYAG


    real function ranu(a, b, iseed)
    ! return on call a random number and updated iseed   
    ! random number is uniformly distributed between a and b
    !  INPUT:
    !        a       real     lower limit of boxcar
    !        b       real     upper limit of boxcar
    !        iseed   integer  seed integer used fot the random number generator
    !  OUTPUT:
    !        ranu    real     uniform random number
    !        iseed   integer  seed used for next call

        implicit none

        real,    intent(IN)    :: a, b
        integer, intent(INOUT) :: iseed
        real :: ran2

        ranu = a + ran2(iseed) * (b - a)
    end function ranu

    real function rang(avg, sigma, iseed)
    ! return on call a random number and updated iseed   
    ! random number is from a gaussian distrbution
    ! used the Marsaglia polar method
    !  INPUT:
    !        avg       real     mean of gaussian dist.
    !        sigma     real     var of gaussian dist
    !        iseed     integer  seed integer used fot the random number generator
    !  OUTPUT:
    !        rang    real     gaussian distributed random number
    !        iseed   integer  seed used for next call

        implicit none

        real,    intent(IN)    :: avg, sigma
        integer, intent(INOUT) :: iseed
        real :: u, s, tmp

        s = 1.d0
        do while(s.ge.1.)
            u = ranu(-1.,1.,iseed)
            s = ranu(-1.,1.,iseed)
            s = s**2. + u**2.
        end do
        tmp = u*sqrt(-2.*log(s)/s)

        rang = avg + sigma*tmp

    end function rang

end MODULE sourceph_mod
