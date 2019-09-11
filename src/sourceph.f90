module sourceph_mod

    implicit none

    contains
    
        subroutine sourcephCO2(xmax,ymax,zmax,xcell,ycell,zcell,iseed)
        ! get phton entry location for C02 laser

            use constants, only : nxg, nyg, nzg, twopi
            use photon_vars

            implicit none


            integer, intent(OUT)   :: xcell, ycell, zcell
            integer, intent(INOUT) :: iseed
            real,    intent(IN)    :: xmax, ymax, zmax
            real                   :: ran2

            real :: spotSize, theta, r

            spotSize    = 250d-4 !um

            !http://mathworld.wolfram.com/DiskPointPicking.html
            !sample circle uniformly
            !sample radius between [0,r^2]
            r = ran2(iseed) * (spotSize/2.)**2
            theta = ran2(iseed) * twopi
            xp = sqrt(r) * cos(theta)
            yp = sqrt(r) * sin(theta)
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
