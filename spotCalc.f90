program spotDist

    implicit none
    !all unit in cm
    real    :: spotSize, spotGap, gridWidth, spotGapCol, spotGapRow, x, y, adjustFactorRow, adjustFactorCol, ran2
    real    :: xp, yp, twopi, r, theta, cenx, ceny
    integer :: i, j, spotsPerRow, spotsPerCol, spotsTotal, iseed, ranx, rany

    twopi = 2. * 4.*atan(1.)
    iseed = -3284984

    gridWidth   = 1.1!cm
    spotsPerRow = 7
    spotsPerCol = 7
    spotsTotal  = spotsPerRow * spotsPerCol
    spotSize    = 150d-4 !um
    spotGapRow  = (gridWidth - (spotsPerRow * spotSize)) / (spotsPerRow + 1.)
    spotGapCol  = (gridWidth - (spotsPercol * spotSize)) / (spotsPerCol + 1.) 

    adjustFactorRow = abs((nint(7./2.) * spotGapRow) + (spotSize/2.) - .55)
    adjustFactorCol = abs((nint(7./2.) * spotGapcol) + (spotSize/2.) - .55)

    do i = 1, 1000000

        ranx = int(ran2(iseed) * 7.) + 1
        rany = int(ran2(iseed) * 7.) + 1

        x = (ranx * spotGapRow) + (spotSize/2.)   
        y = (rany * spotGapCol) + (spotSize/2.)
        cenx = x - .55 + adjustFactorRow
        ceny = y - .55 + adjustFactorCol

        !http://mathworld.wolfram.com/DiskPointPicking.html
        !sample circle uniformly
        !sample radius between [0,r^2]
        r = ran2(iseed) * (spotSize/2.)**2
        theta = ran2(iseed) * twopi
        x = sqrt(r) * cos(theta)
        y = sqrt(r) * sin(theta)
        xp = x - cenx
        yp = y - ceny
        print*,xp,yp

    end do

end program spotDist