MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed)

   use constants, only : nxg,nyg,nzg,twopi
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2

   real :: spotSize, gridWidth, gridHeight, spotGapCol, spotGapRow, x, y, adjustFactorRow, adjustFactorCol, ranx, rany
   real ::  cenx, ceny, theta, r
   integer :: spotsPerRow, spotsPerCol, spotsTotal

   gridWidth   = 2. * xmax!cm
   gridHeight   = 2. * ymax!cm
   spotsPerRow = 7
   spotsPerCol = 7
   spotsTotal  = spotsPerRow * spotsPerCol
   spotSize    = 250d-4 !um
   spotGapRow  = (gridWidth - (spotsPerRow * spotSize)) / (spotsPerRow + 1.)
   spotGapCol  = (gridHeight - (spotsPercol * spotSize)) / (spotsPerCol + 1.) 

   !centre pixels
   adjustFactorRow = abs((nint(7./2.) * spotGapRow) + (spotSize/2.) - xmax)
   adjustFactorCol = abs((nint(7./2.) * spotGapcol) + (spotSize/2.) - ymax)

   !get random integer between 1 and spotsPer
   ranx = int(ran2(iseed) * real(spotsPerRow)) + 1
   rany = int(ran2(iseed) * real(spotsPerCol)) + 1

   !get centre of spot
   x = (ranx * spotGapRow) + (spotSize/2.)   
   y = (rany * spotGapCol) + (spotSize/2.)
   cenx = x - xmax + adjustFactorRow
   ceny = y - ymax + adjustFactorCol

   !http://mathworld.wolfram.com/DiskPointPicking.html
   !sample circle uniformly
   !sample radius between [0,r^2]
   r = ran2(iseed) * (spotSize/2.)**2
   theta = ran2(iseed) * twopi
   x = sqrt(r) * cos(theta)
   y = sqrt(r) * sin(theta)
   !get final point
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
   end subroutine sourceph
end MODULE sourceph_mod
