MODULE sourceph_mod

implicit none
save

CONTAINS
   subroutine sourceph(xmax,ymax,zmax,xcell,ycell,zcell,iseed,j)

   use constants, only : nxg,nyg,nzg,twopi
   use photon_vars

   implicit none


   integer, intent(OUT)   :: xcell, ycell, zcell
   integer, intent(INOUT) :: iseed
   real,    intent(IN)    :: xmax, ymax, zmax
   real                   :: ran2, radius2, offx, offy
   integer :: j


   zp = zmax-(1.e-8*(2.*zmax/nzg))

   xp = xmax*(2.*ran2(iseed)-1.)
   yp = ymax*(2.*ran2(iseed)-1.)
   
   radius2 = .05**2.

   select case(j)
   case(0:100000)!bottom left
      offy =  .5
      offx = -.5
   case(100001:200000)!mid left
      offy = 0.
      offx = -.5
   case(200001:300000)!top left
      offy = -.5
      offx = -.5
   case(300001:400000)!bottom mid
      offy = .5
      offx = 0.
   case(400001:500000)!middle
      offy = 0.
      offx = 0.
   case(500001:600000)!top mid
      offy = -.5
      offx = 0.
   case(600001:700000)!bottom right
      offy = .5
      offx = .5
   case(700001:800000)!mid right
      offy = 0.
      offx = .5
   case(800001:900000)!top right
      offy = -.5
      offx = .5      
   end select  


   do while((xp-offx)**2.+(yp-offy)**2. > radius2)
      xp = xmax*(2.*ran2(iseed)-1.)
      yp = ymax*(2.*ran2(iseed)-1.)
   end do

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
