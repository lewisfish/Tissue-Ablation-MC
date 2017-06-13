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
   case(0:10000)!bottom left
      offy =  .5
      offx = -.5
   case(10001:20000)!mid left
      offy = 0.
      offx = -.5
   case(20001:30000)!top left
      offy = -.5
      offx = -.5
   case(30001:40000)!bottom mid
      offy = .5
      offx = 0.
   case(40001:50000)!middle
      offy = 0.
      offx = 0.
   case(50001:60000)!top mid
      offy = -.5
      offx = 0.
   case(60001:70000)!bottom right
      offy = .5
      offx = .5
   case(70001:80000)!mid right
      offy = 0.
      offx = .5
   case(80001:90000)!top right
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
