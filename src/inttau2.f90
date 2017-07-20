module inttau2

   implicit none
   
CONTAINS

   subroutine tauint1(xmax,ymax,zmax,xcell,ycell,zcell,tflag,iseed,delta)
   !optical depth integration subroutine
   !
   !
      use photon_vars, only : xp, yp, zp
      use iarray,      only : jmean, rhokap
      use opt_prop,    only : mua

      
   
      implicit none
      
      real,    intent(IN)    :: xmax, ymax, zmax, delta
      integer, intent(INOUT) :: xcell, ycell, zcell, iseed
      logical, intent(INOUT) :: tflag
      
      real                   :: tau, taurun, taucell, xcur, ycur, zcur, d, dcell, ran2
      integer                :: celli, cellj, cellk
      logical                :: dir(3)
      
      xcur = xp + xmax
      ycur = yp + ymax
      zcur = zp + zmax
      
      celli = xcell
      cellj = ycell
      cellk = zcell

      taurun = 0.
      d = 0.
      dir = (/.FALSE., .FALSE., .FALSE./)
      
      tau = -log(ran2(iseed))

      do
         dir = (/.FALSE., .FALSE., .FALSE./)
         dcell = wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
         taucell = dcell * rhokap(celli,cellj,cellk)

         if(taurun + taucell < tau)then
            taurun = taurun + taucell
            d = d + dcell
            jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell*mua

            call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .TRUE., dir, delta)
         
         else

            dcell = (tau - taurun) / rhokap(celli,cellj,cellk)
            d = d + dcell
            jmean(celli, cellj, cellk) = jmean(celli, cellj, cellk) + dcell*mua
            call update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, .FALSE., dir, delta)
            exit
         end if
         
         if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            tflag = .true.
            exit
         end if
         
      end do
   
      xp = xcur - xmax
      yp = ycur - ymax
      zp = zcur - zmax
      xcell = celli
      ycell = cellj
      zcell = cellk

   end subroutine tauint1
   
   
   real function wall_dist(celli, cellj, cellk, xcur, ycur, zcur, dir)
   !funtion that returns distant to nearest wall and which wall that is (x,y or z)
   !
   !
      use iarray,      only : xface, yface, zface
      use photon_vars, only : nxp, nyp, nzp
   
      implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur
      logical, intent(INOUT) :: dir(:)
      integer, intent(INOUT) :: celli, cellj, cellk
      real                   :: dx, dy, dz
   
   
      if(nxp > 0.)then
         dx = (xface(celli+1) - xcur)/nxp
      elseif(nxp < 0.)then
         dx = (xface(celli) - xcur)/nxp
      elseif(nxp == 0.)then
         dx = 100000.
      end if
      
      if(nyp > 0.)then
         dy = (yface(cellj+1) - ycur)/nyp
      elseif(nyp < 0.)then
         dy = (yface(cellj) - ycur)/nyp
      elseif(nyp == 0.)then
         dy = 100000.
      end if
      
      if(nzp > 0.)then
         dz = (zface(cellk+1) - zcur)/nzp
      elseif(nzp < 0.)then
         dz = (zface(cellk) - zcur)/nzp
      elseif(nzp == 0.)then
         dz = 100000.
      end if
      
      wall_dist = min(dx, dy, dz)
      if(wall_dist < 0.)print'(A,7F9.5)','dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
      if(wall_dist == dx)dir=(/.TRUE., .FALSE., .FALSE./)
      if(wall_dist == dy)dir=(/.FALSE., .TRUE., .FALSE./)
      if(wall_dist == dz)dir=(/.FALSE., .FALSE., .TRUE./)
      if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
      
   end function wall_dist


   subroutine update_pos(xcur, ycur, zcur, celli, cellj, cellk, dcell, wall_flag, dir, delta)
   !routine that upates postions of photon and calls fresnel routines if photon leaves current voxel
   !
   !
      use photon_vars, only : nxp, nyp, nzp, nxp, nyp, nzp
      use iarray,      only : xface, yface, zface
   
      implicit none
      
      real,    intent(INOUT) :: xcur, ycur, zcur
      real,    intent(IN)    :: dcell, delta
      integer, intent(INOUT) :: celli, cellj, cellk
      logical, intent(IN)    :: wall_flag, dir(:)   
      
      if(wall_flag)then
      
         if(dir(1))then
            if(nxp > 0.)then
               xcur = xface(celli+1) + delta
            elseif(nxp < 0.)then
               xcur = xface(celli) - delta
            else
               print*,'Error in x dir in update_pos', dir, nxp, nyp, nzp
            end if
            ycur = ycur + nyp*dcell 
            zcur = zcur + nzp*dcell
         elseif(dir(2))then
            xcur = xcur + nxp*dcell
            if(nyp > 0.)then
               ycur = yface(cellj+1) + delta
            elseif(nyp < 0.)then
               ycur = yface(cellj) - delta
            else
               print*,'Error in y dir in update_pos', dir, nxp, nyp, nzp
            end if
            zcur = zcur + nzp*dcell
         elseif(dir(3))then
            xcur = xcur + nxp*dcell
            ycur = ycur + nyp*dcell 
            if(nzp > 0.)then
               zcur = zface(cellk+1) + delta
            elseif(nzp < 0.)then
               zcur = zface(cellk) - delta
            else
               print*,'Error in z dir in update_pos', dir, nxp, nyp, nzp
            end if
         else
            print*,'Error in update_pos...',dir
            error stop 0
         end if
      else
      
         xcur = xcur + nxp*dcell
         ycur = ycur + nyp*dcell 
         zcur = zcur + nzp*dcell
      
      end if
      
      
      if(wall_flag)then
         call update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
      end if
      
   end subroutine update_pos


   subroutine update_voxels(xcur, ycur, zcur, celli, cellj, cellk)
   !updates the current voxel based upon position
   !
   !
      use iarray, only : xface, yface, zface
   
      implicit none
      
      real,    intent(IN)    :: xcur, ycur, zcur
      integer, intent(INOUT) :: celli, cellj, cellk
         
      celli = find(xcur, xface) 
      cellj = find(ycur, yface)
      cellk = find(zcur, zface) 

   end subroutine update_voxels


   integer function find(val, a)
   !searchs for bracketing indicies for a value val in an array a
   !
   !
      implicit none
      
      real, intent(IN) :: val, a(:)
      integer          :: n, lo, mid, hi
      logical          :: ascnd

      n = size(a)
      ascnd = (a(n) >= a(1))
      lo = 0
      hi = n+1
      do
         if (hi-lo <= 1) exit
         mid = (hi+lo)/2
         if (ascnd .eqv. (val >= a(mid))) then
            lo = mid
         else
            hi = mid
         end if
      end do

      if (val == a(1)) then
         find = 1
      else if (val == a(n)) then
         find = n-1
      else if(ascnd.and. (val > a(n) .or. val < a(1))) then
         find = -1
      else if(.not.ascnd.and. (val < a(n) .or. val > a(1))) then
         find = -1
      else
         find = lo
      end if

   end function find


   subroutine repeat_bounds(cella, cellb, acur, bcur, amax, bmax, nag, nbg, delta) 
   !if photon leaves grid in a direction a or b, then photon is transported to otherside and continues being simulated
   !
   !
      implicit none
      
      real,    intent(INOUT) :: acur, bcur
      real,    intent(IN)    :: delta, amax, bmax
      integer, intent(IN)    :: nag, nbg
      integer, intent(INOUT) :: cella, cellb

      if(cella == -1)then
         if(acur < delta)then
            acur = 2.*amax  -delta
            cella = nag
         elseif(acur > 2.*amax-delta)then
            acur = delta
            cella = 1
         else
            print*,'Error in Repeat_bounds...'
            error stop 0
         end if
      end if
      if(cellb == -1)then
         if(bcur < delta)then
            bcur = 2.*bmax-delta
            cellb = nbg
         elseif(bcur > 2.*bmax-delta)then
            bcur = delta
            cellb = 1
         else
            print*,'Error in Repeat_bounds...'
            error stop 0
         end if
      endif
      
      
   end subroutine repeat_bounds
   
   
end module inttau2