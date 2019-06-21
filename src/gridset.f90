module gridset_mod

   implicit none

   contains
      subroutine gridset(xmax,ymax,zmax,id)

         use constants, only : nxg,nyg,nzg
         use iarray, only    : rhokap, xface, yface, zface, rhokap, albedo
         use opt_prop, only  : kappa
         use ch_opt

         implicit none

         integer :: i, j, k, id
         real    :: xmax, ymax, zmax, x, y, z, mus, mua, g, n

         if(id == 0)then
            print*, ' '
            print *, 'Setting up density grid....'
         end if
         !**********  Linear Cartesian grid. Set up grid faces ****************
         do i = 1, nxg+1
            xface(i) = (i-1) * 2. * xmax/nxg
         end do
         do i = 1, nyg+1
            yface(i) = (i-1) * 2. * ymax/nyg
         end do
         do i = 1, nzg+1
            zface(i) = (i-1) * 2. * zmax/nzg
         end do
         call init_opt1
         !**************  Loop through x, y, and z to set up grid density and refractive index grid.  ****
         rhokap = 0.
         do i = 1, nxg
               x = xface(i) - xmax + xmax/nxg
            do j = 1, nyg
                  y = yface(j) - ymax + ymax/nyg
               do k = nzg, 1, -1
                  z = zface(k) + zmax/nzg
                  !set density 
                  if(0.1 >= 2.*zmax - z)then
                     !epi
                     mua = 0.78!5.42d0
                     g = 0.87!.75d0
                     mus = 30.9 / (1. - g)!64.3 / (1. - g)
                     n = 1.45d0
                     rhokap(k) = mus + mua
                     albedo(k) = mus / rhokap(k)
                  elseif(0.2 >= 2.*zmax - z)then
                     ! print*,"pap"
                     mua = .66!3.55d0
                     g = .72!.71d0
                     mus = 16.6 / (1. - g)!39.7d0 / (1. - g)
                     n = 1.39d0
                     rhokap(k) = mus + mua
                     albedo(k) = mus / rhokap(k)
                  elseif(0.7 >= 2.*zmax - z)then
                     ! print*,"ret"
                     mua = 0.64!2.90d0
                     g = .72!.71d0
                     mus = 16.6!39.7d0 / (1. - g)
                     n = 1.39
                     rhokap(k) = mus + mua
                     albedo(k) = mus / rhokap(k)
                  else
                     ! print*,"hypo"
                     mua = .18!0.84d0
                     g = .78!.78d0
                     mus = 5.3/ (1. - g)!9.2d0 / (1. - g)
                     n = 1.37d0
                     rhokap(k) = mus + mua
                     albedo(k) = mus / rhokap(k)
                  end if
               end do
            end do
         end do
      end subroutine gridset
end module gridset_mod
