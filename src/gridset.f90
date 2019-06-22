module gridset_mod

   implicit none

   contains
      subroutine gridset(xmax,ymax,zmax,id)

         use constants, only : nxg,nyg,nzg
         use iarray,    only : rhokap, xface, yface, zface, rhokap, albedo, g2, refrac, hgg, muas
         use opt_prop,  only : hgg_epi, hgg_pap, hgg_ret, hgg_hypo, mua_epi, mua_pap, mua_ret, mua_hypo, mus_epi, mus_pap, mus_ret,&
                               mus_hypo, n_epi, n_pap, n_ret, n_hypo
         use ch_opt

         implicit none

         integer :: i, j, k, id
         real    :: xmax, ymax, zmax, x, y, z

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
                     rhokap(k) = mus_epi + mua_epi
                     muas(k) = mua_epi
                     albedo(k) = mus_epi / rhokap(k)
                     refrac(k) = n_epi
                     hgg(k) = hgg_epi
                     g2(k) = hgg_epi**2
                  elseif(0.2 >= 2.*zmax - z)then
                     ! print*,"pap"
                     rhokap(k) = mus_pap + mua_pap
                     muas(k) = mua_pap
                     albedo(k) = mus_pap / rhokap(k)
                     refrac(k) = n_pap
                     hgg(k) = hgg_pap
                     g2(k) = hgg_pap**2
                  elseif(0.7 >= 2.*zmax - z)then
                     ! print*,"ret"
                     rhokap(k) = mus_ret + mua_ret
                     muas(k) = mua_ret
                     albedo(k) = mus_ret / rhokap(k)
                     refrac(k) = n_ret
                     hgg(k) = hgg_ret
                     g2(k) = hgg_ret**2
                  else
                     ! print*,"hypo"
                     rhokap(k) = mus_hypo + mua_hypo
                     muas(k) =   mua_hypo
                     albedo(k) = mus_hypo / rhokap(k)
                     refrac(k) = n_hypo
                     hgg(k) = hgg_hypo
                     g2(k) = hgg_hypo**2
                  end if
               end do
            end do
         end do
      end subroutine gridset
end module gridset_mod
