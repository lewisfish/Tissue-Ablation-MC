module writer_mod

   implicit none

   contains
      subroutine writer(temp, xmax, ymax, zmax)
      !   write out arrays
      !
      !
         use constants, only : nxg, nyg, nzg, fileplace
         use iarray,    only : jmeanGLOBAL, rhokap
         use heat,      only : energyPerPixel, power
         use utils,     only : str

         implicit none

         real, intent(IN) :: temp(0:nxg+1,0:nyg+1,0:nzg+1)
         real, intent(IN) :: xmax, ymax, zmax

         integer :: u


         open(newunit=u,file=trim(fileplace)//"letterjmean-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)jmeanGLOBAL
         close(u)


         open(newunit=u,file=trim(fileplace)//"lettertemp-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)temp - 273.d0
         close(u)

      end subroutine writer
end module writer_mod
