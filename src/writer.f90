module writer_mod

   implicit none

   contains
      subroutine writer(xmax, ymax, zmax, nphotons, numproc, ablateTemp, temp)
      !   write out arrays
      !
      !
         use constants, only : nxg, nyg, nzg, fileplace
         use iarray,    only : jmeanGLOBAL, rhokap
         use heat,      only : watercontent, energyPerPixel
         use utils,     only : str

         implicit none

         integer, intent(IN) :: nphotons, numproc
         real,    intent(IN) :: xmax, ymax, zmax, ablateTemp, temp(0:nxg+1,0:nyg+1,0:nzg+1) 

         integer :: u


         open(newunit=u,file=trim(fileplace)//"ErYAG/jmean-timestep-test2-"&
           //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)jmeanGLOBAL
         close(u)

         open(newunit=u,file=trim(fileplace)//"ErYAG/rhokap-timestep-test2-"&
           //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)rhokap(1:nxg, 1:nyg, 1:nzg)
         close(u)

         open(newunit=u,file=trim(fileplace)//"ErYAG/temp-timestep-test2-"&
         //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)temp - 273.d0
         close(u)

         open(newunit=u,file=trim(fileplace)//"ErYAG/water-timestep-test2-"&
         //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)watercontent
         close(u)
      end subroutine writer
end module writer_mod
