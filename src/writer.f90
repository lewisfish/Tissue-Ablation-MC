module writer_mod

   implicit none

   contains
      subroutine writer(ablateTemp, temp, tissueGlobal, xmax, ymax, zmax, time)
      !   write out arrays
      !
      !
         use constants, only : nxg, nyg, nzg, fileplace
         use iarray,    only : jmeanGLOBAL, rhokap
         use heat,      only : energyPerPixel, power, watercontent
         use utils,     only : str

         implicit none

         real,    intent(IN) :: ablateTemp, temp(0:nxg+1,0:nyg+1,0:nzg+1), tissueGlobal(:,:,:), time(:,:,:,:)
         real, intent(IN) :: xmax, ymax, zmax

         integer :: u

         open(newunit=u,file=trim(fileplace)//"/22-11-18/jmean-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)jmeanGLOBAL
         close(u)

         open(newunit=u,file=trim(fileplace)//"/22-11-18/rhokap-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)rhokap(1:nxg, 1:nyg, 1:nzg)
         close(u)

         open(newunit=u,file=trim(fileplace)//"/22-11-18/temp-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)temp - 273.d0
         close(u)

         open(newunit=u,file=trim(fileplace)//"/22-11-18/water-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)watercontent
         close(u)

         open(newunit=u,file=trim(fileplace)//"/22-11-18/tissue-"//str(int(power))//"w-"&
                              //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
                              //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)tissueGlobal
         close(u)

         ! open(newunit=u,file=trim(fileplace)//"/time-1-"//str(int(power))//"w-"&
         !                      //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
         !                      //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,1)
         ! close(u)


         ! open(newunit=u,file="trim(fileplace)//"time-2-"//str(int(power))//"w-"&
         !                      //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
         !                      //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,2)
         ! close(u)


         ! open(newunit=u,file=trim(fileplace)//"/time-3-"//str(int(power))//"w-"&
         !                      //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"&
         !                      //str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,2)
         ! close(u)
      end subroutine writer
end module writer_mod
