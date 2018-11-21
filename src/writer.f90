module writer_mod

   implicit none

   contains
      subroutine writer(ablateTemp, temp, tissueGlobal, xmax, ymax, zmax)
      !   write out arrays
      !
      !
         use constants, only : nxg, nyg, nzg
         use iarray,    only : jmeanGLOBAL, rhokap
         use heat,      only : energyPerPixel, power
         use utils,     only : str

         implicit none

         real,    intent(IN) :: ablateTemp, temp(0:nxg+1,0:nyg+1,0:nzg+1), tissueGlobal(:,:,:)!, time(:,:,:,:)
         real, intent(IN) :: xmax, ymax, zmax

         integer :: u

     !     open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/jmean-"//str(int(power))//"w-"&
     ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//"-"//str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
     !      ,access="stream",form="unformatted", status="replace")
     !     write(u)jmeanGLOBAL
     !     close(u)

         open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/rhokap4-"//str(int(power))//"w-"&
//str(nzg)//"-"//str(ablateTemp,3)//"-"//str(int(energyPerPixel),3)//"-"//str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
          ,access="stream",form="unformatted", status="replace")
         write(u)rhokap(1:nxg, 1:nyg, 1:nzg)
         close(u)

      !    open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/temp-"//str(int(power))//"w-"&
      ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//"-"//str(xmax,5)//"-"//str(ymax,5)//"-"//str(zmax,5)//".dat" &
      !     ,access="stream",form="unformatted", status="replace")
      !    write(u)temp - 273.d0
      !    close(u)

         ! open(newunit=u,file="/home/lewis/phdshizz/ablation/newdata/water-"&
         ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)watercontent
         ! close(u)

         ! open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/tissue-"//str(int(power))//"w-"&
         ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)tissueGlobal
         ! close(u)

         ! open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/time-1-"//str(int(power))//"w-"&
         ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,1)
         ! close(u)


         ! open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/time-2-"//str(int(power))//"w-"&
         ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,2)
         ! close(u)


         ! open(newunit=u,file="/home/lewis/phdshizz/ablation/testrun/time-3-"//str(int(power))//"w-"&
         ! //str(nzg)//"-"//str(ablateTemp,3)//"-"//str(energyPerPixel,3)//".dat" &
         !  ,access="stream",form="unformatted", status="replace")
         ! write(u)time(:,:,:,3)
         ! close(u)
      end subroutine writer
end module writer_mod
