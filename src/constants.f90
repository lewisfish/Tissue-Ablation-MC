MODULE constants
!
! Module containing constants:
!         PI,TWOPI, the number of grid elements in each direction n%g,
!         the various bin # params,
!         and the vars for the filepaths.
!

implicit none
save

integer, parameter :: nxg=80, nyg=80, nzg=80, spotsPerRow=9, spotsPerCol=9
real,    parameter :: PI = 3.141592, TWOPI = 6.283185, OFFSET=1.e-2*(2.*.5/nxg)
character(len=255) :: cwd, homedir, fileplace, resdir, pulsetype

end MODULE constants
