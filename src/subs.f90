MODULE subs

implicit none
save

CONTAINS
   SUBROUTINE directory
!  subroutine defines vars to hold paths to various folders   
!   
!   
   use constants,only : cwd,homedir,fileplace,resdir
   
!#ifdef intel
!   use ifport
!#endif
   
   implicit none

   integer :: io

   !get current working directory
!#ifdef intel
!   io = getcwd(cwd)
!#else
   CALL getcwd(cwd)
!#endif

   !get 'home' dir from cwd
   homedir=trim(cwd(1:len(trim(cwd))-3))
   !get data dir
   fileplace=trim(homedir)//'data/'
   !checks to see if data folder exists, if not creates it.
!#ifdef intel
!   io = chdir(fileplace)
!#else
   call chdir(fileplace,io)
!#endif
   if(io.ne.0)then
      print*,'data directory does not exist...'
      print*, 'creating directory...'
!#ifdef intel
!      io = system("mkdir "//fileplace)
!      io = chdir(fileplace)
!      io = system("mkdir jmean/")
!      io = system("mkdir deposit/")
!      io = system("mkdir im/")
!#else
      call system("mkdir "//trim(fileplace)) !dosent work with spaces :(
      call chdir(trim(fileplace),io)
      call system("mkdir jmean/")
      call system("mkdir deposit/")
      call system("mkdir im/")
      print*, 'created directory ',trim(fileplace)
!#endif
   end if
   
   !get res dir
   resdir=trim(homedir)//'res/'
   
   end SUBROUTINE directory
   
   SUBROUTINE zarray
   
   use iarray
   
   !sets all arrays to zero
   implicit none
   
   
   jmean=0.
   xface=0.
   yface=0.
   zface=0.
   rhokap=0.
   jmeanGLOBAL = 0.
   
   end SUBROUTINE zarray

   SUBROUTINE alloc_array
!  subroutine allocates allocatable arrays
!   
!   
   use iarray
   use constants,only : nxg,nyg,nzg
   
   implicit none
   
   allocate(xface(nxg+1),yface(nyg+1),zface(nzg+1))
   allocate(rhokap(nxg,nyg,nzg))
   allocate(jmean(nxg,nyg,nzg), jmeanGLOBAL(nxg,nyg,nzg))
   
   
   end SUBROUTINE alloc_array
end MODULE subs
