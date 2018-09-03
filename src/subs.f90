module subs

implicit none

    contains
        subroutine directory
        !  subroutine defines vars to hold paths to various folders   
        !   
        !   
            use constants, only : cwd, homedir, fileplace, resdir

            implicit none

            !get current working directory

            call get_environment_variable('PWD', cwd)

            ! get 'home' dir from cwd
            homedir = trim(cwd(1:len(trim(cwd))-3))
            ! get data dir
            fileplace = trim(homedir)//'data/'
            ! get res dir
            resdir=trim(homedir)//'res/'

        end subroutine directory

        subroutine zarray
        !   sets all arrays to zero
        !
        !
            use iarray

            implicit none

            jmean = 0.
            xface = 0.
            yface = 0.
            zface = 0.
            rhokap = 0.
            jmeanGLOBAL = 0.
        end subroutine zarray


        subroutine alloc_array(numproc)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iarray
            use constants,       only : nxg,nyg,nzg
            use iso_fortran_env, only : int64
            use utils,           only : mem_free

            implicit none

            integer , intent(IN) :: numproc

            integer(int64) :: limit, cnt, i

            limit = mem_free()
            cnt = 0_int64

            allocate(xface(nxg+1))
            inquire(iolength=i)xface(:)
            call chck_mem(cnt, i, limit, 'xface', numproc)

            allocate(yface(nyg+1))
            inquire(iolength=i)yface(:)
            call chck_mem(cnt, i, limit, 'yface', numproc)

            allocate(zface(nzg+1))
            inquire(iolength=i)zface(:)
            call chck_mem(cnt, i, limit, 'zface', numproc)

            allocate(rhokap(0:nxg+1,0:nyg+1,0:nzg+1))
            inquire(iolength=i)rhokap(:,:,:)
            call chck_mem(cnt, i, limit, 'rhokap', numproc)

            allocate(jmean(nxg, nyg, nzg))
            inquire(iolength=i)jmean(:,:,:)
            call chck_mem(cnt, i, limit, 'jmean', numproc)

            allocate(jmeanGLOBAL(nxg, nyg, nzg))
            inquire(iolength=i)jmeanGLOBAL(:,:,:)
            call chck_mem(cnt, i, limit, 'jmeanGLOBAL', numproc)
        end subroutine alloc_array


        subroutine chck_mem(cur, new, limit, name, numproc)
        !   routine to check if the system has enough RAM available in order to run the simulation
        !   cur: current memory assigned, new: new memory to be assigned
        !   limit: the limit of RAM available, name: name of array to be assigned, numproc: processor #

            use iso_fortran_env, only : int64
            use utils,           only : str

            implicit none

            integer(int64), intent(IN)    :: new, limit
            integer(int64), intent(INOUT) :: cur 
            integer,        intent(IN)    :: numproc
            character(*),   intent(IN)    :: name

            integer :: error

            cur = cur + new * numproc
            if(cur > limit)then
                print*,'Need '//str(cur-limit)//' more memory to run. '//name
                call mpi_finalize(error)
                stop
            end if
        end subroutine chck_mem
end module subs
