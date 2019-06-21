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
            albedo = 0.
            jmeanGLOBAL = 0.
        end subroutine zarray


        subroutine alloc_array(numproc, id)
        !  subroutine allocates allocatable arrays
        !   
        !   
            use iarray
            use constants,       only : nxg, nyg, nzg
            use iso_fortran_env, only : int64
            use memorymodule

            implicit none

            integer , intent(IN) :: numproc, id

            integer(int64) :: limit
            integer :: N

            limit = mem_free()
            ! cnt = 0_int64


            call checkallocate(xface, [nxg+1], "xface", numproc)
            call checkallocate(yface, [nyg+1], "yface", numproc)
            call checkallocate(zface, [nzg+1], "zface", numproc)

            call checkallocate(rhokap, [nzg+1], "rhokap", numproc)
            call checkallocate(albedo, [nzg], "albedo", numproc)
            call checkallocate(jmean, [nxg, nyg, nzg], "jmean", numproc)
            call checkallocate(jmeanGLOBAL, [nxg, nyg, nzg], "jmeanGlobal", numproc)

            call checkallocate(tissue, [nxg, nyg, nzg], "tissue", numproc)
            N = nzg ! points for heat sim
            call checkallocate(temp, [N+1, N+1, N+1], "temp", numproc, [0,0,0])

            if(id==0)print'(A,1X,F5.2,A)','allocated:',dble(TotalMem)/dble(limit)*100.d0,' % of total RAM'
            
        end subroutine alloc_array
end module subs
