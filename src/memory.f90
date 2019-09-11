module memoryModule
    
    use iso_fortran_env, only : int64

    implicit none

    integer(int64) :: TotalMem = 0_int64
    
    interface checkallocate
        module procedure checkallocate_1R
        module procedure checkallocate_2R
        module procedure checkallocate_3R
        module procedure checkallocate_4R
    end interface

    private
    public :: TotalMem, checkallocate, mem_free

    contains

        subroutine checkallocate_4R(array, dims, name, numproc, start)
        ! subroutine to allocate array if enough RAM is left
        ! Takes:
        !       n-dimensional arrayto be allocated         : array
        !       1D array length n, contains array dims     : dims
        !       character var, for error messgages         : name
        !       integer var, optional for calc if parallel : numproc

            implicit none

            real, allocatable, intent(INOUT) :: array(:,:,:,:)
            integer,           intent(IN)    :: dims(4)
            integer, optional, intent(IN)    :: numproc, start(4)
            character(*),      intent(IN)    :: name

            integer(int64) :: limit, procs, byte, tot, memCount

            real    :: tmp
            integer :: val, j

            !get current RAM available
            limit = mem_free()
            !set memory used counter to 0
            memCount = 0_int64

            !check if code is parallel
            if(.not.present(numproc))then
                procs = 1_int64
            else
                procs = numproc
            end if

            !get precision
            val = maxexponent(tmp)
            if(val == 1024)then
                !dbl prec
                byte = 8_int64
            elseif(val == 128)then
                !sgle prec
                byte = 4_int64
            else
                error stop
            end if

            !calc space for new array
            tot = dims(1)
            do j = 2, 4
                tot = tot * dims(j)
            end do

            ! calc memory space for array, and add to total counter
            memCount = memCount + (tot*procs*byte)
            TotalMem = TotalMem + memCount

            !check if array can be allocated
            if(memCount > limit)then
                print*,"Not enough ram to allocate: ",name
                error stop
            end if
            !allocate array
            if(present(start))then
                allocate(array(start(1):dims(1), start(2):dims(2), start(3):dims(3), start(4):dims(4)))
            else
                allocate(array(dims(1), dims(2), dims(3), dims(4)))
            end if
        end subroutine checkallocate_4R


        subroutine checkallocate_3R(array, dims, name, numproc, start)
        ! subroutine to allocate array if enough RAM is left
        ! Takes:
        !       n-dimensional arrayto be allocated         : array
        !       1D array length n, contains array dims     : dims
        !       character var, for error messgages         : name
        !       integer var, optional for calc if parallel : numproc

            implicit none

            real, allocatable, intent(INOUT) :: array(:,:,:)
            integer,           intent(IN)    :: dims(3)
            integer, optional, intent(IN)    :: numproc, start(3)
            character(*),      intent(IN)    :: name

            integer(int64) :: limit, procs, byte, tot, memCount

            real    :: tmp
            integer :: val, j

            !get current RAM available
            limit = mem_free()
            !set memory used counter to 0
            memCount = 0_int64

            !check if code is parallel
            if(.not.present(numproc))then
                procs = 1_int64
            else
                procs = numproc
            end if

            !get precision
            val = maxexponent(tmp)
            if(val == 1024)then
                !dbl prec
                byte = 8_int64
            elseif(val == 128)then
                !sgle prec
                byte = 4_int64
            else
                error stop
            end if

            !calc space for new array
            tot = dims(1)
            do j = 2, 3
                tot = tot * dims(j)
            end do

            ! calc memory space for array, and add to total counter
            memCount = memCount + (tot*procs*byte)
            TotalMem = TotalMem + memCount

            !check if array can be allocated
            if(memCount > limit)then
                print*,"Not enough ram to allocate: ",name
                error stop
            end if
            !allocate array
            if(present(start))then
                allocate(array(start(1):dims(1), start(2):dims(2), start(3):dims(3)))
            else
                allocate(array(dims(1), dims(2), dims(3)))
            end if
        end subroutine checkallocate_3R


        subroutine checkallocate_2R(array, dims, name, numproc, start)
        
            implicit none

            real, allocatable, intent(INOUT) :: array(:,:)
            integer,           intent(IN)    :: dims(2)
            integer, optional, intent(IN)    :: numproc, start(2)
            character(*),      intent(IN)    :: name

            integer(int64) :: limit, memCount, procs, byte, tot

            real    :: tmp
            integer :: val, j

            limit = mem_free()
            memCount = 0_int64

            if(.not.present(numproc))then
                procs = 1_int64
            else
                procs = numproc
            end if

            val = maxexponent(tmp)
            if(val == 1024)then
                !dbl prec
                byte = 8_int64
            elseif(val == 128)then
                !sgle prec
                byte = 4_int64
            else
                error stop
            end if

            tot = dims(1)
            do j = 2, 2
                tot = tot * dims(j)
            end do

            memCount = memCount + (tot*procs*byte)
            TotalMem = TotalMem + memCount
            
            if(memCount > limit)then
                print*,"Not enough ram to allocate: ",name
                error stop
            end if
            if(present(start))then
                allocate(array(start(1):dims(1), start(2):dims(2)))
            else
                allocate(array(dims(1), dims(2)))
            end if
        end subroutine checkallocate_2R


        subroutine checkallocate_1R(array, dims, name, numproc)
        
            implicit none

            real, allocatable, intent(INOUT) :: array(:)
            integer,           intent(IN)    :: dims(1)
            integer, optional, intent(IN)    :: numproc
            character(*),      intent(IN)    :: name

            integer(int64) :: limit, memCount, procs, byte, tot

            real    :: tmp
            integer :: val

            limit = mem_free()
            memCount = 0_int64

            if(.not.present(numproc))then
                procs = 1_int64
            else
                procs = numproc
            end if

            val = maxexponent(tmp)
            if(val == 1024)then
                !dbl prec
                byte = 8_int64
            elseif(val == 128)then
                !sgle prec
                byte = 4_int64
            else
                error stop
            end if

            tot = dims(1)

            memCount = memCount + (tot*procs*byte)
            TotalMem = TotalMem + memCount

            if(memCount > limit)then
                print*,"Not enough ram to allocate: ",name
                error stop
            end if
            allocate(array(dims(1)))
        end subroutine checkallocate_1R


        function mem_free()
        ! func that returns amount of free memory
        ! should be portable across most linux systems
        ! returns memory available in b

            implicit none

            integer(int64)    :: mem_free, available, pagecache, active, inactive, sreclaimable, freeram
            integer(int64)    :: i, low
            character(len=15) :: tmp
            integer           :: u, io


            open(newunit=u,file='/proc/zoneinfo',status='old')
            low = 0_int64
            sreclaimable = 0_int64
            inactive = 0_int64
            freeram = 0_int64
            active = 0_int64
            do 
                read(u,*, iostat=io)tmp, i
                if(IS_IOSTAT_END(io))exit
                if(verify(trim(tmp), 'low') == 0)then
                    low = low + i
                end if
            end do
            close(u)

            open(newunit=u,file='/proc/meminfo',status='old')
            mem_free = 0
            do 
                read(u,*,iostat=io)tmp, i
                if(IS_IOSTAT_END(io))exit

                if(verify(trim(tmp), 'MemAvailable:') == 0)then
                    mem_free = i * 1024_int64
                    close(u)
                    return!return early if MemAvailable availble on OS. Else do calculation for MemAvailable
                end if
                if(verify(trim(tmp), 'MemFree:') == 0)then
                    freeram = i
                end if
                if(trim(tmp) == 'Active(file):' )then
                    active = i
                end if
                if(verify(tmp, 'Inactive(file):') == 0)then
                    inactive = i
                end if
                if(trim(tmp) == 'SReclaimable:')then
                    sreclaimable = i
                end if
            end do
            close(u)

            !algorithm from kernal source
        !https://git.kernel.org/pub/scm/linux/kernel/git/torvalds/linux.git/commit/?id=34e431b0ae398fc54ea69ff85ec700722c9da773
            available = freeram - low
            pagecache = active + inactive
            pagecache = pagecache - min(pagecache/2, low)
            available = available + pagecache
            available = available + (sreclaimable - min(sreclaimable/2, low))
            mem_free = available * 1024_int64 !convert from Kib to b
        end function mem_free

end module memoryModule

! program test

!     use iso_fortran_env, only : int64

!     implicit none
    
!     integer(int64) :: limit
!     real,    allocatable :: threedR(:,:,:), twodR(:,:), onedR(:), threedR2(:,:,:),threedR3(:,:,:)

!     limit = mem_free()
    
!     !allocate 77% of RAM (~24.9Gb) on 31.2Gb system with single precision
!     !fails at "three" on double precision
!     call checkallocate(twodR, [100,100], "one", 1, [0,0])
!     twodR = 0.
!     twodR(0,0)=1.
!     print*,twodR(0,0)
!     call checkallocate(threedR2, [1000,1000,1000], "two", 1)
!     threedR2 = 0.
!     call checkallocate(threedR, [2000,2000,1000], "three", 1)
!     threedR=0.

!     call checkallocate(onedR, [1000], "four", 1)
!     onedR = 0.

!     call checkallocate(threedR3, [1000,1000,1000], "five", 1)
!     threedR3 = 0.

!     print'(A,1X,F5.2,A)','allocated:',dble(TotalMem)/dble(limit)*100.d0,' % of total RAM'

! end program test