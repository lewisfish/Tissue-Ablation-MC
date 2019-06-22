MODULE ch_opt

implicit none

CONTAINS
     
   subroutine init_opt2
!
!  subroutine to set tissue optical properties 975nm
!
       use opt_prop
    
       implicit none
    
    
       hgg_epi = 0.87d0
       hgg_pap = 0.72d0
       hgg_ret = 0.72d0
       hgg_hypo = 0.78d0
    
       mua_epi = 0.78d0
       mua_pap = 0.66d0
       mua_ret = 0.64d0
       mua_hypo = 0.18d0
    
       mus_epi = 30.9d0 / (1.d0 - hgg_epi)
       mus_pap = 16.6d0 / (1.d0 - hgg_pap)
       mus_ret = 16.6d0 / (1.d0 - hgg_ret)
       mus_hypo = 5.3d0 / (1.d0 - hgg_hypo)
    
       n_epi = 1.42d0
       n_pap = 1.37d0
       n_ret = 1.37d0
       n_hypo = 1.37d0

   end subroutine init_opt2
   
   subroutine init_opt3
!
!  subroutine to set tissue optical properties 1064nm
!
       use opt_prop

       implicit none

       hgg_epi = 0.89d0
       hgg_pap = 0.72d0
       hgg_ret = 0.72d0
       hgg_hypo = 0.78d0
    
       mua_epi = 0.51d0
       mua_pap = 0.38d0
       mua_ret = 0.38d0
       mua_hypo = 0.08d0
    
       mus_epi = 28.7d0 / (1.d0 - hgg_epi)
       mus_pap = 14.7d0 / (1.d0 - hgg_pap)
       mus_ret = 14.7d0 / (1.d0 - hgg_ret)
       mus_hypo = 5.0d0 / (1.d0 - hgg_hypo)
    
       n_epi = 1.42d0
       n_pap = 1.37d0
       n_ret = 1.37d0
       n_hypo = 1.37d0

   end subroutine init_opt3
   
   subroutine init_opt1
!
!  subroutine to set tissue optical properties 532nm
!
       use opt_prop
       
       implicit none

       hgg_epi = 0.75d0
       hgg_pap = 0.71d0
       hgg_ret = 0.71d0
       hgg_hypo = 0.78d0
    
       mua_epi = 5.42d0
       mua_pap = 3.55d0
       mua_ret = 2.90d0
       mua_hypo = 0.84d0
    
       mus_epi = 64.3d0 / (1.d0 - hgg_epi)
       mus_pap = 39.7d0 / (1.d0 - hgg_pap)
       mus_ret = 39.7d0 / (1.d0 - hgg_ret)
       mus_hypo = 9.2d0 / (1.d0 - hgg_hypo)
    
       n_epi = 1.45d0
       n_pap = 1.39d0
       n_ret = 1.39d0
       n_hypo = 1.37d0

   end subroutine init_opt1

   subroutine sample(array, size_of, cdf, wave, iseed)
!      
!  samples a random value from an array based upon its cdf     
!      
      implicit none
      
      integer, intent(IN)    :: iseed, size_of
      real,    intent(IN)    :: array(size_of, 2), cdf(size_of)
      real,    intent(OUT)   :: wave

      real :: ran2, value
      integer :: nlow
      
      value = ran2(iseed)
      
      call search_1D(size(cdf), cdf, nlow, value)
      call lin_inter_1D(array, cdf, value, size(cdf), nlow, wave)
   
   end subroutine sample
   
   subroutine lin_inter_1D(array, cdf, value, length, nlow, y)
!
!  linear interpolates between values for an array and its cdf
!   
      implicit none
   
      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2),cdf(length-1)
      integer, intent(IN)   :: nlow
   
      y = array(nlow+1,1) + (array(nlow+2,1) - array(nlow+1,1)) * (value - cdf(nlow))/(cdf(nlow+1) - cdf(nlow))
   
   end subroutine lin_inter_1D
   
   subroutine lin_inter_2D(array,value,length,nlow,y)
!
!  linear interpolation for an array
!
      implicit none

      real,    intent(OUT)  :: y
      integer, intent(IN)   :: length
      real,    intent(IN)   :: value,array(length,2)
      integer, intent(IN)   :: nlow
   
      y = array(nlow,2) + (array(nlow+1,2) - array(nlow,2)) * (value - array(nlow,1))/(array(nlow+1,1) - array(nlow,1))
   
   end subroutine lin_inter_2D
   
   subroutine search_1D(length,array,nlow,value)
!
!  search by bisection for 1D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_1D
   
   subroutine search_2D(length,array,nlow,value)
!
!  search by bisection for 2D array
!
      implicit none
      
      integer              :: nup,length,middle
      integer, intent(OUT) :: nlow
      real,    intent(in)  :: array(length,2),value
      
      nup = length
      nlow = 1
      middle = int((nup+nlow)/2.)

      do while((nup - nlow).gt.1)
         middle = int((nup + nlow)/2.)
         if(value.gt.array(middle,1))then
            nlow = middle
         else
            nup = middle   
         end if
      end do
   end subroutine search_2D
   
   subroutine mk_cdf(array,cdf,length)
!
!  subroutine that creates cdf for an array of values.
!
      implicit none

      integer, intent(IN)    :: length
      real,    intent(IN)    :: array(length,2)
      real,    intent(INOUT) :: cdf(length)
      real                   :: summ
      integer                :: i,j
   
      do j=1,length-1
         summ=0.
         do i=1,j   
            summ=summ+0.5*(array(i+1,2)+array(i,2))*(array(i+1,1)-array(i,1))
         end do
         cdf(j)=summ      
      end do
      cdf=cdf/cdf(length-1)
   
   end subroutine mk_cdf
end module ch_opt
