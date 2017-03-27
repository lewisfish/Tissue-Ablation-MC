MODULE iarray
!
!  Contains all array var names.
!
    implicit none
    save

    real, allocatable    :: xface(:),yface(:),zface(:)
    real, allocatable    :: rhokap(:,:,:), absorb(:,:,:)
    real, allocatable    :: jmean(:,:,:),jmeanGLOBAL(:,:,:)
end MODULE iarray
