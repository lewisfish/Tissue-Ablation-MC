MODULE iarray
!
!  Contains all array var names.
!
    implicit none
    save

    real, allocatable :: xface(:),yface(:),zface(:)
    real, allocatable :: rhokap(:,:,:)
    real, allocatable :: jmean(:,:,:),jmeanGLOBAL(:,:,:)
    real, allocatable :: temp(:,:,:), tissue(:,:,:)
    real, allocatable :: ThresTime(:,:,:,:), ThresTimeGLOBAL(:,:,:,:)

end MODULE iarray
