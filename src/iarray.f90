MODULE iarray
!
!  Contains all array var names.
!
    implicit none
    save

    real, allocatable :: xface(:),yface(:),zface(:)
    real, allocatable :: rhokap(:), albedo(:)
    real, allocatable :: jmean(:,:,:),jmeanGLOBAL(:,:,:)
    real, allocatable :: temp(:,:,:), tissue(:,:,:)

end MODULE iarray
