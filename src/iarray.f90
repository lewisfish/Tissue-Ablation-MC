MODULE iarray
!
!  Contains all array var names.
!
    implicit none
    save

    real, allocatable :: xface(:),yface(:),zface(:)
    real, allocatable :: rhokap(:), albedo(:), muas(:)
    real, allocatable :: jmean(:,:,:),jmeanGLOBAL(:,:,:)
    real, allocatable :: temp(:,:,:), tissue(:,:,:)
    real, allocatable :: g2(:), refrac(:), hgg(:)

end MODULE iarray
