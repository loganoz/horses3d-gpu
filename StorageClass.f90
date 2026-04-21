! Updated StorageClass.f90 with mu_artificial allocation

MODULE StorageClass
  IMPLICIT NONE

  ! Declaration of variables
  REAL :: mu_artificial

  ! Constructor
  CONTAINS

  SUBROUTINE allocate_mu_artificial()
    mu_artificial = 0.0
  END SUBROUTINE allocate_mu_artificial

END MODULE StorageClass