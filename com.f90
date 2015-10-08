MODULE com

  IMPLICIT NONE

    integer,parameter :: rsh = 4
    integer,parameter :: rlg = 8

    INTEGER, parameter       :: nread=10000
    integer,parameter :: nmax = 620, lchain=100
    REAL(kind=rsh),parameter :: pi=3.14159265358979, grav=9.810
    real(kind=rsh),allocatable,dimension(:) :: sig,dzu,dz,dzw,sigw,z
    real(kind=rlg) :: t,deltaZ
    integer :: ndtz
    real(kind=rsh) :: latitude, longitude
    REAL(kind=rsh)  :: h0
    INTEGER :: kmax,nbouc
    real(kind=rlg) ,PARAMETER :: tref=59958230400_rlg
    INTEGER, allocatable,dimension(:) :: ncid
!    LOGICAL, PUBLIC :: ctd,col
    CHARACTER(LEN=5), allocatable, dimension(:) :: station
    real(kind=rsh),parameter       :: valmanq=-999.0
    CHARACTER(LEN=19) :: date
    
!    CHARACTER(LEN=5)  :: stat
!    CHARACTER(LEN=130) :: rec1

   
END MODULE com

