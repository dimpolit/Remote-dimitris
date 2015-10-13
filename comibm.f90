MODULE comibm
!---------------------------------------------------------------------------------
  
  USE com
  !USE com, ONLY : rsh,rlg,lchain

  IMPLICIT NONE
  
 
                  ! number of patches
   INTEGER                                      :: nb_patch
                   ! number of particles inside each patch
   
   public

      integer, parameter :: ncohort=1 !nages=6, npreys=3
 !     real, parameter :: fishsteps=3 !(3*1200=1hour)
 !     integer, parameter :: nteggs=2

      real(8),dimension(ncohort)  :: xfi,bl,iage
      real(8),dimension(ncohort)  :: predator, boat        

END MODULE comibm
