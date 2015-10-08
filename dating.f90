! First module ******************************
module funND2TT
PUBLIC ND2TT

CONTAINS

!*************************************
!     exp. 1997,12,31,23,59,59 \rightarrow 6.223158719900000E+10
!*************************************
real(8) function ND2TT(Iyr, Imon, Iday, Ihour, Imin, Isec)
  !
  integer   :: IM2D(12,0:1) = &
       reshape( (/ 0,31,59,90,120,151,181,212,243,273,304,334,  &
       0,31,60,91,121,152,182,213,244,274,305,335  /), (/12,2/) )
  integer   :: Iyr, Imon, Iday, Ihour, Imin, Isec
  integer   :: Iy4, Iy1, Ileap, Im, Itt
  !
  Iy4 = 1461 * ( Iyr / 4 )
  Iy1 = 365 * mod( Iyr, 4 )
  !
  if ( mod( Iyr, 4 ) .ne. 0 ) then
     Ileap = 0
  else
     Ileap = 1
  end if
  Im = IM2D( Imon, Ileap)
  !
  Itt = Iy4 + Iy1 + Im + Iday - Ileap
  !
  ND2TT = Ihour * 3600 + Imin * 60 + Isec
  ND2TT = ND2TT + Itt * 86400.0D0
  !
  return
end function ND2TT

end module funND2TT

! 2nd module ********************************************
module funCD2TT
   
   IMPLICIT NONE
   PUBLIC  CD2TT

CONTAINS

!********************************
!* Utilities for Date Control  Writtien by Yasuhiro Yamanaka (galapen@ees.hokudai.ac.jp)
!*********************************
!     exp. 1997/12/31 23:59:59 \rightarrow 6.223158719900000E+10
!     exp. 0000/01/01 00:00:00 \rightarrow 0.000000000000000E+00
!***********************************
real(8) function CD2TT( Cdate )
  !
  USE funND2TT,only : ND2TT
  integer       :: Iyr, Imon, Iday , Ihour, Imin, Isec
 ! real(8)       :: ND2TT
  character(19) :: Cdate
  !
  if ( len( Cdate ) .ne. 19 ) then
     write(*,*) ' Length of date is no good '
     stop
  end if
  read (Cdate( 1: 4),*)  Iyr
  read (Cdate( 6: 7),*)  Imon
  read (Cdate( 9:10),*)  Iday
  read (Cdate(12:13),*)  Ihour
  read (Cdate(15:16),*)  Imin
  read (Cdate(18:19),*)  Isec
  !
  CD2TT = ND2TT(Iyr, Imon, Iday , Ihour, Imin, Isec)
  !
  return
end function CD2TT

end module funCD2tt


!3rd*************************************
module subTT2ND

PUBLIC TT2ND

CONTAINS
!**********************************
!     exp. 6.223158719900000E+10 \rightarrow 1997,12,31,23,59,59
!***********************************
subroutine TT2ND(                                        &
     Iyr   , Imon  , Iday   , Ihour, Imin, Isec,   & !O & I
     tt    )
  !
  USE funND2TT, only : ND2TT
  integer :: Iyr, Imon, Iday , Ihour, Imin, Isec
  integer :: Itt, Iy, Iy4, Iyd, Iy1, Ileap, Imd, Im, Its
  integer   :: IM2D(12,0:1) = &
       reshape( (/ 0,31,59,90,120,151,181,212,243,273,304,334,  &
       0,31,60,91,121,152,182,213,244,274,305,335  /), (/12,2/) )
  integer :: IY2D(4) = (/0,366,731,1096/)
  real(8) :: tt, tt0 !, ND2TT
  !
  !     ..... ITT [day] .....
  Itt = 1 + tt / 86400.0D0
  !
  Iy4   = (Itt-1) / 1461
  Iyd   = Itt - Iy4 * 1461
  do IY = 1, 4
     if ( IY2D(Iy) + 1 .le. Iyd ) then
        Iy1 = Iy
     end if
  end do
  !
  Iyr   = Iy4 * 4 + Iy1 - 1
  if ( mod(Iyr,4) .ne. 0 ) then
     Ileap = 0
  else
     Ileap = 1
  end if
  IMD = IYD - IY2D(IY1)
  !
  do IM = 1, 12
     if ( IM2D(IM,ILEAP)+1 .le. IMD ) then
        IMON = IM
     end if
  end do
  IDAY = IMD - IM2D(IMON,ILEAP)
  !
  TT0 = ND2TT(IYR, IMON, IDAY ,0,0,0)
  ITS = nint(  TT - TT0 )
  Ihour = ITS / 3600
  Imin  = ( ITS - Ihour * 3600 ) / 60
  Isec  = ITS - Ihour * 3600 - Imin * 60
  !
  return
end subroutine TT2ND

end module subTT2ND


! 4rth module
module funTT2CD

PUBLIC TT2CD

CONTAINS


!*************************************
!     exp. 6.223158719900000E+10 \rightarrow 1997/12/31 23:59:59
!*************************************

 character(19) function TT2CD(tt)
  !
  USE subTT2ND, only: TT2ND
  integer :: Iyr, Imon, Iday , Ihour, Imin, Isec
  real(8) :: tt
  !
  call TT2ND( Iyr, Imon, Iday, Ihour, Imin, Isec , tt )
  !
  write(TT2CD,'(I4.4,5(A,I2.2))') Iyr, '/', Imon, '/', Iday, &
       ' ', Ihour, ':', Imin, ':', Isec
  !
  return
end function TT2CD

end module funTT2CD



