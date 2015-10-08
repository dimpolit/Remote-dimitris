 MODULE bio
  USE funCD2TT
  USE funTT2CD
  USE funND2TT
  USE subTT2ND 
  
 !use countdate, only: CD2TT,TT2CD,ND2TT,TT2ND, larcalkl,juvcalkl,calkl
 IMPLICIT NONE
 PRIVATE
 PUBLIC  Anchovy 
! integer,parameter :: ncc=1
! real(8), public :: xfi(ncc) , bl(ncc)

 CONTAINS

Subroutine Anchovy(p,First,TTime, Tbefore, z45, z55, z65, temp5, z410, z510, z610, temp10, z420, z520, z620, temp20, &
z425, z525, z625, temp25, z435, z535, z635, temp35)
  !
  use comibm
  !use countdate, only: CD2TT,TT2CD,ND2TT,TT2ND,larcalkl,juvcalkl,calkl
  implicit none
   REAL(4), DIMENSION(:), INTENT(IN) :: p
   integer          :: First
  integer,parameter :: nc =1  ! number of fish cohort
  integer,parameter :: nl =6   ! number of fish life stages
  real(8)           :: CTM, CTO,CTOj, CQ, RTM, RTO, RQ ! parameters
  real(8)           :: TTime, Tbefore, z4, z5, z6, temp1,temp
  real(8)           :: z45,z410,z420,z425,z435
  real(8)           :: z55,z510,z520,z525,z535
  real(8)           :: z65,z610,z620,z625,z635
  real(8)           :: temp5,temp10,temp20,temp25,temp35
  real(8),parameter :: d2s= 86400.0d0    ! day ---> sec
  integer           :: Iyr, Imon, Iday, Ihour, Imin, Isec
  character(19)     :: StAge(nl)
  character(19)     :: CTime
  real(8)           :: TAge(nl)
  integer           :: iage(nc)
  integer, save     :: IyrA(nl), ImonA(nl), IdayA(nl) ,IhourA(nl), IminA(nl), IsecA(nl)
  integer           :: JJday,CJJday,JJJday,hourday
  integer           :: ii,i,ic
  real(8)           :: ZooP1, ZooP2, ZooP3, tt1, ttt1
  real(8)           :: t1,t2 !, calkl,larcalkl, juvcalkl
  real(8)           :: TTimes(nc),TAgecur,TAgebef, CTTimes(nc)
  real(8)           ::  xdot(nc),eggprod(nc),blb(nc)
 ! real(8)           :: cd2tt, nd2tt, xfi(nc),bl(nc)
 ! character(19)     :: tt2cd
  real(8)           :: mcAge,mcTTimes
  real(8)           :: v1,x1,y1,z1,v2,x2,y2,z2,act(nc)
  real(8)           :: vul(3)=0.0, k(3,nl)
  real(8)           :: cfac(nc)
  real(8)           :: v, a, u, resp(nc)
  real(8)           :: gctemp,gcmaxfi(nc),grtemp
  real(8)           :: cnum,c1,c2,c3,con1,con2,con3,con(nc)
  real(8)           :: f(nc),e(nc), sda(nc), wtemp(nc), wwtemp(nc)
  real(8)           :: edf(nc) ! energy density of fish; depends partially on length
  real(8)           :: Tm
  real(4)           :: co1,r1,f1
  real(4)           :: zz1,zz2,zz3
!  integer, save     :: First = 1
  character(19)    ::  Cmc = '1990/06/10 00:00:00' !output montecarlo
!  REAL(8), DIMENSION(:), INTENT(IN) :: p
  integer, save     :: Iyrmc, Imonmc, Idaymc ,Ihourmc, Iminmc, Isecmc
  !
  !     ==================

!  print*, 'InitF', First
  if ( First .eq. 1 ) then; 
          First = 0
!         StAge(1) ='0000/06/10 00:00:00' ! embryonic stage=eggs hatching(1.5-3 days)+ Yolk-sac stage ~6 days
!         StAge(2) ='0000/06/16 00:00:00' ! larval stage ~ 2 months
!         StAge(3) ='0000/08/26 00:00:00' ! juvenile until first year
!         StAge(4) ='0001/06/10 00:00:00' ! Adult stage
!         StAge(5) ='0002/06/10 00:00:00'
!         StAge(6) ='0003/12/10 00:00:00' ! End at 3 years+2months
        StAge(1) ='0000/11/01 02:00:00' ! embryonic stage=eggs hatching(1.5-3 days)+ Yolk-sac stage ~6 days
        StAge(2) ='0000/11/05 00:00:00' ! larval stage ~ 2 months
        StAge(3) ='0000/08/16 00:00:00' ! juvenile until first year
        StAge(4) ='0001/06/01 00:00:00' ! Adult stage
        StAge(5) ='0002/06/01 00:00:00'
        StAge(6) ='0003/12/01 00:00:00' ! End at 3 years+2months

       do ii=1,nl   !6 life stages
          TAge(ii) = CD2TT( StAge(ii) )
          call TT2ND(IyrA(ii), ImonA(ii), IdayA(ii) ,IhourA(ii), IminA(ii), IsecA(ii) ,TAge(ii) )
          TAge(ii) = CD2TT( StAge(ii) ) - CD2TT( StAge(1) )
!  pernaei aytes tis xronikes stigmes mesa sto kodika.
!        write(*,*) tage(ii)/60.0/60.0/24.0
        end do

         do i=1,nc      ! gia ta 3 cohort ta mhdenizei ola
          xfi(i)=0.0d0
          bl(i)=0.0d0
          blb(i)=0.0d0
          iage(i)=0
         end do
        
         
 mcAge = CD2TT(Cmc)
 call TT2ND(Iyrmc, Imonmc, Idaymc ,Ihourmc, Iminmc, Isecmc ,mcAge )

     open (21,  file='sardine.txt')
     open (22,  file='zoop1.txt')
     open (223, file='zoop2.txt')
     open (24,  file='zoop3.txt')
     open (30,  file='sardine.bin', form='UNFORMATTED')
     open (101, file='processes.bin', form='UNFORMATTED')
     open (111,  file='processes.txt')
     open (103, file='blen.bin', form='UNFORMATTED')
     open (44,  file='zoo.txt')
     open (55,  file='wtemp.bin',form='UNFORMATTED')
!     print*, 'IyrA, IdayA', IyrA, IdayA
    !
    !  k(1)=p(1)
    !  k(2)=p(2)
    !  k(3)=p(3)

    k(1,2)=p(1);k(2,2)=p(2);k(3,2)=p(3); !Tm=p(4);
    CTO=p(5);k(1,3)=p(6);k(2,3)=p(7);CTOj=p(8)    

  end if ! of First
!     =====================================
  Tm=17.0
  !
  CTime = TT2CD(TTime)                                 ! present time (charactor form)
!  print*, 'ctimeFish', CTime  !ana ora trexei o kodikas
  CALL TT2ND(Iyr, Imon, Iday ,Ihour, Imin, Isec ,TTime)
  JJday = 1 + (TTime - ND2TT(1990,1,1,0,0,0) ) / d2s
! print*, jjday
    do i=1,nc           !gia ta 3 cohorts
    if(iage(i) .eq. 0) then
           TTimes(i)=ND2TT(Iyr, ImonA(1), IdayA(1) ,IhourA(1), IminA(1), IsecA(1))
!         TTimes(2)=ND2TT(Iyr, ImonA(4), IdayA(4) ,IhourA(4), IminA(4), IsecA(4))
!     TTimes(3)=ND2TT(Iyr, ImonA(1), IdayA(1) ,IhourA(1), IminA(1), IsecA(1))
! PROSOXH: To kathe cohort exei diafora 1 etous,me ton parapano diaxorismo
! ksekinaei opote theloume
!         CTTimes(i)=TT2CD(TTimes(i))
!         print*, 'TTimes(i)', TTimes/60.0/60.0/24.0
!           print*, 'CTT', iage(1), Tbefore, TTime,TTimes(1)
!           print*, 'Timing',First,Iyr,CTime,TTime
!     pause
        if( (TTime .ge. TTimes(i)) .and. (Tbefore .lt. TTimes(i))  ) then !.and. (Tbefore .lt. TTimes(i)) 
               iage(i) = 1
               TTimes(i)=TTime
               xfi(i)=2.7962e-05
               bl(i)= larcalkl(xfi(i))
               blb(i)=bl(i)
    print*, 'Aging to 1 ', ctime,bl(i),xfi(i)
 !           goto 998
            endif
         endif
       end do
!998 continue
       do i=1,nc
       if( iage(i) .eq. 0 ) then
          xdot(i)=0.0d0
          goto 999
       end if

JJday= 1 + ( TTime - ND2TT(Iyr, ImonA(1), IdayA(1) ,IhourA(1), IminA(1), IsecA(1) )) / d2s
 mcTTimes=ND2TT(Iyr, Imonmc, Idaymc ,Ihourmc, Iminmc, Isecmc)

               if (bl(i) .lt. 45.) then
                    edf(i)=8000.0
               else if ((bl(i) .ge. 45.) .and. (bl(i) .lt. 115.)) then
                    edf(i)=6000.0
               else 
                    edf(i)=5600.0
               endif
! print*, ctime, edf(1)

!     ..... Aging of Sardine ......
        TAgecur = TTime-TTimes(i)
!        TAgebef = Tbefore-TTimes(i)
!        CTTimes=TT2CD(TTimes)
!   print*, 'CTTimes',  CTTimes
! perasma kai diarkeia tou larval stage
         if (iage(i) .eq. 1)  then
           if( (TAgecur .ge. TAge(iage(i)+1)) .and.  (bl(i) .gt. 3.5)) then
                  iage(i) = iage(i) + 1
          write(*,*)  'Aging +2 ',  CTime,  i, iage(i), xfi(i), bl(i)
!            write(*,*) 'CTime', CTime, xfi(1)
       endif
        end if
! perasma k diarkeia tou juvenile stage
        if (iage(i) .eq. 2)  then
          if (bl(i) .ge. 45.0) then
           iage(i) = iage(i) + 1
         write(*,*)'3Aging3 ', JJday, xfi(1), bl(1)
          endif
        endif
! perasma kia diarkeia tou adult stage

        if (iage(i) .eq. 3) then
          if(bl(i) .gt. 114.8) then 
            iage(i) = iage(i) + 1
       write(*,*) 'Aging to 4 ', CTime, xfi(1), bl(1)
          endif
        endif

  JJJday = 1 + ( TTime - ND2TT(Iyr ,1,1,0,0,0) )/ d2s
!  JJday= 1 + ( TTime - ND2TT(Iyr, ImonA(1), IdayA(1) ,IhourA(1), IminA(1), IsecA(1) )) / d2s
  if (JJday .lt. 0)then
  JJday=JJday+364
  endif
  ! to jjday den einai ayksousa akolouthia
  !------convert Nemuro zoop in uM N/L to g ww/m3
  !------tt1 is conversion from uM N/liter to g ww/m3
  !------14 ug N/uM * 1.0e-6 g/ug * 1 g dw/0.07 g N dw * 1 g ww/0.2 g dw
  !------1.0e6 liters/m3
  !
  hourday=2+( TTime - ND2TT(Iyr ,Imon,Iday,1,0,0) )/ 3600

        tt1=14.0*1.0e-6*(1.0/0.07)*(1.0/0.2)*1.0e3

  !     ..... Temperature Setting .....
   If ((iage(i) .eq. 1) .or. (iage(i) .eq. 2)) then
     zoop1 =(z410*tt1)
     zoop2 =(z510*tt1)
     zoop3 =(z610*tt1)
     wtemp(i)=temp10
! print*, CTime,zoop1
    elseif  (iage(i) .eq. 3) then
       if (hourday .le. 12) then
    zoop1 =(z420*tt1)
    zoop2 =(z520*tt1)
    zoop3 =(z620*tt1)
    wtemp(i)=temp25
        else
    zoop1 =(z410*tt1)
    zoop2 =(z510*tt1)
    zoop3 =(z610*tt1)
    wtemp(i)=temp10
       endif
    elseif ((iage(i) .eq. 4) .or. (iage(i) .eq. 5) .or. (iage(i) .eq. 6)) then
           if (hourday .le. 16) then
   zoop1 =(z435*tt1)
   zoop2 =(z535*tt1)
   zoop3 =(z635*tt1)
   wtemp(i)=temp35
   else
   zoop1 =(z420*tt1)
   zoop2 =(z520*tt1)
   zoop3 =(z620*tt1)
   wtemp(i)=temp20

      endif
endif

  IF(wtemp(i).le.1.0) wtemp(i)=1.0

  !
  !--- Sardine weight state variable = xfi(i)
  !
  !----- set vulnerabilities and k values for 4 zoop groups

  If (iage(i) .eq. 1) then
  vul(1) =   0.0; vul(2) =  0.0; vul(3) =   0.0
  k (1,iage(i)) = 0.1; k (2,iage(i)) = 0.1; 
  k (3,iage(i)) = 0.1;
  elseif (iage(i) .eq. 2) then
  vul(1) =   0.3; vul(2) =  1.0; vul(3) = 0.5
!  k (1,iage(i)) = 0.01*0.25; !k  (2,iage(i)) = 0.045*0.25; k (3,iage(i)) = 0.045
  elseif (iage(i) .eq. 3) then
  vul(1) = 0.7; vul(2) = 0.7; vul(3) = 0.0
  k (3,iage(i)) = 0.1;  ! k (2,iage(i)) = 0.1; k (1,iage(i)) = 0.2;
  else
  vul(1) =  1.0; vul(2) =  0.5; vul(3) =  0.0
!  k (1) = 0.06; k (2) = 0.06; 
  k(3,iage(i)) = 0.15; k(2,iage(i)) = 0.15;k(1,iage(i)) = 0.15
  END IF

  ! --- The 5.258 puts resp (g oxygen/fish) into units of g zoop/g fish/day
  ! --- [13560 joules/gram oxygen]/4.18 joules/cal = 3244 cal/gO2
  ! --- [2580 joules/gram zoop]/4.18 joules/cal = 617 cal/g zoop
  ! --- so respiration in grams/oxygen/g fish/day is multiplied by 3244/617 = 5.258
  ! --- to get food energy equivalents of a gram of oxygen respired

   If  (iage(i) .eq. 2) then
        CTM=25.0;  CQ=2.22;  RTO=28.0; RTM=33.0; RQ=2.22
  elseif  (iage(i) .eq. 3)  then
         CTM=25.0;  CQ=2.22; RTO=25.0; RTM=33.0;  RQ=2.22
  else
        CTO=20.0; CTM=32.0; CQ=2.22; RTO=25.0; RTM=33.0;  RQ=2.22
  endif
!print*, Ctime, wtemp
!optimal and maximum temperature for consumption and respiration
!depends on life stage
 ttt1 = 1.0 / xfi(i)
! print*, tt1

!--- temperature-dependent function for maximum consumption
   z1=(CTM-CTO)*log(CQ)
   y1=(CTM-CTO+2.0)*log(CQ)
   v1=(CTM-wtemp(i))/(CTM-CTO)
   x1=((z1**2.0)*((1.0+sqrt(1.0+40.0/y1))**2.0))/400.0
   gctemp=(v1**x1)*exp(x1*(1.0-v1))
! print*, ctime, gctemp
!---Maximum Consumption
  if (iage(i) .eq. 1) then
       gcmaxfi(i)=0.0d0
     else
       gcmaxfi(i) = 0.24 * ttt1**0.342* gctemp
    end if
 !print*,  ctime,  gcmaxfi(i)
 grtemp=exp((wtemp(i)-Tm)*log(2.2)/10)

 if (ii .eq. 2) then
          u=9.1*xfi(i)**0.24
        elseif (ii .eq. 3) then
          u=5.32*xfi(i)**0.32
        else
          u=10.73*xfi(i)**0.31
        endif


   if (iage(i) .eq. 1) then
    resp(i)=0.0d0
   else
    resp(i)=0.0033*(ttt1**0.227)*grtemp*5.258*exp(0.003*u)
   endif
! print*, ctime, resp(1)
  ! --- multispecies functional response
  ! --- use either this or adjust little p
  ! ---consumption
  cnum=zoop1 * vul(1)/k(1,iage(i))+zoop2*vul(2)/k(2,iage(i))+zoop3*vul(3)/k(3,iage(i))
  c1=gcmaxfi(i)*zoop1*vul(1)/k(1,iage(i))
  c2=gcmaxfi(i)*zoop2*vul(2)/k(2,iage(i))
  c3=gcmaxfi(i)*zoop3*vul(3)/k(3,iage(i))
!
  con1=c1/(1.0+cnum)
  con2=c2/(1.0+cnum)
  con3=c3/(1.0+cnum)
!
  con(i)= con1+con2+con3
  ! --- egestion ---
  !
  wwtemp(i)=1.0/wtemp(i)

  if (iage(i) .eq. 1) then
        f(i)=0.0d0
     elseif (iage(i) .eq. 2) then
        f(i)=0.1*con(i)
  else
        f(i)=0.1*con(i)
  end if
  !
  ! --- excretion
  if (iage(i) .eq. 1) then
    e(i)=0.0d0
  elseif (iage(i) .eq. 2) then
   e(i)=0.16*(con(i)-f(i))
  else
   e(i)=0.16*(con(i)-f(i))
  end if
  !
  ! --- Specific Dynamic Action
  !
   if (iage(i) .eq. 1) then
     sda(i)=0.0d0
   elseif (iage(i) .eq. 2) then
    sda(i)=0.175*(con(i)-f(i))
    else
     sda(i)=0.175*(con(i)-f(i))
   end if

  !C --- use the ratio of calories/g of zoop (2580) to calories/g of fish (edf)
  !C
  !C --- bioenergetics differential equation
  !C

     xdot(i)=(con(i)-resp(i)-f(i)-e(i)-sda(i))*xfi(i)*2580./edf(i)
       write(*,*) 'Bio',con(i),resp(i)

  IF(wtemp(i).le.1.0) xdot(i)=0.0
  !
  !     Time Integration
  !
 999  xfi(i) = xfi(i) + 3600.0d0/d2s * xdot(i)

        if (iage(i) .eq. 2)   then
            bl(i)= larcalkl(xfi(i))
         elseif (iage(i) .eq. 3) then
            bl(i)= juvcalkl(xfi(i))
         else
            bl(i)= 10**calkl(xfi(i))
         end if

      cfac(i)=xfi(i)/bl(i)

        if( bl(i) .gt. blb(i) ) then
          blb(i)=bl(i)
          else
          bl(i)=blb(i)
        end if
      end do

       co1=con(1)*xfi(1)
       r1=resp(1)*xfi(1)
       f1=f(1)*xfi(1)
       zz1=zoop1
       zz2=zoop2
       zz3=zoop3
  !     ..... for Check .....
!  if( (TTime .ge. mcTTimes) .and. (Tbefore .lt. mcTTimes) ) then
  if ( int(TTime/d2s) .ne. int(Tbefore/d2s) ) then
     !!      write(*,*) TZS, zop1(JJday), TZL,zop2(JJday), TZP,zop3(JJday)
     !!      write(*,*) TZP*1.0d6, zop3(JJday)
     write (21,*) JJday, CTime, xfi(1),bl(1)
     write (30)  real(xfi(1))
     write (30)  real(bl(1))
     write(111,*) JJday,co1,r1,f1
     write (101) real(con(1)*xfi(1))
     write (101) real(resp(1)*xfi(1))
     write (101) real(f(1)*xfi(1))
     write (101) real(e(1)*xfi(1))
     write (101) real(sda(1)*xfi(1))
     write (101) real(eggprod(1)*xfi(1))
     write (22,*) CTime,zoop1
     write (223,*) CTime,zoop2
     write (24,*) CTime,zoop3
     write (44,*)  JJday,zz1,zz2,zz3
     write(55)   real(temp1)
 !    write (55)  real(wtemp)
!        print*,  jjday, CTime, xfi(1),bl(1)
  end if
! endif
  return
  !

 ! stop
end Subroutine Anchovy

!***********************************************
real(8) function larcalkl(x)
!
      real(8) :: x, x1, b0, b1
!
!      b0=5.5*(10**6)
!      b1=4.0276
!from NAEG
       b0=45.94
       b1=4.030
      larcalkl=b0*(x**(1.0/b1))
!
      return
      end function
!***********************************************
 real(8) function juvcalkl(x)
!
      real(8) :: x, x1, b0, b1,b2,d1
!
      b0=183823.52
      b1=3.05
      juvcalkl=(b0*x)**(1/b1)
!
      return
      end function
!************************************************
real(8) function calkl(x)
!
      real(8) :: x, x1, b0, b1,b2,b3,d1,d2

      x1=log10(x)
      b0=-6.11580
      b1=3.57642
      b2=-0.616
      b3=0.71370
      d1=1.5798
      d2=2.0212

     calkl=(x1-b0+b2*d1+b3*d2)/(b1+b2+b3)
!
      return
      end function
!*************************************************
!=============================================

END MODULE bio
