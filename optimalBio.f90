program stab

  USE com
  USE optim
  USE fcost_mod
 
 implicit none

 INTERFACE
     FUNCTION mod1dv(p)
       USE com
       IMPLICIT NONE
       REAL(KIND=rsh), DIMENSION(:), INTENT(IN) :: p
       REAL(KIND=rsh) :: mod1dv
     END FUNCTION mod1dv
  END INTERFACE
       REAL(KIND=rsh) :: r
!--------------------------------------------------------------
! Initialisation ----------------------------------------------
!----------------------------------------------------------------

        optimisation=.false.
        file_optim= './optimisation/result_cost_food.txt'
        met_optim='amoeba'
        npara=8
! Define initial parameters for optimization

 allocate(para_init(npara))
  
  para_init(1)=0.02_rsh  !
  para_init(2)=0.001_rsh
  para_init(3)=0.04_rsh
  para_init(4)=0.2_rsh
  para_init(5)=0.1_rsh  ! 

  para_init(6)=17.0_rsh
  para_init(7)=16.0_rsh
  para_init(8)=18.0_rsh
  

 if (optimisation) then   ! boucle optimisation
     allocate(para_ref(npara),l_bound(npara),h_bound(npara))

  para_ref(1)=0.025_rsh  !
  para_ref(2)=0.0011_rsh
  para_ref(3)=0.035_rsh
  para_ref(4)=0.018_rsh
  para_ref(5)=0.12_rsh  ! 

  para_ref(6)=16.5_rsh
  para_ref(7)=15.5_rsh
  para_ref(8)=18.5_rsh


  l_bound(1)=0.000001_rsh
  h_bound(1)=0.8_rsh
  l_bound(2)=0.000001_rsh
  h_bound(2)=0.8_rsh
  l_bound(3)=0.0001_rsh
  h_bound(3)=0.8_rsh
  l_bound(4)=0.0001_rsh
  h_bound(4)=0.8_rsh
  l_bound(5)=0.0001_rsh
  h_bound(5)=0.8_rsh

  l_bound(6)=15.0_rsh
  h_bound(6)=20.0_rsh
  l_bound(7)=15.0_rsh
  h_bound(7)=20.0_rsh
  l_bound(8)=16_rsh
  h_bound(8)=20.0_rsh
     
     CALL do_optim(mod1dv,para_init)
     else ! run classique
     r=mod1dv(para_init)
  endif
      
  stop
end program stab
!**************************************

FUNCTION mod1dv(p)

  USE comibm
  USE optim, only : optimisation, para_ref,h_bound,l_bound,npara
  USE fcost_mod
  !use countdate, only: CD2TT,TT2CD,ND2TT,TT2ND
  USE funCD2TT
  USE funTT2CD
  USE funND2TT
  USE subTT2ND  
  USE bio, only : Anchovy

 !! * Arguments
  REAL(KIND=rsh) :: mod1dv
  REAL(KIND=rsh), DIMENSION(:), INTENT(IN) :: p

 !! * Private variables
  INTEGER :: nb_iter,nb,st,n,status,try,i
  CHARACTER(LEN=19) :: tool_sectodat
  CHARACTER(LEN=1)  :: str1
  CHARACTER(LEN=12) :: STRING

  REAL(KIND=rlg)    :: tool_datosec,date_start,date_start_ibm,date_end
  INTEGER, parameter, dimension(1) :: putseed=1
  INTEGER ::  jj,mm,aaaa,hh,minu,sec,jjs,compt
  REAL(KIND=rsh), allocatable,dimension(:)     :: para_real

 !     ..... Control for Time Integration .....
  character(19)     :: Cstart = '1990/11/01 00:00:00'  ! Starting date\\
  character(19)     :: Cend   = '1991/05/20 00:00:00'  ! Ending date\\
  character(19)     :: Cstep  = '0000/00/00 01:00:00'  ! Time step\\
  character(19)     :: Cmon   = '0000/00/01 00:00:00'  ! Monitor Interval\\
  character(19)     :: CTime,iarTime,jjdayTime,dtTime,StepTime 
  real(8)           :: dt, TTime, Tbefore, Season, Tmon,zero 
  integer           :: Iyr, Imon, Iday, Ihour, Imin, Isec 
 ! real(8)           :: cd2tt, nd2tt 
 ! character(19)     :: tt2cd 
  real(8)           :: Td, GraF, a, b, c 
  integer           :: lt=0 
  integer           :: nt,jjday,jjjday,iar,nouts
  real(8)           :: z45(365),z410(365),z420(365),z425(365),z435(365)
  real(8)           :: z55(365),z510(365),z520(365),z525(365),z535(365)
  real(8)           :: z65(365),z610(365),z620(365),z625(365),z635(365)
  real(8)           :: temp5(365),temp10(365),temp20(365),temp25(365),temp35(365)
  real(8),parameter :: d2s= 86400.0d0    ! day -> sec
  integer           :: hourday,init
  integer, parameter :: dd=200
  real(4)           :: z1(dd),z2(dd),z3(dd),tem(dd)   
!------------------------------------------------------
! executable statement
!---------------------------------------------------------

  allocate(para_real(npara))
  para_real=p

  !
  !     ***** Initial Setting *****
  !     ..... for time control .....
  TTime = cd2tt(Cstart)                                ! Starting Date
! kanei thn metatropi ths hmeromhnias Cstart pou einai se character form
! se pragmatiko arithmo.
!  print *, 'ttime= ', ttime                   ! 48902400.0000000
  CTime = TT2CD(TTime)                                 ! present time (character form)
!  print *, 'ctime = ', ctime                  ! '0001/07/20 00:00:00'
!prokyptei pali to Cstart

  dt    = cd2tt(Cstep) - cd2tt('0000/00/00 00:00:00')   !Time Step of integration (real8)
  Tmon  = cd2tt(Cmon)  - cd2tt('0000/00/00 00:00:00')  ! Monitor Interval (real8 form)
  dtTime=TT2CD(dt)
!  print *, 'dtTime = ', dtTime    ! To dtTime=0000/01/01 01:00:00 diaforo apo to Cstep

   zero=cd2tt('0000/00/00 00:00:00')
!  print *, 'dt = ', dt     !dt=3600
!  print *, 'zero = ', zero !zero=28771200.00000000


  nt = NINT( ( cd2tt(Cend) - cd2tt(Cstart) ) / dt )    ! Total Time Steps
!  print *, 'time steps (nt) = ', nt    !87672
!  pause

if (optimisation) then
     cost=0.0_rsh
     m_size(:)=0.0_rsh
     !   Transformation inverse des parametres
     do i=1,npara 
        if (p(i)>0.0_rsh) then 
           para_real(i)=(p(i)*h_bound(i)+para_ref(i))/(p(i)+1)
        else
           para_real(i)=(p(i)*l_bound(i)-para_ref(i))/(p(i)-1)
        endif
     enddo
  endif

 !  write(*,*) 'Nombre de stations', nstation
 !   do st=1,nstation
 !   write(*,*)'station', station(st)
OPEN(UNIT=112,FILE='environment.dat',STATUS='unknown')
  !-----read the 3 zoop groups and temp from 3d model
    do jjday=1,dd
     READ(112,1998) nouts,z1(jjday),z2(jjday),z3(jjday),tem(jjday)
 1998 format(i3,f7.4,f7.4,f7.4,f7.4)
!      Print*, 'Field', jjday, z1(jjday), z2(jjday), z3(jjday),tem(jjday)
    end do


     if (optimisation) then
        CALL random_seed(put=putseed) ! to make sure we have same simulation and fcost for a common parameter of the loop, and allow correct optimisation
     endif

     init=1
     compt=0
     nbouc=0
    !     write(*,*) 'iar',nt, dt

! MAIN LOOP
  do lt = 1, nt
     Tbefore = TTime                                 ! one step before present time
     TTime   = TTime + dt                            ! present time (real8 form)
     CTime = TT2CD(TTime)                            ! present time (charactor form)
     CALL TT2ND(Iyr, Imon, Iday ,Ihour, Imin, Isec ,TTime)

     hourday=(TTime - ND2TT(Iyr,Imon,Iday,1,0,0) )/ 3600

     jjday = 1 + (TTime - ND2TT(1990,1,1,0,0,0)) / d2s
     iar = mod(jjday,365)+1

     jjjday = 1 + ( TTime - ND2TT(1990,11,1,0,0,0) )/ d2s

!     write(*,*)'Optime', CTime,jjjday,z1(jjjday),tem(jjjday)
  
    call Anchovy(para_real,init,TTime,Tbefore,z1(jjjday),z2(jjjday),z3(jjjday),tem(jjjday) ) 
!   Save

   if (mod(hourday,24)==0 .and. xfi(1)/=0.0) then

   !     write(*,*)'Monitor1', CTime,hourday,xfi(1),bl(1)
   !     write(*,*)'Monitor2', iar,jjday,hourday,nbouc
           
   !        if (compt>=timestep_out) then
              if (optimisation) then
                 m_size(nbouc+1)=bl(1)
              endif
              nbouc=nbouc+1
              compt=dt
       else
              compt=compt+dt
   endif ! end of daily monitoring of model
            
   enddo ! end of main loop
     
    close(112)

! Optimisation
  !   if (optimisation) then
        do n=1,1 !nb_patch
           CALL fcost(n) !station(st)
        enddo
   ! endif

!  enddo ! fin boucle sur stations

 if (optimisation) then
     cost=cost !/nstation
     write(str1,'(i1)') npara
     STRING='('//str1//'(F10.4,x))'
     write(23,STRING) para_real, cost
!     write(*,*) 'parameter set ', para_real
     write(*,*) 'FCost ',JJday, cost
  endif

mod1dv=cost

END FUNCTION


