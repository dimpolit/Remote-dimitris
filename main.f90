program main

  USE com
  USE optim
  USE fcost_mod
  USE initsave

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

!------------------------------------------------------
! executable statement
!---------------------------------------------------------

!--------------------------------------------------------------
! Initialisation ----------------------------------------------
!----------------------------------------------------------------

  CALL namelist

!--------------------------------------------------------------
! optimisation ou run classique ----------------------------------------------
!----------------------------------------------------------------
  
  allocate(para_init(npara))
  para_init(1)=0.787_rsh  ! initialisation des paramètres (set to para_ref values if not twin experiment)
  para_init(2)=2148.0_rsh
  para_init(3)=0.836_rsh
  para_init(4)=2.735_rsh
  para_init(5)=827.0_rsh
  para_init(6)=0.66_rsh

  if (optimisation) then   ! boucle optimisation
     allocate(para_ref(npara),l_bound(npara),h_bound(npara))

     para_ref(1)=0.643_rsh  ! initialisation des paramètres
     para_ref(2)=2483.0_rsh
     para_ref(3)=0.033_rsh
     para_ref(4)=3.3_rsh
     para_ref(5)=832.0_rsh
     para_ref(6)=1.5_rsh

     l_bound(1)=0.0_rsh
     h_bound(1)=1.0_rsh
     l_bound(2)=2000.0_rsh
     h_bound(2)=3000.0_rsh
     l_bound(3)=0.0_rsh
     h_bound(3)=1.0_rsh
     l_bound(4)=0.0_rsh
     h_bound(4)=200.0_rsh
     l_bound(5)=0.0_rsh
     h_bound(5)=100000.0_rsh
     l_bound(6)=0.1_rsh
     h_bound(6)=100.0_rsh

     CALL do_optim(mod1dv,para_init)
  
  else ! run classique
     r=mod1dv(para_init)
  endif

END PROGRAM main

!=============================================
FUNCTION mod1dv(p)

  USE comibm
  USE ibm
  USE netcdf
  USE phys
  USE meteo
  USE initsave
  USE optim, only : optimisation, para_ref,h_bound,l_bound,npara
  USE fcost_mod
  USE DEB

 !! * Arguments
  REAL(KIND=rsh) :: mod1dv
  REAL(KIND=rsh), DIMENSION(:), INTENT(IN) :: p


 !! * Private variables
  INTEGER :: nb_iter,nb,st,n,status,try,i
  CHARACTER(LEN=19) :: tool_sectodat
 CHARACTER(LEN=1) :: str1
 CHARACTER(LEN=12) :: STRING

  REAL(KIND=rlg)    :: tool_datosec,date_start,date_start_ibm,date_end
  INTEGER, parameter, dimension(1) :: putseed=1
  INTEGER ::  jj,mm,aaaa,hh,minu,sec,jjs,compt
  REAL(KIND=rsh), allocatable,dimension(:)     :: para_real
    
!------------------------------------------------------
! executable statement
!---------------------------------------------------------

  allocate(para_real(npara))
  para_real=p

  if (optimisation) then
     cost=0.0_rsh
     m_size(:)=0.0_rsh
     m_weight(:)=0.0_rsh
     !   Transformation inverse des parametres
     do i=1,npara 
        if (p(i)>0.0_rsh) then 
           para_real(i)=(p(i)*h_bound(i)+para_ref(i))/(p(i)+1)
        else
           para_real(i)=(p(i)*l_bound(i)-para_ref(i))/(p(i)-1)
        endif
     enddo
  endif

  write(*,*) 'Nombre de stations', nstation


! boucle sur nombre de stations
  do st=1, nstation

   !  write(*,*) station(st)

     if (optimisation) then
        CALL random_seed(put=putseed) ! to make sure we have same simulation and fcost for a common parameter of the loop, and allow correct optimisation
     endif


! IF PROFIL FROM CTD
     if (ctd)  then
        CALL read_ctd(station(st),date_end)  ! overide h0
        date_start=date_end-duration*24.0_rlg*3600.0_rlg
     else
        date_start=tool_datosec(datec_start)
        date_end=date_start+duration*24.0_rlg*60.0_rlg*60.0_rlg ! en s
     endif

     !Date ibm
     date_start_ibm=date_end-duration_ibm*24.0_rlg*3600.0_rlg


     !Vertical layers
     kmax=int(h0/deltaz)

     !Time
     nb_iter=INT((date_end-date_start)/dt)  
     t=date_start
     compt=0
     nbouc=0
     date=tool_sectodat(t)    ! Hour of the day
     call tool_decompdate(date,jj,mm,aaaa,hh,minu,sec)
     jjs=jj

     call init_phy(kmax, tempctd, salctd,p)
     ! call init_bio(kmax) ! pour eulerien, not ready

     ! IBM
     CALL ibm_init(st,para_real) ! Aytes tis times pairnei h deb_init kai katopin h deb
  !  CALL ibm_init(st)

     !-------------------------------------------------------------------------
     ! Simulation-----------------------------------------------------
     !-------------------------------------------------------------------------
     do nb=1,nb_iter
        
        t=t+dt
        
        date=tool_sectodat(t) 

        if (var_phys) then
! Process
     !      call diffusion(tra,dz,dzw,kz,phis,phif,kmax,h,tetha) ! Eulerien
           call get_wind(p)
           call tool_decompdate(date,jj,mm,aaaa,hh,minu,sec)
           if (jj /= jjs) then
              jjs=jj
              call read_coef
              umaxmaree=(vbtmax*coefmaree)/70.0_rsh
              F0=umaxmaree*(f+omega2m)
           endif
           call calu(tsx)
           call calv(tsy)
           call turbk(u,v,buoy,h0,tsx,tsy,kmax)
           call turbpsi(u,v,buoy,h0,kmax)
           call turblepskznz(buoy,u,v,h0,kmax)
           
        endif    ! var_phys

        !IBM
        if (t >= date_start_ibm) then 
        ! Create file and write initial setup
           if (nbouc == 0) then
              if (optimisation) then
                 m_size(nbouc+1)=part(1,1)%size
                 m_weight(nbouc+1)=part(1,1)%Wdeb
              endif
              if (.not.optimisation) then
                 CALL traj_save(station(st))
              endif
              nbouc=nbouc+1
           endif
        !run IBM   
           call ibm_1d
        !Save
           if (compt>=timestep_out) then
              if (optimisation) then
                 m_size(nbouc+1)=part(1,1)%size
                 m_weight(nbouc+1)=part(1,1)%Wdeb
              endif
              if (.not.optimisation) then
                 write(*,*) 'date=',t,date,nb,nb_iter,station(st)
                 CALL traj_save(station(st))
              endif
              nbouc=nbouc+1
              compt=dt
           else
              compt=compt+dt
           endif
        endif
        
     enddo
     
! Optimisation
     if (optimisation) then
        do n=1,nb_patch
           CALL fcost(n,station(st))
        enddo
     endif

     ! FIN
     if (.not.optimisation) then
        do n=1, nb_patch
           status=NF90_CLOSE(ncid(n))           ! close: save new netCDF dataset
           deallocate(ncid)
        enddo
     endif

     deallocate(temp,u,v,sal,dz,dzw,dzu,sig,sigw,z,kz,psi,ect,eps,nz,buoy,part,nb_part,i_stage,i_density,size_cm,filetraj1dout,name_patch)
     if (var_phys) deallocate(latmet,lonmet,tabtime)
     if (read_environ) deallocate(slope,intercept,tempread,biomzoo)
   

  enddo ! fin boucle sur stations

 if (optimisation) then
     cost=cost/nstation
     write(str1,'(i1)') npara
     STRING='('//str1//'(F10.4,x))'
     write(23,STRING) para_real, cost
     write(*,*) 'parameter set', para_real
     write(*,*) 'FCost' , cost
  endif

mod1dv=cost

END FUNCTION
