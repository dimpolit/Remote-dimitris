MODULE DEB
  !!======================================================================
   !!                   ***  MODULE DEB  ***
   !!             DEB anchois (larves et adultes)
   !!======================================================================
   !! * Modules used

  USE comibm
  USE com, ONLY : rsh,rlg,lchain,dt

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC  deb_init, deb_cycle             ! routine called ibm

   !! * Shared module variables
   LOGICAL, PUBLIC :: debuse,read_environ,T_fix,F_Fix
   REAL(kind=rsh), PUBLIC :: ffix,tempfix
   REAL(kind=rsh), PUBLIC,allocatable, dimension(:) :: tempread,slope, intercept,biomzoo

   !! * Private variables
   !REAL(kind=rsh), allocatable, dimension(:) :: slope, intercept,biomzoo
   logical      :: season
   INTEGER      :: yearspawn


 !! Set anchovy parameters
  REAL(kind=rsh), parameter :: TA = 9800.0_rsh    ! K, Arrhenius temperature ; ref = after Regner 1996;
  REAL(kind=rsh), parameter :: T1 = 286.0_rsh       ! K, Reference temperature ; ; ref = ;
  REAL(kind=rsh), parameter :: TAL = 30000.0_rsh   ! K, Arrhenius temperature at low temp {Kooy2000}
  REAL(kind=rsh), parameter :: TL = 281.0_rsh          ! K, lower boundary temp range
 ! REAL(kind=rsh), parameter :: K = 1.0_rsh          ! mgC.m3, saturation constant 
  REAL(kind=rsh), parameter :: Hb=0.004_rsh         ! d.cm^2, maturity at birth (first feeding)
!  REAL(kind=rsh), parameter :: Hj = 1000.0_rsh         ! d.cm^2, maturity at metamorphosis
  REAL(kind=rsh), parameter :: Hp =8000.0_rsh        ! J, maturity at puberty
!# REAL(kind=rsh), parameter :: g = 6.0_rsh              ! energy investment ratio
!  REAL(kind=rsh), parameter :: vc = 0.033_rsh !       ! cm/j, energy conductance
!  REAL(kind=rsh), parameter :: kap = 0.65_rsh        ! -, Veer
!#  REAL(kind=rsh), parameter :: kM = 0.015_rsh        ! 1/d, somatic maintenance rate coeff
!#  REAL(kind=rsh), parameter :: kJ = kM              ! 1/d, maturity maintenance rate coeff
!  REAL(kind=rsh), parameter :: pAm = 22.0_rsh    ! mol/d.cm2, spec max assim rate {J_EAm} ! 
  REAL(kind=rsh), parameter :: shape = 0.200_rsh  ! -, shape coefficient 
  REAL(kind=rsh), parameter :: shapeb = 0.136_rsh    ! -, shape coefficient at birth
  REAL(kind=rsh), parameter :: Sizeb = 0.4_rsh ! -cm, Length at first feeding (birth in DEB)
  REAL(kind=rsh), parameter :: d_V  = 1.0_rsh          ! g/cm^3, specific density of structure
  REAL(kind=rsh), parameter :: rho_V = 8000.0_rsh   !  J/g, = mu_V / w_V; Energy density of structure
  REAL(kind=rsh), parameter :: rho_E = 8000.0_rsh  !  J/g, = mu_E / w_E; Energy density of reserve

! Reproduction
 REAL(kind=rsh), parameter :: TR = 286.0_rsh      ! 6, K, temperature at start spawning
!# REAL(kind=rsh), parameter :: kapR = 0.95_rsh      ! 14, -, Fraction of reproduction energy fixed in eggs; ref = ;
 REAL(kind=rsh), parameter :: vB = 0.1_rsh             ! 17, cm/d, batch speed (low value gives large batches)
!# REAL(kind=rsh), parameter :: spawn = 3.0_rsh/6.0_rsh  ! 29, -, (approximate) length of spawning season as fraction of year
 REAL(kind=rsh), parameter :: Eb = 80.0_rsh             ! , J/cm3, specific batch energy 
 REAL(kind=rsh), parameter :: E0=0.15_rsh             ! J, Energy for one egg, after Heidi

 ! Ajout Martin
  REAL(kind=rsh), parameter :: sh_max =0.75_rsh     !  fraction, max shrinking rate from L_max before dying

 ! New after Romain
! vc ou Em suffisent
!#  REAL(kind=rsh), parameter :: Em=pAm/vc       ! J/cm3, max Reserve Energy (max. storage density);
!#  REAL(kind=rsh), parameter :: Em=666.0_rsh       ! J/cm3, max Reserve Energy (max. storage density);
! EG ou g suffisent
!#REAL(kind=rsh), parameter :: EG=g*kap*Em     ! J/cm3, volume-spec costs of structure
!!  REAL(kind=rsh), parameter :: EG=2600             ! J/cm3, volume-spec costs of structure, after Heidi
  REAL(kind=rsh), parameter :: pMi= 39.0_rsh       ! J/cm3/d, volume-spec maintenance costs,

  ! pour optimisation
    REAL(kind=rsh) :: kap, EG,vc,pAm,Hj,K
! Lm = k*pAmi/pMi		! cm,=vT/(g*kM), ultimate length;

CONTAINS

!!======================================================================

      subroutine deb_init(p)

  !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE deb_init  ***
   !&E
   !&E ** Purpose : initialize deb values for state variables
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : ibm
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  M. Huret (02-2011)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

        USE optim, only : optimisation


   !! * Arguments
        REAL(KIND=rsh), DIMENSION(:), INTENT(IN) :: p

   !! * Local declarations
        real(kind=rsh) ::  f_init, R_init, e_init
        INTEGER    :: i
        integer jj,mm,aaaa,hh,minu,sec
   !!----------------------------------------------------------------------
   !! * Executable part

! Initialisation if optimisation
   !     if (optimisation) then
           kap=p(1)
           EG=p(2)
           vc=p(3)
           pAm=p(4)
           Hj=p(5)
           K=p(6)
  !      endif

! Initialisation of embryo stage
!!----------------------------------------------------------------------

! food-----------------------
  f_init = ffix  
  part(:,:)%f=f_init
  if (read_environ) then
     CALL read_env
     part(:,:)%temp=tempread(1)
  endif

   if (T_fix) part(:,:)%temp=tempfix

  season=.FALSE.
  call tool_decompdate(date,jj,mm,aaaa,hh,minu,sec)
  yearspawn=aaaa+1
    
! to start for larvae
  part(:,:)%L = Sizeb*shapeb 
  part(:,:)%H = Hb
  part(:,:)%E = E0

! to start for adults
!  part(:,:)%L = part(:,:)%size*shape
!  part(:,:)%H = Hp

  part(:,:)%L_max=part(:,:)%L  ! sauvegarde de L_max atteint au cours du temps

!Reproduction
  R_init = 0.0_rsh
  part(:,:)%R = R_init
  part(:,:)%Gam=0.0_rsh
  part(:,:)%Nbatch=0.0_rsh

!Weight
  part(:,:)%Wdeb=0.0_rsh

  !         rB = kM * g/ (3 * (f + g))  !  von Bert growth rate
  !        Lm = v/ (kM * g)            !  max, ultimate length

END SUBROUTINE deb_init

!!======================================================================

      subroutine deb_cycle(n,m)
 !&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE deb_init  ***
   !&E
   !&E ** Purpose :  Modele D.E.B. generique adapté anchois 
   !&E ** Calculate variation in DEB compartments with ODE 
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : ibm
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  M. Huret (02-2011) from DEB-Larve L. Pecquerie
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used


   !! * Arguments
   INTEGER, INTENT( in )       :: n,m 

   !! * Local declarations
   real(kind=rsh) ::   cor_T, K_food,cor_shape,shape_fun,phys_L,tiers,Hfin,Hdeb
   real(kind=rsh) ::   pA,pC,pG,pJ,pR,pMT,pM,pM2,pM3,pR2,pGam
   real(kind=rsh) ::   kJT, kMT,pAmT,vT,f
   real(kind=rsh) ::   E,L,R,H,L_O,WV,WE,WR,WG,NRJ_V,NRJ_g,Ebatch
   real(kind=rsh) ::   T,X
   real(kind=rsh) ::   dR,dH,dE,dL,dGam
   integer        ::    jj,mm,aaaa,hh,minu,sec
   !!----------------------------------------------------------------------
   !! * Executable part

if (part(n,m)%number>0) then

! Affectation  
   E=part(n,m)%E
   L=part(n,m)%L 
   R=part(n,m)%R
   H=part(n,m)%H

! Temperature
!   if (read_environ) part(n,m)%temp=tempread(m) ! on force avec climato  
   if (read_environ) part(n,m)%temp=tempread(nbouc) ! on force avec climato  
   if (T_fix) part(n,m)%temp=tempfix

   ! Sinusoide de la temperature
    call tool_decompdate(date,jj,mm,aaaa,hh,minu,sec)
 !  part(n,m)%temp=15.0_rsh+5.0_rsh*sin(((mm-1)*30+jj+220.0_rsh)*2.0_rsh*pi/365.0_rsh)  ! entre 10 et 20 avec mini fin fevrier (150=270(dephasage pour avoir min de sinus a t=0) -70 (jour de premier rechauffement)

   T=part(n,m)%temp+273.0_rsh
   cor_T = exp(TA/ T1 - TA/ T) * &   
       (1 + exp(TAL/ T1 - TAL/ TL))/ (1.0_rsh + exp(TAL/ T - TAL/ TL)) ! more elaborated correction for lower temp limit

! Food
  X=0.0_rsh
  K_food = K !  + L * 10; % account for changes in preference
  f=ffix

if (read_environ .and. .not.F_Fix .and. part(n,m)%size>0.0_rsh) then
   X=get_Xdeb(n,m,part(n,m)%size)
   f = X/ (X + K_food)                                                         ! functional response    
endif


  part(n,m)%f=f
  part(n,m)%X=X


!-------------------------------------------------------------------------------
! shape correction function
cor_shape = (MAX( Hb, MIN(H, Hj))/ Hb)**(1.0_rsh/3.0_rsh)  ! cf p161 Laure ? cor_shape=1 pour H=Hb, donc pour jeune larve
!cor_shape=1.0_rsh

! CORRECTED PARAMETERS--------------------------------------
pAmT = pAm*cor_shape*cor_T	
vT = vc*cor_shape*cor_T
!kMT = kM*cor_T   
pMT=pMi*cor_T         
!kJT = kJ*cor_T	

! FLUXES------------------------------------------------------------------------
pA = pAmT*f*L**2                        ! Assimilation
!pM = kMT*EG*L**3		        ! J/d somatic maintenance;
pM = pMT*L**3	                        ! J/d somatic maintenance;

kMT=pMT/EG
kJT=kMT    ! choix ??

pC = (E/L**3)*((vT*EG*L**2)+pM)/(EG+(kap*E/L**3))	! J/d, Reserve Mobilisation, ref=Kooijman (Book:p37); (ok avec eq Laure 3.12)
pG =  (kap*pC)-pM						! J/d, Energy allocated to somatic growth;
pJ = kJT*H							! J/d, allocation to maturity maintenance;
pR = (1-kap)*pC-pJ                                	! J/d, allocation to maturity or reproduction;

! Reproduction and  Emergency maintenance
if (.not.season .and. T>TR .and. mm>3) then
   if (aaaa==yearspawn) then
      season=.TRUE.
      yearspawn=aaaa+1
   endif
endif

pR = (1-kap)*pC-pJ                                	! J/d, allocation to maturity or reproduction;

if (H>=Hp) then
   pR2=R*(vT/L + pMT/EG)*(1-kap*E/(EG*L**3+kap*E))	! J/d, potential allocation to gametogenesis (Pethybridge13 or Bernard11)
   pR2=min(pR2,Eb*L**3.0_rsh/3.0_rsh)                   ! on fixe une limite pour l'energie alloue à Ebatch / 3jours. Faire mieux.
   pM2=0.0_rsh
   pM3=0.0_rsh
   pGam=0.0_rsh
   if (pG<0.0_rsh) then                                                       !Bad condition, need to remobilise energy 
      pM2=min(-pG,pR2)	                                                  ! J/d, first mobilization from R (same speed as for gametogenesis)
      if (season) pM3=max((-pG-pM2),0.0_rsh)	          ! J/d, second mobilization from gametes (/EG cf Heidi ???) 
   endif
   if (season) pGam=pR2-pM2
endif

!ODE---------------------------------------------------------------------
dE = pA-pC				          ! J/d, Dynamics of Reserve Energy;

if (pG>=0.0_rsh) then
   dL = pG/(3.0_rsh*EG*L**2)             ! cm/d, growth increment;
else                                                        ! Test for mortality 
   if (pG+pM2+pM3<0.0_rsh) part(n,m)%number=0                ! maintenance pas assuree then mortality
endif

if (H<Hp) then
   dH = pR				   ! J/d, variation in level of maturity;
   dR = 0.0_rsh			                   ! J/d, variation in reproduction buffer energy;
else
   Ebatch=0.0_rsh
   dH = 0.0_rsh
   if (.not.season) then 
      dR=pR-pM2 ! only necessary fraction is mobilized for emergency
   else
      dR=pR-pR2
   endif
   dGam=pGam-pM3
endif
 
part(n,m)%E=E+dE*dt/86400.0_rsh
part(n,m)%L=L+dL*dt/86400.0_rsh
part(n,m)%H=H+dH*dt/86400.0_rsh
part(n,m)%R=R+dR*dt/86400.0_rsh

if (H>=Hp .and. season) then
   part(n,m)%Gam=part(n,m)%Gam+dGam*dt/86400.0_rsh
   if (part(n,m)%Gam>(2.0_rsh*Eb*L**3)) then        ! enough energy to spawn 1 batch and same amount remains for vitellogenesis (Somarakis Reproduce Del.)
      Ebatch=Eb*L**3
      part(n,m)%Nbatch=part(n,m)%Nbatch+1
   endif 
   if (part(n,m)%R<=0 .or. mm>8) then     ! last batch
      Ebatch=part(n,m)%Gam
      part(n,m)%Nbatch=part(n,m)%Nbatch+1
   endif
   part(n,m)%Gam=part(n,m)%Gam-Ebatch
endif

!----------------------------------------
! Test condition for reproductive season
!if (part(n,m)%R<=0) season=.FALSE. ! Ne marche pas pour l'instant 
if (mm>8 ) season=.FALSE.

!---------------------------------------------------------
! Shrinking and test for mortality
!if (part(n,m)%L > part(n,m)%L_max) part(n,m)%L_max=part(n,m)%L
!if (part(n,m)%L < sh_max*part(n,m)%L_max) part(n,m)%number=0

!---------------------------------------------------------
! Calcul du shape si variable en fonction du stade de developement
if (H<Hj) shape_fun=shapeb
if (H>=Hj) shape_fun=shape
!tiers=1.0_rsh/3.0_rsh
!if (H<Hj) shape_fun=((H**tiers-Hb**tiers)*shape+(Hj**tiers-H**tiers)*shapeb)/(Hj**tiers-Hb**tiers) !(Martin)

!Hdeb=0.01_rsh
!Hfin=20.0_rsh
!if (H<=Hdeb) shape_fun=shapeb
!if (H>Hdeb .and. H<=Hfin) shape_fun=((H**tiers-Hdeb**tiers)*shape+(Hfin**tiers-H**tiers)*shapeb)/(Hfin**tiers-Hdeb**tiers) !(Martin) en considerant debut metamorphose pour Hj quand ralentissement croissance, et Hfin pour reprise croissance avec Von Bertalanffy et shape adulte) Acquisition du shape final en cours de metamorphose
!if (H>=Hfin) shape_fun=shape

phys_L = part(n,m)%L/shape_fun !  in cm
part(n,m)%size=phys_L

! Weight
  WV = d_V*part(n,m)%L**3 
  WE = part(n,m)%E/rho_E ! W of E
  WR = part(n,m)%R/rho_E !  W of E_R
  WG= part(n,m)%Gam/rho_E !  W of Gam
  part(n,m)%Wdeb = WV+WE+WR+WG

! Energy
  NRJ_V = d_V*rho_V*part(n,m)%L**3
  NRJ_g = (NRJ_V + part(n,m)%E + part(n,m)%R+part(n,m)%Gam) !/ part(n,m)%Wdeb 

endif ! test on number of living individuals in Super Individual

END SUBROUTINE deb_cycle


!!======================================================================

      subroutine read_env

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE read_end ***
   !&E
   !&E ** Purpose :  read the slope and intercept of NBSS as well as temp 
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by : DEB_init
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

     !! * Local declarations
     INTEGER,    PARAMETER   :: nline=3000
     INTEGER                            :: nlinef,i
     REAL(kind=rsh) :: fake
     CHARACTER(len=50) :: fileopen,varname

  !---------------------------------------------------------------------
  
   fileopen='../input/trajenv2.txt'
     open (95, file=fileopen, status='old',form='formatted')
     nlinef=0
     read(95,'(a)') varname 
     do i=1,nline
         read(95,*,end=990) fake
         nlinef=nlinef+1
     enddo
990 continue
     close(95)
     open (95, file=fileopen, status='old',form='formatted')
     allocate(slope(nlinef),intercept(nlinef),tempread(nlinef),biomzoo(nlinef))
     read(95,'(a)') varname 
     do i=1,nlinef
         read(95,*) fake,fake,fake,biomzoo(i),slope(i),intercept(i),tempread(i)
         nlinef=0
     enddo
     close(95)

END SUBROUTINE read_env

!!======================================================================
   function get_Xdeb(n,m,size)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE  get_Xdeb ***
   !&E
   !&E ** Purpose : calculate f (standardized food response) for the DEB
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

 !! * Arguments
     INTEGER    :: n,m 
      REAL(kind=rsh) , INTENT (in) :: size
      REAL(kind=rsh) :: get_Xdeb

     !! * Local declarations
      real(kind=rsh) ::  slope_max_preywidth,slope_min_preywidth,maxwidth,esdmax,esdmin,minwidth
      real(kind=rsh) ::  biom,sl,inter,newinter

  !---------------------------------------------------------------------

  slope_max_preywidth=0.035_rsh    ! Delivrable FACTS
  slope_min_preywidth=0.005 ! Test 
  maxwidth=slope_max_preywidth*part(n,m)%size*10.0_rsh ! cm -> mm
  maxwidth=min(maxwidth,3.0_rsh) ! arbitraire ! pour cohérence à maxzoowidth issu du modele ci-dessous
  esdmax=size2ESD(maxwidth)
  minwidth=0.04_rsh+slope_min_preywidth*(part(n,m)%size*10.0_rsh-4.0_rsh)
  minwidth=min(minwidth,0.2_rsh) !toute gamme mesozoo 
  esdmin=size2ESD(minwidth)

!  slope(m)=-1.0_rsh
!  intercept(m)=0.0_rsh
!  biom=NBSS_2biom(slope(m),intercept(m),esdmin,esdmax) ! Martin pour utiliser climato LOPC (plusieurs particules sur chacune une maille de la grille
!  biom=NBSS_2biom(slope(nbouc),intercept(nbouc),esdmin,esdmax) ! Martin pour utiliser climato LOPC (à partir de trajectoire temporelle)
  newinter=biom_sl2inter(slope(nbouc),biomzoo(nbouc))                                            ! on recalcule l'intercept du NBSS à partir de biomasse modele
  biom=NBSS_2biom(slope(nbouc),newinter,esdmin,esdmax)                        ! on recalcule l'abondance dispo en fonction de la taille  
 

   get_Xdeb=biom
  

END function get_Xdeb

!!======================================================================

      function biom_sl2inter(sl,biom)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE  biom_sl2inter ***
   !&E
   !&E ** Purpose : calculate the intercept from NBSS slope and model biomass
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

 !! * Arguments
      REAL(kind=rsh) , INTENT (in) :: sl,biom
      REAL(kind=rsh) ::  biom_sl2inter

     !! * Local declarations
      REAL(kind=rsh)  :: minesd,maxesd,minwidth,maxwidth,slin,xx1,xx2

  !---------------------------------------------------------------------
  
      minwidth=0.02_rsh  !  min width zoo dans modele en mm
      maxwidth=3.0_rsh
      minesd=size2ESD(minwidth) ! transfo en ESD
      maxesd=size2ESD(maxwidth)
      slin=sl
      if (slin==-1) slin= -1.01

      xx1=ESD2weight(minesd)
      xx2=ESD2weight(maxesd)

      biom_sl2inter=log((biom*(slin+1))/(xx2**(slin+1)-xx1**(slin+1)))


    END function biom_sl2inter

!!======================================================================

      function NBSS_2biom(sl,int,size1,size2)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE  NBSS_2biom ***
   !&E
   !&E ** Purpose : calculate the biomasse within a size range given a NBSS (slope + intercept)
   !&E
   !&E ** Description :
   !&E
   !&E ** Called by :
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

 !! * Arguments
      REAL(kind=rsh) , INTENT (in) :: sl,int,size1,size2
      REAL(kind=rsh) :: NBSS_2biom

     !! * Local declarations
      REAL(kind=rsh)  :: S,DW,carb,xmid,xx1,xx2
      REAL(kind=rsh)  :: slin, interin,biom

  !---------------------------------------------------------------------
  
slin=sl
interin=int
if (slin==-1) slin= -1.01

interin = exp(int)  ! log(y) = sl * log(x) + int     eq.    y = exp(int) * x^sl

! convert x1 and x2 from mm esd to mgC:
  S=pi*(size1/2.0_rsh)**2.0_rsh
  DW=43.38_rsh*S**1.5_rsh   !from Lehette & Hernandez-Leon 2009 (general mesozooplankton)
  carb=0.447_rsh*DW     ! from Mauchline 1998 (average carbon content of copepods)
  xx1=carb
  S=pi*(size2/2.0_rsh)**2.0_rsh
  DW=43.38_rsh*S**1.5_rsh   !from Lehette & Hernandez-Leon 2009 (general mesozooplankton)
  carb=0.447_rsh*DW     ! from Mauchline 1998 (average carbon content of copepods)
  xx2=carb


!compute xmid in mugC (mean so as to have biom[x1,xmid]=biom[xmid,x2]

!  xmid = ((xx1_2[2]^(sl+1) + xx1_2[1]^(sl+1))/2)^(1/(sl+1));

! compute biom in mgC.m^{-3} and ab in #.m^{-3}
  biom = (interin*xx2**(slin+1))/(slin+1) - (interin*xx1**(slin+1))/(slin+1)
 ! ab = biom / (xmid/1000)

  NBSS_2biom=biom

END function NBSS_2biom

!=========================================================================

      function size2ESD(width)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE size2esd ***
   !&E
   !&E ** Purpose : calculate the esd from width
   !&E
   !&E ** Description :
   !&E!      ESD=2*sqrt(ab) ou ESD=sqrt(WL), Herman92 avec w et l largeur et longueur. Pour a=3b, ESD=sqrt(3*w²)=w*sqrt(3)
   !&E ** Called by :
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

 !! * Arguments
      REAL(kind=rsh) , INTENT (in) :: width
      REAL(kind=rsh)  :: size2ESD

     !! * Local declarations
    
  !---------------------------------------------------------------------

      size2ESD=width*3.0_rsh**0.5_rsh

END function size2ESD

!=========================================================================

      function ESD2weight(esd)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE size2esd ***
   !&E
   !&E ** Purpose : calculate Carbon weight from ESD
   !&E
   !&E ** Description :
   !&E!   ! from Lehette & Hernandez-Leon 2009 (DW general mesozooplankton)
   !&E    ! from Mauchline 1998 (average carbon content of copepods)
   !&E ** Called by :
   !&E
   !&E ** History :
   !&E       !  M. Huret (07-2012)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

 !! * Arguments
      REAL(kind=rsh) , INTENT (in) :: esd
      REAL(kind=rsh)  :: ESD2weight

     !! * Local declarations
      REAL(kind=rsh)  :: S,DW,carb

  !---------------------------------------------------------------------
  S=pi*(esd/2.0_rsh)**2.0_rsh
  DW=43.38_rsh*S**1.5_rsh   !from Lehette & Hernandez-Leon 2009 (general mesozooplankton)
  ESD2weight=0.447_rsh*DW   !from Mauchline 1998 (average carbon content of copepods)

END function ESD2weight

!=========================================================================


END MODULE DEB
