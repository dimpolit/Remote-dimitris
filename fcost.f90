MODULE fcost_mod

USE com
!USE comibm

IMPLICIT none

  PRIVATE

  PUBLIC  fcost

 !! * Shared module variables
  REAL(kind=rsh), PUBLIC  ::  cost
  REAL(kind=rsh), PUBLIC  :: final_cost
  REAL(kind=rsh), PUBLIC,dimension(280)  :: m_size,m_weight  ! attention durée en dur ici


  !! * Private variables
   REAL(kind=rsh), allocatable,dimension(:)  :: d_size,d_weight,d_age

CONTAINS
  !!======================================================================

   subroutine read_data  !(stat)

  !! * Modules used

  !! * Arguments
  !   CHARACTER(LEN=5) ,INTENT (in) :: stat

     !! * Local declarations
     
     CHARACTER(LEN=100)  :: file_data
     INTEGER,    PARAMETER   :: nline=15000
     REAL(kind=rsh), dimension(nline)  :: age,length
     character(LEN=180) :: rec1
     integer                     :: i,j
  !
  ! executable statement
  !---------------------------------------------------------

     j=1
!     file_data='./optimisation/size.dat'
     file_data='./optimisation/sardinedata.csv'
     open (22, file=file_data, status='old')
  !   read(22, '(a)')
     do i=1,nline
        read(22,'(a)',end=996) rec1
        read(rec1,*) age(j),length(j)
  !      write(*,*) 'openfile ',age(i),length(i)
        j=j+1
     enddo
996     continue
     allocate(d_size(j-1),d_age(j-1))
     d_size(:)=length(1:j-1)
     d_age(:)=age(1:j-1)
   
     close(22)

   END SUBROUTINE read_data

!!================================================================================
 subroutine fcost(n) !stat

!! * Modules used

 !! * Arguments
 !  CHARACTER(LEN=5) ,INTENT (in) :: stat
   INTEGER,INTENT(in) :: n

    !! * Local declarations
   INTEGER                                               ::   k,ndata,n1,n2
   REAL(kind=rsh) :: cost1,cost2
   INTEGER                                               ::  decalage
  !
  ! executable statement
  !---------------------------------------------------------

   decalage=0  ! dans data age=ouverture de la bouche (stries journalières) alors que modèle à partir de ponte
   n1=0
   n2=0
   cost1=0.0_rsh
   cost2=0.0_rsh
   CALL read_data !(stat)
   ndata=size(d_size)
  write(*,*)'Sizedata', ndata
   do k=1,ndata
      if (d_size(k) /= valmanq) then
         cost1=cost1+((m_size(int(d_age(k))+decalage)-d_size(k)))**2
         n1=n1+1
      endif
!      if (d_weight(k) /= valmanq .and. d_weight(k)<4.0_rsh) then
!         cost2=cost2+((m_weight(int(d_age(k))+decalage)-d_weight(k)) )**2
!         n2=n2+1
!      endif
   enddo
   cost=cost1/n1 !+cost2/n2
   
  ! write(*,*) m_size(int(d_age(1))+decalage),d_size(1),int(d_age(1))
  ! write(*,*) m_weight(int(d_age(1))+decalage),d_weight(1),int(d_age(1))


   deallocate(d_size,d_age)

 end subroutine fcost

!!================================================================================

END MODULE fcost_mod
