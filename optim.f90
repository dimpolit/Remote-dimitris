MODULE f1dim_mod  !Used for communication from linmin to f1dim.
USE nrtype
INTEGER(I4B) :: ncom
REAL(SP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS

FUNCTION f1dim(x)
IMPLICIT NONE
REAL(SP), INTENT(IN) :: x
REAL(SP) :: f1dim
!Used by linmin as the one-dimensional function passed to mnbrak and brent.

INTERFACE
FUNCTION mod1dv(x)
USE nrtype
REAL(SP), DIMENSION(:), INTENT(IN) :: x
REAL(SP) :: mod1dv
END FUNCTION mod1dv
END INTERFACE

REAL(SP), DIMENSION(:), ALLOCATABLE :: xt
allocate(xt(ncom))
xt(:)=pcom(:)+x*xicom(:)
f1dim=mod1dv(xt)
deallocate(xt)
END FUNCTION f1dim

END MODULE f1dim_mod
!===============================

MODULE optim

USE com

IMPLICIT none

  PRIVATE

  PUBLIC  do_optim

 !! * Shared module variables
  REAL(KIND=rsh), allocatable,dimension(:),PUBLIC     :: para_ref,h_bound,l_bound,para_init
  CHARACTER(LEN=100), PUBLIC                                :: file_optim,met_optim
  LOGICAL,              PUBLIC                                           :: optimisation
  INTEGER, PUBLIC                                                         :: npara
  ! Amoeba
  REAL(KIND=rsh), allocatable,dimension(:),PUBLIC       :: initval
  REAL(KIND=rsh), allocatable,dimension(:,:),PUBLIC     :: deltap,initvert
  ! Powell
  REAL(kind=rsh), allocatable,dimension(:,:), PUBLIC       :: ni

CONTAINS

!=============================================
SUBROUTINE  do_optim(func,p)

!&E---------------------------------------------------------------------
   !&E                 ***  ROUTINE set_optim  ***
   !&E
   !&E ** Purpose : set parameters of the optimisation and run it
   !&E
   !&E ** Description : 
   !&E
   !&E ** Called by : main
   !&E
   !&E ** External calls :
   !&E
   !&E ** Reference :
   !&E
   !&E ** History :
   !&E       !  M. Huret (10-2013)
   !&E
   !&E---------------------------------------------------------------------
   !! * Modules used

  USE fcost_mod, only : final_cost

   !! * Arguments

  REAL(KIND=rsh), DIMENSION(:), INTENT(IN) :: p

  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE

   !! * Local declarations
  REAL(KIND=rsh), allocatable,dimension(:) :: para_stand,para_optim
  INTEGER                 :: i,j,niter
  REAL(kind=rsh)       :: ftole

   !!----------------------------------------------------------------------
   !! * Executable part
  open (23, file=file_optim, status='new')
  allocate(para_stand(npara),para_optim(npara))

! transformation des para pour bounding
   do i=1,npara 
      if (p(i)>para_ref(i)) para_stand(i)=(p(i)-para_ref(i))/(h_bound(i)-p(i))
      if (p(i)<=para_ref(i)) para_stand(i)=(p(i)-para_ref(i))/(p(i)-l_bound(i))
   enddo

   ftole=0.000001_rlg
   final_cost=0.0_rsh
   niter=0

   if (met_optim=='powell') then
      allocate(ni(npara,npara))
      ni(:,:)=100.0_rsh 
!      do i=1,npara 
!         ni(:,i)=para_ref(:)/5
!      enddo
      CALL powell(para_stand,ni,ftole,niter,final_cost)
   endif

   if (met_optim=='amoeba') then
      allocate(initval(npara+1),initvert(npara+1,npara),deltap(npara+1,npara))
      deltap(:,:)=0.0_rsh
      do i=1,npara+1  ! n+1 vertices (or summit) of the simplex
         do j=1,npara
            if (i==j+1) then 
               deltap(i,j)=para_stand(j)/2  ! first vertices is from para_ref standardized, 
                                                          !then delta is applied successively to parameters one by one
        !       deltap(i,j)=1.0_rsh
            endif
            initvert(i,j)=para_stand(j)+deltap(i,j)
         enddo
         initval(i)=func(initvert(i,:))
      enddo
       !  write(*,*)'optimi', initvert(:,1)
      CALL amoeba(initvert,initval,ftole,func,niter)
      final_cost=initval(1)
      para_stand=initvert(1,:)
   endif

   ! transformation inverse
   do i=1,npara 
      if (para_stand(i)>0.0_rsh) para_optim(i)=(para_stand(i)*h_bound(i)+para_ref(i))/(para_stand(i)+1)
      if (para_stand(i)<=0.0_rsh)  para_optim(i)=(para_stand(i)*l_bound(i)-para_ref(i))/(para_stand(i)-1)
   enddo

     write(*,*) '------------------'
     write(*,*) 'FIN'
     write(*,*) 'Nombre d iterations:', niter
     write(*,*) 'Final cost:', final_cost
     write(*,*) 'Jeu de paramètres:', para_optim
     write(*,*) '------------------'
     
END SUBROUTINE do_optim

!=============================================

SUBROUTINE  powell(p,xi,ftol,iter,fret)

  USE nrtype; USE nrutil, ONLY :assert_eq, nrerror

  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: p
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: xi
  INTEGER(I4B), INTENT(OUT) :: iter
  REAL(SP), INTENT(IN) :: ftol
  REAL(SP), INTENT(OUT) :: fret
  
  INTERFACE
     FUNCTION mod1dv(p)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: p
       REAL(SP) :: mod1dv
     END FUNCTION mod1dv
  END INTERFACE
  
  INTEGER(I4B), PARAMETER :: ITMAX=200
  REAL(SP), PARAMETER :: TINY=1.0e-25_sp
  !Minimization of a function func of N variables. (func is not an argument, it is a fixed
  !function name.) Input consists of an initial starting point p, a vector of length N, an
  !initial N × N matrix xi whose columns contain the initial set of directions (usually the N
  !unit vectors), and ftol, the fractional tolerance in the function value such that failure to
  !decrease by more than this amount on one iteration signals doneness. On output, p is set
  !to the best point found, xi is the then-current direction set, fret is the returned function
  !value at p, and iter is the number of iterations taken. The routine linmin is used.
  !Parameters: Maximum allowed iterations, and a small number.
  INTEGER(I4B) :: i,ibig,n
  REAL(SP) :: del,fp,fptt,t
  REAL(SP), DIMENSION(size(p)) :: pt,ptt,xit

  n=assert_eq(size(p),size(xi,1),size(xi,2),'powell')
  fret=mod1dv(p)
  pt(:)=p(:)              ! Save the initial point.
  iter=0
  do
     iter=iter+1
     fp=fret
     ibig=0
     del=0.0            !Will be the biggest function decrease.
     do i=1,n           !Loop over all directions in the set.
        xit(:)=xi(:,i)  !Copy the direction,
        fptt=fret
        call linmin(p,xit,fret)                !minimize along it,
        if (fptt-fret > del) then              !and record it if it is the largest decrease so far
           del=fptt-fret
           ibig=i
        end if
     end do
     if (2.0_sp*(fp-fret) <= ftol*(abs(fp)+abs(fret))+TINY) RETURN
     !Termination criterion.
     if (iter == ITMAX) call nrerror('powell exceeding maximum iterations')
     ptt(:)=2.0_sp*p(:)-pt(:)             !Construct the extrapolated point and the average direction moved. Save the old starting point.
     xit(:)=p(:)-pt(:)
     pt(:)=p(:)
     fptt=mod1dv(ptt)                     !Function value at extrapolated point.
     if (fptt >= fp) cycle                !One reason not to use new direction.
     t=2.0_sp*(fp-2.0_sp*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
     if (t >= 0.0) cycle                  !Other reason not to use new direction.
     call linmin(p,xit,fret)              !Move to minimum of the new direction,
     xi(:,ibig)=xi(:,n)                   ! and save the new direction.
     xi(:,n)=xit(:)
  end do                                  !Back for another iteration.
END SUBROUTINE powell


!----------------------------------------------


SUBROUTINE linmin(p,xi,fret)

  USE nrtype; USE nrutil, ONLY :assert_eq
  USE f1dim_mod
  IMPLICIT NONE
  REAL(SP), INTENT(OUT) :: fret
  REAL(SP), DIMENSION(:), TARGET, INTENT(INOUT) :: p,xi
  REAL(SP), PARAMETER :: TOL=1.0e-4_sp
  !Given an N-dimensional point p and an N-dimensional direction xi, both vectors of length
  !N, moves and resets p to where the fixed-name function func takes on a minimum along
  !the direction xi from p, and replaces xi by the actual vector displacement that p was
  !moved. Also returns as fret the value of func at the returned location p. This is actually
  !all accomplished by calling the routines mnbrak and brent.
  !Parameter: Tolerance passed to brent.

  REAL(SP) :: ax,bx,fa,fb,fx,xmin,xx
  ncom=assert_eq(size(p),size(xi),'linmin')
  pcom=>p !Communicate the global variables to f1dim.
  xicom=>xi
  ax=0.0 !Initial guess for brackets.
  xx=1.0
  call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
  fret=brent(ax,xx,bx,f1dim,TOL,xmin)
  xi=xmin*xi !Construct the vector results to return.
  p=p+xi
END SUBROUTINE linmin

!--------------------------------------------------------------


SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)

  USE nrtype; USE nrutil, ONLY :swap
  IMPLICIT NONE
  REAL(SP), INTENT(INOUT) :: ax,bx
  REAL(SP), INTENT(OUT) :: cx,fa,fb,fc

  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE

  REAL(SP), PARAMETER :: GOLD=1.618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
  !Given a function func, and given distinct initial points ax and bx, this routine searches
  !in the downhill direction (defined by the function as evaluated at the initial points) and
  !returns new points ax, bx, cx that bracket a minimum of the function. Also returned are
  !the function values at the three points, fa, fb, and fc.
  !Parameters: GOLD is the default ratio by which successive intervals are magnified; GLIMIT
  !is the maximum magnification allowed for a parabolic-fit step.
  REAL(SP) :: fu,q,r,u,ulim
  
  fa=func(ax)
  fb=func(bx)
  if (fb > fa) then !Switch roles of a and b so that we can go downhill in the direction from a to b.
     call swap(ax,bx)
     call swap(fa,fb)
  end if
  cx=bx+GOLD*(bx-ax) !First guess for c.
  fc=func(cx)
  do !Do-while-loop: Keep returning here
     if (fb < fc) RETURN !until we bracket.
     !Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any possible
     !division by zero.
     r=(bx-ax)*(fb-fc)
     q=(bx-cx)*(fb-fa)
     u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
     ulim=bx+GLIMIT*(cx-bx)
     !We won’t go farther than this. Test various possibilities:
     if ((bx-u)*(u-cx) > 0.0) then !Parabolic u is between b and c: try it
        fu=func(u) 
        if (fu < fc) then !Got a minimum between b and c.
           ax=bx
           fa=fb
           bx=u
           fb=fu
           RETURN
        else if (fu > fb) then !Got a minimum between a and u.
           cx=u
           fc=fu
           RETURN
        end if
        u=cx+GOLD*(cx-bx) !Parabolic fit was no use. Use default magnification.
        fu=func(u) 
     else if ((cx-u)*(u-ulim) > 0.0) then !Parabolic fit is between c and its allowed limit.
        fu=func(u) 
        if (fu < fc) then
           bx=cx
           cx=u
           u=cx+GOLD*(cx-bx)
           call shft(fb,fc,fu,func(u))
        end if
     else if ((u-ulim)*(ulim-cx) >= 0.0) then !Limit parabolic u to maximum allowed value.
        u=ulim 
        fu=func(u)
     else !Reject parabolic u, use default magnification
        u=cx+GOLD*(cx-bx) 
        fu=func(u)
     end if
     call shft(ax,bx,cx,u)
     call shft(fa,fb,fc,fu) !Eliminate oldest point and continue.
  end do
CONTAINS
  
  SUBROUTINE shft(a,b,c,d)
    REAL(SP), INTENT(OUT) :: a
    REAL(SP), INTENT(INOUT) :: b,c
    REAL(SP), INTENT(IN) :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft

END SUBROUTINE mnbrak

!------------------------

FUNCTION brent(ax,bx,cx,func,tol,xmin)
  USE nrtype; USE nrutil, ONLY :nrerror
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: ax,bx,cx,tol
  REAL(SP), INTENT(OUT) :: xmin
  REAL(SP) :: brent
  
  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE

  INTEGER(I4B), PARAMETER :: ITMAX=100
  REAL(SP), PARAMETER :: CGOLD=0.3819660_sp,ZEPS=1.0e-3_sp*epsilon(ax)
  !Given a function func, and given a bracketing triplet of abscissas ax, bx, cx (such that bx
  !is between ax and cx, and func(bx) is less than both func(ax) and func(cx)), this
  !routine isolates the minimum to a fractional precision of about tol using Brent’s method.
  !The abscissa of the minimum is returned as xmin, and the minimum function value is
  !returned as brent, the returned function value.
  !Parameters: Maximum allowed number of iterations;g olden ratio;a nd a small number that
  !protects against trying to achieve fractional accuracy for a minimum that happens to be
  !exactly zero.
  INTEGER(I4B) :: iter
  REAL(SP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

  a=min(ax,cx) !a and b must be in ascending order, though the
  b=max(ax,cx) !input abscissas need not be.
  v=bx !Initializations...
  w=v
  x=v
  e=0.0 !This will be the distance moved on the step
  fx=func(x) !before last.
  fv=fx
  fw=fx
  do iter=1,ITMAX !Main program loop.
     xm=0.5_sp*(a+b)
     tol1=tol*abs(x)+ZEPS
     tol2=2.0_sp*tol1
     if (abs(x-xm) <= (tol2-0.5_sp*(b-a))) then !Test for done here.
        xmin=x !Arrive here ready to exit with best values.
        brent=fx
        RETURN
     end if
     if (abs(e) > tol1) then !Construct a trial parabolic fit.
        r=(x-w)*(fx-fv)
        q=(x-v)*(fx-fw)
        p=(x-v)*q-(x-w)*r
        q=2.0_sp*(q-r)
        if (q > 0.0) p=-p
        q=abs(q)
        etemp=e
        e=d
        if (abs(p) >= abs(0.5_sp*q*etemp) .or. &
             p <= q*(a-x) .or. p >= q*(b-x)) then
           !The above conditions determine the acceptability of the parabolic fit. Here it is
           !not o.k., so we take the golden section step into the larger of the two segments.
           e=merge(a-x,b-x, x >= xm )
           d=CGOLD*e
        else !Take the parabolic step.
           d=p/q
           u=x+d
           if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
        end if
     else !Take the golden section step into the larger
        e=merge(a-x,b-x, x >= xm ) !of the two segments.
        d=CGOLD*e
     end if
     u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
     !Arrive here with d computed either from parabolic fit, or else from golden section.
     fu=func(u)
     !This is the one function evaluation per iteration.
     if (fu <= fx) then !Now we have to decide what to do with our
        if (u >= x) then !function evaluation. Housekeeping follows:
           a=x
        else
           b=x
        end if
        call shft(v,w,x,u)
        call shft(fv,fw,fx,fu)
     else
        if (u < x) then
           a=u
        else
           b=u
        end if
        if (fu <= fw .or. w == x) then
           v=w
           fv=fw
           w=u
           fw=fu
        else if (fu <= fv .or. v == x .or. v == w) then
           v=u
           fv=fu
        end if
     end if
  end do !Done with housekeeping. Back for another iteration
  call nrerror('brent:exceed maximum iterations') 
CONTAINS
  SUBROUTINE shft(a,b,c,d)
    REAL(SP), INTENT(OUT) :: a
    REAL(SP), INTENT(INOUT) :: b,c
    REAL(SP), INTENT(IN) :: d
    a=b
    b=c
    c=d
  END SUBROUTINE shft
END FUNCTION brent

!===============================================================

SUBROUTINE amoeba(p,y,ftol,func,iter)

  USE nrtype; USE nrutil, ONLY :assert_eq, imaxloc,iminloc,nrerror,swap
  IMPLICIT NONE
  INTEGER(I4B), INTENT(OUT) :: iter
  REAL(SP), INTENT(IN) :: ftol
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: y
  REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: p

  INTERFACE
     FUNCTION func(x)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), DIMENSION(:), INTENT(IN) :: x
       REAL(SP) :: func
     END FUNCTION func
  END INTERFACE

  INTEGER(I4B), PARAMETER :: ITMAX=5000
  REAL(SP), PARAMETER :: TINY=1.0e-10
  !Minimization of the function func in N dimensions by the downhill simplex method of
  !Nelder and Mead. The (N + 1) × N matrix p is input. Its N + 1 rows are N-dimensional
  !vectors that are the vertices of the starting simplex. Also input is the vector y of length
  !N + 1, whose components must be preinitialized to the values of func evaluated at the
  !N + 1 vertices (rows) of p; and ftol the fractional convergence tolerance to be achieved
  !in the function value (n.b.!). On output, p and y will have been reset to N +1 new points
  !all within ftol of a minimum function value, and iter gives the number of function
  !evaluations taken.
  !Parameters: The maximum allowed number of function evaluations, and a small number.

  INTEGER(I4B) :: ihi,ndim                      !Global variables.
  REAL(SP), DIMENSION(size(p,2)) :: psum

  call amoeba_private

CONTAINS

  SUBROUTINE amoeba_private

    IMPLICIT NONE
    INTEGER(I4B) :: i,ilo,inhi
    REAL(SP) :: rtol,ysave,ytry,ytmp
    
    ndim=assert_eq(size(p,2),size(p,1)-1,size(y)-1,'amoeba')
    iter=0
    psum(:)=sum(p(:,:),dim=1)
    do                                    ! Iteration loop.
       ilo=iminloc(y(:))           !Determine which point is the highest (worst),next-highest,
       ihi=imaxloc(y(:))           !and lowest (best).
       ytmp=y(ihi)
       y(ihi)=y(ilo)
       inhi=imaxloc(y(:))
       y(ihi)=ytmp
       rtol=2.0_sp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
       !Compute the fractional range from highest to lowest and return if satisfactory.
       if (rtol < ftol) then              !If returning, put best point and value in slot 1.
          call swap(y(1),y(ilo))
          call swap(p(1,:),p(ilo,:))
          RETURN
       end if
       if (iter >= ITMAX) call nrerror('ITMAX exceeded in amoeba')
       !Begin a new iteration. First extrapolate by a factor −1 through the face of the simplex
       !across from the high point, i.e., reflect the simplex from the high point.
       ytry=amotry(-1.0_sp)
       iter=iter+1
       if (ytry <= y(ilo)) then !Gives a result better than the best point, sotry an additional extrapolation by a factor of 2.
          ytry=amotry(2.0_sp)
          iter=iter+1
       else if (ytry >= y(inhi)) then !The reflected point is worse than the second highest, so look for an intermediate
          ysave=y(ihi)                      ! lower point, i.e., do a one-dimensional contraction.
          ytry=amotry(0.5_sp)
          iter=iter+1
          if (ytry >= ysave) then
             !Can’t seem to get rid of that high point. Better contract around the lowest (best) point.
             p(:,:)=0.5_sp*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
             do i=1,ndim+1
                if (i /= ilo) y(i)=func(p(i,:))
             end do
             iter=iter+ndim                   !Keep track of function evaluations.
             psum(:)=sum(p(:,:),dim=1)
          end if
       end if
    end do                                        !Go back for the test of doneness and the next iteration
  END SUBROUTINE amoeba_private

  FUNCTION amotry(fac)
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: fac
    REAL(SP) :: amotry
    !Extrapolates by a factor fac through the face of the simplex across from the high point,
    !tries it, and replaces the high point if the new point is better.
    REAL(SP) :: fac1,fac2,ytry
    REAL(SP), DIMENSION(size(p,2)) :: ptry
    fac1=(1.0_sp-fac)/ndim
    fac2=fac1-fac
    ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
    ytry=func(ptry)                         !Evaluate the function at the trial point.
    if (ytry < y(ihi)) then                 !If it’s better than the highest, then replace
       y(ihi)=ytry                             ! the highest.
       psum(:)=psum(:)-p(ihi,:)+ptry(:)
       p(ihi,:)=ptry(:)
    end if
    amotry=ytry
  END FUNCTION amotry

END SUBROUTINE amoeba

END MODULE 

!===============================================================


