Program path
use Defs

implicit none

!  .. Subroutines ..
real (kind=8) :: dnrm2
external dnrm2
external init
external fullpath
external term
external ssolve
external sfun

!  .. Local declarations ..
integer, parameter :: fid = 10
integer :: np1, ifail, i
real (kind=8) :: Y(n+1), lambda(2), Z(n), S(n)

np1 = n+1

call init(lambda,Y(2:N+1))

if (ipar(1) <= 1) then ! Use homotopic method with given Lf
   call fullpath(LAMBDA,Y(2:NP1))
elseif (ipar(1) == 2) then ! Search optimal Lf and then use homotopic method
   call find(LAMBDA,Y(2:NP1))
elseif (ipar(1) == 3) then ! Integrate BVP for drawing
   call term(lambda,Y(2:N+1))
elseif (ipar(1) == 4) then ! Find solution for a given lambda with single shooting
   call ssolve(LAMBDA,Y(2:NP1),Z,ifail)
   call sfun(LAMBDA,Z,S)
   print *, "********************   AFTER SHOOT   *******************"
   print *, "   LAMBDA =", LAMBDA
   print *, "   Ysol   = ", Z/SCALE(FREE0(:))
   print *, "   Ftir(Ysol)   = ", S
   print *, "   || Ftir(Ysol) || = ", dnrm2(N,S,1)
   print *, "********************************************************"
   call term(lambda,Z(:))
end if;

stop
end program Path

!************************************************************
!************************************************************

subroutine Fullpath(lambda,ZI)

Use HOMPACK90, ONLY : fixpqf
Use Defs
implicit none

! Function (lambda,ZI) = fullpath(lambda,ZI), compute a solution
! at lambda = (1,1) (ie, for minimum consumption problem) from a
! solution ZI at input lambda (if lambda(1) > 1) or from nothing
! if lambda = (0,0)

!  .. Arguments ..
!  lambda    homotopic parameters
!  ZI        initial Shooting unknown

real (kind=8), dimension(2), intent(inout) :: lambda
real (kind=8), dimension(n), intent(inout) :: ZI

!  .. Subroutines ..
real (kind=8) :: dnrm2
external dnrm2
external homCI
external sfun
external ssolve
external term
external eltime
external rescale

!  .. Local declarations ..
integer, parameter :: fid = 20
integer :: iflag, nfe, np1, ifail, nfefake
real (kind=8) :: A(1), tci, tho, Y(n+1)
real (kind=8) :: arclen
real (kind=8), dimension(n) :: Z(n), S(n)

!  .. First excecutable statements ..
np1 = n+1

tci = 0.d0
tho = 0.d0

call eltime(tci)

call homCI(lambda,ZI,Z)

call sfun(lambda,Z,S)

if (dnrm2(n,S,1) > 1.d-3) then ! Initial conditions homotopy failed
   print *, 'Initial conditions homotopy failed with |Ftir| = ', dnrm2(n,S,1)
   print *, 'Cannot pursue'
   return
end if

call eltime(tci)

if (trace > 1) then
   print *, '******   Starting point of homotopy   ******'
   print *, '** Ysol = ', Z/SCALE(FREE0)
   print *, '** Sol = ', S
   print *, '||FSol|| = ', dnrm2(n,S,1)
   print *, 'Execution time (s) = ', tci
   print *, '********************************************'
else
   print *, 'Beginning of Differential continuation'
end if

call eltime(tho)

!
! Head of trace file
!

if (trace > 0) then
   WRITE (trace,*) '# path '
   WRITE (trace,*) ' # FREE0 = ', FREE0(:)
   WRITE (trace,'(" # Y0 = ",5(e24.16,1x))') Y0(:)/scale(:)
   WRITE (trace,'(" # T0 = ",(e24.16))') t0
   WRITE (trace,'(" # TF = ",(e24.16))') tmin
   WRITE (trace,'(" # PAR = ",5(e24.16,1x))') par(:)
   WRITE (trace,'(" # IPAR = ",5(i3,1x))') ipar(:)
   WRITE (trace,*) ' # NIT = ', nit
   WRITE (trace,*) ' # Data'
end if

!
! CALL TO HOMPACK ROUTINE
!

l_com = lambda

Y(1) = lambda(2)
Y(2:np1) = Z(:)
iflag = -2
nfe = 0
arclen = 0.d0

print *, ' Step No *** Lambda '
call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,trace,A,SSPAR,nfe,arclen)

do while (iflag == 3)
   lambda(2) = Y(1)
   call rescale(Y(2:NP1))
   call ssolve(lambda,Y(2:NP1),Z,ifail)
   if (trace > 1) then
      call sfun(lambda,Y(2:NP1),S)
      print *, "*************   intermediate result  ************************"
      print *, "   LAMBDA =", lambda
      print *, "   Ysol   = ", Z/SCALE(FREE0(:))
      print *, "   Ftir(Ysol)   = ", S
      print *, "   || Ftir(Ysol) || = ", dnrm2(n,S,1)
      print *, "*************************************************************"
   else
      print *, ' **** Refinement **** '
   end if
   print *, ' Step No *** Lambda '
   iflag = -2
   nfefake = 0
   call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,trace,A,&
        &      SSPAR(1:4),nfefake,arclen)
   nfe = nfe + nfefake
end do

lambda(2) = Y(1)

call sfun(lambda,Y(2:np1),S)

print *, "*******************   BEFORE REFINEMENT   &
     &************************"
print *, "   LAMBDA =", lambda
print *, "   Ysol   = ", Y(2:np1)/SCALE(FREE0(:))
print *, "   Ftir(Ysol)   = ", S
print *, "   || Ftir(Ysol) || = ", dnrm2(n,S,1)
print *, "******************************************&
     &************************"

call rescale(Y(2:NP1))

lambda(2) = 1.d0

! Refine solution
call ssolve(lambda,Y(2:NP1),Z,ifail)

if (trace > 0) then
   WRITE (trace,'(3X,e24.16,3X,(1X,e24.16))') lambda(2), Z/SCALE(FREE0)
end if

call sfun(lambda,Z,S)

print *, "*************************   AFTER REFINEMENT   &
     &************************"
print *, "   LAMBDA =", lambda
print *, "   Ysol   = ", Z/SCALE(FREE0(:))
print *, "   Ftir(Ysol)   = ", S
print *, "   || Ftir(Ysol) || = ", dnrm2(n,S,1)
print *, "***********************************************&
     &************************"

! Print solution in "out.dat" and "next.dat"
call term(lambda,Z)

call eltime(tho)

print *, '** Jacobian evaluations = ', nfe
print *, '** Differential homotopy execution time (secs) = ', tho
print *, '** Total execution time  (secs) = ', tci+tho

return
end subroutine Fullpath

!*************************************************************
!*************************************************************

Subroutine homCI(lambda,ZI,Z)

use Defs
use HOMPACK90, only : fixpqf
implicit none

! Function Z = homCI(lambda,ZI), compute a solution of the energy problem (ie,
! for lambda(1) = 1). The only restriction is that ZI(1) > 0. This solution
! is compute with discrete continuation on initial conditions (and differential
! homotopy if discrete continuation failed)

!  ... Arguments ...
real (kind=8), intent(inout) :: lambda(2)
real (kind=8), intent(in)    :: ZI(n)
real (kind=8), intent(out)   :: Z(n)

!  ... External subroutine ...
real (kind=8) :: dnrm2
external dnrm2
external rescale
external sfun
external ssolve

!  ... Local declarations ...
logical :: failed
integer :: ifail, nfe, iflag, np1
real (kind=8) :: Zaux(n), Y(n+1), S(n)
real (kind=8) :: arclen, step, last_step, lin, lout(2), A(1)

!  ... First executable statements ...

np1 = n+1
step = par(8)
failed = (step <= 0.d0)

if (lambda(1) > 0.d0) then ! Verifying if we have a zero
   call sfun(lambda,ZI,S)
   if (dnrm2(n,S,1)>1.d-2) then ! Not a zero but maybe a good estimation
      call ssolve(lambda,ZI,Z,ifail)
      call sfun(lambda,Z,S)
      if (dnrm2(N,S,1)>1.d-2) then
         print *, 'Bad initial point, restart at lambda(1) = 0'
         lambda(1) = 0.d0
         Z = 0.d0
         scale(free0(:)) = 1.d0
      else
         Z(:) = ZI(:)
      end if
   else
      Z(:) = ZI(:)
   endif
else
   Z(:) = ZI(:)
end if

lin = lambda(1)
lout = lambda
last_step = 0.d0

do while (.not.((lin == 1.d0).or.(failed)))
    lout(1) = lin + step
    
    call ssolve(lout,Z,Zaux,ifail)
    call sfun(lout,Zaux,S)
    if (trace > 1) then
       print *, 'lambda = ', lout(1)
       print *, '|Ftir| = ', dnrm2(n,S,1)
    end if
    
    if ((ifail == 1).and.(dnrm2(n,S,1)<1.d-3)) then
       if (last_step == step) then ! increase step
          step = min(2.d0*step,1.d0-lout(1),par(8))
       else
          last_step = step
       end if
       lin = lout(1)
       Z(:) = Zaux(:)
    else
       if (dnrm2(n,S,1)<1.d-5) then
          if (last_step == step) then ! increase step
             step = min(2.d0*step,1.d0-lout(1),par(8))
          else
             last_step = step
          end if
          lin = lout(1)
          Z(:) = Zaux(:)
          ifail = 1
       else ! Single shooting did not converge 
          step = step/2.d0
          lout(1) = lin+step
          failed = (step < par(9))
       end if
    end if
    
    if (ifail == 1) then
       call rescale(Z(:))
    end if

 end do

lambda(1) = lin

if (failed) then ! discrete continuation did not touch 1
   l_com = lambda
   Y(1) = lambda(1)
   Y(2:np1) = Z(:)
   iflag = -2
   A(:) = 0.d0
   nfe = 0
   arclen = 0.d0
   call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,0,A,SSPAR,nfe,arclen)
   do while (iflag == 3)
      lambda(1) = Y(1)
      call rescale(Y(2:NP1))
      call ssolve(lambda,Y(2:NP1),Z,ifail)
      call sfun(lambda,Y(2:NP1),S)
      print *, "*************   Refinement  *******************"
      iflag = -2
      call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,0,A,&
           &      sspar,nfe,arclen)
   end do
   call rescale(Y(2:NP1))
   call ssolve(lambda,Y(2:NP1),Z,ifail)
   call rescale(Z(:))
end if

return
end subroutine homCI

!************************************************************
!************************************************************

subroutine Find(lambda,ZI)

use Defs
implicit none

!  .. Arguments ..
!  lambda    homotopic parameters
!  ZI        initial Shooting unknown
real (kind=8), dimension(2), intent(inout) :: lambda
real (kind=8), dimension(n), intent(inout) :: ZI

!  .. Subroutines ..
real (kind=8) :: dnrm2
external dnmr2
external sfun
external solve
external eltime

!  .. Local declarations ..
logical :: start, ok
integer :: np1
real (kind=8) :: cl0, cl1, cl2, cl3, cl4, cl5, stepL, ctf
real (kind=8) :: mcl1, mcl2, mcl3, mcl4, mcl5
real (kind=8) :: t1
real (kind=8) :: S(n), tcl(5,2,N+4)

!  .. First executable Statements ..

t1 = 0.d0
call eltime(t1)

np1 = n+1
cl0 = par(3)

stepL = 1.d-1
cl1 = cl0-stepL
cl2 = cl0
cl3 = cl0+stepL

! Invalidation of data in tcl
tcl(:,:,1) = -1.d0

print *, '[SEARCH] Start search for ctf = ', par(1)
ok = .false.
! solve pb for cl1
call solve(ok,cl1,tcl,1)
if (tcl(1,1,1) > 0.d0) then
   mcl1 = tcl(1,1,N+4)
else
   print *, 'cannot solve for this cLf'
endif

! solve pb for cl2
call solve(ok,cl2,tcl,2)
if (tcl(2,1,1) > 0.d0) then
   mcl2 = tcl(2,1,N+4)
else
   print *, 'cannot solve for this cLf'
endif

! solve pb for cl3
call solve(ok,cl3,tcl,3)
if (tcl(3,1,1) > 0.d0) then
   mcl3 = tcl(3,1,N+4)
else
   print *, 'cannot solve for this cLf'
endif

start = .TRUE.

do while (start)
   
   if ((mcl2 >= mcl1).and.(mcl2 >= mcl3)) then
      start = .false.
   elseif ((mcl1 >= mcl2).and.(mcl2 >= mcl3)) then
      cl2 = cl1
      cl1 = cl1 - stepL
      mcl2 = mcl1
      cl3 = cl2
      mcl3 = mcl2
      tcl(3,:,:) = tcl(2,:,:)
      tcl(2,:,:) = tcl(1,:,:)
      call solve(ok,cl1,tcl,1)
      if (tcl(1,1,1) > 0.d0) then
         mcl1 = tcl(1,1,N+4)
      else
         print *, 'cannot solve for this ct'
      endif
   else
      cl2 = cl3
      cl3 = cl3 + stepL
      mcl2 = mcl3
      cl1 = cl2
      mcl1 = mcl2
      tcl(1,:,:) = tcl(2,:,:)
      tcl(2,:,:) = tcl(3,:,:)
      call solve(ok,cl3,tcl,3)
      if (tcl(3,1,1) > 0.d0) then
         mcl3 = tcl(3,1,N+4)
      else
         print *, 'cannot solve for this ct'
      endif
   endif
enddo

! On a maintenant un encadrement de la solution

print *, '  [SEARCH] We have a framing of the solution'

do while (cl3-cl1 > par(2))
   cl4 = (cl1+cl2)/2.d0
   cl5 = (cl2+cl3)/2.d0
   
   call solve(ok,cl4,tcl,4)
   
   call solve(ok,cl5,tcl,5)
   
   mcl4 = tcl(4,1,N+4)
   mcl5 = tcl(5,1,N+4)
   
   if ((mcl4 >= mcl1).and.(mcl4 >=mcl2)) then
      cl3 = cl2
      mcl3 = mcl2
      tcl(3,:,:) = tcl(2,:,:)
      cl2 = cl4
      mcl2 = mcl4
      tcl(2,:,:) = tcl(4,:,:)
   elseif ((mcl5 >= mcl2).and.(mcl5>=mcl3)) then
      cl1 = cl2
      mcl1 = mcl2
      tcl(1,:,:) = tcl(2,:,:)
      cl2 = cl5
      mcl2 = mcl5
      tcl(2,:,:) = tcl(5,:,:)
   else
      cl1 = cl4
      mcl1 = mcl4
      tcl(1,:,:) = tcl(4,:,:)
      cl3 = cl5
      mcl3 = mcl5
      tcl(3,:,:) = tcl(5,:,:)
   endif
   
enddo

ok = .true.
call solve(ok,cl2,tcl,2)

par(3) = cl2

mcl2 = tcl(2,2,N+4)

scale(free0(:)) = 1.d0
call term((/1.d0,1.d0/),tcl(2,2,4:N+3))

print *, '  [SEARCH] For ctf = ', par(1), ' ,best cLf = ', par(3)
print *, '  [SEARCH] Corresponding mass = ', mcl2

call eltime(t1)

print *, '*[SEARCH]* EXECUTION TIME (SECS) = ', t1

RETURN

END SUBROUTINE FIND

!************************************************************
!************************************************************

Subroutine solve(find,cLf,tcl,position)

use HOMPACK90, only : FIXPQF
use Defs
implicit none

logical, intent(IN)   :: find
real (kind=8), intent(IN)   :: cLf
real (kind=8), intent(INOUT):: tcl(5,2,N+4)
integer, intent(IN)   :: position

!  .. Subroutines ..
real (kind=8) :: dnrm2
external dnrm2
external Pfun
external Ssolve
external rescale
external HomCI

!  .. Local declarations ..
integer, parameter :: ndima = 1
real (kind=8) :: Z(N), YF(2*N)
real (kind=8) :: S(N), A(NDIMA), Y(N+1), arclen
real (kind=8) :: norm_S, lambda(2)
integer :: iflag, nfe, np1, np3, J, I, K, ifail

!  .. First executable statements

! Set shooting function parameters

NP1 = N+1
NP3 = N+3

par(3) = cLf

if (.not.find) then
   
   lambda(:) = 0.d0
   Z(:) = 0.d0
   scale(free0(:)) = 1.d0
   ! homotopy on Initial conditions
   call HomCI(lambda,Z,Y(2:NP1))
   lambda(1) = 1.d0

   call sfun(lambda,Y(2:NP1),S)
   norm_S = dnrm2(N,S,1)
   
   ! complete tcl
   if (norm_S < 1.d-2) then
      tcl(position,1,1) = 1.d0
      tcl(position,1,4:NP3) = Y(2:NP1)/scale(free0)
      call PFUN((/1.d0,0.d0/),Y(2:NP1),YF)
      tcl(position,1,N+4) = YF(N)/SCALE(N)
   else ! not a valid solution
      tcl(position,1,1) = -1.d0
   endif
   
   RETURN
   
else
   
   print *, '* Energy solution for ctf = ', par(1), '  and cLf = ', cLf
   print *, '** Ysol = ', tcl(position,1,4:NP3)
   
   !
   ! CALL TO HOMPACK ROUTINE
   !
   
   lambda(1) = 1.d0
   lambda(2) = 0.d0
   l_com = lambda
   Y(1) = 0.d0
   Y(2:NP1) = tcl(position,1,4:NP3)
   scale(free0(1:N)) = 1.d0
   call rescale(Y(2:NP1))

   iflag = -2
   nfe = 0
   arclen = 0.d0

   if (trace > 0) then
      WRITE (trace,*) '# path '
      WRITE (trace,*) ' # FREE0 = ', FREE0(:)
      WRITE (trace,'(" # Y0 = ",5(e12.6,1x))') Y0(:)/scale(:)
      WRITE (trace,'(" # T0 = ",e24.16)') t0
      WRITE (trace,'(" # TF = ",e24.16)') tmin
      WRITE (trace,'(" # PAR = ",5(e24.16,1x))') par(:)
      WRITE (trace,'(" # IPAR = ",5(i3,1x))') ipar(:)
      WRITE (trace,*) ' # NIT = ', nit
      WRITE (trace,*) ' # Data'
   end if

   print *, ' Step No *** Lambda '
   call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,trace,A,SSPAR,nfe,arclen)

   do while (iflag == 3)
      lambda(2) = Y(1)
      call rescale(Y(2:NP1))
      call ssolve(lambda,Y(2:NP1),Z,ifail)
      if (trace > 1) then
         call sfun(lambda,Y(2:NP1),S)
         print *, "*************   intermediate result  ************************"
         print *, "   LAMBDA =", lambda
         print *, "   Ysol   = ", Z/SCALE(FREE0(:))
         print *, "   Ftir(Ysol)   = ", S
         print *, "   || Ftir(Ysol) || = ", dnrm2(n,S,1)
         print *, "*************************************************************"
      else
         print *, ' **** Refinement **** '
      end if
      print *, ' Step No *** Lambda '
      iflag = -2
      call fixpqf(n,Y,iflag,arcre,arcae,ansre,ansae,trace,A,&
           &      SSPAR(1:4),nfe,arclen)
   end do
   
   lambda(2) = 1.d0
   
   call ssolve(lambda,Y(2:NP1),Z,ifail)
   
   call sfun(lambda,Z,S)  
   norm_S = dnrm2(N,S,1)
   
   ! complete tcl
   if (norm_S < 1.d-2) then
      tcl(position,2,1) = 1.d0
      tcl(position,2,2) = par(1)
      tcl(position,2,3) = cLf
      tcl(position,2,4:NP3) = Z/scale(free0(:))
      call PFUN((/1.d0,1.d0/),Z,YF)
      tcl(position,2,N+4) = YF(N)/SCALE(N)
   else ! not a valid solution
      tcl(position,2,1) = -1.d0
   endif
   
endif

return
end subroutine solve

!************************************************************
!************************************************************

subroutine rescale(Y)

use Defs
implicit none

! This subroutine compute scaling for Y according to scaltol, scalrange and scalpL

!  .. Arguments ..
REAL (KIND=8), INTENT (INOUT) :: Y(N)

!  .. Subroutines ..

!  .. Local declarations ..
INTEGER :: k, I

!  .. First executable statements

! Descaling
Y(:) = Y(:)/scale(free0(:))
DO I = 1,N
   if (abs(Y(i)) > scaltol) then
      k = int(log10(abs(Y(i))))+scalrange
      if (abs(Y(i)) < 1.d0) k = k-1
      scale(free0(i)) = 10.d0**(-k)
   else
      scale(free0(i)) = 1.d0
   endif
ENDDO

! Scaling
Y(:) = Y(:)*scale(free0(:))

RETURN

END SUBROUTINE RESCALE
