!  This module provides global allocatable arrays used for the sparse
!  matrix data structures, and by the polynomial system solver.  The
!  MODULE HOMOTOPY uses this module.
!
      MODULE HOMPACK90_GLOBAL
      
      INTEGER, DIMENSION(:), ALLOCATABLE:: COLPOS, IPAR, ROWPOS
      real (kind=8), DIMENSION(:), ALLOCATABLE:: PAR, PP, QRSPARSE
      END MODULE HOMPACK90_GLOBAL


      MODULE HOMOTOPY       ! Interfaces for user written subroutines.
      USE HOMPACK90_GLOBAL
!
C     Interface for subroutine that evaluates RHO(A,LAMBDA,X) and returns it 
C in the vector V.
      INTERFACE
         SUBROUTINE RHO(A,LAMBDA,X,V)
         real (kind=8), INTENT(IN):: A(:),X(:)
         real (kind=8), INTENT(IN OUT):: LAMBDA
         real (kind=8), INTENT(OUT):: V(:)
         END SUBROUTINE RHO
      END INTERFACE
C Interface for subroutine that returns in the vector V the Kth column
C of the Jacobian matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at the
C point (LAMBDA, X).
      INTERFACE
         SUBROUTINE JACFUN(A,LAMBDA,X,Q)
         real (kind=8), INTENT(IN):: A(:),X(:)
         real (kind=8), INTENT(IN OUT):: LAMBDA
         real (kind=8), INTENT(OUT):: Q(:,:)
         END SUBROUTINE JACFUN
      END INTERFACE
      END MODULE HOMOTOPY


      MODULE HOMPACK90
!
!  This MODULE is an encapsulation of the HOMPACK90 drivers, and uses
!  the modules REAL_PRECISION (defines real precision for all
!  routines), HOMPACK90_GLOBAL (allocatable global data structures for
!  sparse matrices), and HOMOTOPY (interfaces to user written routines
!  defining the problem).
!
!  The intended usage is that the calling program would include a
!  statement like, for example,
!     USE HOMPACK90, ONLY : FIXPNF
!
      USE HOMPACK90_GLOBAL             ! Allocated data structures.
      CONTAINS
!


      SUBROUTINE FIXPQF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &     SSPAR,NFE,ARCLEN)
C
C Subroutine  FIXPQF  finds a fixed point or zero of the 
C N-dimensional vector function  F(X), or tracks a zero curve of a 
C general homotopy map  RHO(A,LAMBDA,X).  For the fixed point problem
C F(X) is assumed to be a C2 map of some ball into itself.  The 
C equation  X=F(X)  is solved by following the zero curve of the 
C homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A) ,
C
C starting from  LAMBDA = 0, X = A.   The curve is parameterized
C by arc length  S, and is followed by solving the ordinary 
C differential equation  D(HOMOTOPY MAP)/DS = 0  for  
C Y(S) = (LAMBDA(S), X(S)).  This is done by using a Hermite cubic 
C predictor and a corrector which returns to the zero curve in a 
C hyperplane perpendicular to the tangent to the zero curve at the 
C most recent point.
C
C For the zero finding problem  F(X)  is assumed to be a C2 map
C such that for some  R > 0,  X*F(X) >= 0  whenever  NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve of
C the homotopy map
C
C  LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from  LAMBDA = 0, X = A.
C
C A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to 
C be a C2 map from  E**M X [0,1) X E**N  into  E**N, which for 
C almost all parameter vectors  A  in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX] = N
C
C for all points  (LAMBDA,X)  such that  RHO(A,LAMBDA,X) = 0.  It is
C further assumed that
C
C         rank [ D RHO(A,0,X0)/DX ] = N.
C
C With  A  fixed, the zero curve of  RHO(A,LAMBDA,X)  emanating from
C LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by solving the 
C ordinary differential equation    D RHO(A,LAMBDA(S),X(S))/DS = 0
C for  Y(S) = (LAMBDA(S), X(S)), where  S  is arc length along the
C zero curve.  Also the homotopy map  RHO(A,LAMBDA,X)  is assumed to
C be constructed such that
C
C         D LAMBDA(0)/DS > 0.
C
C For the fixed point and zero finding problems, the user must supply
C a subroutine  F(X,V)  which evaluates  F(X)  at  X  and returns the
C vector F(X) in  V, and a subroutine  FJAC(X,V,K)  which returns in  V
C the Kth column of the Jacobian matrix of F(X) evaluated at X.  For
C the curve tracking problem, the user must supply a subroutine
C RHO(A,LAMBDA,X,V)  which evaluates the homotopy map  RHO at
C (A,LAMBDA,X)  and returns the vector  RHO(A,LAMBDA,X)  in  V, and
C a subroutine  RHOJAC(A,LAMBDA,X,V,K)  which returns in  V
C the Kth column of the  N X (N+1)  Jacobian matrix  
C [D RHO/D LAMBDA, D RHO/DX]  evaluated at  (A,LAMBDA,X).  FIXPQF
C directly or indirectly uses the subroutines  F (or RHO), 
C   FJAC (or RHOJAC),  ROOT,  ROOTQF,  STEPQF,  TANGQF,  UPQRQF,
C the LAPACK routines  DGEQRF,  DORGQR, their auxiliary routines,
C and the BLAS routines  DCOPY,  DDOT,  DGEMM,  DGEMV,
C   DGER,  DNRM2,  DSCAL,  DTPMV,  DTPSV,  DTRMM,  and   DTRMV.
C The module  REAL_PRECISION  specifies 64-bit real arithmetic,
C which the user may want to change.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y(1:N+1)  contains the starting point for tracking the homotopy map.
C    (Y(2),...,Y(N+1)) = A  for the fixed point and zero finding 
C    problems.  (Y(2),...,Y(N+1)) = X0  for the curve tracking problem.
C    Y(1)  need not be defined by the user.
C
C IFLAG can be -2, -1, 0, 2, or 3.  IFLAG should be 0 on the first
C    call to  FIXPQF  for the problem  X=F(X), -1 for the problem
C    F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.   In certain
C    situations  IFLAG  is set to 2 or 3 by  FIXPQF, and  FIXPQF  can
C    be called again without changing  IFLAG.
C
C ARCRE, ARCAE  are the relative and absolute errors, respectively,
C    allowed the quasi-Newton iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input, it is reset to  .5*SQRT(ANS?E).
C    Normally  ARC?E  should be considerably larger than  ANS?E.
C
C ANSRE, ANSAE  are the relative and absolute error values used for 
C    the answer at  LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
C    satisfies
C
C      |Y(1) - 1| .LE. ANSRE + ANSAE      .AND.
C  
C      ||DZ|| .LE. ANSRE*||Y|| + ANSAE      where
C
C      DZ is the quasi-Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems,  A  need not be initialized by the 
C    user, and is assumed to have length  N.  For the curve
C    tracking problem,  A  must be initialized by the user.
C
C SSPAR(1:4) =  (HMIN, HMAX, BMIN, BMAX)  is a vector of parameters 
C    used for the optimal step size estimation.  A default value
C    can be specified for any of these four parameters by setting it
C    .LE. 0.0  on input.  See the comments in  STEPQF  for more
C    information about these parameters.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and  Y  is an approximate
C    zero of the homotopy map.  Normally  LAMBDA = 1  and  X  is a
C    fixed point or zero of  F(X).   In abnormal situations,  LAMBDA
C    may only be near 1 and  X  near a fixed point or zero.
C
C IFLAG =
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE, ARCAE, ANSRE, ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPQF  again
C       without changing any parameters.
C
C   3   STEPQF  has been called 1000 times.  To continue, call
C       FIXPQF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the 
C       homotopy map and is not making progress.  The error 
C       tolerances  ARC?E  and  ANS?E  were too lenient.  The problem 
C       should be restrarted by calling  FIXPQF  with smaller error 
C       tolerances and  IFLAG = 0 (-1, -2).
C
C   6   The quasi-Newton iteration in  STEPQF  or  ROOTQF  failed to
C       converge.  The error tolerances  ANS?E  may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCRE, ARCAE, ANSRE, ANSAE  are unchanged after a normal return
C    (IFLAG = 1).  They are increased to appropriate values on the
C    return  IFLAG = 2.
C
C NFE  is the number of Jacobian evaluations.
C
C ARCLEN  is the approximate length of the zero curve.  
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to the
C    zero curve at the current point  Y.
C
C YOLD(1:N+1) is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1) is a work array containing the tangent vector to
C    the zero curve at  YOLD.
C
C Q(1:N+1,1:N+1), R((N+1)*(N+2)/2), F0(1:N+1), F1(1:N+1), Z0(1:N+1),
C    DZ(1:N+1), W(1:N+1), T(1:N+1), YSAV(1:N+1)  are all work arrays 
C    used by  STEPQF, TANGQF and ROOTQF to calculate the tangent 
C    vectors and quasi-Newton steps.
C
C
C ***** DECLARATIONS *****
      
      use Defs, only : scale, free0, limitd, epstol
C
C LIMITD IS AN UPPER BOUND ON THE NUMBER OF STEPS.
C

C
C     FUNCTION DECLARATIONS 
C 
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        INTEGER:: N,STRIDE
        real (kind=8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE

C
C     LOCAL VARIABLES 
C
      real (kind=8), SAVE:: ABSERR, H, HOLD, RELERR, S, WK 
      INTEGER, SAVE:: IFLAGC, ITER, JW, LIMIT, NP1
      LOGICAL, SAVE:: CRASH, START       
C
C     SCALAR ARGUMENTS 
C
      real (kind=8):: ARCRE, ARCAE, ANSRE, ANSAE, ARCLEN
      INTEGER:: N,IFLAG,TRACE,NFE
C
C     ARRAY DECLARATIONS 
C
      real (kind=8), DIMENSION(:), ALLOCATABLE, SAVE::
     &  R,YOLD,YP,YPOLD
      real (kind=8), DIMENSION(:,:), ALLOCATABLE, SAVE:: Q
      real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1),
     &    SSPAR(4), T(N+1), W(N+1), Y(:), YSAV(N+1), Z0(N+1)
C 
C ***** END OF DECLARATIONS *****
      INTERFACE
        SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
     &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
     &    F0,F1,Z0,DZ,W,T,SSPAR)
        INTEGER:: N, NFE, IFLAG
        LOGICAL:: START, CRASH
        real (kind=8):: HOLD, H, WK, RELERR, ABSERR, S
        real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &    Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
     &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
        END SUBROUTINE STEPQF
        SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &    YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
        real (kind=8):: RELERR, ABSERR
        INTEGER:: N, NFE, IFLAG
        real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &    Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
     &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
        END SUBROUTINE ROOTQF
      END INTERFACE
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C CHECK IFLAG
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
      IF (IFLAG .GE. -2 .AND. IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 50
      IF (IFLAG .EQ. 3) GO TO 40
C
C ONLY VALID INPUT FOR IFLAG IS -2, -1, 0, 2, 3.
C
      IFLAG = 7
      RETURN
C
C ***** INITIALIZATION BLOCK  *****
C
 10   ARCLEN = 0.0
      IF (ARCRE .LE. 0.0) ARCRE = .5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE = .5*SQRT(ANSAE)
      NFE=0
      IFLAGC = IFLAG
      NP1=N+1
      IF (ALLOCATED(Q)) DEALLOCATE(Q)
      IF (ALLOCATED(R)) DEALLOCATE(R)
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(Q(NP1,NP1),R(NP1*(N+2)/2),YOLD(NP1),YP(NP1),YPOLD(NP1),
     &  STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C 
C SET INITIAL CONDITIONS FOR FIRST CALL TO STEPQF.
C
      START=.TRUE.
      CRASH=.FALSE.
      RELERR = ARCRE
      ABSERR = ARCAE
      HOLD=1.0
      H=0.1
      S=0.0
      YPOLD(1) = 1.0
C      Y(1) = 0.0
      YPOLD(2:NP1)=0.0
C
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C
C     MINIMUM STEP SIZE HMIN
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)=(SQRT(N+1.0)+4.0)*epstol
C     MAXIMUM STEP SIZE HMAX
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= 1.0
C     MINIMUM STEP REDUCTION FACTOR BMIN
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= 0.1d0
C     MAXIMUM STEP EXPANSION FACTOR BMAX
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)= 7.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
C
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
C
40    LIMIT=LIMITD
C
C ***** END OF INITIALIZATION BLOCK. *****
C
      IF (TRACE .GT. 0) THEN
          WRITE(TRACE,218) Y(1), Y(2:NP1)/scale(free0(1:N))
 218      FORMAT(3X,e24.16,(1X,e24.16))
      ENDIF

50    DO ITER=1,LIMIT   ! ***** MAIN LOOP. *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN = S
        IFLAG = 5
        CALL CLEANUP ; RETURN
      END IF
C
C TAKE A STEP ALONG THE CURVE.
C
      if (trace .gt. 0 ) then
         print '(2x,i4,4x,"***",1x,f18.16)', iter, y(1)
      end if

      CALL STEPQF(N,NFE,IFLAGC,START,CRASH,HOLD,H,WK,
     &    RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,F0,F1,Z0,DZ,
     &    W,T,SSPAR)
C
C PRINT LATEST POINT ON CURVE IF REQUESTED.
C
      IF (TRACE .GT. 0) THEN
          WRITE(TRACE,217) Y(1), Y(2:NP1)/scale(free0(1:N))
 217      FORMAT(3X,e24.16,(1X,e24.16))
      ENDIF
C
C CHECK IF THE STEP WAS SUCCESSFUL.
C
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUP ; RETURN
      END IF
C
      IF (CRASH) THEN
C
C         RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
C      
        IFLAG=2
C
C         CHANGE ERROR TOLERANCES.
C
        IF (ARCRE .LT. RELERR) THEN
          ARCRE=RELERR
          ANSRE=RELERR
        END IF
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
C
C         CHANGE LIMIT ON NUMBER OF ITERATIONS.
C
        LIMIT = LIMIT - ITER
        RETURN
      END IF
C
C IF LAMBDA >= 1.0, USE ROOTQF TO FIND SOLUTION.
C
      IF (Y(1) .GE. 1.0) GOTO 500
C
      END DO   ! ***** END OF MAIN LOOP *****
C
C DID NOT CONVERGE IN  LIMIT  ITERATIONS, SET  IFLAG  AND RETURN.
C
      ARCLEN = S
      IFLAG = 3
      RETURN
C
C ***** FINAL STEP -- FIND SOLUTION AT LAMBDA=1 *****
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
C
C ***************************************************************
C modification 1: 27th july 2007 Modification by J. Gergaud
C We stop just after we have lambda>1 because we use HYBRD after
 500  RETURN
C 500  YSAV=YOLD
C end of the modification
C ***************************************************************
C
C FIND SOLUTION.
C
      CALL ROOTQF(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,
     &    YPOLD,A,Q,R,DZ,Z0,W,T,F0,F1)
C
C CHECK IF SOLUTION WAS FOUND AND SET  IFLAG  ACCORDINGLY.
C
      IFLAG=1
C
C     SET ERROR FLAG IF ROOTQF COULD NOT GET THE POINT ON THE ZERO
C     CURVE AT  LAMBDA = 1.0.
C
      IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C
C CALCULATE FINAL ARC LENGTH.
C
      DZ = Y - YSAV
      ARCLEN = S - HOLD + DNRM2(NP1,DZ,1)
C
C ***** END OF FINAL STEP *****
C
      CALL CLEANUP ; RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(Q)) DEALLOCATE(Q)
        IF (ALLOCATED(R)) DEALLOCATE(R)
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPQF

      END MODULE HOMPACK90


      SUBROUTINE ROOT(T,FT,B,C,RELERR,ABSERR,IFLAG)
C
C  ROOT COMPUTES A ROOT OF THE NONLINEAR EQUATION F(X)=0
C  WHERE F(X) IS A CONTINOUS REAL FUNCTION OF A SINGLE REAL
C  VARIABLE X.  THE METHOD USED IS A COMBINATION OF BISECTION
C  AND THE SECANT RULE.
C
C  NORMAL INPUT CONSISTS OF A CONTINUOS FUNCTION F AND AN
C  INTERVAL (B,C) SUCH THAT F(B)*F(C).LE.0.0.  EACH ITERATION
C  FINDS NEW VALUES OF B AND C SUCH THAT THE INTERVAL(B,C) IS
C  SHRUNK AND F(B)*F(C).LE.0.0.  THE STOPPING CRITERION IS
C
C          DABS(B-C).LE.2.0*(RELERR*DABS(B)+ABSERR)
C
C  WHERE RELERR=RELATIVE ERROR AND ABSERR=ABSOLUTE ERROR ARE
C  INPUT QUANTITIES.  SET THE FLAG, IFLAG, POSITIVE TO INITIALIZE
C  THE COMPUTATION.  AS B,C AND IFLAG ARE USED FOR BOTH INPUT AND
C  OUTPUT, THEY MUST BE VARIABLES IN THE CALLING PROGRAM.
C
C  IF 0 IS A POSSIBLE ROOT, ONE SHOULD NOT CHOOSE ABSERR=0.0.
C
C  THE OUTPUT VALUE OF B IS THE BETTER APPROXIMATION TO A ROOT
C  AS B AND C ARE ALWAYS REDEFINED SO THAT DABS(F(B)).LE.DABS(F(C)).
C
C  TO SOLVE THE EQUATION, ROOT MUST EVALUATE F(X) REPEATEDLY. THIS
C  IS DONE IN THE CALLING PROGRAM.  WHEN AN EVALUATION OF F IS
C  NEEDED AT T, ROOT RETURNS WITH IFLAG NEGATIVE.  EVALUATE FT=F(T)
C  AND CALL ROOT AGAIN.  DO NOT ALTER IFLAG.
C
C  WHEN THE COMPUTATION IS COMPLETE, ROOT RETURNS TO THE CALLING
C  PROGRAM WITH IFLAG POSITIVE=
C
C     IFLAG=1  IF F(B)*F(C).LT.0 AND THE STOPPING CRITERION IS MET.
C
C          =2  IF A VALUE B IS FOUND SUCH THAT THE COMPUTED VALUE
C              F(B) IS EXACTLY ZERO.  THE INTERVAL (B,C) MAY NOT
C              SATISFY THE STOPPING CRITERION.
C
C          =3  IF DABS(F(B)) EXCEEDS THE INPUT VALUES DABS(F(B)),
C              DABS(F(C)).  IN THIS CASE IT IS LIKELY THAT B IS CLOSE
C              TO A POLE OF F.
C
C          =4  IF NO ODD ORDER ROOT WAS FOUND IN THE INTERVAL.  A
C              LOCAL MINIMUM MAY HAVE BEEN OBTAINED.
C
C          =5  IF TOO MANY FUNCTION EVALUATIONS WERE MADE.
C              (AS PROGRAMMED, 500 ARE ALLOWED.)
C
C  THIS CODE IS A MODIFICATION OF THE CODE ZEROIN WHICH IS COMPLETELY
C  EXPLAINED AND DOCUMENTED IN THE TEXT  NUMERICAL COMPUTING:  AN
C  INTRODUCTION,  BY L. F. SHAMPINE AND R. C. ALLEN.
C
C
      use Defs, only : epstol

      real (kind=8):: A,ABSERR,ACBS,ACMB,AE,B,C,CMB,FA,FB,
     &  FC,FT,FX,P,Q,RE,RELERR,T,TOL,U
      INTEGER IC,IFLAG,KOUNT
      SAVE
C
      IF(IFLAG.GE.0) GO TO 100
      IFLAG=ABS(IFLAG)
      GO TO (200,300,400), IFLAG
  100 U=epstol
      RE=MAX(RELERR,U)
      AE=MAX(ABSERR,0.d0)
      IC=0
      ACBS=ABS(B-C)
      A=C
      T=A
      IFLAG=-1
      RETURN
  200 FA=FT
      T=B
      IFLAG=-2
      RETURN
  300 FB=FT
      FC=FA
      KOUNT=2
      FX=MAX(ABS(FB),ABS(FC))
    1 IF(ABS(FC).GE.ABS(FB))GO TO 2
C
C  INTERCHANGE B AND C SO THAT ABS(F(B)).LE.ABS(F(C)).
C
      A=B
      FA=FB
      B=C
      FB=FC
      C=A
      FC=FA
    2 CMB=0.5*(C-B)
      ACMB=ABS(CMB)
      TOL=RE*ABS(B)+AE
C
C  TEST STOPPING CRITERION AND FUNCTION COUNT.
C
      IF(ACMB.LE.TOL)GO TO 8
      IF(KOUNT.GE.500)GO TO 12
C
C  CALCULATE NEW ITERATE EXPLICITLY AS B+P/Q
C  WHERE WE ARRANGE P.GE.0.  THE IMPLICIT
C  FORM IS USED TO PREVENT OVERFLOW.
C
      P=(B-A)*FB
      Q=FA-FB
      IF(P.GE.0.0)GO TO 3
      P=-P
      Q=-Q
C
C  UPDATE A, CHECK IF REDUCTION IN THE SIZE OF BRACKETING
C  INTERVAL IS SATISFACTORY.  IF NOT BISECT UNTIL IT IS.
C
    3 A=B
      FA=FB
      IC=IC+1
      IF(IC.LT.4)GO TO 4
      IF(8.0*ACMB.GE.ACBS)GO TO 6
      IC=0
      ACBS=ACMB
C
C  TEST FOR TOO SMALL A CHANGE.
C
    4 IF(P.GT.ABS(Q)*TOL)GO TO 5
C
C  INCREMENT BY TOLERANCE
C
      B=B+SIGN(TOL,CMB)
      GO TO 7
C
C  ROOT OUGHT TO BE BETWEEN B AND (C+B)/2
C
    5 IF(P.GE.CMB*Q)GO TO 6
C
C  USE SECANT RULE.
C
      B=B+P/Q
      GO TO 7
C
C  USE BISECTION.
C
    6 B=0.5*(C+B)
C
C  HAVE COMPLETED COMPUTATION FOR NEW ITERATE B.
C
    7 T=B
      IFLAG=-3
      RETURN
  400 FB=FT
      IF(FB.EQ.0.0)GO TO 9
      KOUNT=KOUNT+1
      IF(SIGN(1.d0,FB).NE.SIGN(1.d0,FC))GO TO 1
      C=A
      FC=FA
      GO TO 1
C
C FINISHED.  SET IFLAG.
C
    8 IF(SIGN(1.d0,FB).EQ.SIGN(1.d0,FC))GO TO 11
      IF(ABS(FB).GT.FX)GO TO 10
      IFLAG=1
      RETURN
    9 IFLAG=2
      RETURN
   10 IFLAG=3
      RETURN
   11 IFLAG=4
      RETURN
   12 IFLAG=5
      RETURN
      END SUBROUTINE ROOT


      SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
     &   YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
C
C ROOTQF  FINDS THE POINT  YBAR = (1, XBAR)  ON THE ZERO CURVE OF THE
C HOMOTOPY MAP.  IT STARTS WITH TWO POINTS  YOLD=(LAMBDAOLD,XOLD)  AND
C Y=(LAMBDA,X)  SUCH THAT  LAMBDAOLD < 1 <= LAMBDA, AND ALTERNATES
C BETWEEN USING A SECANT METHOD TO FIND A PREDICTED POINT ON THE 
C HYPERPLANE  LAMBDA=1, AND TAKING A QUASI-NEWTON STEP TO RETURN TO THE 
C ZERO CURVE OF THE HOMOTOPY MAP.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE ROOTQF(N,NFE,IFLAG,RELERR,ABSERR,Y,YP,YOLD,
C    &    YPOLD,A,Q,R,DZ,Z,W,T,F0,F1)
C       
C       real (kind=8):: RELERR, ABSERR
C       INTEGER:: N, NFE, IFLAG
C       real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
C    &    Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
C    &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
C       END SUBROUTINE ROOTQF
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT Y=(LAMBDA,X) IS FOUND 
C    SUCH THAT
C
C    |Y(1) - 1| <= RELERR + ABSERR              AND
C
C    ||DZ|| <= RELERR*||Y|| + ABSERR,           WHERE
C
C    DZ  IS THE QUASI-NEWTON STEP TO Y.
C
C Y(1:N+1) = POINT (LAMBDA(S), X(S)) ON ZERO CURVE OF HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y.
C
C YOLD(1:N+1) = A POINT DIFFERENT FROM  Y  ON THE ZERO CURVE.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C Q(1:N+1,1:N+1) CONTAINS  Q  OF THE QR FACTORIZATION OF
C    THE AUGMENTED JACOBIAN MATRIX EVALUATED AT THE POINT Y.
C
C R((N+1)*(N+2)/2) CONTAINS THE UPPER TRIANGLE OF THE R PART OF
C    THE QR FACTORIZATION, STORED BY COLUMNS.
C
C DZ(1:N+1), Z(1:N+1), W(1:N+1), T(1:N+1), F0(1:N+1), F1(1:N+1)
C    ARE WORK ARRAYS USED FOR THE QUASI-NEWTON STEP AND THE SECANT
C    STEP.
C
C
C ON OUTPUT:
C
C N, RELERR, ABSERR, AND A  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG 
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A SINGULAR JACOBIAN MATRIX OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE.  Y  AND  YOLD  CONTAIN
C        THE LAST TWO POINTS OBTAINED BY QUASI-NEWTON STEPS, AND  YP
C        CONTAINS A POINT OPPOSITE OF THE HYPERPLANE  LAMBDA=1  FROM
C        Y.
C
C Y  IS THE POINT ON THE ZERO CURVE OF THE HOMOTOPY MAP AT  LAMBDA = 1.
C
C YP  AND  YOLD  CONTAIN POINTS NEAR THE SOLUTION.
C
C CALLS  DGEMV, DNRM2, DTPSV, F (OR RHO), ROOT, UPQRQF.
C
C ***** DECLARATIONS ***** 
      USE HOMOTOPY
      Use Defs, only : epstol
      
C
C FUNCTION DECLARATIONS 
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        INTEGER:: N,STRIDE
        real (kind=8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
      real (kind=8):: QOFS
C
C LOCAL VARIABLES 
C
      real (kind=8):: AERR, DD001, DD0011, DD01, DD011, DELS, ETA, 
     &   ONE, P0, P1, PP0, PP1, QSOUT, RERR, S, SA, SB, SOUT,
     &   U, ZERO
      INTEGER:: ISTEP, I, LCODE, LIMIT,NP1
      LOGICAL:: BRACK
C
C SCALAR ARGUMENTS 
C
      real (kind=8):: RELERR, ABSERR
      INTEGER:: N, NFE, IFLAG
C
C ARRAY DECLARATIONS 
C
      real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &   Q(N+1,N+1), R((N+1)*(N+2)/2), T(N+1), W(N+1),
     &   Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z(N+1)
C
C ***** END OF DECLARATIONS *****
C
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(P0,P1,DELS)=(P1-P0)/DELS
      DD001(P0,PP0,P1,DELS)=(DD01(P0,P1,DELS)-PP0)/DELS
      DD011(P0,P1,PP1,DELS)=(PP1-DD01(P0,P1,DELS))/DELS
      DD0011(P0,PP0,P1,PP1,DELS)=(DD011(P0,P1,PP1,DELS) - 
     &                              DD001(P0,PP0,P1,DELS))/DELS
      QOFS(P0,PP0,P1,PP1,DELS,S)=((DD0011(P0,PP0,P1,PP1,DELS)*
     &     (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
      
      
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C ***** INITIALIZATION *****
C
C ETA = PARAMETER FOR BROYDEN'S UPDATE.
C LIMIT = MAXIMUM NUMBER OF ITERATIONS ALLOWED.
C
      ONE=1.0
      ZERO=0.0
      U=epstol
      RERR=MAX(RELERR,U)
      AERR=MAX(ABSERR,ZERO)
      NP1=N+1
      ETA = 100.0*U
      LIMIT = 2*(INT(ABS(LOG10(AERR+RERR)))+1)
C
C F0 = (RHO(Y), YP*Y) TRANSPOSE.
C
      CALL RHO(A,Y(1),Y(2:NP1),F0(1:N))
      F0(NP1) = DOT_PRODUCT(YP,Y)
C
C ***** END OF INITIALIZATION BLOCK *****
C
C ***** COMPUTE FIRST INTERPOLANT WITH A HERMITE CUBIC *****
C
C FIND DISTANCE BETWEEN Y AND YOLD.  DZ=||Y-YOLD||.
C
      DZ = Y - YOLD
      DELS=DNRM2(NP1,DZ,1)
C
C USING TWO POINTS AND TANGENTS ON THE HOMOTOPY ZERO CURVE, CONSTRUCT 
C THE HERMITE CUBIC INTERPOLANT Q(S).  THEN USE  ROOT  TO FIND THE S
C CORRESPONDING TO  LAMBDA = 1.  THE TWO POINTS ON THE ZERO CURVE ARE
C ALWAYS CHOSEN TO BRACKET  LAMBDA=1, WITH THE BRACKETING INTERVAL
C ALWAYS BEING [0, DELS].
C
      SA=0.0
      SB=DELS
      LCODE=1
 40     CALL ROOT(SOUT,QSOUT,SA,SB,RERR,AERR,LCODE)
        IF (LCODE .GT. 0) GO TO 50
        QSOUT=QOFS(YOLD(1),YPOLD(1),Y(1),YP(1),DELS,SOUT) - 1.0
      GO TO 40
C
C IF  LAMBDA = 1  WERE BRACKETED,  ROOT  CANNOT FAIL.
C
 50   IF (LCODE .GT. 2) THEN
        IFLAG=6
        RETURN
      ENDIF
C
C CALCULATE  Q(SA)  AS THE INITIAL POINT FOR A NEWTON ITERATION.
C
      DO I=1,NP1
        Z(I)=QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),DELS,SA)
      END DO
C
C CALCULATE DZ = Z-Y.
C
      DZ = Z - Y
C
C ***** END OF CALCULATION OF CUBIC INTERPOLANT *****
C
C TANGENT INFORMATION  YPOLD  IS NO LONGER NEEDED.  HEREAFTER,  YPOLD
C REPRESENTS THE MOST RECENT POINT WHICH IS ON THE OPPOSITE SIDE OF
C LAMBDA=1  FROM  Y.
C
C ***** PREPARE FOR MAIN LOOP *****
C
      YPOLD = YOLD
C
C INITIALIZE  BRACK  TO INDICATE THAT THE POINTS  Y  AND  YOLD  BRACKET
C LAMBDA=1,  THUS YOLD = YPOLD.
C
      BRACK = .TRUE.
C
      DO ISTEP=1,LIMIT   ! ***** MAIN LOOP *****
C
C UPDATE JACOBIAN MATRIX.
C  
C       F1=(RHO(Z), YP*Z) TRANSPOSE.
C
         CALL RHO(A,Z(1),Z(2:NP1),F1(1:N))
         F1(NP1) = DOT_PRODUCT(YP,Z)
C
C
C PERFORM BROYDEN UPDATE.
C
         CALL UPQRQF(NP1,ETA,DZ,F0,F1,Q,R,W,T)
C
C QUASI-NEWTON STEP.
C
C COMPUTE NEWTON STEP.
C
         T(1:N) = -F1(1:N)
         T(NP1) = 0.0
         CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
         CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE NEWTON STEP.
C
         W = Z
         Z = Z + DZ
C
C CHECK FOR CONVERGENCE. 
C
         IF ((ABS(Z(1)-1.0) .LE. RERR+AERR) .AND. 
     &        (DNRM2(NP1,DZ,1) .LE. RERR*DNRM2(N,Z(2),1)+AERR)) THEN
            Y = Z
            RETURN
         END IF
C
C PREPARE FOR NEXT ITERATION.
C
         F0 = F1  
C
C IF  Z(1) = 1.0  THEN PERFORM QUASI-NEWTON ITERATION AGAIN
C WITHOUT COMPUTING A NEW PREDICTOR.
C
         IF (ABS(Z(1)-1.0) .LE. RERR+AERR) THEN
            DZ = Z - W
            CYCLE
         END IF
C
C     UPDATE  Y  AND  YOLD.
C     
         YOLD = Y
         Y = Z
C     
C UPDATE  YPOLD  SUCH THAT  YPOLD  IS THE MOST RECENT POINT 
C OPPOSITE OF  LAMBDA=1  FROM  Y.  SET  BRACK = .TRUE.  IFF  
C Y & YOLD  BRACKET  LAMBDA=1  SO THAT  YPOLD=YOLD. 
C
         IF ((Y(1)-1.0)*(YOLD(1)-1.0) .GT. 0) THEN
            BRACK = .FALSE.
         ELSE
            BRACK = .TRUE.
            YPOLD = YOLD
         END IF
C
C COMPUTE DELS = ||Y-YPOLD||.
C              
         DZ = Y - YPOLD
         DELS=DNRM2(NP1,DZ,1)
C
C COMPUTE  DZ  FOR THE LINEAR PREDICTOR   Z = Y + DZ,
C           WHERE  DZ = SA*(YOLD-Y).
C
         SA = (1.0-Y(1))/(YOLD(1)-Y(1))
         DZ = SA*(YOLD - Y)
C
C TO INSURE STABILITY, THE LINEAR PREDICTION MUST BE NO FARTHER 
C FROM  Y  THAN  YPOLD  IS.  THIS IS GUARANTEED IF  BRACK = .TRUE.
C IF LINEAR PREDICTION IS TOO FAR AWAY, USE BRACKETING POINTS
C COMPUTE LINEAR PREDICTION.
C
         IF (.NOT. BRACK) THEN
            IF (DNRM2(NP1,DZ,1) .GT. DELS) THEN
C
C COMPUTE  DZ = SA*(YPOLD-Y).
C          
               SA = (1.0-Y(1))/(YPOLD(1)-Y(1))
               DZ = SA*(YPOLD - Y)
            END IF
         END IF
C
C COMPUTE PREDICTOR Z = Y+DZ, AND DZ = NEW Z  - OLD Z (USED FOR
C QUASI-NEWTON UPDATE).
C
         Z = Z + DZ
         DZ = Z - W
      END DO                    ! ***** END OF MAIN LOOP. *****
C
C THE ALTERNATING OSCULATORY LINEAR PREDICTION AND QUASI-NEWTON 
C CORRECTION HAS NOT CONVERGED IN  LIMIT  STEPS.  ERROR RETURN.
      IFLAG=6
      RETURN
C
      END SUBROUTINE ROOTQF

      SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
     &   WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
     &   F0,F1,Z0,DZ,W,T,SSPAR)
C
C SUBROUTINE  STEPQF  TAKES ONE STEP ALONG THE ZERO CURVE OF THE 
C HOMOTOPY MAP  RHO(LAMBDA,X)  USING A PREDICTOR-CORRECTOR ALGORITHM.
C THE PREDICTOR USES LINEAR INTERPOLANT, AND THE CORRECTOR 
C RETURNS TO THE ZERO CURVE USING A QUASI-NEWTON ALGORITHM, REMAINING
C IN A HYPERPLANE PERPENDICULAR TO THE MOST RECENT TANGENT VECTOR.
C  STEPQF  ALSO ESTIMATES A STEP SIZE  H  FOR THE NEXT STEP ALONG THE 
C ZERO CURVE.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE STEPQF(N,NFE,IFLAG,START,CRASH,HOLD,H,
C    &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,
C    &    F0,F1,Z0,DZ,W,T,SSPAR)
C       
C       INTEGER:: N, NFE, IFLAG
C       LOGICAL:: START, CRASH
C       real (kind=8):: HOLD, H, WK, RELERR, ABSERR, S
C       real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
C    &    Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
C    &    Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
C       END SUBROUTINE STEPQF
C     END INTERFACE
C
C
C ON INPUT:
C 
C N = DIMENSION OF  X. 
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPQF, .FALSE. OTHERWISE.
C         SHOULD NOT BE MODIFIED BY THE USER AFTER THE FIRST CALL.
C
C HOLD = ||Y - YOLD|| ; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST
C    BE SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPQF.
C    THEREAFTER,  STEPQF  CALCULATES AN OPTIMAL VALUE FOR  H, AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE LAST STEP (COMPUTED BY PREVIOUS
C    CALL TO  STEPQF).  UNDEFINED ON FIRST CALL.  SHOULD NOT BE
C    MODIFIED BY THE USER.
C  
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION
C    IS CONSIDERED TO HAVE CONVERGED WHEN A POINT  Z=(LAMBDA,X)  IS 
C    FOUND SUCH THAT
C       ||DZ|| .LE. RELERR*||Z|| + ABSERR,
C    WHERE  DZ  IS THE LAST QUASI-NEWTON STEP.
C
C S  = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (LAMBDA(S), X(S)).
C
C Y(1:N+1) = PREVIOUS POINT (LAMBDA(S),X(S)) FOUND ON THE ZERO CURVE
C    OF THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  Y.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL
C    TO  STEPQF.
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY
C    MAP.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL TO 
C    STEPQF.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE 
C    HOMOTOPY MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C Q(1:N+1,1:N+1) =  Q  OF THE QR FACTORIZATION OF
C    THE AUGMENTED JACOBIAN MATRIX AT  Y.
C
C R((N+1)*(N+2)/2) = THE UPPER TRIANGLE  R  OF THE QR 
C    FACTORIZATION, STORED BY COLUMNS.
C
C F0(1:N+1), F1(1:N+1), Z0(1:N+1), DZ(1:N+1), W(1:N+1), T(1:N+1) ARE
C    WORK ARRAYS.  
C 
C SSPAR(1:4) = PARAMETERS USED FOR COMPUTATION OF THE OPTIMAL STEP SIZE.  
C    SSPAR(1) = HMIN, SSPAR(2) = HMAX, SSPAR(3) = BMIN, SSPAR(4) = BMAX.  
C    THE OPTIMAL STEP  H  IS RESTRICTED SUCH THAT 
C       HMIN .LE. H .LE. HMAX, AND  BMIN*HOLD .LE. H .LE. BMAX*HOLD.
C
C
C ON OUTPUT:
C
C NFE HAS BEEN UPDATED.
C
C IFLAG
C
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK <  N  HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE. 
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPQF  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPQF  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y-YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C     NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE STEP TAKEN BY  STEPQF.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY 
C     MAP UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.  THEY ARE POSSIBLY
C     CHANGED IF  CRASH  = .TRUE. (SEE DESCRIPTION OF  CRASH  ABOVE).
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C     VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C Q, R  STORE THE QR FACTORIZATION OF THE AUGMENTED JACOBIAN MATRIX 
C     EVALUATED AT  Y.
C
C
C CALLS  DGEMV, DGEQRF, DNRM2, DORGQR, DTPSV, F (OR RHO),
C     FJAC (OR RHOJAC), TANGQF, UPQRQF.
C
C ***** DECLARATIONS *****
      USE HOMOTOPY
      use Defs, only : epstol

C
C     FUNCTION DECLARATIONS  
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        INTEGER:: N,STRIDE
        real (kind=8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C     LOCAL VARIABLES
C
      real (kind=8):: ALPHA, DELS, ETA, FOURU, GAMMA, HFAIL, HTEMP,
     &  IDLERR, ONE, P0, P1, PP0, PP1, TEMP, TWOU, WKOLD, ZERO         
      INTEGER:: I, ITCNT, LITFH, J, JP1, NP1
      LOGICAL:: FAILED
C
C     SCALAR ARGUMENTS 
C
      INTEGER:: N, NFE, IFLAG
      LOGICAL:: START, CRASH
      real (kind=8):: HOLD, H, WK, RELERR, ABSERR, S
C
C     ARRAY DECLARATIONS
C
      real (kind=8):: A(:), DZ(N+1), F0(N+1), F1(N+1), 
     &   Q(N+1,N+1), R((N+1)*(N+2)/2), SSPAR(4), T(N+1), W(N+1),
     &   Y(:), YOLD(N+1), YP(N+1), YPOLD(N+1), Z0(N+1)
C
      SAVE
C
C ***** END OF DECLARATIONS *****
C
      INTERFACE
        SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
        USE HOMOTOPY
        INTEGER:: N, IFLAG, NFE
        real (kind=8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
     &    S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
        END SUBROUTINE TANGQF
      END INTERFACE
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C
C ***** INITIALIZATION *****
C
C ETA = PARAMETER FOR BROYDEN'S UPDATE.
C LITFH = MAXIMUM NUMBER OF QUASI-NEWTON ITERATIONS ALLOWED.
C
      ONE = 1.0
      ZERO = 0.0
      TWOU = 2.0*epstol
      FOURU = TWOU + TWOU
      NP1 = N+1
      FAILED = .FALSE.
      CRASH = .TRUE.
      ETA = 50.0*TWOU
      LITFH = 4*(INT(ABS(LOG10(ABSERR+RELERR)))+1)
C
C CHECK THAT ALL INPUT PARAMETERS ARE CORRECT.
C
C     THE ARCLENGTH  S MUST BE NONNEGATIVE.
C
      IF (S .LT. 0.0) RETURN
C
C     IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
C   
      IF (H .LT. FOURU*(1.0+S)) THEN
        H=FOURU*(1.0 + S)
        RETURN
      END IF
C
C     IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE 
C     VALUES.
C
      TEMP=DNRM2(NP1,Y,1) + 1.0
      IF (.5*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. 0.0) THEN
          RELERR = FOURU*(1.0+FOURU)
          TEMP = 0.0
          ABSERR = MAX(ABSERR,TEMP)
        ELSE
          ABSERR=FOURU*TEMP
        END IF
        RETURN
      END IF
C
C     INPUT PARAMETERS WERE ALL ACCEPTABLE.
C
      CRASH = .FALSE.
C
C COMPUTE  YP  ON FIRST CALL.
C NOTE:  DZ  IS USED SIMPLY AS A WORK ARRAY HERE.
C
      IF (START) THEN
         CALL TANGQF(Y,YP,YPOLD,A,Q,R,W,DZ,T,N,IFLAG,NFE)
         IF (IFLAG .GT. 0) RETURN
      END IF
C
C F0 = (RHO(Y), YP*Y) TRANSPOSE (DIFFERENT FOR EACH PROBLEM TYPE).
C
      CALL RHO(A,Y(1),Y(2:NP1),F0(1:N))
C
C DEFINE LAST ROW OF F0  =  YP*Y.
C
      F0(NP1) = DOT_PRODUCT(YP,Y)
C
C ***** END OF INITIALIZATION *****
C
C ***** COMPUTE PREDICTOR POINT Z0 *****
C
 20   if (y(1)+h*yp(1) > 1.d0) then ! prediction is going too far
C *****************************************************************
C modification 1: 27th july 2007
C author: J. Gergaud
         print *, 'prediction of lambda more than 1:  ', y(1)+h*yp(1)
C end of the modification
C ******************************************************************	
         h = (1.d0-y(1))/yp(1)
      end if
      Z0 = Y + H*YP
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
      CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
      F1(NP1) = DOT_PRODUCT(YP,Z0)
C
C ***** END OF PREDICTOR SECTION *****
C
C ***** SET-UP FOR QUASI-NEWTON ITERATION *****
C
      IF (FAILED) THEN
C        
C GENERATE Q = AUGMENTED JACOBIAN MATRIX FOR POINT Z0=(LAMBDA,X).
C
         call Jacfun(A,Z0(1),Z0(2:NP1),Q(1:N,1:NP1))
C
C DEFINE LAST ROW OF Q = YP.
C
         Q(NP1,:) = YP
C
C COUNT JACOBIAN EVALUATION.
C
         NFE = NFE+1
C
C DO FIRST QUASI-NEWTON STEP.
C
C FACTOR AUG.
C
         CALL DGEQRF(NP1,NP1,Q,NP1,T,W,NP1,I)
C
C PACK UPPER TRIANGLE INTO ARRAY R.
C
         DO I=1,NP1
            R((I*(I-1))/2 + 1:(I*(I-1))/2 + I) = Q(1:I,I)
         END DO
C
C CHECK FOR SINGULARITY.
C
         J = 1
         DO I = 1, N
            IF( R(J+I-1) .EQ. ZERO ) THEN
               IFLAG = 4
               RETURN
            END IF
            J = J + I
         END DO
C
C EXPAND HOUSEHOLDER REFLECTIONS INTO FULL MATRIX Q .
C
         CALL DORGQR(NP1, NP1, N, Q, NP1, T, W, NP1, I)
C
C COMPUTE NEWTON STEP.
C
         T(1:N) = -F1(1:N)
         T(NP1) = 0.0
         CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
         CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE STEP AND SET F0 = F1.
C
         Z0 = Z0 + DZ
         F0 = F1
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
         CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
         F1(NP1) = DOT_PRODUCT(YP,Z0)
C
      ELSE
C
C IF NOT FAILED THEN DEFINE  DZ=Z0-Y  PRIOR TO MAIN LOOP.
C
         DZ = Z0 - Y
      END IF
C
C ***** END OF PREPARATION FOR QUASI-NEWTON ITERATION *****
C
      DO ITCNT = 1,LITFH  ! ***** QUASI-NEWTON ITERATION *****
C
C PERFORM UPDATE FOR NEWTON STEP JUST TAKEN.
C
        CALL UPQRQF(NP1,ETA,DZ,F0,F1,Q,R,W,T)
C
C COMPUTE NEXT NEWTON STEP.
C
        T(1:N) = -F1(1:N)
        T(NP1) = 0.0
        CALL DGEMV('T',NP1,NP1,ONE,Q,NP1,T,1,ZERO,DZ,1)
        CALL DTPSV('U', 'N', 'N', NP1, R, DZ, 1)
C
C TAKE STEP.
C
        Z0 = Z0 + DZ
        if (z0(1) > 1.1d0) then ! Newton step is DV
           goto 170
        end if
C
C CHECK FOR CONVERGENCE.
C
        IF (DNRM2(NP1,DZ,1) .LE. RELERR*DNRM2(NP1,Z0,1)+ABSERR) THEN
           GO TO 180
        END IF

C
C IF NOT CONVERGED, PREPARE FOR NEXT ITERATION.
C
        F0 = F1
C
C F1 = (RHO(Z0), YP*Z0) TRANSPOSE.
C
        CALL RHO(A,Z0(1),Z0(2:NP1),F1(1:N))
        F1(NP1) = DOT_PRODUCT(YP,Z0)
C
      END DO                    ! ***** END OF QUASI-NEWTON LOOP *****
C
C ***** DIDN'T CONVERGE OR TANGENT AT NEW POINT DID NOT MAKE
C       AN ACUTE ANGLE WITH YPOLD -- TRY AGAIN WITH A SMALLER H *****
C      
 170  FAILED = .TRUE.
      HFAIL = H
      IF (H .LE. FOURU*(1.0 + S)) THEN
        IFLAG = 6
        RETURN
      ELSE
        H = .5 * H
      END IF
      GO TO 20
C
C ***** END OF CONVERGENCE FAILURE SECTION *****
C
C ***** CONVERGED -- MOP UP AND RETURN *****
C
C COMPUTE TANGENT & AUGMENTED JACOBIAN AT  Z0.
C NOTE:  DZ  AND  F1  ARE USED SIMPLY AS WORK ARRAYS HERE.
C
 180  CALL TANGQF(Z0,T,YP,A,Q,R,W,DZ,F1,N,IFLAG,NFE)
      IF (IFLAG .GT. 0) RETURN
C
C CHECK THAT COMPUTED TANGENT  T  MAKES AN ANGLE NO LARGER THAN
C 60 DEGREES WITH CURRENT TANGENT  YP.  (I.E. COS OF ANGLE < .5)
C IF NOT, STEP SIZE WAS TOO LARGE, SO THROW AWAY Z0, AND TRY
C AGAIN WITH A SMALLER STEP.
C
      ALPHA = DOT_PRODUCT(T,YP)
      IF (ALPHA .LT. 0.5) GOTO 170
      ALPHA = ACOS(ALPHA)
C
C SET UP VARIABLES FOR NEXT CALL.
C
      YOLD = Y
      Y = Z0
      YPOLD = YP
      YP = T
C
C UPDATE ARCLENGTH   S = S + ||Y-YOLD||.
C
      HTEMP = HOLD
      Z0 = Z0 - YOLD
      HOLD = DNRM2(NP1,Z0,1)
      S = S+HOLD
C
C COMPUTE OPTIMAL STEP SIZE. 
C   IDLERR = DESIRED ERROR FOR NEXT PREDICTOR STEP.
C   WK = APPROXIMATE CURVATURE = 2*SIN(ALPHA/2)/HOLD  WHERE 
C        ALPHA = ARCCOS(YP*YPOLD).
C   GAMMA = EXPECTED CURVATURE FOR NEXT STEP, COMPUTED BY 
C        EXTRAPOLATING FROM CURRENT CURVATURE  WK, AND LAST 
C        CURVATURE  WKOLD.  GAMMA IS FURTHER REQUIRED TO BE 
C        POSITIVE.
C
      IF (.NOT. START) WKOLD = WK
      IDLERR = SQRT(SQRT(ABSERR + RELERR*DNRM2(NP1,Y,1)))
C
C     IDLERR SHOULD BE NO BIGGER THAN 1/2 PREVIOUS STEP.
C
      IDLERR = MIN(.5*HOLD,IDLERR)
      WK = 2.0*ABS(SIN(.5*ALPHA))/HOLD
      IF (START) THEN
         GAMMA = WK
      ELSE 
         GAMMA = WK + HOLD/(HOLD+HTEMP)*(WK-WKOLD)
      END IF
      GAMMA = MAX(GAMMA, 0.01*ONE)
      H = SQRT(2.0*IDLERR/GAMMA)
C
C     ENFORCE RESTRICTIONS ON STEP SIZE SO AS TO ENSURE STABILITY.
C        HMIN <= H <= HMAX, BMIN*HOLD <= H <= BMAX*HOLD.
C
      H = MIN(MAX(SSPAR(1),SSPAR(3)*HOLD,H),SSPAR(4)*HOLD,SSPAR(2))
      IF (FAILED) H = MIN(HFAIL,H)
      START = .FALSE.
C
C ***** END OF MOP UP SECTION *****
C
      RETURN
C
      END SUBROUTINE STEPQF

      SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
C
C SUBROUTINE  TANGQF  COMPUTES THE UNIT TANGENT VECTOR  YP  TO THE
C ZERO CURVE OF THE HOMOTOPY MAP AT  Y  BY GENERATING THE AUGMENTED 
C JACOBIAN MATRIX  
C
C           --           --
C           |  D(RHO(Y))  |      
C     AUG = |        T    |,   WHERE RHO IS THE HOMOTOPY MAP,
C           |   YPOLD     | 
C           --           --
C
C SOLVING THE SYSTEM
C                                T
C         AUG*YPT = (0,0,...,0,1)    FOR YPT,
C
C AND FINALLY COMPUTING  YP = YPT/||YPT||.
C
C IN ADDITION, THE MATRIX AUG IS UPDATED SO THAT THE LAST ROW IS
C YP  INSTEAD OF  YPOLD  ON RETURN.
C
C THE FOLLOWING INTERFACE BLOCK SHOULD BE INCLUDED IN THE CALLING
C PROGRAM:
C
C     INTERFACE
C       SUBROUTINE TANGQF(Y,YP,YPOLD,A,Q,R,W,S,T,N,IFLAG,NFE)
C       USE HOMOTOPY
C       
C       INTEGER:: N, IFLAG, NFE
C       real (kind=8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
C    &    S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
C       END SUBROUTINE TANGQF
C     END INTERFACE
C
C
C ON INPUT:
C 
C Y(1:N+1) = CURRENT POINT (LAMBDA(S), X(S)).
C
C YP(1:N+1)  IS UNDEFINED ON INPUT.
C 
C YPOLD(1:N+1) = UNIT TANGENT VECTOR AT THE PREVIOUS POINT ON THE 
C    ZERO CURVE OF THE HOMOTOPY MAP.
C
C A(:)  IS THE PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C W(1:N+1), S(1:N+1), T(1:N+1)  ARE WORK ARRAYS.
C
C N  IS THE DIMENSION OF X, WHERE  Y=(LAMBDA(S),X(S)).
C
C IFLAG  IS -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C NFE  IS THE NUMBER OF JACOBIAN EVALUATIONS.
C
C
C ON OUTPUT:
C
C Y, YPOLD, A, N  ARE UNCHANGED.
C
C YP(1:N+1)  CONTAINS THE NEW UNIT TANGENT VECTOR TO THE ZERO
C    CURVE OF THE HOMOTOPY MAP AT  Y(S) = (LAMBDA(S), X(S)).
C
C Q(1:N+1,1:N+1)  CONTAINS  Q  OF THE QR FACTORIZATION OF
C    THE JACOBIAN MATRIX OF RHO EVALUATED AT  Y  AUGMENTED BY  
C    YP TRANSPOSE.
C
C R(1:(N+1)*(N+2)/2)  CONTAINS THE UPPER TRIANGLE (STORED BY COLUMNS)
C    OF THE  R  PART OF THE QR FACTORIZATION OF THE AUGMENTED JACOBIAN
C    MATRIX.
C
C IFLAG  = -2, -1, OR 0, (UNCHANGED) ON A NORMAL RETURN.
C        = 4 IF THE AUGMENTED JACOBIAN MATRIX HAS RANK LESS THAN N+1.
C 
C NFE  HAS BEEN INCREMENTED BY 1.
C
C
C CALLS  DGEQRF, DNRM2, DORGQR, DTPSV, F (OR RHO IF IFLAG = -2),
C FJAC (OR RHOJAC, IF IFLAG = -2),
C R1UPQF (WHICH IS AN ENTRY POINT OF UPQRQF).
C        
C ***** DECLARATIONS *****
      USE HOMOTOPY
      
C
C FUNCTION DECLARATIONS
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        INTEGER:: N,STRIDE
        real (kind=8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C LOCAL VARIABLES
C
      real (kind=8):: LAMBDA, YPNRM
      INTEGER:: I, J, JP1, NP1
C
C SCALAR ARGUMENTS
C
      INTEGER:: N, IFLAG, NFE
C
C ARRAY DECLARATIONS 
C
      real (kind=8):: A(:), Q(N+1,N+1), R((N+1)*(N+2)/2),
     &  S(N+1), T(N+1), W(N+1), Y(:), YP(N+1), YPOLD(N+1)
C
C ***** END OF DECLARATIONS *****
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
      NFE = NFE + 1
      NP1 = N + 1
      LAMBDA = Y(1)
C        
C ***** DEFINE THE AUGMENTED JACOBIAN MATRIX *****
C
C Q = AUG.
C
      call Jacfun(A,LAMBDA,Y(2:NP1),Q(1:N,1:NP1))
C
C DEFINE LAST ROW OF Q = YPOLD.
C
      Q(NP1,:) = YPOLD
C
C ***** END OF DEFINITION OF AUGMENTED JACOBIAN MATRIX *****
C
C                                          T
C ***** SOLVE SYSTEM  AUG*YPT = (0,...,0,1)  *****
C
C FACTOR MATRIX.
C           
      CALL DGEQRF(NP1,NP1,Q,NP1,T,W,NP1,I)
C
C PACK UPPER TRIANGLE INTO ARRAY R .
C
      DO I=1,NP1
        R((I*(I-1))/2 + 1:(I*(I-1))/2 + I) = Q(1:I,I)
      END DO
C
C IF MATRIX IS SINGULAR, THEN RETURN WITH IFLAG = 4,
C ELSE SOLVE SYSTEM  R*YP = QT*(0,...,0,1)  FOR YP.
C
C
C CHECK FOR SINGULARITY.
C
      J = 1
      DO I = 1, N
        IF( R( J+I-1 ).EQ. 0.0 ) THEN
          IFLAG = 4
          RETURN
        END IF
        J = J + I
      END DO
C
C EXPAND HOUSEHOLDER REFLECTIONS INTO FULL MATRIX Q . 
C
      CALL DORGQR(NP1, NP1, N, Q, NP1, T, W, NP1, I)
C
      YP = Q(NP1,:)
      CALL DTPSV('U', 'N', 'N', NP1, R, YP, 1)
C
C COMPUTE UNIT VECTOR.
C
      YPNRM = 1.0/DNRM2(NP1,YP,1)
      YP = YPNRM*YP
      if (YP(1).lt.0.d0) YP = - YP
C
C ***** SYSTEM SOLVED *****
C
C ***** UPDATE AUGMENTED SYSTEM SO THAT LAST ROW IS YP *****
C                        
C S=YP-YPOLD,  T = E(NP1)T*Q.
C      
      S = YP - YPOLD
      T = Q(NP1,:)
      CALL R1UPQF(NP1,S,T,Q,R,W)        
C
      RETURN
C
      END SUBROUTINE TANGQF

      SUBROUTINE UPQRQF(N,ETA,S,F0,F1,Q,R,W,T)
C
C SUBROUTINE  UPQRQF  PERFORMS A BROYDEN UPDATE ON THE  Q R  
C FACTORIZATION OF A MATRIX  A, (AN APPROXIMATION TO J(X0)), 
C RESULTING IN THE FACTORIZATION  Q+ R+ OF
C
C       A+  =  A  +  (Y - A*S) (ST)/(ST * S),
C
C (AN APPROXIMATION TO J(X1))
C WHERE S = X1 - X0, ST = S TRANSPOSE,  Y = F(X1) - F(X0).
C
C THE ENTRY POINT  R1UPQF  PERFORMS THE RANK ONE UPDATE ON THE QR
C FACTORIZATION OF 
C
C       A+ =  A + Q*(T*ST).
C
C
C ON INPUT:
C
C N  IS THE DIMENSION OF X AND F(X).
C
C ETA  IS A NOISE PARAMETER.  IF (Y-A*S)(I) .LE. ETA*(|F1(I)|+|F0(I)|)
C    FOR 1 .LE. I .LE. N, THEN NO UPDATE IS PERFORMED.
C
C S(1:N) = X1 - X0   (OR S FOR THE ENTRY POINT R1UPQF).
C
C F0(1:N) = F(X0).
C
C F1(1:N) = F(X1).
C
C Q(1:N,1:N)  CONTAINS THE OLD Q , WHERE  A = Q*R .
C
C R(1:N*(N+1)/2)  CONTAINS THE OLD R, STORED BY COLUMNS.
C
C W(1:N), T(1:N)  ARE WORK ARRAYS ( T  CONTAINS THE VECTOR T FOR THE
C    ENTRY POINT  R1UPQF ).
C
C 
C ON OUTPUT:
C
C N  AND  ETA  ARE UNCHANGED.
C
C Q  CONTAINS Q+ .
C
C R   CONTAINS R+, STORED BY COLUMNS.
C
C S, F0, F1, W, AND T  HAVE ALL BEEN CHANGED.
C
C
C CALLS   DGEMV, DNRM2, DTPMV.
C
C ***** DECLARATIONS *****
      
C
C FUNCTION DECLARATIONS 
C
      INTERFACE
        FUNCTION DNRM2(N,X,STRIDE)
        INTEGER:: N,STRIDE
        real (kind=8):: DNRM2,X(N)
        END FUNCTION DNRM2
      END INTERFACE
C
C LOCAL VARIABLES 
C
      real (kind=8):: C, DEN, ONE, SS, WW, YY, ZERO
      INTEGER:: I, INDEXC, INDEXD, INDXC2, J, K
      LOGICAL:: SKIPUP
C
C SCALAR ARGUMENTS 
C
      real (kind=8):: ETA
      INTEGER:: N
C
C ARRAY DECLARATIONS  
C
      real (kind=8)::  S(N), F0(N), F1(N), Q(N,N), R(N*(N+1)/2),
     &    W(N), T(N), TT(2)
C
C ***** END OF DECLARATIONS *****  
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
      ONE = 1.0
      ZERO = 0.0
      SKIPUP = .TRUE.
C
C ***** DEFINE T AND S SUCH THAT *****
C
C           A+ = Q*(R + T*ST). 
C
C T = R*S.
C
      T = S
      CALL DTPMV('U','N','N',N,R,T,1)
C
C W = Y - Q*T  = Y - A*S.
C
      W = F1 - F0 - MATMUL(Q,T)
C
C IF W(I) IS NOT SMALL, THEN UPDATE MUST BE PERFORMED,
C OTHERWISE SET W(I) TO 0.
C
      WHERE (ABS(W) .LE. ETA*(ABS(F1) + ABS(F0))) W = 0.0
      IF (ANY(ABS(W) .GT. ETA*(ABS(F1) + ABS(F0)))) SKIPUP = .FALSE.
C
C IF NO UPDATE IS NECESSARY, THEN RETURN.
C
      IF (SKIPUP) RETURN
C
C T = QT*W = QT*Y - R*S.
C
      CALL DGEMV('T',N,N,ONE,Q,N,W,1,ZERO,T,1)
C
C S = S/(ST*S).
C
      S = (1.0/DOT_PRODUCT(S,S))*S
C
C ***** END OF COMPUTATION OF  T & S      *****
C       AT THIS POINT,  A+ = Q*(R + T*ST). 
C
      ENTRY R1UPQF(N,S,T,Q,R,W)
C
C ***** COMPUTE THE QR FACTORIZATION Q- R- OF (R + T*S).  THEN,  *****
C       Q+ = Q*Q-,  AND  R+ = R-.
C
C FIND THE LARGEST  K  SUCH THAT  T(K) .NE. 0.
C
      K = N
      DO
        IF (T(K) .NE. 0.0 .OR. K .LE. 1) EXIT
        K=K-1
      END DO
C
C COMPUTE THE INDEX OF R(K-1,K-1).
C         
      INDEXD = (K*(K-1))/2
C
C ***** TRANSFORM R+T*ST INTO AN UPPER HESSENBERG MATRIX *****
C
C DETERMINE JACOBI ROTATIONS WHICH WILL ZERO OUT ROWS 
C N, N-1,...,2  OF THE MATRIX  T*ST,  AND APPLY THESE
C ROTATIONS TO  R.  (THIS IS EQUIVALENT TO APPLYING THE
C SAME ROTATIONS TO  R+T*ST, EXCEPT FOR THE FIRST ROW.
C THUS, AFTER AN ADJUSTMENT FOR THE FIRST ROW, THE 
C RESULT IS AN UPPER HESSENBERG MATRIX.  THE
C SUBDIAGONAL ELEMENTS OF WHICH WILL BE STORED IN  W.
C
C NOTE:  ROWS N,N-1,...,K+1 ARE ALREADY ALL ZERO.
C
      JACOBI: DO I=K-1,1,-1
C
C         DETERMINE THE JACOBI ROTATION WHICH WILL ZERO OUT
C         ROW  I+1  OF THE  T*ST  MATRIX.
C
        IF (T(I) .EQ. 0.0) THEN
          C = 0.0
C         SS = SIGN(-T(I+1))= -T(I+1)/|T(I+1)|
          SS = -SIGN(ONE,T(I+1))
        ELSE
          DEN = DNRM2(2,T(I),1)
          C = T(I) / DEN
          SS = -T(I+1)/DEN
        END IF
C
C         PREMULTIPLY  R  BY THE JACOBI ROTATION.
C
        YY = R(INDEXD)
        WW = 0.0
        R(INDEXD) = C*YY - SS*WW
        W(I+1) = SS*YY + C*WW
        DO J= I+1,N
C           YY = R(I,J)
C           WW = R(I+1,J)
            INDEXC = ((J-1)*J)/2 + I 
            INDXC2 = INDEXC + 1
            YY = R(INDEXC)
            WW = R(INDXC2)
C           R(I,J) = C*YY - SS*WW
C           R(I+1,J) = SS*YY + C*WW
            R(INDEXC) = C*YY - SS*WW
            R(INDXC2) = SS*YY + C*WW
        END DO
C
C         MULTIPLY  Q  BY THE JACOBI ROTATION.
C
        DO J=1,N
          YY = Q(J,I)
          WW = Q(J,I+1)
          Q(J,I) = C*YY - SS*WW
          Q(J,I+1) = SS*YY + C*WW
        END DO
C
C         UPDATE  T(I)  SO THAT  T(I)*ST(J)  IS THE  (I,J)TH  COMPONENT
C         OF  T*ST, PREMULTIPLIED BY ALL OF THE JACOBI ROTATIONS SO
C         FAR.
C
        IF (T(I) .EQ. 0.0) THEN
          T(I) = ABS(T(I+1))
        ELSE
          T(I) = DNRM2(2,T(I),1)
        END IF
C
C         LET INDEXD = THE INDEX OF R(I-1,I-1).
C
        INDEXD = INDEXD - I
C
      END DO JACOBI
C
C UPDATE THE FIRST ROW OF  R  SO THAT  R  HOLDS  (R+T*ST) 
C PREMULTIPLIED BY ALL OF THE ABOVE JACOBI ROTATIONS.
C
      J=1
      DO I=1,N 
        R(J) = T(1)*S(I) + R(J)
        J=I+J
      END DO
C
C ***** END OF TRANSFORMATION TO UPPER HESSENBERG *****
C
C
C ***** TRANSFORM UPPER HESSENBERG MATRIX INTO UPPER *****
C       TRIANGULAR MATRIX. 
C
C       INDEXD = INDEX OF R(I,I).
C        
      INDEXD = 1
      HESSEN: DO I=1,K-1
C
C         DETERMINE APPROPRIATE JACOBI ROTATION TO ZERO OUT
C         R(I+1,I).
C
        IF (R(INDEXD) .EQ. 0.0) THEN
          C = 0.0
          SS = -SIGN(ONE,W(I+1))
        ELSE
          TT(1) = R(INDEXD)
          TT(2) = W(I+1)
          DEN = DNRM2(2,TT,1)
          C = R(INDEXD) / DEN
          SS = -W(I+1)/DEN
        END IF
C
C         PREMULTIPLY  R  BY JACOBI ROTATION.
C
        YY = R(INDEXD)
        WW = W(I+1)
        R(INDEXD) = C*YY - SS*WW
        W(I+1) = 0.0
        DO J= I+1,N
C           YY = R(I,J)
C           WW = R(I+1,J)  
          INDEXC = ((J-1)*J)/2 + I
          INDXC2 = INDEXC + 1 
          YY = R(INDEXC)
          WW = R(INDXC2)
C           R(I,J) = C*YY -SS*WW
C           R(I+1,J) = SS*YY + C*WW
          R(INDEXC) = C*YY - SS*WW
          R(INDXC2) = SS*YY + C*WW
        END DO
        INDEXD = INDEXD + I + 1
C
C         MULTIPLY  Q  BY JACOBI ROTATION.
C
        DO J=1,N
          YY = Q(J,I)
          WW = Q(J,I+1)
          Q(J,I) = C*YY - SS*WW
          Q(J,I+1) = SS*YY + C*WW
        END DO
      END DO HESSEN
C
C ***** END OF TRANSFORMATION TO UPPER TRIANGULAR *****
C
C
C ***** END OF UPDATE *****
C
C
      RETURN
      END SUBROUTINE UPQRQF
