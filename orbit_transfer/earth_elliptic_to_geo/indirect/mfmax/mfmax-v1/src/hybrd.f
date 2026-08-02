      subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag,
     *                 mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,
     *                 qtf,wa1,wa2,wa3,wa4)
      integer n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr
      double precision xtol,epsfcn,factor
      double precision x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr),
     *                 qtf(n),wa1(n),wa2(n),wa3(n),wa4(n)
      external fcn
C     **********
c
C     subroutine hybrd
c
C     the purpose of hybrd is to find a zero of a system of
C     n nonlinear functions in n variables by a modification
C     of the powell hybrid method. the user must provide a
C     subroutine which calculates the functions. the jacobian is
C     then calculated by a forward-difference approximation.
c
C     the subroutine statement is
c
C       subroutine hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
C                        diag,mode,factor,nprint,info,nfev,fjac,
C                        ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
c
C     where
c
C       fcn is the name of the user-supplied subroutine which
C         calculates the functions. fcn must be declared
C         in an external statement in the user calling
C         program, and should be written as follows.
c
C         subroutine fcn(n,x,fvec,iflag)
C         integer n,iflag
C         double precision x(n),fvec(n)
C         ----------
C         calculate the functions at x and
C         return this vector in fvec.
C         ---------
C         return
C         end
c
C         the value of iflag should not be changed by fcn unless
C         the user wants to terminate execution of hybrd.
C         in this case set iflag to a negative integer.
c
C       n is a positive integer input variable set to the number
C         of functions and variables.
c
C       x is an array of length n. on input x must contain
C         an initial estimate of the solution vector. on output x
C         contains the final estimate of the solution vector.
c
C       fvec is an output array of length n which contains
C         the functions evaluated at the output x.
c
C       xtol is a nonnegative input variable. termination
C         occurs when the relative error between two consecutive
C         iterates is at most xtol.
c
C       maxfev is a positive integer input variable. termination
C         occurs when the number of calls to fcn is at least maxfev
C         by the end of an iteration.
c
C       ml is a nonnegative integer input variable which specifies
C         the number of subdiagonals within the band of the
C         jacobian matrix. if the jacobian is not banded, set
C         ml to at least n - 1.
c
C       mu is a nonnegative integer input variable which specifies
C         the number of superdiagonals within the band of the
C         jacobian matrix. if the jacobian is not banded, set
C         mu to at least n - 1.
c
C       epsfcn is an input variable used in determining a suitable
C         step length for the forward-difference approximation. this
C         approximation assumes that the relative errors in the
C         functions are of the order of epsfcn. if epsfcn is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
c
C       diag is an array of length n. if mode = 1 (see
C         below), diag is internally set. if mode = 2, diag
C         must contain positive entries that serve as
C         multiplicative scale factors for the variables.
c
C       mode is an integer input variable. if mode = 1, the
C         variables will be scaled internally. if mode = 2,
C         the scaling is specified by the input diag. other
C         values of mode are equivalent to mode = 1.
c
C       factor is a positive input variable used in determining the
C         initial step bound. this bound is set to the product of
C         factor and the euclidean norm of diag*x if nonzero, or else
C         to factor itself. in most cases factor should lie in the
C         interval (.1,100.). 100. is a generally recommended value.
c
C       nprint is an integer input variable that enables controlled
C         printing of iterates if it is positive. in this case,
C         fcn is called with iflag = 0 at the beginning of the first
C         iteration and every nprint iterations thereafter and
C         immediately prior to return, with x and fvec available
C         for printing. if nprint is not positive, no special calls
C         of fcn with iflag = 0 are made.
c
C       info is an integer output variable. if the user has
C         terminated execution, info is set to the (negative)
C         value of iflag. see description of fcn. otherwise,
C         info is set as follows.
c
C         info = 0   improper input parameters.
c
C         info = 1   relative error between two consecutive iterates
C                    is at most xtol.
c
C         info = 2   number of calls to fcn has reached or exceeded
C                    maxfev.
c
C         info = 3   xtol is too small. no further improvement in
C                    the approximate solution x is possible.
c
C         info = 4   iteration is not making good progress, as
C                    measured by the improvement from the last
C                    five jacobian evaluations.
c
C         info = 5   iteration is not making good progress, as
C                    measured by the improvement from the last
C                    ten iterations.
c
C       nfev is an integer output variable set to the number of
C         calls to fcn.
c
C       fjac is an output n by n array which contains the
C         orthogonal matrix q produced by the qr factorization
C         of the final approximate jacobian.
c
C       ldfjac is a positive integer input variable not less than n
C         which specifies the leading dimension of the array fjac.
c
C       r is an output array of length lr which contains the
C         upper triangular matrix produced by the qr factorization
C         of the final approximate jacobian, stored rowwise.
c
C       lr is a positive integer input variable not less than
C         (n*(n+1))/2.
c
C       qtf is an output array of length n which contains
C         the vector (q transpose)*fvec.
c
C       wa1, wa2, wa3, and wa4 are work arrays of length n.
c
C     subprograms called
c
C       user-supplied ...... fcn
c
C       minpack-supplied ... dogleg,dpmpar,enorm,fdjac1,
C                            qform,qrfac,r1mpyq,r1updt
c
C       fortran-supplied ... dabs,dmax1,dmin1,min0,mod
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2
      integer iwa(1)
      logical jeval,sing
      double precision actred,delta,epsmch,fnorm,fnorm1,one,pnorm,
     *                 prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,
     *                 zero
      double precision dpmpar,enorm
      data one,p1,p5,p001,p0001,zero
     *     /1.0d0,1.0d-1,5.0d-1,1.0d-3,1.0d-4,0.0d0/
c
C     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      info = 0
      iflag = 0
      nfev = 0
c
C     check the input parameters for errors.
c
      if (n .le. 0 .or. xtol .lt. zero .or. maxfev .le. 0
     *    .or. ml .lt. 0 .or. mu .lt. 0 .or. factor .le. zero
     *    .or. ldfjac .lt. n .or. lr .lt. (n*(n + 1))/2) go to 300
      if (mode .ne. 2) go to 20
      do 10 j = 1, n
         if (diag(j) .le. zero) go to 300
   10    continue
   20 continue
c
C     evaluate the function at the starting point
C     and calculate its norm.
c
      iflag = 1
      call fcn(n,x,fvec,iflag)
      nfev = 1
      if (iflag .lt. 0) go to 300
      fnorm = enorm(n,fvec)
c
C     determine the number of calls to fcn needed to compute
C     the jacobian matrix.
c
      msum = min0(ml+mu+1,n)
c
C     initialize iteration counter and monitors.
c
      iter = 1
      ncsuc = 0
      ncfail = 0
      nslow1 = 0
      nslow2 = 0
c
C     beginning of the outer loop.
c
   30 continue
         jeval = .true.
c
C        calculate the jacobian matrix.
c
         iflag = 2
         call fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,
     *               wa2)
         nfev = nfev + msum
         if (iflag .lt. 0) go to 300
c
C        compute the qr factorization of the jacobian.
c
         call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)
c
C        on the first iteration and if mode is 1, scale according
C        to the norms of the columns of the initial jacobian.
c
         if (iter .ne. 1) go to 70
         if (mode .eq. 2) go to 50
         do 40 j = 1, n
            diag(j) = wa2(j)
            if (wa2(j) .eq. zero) diag(j) = one
   40       continue
   50    continue
c
C        on the first iteration, calculate the norm of the scaled x
C        and initialize the step bound delta.
c
         do 60 j = 1, n
            wa3(j) = diag(j)*x(j)
   60       continue
         xnorm = enorm(n,wa3)
         delta = factor*xnorm
         if (delta .eq. zero) delta = factor
   70    continue
c
C        form (q transpose)*fvec and store in qtf.
c
         do 80 i = 1, n
            qtf(i) = fvec(i)
   80       continue
         do 120 j = 1, n
            if (fjac(j,j) .eq. zero) go to 110
            sum = zero
            do 90 i = j, n
               sum = sum + fjac(i,j)*qtf(i)
   90          continue
            temp = -sum/fjac(j,j)
            do 100 i = j, n
               qtf(i) = qtf(i) + fjac(i,j)*temp
  100          continue
  110       continue
  120       continue
c
C        copy the triangular factor of the qr factorization into r.
c
         sing = .false.
         do 150 j = 1, n
            l = j
            jm1 = j - 1
            if (jm1 .lt. 1) go to 140
            do 130 i = 1, jm1
               r(l) = fjac(i,j)
               l = l + n - i
  130          continue
  140       continue
            r(l) = wa1(j)
            if (wa1(j) .eq. zero) sing = .true.
  150       continue
c
C        accumulate the orthogonal factor in fjac.
c
         call qform(n,n,fjac,ldfjac,wa1)
c
C        rescale if necessary.
c
         if (mode .eq. 2) go to 170
         do 160 j = 1, n
            diag(j) = dmax1(diag(j),wa2(j))
  160       continue
  170    continue
c
C        beginning of the inner loop.
c
  180    continue
c
C           if requested, call fcn to enable printing of iterates.
c
            if (nprint .le. 0) go to 190
            iflag = 0
            if (mod(iter-1,nprint) .eq. 0) call fcn(n,x,fvec,iflag)
            if (iflag .lt. 0) go to 300
  190       continue
c
C           determine the direction p.
c
            call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)
c
C           store the direction p and x + p. calculate the norm of p.
c
            do 200 j = 1, n
               wa1(j) = -wa1(j)
               wa2(j) = x(j) + wa1(j)
               wa3(j) = diag(j)*wa1(j)
  200          continue
            pnorm = enorm(n,wa3)
c
C           on the first iteration, adjust the initial step bound.
c
            if (iter .eq. 1) delta = dmin1(delta,pnorm)
c
C           evaluate the function at x + p and calculate its norm.
c
            iflag = 1
            call fcn(n,wa2,wa4,iflag)
            nfev = nfev + 1
            if (iflag .lt. 0) go to 300
            fnorm1 = enorm(n,wa4)
c
C           compute the scaled actual reduction.
c
            actred = -one
            if (fnorm1 .lt. fnorm) actred = one - (fnorm1/fnorm)**2
c
C           compute the scaled predicted reduction.
c
            l = 1
            do 220 i = 1, n
               sum = zero
               do 210 j = i, n
                  sum = sum + r(l)*wa1(j)
                  l = l + 1
  210             continue
               wa3(i) = qtf(i) + sum
  220          continue
            temp = enorm(n,wa3)
            prered = zero
            if (temp .lt. fnorm) prered = one - (temp/fnorm)**2
c
C           compute the ratio of the actual to the predicted
C           reduction.
c
            ratio = zero
            if (prered .gt. zero) ratio = actred/prered
c
C           update the step bound.
c
            if (ratio .ge. p1) go to 230
               ncsuc = 0
               ncfail = ncfail + 1
               delta = p5*delta
               go to 240
  230       continue
               ncfail = 0
               ncsuc = ncsuc + 1
               if (ratio .ge. p5 .or. ncsuc .gt. 1)
     *            delta = dmax1(delta,pnorm/p5)
               if (dabs(ratio-one) .le. p1) delta = pnorm/p5
  240       continue
c
C           test for successful iteration.
c
            if (ratio .lt. p0001) go to 260
c
C           successful iteration. update x, fvec, and their norms.
c
            do 250 j = 1, n
               x(j) = wa2(j)
               wa2(j) = diag(j)*x(j)
               fvec(j) = wa4(j)
  250          continue
            xnorm = enorm(n,wa2)
            fnorm = fnorm1
            iter = iter + 1
  260       continue
c
C           determine the progress of the iteration.
c
            nslow1 = nslow1 + 1
            if (actred .ge. p001) nslow1 = 0
            if (jeval) nslow2 = nslow2 + 1
            if (actred .ge. p1) nslow2 = 0
c
C           test for convergence.
c
            if (delta .le. xtol*xnorm .or. fnorm .eq. zero) info = 1
            if (info .ne. 0) go to 300
c
C           tests for termination and stringent tolerances.
c
            if (nfev .ge. maxfev) info = 2
            if (p1*dmax1(p1*delta,pnorm) .le. epsmch*xnorm) info = 3
            if (nslow2 .eq. 5) info = 4
            if (nslow1 .eq. 10) info = 5
            if (info .ne. 0) go to 300
c
C           criterion for recalculating jacobian approximation
C           by forward differences.
c
            if (ncfail .eq. 2) go to 290
c
C           calculate the rank one modification to the jacobian
C           and update qtf if necessary.
c
            do 280 j = 1, n
               sum = zero
               do 270 i = 1, n
                  sum = sum + fjac(i,j)*wa4(i)
  270             continue
               wa2(j) = (sum - wa3(j))/pnorm
               wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
               if (ratio .ge. p0001) qtf(j) = sum
  280          continue
c
C           compute the qr factorization of the updated jacobian.
c
            call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
            call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)
            call r1mpyq(1,n,qtf,1,wa2,wa3)
c
C           end of the inner loop.
c
            jeval = .false.
            go to 180
  290    continue
c
C        end of the outer loop.
c
         go to 30
  300 continue
c
C     termination, either normal or user imposed.
c
      if (iflag .lt. 0) info = iflag
      iflag = 0
      if (nprint .gt. 0) call fcn(n,x,fvec,iflag)
      return
c
C     last card of subroutine hybrd.
c
      end
      subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
      integer n,lr
      double precision delta
      double precision r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
C     **********
c
C     subroutine dogleg
c
C     given an m by n matrix a, an n by n nonsingular diagonal
C     matrix d, an m-vector b, and a positive number delta, the
C     problem is to determine the convex combination x of the
C     gauss-newton and scaled gradient directions that minimizes
C     (a*x - b) in the least squares sense, subject to the
C     restriction that the euclidean norm of d*x be at most delta.
c
C     this subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     qr factorization of a. that is, if a = q*r, where q has
C     orthogonal columns and r is an upper triangular matrix,
C     then dogleg expects the full upper triangle of r and
C     the first n components of (q transpose)*b.
c
C     the subroutine statement is
c
C       subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
c
C     where
c
C       n is a positive integer input variable set to the order of r.
c
C       r is an input array of length lr which must contain the upper
C         triangular matrix r stored by rows.
c
C       lr is a positive integer input variable not less than
C         (n*(n+1))/2.
c
C       diag is an input array of length n which must contain the
C         diagonal elements of the matrix d.
c
C       qtb is an input array of length n which must contain the first
C         n elements of the vector (q transpose)*b.
c
C       delta is a positive input variable which specifies an upper
C         bound on the euclidean norm of d*x.
c
C       x is an output array of length n which contains the desired
C         convex combination of the gauss-newton direction and the
C         scaled gradient direction.
c
C       wa1 and wa2 are work arrays of length n.
c
C     subprograms called
c
C       minpack-supplied ... dpmpar,enorm
c
C       fortran-supplied ... dabs,dmax1,dmin1,dsqrt
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,j,jj,jp1,k,l
      double precision alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm,sum,
     *                 temp,zero
      double precision dpmpar,enorm
      data one,zero /1.0d0,0.0d0/
c
C     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
C     first, calculate the gauss-newton direction.
c
      jj = (n*(n + 1))/2 + 1
      do 50 k = 1, n
         j = n - k + 1
         jp1 = j + 1
         jj = jj - k
         l = jj + 1
         sum = zero
         if (n .lt. jp1) go to 20
         do 10 i = jp1, n
            sum = sum + r(l)*x(i)
            l = l + 1
   10       continue
   20    continue
         temp = r(jj)
         if (temp .ne. zero) go to 40
         l = j
         do 30 i = 1, j
            temp = dmax1(temp,dabs(r(l)))
            l = l + n - i
   30       continue
         temp = epsmch*temp
         if (temp .eq. zero) temp = epsmch
   40    continue
         x(j) = (qtb(j) - sum)/temp
   50    continue
c
C     test whether the gauss-newton direction is acceptable.
c
      do 60 j = 1, n
         wa1(j) = zero
         wa2(j) = diag(j)*x(j)
   60    continue
      qnorm = enorm(n,wa2)
      if (qnorm .le. delta) go to 140
c
C     the gauss-newton direction is not acceptable.
C     next, calculate the scaled gradient direction.
c
      l = 1
      do 80 j = 1, n
         temp = qtb(j)
         do 70 i = j, n
            wa1(i) = wa1(i) + r(l)*temp
            l = l + 1
   70       continue
         wa1(j) = wa1(j)/diag(j)
   80    continue
c
C     calculate the norm of the scaled gradient and test for
C     the special case in which the scaled gradient is zero.
c
      gnorm = enorm(n,wa1)
      sgnorm = zero
      alpha = delta/qnorm
      if (gnorm .eq. zero) go to 120
c
C     calculate the point along the scaled gradient
C     at which the quadratic is minimized.
c
      do 90 j = 1, n
         wa1(j) = (wa1(j)/gnorm)/diag(j)
   90    continue
      l = 1
      do 110 j = 1, n
         sum = zero
         do 100 i = j, n
            sum = sum + r(l)*wa1(i)
            l = l + 1
  100       continue
         wa2(j) = sum
  110    continue
      temp = enorm(n,wa2)
      sgnorm = (gnorm/temp)/temp
c
C     test whether the scaled gradient direction is acceptable.
c
      alpha = zero
      if (sgnorm .ge. delta) go to 120
c
C     the scaled gradient direction is not acceptable.
C     finally, calculate the point along the dogleg
C     at which the quadratic is minimized.
c
      bnorm = enorm(n,qtb)
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
      temp = temp - (delta/qnorm)*(sgnorm/delta)**2
     *       + dsqrt((temp-(delta/qnorm))**2
     *               +(one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
      alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp
  120 continue
c
C     form appropriate convex combination of the gauss-newton
C     direction and the scaled gradient direction.
c
      temp = (one - alpha)*dmin1(sgnorm,delta)
      do 130 j = 1, n
         x(j) = temp*wa1(j) + alpha*x(j)
  130    continue
  140 continue
      return
c
C     last card of subroutine dogleg.
c
      end
      double precision function dpmpar(i)
      integer i
C     **********
c
C     Function dpmpar
c
C     This function provides double precision machine parameters
C     when the appropriate set of data statements is activated (by
C     removing the c from column 1) and all other data statements are
C     rendered inactive. Most of the parameter values were obtained
C     from the corresponding Bell Laboratories Port Library function.
c
C     The function statement is
c
C       double precision function dpmpar(i)
c
C     where
c
C       i is an integer input variable set to 1, 2, or 3 which
C         selects the desired machine parameter. If the machine has
C         t base b digits and its smallest and largest exponents are
C         emin and emax, respectively, then these parameters are
c
C         dpmpar(1) = b**(1 - t), the machine precision,
c
C         dpmpar(2) = b**(emin - 1), the smallest magnitude,
c
C         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
c
C     Argonne National Laboratory. MINPACK Project. November 1996.
C     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More'
c
C     **********
      integer mcheps(4)
      integer minmag(4)
      integer maxmag(4)
      double precision dmach(3)
      equivalence (dmach(1),mcheps(1))
      equivalence (dmach(2),minmag(1))
      equivalence (dmach(3),maxmag(1))
c
C     Machine constants for the IBM 360/370 series,
C     the Amdahl 470/V6, the ICL 2900, the Itel AS/6,
C     the Xerox Sigma 5/7/9 and the Sel systems 85/86.
c
C     data mcheps(1),mcheps(2) / z34100000, z00000000 /
C     data minmag(1),minmag(2) / z00100000, z00000000 /
C     data maxmag(1),maxmag(2) / z7fffffff, zffffffff /
c
C     Machine constants for the Honeywell 600/6000 series.
c
C     data mcheps(1),mcheps(2) / o606400000000, o000000000000 /
C     data minmag(1),minmag(2) / o402400000000, o000000000000 /
C     data maxmag(1),maxmag(2) / o376777777777, o777777777777 /
c
C     Machine constants for the CDC 6000/7000 series.
c
C     data mcheps(1) / 15614000000000000000b /
C     data mcheps(2) / 15010000000000000000b /
c
C     data minmag(1) / 00604000000000000000b /
C     data minmag(2) / 00000000000000000000b /
c
C     data maxmag(1) / 37767777777777777777b /
C     data maxmag(2) / 37167777777777777777b /
c
C     Machine constants for the PDP-10 (KA processor).
c
C     data mcheps(1),mcheps(2) / "114400000000, "000000000000 /
C     data minmag(1),minmag(2) / "033400000000, "000000000000 /
C     data maxmag(1),maxmag(2) / "377777777777, "344777777777 /
c
C     Machine constants for the PDP-10 (KI processor).
c
C     data mcheps(1),mcheps(2) / "104400000000, "000000000000 /
C     data minmag(1),minmag(2) / "000400000000, "000000000000 /
C     data maxmag(1),maxmag(2) / "377777777777, "377777777777 /
c
C     Machine constants for the PDP-11. 
c
C     data mcheps(1),mcheps(2) /   9472,      0 /
C     data mcheps(3),mcheps(4) /      0,      0 /
c
C     data minmag(1),minmag(2) /    128,      0 /
C     data minmag(3),minmag(4) /      0,      0 /
c
C     data maxmag(1),maxmag(2) /  32767,     -1 /
C     data maxmag(3),maxmag(4) /     -1,     -1 /
c
C     Machine constants for the Burroughs 6700/7700 systems.
c
C     data mcheps(1) / o1451000000000000 /
C     data mcheps(2) / o0000000000000000 /
c
C     data minmag(1) / o1771000000000000 /
C     data minmag(2) / o7770000000000000 /
c
C     data maxmag(1) / o0777777777777777 /
C     data maxmag(2) / o7777777777777777 /
c
C     Machine constants for the Burroughs 5700 system.
c
C     data mcheps(1) / o1451000000000000 /
C     data mcheps(2) / o0000000000000000 /
c
C     data minmag(1) / o1771000000000000 /
C     data minmag(2) / o0000000000000000 /
c
C     data maxmag(1) / o0777777777777777 /
C     data maxmag(2) / o0007777777777777 /
c
C     Machine constants for the Burroughs 1700 system.
c
C     data mcheps(1) / zcc6800000 /
C     data mcheps(2) / z000000000 /
c
C     data minmag(1) / zc00800000 /
C     data minmag(2) / z000000000 /
c
C     data maxmag(1) / zdffffffff /
C     data maxmag(2) / zfffffffff /
c
C     Machine constants for the Univac 1100 series.
c
C     data mcheps(1),mcheps(2) / o170640000000, o000000000000 /
C     data minmag(1),minmag(2) / o000040000000, o000000000000 /
C     data maxmag(1),maxmag(2) / o377777777777, o777777777777 /
c
C     Machine constants for the Data General Eclipse S/200.
c
C     Note - it may be appropriate to include the following card -
C     static dmach(3)
c
C     data minmag/20k,3*0/,maxmag/77777k,3*177777k/
C     data mcheps/32020k,3*0/
c
C     Machine constants for the Harris 220.
c
C     data mcheps(1),mcheps(2) / '20000000, '00000334 /
C     data minmag(1),minmag(2) / '20000000, '00000201 /
C     data maxmag(1),maxmag(2) / '37777777, '37777577 /
c
C     Machine constants for the Cray-1.
c
C     data mcheps(1) / 0376424000000000000000b /
C     data mcheps(2) / 0000000000000000000000b /
c
C     data minmag(1) / 0200034000000000000000b /
C     data minmag(2) / 0000000000000000000000b /
c
C     data maxmag(1) / 0577777777777777777777b /
C     data maxmag(2) / 0000007777777777777776b /
c
C     Machine constants for the Prime 400.
c
C     data mcheps(1),mcheps(2) / :10000000000, :00000000123 /
C     data minmag(1),minmag(2) / :10000000000, :00000100000 /
C     data maxmag(1),maxmag(2) / :17777777777, :37777677776 /
c
C     Machine constants for the VAX-11.
c
C     data mcheps(1),mcheps(2) /   9472,  0 /
C     data minmag(1),minmag(2) /    128,  0 /
C     data maxmag(1),maxmag(2) / -32769, -1 /
c
C     Machine constants for IEEE machines.
c
      data dmach(1) /2.22044604926d-16/
      data dmach(2) /2.22507385852d-308/
      data dmach(3) /1.79769313485d+308/
c
      dpmpar = dmach(i)
      return
c
C     Last card of function dpmpar.
c
      end
      double precision function enorm(n,x)
      integer n
      double precision x(n)
C     **********
c
C     function enorm
c
C     given an n-vector x, this function calculates the
C     euclidean norm of x.
c
C     the euclidean norm is computed by accumulating the sum of
C     squares in three different sums. the sums of squares for the
C     small and large components are scaled so that no overflows
C     occur. non-destructive underflows are permitted. underflows
C     and overflows do not occur in the computation of the unscaled
C     sum of squares for the intermediate components.
C     the definitions of small, intermediate and large components
C     depend on two constants, rdwarf and rgiant. the main
C     restrictions on these constants are that rdwarf**2 not
C     underflow and rgiant**2 not overflow. the constants
C     given here are suitable for every known computer.
c
C     the function statement is
c
C       double precision function enorm(n,x)
c
C     where
c
C       n is a positive integer input variable.
c
C       x is an input array of length n.
c
C     subprograms called
c
C       fortran-supplied ... dabs,dsqrt
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,3.834d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
C              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
C              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
C           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
C     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     *         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
c
C     last card of function enorm.
c
      end
      subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
     *                  wa1,wa2)
      integer n,ldfjac,iflag,ml,mu
      double precision epsfcn
      double precision x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
C     **********
c
C     subroutine fdjac1
c
C     this subroutine computes a forward-difference approximation
C     to the n by n jacobian matrix associated with a specified
C     problem of n functions in n variables. if the jacobian has
C     a banded form, then function evaluations are saved by only
C     approximating the nonzero terms.
c
C     the subroutine statement is
c
C       subroutine fdjac1(fcn,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
C                         wa1,wa2)
c
C     where
c
C       fcn is the name of the user-supplied subroutine which
C         calculates the functions. fcn must be declared
C         in an external statement in the user calling
C         program, and should be written as follows.
c
C         subroutine fcn(n,x,fvec,iflag)
C         integer n,iflag
C         double precision x(n),fvec(n)
C         ----------
C         calculate the functions at x and
C         return this vector in fvec.
C         ----------
C         return
C         end
c
C         the value of iflag should not be changed by fcn unless
C         the user wants to terminate execution of fdjac1.
C         in this case set iflag to a negative integer.
c
C       n is a positive integer input variable set to the number
C         of functions and variables.
c
C       x is an input array of length n.
c
C       fvec is an input array of length n which must contain the
C         functions evaluated at x.
c
C       fjac is an output n by n array which contains the
C         approximation to the jacobian matrix evaluated at x.
c
C       ldfjac is a positive integer input variable not less than n
C         which specifies the leading dimension of the array fjac.
c
C       iflag is an integer variable which can be used to terminate
C         the execution of fdjac1. see description of fcn.
c
C       ml is a nonnegative integer input variable which specifies
C         the number of subdiagonals within the band of the
C         jacobian matrix. if the jacobian is not banded, set
C         ml to at least n - 1.
c
C       epsfcn is an input variable used in determining a suitable
C         step length for the forward-difference approximation. this
C         approximation assumes that the relative errors in the
C         functions are of the order of epsfcn. if epsfcn is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
c
C       mu is a nonnegative integer input variable which specifies
C         the number of superdiagonals within the band of the
C         jacobian matrix. if the jacobian is not banded, set
C         mu to at least n - 1.
c
C       wa1 and wa2 are work arrays of length n. if ml + mu + 1 is at
C         least n, then the jacobian is considered dense, and wa2 is
C         not referenced.
c
C     subprograms called
c
C       minpack-supplied ... dpmpar
c
C       fortran-supplied ... dabs,dmax1,dsqrt
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,j,k,msum
      double precision eps,epsmch,h,temp,zero
      double precision dpmpar
      data zero /0.0d0/
c
C     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = dsqrt(dmax1(epsfcn,epsmch))
      msum = ml + mu + 1
      if (msum .lt. n) go to 40
c
C        computation of dense approximate jacobian.
c
         do 20 j = 1, n
            temp = x(j)
            h = eps*dabs(temp)
            if (h .eq. zero) h = eps
            x(j) = temp + h
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 30
            x(j) = temp
            do 10 i = 1, n
               fjac(i,j) = (wa1(i) - fvec(i))/h
   10          continue
   20       continue
   30    continue
         go to 110
   40 continue
c
C        computation of banded approximate jacobian.
c
         do 90 k = 1, msum
            do 60 j = k, n, msum
               wa2(j) = x(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               x(j) = wa2(j) + h
   60          continue
            call fcn(n,x,wa1,iflag)
            if (iflag .lt. 0) go to 100
            do 80 j = k, n, msum
               x(j) = wa2(j)
               h = eps*dabs(wa2(j))
               if (h .eq. zero) h = eps
               do 70 i = 1, n
                  fjac(i,j) = zero
                  if (i .ge. j - mu .and. i .le. j + ml)
     *               fjac(i,j) = (wa1(i) - fvec(i))/h
   70             continue
   80          continue
   90       continue
  100    continue
  110 continue
      return
c
C     last card of subroutine fdjac1.
c
      end

      subroutine qform(m,n,q,ldq,wa)
      integer m,n,ldq
      double precision q(ldq,m),wa(m)
C     **********
c
C     subroutine qform
c
C     this subroutine proceeds from the computed qr factorization of
C     an m by n matrix a to accumulate the m by m orthogonal matrix
C     q from its factored form.
c
C     the subroutine statement is
c
C       subroutine qform(m,n,q,ldq,wa)
c
C     where
c
C       m is a positive integer input variable set to the number
C         of rows of a and the order of q.
c
C       n is a positive integer input variable set to the number
C         of columns of a.
c
C       q is an m by m array. on input the full lower trapezoid in
C         the first min(m,n) columns of q contains the factored form.
C         on output q has been accumulated into a square matrix.
c
C       ldq is a positive integer input variable not less than m
C         which specifies the leading dimension of the array q.
c
C       wa is a work array of length m.
c
C     subprograms called
c
C       fortran-supplied ... min0
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,j,jm1,k,l,minmn,np1
      double precision one,sum,temp,zero
      data one,zero /1.0d0,0.0d0/
c
C     zero out upper triangle of q in the first min(m,n) columns.
c
      minmn = min0(m,n)
      if (minmn .lt. 2) go to 30
      do 20 j = 2, minmn
         jm1 = j - 1
         do 10 i = 1, jm1
            q(i,j) = zero
   10       continue
   20    continue
   30 continue
c
C     initialize remaining columns to those of the identity matrix.
c
      np1 = n + 1
      if (m .lt. np1) go to 60
      do 50 j = np1, m
         do 40 i = 1, m
            q(i,j) = zero
   40       continue
         q(j,j) = one
   50    continue
   60 continue
c
C     accumulate q from its factored form.
c
      do 120 l = 1, minmn
         k = minmn - l + 1
         do 70 i = k, m
            wa(i) = q(i,k)
            q(i,k) = zero
   70       continue
         q(k,k) = one
         if (wa(k) .eq. zero) go to 110
         do 100 j = k, m
            sum = zero
            do 80 i = k, m
               sum = sum + q(i,j)*wa(i)
   80          continue
            temp = sum/wa(k)
            do 90 i = k, m
               q(i,j) = q(i,j) - temp*wa(i)
   90          continue
  100       continue
  110    continue
  120    continue
      return
c
C     last card of subroutine qform.
c
      end
      subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
      integer m,n,lda,lipvt
      integer ipvt(lipvt)
      logical pivot
      double precision a(lda,n),rdiag(n),acnorm(n),wa(n)
C     **********
c
C     subroutine qrfac
c
C     this subroutine uses householder transformations with column
C     pivoting (optional) to compute a qr factorization of the
C     m by n matrix a. that is, qrfac determines an orthogonal
C     matrix q, a permutation matrix p, and an upper trapezoidal
C     matrix r with diagonal elements of nonincreasing magnitude,
C     such that a*p = q*r. the householder transformation for
C     column k, k = 1,2,...,min(m,n), is of the form
c
C                           t
C           i - (1/u(k))*u*u
c
C     where u has zeros in the first k-1 positions. the form of
C     this transformation and the method of pivoting first
C     appeared in the corresponding linpack subroutine.
c
C     the subroutine statement is
c
C       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
c
C     where
c
C       m is a positive integer input variable set to the number
C         of rows of a.
c
C       n is a positive integer input variable set to the number
C         of columns of a.
c
C       a is an m by n array. on input a contains the matrix for
C         which the qr factorization is to be computed. on output
C         the strict upper trapezoidal part of a contains the strict
C         upper trapezoidal part of r, and the lower trapezoidal
C         part of a contains a factored form of q (the non-trivial
C         elements of the u vectors described above).
c
C       lda is a positive integer input variable not less than m
C         which specifies the leading dimension of the array a.
c
C       pivot is a logical input variable. if pivot is set true,
C         then column pivoting is enforced. if pivot is set false,
C         then no column pivoting is done.
c
C       ipvt is an integer output array of length lipvt. ipvt
C         defines the permutation matrix p such that a*p = q*r.
C         column j of p is column ipvt(j) of the identity matrix.
C         if pivot is false, ipvt is not referenced.
c
C       lipvt is a positive integer input variable. if pivot is false,
C         then lipvt may be as small as 1. if pivot is true, then
C         lipvt must be at least n.
c
C       rdiag is an output array of length n which contains the
C         diagonal elements of r.
c
C       acnorm is an output array of length n which contains the
C         norms of the corresponding columns of the input matrix a.
C         if this information is not needed, then acnorm can coincide
C         with rdiag.
c
C       wa is a work array of length n. if pivot is false, then wa
C         can coincide with rdiag.
c
C     subprograms called
c
C       minpack-supplied ... dpmpar,enorm
c
C       fortran-supplied ... dmax1,dsqrt,min0
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,j,jp1,k,kmax,minmn
      double precision ajnorm,epsmch,one,p05,sum,temp,zero
      double precision dpmpar,enorm
      data one,p05,zero /1.0d0,5.0d-2,0.0d0/
c
C     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
C     compute the initial column norms and initialize several arrays.
c
      do 10 j = 1, n
         acnorm(j) = enorm(m,a(1,j))
         rdiag(j) = acnorm(j)
         wa(j) = rdiag(j)
         if (pivot) ipvt(j) = j
   10    continue
c
C     reduce a to r with householder transformations.
c
      minmn = min0(m,n)
      do 110 j = 1, minmn
         if (.not.pivot) go to 40
c
C        bring the column of largest norm into the pivot position.
c
         kmax = j
         do 20 k = j, n
            if (rdiag(k) .gt. rdiag(kmax)) kmax = k
   20       continue
         if (kmax .eq. j) go to 40
         do 30 i = 1, m
            temp = a(i,j)
            a(i,j) = a(i,kmax)
            a(i,kmax) = temp
   30       continue
         rdiag(kmax) = rdiag(j)
         wa(kmax) = wa(j)
         k = ipvt(j)
         ipvt(j) = ipvt(kmax)
         ipvt(kmax) = k
   40    continue
c
C        compute the householder transformation to reduce the
C        j-th column of a to a multiple of the j-th unit vector.
c
         ajnorm = enorm(m-j+1,a(j,j))
         if (ajnorm .eq. zero) go to 100
         if (a(j,j) .lt. zero) ajnorm = -ajnorm
         do 50 i = j, m
            a(i,j) = a(i,j)/ajnorm
   50       continue
         a(j,j) = a(j,j) + one
c
C        apply the transformation to the remaining columns
C        and update the norms.
c
         jp1 = j + 1
         if (n .lt. jp1) go to 100
         do 90 k = jp1, n
            sum = zero
            do 60 i = j, m
               sum = sum + a(i,j)*a(i,k)
   60          continue
            temp = sum/a(j,j)
            do 70 i = j, m
               a(i,k) = a(i,k) - temp*a(i,j)
   70          continue
            if (.not.pivot .or. rdiag(k) .eq. zero) go to 80
            temp = a(j,k)/rdiag(k)
            rdiag(k) = rdiag(k)*dsqrt(dmax1(zero,one-temp**2))
            if (p05*(rdiag(k)/wa(k))**2 .gt. epsmch) go to 80
            rdiag(k) = enorm(m-j,a(jp1,k))
            wa(k) = rdiag(k)
   80       continue
   90       continue
  100    continue
         rdiag(j) = -ajnorm
  110    continue
      return
c
C     last card of subroutine qrfac.
c
      end
      subroutine r1mpyq(m,n,a,lda,v,w)
      integer m,n,lda
      double precision a(lda,n),v(n),w(n)
C     **********
c
C     subroutine r1mpyq
c
C     given an m by n matrix a, this subroutine computes a*q where
C     q is the product of 2*(n - 1) transformations
c
C           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
C     and gv(i), gw(i) are givens rotations in the (i,n) plane which
C     eliminate elements in the i-th and n-th planes, respectively.
C     q itself is not given, rather the information to recover the
C     gv, gw rotations is supplied.
c
C     the subroutine statement is
c
C       subroutine r1mpyq(m,n,a,lda,v,w)
c
C     where
c
C       m is a positive integer input variable set to the number
C         of rows of a.
c
C       n is a positive integer input variable set to the number
C         of columns of a.
c
C       a is an m by n array. on input a must contain the matrix
C         to be postmultiplied by the orthogonal matrix q
C         described above. on output a*q has replaced a.
c
C       lda is a positive integer input variable not less than m
C         which specifies the leading dimension of the array a.
c
C       v is an input array of length n. v(i) must contain the
C         information necessary to recover the givens rotation gv(i)
C         described above.
c
C       w is an input array of length n. w(i) must contain the
C         information necessary to recover the givens rotation gw(i)
C         described above.
c
C     subroutines called
c
C       fortran-supplied ... dabs,dsqrt
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
C     **********
      integer i,j,nmj,nm1
      double precision cos,one,sin,temp
      data one /1.0d0/
c
C     apply the first set of givens rotations to a.
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 50
      do 20 nmj = 1, nm1
         j = n - nmj
         if (dabs(v(j)) .gt. one) cos = one/v(j)
         if (dabs(v(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(v(j)) .le. one) sin = v(j)
         if (dabs(v(j)) .le. one) cos = dsqrt(one-sin**2)
         do 10 i = 1, m
            temp = cos*a(i,j) - sin*a(i,n)
            a(i,n) = sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   10       continue
   20    continue
c
C     apply the second set of givens rotations to a.
c
      do 40 j = 1, nm1
         if (dabs(w(j)) .gt. one) cos = one/w(j)
         if (dabs(w(j)) .gt. one) sin = dsqrt(one-cos**2)
         if (dabs(w(j)) .le. one) sin = w(j)
         if (dabs(w(j)) .le. one) cos = dsqrt(one-sin**2)
         do 30 i = 1, m
            temp = cos*a(i,j) + sin*a(i,n)
            a(i,n) = -sin*a(i,j) + cos*a(i,n)
            a(i,j) = temp
   30       continue
   40    continue
   50 continue
      return
c
C     last card of subroutine r1mpyq.
c
      end
      subroutine r1updt(m,n,s,ls,u,v,w,sing)
      integer m,n,ls
      logical sing
      double precision s(ls),u(m),v(n),w(m)
C     **********
c
C     subroutine r1updt
c
C     given an m by n lower trapezoidal matrix s, an m-vector u,
C     and an n-vector v, the problem is to determine an
C     orthogonal matrix q such that
c
C                   t
C           (s + u*v )*q
c
C     is again lower trapezoidal.
c
C     this subroutine determines q as the product of 2*(n - 1)
C     transformations
c
C           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
c
C     where gv(i), gw(i) are givens rotations in the (i,n) plane
C     which eliminate elements in the i-th and n-th planes,
C     respectively. q itself is not accumulated, rather the
C     information to recover the gv, gw rotations is returned.
c
C     the subroutine statement is
c
C       subroutine r1updt(m,n,s,ls,u,v,w,sing)
c
C     where
c
C       m is a positive integer input variable set to the number
C         of rows of s.
c
C       n is a positive integer input variable set to the number
C         of columns of s. n must not exceed m.
c
C       s is an array of length ls. on input s must contain the lower
C         trapezoidal matrix s stored by columns. on output s contains
C         the lower trapezoidal matrix produced as described above.
c
C       ls is a positive integer input variable not less than
C         (n*(2*m-n+1))/2.
c
C       u is an input array of length m which must contain the
C         vector u.
c
C       v is an array of length n. on input v must contain the vector
C         v. on output v(i) contains the information necessary to
C         recover the givens rotation gv(i) described above.
c
C       w is an output array of length m. w(i) contains information
C         necessary to recover the givens rotation gw(i) described
C         above.
c
C       sing is a logical output variable. sing is set true if any
C         of the diagonal elements of the output s are zero. otherwise
C         sing is set false.
c
C     subprograms called
c
C       minpack-supplied ... dpmpar
c
C       fortran-supplied ... dabs,dsqrt
c
C     argonne national laboratory. minpack project. march 1980.
C     burton s. garbow, kenneth e. hillstrom, jorge j. more,
C     john l. nazareth
c
C     **********
      integer i,j,jj,l,nmj,nm1
      double precision cos,cotan,giant,one,p5,p25,sin,tan,tau,temp,
     *                 zero
      double precision dpmpar
      data one,p5,p25,zero /1.0d0,5.0d-1,2.5d-1,0.0d0/
c
C     giant is the largest magnitude.
c
      giant = dpmpar(3)
c
C     initialize the diagonal element pointer.
c
      jj = (n*(2*m - n + 1))/2 - (m - n)
c
C     move the nontrivial part of the last column of s into w.
c
      l = jj
      do 10 i = n, m
         w(i) = s(l)
         l = l + 1
   10    continue
c
C     rotate the vector v into a multiple of the n-th unit vector
C     in such a way that a spike is introduced into w.
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 nmj = 1, nm1
         j = n - nmj
         jj = jj - (m - j + 1)
         w(j) = zero
         if (v(j) .eq. zero) go to 50
c
C        determine a givens rotation which eliminates the
C        j-th element of v.
c
         if (dabs(v(n)) .ge. dabs(v(j))) go to 20
            cotan = v(n)/v(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 30
   20    continue
            tan = v(j)/v(n)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
   30    continue
c
C        apply the transformation to v and store the information
C        necessary to recover the givens rotation.
c
         v(n) = sin*v(j) + cos*v(n)
         v(j) = tau
c
C        apply the transformation to s and extend the spike in w.
c
         l = jj
         do 40 i = j, m
            temp = cos*s(l) - sin*w(i)
            w(i) = sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
   40       continue
   50    continue
   60    continue
   70 continue
c
C     add the spike from the rank 1 update to w.
c
      do 80 i = 1, m
         w(i) = w(i) + v(n)*u(i)
   80    continue
c
C     eliminate the spike.
c
      sing = .false.
      if (nm1 .lt. 1) go to 140
      do 130 j = 1, nm1
         if (w(j) .eq. zero) go to 120
c
C        determine a givens rotation which eliminates the
C        j-th element of the spike.
c
         if (dabs(s(jj)) .ge. dabs(w(j))) go to 90
            cotan = s(jj)/w(j)
            sin = p5/dsqrt(p25+p25*cotan**2)
            cos = sin*cotan
            tau = one
            if (dabs(cos)*giant .gt. one) tau = one/cos
            go to 100
   90    continue
            tan = w(j)/s(jj)
            cos = p5/dsqrt(p25+p25*tan**2)
            sin = cos*tan
            tau = sin
  100    continue
c
C        apply the transformation to s and reduce the spike in w.
c
         l = jj
         do 110 i = j, m
            temp = cos*s(l) + sin*w(i)
            w(i) = -sin*s(l) + cos*w(i)
            s(l) = temp
            l = l + 1
  110       continue
c
C        store the information necessary to recover the
C        givens rotation.
c
         w(j) = tau
  120    continue
c
C        test for zero diagonal elements in the output s.
c
         if (s(jj) .eq. zero) sing = .true.
         jj = jj + (m - j + 1)
  130    continue
  140 continue
c
C     move w back into the last column of the output s.
c
      l = jj
      do 150 i = n, m
         s(l) = w(i)
         l = l + 1
  150    continue
      if (s(jj) .eq. zero) sing = .true.
      return
c
C     last card of subroutine r1updt.
c
      end
