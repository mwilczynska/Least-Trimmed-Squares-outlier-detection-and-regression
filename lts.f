c
c gfortran -o ltswm ltswm.f
c
      PROGRAM LTSWM
      integer npts,n,k,nmax,kmin
      parameter(nmax=500)
      character*30 filename
      character(len=20) fmt
      real*8 xval(nmax),err(nmax),totvar(nmax),errnew(nmax)
      real*8 wm,errwm,cs,csorig,s,ds,s2,sprev,wmorig,errwmorig
      real*8 uwm,siguwm,csprevious,errnewfull(nmax)
      real*8 x2(nmax),err2(nmax),uwm2,siguwm2,err3(nmax)
      real*8 wm2,errwm2,cs1,cs2
c      real*8 alo,ahi,result
      real*4 alo,ahi,result,csn,ecs
      real*4 XCHSQ,x
      real*8 p,xtest,PPND16
      real*8 sigwmnew,wmnew,wheightsumnew,wsumnew,sigwm,weightsum
      integer in,out,a(nmax),ifault,astore(nmax),i,t
      logical mtc
      external XCHSQ

c Initialise sigma_rand and its step size
      s=0.0
      ds=0.0001
c Read in the data
      write(6,*)' File name for data?'
      read(5,'(a30)')filename
      open(unit=1,file=filename,status='old')
      npts=0
      i=1
 1    read(1,*,err=1,end=2)xval(i),err(i)
      i=i+1
      npts=npts+1
      if(npts.gt.nmax)then
        write(6,*)' Max number of points (',nmax,') exceeded. Stop'
        stop
      endif
      goto 1
c Copy original error array into new array
 2    do i = 1,npts
        errnew(i) = err(i)
      end do

c ------------ WEIGHTED MEAN LOOP BEGINS HERE  ----------------
c Calculate weighted mean and chi-squared for the current iteration
      s=0.0
 3    call wmean(npts,xval,errnew,wm,errwm)
      call chisq(npts,xval,errnew,wm,cs)
c If normalised chi-squared > 1, increase sigma_rand, update error array:
      csn = cs/dble(npts)
      if(csn.gt.1.0)then
        s = s + ds
        do j = 1,npts
          errnew(j) = ( err(j)**2 + s**2 )**0.5
        end do
        goto 3
      endif
c ------------ WEIGHTED MEAN LOOP ENDS HERE  ----------------

c Weighted mean, chisq using original errors, i.e. s=0, and ordinary mean
      call wmean(npts,xval,err,wmorig,errwmorig)
      call chisq(npts,xval,err,wmorig,csorig)
      call mean(npts,xval,uwm,siguwm)
      write(6,*)' '
      write(6,*)' ----------------------------------------------------'
      write(6,*)' WEIGHTED AND ORDINARY MEANS:'
      write(6,*)' Total number of points used, n = ',npts
      write(6,*)' Sigma_rand for chisq_n=1 is ',s
      write(6,*)' Wmean (modified errors) = ',wm,' +/- ',errwm
      write(6,*)' Chi-squared (modified errors) = ',cs
      write(6,*)' '
      write(6,*)' Wmean (orig errors) = ',wmorig,' +/- ',errwmorig
      write(6,*)' Chi-squared (orig errors) = ',csorig
      write(6,*)' '
      write(6,*)' Ordinary mean = ',uwm,' +/- ',siguwm
      write(6,*)' ----------------------------------------------------'

c -------------------- START OF LTS CALCULATION -------------------------
c Iteratively repeat whole procedure using LTS for diminishing k,
c sigma_rand always zero. Stop when k = nint(0.85*n).
      csprevious=csorig
      n=npts
      k=nint(0.85*n)

c Compute the chisq_n expectation value for trimmed set
c First get the integration limits for chisq_n expectation value    
      p = (real(k)/real(npts) + 1.)/2.
c      p = 0.739
      ifault=0
      ahi=PPND16(p,IFAULT)
      alo=-abs(ahi)
      ahi=abs(ahi)
c      alo=-2.
c      ahi=2.
      write(6,*)' ----------------------------------------------------'
      write(6,*)' LTS trimming fraction:'
      write(6,*)' n, k = ',n,k
      write(6,*)' alo, ahi =',alo, ahi
      call QSIMP(XCHSQ,alo,ahi,result)
      write(6,*)' Expected Chi^2 =',result
      write(6,*)' ----------------------------------------------------'
      write(6,*)' Calculating s for all combinations of k points'
      t=0
      sprev=1000
11    call nxksrd(n,k,a,mtc,in,out)
      t=t+1
c DEBUG: Write out all possible combinations.  Variable format length:
c http://gcc.gnu.org/onlinedocs/gfortran/Variable-FORMAT-expressions.html
c      write(fmt,*)k
c      write(6,"(1x,"//adjustl(fmt)//"(i2,1x))")a(1:k)
c      write(6,*)' Sigma_rand = ',s
c      write(6,*)' s,csn,wm,errwm,t,',s,csn,wm,errwm,t
c      write(6,*)' Just after calling NXKSRD'

      do j=1,k
        x2(j)=xval(a(j))
        err2(j)=err(a(j))
      end do

ccccccccccccccccccccccccccccc
      s=0.0
 34   call wmean(k,x2,err2,wm,errwm)
      call chisq(k,x2,err2,wm,cs)
c If normalised chi-squared > "result", increase sigma_rand, update error array:
      csn = cs/(real(k)-1.0)
      if(csn.gt.result)then
        s = s + ds
        do j = 1,k
          err2(j) = ( err(a(j))**2 + s**2 )**0.5
        end do
c      write(6,*)' csn',csn
c      write(6,*)' result,', result
c      write(6,*)' Sigma_rand = ',s
c      write(6,*)' csn = ',csn
c      write(6,*)' Wmean = ',wm,' +/- ',errwm

c      do j=1,k
c      write(6,*),x2(j),err2(j)
c      end do

c      write(6,*)' s, csn, wm, errwm,', s, csn, wm, errwm

c DEBUG
c        if(abs(csn-1.).le.0.001)then
c        write(6,*)' s, csn, wm, errwm,', s, csn, wm, errwm
c        do j=1,k
c        write(6,*)' x2, err2 = ',x2(j),err2(j)
c        end do
c        endif
c END DEBUG

        goto 34
      endif

ccccccccccccccccccccccccccccc
c For current combination of k from n, identify s_min and calculate weighted mean
      call wmean(j,x2,err2,wm2,errwm2)
      call chisq(j,x2,err2,wm2,cs2)
c Store the set of points giving smallest chisq. (csorig is chisq for all n points)
      if(s.lt.sprev)s2=s

	  if(s2.eq.s)then
        do i=1,k
          astore(i)=a(i)
        end do
      write(6,*)' ----------------------------------------------------'
      write(fmt,*)k
      write(6,"(1x,"//adjustl(fmt)//"(i2,1x))")astore(1:k)
      endif
      write(6,*)' s, s2, sprev, t = ',s, s2, sprev, t
      write(6,*)' csn, wm, errwm,', csn, wm, errwm
      write(6,*)' ----------------------------------------------------'

      sprev=s2
      if(mtc)goto 11

      write(6,*)' k points with lowest sigma_rand:'
      write(fmt,*)k
      write(6,"(1x,"//adjustl(fmt)//"(i2,1x))")astore(1:k)
      do j=1,k
        x2(j)=xval(astore(j))
        err2(j)=err(astore(j))
        write(6,*)' x2(j), err2(j): ',x2(j),err2(j)
      end do

c      call wmean(j,x2,err2,wm2,errwm2)
c      call chisq(j,x2,err2,wm2,cs2)
c      write(6,*)' s, csn, wm, errwm, t,', s2, csn, wm, errwm, t


      write(6,*)' ----------------------------------------------------'
      write(6,*)' Sigma_rand applied to enire (n) dataset:'
      do i = 1,npts
        errnewfull(i) = ( err(i)**2 + s2**2 )**0.5
      write(6,*),xval(i), errnewfull(i)
      end do
      write(6,*)' ----------------------------------------------------'

      call wmean(npts,xval,errnewfull,wm,errwm)
      call chisq(npts,xval,errnewfull,wm,cs)
      write(6,*)' Values for whole untrimmed sample, '
      write(6,*)' using sigma_rand from trimmed sample'
      write(6,*)' Sigma_rand = ',s2
      write(6,*)' Chisq = ',cs
      write(6,*)' Wmean = ',wm,' +/- ',errwm

      write(6,*)' 3-sigma deviations relative to weighted mean:'
      do i = 1,npts
        x = abs((xval(i) - wm)/errnewfull(i))
c        if(x.ge.3.0)
c     $ write(6,*)' i, xval, err, nsig =',i,xval(i),errnewfull(i),x
      write(6,*)' i, xval, err, nsig =',i,xval(i),errnewfull(i),x
      end do

 40   stop
      end


c-------------------------------------------------------------
c Subroutine CHISQ
c Calculate chi-squared using a given mean
c Inputs:  n - # of points, x - array, sigx - errors on x, m - a mean
c Returns: cs - chi-squared (not normalised).
c-------------------------------------------------------------
      subroutine chisq(n,x,sigx,m,cs)
      real*8 cs,x(n),sigx(n),m
      integer i,n
      cs=0.0
      do i=1,n
        cs = cs + ( (x(i) - m)/sigx(i) )**2
      end do
      return
      end

c-------------------------------------------------------------
c Subroutine WMEAN
c Calculate a simple weighted mean and its error
c Inputs:  n - # of points, x - array, sigx - errors on x
c Returns: wm - weighted mean, sigwm - error on weighted mean.
c-------------------------------------------------------------
      subroutine wmean(n,x,sigx,wm,sigwm)
      real*8 wm,wsum,weightsum,x(n),sigx(n),sigwm
      integer i,n
      wsum=0.0
      weightsum=0.0
      do i=1,n
        wsum=wsum + ( x(i)/(sigx(i)**2) )
        weightsum=weightsum + (1.0/sigx(i)**2)
      end do
      wm=wsum/weightsum
      sigwm = (1.0/weightsum)**0.5
      return
      end

c--------------------------------------------------
c Subroutine MEAN
c Calculate a simple unweighted mean and its error.
c--------------------------------------------------
      subroutine mean(n,x,uwm,siguwm)
      real*8 uwm,siguwm,x(n),sum,sd
      integer n
      sum=0.0
      do i=1,n
        sum=sum+x(i)
      end do
      uwm=sum/dble(n)
      sum=0.0
      do i=1,n
        sum=sum + (x(i) - uwm)**2
      end do
      sd=(sum/dble(n))**0.5
      siguwm=sd/dble(n)**0.5

      return
      end

c-------------------------------------------------------------------------------
c http://www.cs.sunysb.edu/~algorith/implement/wilf/distrib/processed/nxksrd_2.f
c http://www.cs.sunysb.edu/~algorith/
c-------------------------------------------------------------------------------
c Chapter 3: Next k-subset of an n-set(p33)
c-------------------------------------------------------------------------------
c Name of subroutine: NXKSRD
c
c Algorithm:Generating subsets of {1,2,...,n} which succeed input set, 
c in lexicographical order, with optional jumps over supersets.     
c
c Inputs: n and k
c-------------------------------------------------------------------------------
      subroutine nxksrd(n,k,a,mtc,in,out)
      integer a(k),out
      logical mtc
      if(mtc) goto 10
      do 1 i=1,k
    1 a(i)=i
      mtc=k.ne.n
      return
   10 j=0
   20 if(mod(k,2).ne.0)goto 100
   30 j=j+1
c 40 if(j.le.k)goto 60
c 50 a(k)=k
c in=k
c out=n
c return
   60 if(a(j).eq.j)goto 100
   70 out=a(j)
      in=out-1
      a(j)=in
   80 if(j.eq.1)goto 200
   90 in=j-1
      a(j-1)=in
      goto 200
  100 j=j+1
  130 m=n
  110 if(j.lt.k)m=a(j+1)-1
  140 if(m.eq.a(j)) goto 30
  150 in=a(j)+1
      a(j)=in
      out=in-1
      if(j.eq.1)goto 200
      a(j-1)=out
      out=j-1
  200 if(k.eq.1)goto 201
      mtc=a(k-1).eq.k-1
  201 mtc=(.not.mtc).or.a(k).ne.n
      return
      end

c ------------------------------------------------------------------------
c Function to compute the inverse normal cumulative distribution function.
c (see e.g. pages 118-119 Julian King's PhD thesis).
c A description and links to various implementations is here:
c http://home.online.no/~pjacklam/notes/invnorm/
c The code below (claimed accuracy is 10^{-16}) was obtained from:
c http://lib.stat.cmu.edu/apstat/241
c JKW 15/9/12
C ------------------------------------------------------------------------
C	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
C
C	Produces the normal deviate Z corresponding to a given lower
C	tail area of P; Z is accurate to about 1 part in 10**16.
C
C	The hash sums below are the sums of the mantissas of the
C	coefficients.   They are included for use in checking transcription.
C ------------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION PPND16 (P, IFAULT)
	DOUBLE PRECISION ZERO, ONE, HALF, SPLIT1, SPLIT2, CONST1,
     *		CONST2, A0, A1,	A2, A3, A4, A5, A6, A7, B1, B2, B3,
     *          B4, B5, B6, B7,
     *		C0, C1, C2, C3, C4, C5, C6, C7,	D1, D2, D3, D4, D5,
     *		D6, D7, E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3,
     *		F4, F5, F6, F7, P, Q, R
	PARAMETER (ZERO = 0.D0, ONE = 1.D0, HALF = 0.5D0,
     *		SPLIT1 = 0.425D0, SPLIT2 = 5.D0,
     *		CONST1 = 0.180625D0, CONST2 = 1.6D0)
C
C	Coefficients for P close to 0.5
C
	PARAMETER (A0 = 3.38713 28727 96366 6080D0,
     *		   A1 = 1.33141 66789 17843 7745D+2,
     *		   A2 = 1.97159 09503 06551 4427D+3,
     *		   A3 = 1.37316 93765 50946 1125D+4,
     *		   A4 = 4.59219 53931 54987 1457D+4,
     *		   A5 = 6.72657 70927 00870 0853D+4,
     *		   A6 = 3.34305 75583 58812 8105D+4,
     *		   A7 = 2.50908 09287 30122 6727D+3,
     *		   B1 = 4.23133 30701 60091 1252D+1,
     *		   B2 = 6.87187 00749 20579 0830D+2,
     *		   B3 = 5.39419 60214 24751 1077D+3,
     *		   B4 = 2.12137 94301 58659 5867D+4,
     *		   B5 = 3.93078 95800 09271 0610D+4,
     *		   B6 = 2.87290 85735 72194 2674D+4,
     *		   B7 = 5.22649 52788 52854 5610D+3)
C	HASH SUM AB    55.88319 28806 14901 4439
C
C	Coefficients for P not close to 0, 0.5 or 1.
C
	PARAMETER (C0 = 1.42343 71107 49683 57734D0,
     *		   C1 = 4.63033 78461 56545 29590D0,
     *		   C2 = 5.76949 72214 60691 40550D0,
     *		   C3 = 3.64784 83247 63204 60504D0,
     *		   C4 = 1.27045 82524 52368 38258D0,
     *		   C5 = 2.41780 72517 74506 11770D-1,
     *         C6 = 2.27238 44989 26918 45833D-2,
     *		   C7 = 7.74545 01427 83414 07640D-4,
     *		   D1 = 2.05319 16266 37758 82187D0,
     *		   D2 = 1.67638 48301 83803 84940D0,
     *		   D3 = 6.89767 33498 51000 04550D-1,
     *		   D4 = 1.48103 97642 74800 74590D-1,
     *		   D5 = 1.51986 66563 61645 71966D-2,
     *		   D6 = 5.47593 80849 95344 94600D-4,
     *		   D7 = 1.05075 00716 44416 84324D-9)
C	HASH SUM CD    49.33206 50330 16102 89036
C
C	Coefficients for P near 0 or 1.
C
	PARAMETER (E0 = 6.65790 46435 01103 77720D0,
     *		   E1 = 5.46378 49111 64114 36990D0,
     *		   E2 = 1.78482 65399 17291 33580D0,
     *		   E3 = 2.96560 57182 85048 91230D-1,
     *		   E4 = 2.65321 89526 57612 30930D-2,
     *		   E5 = 1.24266 09473 88078 43860D-3,
     *		   E6 = 2.71155 55687 43487 57815D-5,
     *		   E7 = 2.01033 43992 92288 13265D-7,
     *		   F1 = 5.99832 20655 58879 37690D-1,
     *		   F2 = 1.36929 88092 27358 05310D-1,
     *		   F3 = 1.48753 61290 85061 48525D-2,
     *		   F4 = 7.86869 13114 56132 59100D-4,
     *		   F5 = 1.84631 83175 10054 68180D-5,
     *		   F6 = 1.42151 17583 16445 88870D-7,
     *		   F7 = 2.04426 31033 89939 78564D-15)
C	HASH SUM EF    47.52583 31754 92896 71629
C
	IFAULT = 0
	Q = P - HALF
	IF (ABS(Q) .LE. SPLIT1) THEN
	  R = CONST1 - Q * Q
	  PPND16 = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3)
     *			* R + A2) * R + A1) * R + A0) /
     *		      (((((((B7 * R + B6) * R + B5) * R + B4) * R + B3)
     *			* R + B2) * R + B1) * R + ONE)
	  RETURN
	ELSE
	  IF (Q .LT. ZERO) THEN
	    R = P
	  ELSE
	    R = ONE - P
	  END IF
	  IF (R .LE. ZERO) THEN
	    IFAULT = 1
	    PPND16 = ZERO
	    RETURN
	  END IF
	  R = SQRT(-LOG(R))
	  IF (R .LE. SPLIT2) THEN
	    R = R - CONST2
	    PPND16 = (((((((C7 * R + C6) * R + C5) * R + C4) * R + C3)
     *			* R + C2) * R + C1) * R + C0) /
     *		     (((((((D7 * R + D6) * R + D5) * R + D4) * R + D3)
     *			* R + D2) * R + D1) * R + ONE)
	  ELSE
	    R = R - SPLIT2
	    PPND16 = (((((((E7 * R + E6) * R + E5) * R + E4) * R + E3)
     *			* R + E2) * R + E1) * R + E0) /
     *		     (((((((F7 * R + F6) * R + F5) * R + F4) * R + F3)
     *			* R + F2) * R + F1) * R + ONE)
	  END IF
	  IF (Q .LT. ZERO) PPND16 = - PPND16
	  RETURN
	END IF
	END

c ----------------------------------------------------------------------
c Numerical integration. From Numerical Recipes.
c Returns as s the integral of the function func from a to b. 
c The parameters EPS can be set to the desired fractional accuracy 
c and JMAX so that 2 to the power JMAX-1 is the maximum allowed number 
c of steps. Integration is performed by Simpson’s rule.
c USES trapzd
c Converted to real*8 and obsolete PAUSE statement replaced, JKW 15/9/12
c ----------------------------------------------------------------------
      SUBROUTINE QSIMP(XCHSQ,A,B,S)
      real a,b,XCHSQ,s,eps,os,ost,st
      integer jmax,j
      external XCHSQ
      PARAMETER (EPS=1.E-6, JMAX=20)
      OST=-1.E30
      OS= -1.E30
      DO J=1,JMAX
        CALL TRAPZD(XCHSQ,A,B,ST,J)
        S=(4.*ST-OST)/3.

        if(j.gt.5)then
          if (abs(s-os).lt.EPS*abs(os).or.(s.eq.0..and.os.eq.0.)) return
        endif
        os=s
        ost=st
      end do
      WRITE(6,*)' Too many steps'
      end

c          IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
c        OS=S
c        OST=ST
c11    CONTINUE
c      PAUSE 'Too many steps.'
c      WRITE (*,*) 'Too many steps.'
c      READ (*,'()')
c      END

c ----------------------------------------------------------------------------
c From Numerical Recipes
c This routine computes the nth stage of refinement of an extended trapezoidal 
c rule. func is input as the name of the function to be integrated between 
c limits a and b, also input. When called with n=1, the routine returns as s 
c the crudest estimate of Int^b_a f(x)dx. Subsequent calls with n=2,3,... (in 
c that sequential order) will improve the accuracy of s by adding 2n-2
c additional interior points. s should not be modified between sequential calls.
c Converted to real*8, JKW 15/9/12
c -----------------------------------------------------------------------------
      SUBROUTINE TRAPZD(XCHSQ,A,B,S,N)
c JKW change:
      real a,b,s,XCHSQ
      integer n,it,j
      external XCHSQ
c      real*8 func,a,b,s,tnm,del,x,sum
      IF (N.EQ.1) THEN
        S=0.5*(B-A)*(XCHSQ(A)+XCHSQ(B))
        IT=1
      ELSE
c Original:
        TNM=IT
c JKW change:
c        TNM=REAL(IT)
        DEL=(B-A)/TNM
        X=A+0.5*DEL
        SUM=0.
        DO 11 J=1,IT
          SUM=SUM+XCHSQ(X)
          X=X+DEL
11      CONTINUE
        S=0.5*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END

c ------------------------------------------------------------------
c Function to be integrated to get expectation value for normally
c distributed variable. See e.g. equation 4.35, Julian King's thesis
c JKW 15/9/12
c ------------------------------------------------------------------
c      DOUBLE PRECISION FUNCTION XCHSQ(x)
      function XCHSQ(x)
c      real*8 x,pi
      real*4 x,pi
c      parameter(pi=3.141592653589793238462643d0)
      parameter (pi=3.14159265)
c      write(6,*)' Inside XCHSQ'
      xchsq=exp((-x**2)/2.)*(x**2)/sqrt(2*pi)
c      write(6,*)' xchsq = ',xchsq
      return
      end
