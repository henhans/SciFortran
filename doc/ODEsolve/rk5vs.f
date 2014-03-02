      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     &  derivs,rkqs)
c
c .................................................................
c .   Subroutine odeint 
c .   Runge-Kutta driver with adaptive stepsize control.  Integrate
c .   the starting values ystart(1:nvar) from x1 to x2 with accuracy 
c .   eps, storing intermediate results in the common block /path/.
c .   h1 should be set as a guessed first stepsize, hmin as the 
c .   minimum allowed stepsize (can be zero).
c .   On output, nok and nbad are the number of good and bad (but
c .   retried and fixed) steps taken, and ystart is replaced by 
c .   values at the end of the integration interval.
c .   derivs is the user-supplied subroutine for calculating the
c .   right-hand side derivatives, while rkqs is the name of 
c .   the stepper routine to be used.  
c .   /path/ contains its own information about how often an
c .   intermediate value is to be stored.
c .................................................................
c
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL*8 eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,rkqs
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.0d-30)
      INTEGER i,kmax,kount,nstp
      REAL*8 dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     & yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do i=1,nvar
        y(i)=ystart(i)
      end do
      if (kmax.gt.0) xsav=x-2.*dxsav
      do nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
        end do
        if(kmax.gt.0) then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1) then
              kount=kount+1
              xp(kount)=x
              do i=1,nvar
                yp(i,kount)=y(i)
              end do
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1) .gt. 0.0d0) h=x2-x
        call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h) then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1) .ge. 0.0d0)then
          do i=1,nvar
            ystart(i)=y(i)
          end do
          if(kmax.ne.0) then
            kount=kount+1
            xp(kount)=x
            do i=1,nvar
              yp(i,kount)=y(i)
            end do
          endif
          return
        endif
        if(abs(hnext) .lt. hmin) pause
     &     'stepsize smaller than minimum in odeint'
        h=hnext
        end do
      pause 'too many steps in odeint'
      return
      END
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs)
c
c ................................................................
c .   Subroutine rkqs is a fifth-order Runge-Kutta step with 
c .   monitoring of local truncation error to ensure accuracy
c .   and adjust stepsize.   
c .
c .   Inputs:
c .   y(1:n) = the dependent variable vector
c .   dydx(1:n) = deriv of y wrt x at starting value of indep var x
c .   htry = attempted step size
c .   eps = required accuracy
c .   yscal(1:n) = scaling vector for the error
c .   
c .   Outputs:
c .   y = new value of y
c .   x = new value of x
c .   hdid = actual stepsize accomplished
c .   hnext = estimated next stepsize   
c .
c .   derivs is the user-supplied subroutine that computes the 
c .    right-hand side derivatives
c ................................................................
c
      INTEGER n,NMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
c     Maximum number of equations
      INTEGER i
      REAL*8 errmax,h,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,PGROW,
     &  PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      h=htry
    1 call rkck(y,dydx,n,x,h,ytemp,yerr,derivs)
      errmax=0.0d0
      do i=1,n
        errmax=max(errmax,abs(yerr(i)/yscal(i)))
      end do
      errmax=errmax/eps
      if(errmax.gt.1.0d0) then
        htemp=SAFETY*h*(errmax**PSHRNK)
        h=sign(max(abs(htemp), 0.1*abs(h)), h)
        xnew=x+h
        if(xnew.eq.x) pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON) then
          hnext=SAFETY*h*(errmax**PGROW)
        else
          hnext=5.0d0*h
        endif
        hdid=h
        x=x+h
        do i=1,n
          y(i)=ytemp(i)
        end do
        return
      endif
      END
      SUBROUTINE rkck(y,dydx,n,x,h,yout,yerr,derivs)
c
c ...............................................................
c .   Subroutine rkck.
c .   Given values for n variables y and their derivatives dydx
c .   known at x, use the fifth-order Cash-Karp Runge-Kutta method
c .   to advance the solution over an interval h and return the 
c .   incremented variables as yout.  Also return an estimate of
c .   the local truncation error in yout using the embedded 
c .   fourth-order method.  The user supplies the subroutine
c .   derivs (x,y,dydx) which returns derivatives dydx at x.
c ................................................................
c
      INTEGER n,NMAX
      REAL*8 h,x,dydx(n),y(n),yerr(n),yout(n)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
C     USES derivs
      INTEGER i
      REAL*8 ak2(NMAX),ak3(NMAX),ak4(NMAX),ak5(NMAX),ak6(NMAX),
     &  ytemp(NMAX),A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     &  B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,
     &  DC5,DC6
      PARAMETER (A2=0.2,A3=0.3,A4=0.6,A5=1.0,A6=0.875,B21=0.2,
     &  B31=3.0/40.0,B32=9.0/40.0,B41=0.3,B42=-0.9,B43=1.2,
     &  B51=-11.0/54.0,B52=2.5,B53=-70.0/27.0,B54=35.0/27.0,
     &  B61=1631.0/55296.0,B62=175.0/512.0,B63=575.0/13824.0,
     &  B64=44275.0/110592.0,B65=253.0/4096.0,C1=37.0/378.0,
     &  C3=250.0/621.0,C4=125.0/594.0,C6=512.0/1771.0,
     &  DC1=C1-2825.0/27648.0,DC3=C3-18575.0/48384.0,
     &  DC4=C4-13525.0/55296.0,DC5=-277.0/14336.0,DC6=C6-0.25)
      do i=1,n
        ytemp(i)=y(i)+B21*h*dydx(i)
      end do 
      call derivs(x+A2*h,ytemp,ak2)
      do i=1,n
        ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      end do
      call derivs(x+A3*h,ytemp,ak3)
      do i=1,n
        ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      end do
      call derivs(x+A4*h,ytemp,ak4)
      do i=1,n
        ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+B54*ak4(i))
      end do
      call derivs(x+A5*h,ytemp,ak5)
      do i=1,n
        ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+B64*ak4(i)+
     &     B65*ak5(i))
      end do
      call derivs(x+A6*h,ytemp,ak6)
      do i=1,n
        yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+C6*ak6(i))
      end do
c     estimate error as difference between fourth and fifth order methods
      do i=1,n
        yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)+
     &   DC6*ak6(i))
      end do
      return
      END
