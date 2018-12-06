       IMPLICIT NONE

       INTEGER i,j,k,l,m,n,sig
       REAL*8 semi1,ecc1,incl1,aper1,lasc1
       REAL*8 semi2,ecc2,incl2,aper2,lasc2,moid(16)
       REAL*8 p1(3),p2(3),q1(3),q2(3),s1(3),s2(3)
       REAL*8 aa,bb,cc,mm,nn,kk,alpha1,alpha2,u1(16),u2(16),gg,u
       REAL*8 su,cu,ac,ab,bc,a2,b2,c2
       REAL*8 pi,dot,cons,derv,func,eps,x,y,dist,carg,dd
       REAL*8 csu,ssu,mnsc,kksc,epssig
       COMPLEX*16 m1,m2,m3,m4,c8,eps2,roots(16),nroots(16)
       COMPLEX*16 tt,ii,c_exp,c(-8:8),coeff(17),sumgg
       EXTERNAL dot,c_exp,derv,func,carg
       LOGICAL changed 

       open(20,file='gvsu.dat',status='unknown')
       pi = 4.0*atan(1.d0)
       ii = (0D0,1D0)
      
       semi1 = 1.00000011
       ecc1  = 1.672293696353922E-02
       incl1 = 5.E-05
       aper1 = 1.1420783E+02 
       lasc1 = 3.4873936E+02

       semi2 = 9.235960944E+00
       ecc2  = 9.189810923E-01	
       incl2 = 2.896687E+01
       aper2 = 1.960254E+02
       lasc2 = 2.506264E+02

       incl1 = incl1 * pi / 180d0
       aper1 = aper1 * pi / 180d0
       lasc1 = lasc1 * pi / 180d0

       incl2 = incl2 * pi / 180d0
       aper2 = aper2 * pi / 180d0
       lasc2 = lasc2 * pi / 180d0

       p1(1) = cos(aper1)*cos(lasc1) - cos(incl1)*sin(aper1)*sin(lasc1)
       p1(2) = cos(aper1)*sin(lasc1) + cos(incl1)*sin(aper1)*cos(lasc1)
       p1(3) = sin(incl1)*sin(aper1)

       q1(1) = -sin(aper1)*cos(lasc1) - cos(incl1)*cos(aper1)*sin(lasc1)
       q1(2) = -sin(aper1)*sin(lasc1) + cos(incl1)*cos(aper1)*cos(lasc1)
       q1(3) =  sin(incl1)*cos(aper1)

       s1(1) = sqrt(1. - ecc1**2)*q1(1)
       s1(2) = sqrt(1. - ecc1**2)*q1(2)
       s1(3) = sqrt(1. - ecc1**2)*q1(3)
       alpha1 = semi1/semi2

       p2(1) = cos(aper2)*cos(lasc2) - cos(incl2)*sin(aper2)*sin(lasc2)
       p2(2) = cos(aper2)*sin(lasc2) + cos(incl2)*sin(aper2)*cos(lasc2)
       p2(3) = sin(incl2)*sin(aper2)

       q2(1) = -sin(aper2)*cos(lasc2) - cos(incl2)*cos(aper2)*sin(lasc2)
       q2(2) = -sin(aper2)*sin(lasc2) + cos(incl2)*cos(aper2)*cos(lasc2)
       q2(3) =  sin(incl2)*cos(aper2)

       s2(1) = sqrt(1. - ecc2**2)*q2(1)
       s2(2) = sqrt(1. - ecc2**2)*q2(2)
       s2(3) = sqrt(1. - ecc2**2)*q2(3)
       alpha2 = semi2/semi1

! Explicit formulas for C_8:
      
      m1=dot(p1,p2)-dot(s1,s2)-ecc1*ecc2-ii*(dot(s1,p2)+dot(p1,s2))
      m2=dot(p1,p2)-dot(s1,s2)+ecc1*ecc2-ii*(dot(s1,p2)+dot(p1,s2))
      m3=dot(p1,p2)+dot(s1,s2)-ecc1*ecc2-ii*(dot(s1,p2)-dot(p1,s2))
      m4=dot(p1,p2)+dot(s1,s2)+ecc1*ecc2-ii*(dot(s1,p2)-dot(p1,s2))
      c8=m1*m2*m3*m4*(alpha1*ecc1**2/16.)**2

! Calculate the coefficients C_k.
      print*,'********* Printing calculated Coefficientc C_k:'
       do k = -8, 8
         sumgg=(0.d0,0.d0)
         cons=1./21.
         do j=0,20
            u=2.*pi*j*cons
            su = sin(u)
            cu = cos(u)
            aa=dot(p1,s2)*su - dot(s1,s2)*cu
            bb=dot(p1,p2)*su - dot(s1,p2)*cu
            cc=(ecc2*bb) - (alpha1*ecc1*su)*(1 - ecc1*cu)
            mm=dot(p1,p2)*(cu-ecc1)+dot(s1,p2)*su+alpha1*ecc2
            nn=dot(p1,s2)*(ecc1-cu) - dot(s1,s2)*su
            kk=alpha2*ecc2**ecc2
            a2=aa*aa
            b2=bb*bb
            c2=cc*cc
            ac=a2-c2
            ab=a2+b2
            bc=b2-c2
            gg=ac*bc*kk*kk + 2.*kk*cc*(nn*aa*ac + mm*bb*bc)
     %           - ab*(ac*nn*nn +bc*mm*mm - 2.*nn*mm*aa*bb)
            sumgg = sumgg + gg*c_exp(real(k)*u)
         end do
         c(k) = sumgg*cons
       end do 

        k = 0
        write(6,"(/,15x,a13,15x,a21)")'EXPLICIT C_8:','DFT CALC. C_8:'
        do i = -8,8
          k = k + 1
          coeff(k) = c(i) 
          write(6,"(1p,'[',e14.7,',',1x,e14.7,'i',']',
     %   5x,1p,'[',e14.7,',',1x,e14.7,'i',']')")c8,coeff(k)
      end do
        print*,'c8-c_8=',c8-coeff(17)

c Given K complex coefficients C_k = Coeff, find the roots of the 
c polynomial P(Z), where Z is a complex number. This uses the 
c LAGUERRE's method.
       call zroots(coeff,16,roots,.true.)

c Select only complex roots which lie on the unit circle.
c Here you can adjust the epsilon EPS:
        eps = 1.E-7
        j = 0
        do i = 1, size(roots)
         x = real(roots(i))
         y = aimag(roots(i))
         dist = sqrt(x**2 + y**2)
         if (abs(1.d0-dist).le.eps) then
           j = j + 1
           nroots(j)=roots(i)
         end if 
        end do

c Check the sign of m
        write(6, "('-------- PERFORMING SIGN CHECK FOR (m) ------',/) ")
        epssig=1.E-7
	do i = 1, j
          u1(i) = carg(roots(i))    ! calculate the argument(phase) of the roots
	  aa = dot(p1,s2) * sin(u1(i)) - dot(s1,s2) * cos(u1(i))
          bb = dot(p1,p2) * sin(u1(i)) - dot(s1,p2) * cos(u1(i))
          cc = (ecc2 * bb)- alpha1 * ecc1 * sin(u1(i)) * 
     %         (1 -ecc1*cos(u1(i)))
          mm = dot(p1,p2) * cos(u1(i)) + dot(s1,p2)*sin(u1(i)) + 
     %    alpha1*ecc2 - dot(p1,p2)*ecc1 
          nn = dot(p1,s2) * ecc1 - dot(s1,s2)*sin(u1(i)) -
     %         dot(p1,s2)*cos(u1(i))
          kk = alpha2 * (ecc2**2)
          dd = aa*aa + bb*bb - cc*cc
          sig = 1.0
          changed = .false.
  15      csu = (bb*cc + sig*aa*sqrt(dd)) / (aa*aa + bb*bb)
          ssu = (aa*cc - sig*bb*sqrt(dd)) / (aa*aa + bb*bb)
          u2(i) = atan2(ssu,csu)
          mnsc = mm*ssu + nn*csu
          kksc = kk*ssu*csu
          if ((dabs(mnsc-kksc).gt.epssig) .and. (.not.changed)) then
            sig = -1.0
            changed = .true.
            goto 15
          end if
          write(6,"('ABS(mnsc-kksc)= ', 1p,e18.10,10x, 'u2= ',
     %           1p,e13.6,/)")abs(mnsc-kksc),u2(i)
	end do


        do k = 1,j 
          moid(i) = (alpha1+alpha2)/2. + 
     %       (alpha1*ecc1**2+alpha2*ecc2**2)/4. -
     %       dot(p1,p2)*ecc1*ecc2 + 
     %       (dot(p1,p2)*ecc2 -alpha1*ecc1)*cos(u1(k)) + 
     %       dot(s1,p2)*ecc2*sin(u1(k)) +
     %       (dot(p1,p2)*ecc1-alpha2*ecc2)*cos(u2(k)) +
     %       dot(p1,s2)*ecc1*sin(u2(k)) - 
     %       dot(p1,p2)*cos(u1(k))*cos(u2(k)) - 
     %       dot(p1,s2)*cos(u1(k))*sin(u2(k)) -
     %       dot(s1,p2)*sin(u1(k))*cos(u2(k)) -
     %       dot(s1,s2)*sin(u1(k))*sin(u2(k)) +
     %       (alpha1*ecc1**2)*cos(2*u1(k))/4. + 
     %       (alpha2*ecc2**2)*cos(2*u2(k))/4.
          write(6,"(1p,E16.7)")moid(i)
        end do
        
        END
c====================== END OF MAIN PROGRAM ===========================
C 
C======================================================================
c =====================	Begin FUNCTIONS ===============================
C======================================================================
	FUNCTION dot(v1,v2)
        REAL*8 dot,v1(3),v2(3)
	dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
	return
	END
C======================================================================
	FUNCTION carg(z)
	REAL*8 CARG
        DOUBLE COMPLEX Z

	CARG = ATAN2(AIMAG(Z),REAL(Z))
	RETURN
	END
C======================================================================
       FUNCTION c_exp(angle)
       DOUBLE COMPLEX c_exp
       real*8 angle   ! The angle must be in radians.

       c_exp = cmplx(dcos(angle),dsin(angle))
       return
	 END
C======================================================================
       SUBROUTINE laguer(a,m,x,its)
       INTEGER m,its,MAXIT,MR,MT
       REAL*8 EPSS
       DOUBLE COMPLEX a(m+1),x
       PARAMETER (EPSS=6.e-8,MR=8,MT=200,MAXIT=MT*MR)
       INTEGER iter,j
       REAL*8 abx,abp,abm,err,frac(MR)
       DOUBLE COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
       SAVE frac
       DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
       do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
  11   continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.+abx),float(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
  12   continue
       pause 'too many iterations in laguer'
       return
       END

        SUBROUTINE zroots(a,m,roots,polish)  
        INTEGER m,MAXM  
        REAL*8 EPS  
        DOUBLE COMPLEX a(m+1),roots(m)  
        LOGICAL polish  
        PARAMETER (EPS=6.e-8,MAXM=201)  
CU    USES laguer
      	
        INTEGER i,j,jj,its  
        DOUBLE COMPLEX ad(MAXM),x,b,c  
        do 11 j=1,m+1  
           ad(j)=a(j)  
  11    continue  
        do 13 j=m,1,-1  
           x=cmplx(0.,0.)  
           call laguer(ad,j,x,its)  
           if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),
     %          0.) 
           roots(j)=x  
           b=ad(j+1)  
           do 12 jj=j,1,-1  
             c=ad(jj)  
             ad(jj)=b  
             b=x*b+c  
  12     continue  
  13     continue  
      if (polish) then  
        do 14 j=1,m  
          call laguer(a,m,roots(j),its)  
  14    continue  
      endif  
      do 16 j=2,m  
        x=roots(j)  
        do 15 i=j-1,1,-1  
          if(real(roots(i)).le.real(x))goto 10  
          roots(i+1)=roots(i)  
  15    continue  
        i=0  
  10      roots(i+1)=x  
  16    continue  
      return  
      END  

C=======================================================================
	FUNCTION derv(func,x0,h)
      INTEGER i
      REAL*8 h,derv,x0,x(4),y(4),func
      EXTERNAL func
      do i = 1, 4
         if (i .le. 2) then
           x(i) = x0-(3.-i)*h
         else
           x(i) = x0 + (i-2.)*h
         end if
         y(i) = func(x(i))
      end do
      derv = (1./(12d0*h))*(y(1) - 8d0*y(2) + 8d0*y(3) - y(4))
      return
	END
C=======================================================================
      FUNCTION func(x)
      REAL*8 func,x

c Put your custom function "func" here! Current defauls "sin(x)/x"

      func = sin(x)/x
	return
	END

c      SUBROUTINE ROOTS(KNOWN,NMORE,KM,REALRT,EP1,EP2,GUESS,MAXITS,RTS,NFOUND, FNV, IFAIL)

c      COMPLEX*16 RTS(KM), FNV(KM), X1, X2, X3, F1, F2, F3, &
c      X21, X32, F21, F32, F31, F321, B, XNEW, FNEW, DENOM, &
c      FSAVE, CSQRT, RADICL, RT, DIF, CSTEP, FSLAST
c      REAL*8 EP1, EP2, ONE, HUNDRD, TWO, FOUR, HALF, STEP, &
c      ZERO, EP1DEF, EP2DEF
c      INTEGER GUESS
c      LOGICAL REALRT
c      DATA IZERO,IONE,ITWO,ITHREE,ZERO,HALF,ONE,TWO,FOUR,HUNDRD,ITSDEF,EP1DEF,EP2DEF &
c      /0,1,2,3,0.0d0,0.5d0,1.0d0,2.0d0,4.0d0,100.0d0,100,0.5d-5,1.0d-6/
c      EP1 = AMAX1(EP1,EP1DEF)
c      EP2 = AMAX1(EP2,EP2DEF)
c      MAXITS = MIN0(MAXITS,ITSDEF)
c      IFAIL = IZERO
c      NFOUND = IZERO
c      IF (KNOWN.LT.IONE) GO TO 30
c      DO 10 I=1,KNOWN
c      II = I
c      CALL CALCF(RTS(I), FNV(I))
c      IF (ABS(FNV(II)).GE.EP2) GO TO 20
c      10 CONTINUE
c      20 IFAIL =IONE
c      GUESS = GUESS + KNOWN - II +IONE
c      NMORE = NMORE + KNOWN - II +IONE
c      KNOWN = II - IONE
c      30 CONTINUE
c      LOOP1 = KNOWN +IONE
c      LOOP2 = KNOWN + NMORE
c      IF (LOOP1.GT.LOOP2 .OR. LOOP1.LT.IONE) GO TO 180
c      IF (GUESS-NMORE) 40, 70, 60
c      40 ILO = GUESS + IONE
c      DO 50 I=ILO,LOOP2
c      RTS(I) = ZERO
c      50 CONTINUE
c      GO TO 70
c      60 CONTINUE
c      GUESS = NMORE
c      70 CONTINUE
c      STEP = HALF
c      DO 160 NEW=LOOP1,LOOP2
c      KOUNT = ITHREE
c      NEWM1 = NEW - IONE
c      RT = RTS(NEW)
c      X1 = RT - STEP
c      X2 = RT + STEP
c      X3 = RT
c      CALL CALCF(X1, F1)
c      CALL CALCF(X2, F2)
c      CALL CALCF(X3, F3)
c      FSLAST = F3
c      IF (NEW.GT.IONE) CALL TEST(X1, F1, FSAVE, RTS,NEWM1, EP2, KOUNT)
c      IF (NEW.GT.IONE) CALL TEST(X2, F2, FSAVE, RTS,NEWM1, EP2, KOUNT)
c      IF (NEW.GT.IONE) CALL TEST(X3, F3, FSLAST, RTS,NEWM1, EP2, KOUNT)
c      F21 = (F2-F1)/(X2-X1)
c      80 X32 = X3 - X2
c      F32 = (F3-F2)/X32
c      F321 = (F32-F21)/(X3-X1)
c      B = F32 + X32*F321
c      RADICL = B**ITWO - FOUR*F321*F3
c      IF (REALRT .AND. REAL(RADICL).LT.ZERO) RADICL = ZERO
c      RADICL = SQRT(RADICL)
c      IF (REAL(B)*REAL(RADICL)+AIMAG(B)*AIMAG(RADICL).LT.ZERO) RADICL = -RADICL
c      DENOM = B + RADICL
c      IF (ABS(DENOM).NE.ZERO) GO TO 100
c      IF (ABS(F3).GE.EP2) GO TO 90
c      XNEW = X3
c      GO TO 120
c      90 XNEW = X3 + X32
c      GO TO 120
c      100 CSTEP = TWO*F3/DENOM
c      IF (.NOT.REALRT .OR. ABS(F3).EQ.ZERO .OR. ABS(F32).EQ.ZERO) GO TO 110
c      CSTEP = F32/ABS(F32)*F3/ABS(F3)*ABS(CSTEP)
c      110 XNEW = X3 - CSTEP
c      120 CALL CALCF(XNEW, FNEW)
c      FSAVE = FNEW
c      IF (NEW.LE.IONE) GO TO 130
c      CALL TEST(XNEW, FNEW, FSAVE, RTS, NEWM1, EP2,KOUNT)
c      130 KOUNT = KOUNT +IONE
c      IF (KOUNT.GT.MAXITS) GO TO 170
c      DIF = XNEW - X3
c      IF (ABS(DIF).LT.EP1*ABS(XNEW)) GO TO 150
c      IF (ABS(FSAVE).LT.EP2) GO TO 150
c      IF (REALRT .OR. ABS(FSAVE).LT.HUNDRD*ABS(FSLAST)) GO TO 140
c      CSTEP = CSTEP*HALF
c      XNEW = XNEM + CSTEP
c      GO TO 120
c      140 X1 = X2
c      X2 = X3
c      X3 = XNEW
c      F1 = F2
c      F2 = F3
c      F3 = FNEW
c      FSLAST = FSAVE
c      F21 = F32
c      GO TO 80
c      150 CONTINUE
c      RTS(NEW) = XNEW
c      FNV(NEW) = FSAVE
c      NFOUND = NEW
c      160 CONTINUE
c      GO TO 190
c      170 CONTINUE
c      RTS(NEW) = XNEW
c      FNV(NEW) = FSAVE
c      IFAIL = ITHREE
c      GO TO 190
c      180 CONTINUE
c      IFAIL = ITWO
c      190 CONTINUE
c      RETURN
c      END

c      SUBROUTINE TEST(X, F, FSAVE, RTS, K, EPS, KOUNT)
c      COMPLEX*16 X, F, RTS(K), D, FSAVE
c      REAL*8 EPS!, CABS
c      DATA IONE, PERTB /1,0.01d0/
c      10 CONTINUE
c      DO 20 I=1,K
c      D = X - RTS(I)
c      IF (ABS(D).LE.EPS) GO TO 30
c      F = F/D
c      20 CONTINUE
c      GO TO 40
c      30 CONTINUE
c      X = X + PERTB
c      CALL CALCF(X, F)
c      FSAVE = F
c      KOUNT = KOUNT +IONE
c      GO TO 10
c      40 RETURN
c      END
