       IMPLICIT NONE 
       INTEGER i,j,k,l,m,n,sig
       REAL*8 semi1,ecc1,incl1,aper1,lasc1
       REAL*8 semi2,ecc2,incl2,aper2,lasc2,moid(16)
       REAL*8 p1(3),p2(3),q1(3),q2(3),s1(3),s2(3)
       REAL*8 aa,bb,cc,mm,nn,kk,alpha1,alpha2,u1(16),u2(16),gg,u
       REAL*8 su,cu,ac,ab,bc,a2,b2,c2,err2,cf,e2t
       REAL*8 sumgi,sumgr,t0,t1,t2
       REAL*8 pi,dot,cons,derv,func,eps,x,y,dist,carg,dd
       REAL*8 csu,ssu,mnsc,kksc,epssig,p1p2,s1s2,p1s2,s1p2
       REAL*8 cw1,co1,ci1,so1,sw1,si1,cw2,co2,ci2,so2,sw2,si2
       COMPLEX*8 m1,m2,m3,m4,c8,eps2,roots(16),nroots(16)
       COMPLEX*8 tt,c_exp,coeff(17),zz
       COMPLEX*8 sumgg,c(-8:8)
       EXTERNAL dot,c_exp,derv,func,carg,e2t
       LOGICAL changed 

       CALL CPU_TIME(T0)
       pi = 4.d0*atan(1.d0)
       
       semi1 = 1.000448939654d0
       ecc1  = 1.7118738644d-2
       incl1 = 4.18817012d-4
       aper1 = 3.26726215710d+2 
       lasc1 = 1.35082830697d+2


c       semi2 = 3.620686838118346d0
c       ecc2  = 4.328758383063912d-1
c       incl2 = 1.911426464401879d+1
c       aper2 = 2.435185259322855d+1
c       lasc2 = 3.268532150891524d2

       semi2 = 1.44270761033976d0
       ecc2  = 4.839260074970221d-1	
       incl2 = 4.098859089906774d0
       aper2 = 1.029160974050444d+2
       lasc2 = 2.1374973899383d+2

       incl1 = incl1 * pi / 180.d0
       aper1 = aper1 * pi / 180.d0
       lasc1 = lasc1 * pi / 180.d0

       incl2 = incl2 * pi / 180.d0
       aper2 = aper2 * pi / 180.d0
       lasc2 = lasc2 * pi / 180.d0

       p1(1) = cos(aper1)*cos(lasc1)-cos(incl1)*sin(aper1)*sin(lasc1)
       p1(2) = cos(aper1)*sin(lasc1)+cos(incl1)*sin(aper1)*cos(lasc1)
       p1(3) = sin(incl1)*sin(aper1)

       q1(1) = -sin(aper1)*cos(lasc1)-cos(incl1)*cos(aper1)*sin(lasc1)
       q1(2) = -sin(aper1)*sin(lasc1)+cos(incl1)*cos(aper1)*cos(lasc1)
       q1(3) =  sin(incl1)*cos(aper1)

       s1(1) = sqrt(1 - ecc1**2)*q1(1)
       s1(2) = sqrt(1 - ecc1**2)*q1(2)
       s1(3) = sqrt(1 - ecc1**2)*q1(3)
       alpha1 = semi1/semi2

       p2(1) = cos(aper2)*cos(lasc2)-cos(incl2)*sin(aper2)*sin(lasc2)
       p2(2) = cos(aper2)*sin(lasc2)+cos(incl2)*sin(aper2)*cos(lasc2)
       p2(3) = sin(incl2)*sin(aper2)

       q2(1) = -sin(aper2)*cos(lasc2)-cos(incl2)*cos(aper2)*sin(lasc2)
       q2(2) = -sin(aper2)*sin(lasc2)+cos(incl2)*cos(aper2)*cos(lasc2)
       q2(3) =  sin(incl2)*cos(aper2)

       s2(1) = sqrt(1 - ecc2**2)*q2(1)
       s2(2) = sqrt(1 - ecc2**2)*q2(2)
       s2(3) = sqrt(1 - ecc2**2)*q2(3)
       alpha2 = semi2/semi1

! Explicit formulas for C_8: 
      p1p2=dot(p1,p2)
      s1s2=dot(s1,s2)
      s1p2=dot(s1,p2)
      p1s2=dot(p1,s2)
      m1=complex(p1p2-s1s2-ecc1*ecc2,-(s1p2+p1s2))
      m2=complex(p1p2-s1s2+ecc1*ecc2,-(s1p2+p1s2))
      m3=complex(p1p2+s1s2-ecc1*ecc2,-(s1p2-p1s2))
      m4=complex(p1p2+s1s2+ecc1*ecc2,-(s1p2-p1s2))
      c8=m1*m2*m3*m4*(alpha1*ecc1**2/16.d0)**2
  
! Calculate the coefficients C_k via DFT
c      print*,'********* Printing calculated Coefficientc C_k:'
       do k = -8,8
c         sumgg=cmplx(0.d0,0.d0)
         cons=2.d0*pi/21.d0
         sumgr=0.d0
         sumgi=0.d0
         do j=-10,10
            u=real(j)*cons
            su = sin(u)
            cu = cos(u)
            aa=p1s2*su - s1s2*cu
            bb=p1p2*su - s1p2*cu
            cc=ecc2*bb - alpha1*ecc1*su*(1 - ecc1*cu)
            mm=p1p2*(cu-ecc1)+s1p2*su+alpha2*ecc2
            nn=p1s2*(ecc1-cu)-s1s2*su
            kk=alpha2*ecc2*ecc2
            a2=aa*aa
            b2=bb*bb
            c2=cc*cc
            ac=a2-c2
            ab=a2+b2
            bc=b2-c2
            gg=ac*bc*kk*kk + 2*kk*cc*(nn*aa*ac + mm*bb*bc)
     %           - ab*(ac*nn*nn + bc*mm*mm - 2*nn*mm*aa*bb)
c            sumgg = sumgg+gg*c_exp(-ii*k*u)
            sumgr = sumgr + gg*cos(real(k)*u)
            sumgi = sumgi - gg*sin(real(k)*u)     
         end do
         sumgr=sumgr
         sumgi=sumgi
         sumgg=cmplx(sumgr,sumgi)
         c(k) = sumgg*cons*3.342257109/21.d0
       end do
c       call cpu_time(t2)
c       print*,t2-t0
      k = 0
c        write(6,"(/,15x,a13,15x,a21)")'EXPLICIT C_8:','DFT CALC. C_8:'
        do i = -8,8
          k = k + 1
          coeff(k) = c(i) 
c          write(6,"(1p,'[',e14.7,',',1x,e14.7,'i',']',
c     %   5x,1p,'[',e14.7,',',1x,e14.7,'i',']')")c8,coeff(k)
      end do
c        print*,'c8-c_8=',c8-coeff(17)
      err2=(abs(c8-coeff(17))**2)
c      print*,'err=',err2,sqrt(err2)
      zz=cmplx(0.0d0,0.0d0)

c Given K complex coefficients C_k = Coeff, find the roots of the 
c polynomial P(Z), where Z is a complex number. This uses the 
c LAGUERRE's method.
       call zroots(coeff,16,roots,.true.)
c        call croots(16,coeff,zz,roots)
c        print*,abs(roots)
c Select only complex roots which lie on the unit circle.
c Here you can adjust the epsilon EPS:
        eps = 1.d-6
        j = 0
        do i=1,size(roots)
c          print*,abs(roots(i))
          if (abs(1.d0-abs(roots(i))).le.eps) then
            j = j + 1
            nroots(j)=roots(i)
          end if 
        end do

c Check the sign of m
c        write(6, "('-------- PERFORMING SIGN CHECK FOR (m) ------',/) ")
        epssig=1.d-6
	do i = 1, j
          u1(i) = carg(nroots(i))    ! calculate the (phase) of root
          aa=dot(p1,s2)*sin(u1(i)) - dot(s1,s2)*cos(u1(i))
          bb=dot(p1,p2)*sin(u1(i)) - dot(s1,p2)*cos(u1(i))
          cc=(ecc2*bb)- alpha1*ecc1*sin(u1(i))* 
     %         (1 -ecc1*cos(u1(i)))
          mm=dot(p1,p2)*cos(u1(i)) + dot(s1,p2)*sin(u1(i)) + 
     %    alpha2*ecc2 - dot(p1,p2)*ecc1 
          nn=dot(p1,s2) * ecc1 - dot(s1,s2)*sin(u1(i)) -
     %         dot(p1,s2)*cos(u1(i))
          kk=alpha2 * (ecc2**2)
          dd=aa*aa + bb*bb - cc*cc
          sig=1.d0
          changed = .false.
  15      csu=(bb*cc + sig*aa*sqrt(dd)) / (aa*aa + bb*bb)
          ssu=(aa*cc - sig*bb*sqrt(dd)) / (aa*aa + bb*bb)
          mnsc=mm*ssu + nn*csu
          kksc=kk*ssu*csu
          if ((abs(mnsc-kksc).gt.epssig) .and. (.not.changed)) then
            sig= -1.d0
            changed=.true.
            goto 15
          end if
c          u2(i)=atan2(ssu,csu)
          u2(i)=asin(ssu)
c          write(6,"('ABS(mnsc-kksc)= ', 1p,e18.10,10x, 'u2= ',
c     %           1p,e13.6,/)")abs(mnsc-kksc),u2(i)
	 end do
       write(6,"(/,'#########################################',/)")
       write(6,"('#',5x,'TOTAL OF ',i1,1x,' Proximities Detected')") J
       write(6,"(/,'#########################################',/)")
        do k = 1,j 
          moid(i) = (alpha1+alpha2)/2.d0 + 
     %       (alpha1*ecc1**2+alpha2*ecc2**2)/4.d0 -
     %       dot(p1,p2)*ecc1*ecc2 + 
     %       (dot(p1,p2)*ecc2 -alpha1*ecc1)*cos(u1(k)) + 
     %       dot(s1,p2)*ecc2*sin(u1(k)) +
     %       (dot(p1,p2)*ecc1-alpha2*ecc2)*cos(u2(k)) +
     %       dot(p1,s2)*ecc1*sin(u2(k)) - 
     %       dot(p1,p2)*cos(u1(k))*cos(u2(k)) - 
     %       dot(p1,s2)*cos(u1(k))*sin(u2(k)) -
     %       dot(s1,p2)*sin(u1(k))*cos(u2(k)) -
     %       dot(s1,s2)*sin(u1(k))*sin(u2(k)) +
     %       (alpha1*ecc1**2)*cos(2*u1(k))/4.d0 + 
     %       (alpha2*ecc2**2)*cos(2*u2(k))/4.d0
          write(6,"(1p,'MOID=',E14.7,2X,'EAN1=',E14.7,2X,
     %          'EAN2=',E16.7,/)")sqrt(2.d0*semi1*semi2*moid(i)),
     %           e2t(u1(k),ecc1)*180.d0/pi,e2t(U2(K),ecc2)*180.d0/pi
        end do
        call cpu_time(t1)
        print*,"PROGRAM RAN FOR:",T1-T0
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
       COMPLEX*8 Z

	 CARG = ATAN2(AIMAG(Z),REAL(Z))
	 END
C======================================================================
       FUNCTION c_exp(cpx)
       COMPLEX*8 c_exp
       COMPLEX*16 cpx
       real*8 a,b
       a=real(cpx)
       b=aimag(cpx)
       c_exp=cmplx(exp(a)*cos(b),exp(a)*sin(b))

	 END
C======================================================================
       SUBROUTINE laguer(a,m,x,its)
       INTEGER m,its,MAXIT,MR,MT
       REAL*8 EPSS
       COMPLEX*8 a(m+1),x
       PARAMETER (EPSS=1.d-7,MR=8,MT=100,MAXIT=MT*MR)
       INTEGER iter,j
       REAL*8 abx,abp,abm,err,frac(MR)
       COMPLEX*8 dx,x1,b,d,f,g,h,sq,gp,gm,g2
       SAVE frac
       DATA frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1.d0/
       do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.d0,0.d0)
        f=cmplx(0.d0,0.d0)
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
          if (max(abp,abm).gt.0.d0) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.d0+abx),float(iter)))
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
       print*,'too many iterations in laguer'
       return
       END

        SUBROUTINE zroots(a,m,roots,polish)  
        INTEGER m,MAXM  
        REAL*8 EPS  
        COMPLEX*8 a(m+1),roots(m)  
        LOGICAL polish  
        PARAMETER (EPS=1.d-7,MAXM=100)  
CU    USES laguer
      	
        INTEGER i,j,jj,its  
        COMPLEX*8 ad(MAXM),x,b,c  
        do 11 j=1,m+1  
           ad(j)=a(j)  
  11    continue  
        do 13 j=m,1,-1
           x=cmplx(0.d0,0.d0)  
           call laguer(ad,j,x,its)  
           if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),
     %          0.d0) 
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


      SUBROUTINE CROOTS(DG,C,Z,ZR)
      implicit none

C   'DG'  - THE DEGREE OF THE COMPLEX POLYNOMIAL
c   'K '  - HOW MANY ITERATIONS BEFORE ROOT CONVERGES
C   'ITER'- MAXIMUM NUMBER OF ITERATIONS FOR CONVERGING
C   'DM1' - (DG-1)
C   'P'   - COUNTER FOR THE # OF ROOTS
C   'Q'   - INT.USED TO DECREASE THE # OF COEFF. ONCE A ROOT HAS BEEN
C           ISOLATED, THE # INDEX OF COEFF. IS REDUCED BY 1.

C    'EPS' - MAXIMUM PRECISION OF THE ROOTS 


C    'ZR' - STORES THE COMPLEX ROOTS. HAS DIMENSION OF 'DG'.
C    'Z ' - INIT. GUESS AND TEMP. UPDATE OF THE SOUGHT ROOT.
C    'C ' - COMPLEX COEFFICIENTS OF DEGREE (DG+1). THE POLYNOMIAL IS
C           OF THE FORM: C_1*Z^N + C_2*Z^(N-1)+...+C_(N+1)*Z^(N-N)=0 .
C           E.G.: c1*Z^4 + c2*Z^3 + c3*Z^2 + c4*Z + c5=0 .
c    'F'  - THE COMPLEX POLYNOMIAL OF THE FORM AS INDICATED ABOVE.
C    'DF1'- FIRST DERIVATIVE OF THE POLYNOMIAL
C    'DF2'- SECOND DERIVATIVE OF THE POLYNOMIAL
      INTEGER b,d,dg,dm1,dp1,i,j,k,l,m,n,p,q,iter
      parameter(iter=100)
      REAL*16 im1,im2,re1,re2,i12,r12,eps,delta
      PARAMETER(EPS=1.D-8)
      COMPLEX*16 f,df1,df2,zo,z,zn,c(DG+1),zr(DG)
c=======================================================================

      dm1=dg-1 ! DEGREE - 1
      dp1=dg+1 ! DEGREE + 1
      
c      c(1)=cmplx(5.d0,1.d0)
c      c(2)=cmplx(-1.d0,3.d0)
c      c(3)=cmplx(1.d0,1.d0)
c      c(4)=cmplx(-1.d0,0.d0)

c      dg=3
      p=0
      q=0
      do n=dg,1,-1    ! loop backwards over the degree of polynomial
        p=p+1
        k=0
        if (n.eq.1) then              ! THIS IS FOR THE LAST ROOT ONLY,
           zr(p)=-(c(2)+zr(p-1)*c(1))/c(1) ! IT IS FOUND ANALYTICALY
           print*,'root',p,zr(p)
           exit
        end if
        do i = 1, iter
          f=cmplx(0.d0,0.d0)
          df1=cmplx(0.d0,0.d0)
          df2=cmplx(0.d0,0.d0)
          k=k+1
          do j=n,0,-1     ! calculate 'f' and its derivatives 
            l=n+1-j
            f=f+c(l)*z**j       ! e.g. c1*z^3 + c2*z^2 + c3*z
            if (j.ne.0) df1=df1+j*c(l)*z**(j-1)
            if (j.gt.1) df2=df2+j*c(l)*z**(j-2)
          end do
          zn=z-(f/df1)
c          zn=z-2.d0*f*df1/(2.d0*df1**2-f*df2)  ! Halley's method
          delta=(abs(z)**2-abs(zn)**2)/abs(zn)**2
          z=zn
          if(abs(delta).le.eps) then
            exit
          else if (i.eq.iter.and.(abs(delta).gt.eps)) then
             print*,"Solution exceeded maximum iterations!"
             exit
          end if
c          z=zn
c          print*,i
        end do
        zr(p)=z         ! get current best root and store it in 'zr'.
c        print*,p,zr(p)
c       z=cmplx(0.d0,0.d0)
        z=1.d0/conjg(zr(p)) ! initiate guess for the next root, z=1/z* 
        q=q+1               ! is a suggestion by Balyuev (2018)
        if (p.lt.dm1) then
          do j=2,dp1-p
            c(j)=c(j)+c(j-1)*zr(p) ! UPDATE THE COEFFICIENTS ACCORDING
          end do                   ! TO THE HORNER'S SCHEME
        end if
      end do
      END

      SUBROUTINE HORNER(D,RTIN,CFIN,CFOUT)

              INTEGER I,J,K,D
              REAL*16 EPS
              PARAMETER(EPS=1.D-10)
              COMPLEX*16 RTIN,CFIN(D),CFOUT(D-1),TEMP(D)

              TEMP(1)=CFIN(1)

              DO I=2,D
                 TEMP(I)=CFIN(I)+RTIN*CFIN(I-1)
                 IF (I.EQ.D) THEN
                    IF (ABS(CFIN(D)-TEMP(D)).GT.EPS) THEN
                      WRITE(6,"('WARNING:IN HORNER:',1P,2(E16.7,E16.7)
     %              ,2X, 'MAY NOT BE A ROOT OF THE PLYNOMIAL',/)")RTIN
                      EXIT
                    END IF
                    EXIT
                 END IF
              END DO
              DO J=1,D-1
                CFOUT(J)=TEMP(J)
              END DO
              END
      FUNCTION E2T(EAN,ECC)
            REAL*8 EAN,E2T,ECC

            E2T=2.D0*ATAN(SQRT((1.d0+ECC)/(1.d0-ECC))*TAN(EAN/2.D0))

            RETURN
            END
