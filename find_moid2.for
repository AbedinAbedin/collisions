      implicit none
c
      INTEGER I,J,K,D,P
      REAL*8 A1,E1,I1,L1,W1
      REAL*8 A2,E2,I2,L2,W2,moid
      REAL*8, DIMENSION(1000)::W,OM,INC,EC,A
      REAL*8 D2R,PI,emoid,t0,t1

      PI=4.D0*ATAN(1.D0)
      D2R=PI/180.D0
c
      OPEN(20,FILE='PHAs.txt',STATUS='OLD')
c      call cpu_time(t0)
      READ(20,*)
      K=0
      DO
        K=K+1
        READ(20,*,END=666)EMoid,W(K),OM(K),INC(K),EC(K),A(K)
      END DO

  666 K=K-1
      print*,k
      D=1
      P=0
  10    A1=A(D)
        E1=EC(D)
        I1=INC(D)*D2R
        L1=OM(D)*D2R
        W1=W(D)*D2R
        DO J=D+1,K
          A2=A(J)
          E2=EC(J)
          I2=INC(J)*D2R
          L2=OM(J)*D2R
          W2=W(J)*D2R
          CALL MINDIST(A1,E1,I1,L1,W1,
     %                 A2,E2,I2,L2,W2,MOID)
          P=P+1
c          if (isnan(moid)) print*,a1,e1,a2,e2
          write(6,"(1P,'MOID=',E14.7,' [AU]')")MOID
        END DO
        D=D+1
        IF (D.LT.K) GOTO 10
      PRINT*,P
c      CALL CPU_TIME(T1)
c      PRINT*,(T1-T0)
      END

       SUBROUTINE MINDIST(SEMI1,ECC1,INCL1,LASC1,APER1,
     %                    SEMI2,ECC2,INCL2,LASC2,APER2,MOID)

c This code still needs to be fixed
c Last added calculation of the true anomaly U1, at the node of
c the first orbit. Have to check if the computed moid is accurate.

       IMPLICIT NONE
       INTEGER i,j,k,l,m,n,sig
       REAL*8 semi1,ecc1,incl1,aper1,lasc1
       REAL*8 semi2,ecc2,incl2,aper2,lasc2,moid,dist(16)
       REAL*8 p1(3),p2(3),q1(3),q2(3),s1(3),s2(3)
       REAL*8 aa,bb,cc,mm,nn,kk,alpha1,alpha2,u1(16),u2(16),gg,u,vv
       REAL*8 su,cu,ac,ab,bc,a2,b2,c2,err2,sqerr,cf,e2t
       REAL*8 sumgi,sumgr,t0,t1,t2,x_,y_,x1,x2,y1,y2
       REAL*8 pi,dot,cons,eps,x,y,carg,dd,err,err_
       REAL*8 csu,ssu,mnsc,kksc,epssig,p1p2,s1s2,p1s2,s1p2
       REAL*8 cw1,co1,ci1,so1,sw1,si1,cw2,co2,ci2,so2,sw2,si2
       REAL*8 ac_,bd_,ad_,bc_,elm1(5),elm2(5)
       COMPLEX*16 m1,m2,m3,m4,c8,eps2,roots(16),nroots(16)
       COMPLEX*16 tt,c_exp,coeff(17),zz
       COMPLEX*16 sumgg,c(-8:8)
       EXTERNAL dot,c_exp,derv,func,carg,e2t
       LOGICAL changed
c-----------------------------------------------------------------------

       pi = 4.d0*atan(1.d0)

       elm1(1)=semi1; elm1(2)=ecc1; elm1(3)=incl1
       elm1(4)=lasc1; elm1(5)=aper1;

       elm2(1)=semi2; elm2(2)=ecc2; elm2(3)=incl2
       elm2(4)=lasc2; elm2(5)=aper2;

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
       do k = -8,8
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
            sumgr = sumgr + gg*cos(real(k)*u)
            sumgi = sumgi - gg*sin(real(k)*u)
         end do
         sumgr=sumgr
         sumgi=sumgi
         sumgg=cmplx(sumgr,sumgi)
         c(k) = sumgg*cons*3.342257109/21.d0
       end do

      k = 0
      do i = -8,8
         k = k + 1
         coeff(k) = c(i)
      end do
      err2=(abs(c8-coeff(17))**2)
      sqerr=sqrt(sqrt(err2))

      call croots(p1,q1,elm1,elm2,16,coeff,roots)

c Select only complex roots which lie on the unit circle.
c Here you can adjust the epsilon EPS:
        eps = 1.d-6
        j = 0
        do i=1,size(roots)
          if (abs(1.d0-abs(roots(i))).le.eps) then
            j = j + 1
            nroots(j)=roots(i)
          end if
        end do

c Check the sign of m
      epssig=1.d-6
      do i = 1, j
          u1(i) = carg(nroots(i))    ! calculate the (phase) of root
          x=cos(u1(i))
          y=sin(u1(i))
          aa=dot(p1,s2)*y - dot(s1,s2)*x
          bb=dot(p1,p2)*y - dot(s1,p2)*x
          cc=(ecc2*bb)- alpha1*ecc1*y*(1. -ecc1*x)
          mm=dot(p1,p2)*(x-ecc1) + dot(s1,p2)*y + alpha2*ecc2
          nn=dot(p1,s2)*(ecc1-x) - dot(s1,s2)*y
          kk=alpha2*ecc2*ecc2
          vv=aa*aa + bb*bb
          dd=vv-cc*cc
          if (dd.lt.0.d0) then
            dd=0.d0
          else
            dd=sqrt(dd)
          end if
          ac_=aa*cc
          bd_=bb*dd
          ad_=aa*dd
          bc_=bb*cc
  15      x1=(bc_ + ad_)/vv
          y1=(ac_ - bd_)/vv
          x2=(bc_ - ad_)/vv
          y2=(ac_ + bd_)/vv
          err =mm*y1 + nn*x1 - kk*y1*x1
          err_=mm*y2 + nn*x2 - kk*y2*x2
c          print*,"err,err_=",err
          if (abs(err-err_).lt.epssig) then
            u2(i)=atan2(y1,x1)
          else
            u2(i)=atan2(y2,x2)
          end if
      end do
c
      do k = 1,j
          dist(k) = (alpha1+alpha2)/2.d0 +
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
        end do
        moid=sqrt(minval(dist(1:j))*2.*semi1*semi2)
c        write(6,"(/,1p,'MOID=',E13.6,/)")moid
        END
c====================== END OF MAIN PROGRAM ===========================
C
c======================================================================
c ===================== Begin FUNCTIONS ===============================
C======================================================================
      FUNCTION dot(v1,v2)
       REAL*8 dot,v1(3),v2(3)
      dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
      return
      END
C======================================================================
      FUNCTION carg(z)
      REAL*8 CARG
        COMPLEX*16 Z
        CARG = ATAN2(AIMAG(Z),REAL(Z))
      END
C======================================================================
       FUNCTION POLAR(RHO,THETA)
       REAL*8 RHO,THETA,RP,IP
       COMPLEX*16 POLAR

       RP=RHO*COS(THETA)
       IP=RHO*SIN(THETA)
       POLAR=CMPLX(RP,IP)
       RETURN
       END
C=====================================================================

      subroutine croots(p1,q1,el1,el2,dg,c,zr)
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
      INTEGER dg,i,j,k,l,n,maxit
      REAL*8 eps,dmax,t0,t1,t2,p1(3),q1(3),el1(5),el2(5)
      parameter(maxit=100,dmax=1.d-6)
      COMPLEX*16 p(dg+1),q(dg+1),dp,f,df,r,z,zn,discr
      COMPLEX*16 c(DG+1),zr(DG),nz(4)

c======================================================================

      z=cmplx(0.d0,0.d0)
      call root_at_node(q1,p1,el1,el2,nz)
      k=1
      l=1
      do n=dg,3,-1
c Set the polynomial coefficients.
         do j=1,n+1
            p(j)=c(j)
         end do
c The starting values of the 5th,6th,7th and 8th roots
c are chosen so as they lie at the true anomaly of the mutual nodes
c of Orb1 and Orb2. The 9th root is restarted from z=(0.,0.)
         if (l.gt.4.and.l.lt.9) then
             z=nz(k)
             k=k+1
         else if (l.eq.9) then
             z=cmplx(0.d0,0.d0)
         end if
c         print*,n,l,z
c The root-finding iterations start here. The maximum number
c of iterations "MAXIT" is set in the header of the subroutine
  17     do i=1,maxit
            f=c(n+1)
            df=0.d0
c Evaluate simulatneously the polynomial 'f' and its
c derivative 'df'.
            do j=n,1,-1
               df=df*z+f
               f =f *z+c(j)
            end do
c Calculate new apprx. values for 'z', based on the coefficients
c 'C_k', using the Newton's method.
            zn=z-(f/df)
            eps=abs(zn)*dmax       ! Set 'eps' for root convergence
            if (abs(zn-z).le.eps) then
               zr(l)=zn
c              print*,abs(zr(l)),i,n
               exit
            else if(abs(zn-z).gt.eps.and.i.eq.maxit) then
c               z=1./zr(n)
c               go to 17
c               write(6,"('Solution Exceeded Max. Iterations')")
               exit
            end if
c The new 'z' value for continuing iterations
            z=zn
         end do
c Divide the polynomial P(z)=C_(k+1)*z^k, by the approximated root and obtain
c a new polynomial Q(z) of degree (k-1). Here the Horner's scheme and synthetic
c division have been used. By extracting a root in this way, we essentially have
c P(z)=(z-root)*Q(z). Q(z) is the quotinent of the division and if 'z' is indeed
c a root, the the remainder must be 0.000.
         call horner (p,n+1,zr(l),q)
         do j=1,n
            c(j)=q(j+1)     ! Adjust the new coefficients
         end do
c Set the new starting point to z=1./z^*. This has been suggested by
c (Baluev and Kholshevnikov, 2018). This does not apply to roots 5,6,7 and 8.
c Their apprx. values are started at the mutual nodes of the two orbits (see above).
         z=1.d0/conjg(zr(l))
         l=l+1
      end do
c The last two roots are found via the quadratic formula. Note the reverse order
c of the roots, starting from 1 and ending at 16, rather than other way around.
      discr=c(2)*c(2) - 4.*c(1)*c(3)
      zr(16)=(-c(2)+sqrt(discr))/2./c(3)
      zr(15)=(-c(2)-sqrt(discr))/2./c(3)
      END
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c     HORNER.FOR
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c This subroutine uses the Horner's method of deflation the degree
c of a polynomial via SYNTHETIC division.
      subroutine horner(p,d,root,q)
      integer i,j,k,d
      complex*16 p(d),root,q(d)

      do i=1,d
        q(i)=p(i)
      end do

      if (abs(root).lt.1.d0) then
         do k=d-1,1,-1
            q(k)=q(k+1)*root + p(k)
         end do
      else if (abs(root).ge.1.d0) then
         do k=1,d-1,1
           q(k+1)=q(k)/root + p(k+1)
         end do
         q(d:1:-1)=q(1:d:1)
      endif
      end
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FUNCTION truan(EAN,ECC)
      REAL*8 EAN,truan,ECC,sqec1,sqec2,te2

      sqec1=sqrt(1.d0+ecc)
      sqec2=sqrt(1.d0-ecc)
      te2  =tan(ean/2.d0)
      truan=2.d0*atan((sqec1/sqec2)*te2)
      RETURN
      END

      subroutine root_at_node(q1,p1,elm1,elm2,root)

      implicit none
      integer i
      real*8 dot
      real*8 elm1(5),elm2(5),q1(3),p1(3),w(3),sc(2)
      real*8 ci1,ci2,si1,si2,co1,co2,so1,so2,p1w,q1w
      real*8 tnode,eta
      complex*16 polar,root(4)
      external dot,polar

      ci1=cos(elm1(3))
      ci2=cos(elm2(3))
      si1=sin(elm1(3))
      si2=sin(elm2(3))
      co1=cos(elm1(4))
      co2=cos(elm2(4))
      so1=sin(elm1(4))
      so2=sin(elm2(4))
      eta=1.-elm1(2)*elm1(2)

      w(1)=ci1*si2*co2-si1*ci2*co1
      w(2)=ci1*si2*so2-si1*ci2*so1
      w(3)=si1*si2*sin(elm2(4)-elm1(4))

      q1w=dot(q1,w)
      p1w=dot(p1,w)

      tnode=atan2(q1w,p1w)
      sc(1)=sin(tnode)
      sc(2)=cos(tnode)
      root(1)=polar(1.d0, atan2( sc(1)*eta,elm1(2)+sc(2)))
      root(2)=polar(1.d0, atan2(-sc(1)*eta,elm1(2)-sc(2)))
      root(3)=polar(1.d0, atan2(-sc(2)*eta,elm1(2)+sc(1)))
      root(4)=polar(1.d0, atan2( sc(2)*eta,elm1(2)-sc(1)))

      end
