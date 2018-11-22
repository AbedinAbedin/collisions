       IMPLICIT NONE

       INTEGER i,j,k,l,m,n,sig
       REAL*8 semi1,ecc1,incl1,aper1,lasc1
       REAL*8 semi2,ecc2,incl2,aper2,lasc2,moid(16)
       REAL*8 p1(3),p2(3),q1(3),q2(3),s1(3),s2(3)
       REAL aa,bb,cc,mm,nn,kk,alpha1,alpha2,u1(16),u2(16),gg,u(0:20)
       REAL*8 pi,dot,cons,derv,func,eps,x,y,dist,carg,dd
       REAL*8 csu,ssu,mnsc,kksc
       COMPLEX m1,m2,m3,m4,c8,eps2,roots(16),nroots(16)
       COMPLEX tt,ii,c_exp,c(-8:8),coeff(17),sumgg
       EXTERNAL dot,c_exp,derv,func,carg
       LOGICAL changed 
	
       pi = 4.0*atan(1.d0)
       ii = (0,1)

       semi1 = 44.63
       ecc1  = 0.05
       incl1 = 2.45
       aper1 = 181.148
       lasc1 = 159.044

       semi2 = 5.2401
       ecc2  = 0.110
       incl2 = 0.9
       aper2 = 10.0
       lasc2 = 10.0

        print*,'pi=',pi

       incl1 = incl1 * pi / 180.
       aper1 = aper1 * pi / 180.
       lasc1 = lasc1 * pi / 180.

       incl2 = incl2 * pi / 180.
       aper2 = aper2 * pi / 180.
       lasc2 = lasc2 * pi / 180.

       p1(1) = cos(aper1)*cos(lasc1) - cos(incl1)*sin(aper1)*sin(lasc1)
       p1(2) = cos(aper1)*sin(lasc1) + cos(incl1)*sin(aper1)*cos(lasc1)
       p1(3) = sin(incl1)*sin(aper1)

       q1(1) = -sin(aper1)*cos(lasc1) - cos(incl1)*cos(aper1)*sin(lasc1)
       q1(2) = -sin(aper1)*sin(lasc1) + cos(incl1)*cos(aper1)*cos(lasc1)
       q1(3) =  sin(incl1)*cos(aper1)

       s1(1) = sqrt(1 - ecc1**2)*q1(1)
       s1(2) = sqrt(1 - ecc1**2)*q1(2)
       s1(3) = sqrt(1 - ecc1**2)*q1(3)
       alpha1 = semi1/semi2

       p2(1) = cos(aper2)*cos(lasc2) - cos(incl2)*sin(aper2)*sin(lasc2)
       p2(2) = cos(aper2)*sin(lasc2) + cos(incl2)*sin(aper2)*cos(lasc2)
       p2(3) = sin(incl2)*sin(aper2)

       q2(1) = -sin(aper2)*cos(lasc2) - cos(incl2)*cos(aper2)*sin(lasc2)
       q2(2) = -sin(aper2)*sin(lasc2) + cos(incl2)*cos(aper2)*cos(lasc2)
       q2(3) =  sin(incl2)*cos(aper2)

       s2(1) = sqrt(1 - ecc2**2)*q2(1)
       s2(2) = sqrt(1 - ecc2**2)*q2(2)
       s2(3) = sqrt(1 - ecc2**2)*q2(3)
       alpha2 = 1./alpha1

! Explicit formulas for C_8:
      m1=dot(p1,p2)-dot(s1,s2)-ecc1*ecc2-ii*(dot(s1,p2)+dot(p1,s2))
      m2=dot(p1,p2)-dot(s1,s2)+ecc1*ecc2-ii*(dot(s1,p2)+dot(p1,s2))
      m3=dot(p1,p2)+dot(s1,s2)-ecc1*ecc2-ii*(dot(s1,p2)-dot(p1,s2))
      m4=dot(p1,p2)+dot(s1,s2)+ecc1*ecc2-ii*(dot(s1,p2)-dot(p1,s2))
      c8=m1*m2*m3*m4*(alpha1*(ecc1**2)/16.)**2

! Calculate the coefficients C_k.
      print*,'********* Printing calculated Coefficientc C_k:'
       do k = -8, 8
         sumgg = 0.0
         cons = 1./21.
         do j = 0, 20
            u(j) = 2.*pi*j*cons
            aa = (dot(p1,s2) * sin(u(j))) - (dot(s1,s2) * cos(u(j)))
            bb = dot(p1,p2) * sin(u(j)) - dot(s1,p2) * cos(u(j))
            cc = (ecc2 * bb)- alpha1*ecc1*sin(u(j))
     %           * (1 - ecc1*cos(u(j)))
            mm = dot(p1,p2) * cos(u(j)) + dot(s1,p2)*sin(u(j))
     %           + alpha1*ecc2 - dot(p1,p2)*ecc1
            nn = dot(p1,s2) * ecc1 - dot(s1,s2)*sin(u(j))
     %           - dot(p1,s2)*cos(u(j))
            kk = alpha2 * (ecc2**2)
            gg = ((aa**2 - cc**2)*(bb**2 - cc**2)*kk**2)
     %           + 2.*kk*cc*(nn*aa*(aa**2-cc**2) + mm*bb*(bb**2-cc**2))
     %           - (aa**2 + bb**2)*((aa**2-cc**2)*nn**2+((bb**2-cc**2)
     %           * mm**2) - (2.*nn*mm*aa*bb))
            print*,'aa=',aa
            print*,'bb=',bb
            print*,'cc=',cc
            print*,'mm=',mm
            print*,'nn=',nn
            print*,'kk=',kk
            print*,'gg=',gg
            sumgg = sumgg + gg*cexp(ii*k*u(j))
         end do
         c(k) = sumgg*cons
         write(6,"(a2,i2,a2,'(',ES16.9E2,',',ES16.9E2,')')")'C_',k,'=',
     %         c(k)
       end do 

        k = 0
        do i = -8,8
          k = k + 1
          coeff(k) = c(i)
        end do 
        write(6,"(/,15x,a13,15x,a21)")'EXPLICIT C_8:','DFT CALC. C_8:'
        print*,k,c8,coeff(17)
        print*,''

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

        write(6,"(/,a60)")'*'
        write(6, "('         PERFORMING SIGN CHECK FOR (m)    ',/) ")
	do i = 1, j
          u1(i) = carg(nroots(i))    ! calculate the argument(phase) of the roots
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
c          write(6,"(/,'MNSC =',1p,E18.10,10x,'KKSC = ',
c     %          1p,E18.10 )")mnsc,kksc
          if ((abs(mnsc-kksc).gt.1E-5) .and. (.not.changed)) then
            sig = -1
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
          print*,moid(i)
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
	COMPLEX Z

	CARG = ATAN2(AIMAG(Z),REAL(Z))
	RETURN
	END
C======================================================================
       FUNCTION c_exp(angle)
	 COMPLEX c_exp
       real*8 angle   ! The angle must be in radians.

       c_exp = cmplx(cos(angle),sin(angle))
       return
	 END
C======================================================================
       SUBROUTINE laguer(a,m,x,its)
       INTEGER m,its,MAXIT,MR,MT
       REAL*8 EPSS
       COMPLEX a(m+1),x
       PARAMETER (EPSS=6.e-8,MR=8,MT=100,MAXIT=MT*MR)
       INTEGER iter,j
       REAL*8 abx,abp,abm,err,frac(MR)
       COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
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
        COMPLEX a(m+1),roots(m)  
        LOGICAL polish  
        PARAMETER (EPS=6.e-8,MAXM=101)  
CU    USES laguer
      	
        INTEGER i,j,jj,its  
        COMPLEX ad(MAXM),x,b,c  
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
