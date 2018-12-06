        implicit none

        real*8 pi,r2d,d2r
        real*8 semi1,ecc1,incl1,aper1,lasc1
        real*8 semi2,ecc2,incl2,aper2,lasc2
        real*8 meanan1,meanan2,truan1,truan2,moid

        pi = 4.0*atan(1.d0)
        r2d = 180./pi
        d2r = 1./r2d

        semi1 = 1.000525749372538E+00
        ecc1  = 1.731749097026299E-02 
        incl1 = 0.0 *d2r
        aper1 = 2.354347067156837E+02 *d2r
        lasc1 = 0.0 *d2r

        semi2 = 9.223077829454267E-01
        ecc2  = 1.911162084383658E-01
        incl2 = 3.331864522652364E+00 *d2r
        aper2 = 1.264302212363639E+02 *d2r
        lasc2 = 2.044313460672734E+02 *d2r
       
        
        call FIND_MOID(semi1,ecc1,incl1,lasc1,aper1,
     %                 semi2,ecc2,incl2,lasc2,aper2)
        end


        SUBROUTINE FIND_MOID(semi1,ecc1,incl1,lasc1,aper1,
     %                       semi2,ecc2,incl2,lasc2,aper2) 
c Direct search of MOID between two orbits O1 = <a,e,i,w,Om> and 
c               O2 = <a',e',i',w',Om'>. 
c NOTE!!! In order to avoid solving the Kepler's equation for 
c         TruAn or MeanAn, we work in EccAn.
        implicit none

        integer i,j,k,l,m,n,indmin,p,i1,i2,invmod,b,b2,inu1,inu2
        parameter(b=50,b2=b*b)
        real*8 semi1,ecc1,incl1,aper1,lasc1,semi2,ecc2,incl2,aper2,lasc2
        real*8 x1,y1,z1,x2,y2,z2,r1,r2,minu1,maxu1,minu2,maxu2,pi
        real*8 p1(3),q1(3),s1(3),p2(3),q2(3),s2(3),cu1,su1,cu2,su2
        real*8 dist(b2),u1(b2),u2(b),c1,c2
        external invmod
       
        open(20,file='mindist.dat')
        pi = 4.d0*atan(1.d0)

        p1(1) = cos(aper1)*cos(lasc1)-cos(incl1)*sin(aper1)*sin(lasc1)
        p1(2) = cos(aper1)*sin(lasc1)+cos(incl1)*sin(aper1)*cos(lasc1)
        p1(3) = sin(incl1)*sin(aper1)

        q1(1) = -sin(aper1)*cos(lasc1)-cos(incl1)*cos(aper1)*sin(lasc1)
        q1(2) = -sin(aper1)*sin(lasc1)+cos(incl1)*cos(aper1)*cos(lasc1)
        q1(3) =  sin(incl1)*cos(aper1)

        s1(1) = sqrt(1.d0 - ecc1**2)*q1(1)
        s1(2) = sqrt(1.d0 - ecc1**2)*q1(2)
        s1(3) = sqrt(1.d0 - ecc1**2)*q1(3)

        p2(1) = cos(aper2)*cos(lasc2)-cos(incl2)*sin(aper2)*sin(lasc2)
        p2(2) = cos(aper2)*sin(lasc2)+cos(incl2)*sin(aper2)*cos(lasc2)
        p2(3) = sin(incl2)*sin(aper2)

        q2(1) = -sin(aper2)*cos(lasc2)-cos(incl2)*cos(aper2)*sin(lasc2)
        q2(2) = -sin(aper2)*sin(lasc2)+cos(incl2)*cos(aper2)*cos(lasc2)
        q2(3) =  sin(incl2)*cos(aper2)

        s2(1) = sqrt(1.d0 - ecc2**2)*q2(1)
        s2(2) = sqrt(1.d0 - ecc2**2)*q2(2)
        s2(3) = sqrt(1.d0 - ecc2**2)*q2(3)
     
         minu1=0.0
         maxu1=2.*pi
         minu2=0.0
         maxu2=2.*pi

        do i=1,2
          l=0
          dist=0d0
c          print*,minu1,maxu1,minu2,maxu2
          do j=1,b
            u2(j) = minu2 + (j-1)*abs(maxu2-minu2)/b
            cu2=cos(u2(j))
            su2=sin(u2(j))
            x2 = semi2*p2(1)*(cu2-ecc2) + s2(1)*su2
            y2 = semi2*p2(2)*(cu2-ecc2) + s2(2)*su2
            z2 = semi2*p2(3)*(cu2-ecc2) + s2(3)*su2
            do k=1,b
              l=l+1
              u1(l) = minu1 + (k-1)*abs(maxu1-minu1)/b
              cu1=cos(u1(l))
              su1=sin(u1(l))
              x1 = semi1*p1(1)*(cu1-ecc1) + s1(1)*su1
              y1 = semi1*p1(2)*(cu1-ecc1) + s1(2)*su1
              z1 = semi1*p1(3)*(cu1-ecc1) + s1(3)*su1
              dist(l)=sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
c              print*,j,k,l,u2(j),u1(k),dist(l)
            end do
          end do

          call FindInVector(l,(dist(1:l).eq.minval(dist(1:l))),inu1)
          if (mod(inu1,b).eq.0.0) then
            i1=b
          else
            i1=mod(inu1,b)
          end if
          minu1=u1(i1-1)
          maxu1=u1(i1+1)

c          print*,'MOD(200,200)=',mod(200,b)
c          print*,'MOD(1,200)=',mod(1,b)
c          print*,'MOD(100,200)=',mod(100,b)
c          print*,'MOD(400,200)=',mod(200,b)
c          print*,'MOD(500,200)=',mod(500,b)
c          print*,'MOD(513,200)=',mod(513,b)
          i2=int(inu1/b)+1
          
c          print*,i1,i2
          minu2=u2(i2-2)
          maxu2=u2(i2+2)
          print*,i1,u1(i1),i2,u2(i2),dist(inu1)
        end do
        end

c        SUBROUTINE KEP2CART(semi,ecc,incl,aper,lasc,man,x,y,z,tanom)
c
c        implicit none
c        real*8 pi,rad,k
c        real*8 x,y,z,rhc,p
c        real*8 semi,ecc,incl,lasc,aper,tanom,man

c        pi = 4.d0 * atan(1.d0)
c        rad = pi / 180.

c        call mean2tan(man,ecc,tanom)
      
c        p  = semi * (1. - ecc**2)
c        rhc = p/(1.+ecc * cos(tanom))
c        x = rhc * (cos(lasc) * cos(aper + tanom)  -  sin(lasc) *
c     %            sin(aper + tanom) * cos(incl))
c        y = rhc * (sin(lasc) * cos(aper + tanom)  +  cos(lasc) *
c     %           sin(aper + tanom) * cos(incl))
c        z = rhc * sin(aper + tanom) * sin(incl)

c        END

c        SUBROUTINE MEAN2TAN(manom,ecc,tanom)

c        implicit none
c        integer k
c        REAL*8 manom,tanom,eccan,eccan0,ecc

c        eccan0 = manom

c        do k = 1, 300
c          eccan = manom + ecc*sin(eccan0)
c          if (abs(eccan-eccan0).le.1.e-8) then
c              exit
c          else if(abs(eccan-eccan0).gt.1.e-7 .and. k.eq.500) then
c              print*,"KEPLER'S EQUATION FAILD TO CONVERGE:"
c          end if
c          eccan0=eccan
c        end do
        
c        eccan = eccan0
c        tanom = 2.d0*atan(tan(eccan/2.0) * sqrt((1+ecc)/(1-ecc)))

c        end

c        subroutine t2m(tanom,ecc,manom)
c
c        real*8 tanom,ecc,manom,ecanom
c
c        ecanom = 2.d0*atan(tan(tanom/2.0) * sqrt((1-ecc)/(1+ecc)))
c        manom = ecanom - ecc*sin(ecanom)

c        end



        function invmod(p,m,n)
        integer m,n,p

        invmod = ((m - mod(m,n))/p) + 1
        end

        SUBROUTINE FindInVector(n,TF,indx)
        ! Inlet variables
        INTEGER,INTENT(IN):: n      ! Dimension of logical vector
        LOGICAL,INTENT(IN):: TF(n)  ! Logical vector (True or False)
        ! Outlet variables
        INTEGER npos                ! number of "true" conditions
        INTEGER pos(n)              ! position of "true" conditions
    ! Internal variables
        INTEGER i                   ! counter
        INTEGER v(n)                ! vector of all positions
        integer indx

        pos = 0                     ! Initialize pos
        FORALL(i=1:n)   v(i) = i    ! Enumerate all positions
        npos  = COUNT(TF)           ! Count the elements of TF that are .True.
        pos(1:npos)= pack(v, TF)    ! With Pack function, verify position of true conditions
        
        do j = 1, npos
          if (pos(j).ne.0) indx = pos(j)
        end do
        END
