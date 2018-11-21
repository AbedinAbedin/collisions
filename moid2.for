        implicit none

        integer i,j,k,l,indmin,p,i1,i2,invmod
        real*8, allocatable, dimension(:):: truan1,truan2
        real*8 semi1,ecc1,incl1,aper1,lasc1,semi2,ecc2,incl2,aper2,
     %   lasc2
        real*8, allocatable,dimension(:):: moid,x1,y1,z1,x2,y2,z2
        real*8 meanan1, meanan2,x,y,z
        real*8, allocatable, dimension(:) :: dist
        real*8 pi,d2r,r2d,truan,m1,m2,m1min,m1max,m2min,m2max,dm1,dm2
        external invmod

        pi = 4.0*atan(1.d0)
        r2d = 180./pi
        d2r = 1./r2d

        semi1 = 1.44270708d0
        ecc1 = 0.483925457d0
        incl1 = 4.098859d0*d2r
        lasc1 = 213.749769d0*d2r
        aper1 = 102.916097d0*d2r
        meanan1 = 214.3092098d0*d2r

        semi2 = 0.99943376943d0
        ecc2 = 0.0161742442d0
        incl2 = 0.002812662230d0*d2r
        lasc2 = 138.092792634d0*d2r
        aper2 = 326.1612866d0*d2r
        meanan2 = 355.74d0*d2r

        allocate(truan1(10))
        allocate(truan2(10))
        allocate(x1(10))
        allocate(y1(10))
        allocate(z1(10))
        allocate(x2(10))
        allocate(y2(10))
        allocate(z2(10))
        allocate(dist(100))

        m1min = meanan1
        m1max = meanan1 + 2.*pi
        dm1 = abs(m1max - m1min)/10.

        m2min = meanan2
        m2max = meanan2 + 2.*pi
        dm2 = abs(m2max - m2min)/10.

        do l = 1,7
          if(m1min.gt.2*pi) m1min=mod(m1min,2*pi)
          if(m1max.gt.2*pi) m1max=mod(m1max,2*pi)
          if(m2min.gt.2*pi) m2min=mod(m2min,2*pi)
          if(m2max.gt.2*pi) m2max=mod(m2max,2*pi)
          if(m1min.lt.0.0) m1min=m1min+2*pi
          if(m1max.lt.0.0) m1max=m1max+2*pi
          if(m2min.lt.0.0) m2min=m2min+2*pi
          if(m2max.lt.0.0) m2max=m2max+2*pi
          k = 0
          do i = 1, 10
             m1 = m1min + dm1*(i-1)
            call kep2cart(semi1,ecc1,incl1,aper1,lasc1,m1,
     %          x1(i),y1(i),z1(i),truan1(i))
            if (truan1(i).lt.0.0) truan1(i) = truan1(i) + 2*pi
            do j = 1, 10
              k = k + 1
              m2 = m2min + (j-1)*dm2
              call kep2cart(semi2,ecc2,incl2,aper2,lasc2,m2,
     %           x2(j),y2(j),z2(j),truan2(j))
              if (truan2(j).lt.0.0) truan2(j) = truan2(j) + 2*pi
              dist(k) = sqrt((x1(i)-x2(j))**2 + (y1(i)-y2(j))**2 +
     %                  (z1(i)-z2(j))**2)
c              print*,k,i,j,truan1(i)*r2d,truan2(j)*r2d,dist(k)
            end do
          end do
          call FindInVector(100,dist.eq.minval(dist),p)
          i1 = invmod(p,10)
          i2 = mod(p,10)
          
c          print*,m1min*r2d,m1max*r2d,dm1*r2d,m2min*r2d,m2max*r2d,dm2*r2d
          call t2m(truan1(i1-1),ecc1,m1min)
          call t2m(truan1(i1+1),ecc1,m1max)
          call t2m(truan2(i2-1),ecc2,m2min)
          call t2m(truan2(i2+1),ecc2,m2max)
          if (m1min.gt.m1max) then
            dm1=(2*pi-m1min+m1max)/10.
          else
            dm1=abs(m1max-m1min)/10.
          end if
          if (m2min.gt.m2max) then
            dm2=(2*pi-m2min+m2max)/10.
          else
            dm2=abs(m2max-m2min)/10.
          end if
c          print*,i1,i2,truan1(i1)*r2d,truan2(i2)*r2d,minval(dist)
          print*,truan1(i1-1)*r2d,truan1(i1+1)*r2d,truan2(i2-1)*r2d,
     %           truan2(i2+1)*r2d,minval(dist),dm1,dm2
        end do

        end

        SUBROUTINE KEP2CART(semi,ecc,incl,aper,lasc,man,x,y,z,tanom)

        implicit none
        real*8 pi,rad,k
        real*8 x,y,z,rhc,p
        real*8 semi,ecc,incl,lasc,aper,tanom,man

        pi = 4.d0 * atan(1.d0)
        rad = pi / 180.

        call mean2tan(man,ecc,tanom)
        if (tanom .lt. 0d0) tanom = tanom + 2.*pi
        p  = semi * (1. - ecc**2)
        rhc = p/(1.+ecc * cos(tanom))
        x = rhc * (cos(lasc) * cos(aper + tanom)  -  sin(lasc) *
     %            sin(aper + tanom) * cos(incl))
        y = rhc * (sin(lasc) * cos(aper + tanom)  +  cos(lasc) *
     %           sin(aper + tanom) * cos(incl))
        z = rhc * sin(aper + tanom) * sin(incl)

        END

        SUBROUTINE MEAN2TAN(manom,ecc,tanom)

        implicit none
        integer k
        REAL*8 manom,tanom,eccan,eccan0,ecc

        eccan0 = manom

        do k = 1, 5
          eccan = manom + ecc*sin(eccan0)
          if (abs(eccan-eccan0).le.1.e-7) then
              exit
          else if(abs(eccan-eccan0).gt.1.e-7 .and. k.eq.500) then
              print*,"KEPLER'S EQUATION FAILD TO CONVERGE:"
          end if
          eccan0=eccan
        end do
        
        eccan = eccan0
        tanom = 2.d0*atan(tan(eccan/2.0) * sqrt((1+ecc)/(1-ecc)))

        end

        subroutine t2m(tanom,ecc,manom)

        real*8 tanom,ecc,manom,ecanom

        ecanom = 2.d0*atan(tan(tanom/2.0) * sqrt((1-ecc)/(1+ecc)))
        manom = ecanom - ecc*sin(ecanom)

        end



        function invmod(m,n)
        integer m,n

        invmod = ((m - mod(m,n))/10) + 1
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
