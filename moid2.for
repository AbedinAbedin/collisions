        implicit none

        integer i,j,k,indmin,p,i1,i2,invmod
        real*8, allocatable, dimension(:):: truan1,truan2
        real*8 semi1,ecc1,incl1,aper1,lasc1,semi2,ecc2,incl2,aper2,lasc2
        real*8, allocatable,dimension(:):: moid,x1,y1,z1,x2,y2,z2
        real*8 meanan1, meanan2,x,y,z
        real*8, allocatable, dimension(:) :: dist
        real*8 pi,d2r,r2d,truan,m1,m2
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

c        call kep2cart(semi1,ecc1,incl1,aper1,lasc1,meanan1,
c     %       x,y,z,truan1(1))
c        x1(1) = x
c        y1(1) = y
c        z1(1) = z
c        write(6,"(/,'X = ',1p,E18.10,/'Y = ',1p,E18.10,/,'Z = ',1p,
c     %        E18.10)")x,y,z
c        call kep2cart(semi2,ecc2,incl2,aper2,lasc2,meanan2,
c     %       x,y,z,truan2(1))
c        x2(1) = x
c        y2(1) = y
c        z2(1) = z
c        write(6,"(/,'X = ',1p,E18.10,/'Y = ',1p,E18.10,/,'Z = ',1p,
c     %        E18.10)")x,y,z

        allocate(truan1(10))
        allocate(truan2(10))
        allocate(x1(10))
        allocate(y1(10))
        allocate(z1(10))
        allocate(x2(10))
        allocate(y2(10))
        allocate(z2(10))
        allocate(dist(100))

        k = 0
        do i = 1, 10
          m1 = meanan1 + (i-1)*2.*pi/10.
          call kep2cart(semi1,ecc1,incl1,aper1,lasc1,m1,
     %          x1(i),y1(i),z1(i),truan1(i))
          do j = 1, 10
            k = k + 1
            m2 = meanan2 + (j-1)*2.*pi/10.
            call kep2cart(semi2,ecc2,incl2,aper2,lasc2,m2,
     %         x2(j),y2(j),z2(j),truan2(j))
            dist(k) = sqrt((x1(i)-x2(j))**2 + (y1(i)-y2(j))**2 +
     %                  (z1(i)-z2(j))**2)
c            print*,k,i,j,truan1(i)*r2d,truan2(j)*r2d,dist(k)
          end do
        end do
        call FindInVector(100,dist.eq.minval(dist),p)
c        print*,p,minval(dist),truan1(p)*r2d,truan2(p)*r2d
        i1 = invmod(p,10)
        i2 = mod(p,10)
        print*,i1,i2
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
