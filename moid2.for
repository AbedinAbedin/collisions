        implicit none

        integer i,j,k,l,indmin,p,i1,i2,invmod,ddm,ddm2
        real*8, allocatable, dimension(:):: truan1,truan2
        real*8 semi1,ecc1,incl1,aper1,lasc1,semi2,ecc2,incl2,aper2,lasc2
        real*8, allocatable,dimension(:):: moid,x1,y1,z1,x2,y2,z2
        real*8 meanan1, meanan2,x,y,z,r1,r2
        real*8, allocatable, dimension(:) :: dist
        real*8 pi,d2r,r2d,truan,m1,m2,m1min,m1max,m2min,m2max,dm1,dm2
        external invmod

        ddm = 100
        ddm2=ddm*ddm

        pi = 4.0*atan(1.d0)
        r2d = 180./pi
        d2r = 1./r2d

        semi1 = 1.44270708d0
        ecc1 = 0.483925457d0
        incl1 = 4.098859d0*d2r
        lasc1 = 213.749769d0*d2r
        aper1 = 102.916097d0*d2r
        meanan1 = 214.3092098d0*d2r

        semi2 = 5.2034d0
        ecc2 = 0.0483d0
        incl2 = 1.305d0*d2r
        lasc2 = 100.556d0*d2r
        aper2 = 274d0*d2r
        meanan2 = 2.55767d0*d2r

        allocate(truan1(ddm))
        allocate(truan2(ddm))
        allocate(x1(ddm))
        allocate(y1(ddm))
        allocate(z1(ddm))
        allocate(x2(ddm))
        allocate(y2(ddm))
        allocate(z2(ddm))
        allocate(moid(4))
        allocate(dist(ddm2))

        m1min=meanan1
        m2min=meanan2
        dm1 = 2*pi/(ddm-1)
        dm2 = dm1

        do l = 1,4
          k = 0
          do i = 1, ddm
             m1 = m1min + dm1*(i-1)
            call kep2cart(semi1,ecc1,incl1,aper1,lasc1,m1,
     %          x1(i),y1(i),z1(i),truan1(i))

            do j = 1, ddm
              k = k + 1
              m2 = m2min + dm2*(j-1)
              call kep2cart(semi2,ecc2,incl2,aper2,lasc2,m2,
     %           x2(j),y2(j),z2(j),truan2(j))

              dist(k) = sqrt((x1(i)-x2(j))**2 + (y1(i)-y2(j))**2 +
     %                  (z1(i)-z2(j))**2)
c              print*,k,r1,r2
c     %           truan2(j)*r2d,dist(k)
            end do
          end do
          call FindInVector(ddm2,dist.eq.minval(dist),p)
          i1 = invmod(ddm,p,ddm)
          i2 = mod(p,ddm)
          moid(l)=minval(dist)
          print*,truan1(i1),truan2(i2),moid(l)
          call t2m(truan1(i1-1),ecc1,m1min)
          call t2m(truan1(i1+1),ecc1,m1max)
          call t2m(truan2(i2-1),ecc2,m2min)
          call t2m(truan2(i2+1),ecc2,m2max)
c          print*,m1min,m1max
          if (m1min.gt.m1max) then
            dm1 = (-m1max + m1min)/ddm
          else
            dm1 = (m1max-m1min)/ddm
          end if
          if (m2min.gt.m2max) then
            dm2 = (-m2max + m2min)/ddm
          else
            dm2 = (m2max-m2min)/ddm
          end if

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
