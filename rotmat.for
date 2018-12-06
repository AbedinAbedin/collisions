       IMPLICIT NONE

       real*8 x,y,z,xn,yn,zn,incl,lasc,aper
       real*8 semi,ecc,rhel,ppp,qqq,pi,si,sw,so,ci,cw,co

       pi=4.d0*atan(1.d0)

       semi=9.223077829454267E-01
       ecc=1.911162084383658E-01
       lasc=2.044313460672734E+02*pi/180.
       incl=3.331864522652364E+00*pi/180.
       aper=1.264302212363639E+02*pi/180.

       ppp=semi*(1.-ecc*ecc)
       qqq=semi*(1.-ecc)
       
       rhel = ppp/(1+ecc)
       si=sin(incl)
       sw=sin(aper)
       so=sin(lasc)
       ci=cos(incl)
       cw=cos(aper)
       co=cos(lasc)
       

       x=rhel*(cw*co-ci*so*sw)
       y=rhel*(cw*so+ci*co*sw)
       z=rhel*sw*si
       
c       x = 7.640651955680458E-01
c       y = 5.870095330756822E-02
c       z = 1.528647623448541E-02

c 3D transformation of heliocentric ecliptic coordinates
c of a body to the ecliptic ref. plane in performed via 3 
c succesive rotations of the orbit about the axes Z, X, Z
c to angles -LASC 'Om', -INCL 'I' and -APER 'W' respectively.
c NOTE the negative signs, corresponding to clockwise rotation.
       call rotmat('Z',-lasc,x,y,z,xn,yn,zn)
       call rotmat('X',-incl,xn,yn,zn,x,y,z)
       call rotmat('Z',-aper,x,y,z,xn,yn,zn)
       print*,xn,yn,zn
       END


       SUBROUTINE ROTMAT(AX,PHI,XOLD,YOLD,ZOLD,XNEW,YNEW,ZNEW)
c
c Subroutine for performing transformation (rotation) around a
c given Axis.
c              
c INPUT:
c--------------------------------------------------------------
c AX - ('X', 'Y' or 'Z') the axis around the which the rotation
c      is performed. Lower case characters (x,y,z) are also 
c      allowed.
c PHI - the angle of rotatation around given axis. This angle
c       can be either (+) for counter-clockwise or (-) for
c       clockwise rotations, respectively.
c
c XOLD, YOLD, ZOLD - X,Y,Z (old) coordinates which will be 
c                    transformed
c--------------------------------------------------------------
c OUTPUT:
c
c XNEW, YNEW, ZNEW - The new X',Y',Z' after the transformation
c
       IMPLICIT NONE

       REAL*8 XOLD,YOLD,ZOLD,XNEW,YNEW,ZNEW,PHI,CP,SP       
       CHARACTER*1 AX


       CP=COS(PHI)
       SP=SIN(PHI)

       IF (AX.EQ.'X' .OR. AX.EQ.'x') THEN
         XNEW=XOLD
         YNEW=YOLD*CP-ZOLD*SP
         ZNEW=YOLD*SP+ZOLD*CP
       ELSE IF (AX.EQ.'Y' .OR. AX.EQ.'y') THEN
         XNEW=XOLD*CP+ZOLD*SP
         YNEW=YOLD
         ZNEW=-XOLD*SP+ZOLD*CP
       ELSE IF (AX.EQ.'Z' .OR. AX.EQ.'z') THEN
         XNEW=XOLD*CP-YOLD*SP
         YNEW=XOLD*SP+YOLD*CP
         ZNEW=ZOLD
       ELSE IF (AX.NE.'X'.OR.AX.NE.'x'.OR.AX.NE.'Y'.OR.AX.NE.'y'
     %          .or.AX.NE.'Z'.OR.AX.NE.'z') THEN
          WRITE(6,*)'ERROR IN TRANSMAT: Your choice for //
     %               Axis rotation does not match X,Y,Z or x,y,z'
          STOP
       END IF

       END
