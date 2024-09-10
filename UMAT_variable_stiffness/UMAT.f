      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)

      INTEGER :: I, J, MEMBRANE
      REAL(8) :: TOT, D, ANG
      REAL(8) :: E11, E22, NU12, G12
      REAL(8) :: Q11, Q22, Q12, Q66
      REAL(8) :: QBAR(3,3)
      REAL(8) :: AMAT(3,3), DMAT(3,3)
      REAL(8) :: A, T0, T1, TH, PHI
      REAL(8) :: PI
      
C     NUMBER PI
      PI = 3.14159265358979323846D+0

C      WRITE(*,*) 'COORDS:'
C      WRITE(*,*) COORDS
C      WRITE(*,*) 'LAYER: ', LAYER
C      WRITE(*,*) 'NTENS: ', NTENS
C      WRITE(*,*) 'NDI: ', NDI
C      WRITE(*,*) 'NSHR: ', NSHR

      
C     READ INPUTS
C     ELASTIC CONSTANTS
      E11=PROPS(1)
      E22=PROPS(2)
      NU12=PROPS(3)
      G12=PROPS(4)
C     GEOMETRY RELATED
      T0=PROPS(5)
      T1=PROPS(6)
      PHI=PROPS(7)
      A=PROPS(8)
      TH=PROPS(9)
C     FLAG FOR THE MEMBRANE ANALYSIS
C     IN-PLANE ANALYSIS: 1
C     OUT-OF-PLANE ANALSIS: 0
      MEMBRANE=INT(PROPS(10))
      
      
C     CALCULATE POSITION
C     ASSUMING PLATE IS CENTERED AT THE ORIGIN
      D = SQRT(COORDS(1)**2. + COORDS(2)**2.)
      
CC     ANGLE PHI
C      PHI = ATAN2(COORDS(2),COORDS(1))*180./PI
      
C     CALCULATE THE ANGLE (DEG)
      ANG = PHI+ 2.*(T1-T0)/A*D + T0 

CC     MANUALLY SET THE ANGLE
C      ANG=90.
      
C     OUTPUT THE ANGLES TO STATE VARIABLES
      STATEV(1)=ANG
      
C     CONVERT TO RADIANS
      ANG = ANG * PI / 180.
      

      
C     LAMINA STIFFNESS      
      
      
C     CALCULATION OF STIFFNESS PARAMETERS
      Q11 = E11**2./(E11-(NU12**2.*E22))
      Q12 = NU12*E11*E22/(E11-(NU12**2.*E22))
      Q22 = E11*E22/(E11-(NU12**2.*E22))
      Q66 = G12
      
      QBAR=0.
C     TRANSFORMED STIFFNESS
      QBAR(1,1) = COS(ANG)**4.*Q11 + SIN(ANG)**4.*Q22 +
     + COS(ANG)**2.*SIN(ANG)**2.*(2.*Q12+4.*Q66)
      QBAR(1,2) = COS(ANG)**2.*SIN(ANG)**2.*(Q11+Q22-4.*Q66) +
     + (COS(ANG)**4. + SIN(ANG)**4.)*Q12
      QBAR(1,3) = COS(ANG)**3.*SIN(ANG)*(Q11-Q12-2.*Q66) +
     + COS(ANG)*SIN(ANG)**3.*(Q12-Q22+2.*Q66)
      QBAR(2,2) = SIN(ANG)**4.*Q11+
     + COS(ANG)**2.*SIN(ANG)**2.*(2.*Q12+4.*Q66) + COS(ANG)**4.*Q22
      QBAR(2,3) = SIN(ANG)**3.*COS(ANG)*(Q11-Q12-2.*Q66) +
     + COS(ANG)**3.*SIN(ANG)*(Q12-Q22+2.*Q66)
      QBAR(3,3) = COS(ANG)**2.*SIN(ANG)**2.*(Q11+Q22-2.*Q12-2.*Q66) +
     + (SIN(ANG)**4. + COS(ANG)**4.)*Q66

      QBAR(2,1) = QBAR(1,2)
      QBAR(3,1) = QBAR(1,3)
      QBAR(3,2) = QBAR(2,3)
      
C     SHEAR CORRECTION
      QBAR(3,3) = 2.*QBAR(3,3)
      
      
C     BENDING STIFFNESS
C     FOR A SINGLE LAYER ONLY
      DMAT=1./3. * QBAR * TH**3.
      
C     IN-PLANE STIFFNESS
C     FOR A SINGLE LAYER ONLY
      AMAT = QBAR * TH
      

C     AMAT FOR IN-PLANE ANALYSIS
      IF (MEMBRANE.EQ.1) THEN
          
          DO I=1,NTENS
              DO J=1,NTENS
                  DDSDDE(I,J)=AMAT(I,J)
              END DO
          END DO
          
C     DMAT FOR OUT-OF-PLANE ANALYSIS
      ELSEIF (MEMBRANE.EQ.0) THEN
          
          DO I=1,NTENS
              DO J=1,NTENS
                  DDSDDE(I,J)=DMAT(I,J)
              END DO
          END DO
      
      ENDIF
      
            
      DO I=1,NTENS
          TOT =0.0
          DO J=1,NTENS
              TOT = TOT + DDSDDE(I,J) * DSTRAN(J)     
C              WRITE(6,*) ELAS(I,J)     
          ENDDO
          STRESS(I) = STRESS(I) +  TOT
      ENDDO
      

      RETURN
      END