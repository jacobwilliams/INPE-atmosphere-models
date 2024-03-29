      FUNCTION TEMLO(ALTU,C)
C
C------
C
C PURPOSE:
C
C     THE FUNCTION TEMLO  EVALUATES  THE  TEMPERATURE AT THE
C     LOCAL REGION.
C
C INPUTS:
C       ALTU     ALTITUDE OF THE POINT, IN KM (90 TO 2000).
C       C        ARRAY THAT CONTAINS THE PARAMETERS USED  IN
C                THE EVALUATION OF EQUATIONS 3 AND 4, AS GI-
C                VEN BY JACCHIA.
C
C OUTPUTS:
C       TEMLO    LOCAL TEMPERATURE AS  DEFINED  IN EQUATIONS
C                3 AND 4, IN KELVIN.
C
C REFERENCES:
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA  APR/87    V 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION C(7)
C
C------
C  
      DATA ZX   /125.D0/
      DATA ZO   /90.0D0/
      DATA TO   /188.D0/
  
      HIGX   = ALTU - ZX
      HIGO   = ALTU - ZO
      TEMLO  = TO

      IF(HIGO.EQ.0.D0) RETURN
  
      IF(HIGX.GT.0.D0) THEN
       TEMLO = C(7) +
     1         C(4)*DATAN(C(5)*HIGX + C(6)*HIGX*HIGX*HIGX)
      ELSE
       AUXI  = HIGX/HIGO
       TEMLO = C(7) +
     1         C(1)*DATAN(C(2)*HIGX + C(3)*HIGX*AUXI*AUXI)
      ENDIF

      RETURN
      END
  
  
      SUBROUTINE SEALAT(TYFR,SUDC,RLAT,ALTU,AL)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE SEALAT  OBTAINS  THE VARIATIONS  ON THE
C
C     NUMBER DENSITY OF THE ATMOSPHERE, DUE TO THE SEAZONAL-
C                                                  ---
C     LATITUDINAL EFFECT.
C     ---
C
C INPUTS:
C
C       TYFR     FRACTION OF THE TROPIC YEAR, IN  THE  RANGE
C                0. TO 1., STARTING ON JAN. 1ST.
C       SUDC     SUN DECLINATION IN RADIANS (-PI TO PI).
C       RLAT     DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
C                POINT, IN RADIANS (-PI TO PI).
C       ALTU     GEOCENTRIC ALTITUDE  IN KM,  IN  THE  RANGE
C                90. TO 2000..
C
C OUTPUTS:
C
C       AL       ARRAY CONTAINING  THE  SEAZONAL-LATITUDINAL
C                VARIATIONS FOR THE HE (HELIUM), O2 (MOLECU-
C                LAR OXYGEN), N2  (MOLECULAR  NITROGEN),  AR
C                (ARGON), O (ATOMIC OXYGEN)  AND  H  (ATOMIC
C                HYDROGEN) NUMBER DENSITY, RESPECTIVELY.
C
C
C REFERENCES:
C
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA       APR/87             V 1.0
C            BENTO SILVA DE MATOS
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION AL(6)
      DIMENSION CR(6)
C
C------
C  
      DATA PITW /6.28318530718D0/
  
      DATA CR   /-0.79D0,0.D0,0.D0,0.D0,-.16D0,0.D0/
  
      SILA   = DSIN(RLAT)
      DSLT   = SUDC*SILA/0.409157536545D0
      DELZ   = ALTU - 91.D0
      ESSE   = 0.014D0*DELZ*DEXP(-0.0013D0*DELZ*DELZ)
      PCAP   = DSIN(PITW*TYFR + 1.72D0)
      DSLM   = DSIGN(SILA*SILA*ESSE*PCAP,RLAT)
  
      DO   I = 1, 6
       AL(I) = DSLT*CR(I) + DSLM
      ENDDO
  
      RETURN
      END
  
  
      SUBROUTINE SEMIAN(TYFR,ALTU,ALCO)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE SEMIAN GIVES THE CORRECTION FACTOR ALCO
C
C     FOR   THE  ATMOSPHERE   NUMBER  DENSITY,  DUE  TO  THE
C
C     SEMIANNUAL EFFECT.
C     ------
C
C INPUTS:
C
C       TYFR     FRACTION OF THE TROPIC  YEAR,  IN THE RANGE
C                0. TO 1., STARTING ON JAN. 1ST.
C       ALTU     GEOCENTRIC ALTITUDE  IN KM,  IN  THE  RANGE
C                90. TO 2000..
C
C OUTPUTS:
C
C       ALCO     THE SEMIANNUAL VARIATION OF THE  ATMOSPHERE
C                NUMBER DENSITY.
C
C REFERENCES:
C
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA        APR/87            V 1.0
C            BENTO SILVA DE MATOS
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C------
C  
      DATA PITW /6.28318530718D0/
  
      AUXI   = 0.04D0*ALTU*ALTU/1.D+04 + 0.05D0
      FOFT   = AUXI*DEXP(-0.25D-02*ALTU)
      AUXI   = (0.5D0+0.5D0*DSIN(PITW*TYFR+6.04D0))**1.65D0
      TAUC   = 0.0954D0*(AUXI - 0.5D0) + TYFR
      AUXI   = DSIN(2.D0*PITW*TAUC + 4.26D0)*(1.D0 +
     1         0.467D0*DSIN(PITW*TAUC + 4.14D0))
      GOFT   = AUXI*0.382D0 + 0.0284D0
      ALCO   = FOFT*GOFT
  
      RETURN
      END

      FUNCTION DATANH(X)
C
C------
C
C PURPOSE:
C
C     THE DATANH FUNCTION CALCULATES THE HYPERBOLIC  TANGENT
C     ARC OF AN ARGUMENT X, IN DOUBLE PRECISION.
C
C INPUTS:
C       X        HYPERBOLIC TANGENT VALUE (-1<X<1).
C
C OUTPUTS:
C       DATANH   HYPERBOLIC TANGENT ARC.
C
C AUTHOR:    VALDEMIR CARRARA  APR/87    V 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C------
C
      IF(DABS(X).GT.1.D0) THEN
        WRITE(*) 'INVALID HYPERBOLIC TANGENT ARC ARGUMENT'
        STOP
      END IF
  
      ARGU   = (1.D0 + X)/(1.D0 - X)

      DATANH = LOG(ARGU)/2.D0

      RETURN
      END
