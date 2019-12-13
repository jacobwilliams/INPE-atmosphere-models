C     LIBRARY JAC70
C
C----
C
C     The library JAC70 includes several routines to
C     compute the high atmospheric properties, using the
C     Jacchia 1970 model.
C
C----
C
      SUBROUTINE JDYMOS(SA,SU,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C------
C
C
C PURPOSE:
C
C     THE SUBROUTINE JDYMOS GIVES THE  TEMPERATURE,  DENSITY
C
C     AND MOLECULAR WEIGHT OF ATMOSPHERE, USING THE  JACCHIA
C                                                    -
C     70 INTEGRATED  STATIC  AND  DYNAMIC  MODEL,  WITH  THE  
C                                 --       -- 
C     SOLAR  FLUX DATA FILE.
C     -
C
C INPUTS:
C
C     SA(1) RIGHT ASCENSION OF THE  POINT  IN  QUESTION,  IN
C           RADIANS (0 TO 2.*PI)
C     SA(2) DECLINATION (GEOCENTRIC LATITUDE) OF THE  POINT,
C           IN RADIANS (-PI TO PI).
C     SA(3) GEOCENTRIC ALTITUDE OF THE POINT IN  METERS,  IN
C           THE RANGE 110000. TO 2000000.M.
C     SU(1) RIGHT ASCENSION OF THE SUN AT THE DATE, IN RADI-
C           ANS (O TO 2.*PI).
C     SU(2) SUN DECLINATION IN RADIANS (-PI TO PI)
C     RJUD  MODIFIED JULIAN DATE (IF OUT OF RANGE, THE  SUB-
C           ROUTINE WILL PRINT A MESSAGE AND STOP).
C     DAFR  TIME (UT) OF THE DAY, IN SECONDS.
C     GSTI  GREENWICH SIDERAL TIME, IN RADIANS (0 TO 2.*PI),
C           AT THE TIME DAFR OF THE DATE RJUD. (NOT USED. FOR
C           COMPATIBILITY PURPOSE WITH OTHER MODELS ONLY)
C
C OUTPUTS:
C
C     TE(1) EXOSPHERIC  TEMPERATURE  ABOVE  THE   POINT   IN 
C           QUESTION, AS DEFINED IN REFERENCE (1), IN KELVIN
C     TE(2) LOCAL TEMPERATURE AROUND THE POINT, IN KELVIN.
C     AD(1) LOGARITHM IN BASE 10 OF THE HE NUMBER-DENSITY.
C     AD(2) LOGARITHM IN BASE 10 OF THE O2 NUMBER-DENSITY.
C     AD(3) LOGARITHM IN BASE 10 OF THE N2 NUMBER-DENSITY.
C     AD(4) LOGARITHM IN BASE 10 OF THE AR NUMBER-DENSITY.
C     AD(5) LOGARITHM IN BASE 10 OF THE O  NUMBER-DENSITY.
C     AD(6) LOGARITHM IN BASE 10 OF THE H  NUMBER-DENSITY.
C     WMOL  MEAN-MOLECULAR-WEIGHT OF THE ATMOSPHERE  AT  THE
C           POINT IN KG/KGMOL.
C     RHOD  MEAN-MASS-DENSITY OF THE ATMOSPHERE, IN KG/M/M/M
C
C OBS:
C
C     HE    HELIUM
C     O2    MOLECULAR OXYGEN
C     N2    MOLECULAR NITROGEN
C     AR    ARGON
C     O     ATOMIC OXYGEN
C     H     ATOMIC HYDROGEN
C
C SUBCALLS:
C
C     SOFLUD
C     JSDAMO
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1989              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SA(3),SU(2),TE(2),AD(6)
      DIMENSION SF(3),SD(15)
C
C------
C

      RJFL    = RJFL - 1.
      IF(DAFR.LT.61344.) THEN
       RJFL   = RJFL - 1.
       DAFL   = DAFR + 25056.
      ELSE
       DAFL   = DAFR - 61344.
      ENDIF

      CALL SOFLUD(RJFL,DAFL,SD,OUTR)

      IF(OUTR.NE.0.) THEN
       WRITE(6,*) ' ERROR IN ROUTINE JDYMOS: OUTR = ',
     1            INT(OUTR)
       STOP
      ENDIF

      TAUO    = 6.696
      ND      = (DAFR/3600. - SD(6) + 12. - TAUO)/3.
      SF(1)   = SD(9)
      SF(2)   = SD(11)
      SF(3)   = SD(ND)

      CALL JSDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,RHOD)

      RETURN
      END

      SUBROUTINE JSMODS(ALTU,RJUD,DAFR,TE,AL,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE  SUBROUTINE  JSMODS  USES  THE  JACCHIA 70  STATIC 
C                                         -           -
C     MODEL (1) AND  THE  SOLAR  FLUX  DATA  TO  OBTAIN  THE 
C     ---                 -
C     TEMPERATURE, MOLECULAR WEIGHT  AND  DENSITY  OF  LOCAL
C
C     ATMOSPHERE.
C
C INPUTS:
C
C     ALTU  ALTITUDE  OF  THE  POINT  IN  METERS  (110000 TO
C           2000000).
C     RJUD  MODIFIED  JULIAN  DATE  (IF  OUT  OF  RANGE  THE 
C           ROUTINE WILL PRINT A MESSAGE AND STOP).
C     DAFR  TIME (UT) OF THE DAY, IN SECONDS (0  TO  86400).
C
C OUTPUTS:
C
C     TE(1) EXOSPHERIC  TEMPERATURE  ABOVE  THE   POINT   IN
C           QUESTION, AS DEFINED IN REFERENCE (1), IN KELVIN
C     TE(2) LOCAL TEMPERATURE AROUND THE POINT, IN KELVIN.
C     AD(1) LOGARITHM IN BASE 10 OF THE HE NUMBER-DENSITY.
C     AD(2) LOGARITHM IN BASE 10 OF THE O2 NUMBER-DENSITY.
C     AD(3) LOGARITHM IN BASE 10 OF THE N2 NUMBER-DENSITY.
C     AD(4) LOGARITHM IN BASE 10 OF THE AR NUMBER-DENSITY.
C     AD(5) LOGARITHM IN BASE 10 OF THE O  NUMBER-DENSITY.
C     AD(6) LOGARITHM IN BASE 10 OF THE H  NUMBER-DENSITY.
C     WMOL  MEAN-MOLECULAR-WEIGHT OF THE ATMOSPHERE  AT  THE 
C           POINT IN KG/KGMOL.
C     RHOD  MEAN-MASS-DENSITY OF THE ATMOSPHERE, IN KG/M/M/M
C
C SUBCALLS:
C
C     SOFLUD
C     JSMADE
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1989              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TE(2),AL(6)
      DIMENSION SF(3),SD(15)
C
C------
C
      CALL SOFLUD(RJUD,DAFR,SD,OUTR)

      IF(OUTR.NE.0.) THEN
       WRITE(6,*) ' ERROR IN ROUTINE ISMODS: OUTR = ',
     1            INT(OUTR)
       STOP
      ENDIF

      SF(1)   = SD(9)
      SF(2)   = SD(11)
      VARI    = .154*SD(7)
      SF(3)   = 1.89*DLOG(VARI + DSQRT(VARI*VARI + 1.D0))

      CALL JSMADE(ALTU,SF,TE,AL,WMOL,RHOD)

      RETURN
      END


      SUBROUTINE JSDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  JSDAMO  GIVES  THE  DENSITY, MOLECULAR
C
C     WEIGHT AND TEMPERATURE OF THE UPPER ATMOSPHERE,  USING
C
C     THE JACCHIA 70 STATIC  MODEL  (JACCHIA 1972)  AND  THE
C         -          -
C     DYNAMIC ATMOSPHERIC MODEL.
C     -       -           --
C
C INPUTS:
C
C       SA(1)    RIGHT ASCENTION OF THE POINT, IN RADIANS.
C       SA(2)    DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
C                POINT, IN RADIANS (-PI TO PI).
C       SA(3)    GEOCENTRIC ALTITUDE IN METERS, BETWEEN  THE
C                RANGE 110,000-2,000,000.
C       SU(1)    RIGHT ASCENTION OF THE SUN AT THE DATE,  IN
C                RADIANS (0 TO 2*PI).
C       SU(2)    SUN DECLINATION IN RADIANS (-PI TO PI).
C       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
C                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
C                W/M/M/HZ.
C       SF(2)    AVERAGED DAILY OBSERVED FLUX  AS DEFINED BY
C                JACCHIA, IN 1E-22 W/M/M/HZ.
C       SF(3)    3-HOURLY PLANETARY GEOMAGNETIC INDEX KP, AT
C                THE TIME 0.279 DAYS EARLIER.
C       RJUD     MODIFIED JULIAN  DATE,  REFERED  TO  1950.0
C                (JULIAN DATE-2433282.5).
C       DAFR     TIME (UT) OF THE DAY, IN SECONDS.
C       GSTI     GREENWICH SIDEREAL TIME, IN RADIANS, AT THE
C                TIME DAFR OF THE DATE RJUD (0 TO 2*PI).(NOT
C                USED. FOR COMPATIBILITY PURPOSE WITH  OTHER
C                MODELS ONLY)
C
C OUTPUTS:
C
C       TE(1)    EXOSPHERIC TEMPERATURE ABOVE THE  POINT  AS
C                DEFINED BY JACCHIA'S 70 MODEL, IN KELVIN.
C       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
C                KELVIN.
C       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
C       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
C                THE POINT, IN KG/KGMOL.
C       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
C                POINT, IN KG/M/M/M.
C OBS:
C       HE       HELIUM
C       O2       MOLECULAR OXYGEN
C       N2       MOLECULAR NITROGEN
C       AR       ARGON
C       O        ATOMIC OXYGEN
C       H        ATOMIC HYDROGEN
C
C SUBCALLS:
C
C       ADEN
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1989              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION SA(3),SU(2),SF(3),TE(2),AD(6)
      DIMENSION AL(6),ST(3)
  
C
C------
C  
      DATA CONV /2.302585093/

      AMJD   = RJUD + 33282. + DAFR/86400.
      ST(1)  = SA(1)
      ST(2)  = SA(2)
      ST(3)  = SA(3)/1000.

      CALL ADEN (AMJD,SU,ST,SF,TE,AL,WMOL,RHOD)
  
      AD(1)  = AL(5)/CONV
      AD(2)  = AL(2)/CONV
      AD(3)  = AL(1)/CONV
      AD(4)  = AL(4)/CONV
      AD(5)  = AL(3)/CONV
      AD(6)  = AL(6)/CONV
  
      RETURN
      END
  
  
      SUBROUTINE JSMADE(ALTU,SF,TE,AD,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
C
C     HEIGHTS FROM 100 TO 2000 KM, USING THE JACCHIA 70 STA-
C                                            -          -
C     TIC MODEL TO COMPUTE THE ATMOSPHERIC DENSITY (1).
C         -                    -           --
C
C INPUTS:
C
C       ALTU     ALTITUDE OF THE POINT IN METERS.
C       SF(1)    DAILY OBSERVED SOLAR FLUX AT  10.7  CM,  AT
C                THE  TIME  1.71  DAYS  EARLIER,  IN   1E-22
C                W/M/M/HZ.
C       SF(2)    AVERAGED DAILY OBSERVED FLUX  AS DEFINED BY
C                JACCHIA, IN 1E-22 W/M/M/HZ.
C
C OUTPUTS:
C
C       TE(1)    EXOSPHERIC TEMPERATURE ABOVE THE  POINT  AS
C                DEFINED BY JACCHIA'S 70 MODEL, IN KELVIN.
C       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
C                KELVIN.
C       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
C       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
C                THE POINT, IN KG/KGMOL.
C       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
C                POINT, IN KG/M/M/M.
C
C OBS:
C       HE       HELIUM
C       O2       MOLECULAR OXYGEN
C       N2       MOLECULAR NITROGEN
C       AR       ARGON
C       O        ATOMIC OXYGEN
C       H        ATOMIC HYDROGEN
C
C SUBCALLS:
C
C       JMOWEI
C
C OUTPUTS:
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1987              V. 1.0
C     AUG  2011			   V. 1.1 (INCLUDED MISSING TEMPERATURE)
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION SF(3),TE(2),AD(6)
      DIMENSION WM(6),AL(6)
C
C------
C
      DATA CONV /2.302585093/
      DATA AVOG /6.02217D+26/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/
  
      FLUX   = SF(1)
      FBAR   = SF(2)
      THAF   = 379.0 + 3.24*FBAR + 1.3*(FLUX - FBAR)
      HEIG   = ALTU/1.D3
  
      CALL SMADEN (THAF, HEIG, TZ, AL)

      AD(1)  = AL(5)/CONV
      AD(2)  = AL(2)/CONV
      AD(3)  = AL(1)/CONV
      AD(4)  = AL(4)/CONV
      AD(5)  = AL(3)/CONV
      AD(6)  = AL(6)/CONV
      ANUT   = 0.
      WEIG   = 0.

      DO IC = 1 , 6
       ANAC   = 10.**AD(IC)
       WEIG   = WEIG + WM(IC) * ANAC
       ANUT   = ANUT + ANAC
      ENDDO

      WMOL = WEIG/ANUT
      RHOD = WEIG/AVOG
	TE(1) = THAF
	TE(2) = TZ

      RETURN
      END

      SUBROUTINE JMOWEI (TINF, HEIG, AD, WMOL, RHOD)
C
C------
C
C PURPOSE:
C  
C     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
C
C     HEIGHTS FROM 100 TO 2000 KM, USING THE JACCHIA 70 STA-
C                                            -          -
C     TIC MODEL TO COMPUTE THE ATMOSPHERIC DENSITY AND MOLE-
C         -
C     CULAR WEIGHT (1).
C           ---
C
C INPUTS:
C
C       TINF     EXOSPHERIC  TEMPERATURE   AS   DEFINED   BY
C                JACCHIA'S 1970 MODEL, IN KELVIN.
C       HEIG     ALTITUDE OF THE POINT IN KM.
C
C OUTPUTS:
C
C       TLOC     LOCAL TEMPERATURE, IN KELVIN
C       AD(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AD(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AD(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AD(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AD(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AD(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
C       WMOL     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE  AT
C                THE POINT, IN KG/KGMOL.
C       RHOD     MEAN MASS DENSITY OF THE ATMOSPHERE AT  THE
C                POINT, IN KG/M/M/M.
C
C OBS:
C       HE       HELIUM
C       O2       MOLECULAR OXYGEN
C       N2       MOLECULAR NITROGEN
C       AR       ARGON
C       O        ATOMIC OXYGEN
C       H        ATOMIC HYDROGEN
C
C SUBCALLS:
C
C       SMADEN
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1989              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION AD(6)
      DIMENSION AL(6),WM(6)
C
C------
C
      DATA CONV /2.302585093/
      DATA AVOG /6.02217D+26/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/

      CALL SMADEN (TINF, HEIG, TZ, AL)

      AD(1)  = AL(5)/CONV
      AD(2)  = AL(2)/CONV
      AD(3)  = AL(1)/CONV
      AD(4)  = AL(4)/CONV
      AD(5)  = AL(3)/CONV
      AD(6)  = AL(6)/CONV
      ANUT   = 0.
      WEIG   = 0.

      DO IC = 1 , 6
       ANAC   = 10.**AD(IC)
       WEIG   = WEIG + WM(IC) * ANAC
       ANUT   = ANUT + ANAC
      ENDDO

      WMOL = WEIG/ANUT
      RHOD = WEIG/AVOG

      END


      SUBROUTINE ADEN (AMJD,SUN,SAT,GEO,TEMP,ALIN,AMMW,RHO)

C
C------
C
C PURPOSE:
C
C     THE ROUTINE ADEN COMPUTES THE ATMOSPHERIC  PROPERTIES:
C                                   -
C     TEMPERATURE, NUMBER OF DENSITY,  MEAN  MOLECULAR  MASS
C 
C     DENSITY. IT WAS DEVELOPED INITIALLY BY JACCHIA (1).
C     ---
C
C INPUTS:
C
C     AMJD     DATE AND TIME, IN MODIFIED  JULIAN  DAYS  AND
C              FRACTION THEREOF (MJD=JD-2400000.5).
C     SUN(1)   RIGHT ASCENSION OF THE SUN, IN RADIANS
C     SUN(2)   DECLINATION OF THE SUN, IN RADIANS
C     SAT(1)   RIGHT ASCENSION OF  THE  POINT  IN  QUESTION,
C              IN RADIANS
C     SAT(2)   DECLINATION  (GEOCENTRIC  LATITUDE)  OF   THE
C              POINT IN QUESTION, IN RADIANS
C     SAT(3)   HEIGHT OF THE POINT IN QUESTION, IN KM.
C     GEO(1)   10.7 CM  SOLAR  FLUX,  IN  UNITS  OF  1.0E-22
C              WATT*(M**-2)*(HERTZ**-1), FOR A TABULAR  TIME
C              1.71 DAYS EARLIER.
C     GEO(2)   10.7 CM SOLAR FLUX, AVERAGED OVER FOUR  SOLAR
C              ROTATIONS, CENTERED ON THE TIME IN QUESTION
C     GEO(3)   THE  GEOMAGNETIC  PLANETARY  THREE-HOUR-RANGE
C              INDEX KP (3- = 2.667, 30 = 3.000, 3+ = 3.333, 
C              ETC.), FOR A TABULAR TIME 0.279 DAYS EARLIER.
C
C OUTPUTS:
C
C     TEMP(1)  EXOSPHERIC TEMPERATURE  ABOVE  THE  POINT  IN
C              QUESTION, KELVIN.
C     TEMP(2)  TEMPERATURE AT THE POINT IN QUESTION, KELVIN.
C     ALIN(1)  COMMON LOGARITHM OF THE N2 NUMBER-DENSITY, IN
C              M**-3
C     ALIN(2)  COMMON LOGARITHM OF THE 02 NUMBER-DENSITY
C     ALIN(3)  COMMON LOGARITHM OF THE 0  NUMBER-DENSITY
C     ALIN(4)  COMMON LOGARITHM OF THE A  NUMBER-DENSITY
C     ALIN(5)  COMMON LOGARITHM OF THE HE NUMBER-DENSITY
C     ALIN(6)  COMMON LOGARITHM OF THE H  NUMBER-DENSITY
C     AMMW     MEAN-MOLECULAR-WEIGHT
C     RHO      TOTAL MASS-DENSITY, IN KG*(M**-3)
C
C OBS:
C       N2       MOLECULAR NITROGEN
C       O2       MOLECULAR OXYGEN
C       O        ATOMIC OXYGEN
C       A        ARGON
C       HE       HELIUM
C       H        ATOMIC HYDROGEN
C
C SUBCALLS:
C
C       SMADEN
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHORS:
C
C     VALDEMIR CARRARA (BASED ON THE ADEN SUBRROUTINE).
C
C DATE:
C
C     APR. 1989              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SUN(2), SAT(3), GEO(3), TEMP(2), ALIN(6)
      DIMENSION ALN(6), AMW(6)
C
C------
C
C          AL10 IS ALOG(1.0)
C          THE-AMW ARE THE MOLECULAR WEIGHTS IN THE ORDER
C          N2,02,0,A,HE, AND H
C          AVOGAD IS AVOGADROS NUMBER IN MKS UNITS
C          CONS25 IS SIN(PI/4.0)**3, USED IN EQUATION (25)
C          FOURPI IS 4.0*PI
C          TWOPI IS 2.0*PI
C          PIOV4 IS PI/4.0

      DATA AL10   /2.302585093/
      DATA AMW    /28.0134, 31.9988, 15.9994,
     1              39.948,  4.0026, 1.00797/
      DATA AVOGAD /6.02257E26/
      DATA CONS25 /0.35355339/
      DATA FOURPI /12.56637062/
      DATA TWOPI  /6.283185308/
      DATA PIOV4  /0.7853981635/

      TSUBC   = 379.0 + 3.24*GEO(2) + 1.3*(GEO(1)-GEO(2))
      ETA     = 0.5*ABS(SAT(2)-SUN(2))
      THETA   = 0.5*ABS(SAT(2)+SUN(2))
      H       = SAT(1) - SUN(1)
      TAU     = H - 0.64577182 + 0.10471976*SIN(H+0.75049158)
      C       = COS(ETA)**2.2
      S       = SIN(THETA)**2.2
      DF      = S + (C-S)*ABS(COS(0.5*TAU))**3
      TSUBL   = TSUBC*(1.0+0.3*DF)
      EXPKP   = EXP(GEO(3))
      DTG18   = 28.0*GEO(3) + 0.03*EXPKP
      DTG20   = 14.0*GEO(3) + 0.02*EXPKP
      DLR20   = 0.012*GEO(3) + 1.2E-5*EXPKP
      F       = 0.5*(TANH(0.04*(SAT(3)-350.0))+1.0)
      DLRGM   = DLR20*(1.0-F)
      DTG     = DTG20*(1.0-F) + DTG18*F
      TINF    = TSUBL+DTG

      CALL SMADEN (TINF, SAT(3), TZ, ALN)

      TEMP(1) = TINF
      TEMP(2) = TZ
     
      CAPPHI  = MOD((AMJD-36204.0D0)/365.2422D0,1.0D0)
      TAU     = CAPPHI + 0.09544*((0.5 +
     1          0.5*SIN(TWOPI*CAPPHI+6.035))**1.650-0.5)
      GOFT    = 0.02835 + 0.3817*(1.0 +
     1          0.4671*SIN(TWOPI*TAU+4.137))*
     2          SIN(FOURPI*TAU+4.259)
      FOFZ    = (5.876E-7*SAT(3)**2.331 + 0.06328)*
     1          EXP(-2.868E-3*SAT(3))
      DLRSA   = FOFZ*GOFT
      DLRSL   = 0.014*(SAT(3)-90.0)*
     1          EXP(-0.0013*(SAT(3)-90.0)**2)*
     2          SIGN(1.0D0,SAT(2))*SIN(TWOPI*CAPPHI+1.72)*
     3          SIN(SAT(2))**2
      DLR     = AL10*(DLRGM + DLRSA + DLRSL)

      DO I = 1, 6
       ALN(I) = ALN(I) + DLR
      ENDDO

      DLNHE   = 0.65*ABS(SUN(2)/0.4091609D0)*
     1          (SIN(PIOV4-0.5*SAT(2)*
     2          SIGN(1.0D0,SUN(2)))**3-CONS25)
      ALN(5)  = ALN(5) + AL10*DLNHE
      SUMN    = 0.0
      SUMNM   = 0.0

      DO I = 1, 6
       AN     = EXP(ALN(I))
       SUMN   = SUMN + AN
       SUMNM  = SUMNM + AN*AMW(I)
       ALIN(I)= ALN(I)/AL10
      ENDDO

      AMMW    = SUMNM/SUMN
      RHO     = SUMNM/AVOGAD

      RETURN
      END

      SUBROUTINE SMADEN (TINF, SAT3, TZ, ALN)
C
C-----
C
C PURPOSE:
C
C     THE SUBROUTINE SMADEN OBTAINS THE STATIC CONDITIONS  OF
C                                       -
C     THE UPPER ATMOSPHERE USING THE JACCHIA 70 MODEL (ADEN).
C                                               -      ----
C
C INPUTS:
C
C     TINF  EXOSPHERIC  TEMPERATURE  ABOVE  THE   POINT   IN
C           QUESTION, IN KELVIN.
C     SAT3  ALTITUDE OF THE POINT IN QUESTION, IN KILOMETERS
C
C OUTPUTS:
C
C     TZ    TEMPERATURE AT THE POINT IN QUESTION, IN KELVIN.
C     ALN(1) COMMON LOG OF THE N2 NUMBER DENSITY, IN M**-3.
C     ALN(2) COMMON LOG OF THE O2 NUMBER DENSITY, IN M**-3.
C     ALN(3) COMMON LOG OF THE  O NUMBER DENSITY, IN M**-3.
C     ALN(4) COMMON LOG OF THE  A NUMBER DENSITY, IN M**-3.
C     ALN(5) COMMON LOG OF THE HE NUMBER DENSITY, IN M**-3.
C     ALN(6) COMMON LOG OF THE  H NUMBER DENSITY, IN M**-3.
C
C SUBCALLS:
C
C     AMBAR
C     TLOCAL
C     GRAV
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHOR:
C
C     VALDEMIR CARRARA (BASED ON THE ADEN SUBROUTINE)
C
C DATE:
C
C     ABR. 89             V 1.0
C

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALN(6), ALPHA(5), AMW(6), FRAC(4)
      DIMENSION TC(4), WT(5)

C
C-----
C
C          THE ALPHA ARE THE THERMAL DIFFUSION COEFFICIENTS IN
C          EQUATION (6)
C          AL10 IS ALOG(1.0)
C          THE-AMW ARE THE MOLECULAR WEIGHTS IN THE ORDER
C          N2,02,0,A,HE, AND H
C          THE FRAC ARE THE ASSUMED SEA-LEVEL VOLUME-FRACTIONS IN THE
C          ORDER N2,02,A, AND HE
C          AVOGAD IS AVOGADROS NUMBER IN MKS UNITS
C          PIOV2 IS PI/2.0
C          RSTAR IS THE UNIVERSAL GAS-CONSTANT IN MKS UNITS
C          THE R ARE VALUES USED TO ESTABLISH HEIGHT-STEP-SIZES IN
C          THE REGIMES 90 KM TO 100 KM; 100 KM TO 500 KM, AND
C          500 KM UPWARDS
C          THE WT ARE WEIGHTS FOR THE NEWTON-COTES FIVE-POINT
C          QUADRATURE FORMULA

      DATA ALPHA  /0.0, 0.0, 0.0, 0.0, -0.38/
      DATA AL10   /2.3025851/
      DATA AMW    /28.0134, 31.9988, 15.9994,
     1              39.948,  4.0026, 1.00797/
      DATA FRAC   /0.78110, 0.20955, 9.3432E-3, 6.1471E-6/
      DATA AVOGAD /6.02257E26/
      DATA PIOV2  /1.5707963/
      DATA RSTAR  /8314.32/
      DATA R1     /0.010/
      DATA R2     /0.025/
      DATA R3     /0.075/
      DATA WT     /0.3111111111, 1.4222222222, 0.5333333333, 
     *             1.4222222222, 0.3111111111/


      TSUBX   = 371.6678 + 0.0518806*TINF - 
     1          294.3503*EXP(-0.00216222*TINF)
      GSUBX   = 0.054285714*(TSUBX-183.0)
      TC(1)   = TSUBX
      TC(2)   = GSUBX
      TC(3)   = (TINF-TSUBX)/PIOV2
      TC(4)   = GSUBX/TC(3)
      Z1      = 90.0D0
      Z2      = MIN(SAT3,100.0D0)
      AL      = LOG(Z2/Z1)
      N       = INT(AL/R1) + 1
      ZR      = EXP(AL/FLOAT(N))
      AMBAR1  = AMBAR(Z1)
      TLOC1   = TLOCAL(Z1,TC)
      ZEND    = Z1
      SUM2    = 0.0
      AIN     = AMBAR1*GRAV(Z1)/TLOC1

      DO I = 1, N
       Z      = ZEND
       ZEND   = ZR*Z
       DZ     = 0.25*(ZEND-Z)
       SUM1   = WT(1)*AIN

       DO J = 2, 5
        Z     = Z + DZ
        AMBAR2= AMBAR(Z)
        TLOC2 = TLOCAL(Z,TC)
        GRAVL = GRAV(Z)
        AIN   = AMBAR2*GRAVL/TLOC2
        SUM1  = SUM1 + WT(J)*AIN
       ENDDO

       SUM2   = SUM2 + DZ*SUM1
      ENDDO
      FACT1   = 1000.0/RSTAR
      RHO     = 3.46E - 6*AMBAR2*TLOC1*
     1          EXP(-FACT1*SUM2)/AMBAR1/TLOC2
      ANM     = AVOGAD*RHO
      AN      = ANM/AMBAR2
      FACT2   = ANM/28.960
      ALN(1)  = LOG(FRAC(1)*FACT2)
      ALN(4)  = LOG(FRAC(3)*FACT2)
      ALN(5)  = LOG(FRAC(4)*FACT2)
      ALN(2)  = LOG(FACT2*(1.0+FRAC(2))-AN)
      ALN(3)  = LOG(2.*(AN-FACT2))

      IF (SAT3.LE.100.0) THEN
       TZ     = TLOC2
       ALN(6) = ALN(5) - 25.0
      ELSE
       Z3     = MIN(SAT3,500.0D0)
       AL     = LOG(Z3/Z)
       N      = INT(AL/R2) + 1
       ZR     = EXP(AL/FLOAT(N))
       SUM2   = 0.0
       AIN    = GRAVL/TLOC2

       DO I = 1, N
        Z     = ZEND
        ZEND  = ZR*Z
        DZ    = 0.25*(ZEND-Z)
        SUM1  = WT(1)*AIN

        DO J = 2, 5
         Z    = Z + DZ
         TLOC3= TLOCAL(Z,TC)
         GRAVL= GRAV(Z)
         AIN  = GRAVL/TLOC3
         SUM1 = SUM1 + WT(J)*AIN
        ENDDO

        SUM2  = SUM2 + DZ*SUM1
       ENDDO

       Z4     = MAX(SAT3,500.0D0)
       AL     = LOG(Z4/Z)
       R      = R2
       IF (SAT3.GT.500.0) R = R3
       N      = INT(AL/R) + 1
       ZR     = EXP(AL/FLOAT(N))
       SUM3   = 0.0

       DO I = 1, N
        Z     = ZEND
        ZEND  = ZR*Z
        DZ    = 0.25*(ZEND-Z)
        SUM1  = WT(1)*AIN
 
        DO J = 2, 5
         Z    = Z + DZ
         TLOC4= TLOCAL(Z,TC)
         GRAVL= GRAV(Z)
         AIN  = GRAVL/TLOC4
         SUM1 = SUM1 + WT(J)*AIN
        ENDDO

        SUM3  = SUM3 + DZ*SUM1
       ENDDO

       IF (SAT3.LE.500.0) THEN
        T500  = TLOC4
        TZ    = TLOC3
        ALTR  = LOG(TLOC3/TLOC2)
        FACT2 = FACT1*SUM2
        HSIGN = 1.0
       ELSE
        T500  = TLOC3
        TZ    = TLOC4
        ALTR  = LOG(TLOC4/TLOC2)
        FACT2 = FACT1*(SUM2+SUM3)
        HSIGN = -1.0
       ENDIF

       DO I = 1, 5
        ALN(I)= ALN(I) - (1.0+ALPHA(I))*ALTR - FACT2*AMW(I)
       ENDDO
       AL10T5 = LOG10(T500)
       ALNH5  = (5.5*AL10T5-39.40)*AL10T5 + 73.13
       ALN(6) = AL10*(ALNH5+6.0) + HSIGN*(LOG(TLOC4/TLOC3) +
     *          FACT1*SUM3*AMW(6))
      ENDIF

      RETURN
      END

      FUNCTION AMBAR(Z)
C
C-----
C
C PURPOSE:
C
C     TO EVALUATE THE EQUATION (1) OF REFERENCE (1).
C
C INPUTS:
C
C     Z     HEIGHT IN KM, IN THE RANGE 90 TO 100 KM.
C
C OUTPUTS:
C
C     AMBAR MEAN  MOLECULAR  MASS  OF  THE  ATMOSPHERE,  IN
C           KG/KGMOL.
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHOR:
C
C     VALDEMIR CARRARA (BASED ON THE ADEN SUBROUTINE)
C
C DATE:
C
C     ABR. 89             V 1.0
C

      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION C(7)
C
C------
C
      DATA C /28.82678E+0, -7.40066E-2, -1.19407E-2,
     1         4.51103E-4, -8.21895E-6,  1.07561E-5, 
     2        -6.97444E-7/

      DZ      = Z - 90.0
      AMB     = C(7)

      DO I = 1, 6
       J      = 7 - I
       AMB    = DZ*AMB + C(J)
      ENDDO

      AMBAR   = AMB

      RETURN
      END


      FUNCTION GRAV(Z)
C
C-----
C
C PURPOSE:
C
C     TO EVALUATE THE EQUATION (8) OF REFERENCE (1).
C
C INPUTS:
C
C     Z     HEIGHT IN KM.
C
C OUTPUTS:
C
C     GRAV  ACCELERATION DUE TO GRAVITY, IN M/S/S
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHOR:
C
C     VALDEMIR CARRARA (BASED ON THE ADEN SUBROUTINE)
C
C DATE:
C
C     ABR. 89             V 1.0
C

      IMPLICIT REAL*8 (A-H,O-Z)

C
C------
C
      DENO    = (1.0+Z/6356.766)
      GRAV    = 9.80665/DENO/DENO

      RETURN
      END


      FUNCTION TLOCAL(Z,TC)
C
C-----
C
C PURPOSE:
C
C     TO EVALUATE EQUATION (10) OR  (13)  OF  REFERENCE (1),
C     DEPENDING ON Z.
C
C INPUTS:
C
C     Z     HEIGHT IN KM.
C
C OUTPUTS:
C
C     TLOCAL  LOCAL TEMPERATURE OF THE ATMOSPHERE, IN KELVIN
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C AUTHOR:
C
C     VALDEMIR CARRARA (BASED ON THE ADEN SUBROUTINE)
C
C DATE:
C
C     ABR. 89             V 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TC(4)

C
C------
C
      DZ      = Z - 125.0
      IF(DZ.LE.0.0) THEN
       TLOCAL = ((-9.8204695E-6*DZ-7.3039742E-4)*
     1          DZ**2+1.0)*DZ*TC(2) + TC(1)
      ELSE
       TLOCAL = TC(1) + 
     1          TC(3)*ATAN(TC(4)*DZ*(1.0+4.5E-6*DZ**2.5))
      ENDIF

      RETURN
      END
