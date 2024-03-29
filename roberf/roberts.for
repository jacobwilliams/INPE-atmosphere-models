C     LIBRARY DENSITY
C
C----
C
C     The library DENSITY includes several routines to
C     compute the high atmospheric properties, using the
C     Robert's version of the Jacchia 1970 model.
C
C----
C
      SUBROUTINE RDYMOS (SA, SU, RJUD, DAFR, GSTI, TE, AD,
     1                   WMOL, RHOD)
C
C------
C
C
C PURPOSE:
C
C     THE SUBROUTINE RDYMOS GIVES THE  TEMPERATURE,  DENSITY
C
C     AND MOLECULAR WEIGHT OF ATMOSPHERE, USING THE  ROBERTS
C                                                    -
C     VERSION (2) OF THE JACCHIA 70 DYNAMIC MODEL, WITH  THE  
C                                   --      -- 
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
C           AT THE TIME DAFR OF THE DATE RJUD.(NOT USED. FOR
C           COMPATIBILITY PURPOSE WITH  OTHER MODELS ONLY)
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
C     RSDAMO
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
C         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
C         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
C         DEC. 1971.
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

      CALL RSDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,RHOD)

      RETURN
      END

      SUBROUTINE RSMODS (ALTU, RJUD, DAFR, TE, AL, WMOL, RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE RSMODS USES THE ROBERTS VERSION (2)  OF 
C                                    -
C     THE JACCHIA 70 STATIC AND DYNAMIC MODEL (1), WITH  THE
C                    -                  ---
C     SOLAR FLUX DATA TO OBTAIN THE  TEMPERATURE,  MOLECULAR
C     -
C     WEIGHT AND DENSITY OF LOCAL ATMOSPHERE.
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
C     RSMADE
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
C         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
C         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
C         DEC. 1971.
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

      CALL RSMADE(ALTU,SF,TE,AL,WMOL,RHOD)

      RETURN
      END


      SUBROUTINE RSDAMO (SA, SU, SF, RJUD, DAFR, GSTI, TE,
     1                   AD, WMOL, RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  RSDAMO  GIVES  THE  DENSITY, MOLECULAR
C
C     WEIGHT AND TEMPERATURE OF THE UPPER ATMOSPHERE,  USING
C
C     THE ROBERTS VERSION (2) OF THE JACCHIA 70  STATIC  AND 
C         -                                      -
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
C       DYJRMO
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
C         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
C         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
C         DEC. 1971.
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
      DIMENSION AL(6)
  
C
C------
C  

      AMJD   = RJUD + 33282. + DAFR/86400.

      CALL DYJRMO (AMJD,SU,SA,SF,TE,AL,WMOL,RHOD)
  
      AD(1)  = AL(3)
      AD(2)  = AL(4)
      AD(3)  = AL(1)
      AD(4)  = AL(2)
      AD(5)  = AL(5)
      AD(6)  = AL(6)
  
      RETURN
      END
  
  
      SUBROUTINE RSMADE (ALTU, SF, TE, AD, WMOL, RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
C
C     HEIGHTS FROM 110 TO 2000 KM, USING THE ROBERTS VERSION
C                                            -
C     (2) OF THE JACCHIA 70  STATIC  MODEL  TO  COMPUTE  THE
C                            -       -
C     ATMOSPHERIC DENSITY (1).
C     -           --
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
C       RMOWEI
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
C         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
C         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
C         DEC. 1971.
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     APR. 1987              V. 1.0
C     AUG. 2011              V. 1.1 (INCLUDED MISSING TEMPERATURE)
C
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION SF(3),TE(2),AD(6)
      DIMENSION WM(6),AL(6)
C
C------
C
      DATA AVOG /6.02217D+26/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/
  
      FLUX   = SF(1)
      FBAR   = SF(2)
      THAF   = 379.0 + 3.24*FBAR + 1.3*(FLUX - FBAR)
      HEIG   = ALTU/1.D3

      CALL STJRMO (THAF, HEIG, TZ, AL)

      AD(1)  = AL(3)
      AD(2)  = AL(4)
      AD(3)  = AL(1)
      AD(4)  = AL(2)
      AD(5)  = AL(5)
      AD(6)  = AL(6)
      ANUT   = 0.
      WEIG   = 0.
  
      DO IC = 1 , 6
       ANAC   = 10.**AD(IC)
       WEIG   = WEIG + WM(IC) * ANAC
       ANUT   = ANUT + ANAC
      ENDDO

      WMOL = WEIG/ANUT
      RHOD = WEIG/AVOG
	TE(1)= THAF
	TE(2)= TZ
  
      RETURN
      END

      SUBROUTINE RMOWEI (TINF, HEIG, AD, WMOL, RHOD)
C
C------
C
C PURPOSE:
C  
C     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
C
C     HEIGHTS FROM 110 TO 2000 KM, USING THE ROBERTS VERSION
C                                            -
C     (2) OF THE JACCHIA'S 70 STATIC MODEL  TO  COMPUTE  THE
C                                    --
C     ATMOSPHERIC DENSITY AND MOLECULAR WEIGHT (1).
C                                       ---
C
C INPUTS:
C
C       TINF     EXOSPHERIC  TEMPERATURE   AS   DEFINED   BY
C                JACCHIA'S 1970 MODEL, IN KELVIN.
C       HEIG     ALTITUDE OF THE POINT IN KM.
C
C OUTPUTS:
C
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
C       STJRMO
C
C REFERENCES:
C
C     (1) JACCHIA, L. G.  ATMOSPHERIC MODELS IN  THE  REGION
C         FROM  110  TO  2000 KM.  IN:  COMMITTEE  ON  SPACE 
C         RESEARCH (COSPAR) "CIRA 1972".  BERLIM,  AKADEMIK-
C         VERLAG, 1972. PART 3, P. 227-338.
C
C     (2) ROBERTS JR., C. E. AN ANALYTICAL MODEL  FOR  UPPER
C         ATMOSPHERIC DENSITIES BASED  UPON  JACCHIA'S  1970
C         MODELS.  "CELESTIAL  MECHANICS",   4(3/4):368-377,
C         DEC. 1971.
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
      DATA AVOG /6.02217D+26/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/

      CALL STJRMO (TINF, HEIG, TZ, AL)

      AD(1)  = AL(3)
      AD(2)  = AL(4)
      AD(3)  = AL(1)
      AD(4)  = AL(2)
      AD(5)  = AL(5)
      AD(6)  = AL(6)
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


      SUBROUTINE DYJRMO(DJM,SUN,SAT,GEO,TEMP,DN,AMW,DENS)
C
C-----
C
C  PURPOSE : COMPUTATION OF THE ATMOSPHERIC PROPERTIES
C            ACCORDING TO THE ANALYTICAL ROBERTS(1972)
C            METHOD APPLIED TO THE JACCHIA(1971) MODEL
C
C INPUTS : DJM ... MODIFIED JULIAN DATE DJM=JD-2400000.5
C          SUN(1). RIGHT ASCENSION OF SUN (RAD)
C          SUN(2). DECLINATION     OF SUN (RAD)
C          SAT(1). RIGHT ASCENSION OF THE POINT (RAD)
C          SAT(2). DECLINATION     OF THE POINT (RAD)
C          SAT(3). ALTITUDE        OF THE POINT (M)
C          GEO(1). 10.7 CM SOLAR FLUX,IN UNITS OF
C                  1.E-22 WATTS M**2/HERTZ , FOR A
C                  TABULAR TIME 1.71 DAYS EARLIER
C          GEO(2). 10.7 CM SOLAR FLUX AVERAGED OVER
C                  FOUR SOLAR ROTATIONS,CENTERED ON
C                  THE PRESENT TIME
C          GEO(3). GEOMAGNETIC PLANETARY THREE HOUR
C                  RANGE  INDEX "KP" FOR A TABULAR
C                  TIME 0.279 DAYS EARLIER
C
C OUTPUTS : TEMP(1)...EXOSPHERIC TEMPERATURE (KELVIN)
C           TEMP(2)...LOCAL TEMPERATURE      (KELVIN)
C           DN(1) ... LOG10 OF N2 DENSITY NUMBER (M**-3)
C           DN(2) ... LOG10 OF A  DENSITY NUMBER
C           DN(3) ... LOG10 OF HE DENSITY NUMBER
C           DN(4) ... LOG10 OF O2 DENSITY NUMBER
C           DN(5) ... LOG10 OF O  DENSITY NUMBER
C           DN(6) ... LOG10 OF H  DENSITY NUMBER
C           AMW   ... MEAN MOLECULAR WEIGHT (KG/KGMOL)
C           DENS  ... ATMOSPHERIC DENSITY (KG/M**3)
C
C REF. KUGA,H.K. "REFORMULACAO COMPUTACIONAL DO
C         MODELO DE JACCHIA-ROBERTS PARA A DEN-
C         SIDADE ATMOSFERICA".INPE,SAO JOSE DOS
C         CAMPOS, OCT.1985 (INPE-3691-RPE/493).
C
C      JACCHIA,L.G."REVISED STATIC MODELS OF THE
C         THERMOSPHERE AND EXOSPHERE WITH EMPIRICAL
C         TEMPERATURE PROFILES."SAO,CAMBRIDGE,MA,
C         1971.SAO SPECIAL REPORT NO. 332
C
C      ROBERTS JR,C.E."AN ANALYTIC MODEL FOR UPPER
C         ATMOSPHERE DENSITIES BASED UPON JACCHIA'S
C         1970 MODELS."CELESTIAL MECHANICS 4:368-
C         377,1971
C
C AUTHOR : HELIO KOITI KUGA - JUNE 1985 -INPE-DMC/DDO
C
C DIMENSION ARRAY'S, VARIABLES AND CONSTANTS
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SUN(2),SAT(3),GEO(3),TEMP(2),DN(6)
C-----
      DATA PIV2,PIV4,PID4,CONS25 /6.2831853D0,
     *     12.566371D0,0.78539816D0,0.35355339D0/
C
C      PIV2 = 2 * PI
C      PIV4 = 4 * PI
C      PID4 = PI / 4
C      CONS25 = SIN (PI/4) **3
C
      SUN1 = SUN(1)
      SUN2 = SUN(2)
      SAT1 = SAT(1)
      SAT2 = SAT(2)
      SAT3 = SAT(3)/1000.
      FS   = GEO(1)
      FSM  = GEO(2)
      PK   = GEO(3)
C
C       MINIMUM NIGHT-TIME TEMPERATURE OF THE GLOBAL
C       EXOSPHERIC TEMPERATURE DISTRIBUTION WHEN THE
C       GEOMAGNETIC ACTIVITY INDEX KP = 0
C       EQUATION 14J
C
      TSUBC = 379.+3.24*FSM+1.3*(FS-FSM)
C
C       EQUATION 15J
C
      ETA  = 0.5*ABS(SAT2-SUN2)
      THETA= 0.5*ABS(SAT2+SUN2)
C
C       EQUATION 16J
C
      H  = SAT1 - SUN1
      TAU = H - 0.64577182 + 0.10471976*SIN(H+0.75049158)
C
C       EXOSPHERIC TEMPERATURE TSUBL WITHOUT CORRECTION
C       FOR GEOMAGNETIC ACTIVITY
C       EQUATION 17J
C
      S    = SIN(THETA)**2.2
      DF   = S + ( COS(ETA)**2.2D0-S ) * ABS( COS(0.5*TAU) )**3
      TSUBL= TSUBC*(1.+0.3*DF)
C
C       EQUATION 18J
C
      EXPKP  = EXP(PK)
      DTG18  = 28.*PK + 0.03*EXPKP
C
C       EQUATION 20J
C
      DTG20 = 14.*PK + 0.02*EXPKP
      DLR20 = 0.012*PK + 1.2D-05*EXPKP
C
C       THE FOLLOWING STATEMENTS EFFECT A CONTINUOUS
C       TRANSITION FROM EQ. 20J AT HEIGHTS WELL BELOW
C       350 KM TO EQ. 18J AT HEIGHTS WELL ABOVE
C       350 KM .
C
      F = 0.5*(TANH(0.04*(SAT3-350.))+1.)
      DLRGM = DLR20 * (1.-F)
      DTG   = DTG20 * (1.-F) + DTG18 * F
C
C       EXOSPHERIC TEMPERATURE
C
      TINF  = TSUBL + DTG
C
C   STATIC MODEL OF JACCHIA-ROBERTS FOR THE
C   ATMOSPHERIC DENSITY
C
      CALL STJRMO(TINF,SAT3,TZ,DN)
C
C   EQ. 23J   PHASE OF THE SEMI-ANNUAL VARIATION
C
      CAPPHI = MOD((DJM-36204.)/365.2422,1.D0)
C
C   EQ. 22J
C
      TAU = CAPPHI + 0.09544*(
     *      (0.5+0.5*SIN(PIV2*CAPPHI+6.035))**1.650-0.5)
      GDFT = 0.02835+0.3817*(1.+0.4671*SIN(PIV2*TAU+4.137))
     *              *SIN(PIV4*TAU+4.259)
      FDFZ = (5.876D-07*SAT3**2.331D0+0.06328)*
     *       EXP(-2.868D-03*SAT3)
C
C   EQ. 21J  SEMI-ANNUAL VARIATION
C
      DLRSA = FDFZ*GDFT
C
C   EQ. 24J  SEASONAL-LATITUDINAL VARIATION OF THE
C            LOWER THERMOSPHERE
C
      DLRSL = 0.014*(SAT3-90.)*EXP(-0.0013*(SAT3-90.)**2)
     *       *SIGN(1.D0,SAT2)*SIN(PIV2*CAPPHI+1.72)
     *       *SIN(SAT2)**2
C
C   SUM THE CORRECTIONS AND APPLY TO THE
C   NUMBER DENSITIES
C
      DLR = DLRGM + DLRSA + DLRSL
      DN(1) = DN(1) +  DLR
      DN(2) = DN(2) +  DLR
      DN(3) = DN(3) +  DLR
      DN(4) = DN(4) +  DLR
      DN(5) = DN(5) +  DLR
      DN(6) = DN(6) +  DLR
C
C   EQ. 25J  SEASONAL-LATITUDINAL VARIATION
C            OF HELIUM
C
      DLHE = 0.65*ABS(SUN2/0.4091609)*
     *       (SIN(PID4-0.5*SAT2*SIGN(1.D0,SUN2))**3
     *      -CONS25)
      DN(3)   = DN(3) + DLHE
C
C  COMPUTE DENSITY AND MEAN MOLECULAR WEIGHT
C
      D1 = 10.**DN(1)
      D2 = 10.**DN(2)
      D3 = 10.**DN(3)
      D4 = 10.**DN(4)
      D5 = 10.**DN(5)
      D6 = 10.**DN(6)
      SUMN  = D1+D2+D3+D4+D5+D6
      SUMNM = 28.0134 * D1
     1       +39.9480 * D2
     2       + 4.0026 * D3
     3       +31.9988 * D4
     4       +15.9994 * D5
     5       + 1.00797* D6
      AMW   = SUMNM/SUMN
      DENS  = SUMNM/6.02257D+26
      TEMP(1) = TINF
      TEMP(2) = TZ
C
      RETURN
      END

      SUBROUTINE STJRMO(TINF,SAT3,TZ,DN)
C
C-----
C
C  PURPOSE : STATIC MODEL FOR CALCULATION OF
C            ATMOSPHERIC PROPERTIES AT A GIVEN
C            ALTITUDE.
C
C  INPUTS : TINF...EXOSPHERIC TEMPERATURE (KELVIN)
C           SAT3...ALTITUDE               (KM)
C
C  OUTPUTS : TZ  ... LOCAL TEMPERATURE    (KELVIN)
C            DN  ... LOG10 OF DENSITY NUMBERS (M**-3)
C            DN(1) ... N2
C            DN(2) ... A
C            DN(3) ... HE
C            DN(4) ... O2
C            DN(5) ... O
C            DN(6) ... H
C
C  AUTHOR : HELIO KOITI KUGA - JUNE 1985
C           INPE - DMC/DDO
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DN(6)
C-----
      IF(SAT3.GT.125.) GO TO 20
      IF(SAT3.GT.100.) GO TO 10
      IF(SAT3.LT.90.)  GO TO 30
      CALL STJR01(TINF,SAT3,TZ,DN)
      RETURN
C
   10 CONTINUE
      CALL STJR02(TINF,SAT3,TZ,DN)
      RETURN
C
   20 CONTINUE
      CALL STJR03(TINF,SAT3,TZ,DN)
      RETURN
C
   30 CONTINUE
      PRINT 1000
 1000 FORMAT(1X,'ATENCA0 : MENSAGEM DA ROTINA DE ',
     1     /,1X,'*******   CALCULO DA DENSIDADE AT-',
     2     /,1X,'          MOSFERICA.',
     3    //,3X,'ALTITUDE DO SATELITE MENOR QUE 90 KM')
C
      RETURN
      END


      SUBROUTINE STJR01(TINF,SAT3,TL2,AL10N)
C
C-----
C
C  PURPOSE : STATIC MODEL FOR CALCULATION OF
C            ATMOSPHERIC PROPERTIES FOR THE
C            BAND FROM 90 TO 100 KM.
C
C  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
C           SAT3 ... ALTITUDE               (KM)
C
C  OUTPUTS : TL2 ... LOCAL TEMPERATURE (KELVIN)
C            AL10N . ALOG10 OF DENSITY NUMBERS (M**-3)
C
C            AL10N(1) ... N2
C            AL10N(2) ... A
C            AL10N(3) ... HE
C            AL10N(4) ... O2
C            AL10N(5) ... O
C            AL10N(6) ... H
C
C  AUTHOR : HELIO KOITI KUGA - JUNE 1985
C           INPE - DMC/DDO
C
C  REF. : JACCHIA,L.G.-"ATMOSPHERIC MODELS IN THE
C              REGION FROM 110 TO 2000 KM".IN :
C              CIRA 1972,PART.3,PP. 225-338
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WT(5),AL10N(6)
C
C-----
C
C   RA = POLAR EARTH RADIUS (KM)
C   WT = WEIGHTS FOR THE NEWTON-COTES
C        FIVE POINT QUADRATURE FORMULAE
C
      DATA RA,WT/6356.766,0.31111111,1.4222222,
     *           0.53333333,1.4222222,0.31111111/
C
      TX = 371.668+0.0518806*TINF
     *    -294.3503*EXP(-0.00216222*TINF)
      GX = 0.054285714*(TX-183.)
      AL = LOG(SAT3/90.)
      N  = INT(AL/0.050) + 1
      ZR = EXP(AL/DFLOAT(N))
      AM1 = 28.82678
      TL1 = 183.
      ZEND  = 90.
      SUM2  = 0.
      AIN   = AM1*9.534750028/TL1
C
       DO 2 I=1,N
       Z = ZEND
       ZEND = ZR*Z
       DZ = 0.25*(ZEND-Z)
       SUM1 = 0.31111111*AIN
        DO 1 J=2,5
         Z = Z + DZ
C
C       MOLECULAR WEIGHT FOR Z BETWEEN 90 KM AND
C       100 KM . ACCORDING TO JACCHIA 1971,EQ.1J
C
      ZD  = Z-90.
      AM2 = 28.82678
     1     -7.40066 D-02 * ZD
     2     +ZD * (-1.19407 D-02 * ZD
     3     +ZD * ( 4.51103 D-04 * ZD
     4     +ZD * (-8.21895 D-06 * ZD
     5     +ZD * ( 1.07561 D-05 * ZD
     6            -6.97444 D-07 * ZD * ZD ))))
C
C       TEMPERATURE FOR Z BETWEEN 90 AND 100 KM
C       EQ. 5R
C
      DZX = Z - 125.
      TL2 = TX
     1     +((-9.8204695 D-06*DZX
     2        -7.3039742 D-04)*DZX*DZX
     3        +1.)*DZX*GX
         GZ  = 9.80665*(RA/(RA+Z))**2
         AIN = AM2*GZ/TL2
    1   SUM1 = SUM1 + WT(J)*AIN
    2 SUM2 = SUM2 + DZ*SUM1
C
      FACT1 = 0.12027444181
      DENS  = 3.46D-06*AM2*TL1*
     *        EXP(-FACT1*SUM2)/AM1/TL2
      ANM   = 6.02257D+26*DENS
      AN    = ANM/AM2
      FACT2 = ANM/28.960
C
      AL10N(1) = LOG10(0.78110*FACT2)
      AL10N(2) = LOG10(9.3432D-03*FACT2)
      AL10N(3) = LOG10(6.1471D-06*FACT2)
      AL10N(4) = LOG10(1.20955*FACT2-AN)
      AL10N(5) = LOG10(2.*(AN-FACT2))
      AL10N(6) = AL10N(5)-15.
C
      RETURN
      END


      SUBROUTINE STJR02(TINF,SAT3,TZ,DN)
C
C-----
C
C  PURPOSE : STATIC MODEL FOR CALCULATION OF
C            ATMOSPHERIC PROPERTIES FOR THE
C            BAND FROM 100 TO 125 KM.
C
C  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
C           SAT3 ... ALTITUDE               (KM)
C
C  OUTPUTS : TZ ... LOCAL TEMPERATURE (KELVIN)
C            DN ... ALOG10 OF DENSITY NUMBERS (M**-3)
C
C            DN(1) ... N2
C            DN(2) ... A
C            DN(3) ... HE
C            DN(4) ... O2
C            DN(5) ... O
C            DN(6) ... H
C
C  AUTHOR : HELIO KOITI KUGA - JUNE 1985
C           INPE - DMC/DDO
C
C  REF. : KUGA,H.K. "REFORMULACAO COMPUTACIONAL DO
C              MODELO DE JACCHIA-ROBERTS PARA A
C              DENSIDADE ATMOSFERICA".INPE,SAO JOSE
C              DOS CAMPOS,1985.A SER PUBLICADO
C
C         ROBERTS JR,C.E."AN ANALYTICAL MODEL FOR
C              UPPER ATMOSPHERE DENSITIES BASED UPON
C              JACCHIA'S 1970 MODELS."CELESTIAL
C              MECHANICS 4:368-377,1971.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION DN(6)
C
C-----
C       R ... UNIVERSAL GAS CONSTANT(JOULES/K MOLE)
C       RA... POLAR EARTH RADIUS    (KM)
C       RAS.. RA**2                 (KM**2)
C
      DATA R / 8.31432 D0/
      DATA RA,RAS / 6356.766 D0,4.04084739788 D+07 /
C
C       DENSITY ANALYTICALLY CALCULATED
C
C         EQ. 9J = EQ. 2R
C
      TSUBX = 371.6678 + 0.0518806*TINF
     *       -294.3503*EXP(-0.00216222*TINF)
C
C       EQ. 11J
C
      TXMT0 = TSUBX - 183.
      GSUBX = 0.054285714*TXMT0
C
C       VALUE OF SMALL K <= SK  AND SMALL F <= SF
C
       SKSF  = 9.80665/(R*TXMT0)*1500625.*RAS/0.8
C
C       VALUE OF  C0* <= C0A FOR COMPOSING THE
C       FOURTH DEGREE POLYNOMIAL
C
      C0A   =  -87783750. + 274614375./TXMT0
C
C       NEWTON-RAPHSON PROCEDURE FOR OBTAINING
C       THE TWO REAL ROOTS OF THE QUARTIC
C       POLYNOMIAL C4*P(Z) , EQ. 10R
C
C               INITIAL GUESSES
C
      TEMP = (TSUBX-300.)/200.
      R1   = 167.77 - 3.35*TEMP
      R2   =  57.34 + 7.95*TEMP
C
        DO 10 I=1,7
         PZ1 = C0A + 3542400.*R1 +
     *         R1*(R1*(-52687.5+340.5*R1-0.8*R1*R1))
         PZ2 = C0A + 3542400.*R2 +
     *         R2*(R2*(-52687.5+340.5*R2-0.8*R2*R2))
         DPZ1 = 3542400.-105375.*R1+
     *          R1*(1021.5*R1-3.2*R1*R1)
         DPZ2 = 3542400.-105375.*R2+
     *          R2*(1021.5*R2-3.2*R2*R2)
         R1N = R1 - PZ1/DPZ1
         R2N = R2 - PZ2/DPZ2
      IF(ABS(R1N-R1).LT.1.D-07.AND.
     *   ABS(R2N-R2).LT.1.D-07     ) GO TO 20
         R1 = R1N
         R2 = R2N
   10   CONTINUE
   20   CONTINUE
      R1 = R1N
      R2 = R2N
C
C       COMPLEX ROOTS OR X & X**2+Y**2
C
        SOMA = R1 + R2
        PROD = R1 * R2
        DIFE = R1 - R2
        X      = -0.5*( SOMA-425.625 )
        X2Y2   = -C0A/(0.8*PROD)
C
C       CALCULATE U(R1),U(R2),W(R1),W(R2),CX(CAPITAL X),
C                 AND V(-RA)
C
C  EXPRESSION OF W CORRECTED ACCORDING TO GSFC (NASA,1976)
C
      H2 = R1 + RA
      H3 = R2 + RA
      H4 = RAS + 2.*X*RA + X2Y2
C
      UR1H2 = H2*(R1*R1-2.*X*R1+X2Y2)*DIFE
      UR2H3 = H3*(R2*R2-2.*X*R2+X2Y2)*DIFE
      WR1 = RA+X2Y2/R1
      WR2 = RA+X2Y2/R2
      VRA = H4 * H2 * H3
      CX  = -H4 - H4
C
      DE100 = (((((0.7026942D-32*TINF*TINF
     2            -0.7734110D-28*TINF)*TINF
     3            +0.3727894D-24*TINF)*TINF
     4            -0.1021474D-20*TINF)*TINF
     5            +0.1711735D-17*TINF)*TINF
     6            -0.1833490D-14*TINF
     7            +0.1985549D-10 )
      T100   = TSUBX - 0.94585589*TXMT0
      Z      = SAT3
C
C      NUMBER OF PARTICLES PER M**3 AT 100 KM
C
C         D  ... TOTAL NUMBER
C         D1 ... N2 NITROGEN
C         D2 ... AR ARGON
C         D3 ... HE HELIUM
C         D4 ... O2 DIATOMYC OXYGEN
C         D5 ... O  MONOATOMYC OXYGEN
C
      AM100 = 27.6396281382
      DEAVOG = DE100*6.02257D+29
      D1 = 0.78110    * DEAVOG
      D2 = 0.0093432  * DEAVOG
      D3 = 6.1471D-06 * DEAVOG
      D4 = (1.20955-28.96/AM100) *DEAVOG
      D5 = 2.*(28.96-AM100)/AM100*DEAVOG
C
C      Q(I) PARAMETERS
C
      UR1 = H2*UR1H2
      UR2 = H3*UR2H3
      Q2 = 1/UR1
      Q3 =-1/UR2
      Q5 = 1/VRA
      Q4 =( 1./(PROD*RA)
     1     +(RA-X2Y2/RA)/VRA
     2     +WR1/UR1H2
     3     -WR2/UR2H3)/CX
      Q6 = -Q5
     1     -2.*(X+RA)*Q4
     2     +1./UR2H3
     3     -1./UR1H2
      Q1 =-Q4-Q4-Q3-Q2
C
C       TEMPERATURE FOR Z BETWEEN 100 AND 125 KM
C       EQ. 5R
C
      DZX = Z - 125.
      TZ  = TSUBX
     1     +((-9.8204695 D-06*DZX
     2        -7.3039742 D-04)*DZX*DZX
     3        +1.)*DZX*GSUBX
C
      AUX = Z - 100.
      Y   = SQRT(X2Y2-X*X)
      AUX1 = Z + RA
      AUX2 = RA + 100.
C
      F3 = DLOG(AUX1/AUX2)*Q1
     1    +LOG( (Z-R1)/(100.-R1) )*Q2
     2    +LOG( (Z-R2)/(100.-R2) )*Q3
     3    +LOG( (Z*Z-2.*X*Z+X2Y2)
     4      /(10000.-200.*X+X2Y2) )*Q4
C
      F4 = Q5*AUX/(AUX1*AUX2)
     1    +Q6/Y*ATAN(Y*AUX/(X2Y2+100.*Z-(100.+Z)*X))
C
C      DENSITY NUMBERS D(I) : N2,AR,HE,O2,O,H
C      EQ. 20R
C
      T100TZ = T100/TZ
      SKSF34 = SKSF*(F3+F4)
C
      DN(1) = LOG10(D1*T100TZ*EXP(28.0134*SKSF34))
      DN(2) = LOG10(D2*T100TZ*EXP(39.9480*SKSF34))
      DN(3) = LOG10(D3*T100TZ**0.62*EXP(4.0026*SKSF34))
      DN(4) = LOG10(D4*T100TZ*EXP(31.9988*SKSF34))
      DN(5) = LOG10(D5*T100TZ*EXP(15.9994*SKSF34))
C
      RETURN
      END


      SUBROUTINE STJR03(TINF,SAT3,TZ,DN)
C
C-----
C
C  PURPOSE : STATIC MODEL FOR CALCULATION OF
C            ATMOSPHERIC PROPERTIES FOR THE
C            BAND ABOVE 125 KM.
C
C  INPUTS : TINF ... EXOSPHERIC TEMPERATURE (KELVIN)
C           SAT3 ... ALTITUDE               (KM)
C
C  OUTPUTS : TZ ... LOCAL TEMPERATURE (KELVIN)
C            DN ... ALOG10 OF DENSITY NUMBERS (M**-3)
C
C            DN(1) ... N2
C            DN(2) ... A
C            DN(3) ... HE
C            DN(4) ... O2
C            DN(5) ... O
C            DN(6) ... H
C
C  AUTHOR : HELIO KOITI KUGA - JUNE 1985
C           INPE - DMC/DDO
C
C  REFS. KUGA,H.K."REFORMULACAO COMPUTACIONAL DO
C              MODELO ATMOSFERICO DE JACCHIA-
C              ROBERTS PARA A DENSIDADE ATMOSFERICA".
C              INPE,SAO JOSE DOS CAMPOS,1985.A SER
C              PUBLICADO.
C
C        ROBERTS JR,C.E."AN ANALYTIC MODEL FOR UPPER
C              ATMOSPHERE DENSITIES BASED UPON
C              JACCHIA'S 1970 MODELS."CELESTIAL
C              MECHANICS 4:368-377,1971
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/NCALL/ICOUNT
      DIMENSION DN(6)
C
C-----
C
C       R....UNIVERSAL GAS CONSTANT (JOULES/K MOLE)
C       RA...POLAR EARTH RADIUS     (KM)
C       RAS..RA**2                  (KM**2)
C
      DATA R / 8.31432 D0 /
      DATA RA,RAS / 6356.766 D0 , 4.04084739788 D+07 /
      ICOUNT = ICOUNT+1
C
C       DENSITY ANALYTICALLY CALCULATED
C
C         EQ. 9J = EQ. 2R
C
      TSUBX = 371.6678 + 0.0518806*TINF
     *       -294.3503*EXP(-0.00216222*TINF)
C
      D1 = ((((-0.2296182D-19*TINF*TINF
     1         +0.1969715D-15*TINF)*TINF
     2         -0.7139785D-12*TINF)*TINF
     3         +0.1420228D-08*TINF)*TINF
     4         -0.1677341D-05*TINF)*TINF
     5         +0.1186783D-02*TINF
     6         +0.1093155D+02
      D2 = ((((-0.4837461D-19*TINF*TINF
     1         +0.4127600D-15*TINF)*TINF
     2         -0.1481702D-11*TINF)*TINF
     3         +0.2909714D-08*TINF)*TINF
     4         -0.3391366D-05*TINF)*TINF
     5         +0.2382822D-02*TINF
     6         +0.8049405D+01
      D3 = (((-0.1270838D-16*TINF*TINF
     1        +0.9451989D-13*TINF)*TINF
     2        -0.2894886D-09*TINF)*TINF
     3        +0.4694319D-06*TINF)*TINF
     4        -0.4383486D-03*TINF
     5        +0.7646886D+01
      D4 = ((((-0.3131808D-19*TINF*TINF
     1         +0.2698450D-15*TINF)*TINF
     2         -0.9782183D-12*TINF)*TINF
     3         +0.1938454D-08*TINF)*TINF
     4         -0.2274761D-05*TINF)*TINF
     5         +0.1600311D-02*TINF
     6         +0.9924237D+01
      D5 = ((( 0.5116298D-17*TINF*TINF
     1        -0.3490739D-13*TINF)*TINF
     2        +0.9239354D-10*TINF)*TINF
     3        -0.1165003D-06*TINF)*TINF
     4        +0.6118742D-04*TINF
     5        +0.1097083D+02
      Z   = SAT3
      AL = (( 0.2462708D-09*TINF*TINF
     1       -0.1252487D-05*TINF)*TINF
     2       +0.1579202D-02*TINF)*TINF
     3       +0.2341230D+01*TINF
     4       +0.1031445D+05
C
C      TEMPERATURE PROFILE EQ. 23R
C
      TXMT0 = TSUBX - 183.
      TETX = TINF - TSUBX
      TZ   = TETX*EXP(-TXMT0/TETX*(Z-125.)/35.*AL/(Z+RA))
C
C      PARAMETERS G(I) : N2,AR,HE,O2,O EQ. 25'R
C
      AUX = 9.80665*RAS/(R*AL*TINF)
     *     *TETX/TXMT0*35./6481.766
C
      AUX1 = TSUBX/(TINF - TZ)
      AUX2 = TZ / TETX
      A1   = LOG10(AUX1)
      A2   = LOG10(AUX2)
      A1A2A= (A1+A2)*AUX
C
C        DENSITY NUMBERS D(I) EQ.25R
C
      D1 = D1+     A1+28.0134*A1A2A
      D2 = D2+     A1+39.9480*A1A2A
      D3 = D3+0.62*A1+ 4.0026*A1A2A
      D4 = D4+     A1+31.9988*A1A2A
      D5 = D5+     A1+15.9994*A1A2A
C
C       CALCULATE TZ(500) EQ.23R
C       T(Z) ALREADY CALCULATED
C
      TZ500 = TINF-TETX*EXP(
     *       -TXMT0/TETX*(375./35.*AL/(500.+RA)))
C
C       INCLUSION OF HYDROGEN
C
C       DENSITY NUMBER FROM EQS. 26R,27R
C
      AUX1 = LOG10(TINF)
      H500 = 73.13-(39.4-5.5*AUX1)*AUX1
      A1   = LOG10(TZ500/(TINF - TZ))
      A2   = LOG10(TZ/(TINF-TZ500))
      D6   = H500+A1+1.00797*AUX*(A1+A2)
C
C        LOAD ALOG10 OF DENSITY NUMBERS
C        IN M**-3
C
      DN(1) = D1 + 6.
      DN(2) = D2 + 6.
      DN(3) = D3 + 6.
      DN(4) = D4 + 6.
      DN(5) = D5 + 6.
      DN(6) = D6 + 6.
	TZ    = TINF - TZ
C
      RETURN
      END
