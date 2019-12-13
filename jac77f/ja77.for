      SUBROUTINE IDYMOS(SA,SU,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C------
C
C
C PURPOSE:
C
C     THE SUBROUTINE IDYMOS GIVES THE  TEMPERATURE,  DENSITY
C
C     AND MOLECULAR WEIGHT OF ATMOSPHERE, USING THE INTEGRA-
C                                                   -
C     TED STATIC AND DYNAMIC  MODEL,  WITH  THE  SOLAR  FLUX
C                    --       --                 -
C     DATA FILE.
C
C INPUTS:
C
C     SA(1) RIGHT ASCENSION OF THE  POINT  IN  QUESTION,  IN
C           RADIANS (0 TO 2.*PI)
C     SA(2) DECLINATION (GEOCENTRIC LATITUDE) OF THE  POINT,
C           IN RADIANS (-PI TO PI).
C     SA(3) GEOCENTRIC ALTITUDE OF THE POINT IN  METERS,  IN
C           THE RANGE 90000. TO 2000000.M.
C     SU(1) RIGHT ASCENSION OF THE SUN AT THE DATE, IN RADI-
C           ANS (O TO 2.*PI).
C     SU(2) SUN DECLINATION IN RADIANS (-PI TO PI)
C     RJUD  MODIFIED JULIAN DATE (IF OUT OF RANGE, THE  SUB-
C           ROUTINE WILL PRINT A MESSAGE AND STOP).
C     DAFR  TIME (UT) OF THE DAY, IN SECONDS.
C     GSTI  GREENWICH SIDERAL TIME, IN RADIANS (0 TO 2.*PI),
C           AT THE TIME DAFR OF THE DATE RJUD.
C
C OUTPUTS:
C
C     TE(1) MEAN EXOSPHERIC TEMPERATURE ABOVE THE  POINT  IN 
C           QUESTION, AS DEFINED BY EQUATION  20,  IN  REFE-
C           RENCE (1), IN KELVIN.
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
C     ISDAMO
C
C REFERENCES:
C
C     (1) JACCHIA, L. G. "THERMOSPHERIC TEMPERATURE, DENSITY
C         AND COMPOSTION: NEW MODELS." CAMBRIDGE,  MA,  SAO,
C         1977. (SAO SPECIAL REPORT NO 375).
C
C AUTHORS:
C
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     JAN. 1985              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SA(3),SU(2),TE(2),AD(6)
      DIMENSION SF(3),SD(15)
C
C------
C

      CALL SOFLUD(RJUD,DAFR,SD,OUTR)

      IF(OUTR.NE.0.) THEN
       WRITE(6,*) ' ERROR IN ROUTINE IDYMOS: OUTR = ',
     1            INT(OUTR)
       STOP
      ENDIF

      RLAT    = SA(2)
      RLON    = SA(1) - GSTI
      SIFI    = 0.9792*SIN(RLAT) + 
     1          0.2028*COS(RLAT)*COS(RLON - 5.07891)
      COFI    = SQRT(1 - SIFI*SIFI)
      TAUO    = 2.4 + 4.8*COFI*COFI

      ND      = (DAFR/3600. - SD(6) + 12. - TAUO)/3.
      SF(1)   = SD(9)
      SF(2)   = SD(10)
      SF(3)   = SD(ND)

      CALL ISDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,RHOD)

      RETURN
      END

      SUBROUTINE ISMODS(ALTU,RJUD,DAFR,TE,AL,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE  SUBROUTINE  ISMODS  USES  THE  INTEGRATED  STATIC 
C                                         -           -
C     MODEL (1) AND  THE  SOLAR  FLUX  DATA  TO  OBTAIN  THE 
C     ---                 -
C     TEMPERATURE, MOLECULAR WEIGHT  AND  DENSITY  OF  LOCAL
C
C     ATMOSPHERE.
C
C INPUTS:
C
C     ALTU  ALTITUDE  OF  THE  POINT  IN  METERS  (90000  TO
C           2000000).
C     RJUD  MODIFIED  JULIAN  DATE  (IF  OUT  OF  RANGE  THE 
C           ROUTINE WILL PRINT A MESSAGE AND STOP).
C     DAFR  TIME (UT) OF THE DAY, IN SECONDS (0  TO  86400).
C
C OUTPUTS:
C
C     TE(1) MEAN EXOSPHERIC TEMPERATURE ABOVE THE  POINT  IN
C           QUESTION, AS DEFINED BY EQUATION  20,  IN  REFE-
C           RENCE (1), IN KELVIN.
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
C     ISMADE
C
C REFERENCES:
C
C   (1) JACCHIA, L. G. "THERMOSPHERIC TEMPERATURE,   DENSITY
C       AND COMPOSTION: NEW MODELS." CAMBRIDGE,   MA,   SAO,
C       1977. (SAO SPECIAL REPORT NO 375).
C
C AUTHORS:
C
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     JAN. 1985              V. 1.0
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
      SF(2)   = SD(10)
      VARI    = .154*SD(7)
      SF(3)   = 1.89*DLOG(VARI + DSQRT(VARI*VARI + 1.D0))

      CALL ISMADE(ALTU,SF,TE,AL,WMOL,RHOD)

      RETURN
      END


      SUBROUTINE ISDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  ISDAMO  GIVES  THE  DENSITY, MOLECULAR
C
C     WEIGHT AND TEMPERATURE OF THE UPPER ATMOSPHERE,  USING
C
C     THE INTEGRATED STATIC VERSION (JACCHIA 1977)  AND  THE
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
C                RANGE 90,000-2,000,000.
C       SU(1)    RIGHT ASCENTION OF THE SUN AT THE DATE,  IN
C                RADIANS (0 TO 2*PI).
C       SU(2)    SUN DECLINATION IN RADIANS (-PI TO PI).
C       SF(1)    DAILY SOLAR FLUX AT  10.7 CM,  ADJUSTED FOR
C                THE DISTANCE BETWEEN THE EARTH AND THE  SUN
C                AT THE DATE, IN 1E-22 W/M/M/HZ.
C       SF(2)    AVERAGED DAILY FLUX  AS DEFINED BY JACCHIA,
C                IN 1E-22 W/M/M/HZ, ADJUSTED FOR THE  EARTH-
C                SUN DISTANCE.
C       SF(3)    3-HOURLY PLANETARY GEOMAGNETIC INDEX KP, AT
C                THE TIME DAFR - TAU, WHERE TAU IS GIVEN  BY
C                JACCHIA'S 1977 MODEL.
C       RJUD     MODIFIED JULIAN  DATE,  REFERED  TO  1950.0
C                (JULIAN DATE-2433282.5).
C       DAFR     TIME (UT) OF THE DAY, IN SECONDS.
C       GSTI     GREENWICH SIDERAL TIME, IN RADIANS, AT  THE
C                TIME DAFR OF THE DATE RJUD (0 TO 2*PI).
C
C OUTPUTS:
C
C       TE(1)    MEAN EXOSPHERIC TEMPERATURE ABOVE THE POINT
C                AS DEFINED BY EQUATION 20 IN THE  JACCHIA'S
C                1977 MODEL, IN KELVIN.
C       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
C                KELVIN
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
C       ISMADE
C       DIVARI
C       GEOACI
C       SEALAT
C       SEMIAN
C
C REFERENCES:
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     ABR. 1987              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION SA(3),SU(2),SF(3),TE(2),AD(6)
      DIMENSION AL(6),WM(6)
  
      COMMON /STADE/ AC(6),TH(6)
C
C------
C  
      DATA AVOG /6.02217D+26/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/
  
      YEFR   = (RJUD + DAFR/.86400D05)/.3652422D03
      TYFR   = YEFR - INT(YEFR)
      GEAC   = SF(3)
      RLAT   = SA(2)
      RLON   = SA(1) - GSTI
      HEIG   = SA(3)
      ALTU   = HEIG/1.D03
      SUDC   = SU(2)
  
      CALL ISMADE(HEIG,SF,TE,AL,WMOL,RHOD)
  
      TEXO   = TE(1)
  
      CALL DIVARI(TEXO,SA,SU,WMOL,AL)
  
      DO I   = 1, 6
       AD(I) = AL(I)
      ENDDO
  
      TQUT   = TH(6)
  
      CALL GEOACI(TQUT,GEAC,RLAT,RLON,ALTU,AL)
  
      DO I   = 1, 6
       AD(I) = AD(I) + AL(I) - AC(I)
      ENDDO
  
      CALL SEALAT(TYFR,SUDC,RLAT,ALTU,AL)
  
      DO I   = 1, 6
       AD(I) = AD(I) + AL(I)
      ENDDO
  
      CALL SEMIAN(TYFR,ALTU,ALCO)
  
      DO I   = 1, 6
       AD(I) = AD(I) + ALCO
      ENDDO
  
      WEIG   = 0.D0
      ANUT   = 0.D0
  
      DO   I = 1, 6
       ANAC  = 10.D0**AD(I)
       ANUT  = ANUT + ANAC
       WEIG  = WEIG + WM(I)*ANAC
      ENDDO
  
      WMOL   = WEIG/ANUT
      RHOD   = WEIG/AVOG
  
      RETURN
      END
  
  
      SUBROUTINE ISMADE(ALTU,SF,TE,AL,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE CALCULATES THE ATMOSPHERIC DENSITY  FOR
C
C     HEIGHTS FROM 100 TO 2000 KM, USING THE INTEGRATED STA-
C                                            -          -
C     TIC MODEL TO COMPUTE THE ATMOSPHERIC DENSITY (1).
C         -                    -           --
C
C INPUTS:
C
C       ALTU     ALTITUDE OF THE POINT IN METERS.
C       SF(1)    DAILY OBSERVED  SOLAR  FLUX AT 10.7 CM,  IN
C                1E-22 W/M/M/HZ.
C       SF(2)    AVERAGED DAILY FLUX, AS DEFINED BY JACCHIA,
C                IN 1E-22 W/M/M/HZ.
C
C OUTPUTS:
C
C       TE(1)    MEAN EXOSPHERIC TEMPERATURE ABOVE THE POINT
C                AS DEFINED BY EQUATION 20 IN THE  JACCHIA'S
C                1977 MODEL, IN KELVIN.
C       TE(2)    LOCAL  TEMPERATURE  AROUND  THE  POINT,  IN
C                KELVIN
C       AL(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AL(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AL(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AL(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AL(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AL(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
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
C       IMOWEI
C       TEMLO
C
C REFERENCES:
C
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     ABR. 1987              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION SF(3),TE(2),AL(6)
      DIMENSION WM(6),CP(7)
C
C------
C  
      DATA AVOG /6.02217D+26/
      DATA PI   /3.14159265359D0/
  
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/
      DATA TO   /188.D0/
      DATA ZX   /125.D0/
      DATA ZO   /90.0D0/
      DATA BE   /5.5D-5/
  
      FLUX   = SF(1)
      FBAR   = SF(2)
      THAF   = 5.48D0*(FBAR**.8D0) + 101.8D0*(FLUX**.4D0)
      HEIG   = ALTU/1.D3
  
      CALL IMOWEI(THAF,HEIG,AL,WMOL,RHOD)
  
      POTE   = 0.0045D0*(THAF - TO)
      ARGU   = POTE/DSQRT(1.D0 + POTE*POTE)
      TX     = TO + 110.5D0*DATANH(ARGU)
      GX     = 1.9D0*(TX - TO)/(ZX - ZO)
  
      CP(1)  = 2.D0*(TX - TO)/PI
      CP(2)  = GX/CP(1)
      CP(3)  = 1.7D0*CP(2)
      CP(4)  = 2.D0*(THAF - TX)/PI
      CP(5)  = GX/CP(4)
      CP(6)  = BE*CP(5)
      CP(7)  = TX
  
      TE(1)  = THAF
      TE(2)  = TEMLO(HEIG,CP)
  
      RETURN
      END
  
  
      SUBROUTINE DIVARI(TEXO,SA,SU,WMOL,AL)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  CALCULATES THE  DIURNAL VARIATIONS  IN
C                                     --      ---
C     THE ATMOSPHERE DENSITY, AS PROPOSED BY JACCHIA,  USING
C
C     THE INTEGRATED STATIC MODEL OF JACCHIA (1).
C         -
C
C INPUTS:
C
C       TEXO     MEAN EXOSPHERIC TEMPERATURE, IN KELVIN.
C       SA(1)    RIGHT ASCENTION OF THE POINT, IN RADIANS.
C       SA(2)    DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
C                POINT, IN RADIANS (-PI TO PI).
C       SA(3)    GEOCENTRIC ALTITUDE IN METERS, BETWEEN  THE
C                RANGE 90,000-2,000,000.
C       SU(1)    RIGHT ASCENTION OF THE SUN AT THE DATE,  IN
C                RADIANS (0 TO 2*PI).
C       SU(2)    SUN DECLINATION IN RADIANS (-PI TO PI).
C       WMOL     MEAN MOLECULAR WEIGHT AT THE LOCAL  REGION,
C                IN KG/KGMOL.
C
C OUTPUTS:
C
C       AL(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AL(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AL(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AL(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AL(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AL(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
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
C COMMON OUTPUTS:
C       COMMON /STADE/ AC(6),TH(6)
C                AC     SAME AS  AL,  BUT  USING THE  PSEUDO
C                       EXOSPHERIC   TEMPERATURE   FOR   THE
C                       HYDROGEN.
C                TH(1)  EXOSPHERIC TEMPERATURE FOR HE, IN K.
C                TH(2)  EXOSPHERIC TEMPERATURE FOR O2, IN K.
C                TH(3)  EXOSPHERIC TEMPERATURE FOR N2, IN K.
C                TH(4)  EXOSPHERIC TEMPERATURE FOR AR, IN K.
C                TH(5)  EXOSPHERIC TEMPERATURE FOR  O, IN K.
C                TH(6)  EXOSPHERIC TEMPERATURE FOR  H, IN K.
C
C SUBCALLS:
C
C       IMOWEI
C
C REFERENCES:
C
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     ABR. 1987              V. 1.0
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION SA(3),SU(2),AL(6)
      DIMENSION AN(6),WM(6),CP(7)
  
      COMMON /STADE/ AC(6),TH(6)
C
C------
C  
      DATA AVOG /6.02217D+26/
      DATA PI   /3.14159265359D0/
      DATA RAD  /1.7453292778D-2/
  
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1           1.00797/
  
      RLAT   = SA(2)
      DELT   = SU(2)
      SOAN   = SA(1) - SU(1)
      ALTU   = SA(3)/1.D3
      ARGU   = RLAT*RLAT/1.5707963D0
      COAR   = DCOS(ARGU)
      COLA   = DCOS(RLAT)
      E1     = 2.D0 + COAR*COAR
      TEAU   = 1.D0 + 0.3666069D0*DELT*DSIN(RLAT)
  
      DO   I = 1, 5
       BISH  = RAD*(27.D0*(WMOL/WM(I) - 1.D0) - 35.D0)
       ARGU  = SOAN + BISH
       FHAG  = 0.08D0*DCOS(3.D0*ARGU - 1.3089969D0) +
     1         DABS(DCOS(0.5D0*ARGU))**E1
       TH(I) = (0.24D0*COLA*(FHAG - 0.5D0) + TEAU)*TEXO
      ENDDO
  
      ARGU   = SOAN - 1.0471976D0
      FHAG   = 0.08D0*DCOS(3.D0*ARGU - 1.3089969D0) +
     1         DABS(DCOS(0.5D0*ARGU))**E1
      TH(6)  = TEXO*(TEAU + 0.24*COLA*(FHAG - 0.5D0))
  
      DO   I = 1, 6
       TEMP   = TH(I)
       CALL IMOWEI(TEMP,ALTU,AN,WEIG,RHOD)
       AL(I)  = AN(I)
      ENDDO
  
      DO   I = 1, 6
       AC(I)  = AN(I)
      ENDDO
  
      RETURN
      END
  
  
      SUBROUTINE GEOACI(TQUT,GEAC,RLAT,RLON,ALTU,DN)
C
C------
C
C PURPOSE:
C
C     THE  SUBROUTINE  GEOACA  OBTAINS  THE VARIATION OF THE
C
C     ATMOSPHERE  NUMBER-DENSITY,  DUE  TO  THE  GEOMAGNETIC
C                                                ---
C     ACTIVITY, USING THE INTEGRATED STATIC MODEL.
C     --                  -
C
C INPUTS:
C
C       TQUT     EXOSPHERIC  TEMPERATURE  FOR  THE  HYDROGEN
C                ABOVE THE POINT, IN KELVIN.
C       GEAC     THE PLANETARY THREE-HOUR INDEX, KP, IN  THE
C                RANGE 0. TO 9.0
C       RLAT     DECLINATION (GEOCENTRIC  LATITUDE)  OF  THE
C                POINT, IN RADIANS (-PI TO PI).
C       RLON     LOCAL LONGITUDE, IN RADIANS (-PI TO PI).
C       ALTU     GEOCENTRIC  ALTITUDE  IN  KM,  IN THE RANGE
C                90. TO 2000..
C
C OUTPUTS:
C
C       DN       ARRAY WITH THE COMMON LOGARITHM  (BASE 10),
C                OF THE VARIATION IN THE  NUMBER DENSITY  OF
C                THE HE (HELIUM), O2 (MOLECULAR OXYGEN),  N2
C                (MOLECULAR NITROGEN), AR (ARGON), O (ATOMIC
C                OXYGEN) AND H (ATOMIC HYDROGEN),  RESPECTI-
C                VELY, DUE TO THE GEOMAGNETIC ACTIVITY.
C
C SUBCALLS:
C
C       IMOWEI
C
C REFERENCES:
C
C       [1]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C AUTHORS:
C
C     VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C     BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C
C DATE:
C
C     ABR. 1987              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION DN(6)
      DIMENSION GM(6)
C
C------
C  
      DATA GM   /-6.30D-05,1.03D-05,0.D0,3.07D-05,-4.85D-05,
     1          0.D0/
  
      AUXI   = 57.5D0*GEAC*(1.D0 + 0.027D0*DEXP(0.4D0*GEAC))
      SILA   = DSIN(RLAT)
      SENO   = 0.9792*SILA +
     1         0.2028*DCOS(RLAT)*DCOS(RLON - 5.0789081D0)
      SEN2   = SENO*SENO
      COL2   = 1.D0 - SEN2
      DELT   = AUXI*SEN2*SEN2
      TEMP   = TQUT + DELT
  
      CALL IMOWEI(TEMP,ALTU,DN,WMOL,RHOD)
  
      DELT   = 0.01D0*DELT
      DEZH   = 5000.D0*DLOG(DSQRT(1.D0 + DELT*DELT) + DELT)
      DELE   = 5.2D-04*AUXI*COL2*COL2
  
      DO   I = 1, 6
       DN(I) = DN(I) + GM(I)*DEZH + DELE
      ENDDO
  
      RETURN
      END
  
  
      SUBROUTINE IMOWEI(TEXO,ALTU,AN,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THIS SUBROUTINE USES THE  NUMERICALLY  INTEGRATED  MO-
C                                            -
C     DEL (JACCHIA, 1977-1) FOR THE MOLECULAR WEIGHT OF  THE 
C                                   --        ---
C     UPPER  ATMOSPHERE.
C
C INPUTS:
C
C     TEXO  EXOSPHERIC TEMPERATURE, AS DESCRIBED IN  JACCHIA
C           (1), IN KELVIN.
C     ALTU  GEOMETRICAL HEIGHT WHERE  THE  MOLECULAR  WEIGHT
C           WILL BE CALCULATED, IN KM. (90 TO 2000KM.).
C
C OUTPUTS:
C
C     AN(1:6) COMMON LOGARITHM (BASE 10) OF THE HE (HELIUM),
C             02 (MOLECULAR OXYGEN), N2 (MOLECULAR    NITRO-
C             GEN), AR (ARGON), O (ATOMIC OXYGEN) AND H  (A-
C             TOMIC HIDROGEN) NUMBER-DENSITY, RESPECTIVELY.
C
C     WMOL  MEAN MOLECULAR WEIGHT, IN KG/KGMOL.
C
C     RHOD  ATMOSPHERIC DENSITY IN KG/M/M/M.
C
C REFERENCES:
C
C    (1)  JACCHIA, L. G. "THERMOSPHERIC TEMPERATURE, DENSITY
C         AND COMPOSITION: NEW MODELS." CAMBRIDGE,  MA  SAO,
C         1977. (SAO SPECIAL REPORT NO 375).
C
C
C AUTHORS:
C
C         BENTO SILVA DE MATOS   - INPE - S.J.CAMPOS - BR
C         VALDEMIR CARRARA       - INPE - S.J.CAMPOS - BR
C
C DATE:
C         SET. 1984              V. 1.0
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AN(6)
      DIMENSION WT(5),AP(6),WM(6),BN(6),AA(7),AL(6),CI(6)
C
C------
C

      DATA WT   /0.31111111111,1.4222222222,0.53333333333,
     1           1.4222222222,0.31111111111/
      DATA CI   /28.89122,-2.83071E-02,-6.59924E-03,
     1           -3.39574E-04,6.19256E-05,-1.84796E-06/

      DATA CONV /2.302585092994/
      DATA PI   /3.14159265359/

      DATA AVOG /6.02217E+26/
      DATA RGAS /8.31432  /
      DATA AP   /-0.38,0.0,0.0,0.0,0.0,-0.25/
      DATA WM   /4.0026,31.9988,28.0134,39.948,15.9994,
     1          1.00797/

      DATA FRC1 /5.242E-06/
      DATA FRC2 /0.20955  /
      DATA FRC3 /0.78110  /
      DATA FRC4 /9.343E-03/
      DATA RT   /6356.766 /
      DATA GO   /9.80665  /
      DATA RHOI /3.43E-06 /
      DATA RMLI /28.960   /
      DATA BE   /5.5E-05  /
      DATA ZO   /90.0     /
      DATA ZX   /125.0    /
      DATA TO   /188.0    /

      POTE    = 0.0045*(TEXO - TO)
      AUXI    = POTE/SQRT(1. + POTE*POTE)
      TX      = TO + 110.5*DATANH(AUXI)
      AA(1)   = 2.*(TX - TO)/PI
      AA(4)   = 2.*(TEXO - TX)/PI
      GX      = 1.9*(TX - TO)/(ZX - ZO)
      AA(2)   = GX/AA(1)
      AA(3)   = 1.7*AA(2)
      AA(5)   = GX/AA(4)
      AUXI    = 28.9/TEXO**0.25
      PHIS    = 10.**(6.90 + AUXI)/2.0E+20
      HND5    = (5.94 + AUXI)*CONV
      AA(6)   = BE*AA(5)
      AA(7)   = TX
      NC      = 5

      STEP    = 0.05
      RSLT    = 0.
      ZINI    = ZO
      ZEND    = 100.

      IF(ALTU.LT.ZEND) ZEND  = ALTU

      ZAUX    = (ZEND - ZINI)/4.
      AUXI    = INT(ZAUX/STEP)
      IF(AUXI.LE.0.) AUXI  = 1.
      STEP    = ZAUX/AUXI

      DO WHILE (ABS(ZEND-ZINI).GT.0.0001)
       SUMM   = 0.
       ZAUX   = ZINI + 2.*STEP
       AUXI   = 1. + ZAUX/RT
       GRAV   = GO/AUXI/AUXI
       ZAUX   = ZINI

       DO   I = 1, 5
        H     = ZAUX - ZO
        WEIG  = ((((CI(6)*H + CI(5))*H + CI(4))*H + 
     1          CI(3))*H + CI(2))*H + CI(1)
        SUMM  = SUMM + GRAV*WEIG/TEMLO(ZAUX,AA)*WT(I)
        ZAUX  = ZAUX + STEP
       ENDDO

       RSLT   = RSLT + STEP*SUMM
       ZINI   = ZINI + 4.*STEP
      ENDDO

      RHOS    = RHOI*EXP(-RSLT/RGAS)
      H       = ZEND - ZO
      WEIG    = ((((CI(6)*H + CI(5))*H + CI(4))*H + 
     1          CI(3))*H + CI(2))*H + CI(1)
      TEMF    = TEMLO(ZEND,AA)
      DENU    = AVOG*RHOS/RMLI*TO/TEMF
      RHOS    = DENU*WEIG
      FAT2    = RHOS/RMLI
      AN(1)   = LOG(FRC1*FAT2)
      AN(2)   = LOG(FAT2*(1. + FRC2) - DENU)
      AN(3)   = LOG(FRC3*FAT2)
      AN(4)   = LOG(FRC4*FAT2)
      AN(5)   = LOG(2.*(DENU - FAT2))
      AN(6)   = 0.

      IF(ALTU.LE.100.) GOTO 800

      STEP    = 0.05
      ZINI    = 100.
      TEMI    = TEMF
      RSLT    = 0.
      ZEND    = 140.
      IF(ALTU.LT.ZEND) ZEND  = ALTU

      ZAUX    = (ZEND - ZINI)/4.
      AUXI    = INT(ZAUX/STEP)
      IF(AUXI.LE.0.) AUXI  = 1.
      STEP    = ZAUX/AUXI

      DO WHILE (ABS(ZEND-ZINI).GT.0.0001)
       ZAUX   = ZINI + 2.*STEP
       AUXI   = 1. + ZAUX/RT
       GRAV   = GO/AUXI/AUXI/RGAS
       SUMM   = 0.
       ZAUX   = ZINI

       DO   I = 1, 5
        SUMM  = SUMM + WT(I)*GRAV/TEMLO(ZAUX,AA)
        ZAUX  = ZAUX + STEP
       ENDDO

       RSLT   = RSLT + STEP*SUMM
       ZINI   = ZINI + 4.*STEP
      ENDDO

      TEMF    = TEMLO(ZEND,AA)
      AUXI    = LOG(TEMI/TEMF)

      DO    I = 1, 5
       AN(I)  = AN(I) - RSLT*WM(I) + AUXI*(1 + AP(I))
      ENDDO

      IF(ALTU.LE.140.) GOTO 800

      NC      = 6
      ZINI    = 140.
      STEP    = 5.
      RSLT    = 0.
      TEMI    = TEMF
      ZEND    = 500.

      DO WHILE (ABS(ZEND-ZINI).GT.0.0001)
       ZAUX   = ZINI + 2.*STEP
       AUXI   = 1. + ZAUX/RT
       GRAV   = GO/AUXI/AUXI/RGAS
       SUMM   = 0.
       ZAUX   = ZINI

       DO   I = 1, 5
        SUMM  = SUMM + GRAV*WT(I)/TEMLO(ZAUX,AA)
        ZAUX  = ZAUX + STEP
       ENDDO

       RSLT   = RSLT + STEP*SUMM
       ZINI   = ZINI + 4.*STEP
      ENDDO

      TEMF    = TEMLO(ZEND,AA)
      AUXI    = LOG(TEMI/TEMF)
      AN(6)   = HND5

      DO    I = 1, 5
       AN(I)  = AN(I) - RSLT*WM(I) + AUXI*(1 + AP(I))
      ENDDO

      IF(ALTU.EQ.500.) GOTO 800

      IF(ALTU.LT.500.) THEN
       TEMI   = TEMF
       STEP   = -5.
       ZINI   = 500.
       ZEND   = ALTU

       ZAUX   = (ZEND - ZINI)/4.
       AUXI   = INT(ZAUX/STEP)
       IF(AUXI.LE.0.) AUXI  = 1.
       STEP   = ZAUX/AUXI

       DO WHILE (ABS(ZEND-ZINI).GT.0.0001)
        AL(1) = AN(1)
        AL(2) = AN(2) - 
     1        0.07*(1. + TANH(0.18*(ZINI - 111.)))*CONV
        AL(3) = AN(3)
        AL(4) = AN(4)
        AL(5) = AN(5) - 
     1        0.24*EXP(-0.009*(ZINI-97.7)*(ZINI-97.7))*CONV

        ZAUX  = ZINI + 2.*STEP
        AUXI  = 1. + ZAUX/RT
        GRAV  = GO/AUXI/AUXI/RGAS
        SUMM  = 0.
        SUMO  = 0.
        ZAUX  = ZINI

        DO  I = 1, 5
         SUMM = SUMM + GRAV*WT(I)/TEMLO(ZAUX,AA)
         SUMO = SUMO + EXP(AL(I))
         ZAUX = ZAUX + STEP
        ENDDO

        ZINI  = ZINI + 4.*STEP
        RSLT  = STEP*SUMM
        TEMF  = TEMLO(ZINI,AA)
        TITF  = LOG(TEMI/TEMF)
        TEMI  = TEMF

        DO  I = 1, 5
         AN(I)= AN(I) - RSLT*WM(I) + TITF*(1. + AP(I))
        ENDDO

        SOMA  = SUMO/EXP(AN(6))*PHIS
        ZAUX  = ZINI + 2.*STEP
        AUXI  = 1. + ZAUX/RT
        GRAV  = GO/AUXI/AUXI/RGAS
        ZAUX  = ZINI
        SUMM  = 0.
        SUMO  = 0.
        SUMH  = 0.

        DO  I = 1, 5
         TEMP = TEMLO(ZAUX,AA)
         SUMM = SUMM + WT(I)/TEMP
         AUXI = WT(I)/SQRT(TEMP)
         SUMO = SUMO + SOMA*AUXI
         SUMH = SUMH + AUXI 
         ZAUX = ZAUX + STEP
        ENDDO

        AUXI  = GRAV*SUMM*STEP*WM(6)
        RSLT  = AUXI - TITF*(1. + AP(6))
        AUXI  = 1000.*STEP
        RSLS  = SUMO*AUXI
        AUXI  = SUMH*PHIS*AUXI
        AN(6) = AN(6) - RSLT - RSLS - AUXI
       ENDDO
      ELSE
       ZINI   = 500.
       STEP   = 2.5
       RSLT   = 0.
       RSLS   = 0.
       TEMI   = TEMF
       ZEND   = ALTU

       ZAUX   = (ZEND - ZINI)/4.
       AUXI   = INT(ZAUX/STEP)
       IF(AUXI.LE.0.) AUXI  = 1.
       STEP   = ZAUX/AUXI

       DO WHILE (ABS(ZEND-ZINI).GT.0.0001)
        ZAUX  = ZINI + 2.*STEP
        AUXI  = 1. + ZAUX/RT
        GRAV  = GO/AUXI/AUXI/RGAS
        SUMM  = 0.
        SUMO  = 0.
        ZAUX  = ZINI

        DO  I = 1, 5
         TEMP = TEMLO(ZAUX,AA)
         SUMM = SUMM + GRAV/TEMP*WT(I)
         SUMO = SUMO + WT(I)/SQRT(TEMP)
         ZAUX = ZAUX + STEP
        ENDDO

        RSLT  = RSLT + STEP*SUMM
        RSLS  = RSLS + 1000.*STEP*PHIS*SUMO
        ZINI  = ZINI + 4.*STEP
       ENDDO

       TEMF   = TEMLO(ZEND,AA)
       AUXI   = LOG(TEMI/TEMF)
  
       DO   I = 1, 6
        AN(I) = AN(I) - RSLT*WM(I) + AUXI*(1. + AP(I))
       ENDDO

       AN(6)  = AN(6) - RSLS
      ENDIF

  800 CONTINUE
      AN(2)   = AN(2) - 
     1        0.07*(1. + TANH(0.18*(ZEND - 111.)))*CONV
      AUXI    = ZEND - 97.7
      AN(5)   = AN(5) - 0.24*EXP(-0.009*AUXI*AUXI)*CONV
      WEIG    = 0.
      SUMM    = 0.

      DO    I = 1, NC
       AUXI   = EXP(AN(I))
       AN(I)  = AN(I)/CONV
       IF(AN(I).LT.0.) AN(I)   = 0.
       WEIG   = WEIG + AUXI*WM(I)
       SUMM   = SUMM + AUXI
      ENDDO

      WMOL    = WEIG/SUMM
      RHOD    = WEIG/AVOG

      RETURN
      END
