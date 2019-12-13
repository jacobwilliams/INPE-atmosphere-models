C     LIBRARY LAFON
C
C----
C
C     Purpose: To compute the atmospheric density, tempe-
C              perature and molecular weigh using the
C              Lafontaine & Hughes analitical model.
C
C----
C
      SUBROUTINE ADYMOS(SA,SU,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C----
C
C
C PURPOSE:
C
C     THE SUBROUTINE ADYMOS GIVES THE  TEMPERATURE,  DENSITY
C
C     AND MOLECULAR WEIGHT OF ATMOSPHERE, USING THE ANALITI-
C                                                   -
C     CAL STATIC AND DYNAMIC  MODEL,  WITH  THE  SOLAR  FLUX
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
C           RENCE (2), IN KELVIN.
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
C     ASDAMO
C
C REFERENCES:
C
C     (1)   LAFONTAINE, J.; HUGHES, P.  AN ANALYTIC VERSION   OF
C           JACCHIA'S 1977 MODEL ATMOSPHERE. "CELESTIAL   MECHA-
C           NICS" 29(3-26) 1983.
C
C     (2)   JACCHIA, L. G. "THERMOSPHERIC TEMPERATURE,   DENSITY
C           AND COMPOSTION: NEW MODELS." CAMBRIDGE,   MA,   SAO,
C           1977. (SAO SPECIAL REPORT NO 375).
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
C----
C
      RLAT    = SA(2)
      RLON    = SA(1) - GSTI
      SIFI    = 0.9792*SIN(RLAT) + 
     1          0.2028*COS(RLAT)*COS(RLON - 5.07891)
      COFI    = SQRT(1 - SIFI*SIFI)
      TAUO    = 2.4 + 4.8*COFI*COFI

      CALL SOFLUD(RJUD,DAFR,SD,OUTR)

      IF(OUTR.NE.0.) THEN
       WRITE(6,*) ' ERROR IN ROUTINE ADYMOS: OUTR = ',
     1            INT(OUTR)
       STOP
      ENDIF

      ND      = (DAFR/3600. - SD(6) + 12. - TAUO)/3.
      SF(1)   = SD(9)
      SF(2)   = SD(10)
      SF(3)   = SD(ND)

      CALL ASDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,RHOD)

      RETURN
      END

      SUBROUTINE ASMODS(ALTU,RJUD,DAFR,TE,AL,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE  SUBROUTINE  ASMODS  USES  THE  ANALITICAL  STATIC 
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
C           RENCE (2), IN KELVIN.
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
C     ASMADE
C
C REFERENCES:
C
C   (1) LAFONTAINE, J.; HUGHES, P.  AN ANALYTIC VERSION   OF
C       JACCHIA'S 1977 MODEL ATMOSPHERE. "CELESTIAL   MECHA-
C       NICS" 29(3-26) 1983.
C
C   (2) JACCHIA, L. G. "THERMOSPHERIC TEMPERATURE,   DENSITY
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
       WRITE(6,*) ' ERROR IN ROUTINE ASMODS: OUTR = ',
     1            INT(OUTR)
       STOP
      ENDIF

      SF(1)   = SD(9)
      SF(2)   = SD(10)
      VARI    = .154*SD(7)
      SF(3)   = 1.89*DLOG(VARI + DSQRT(VARI*VARI + 1.D0))

      CALL ASMADE(ALTU,SF,TE,AL,WMOL,RHOD)

      RETURN
      END


      SUBROUTINE ASDAMO(SA,SU,SF,RJUD,DAFR,GSTI,TE,AD,WMOL,
     1                  RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  ASDAMO  GIVES  THE  DENSITY, MOLECULAR
C
C     WEIGHT AND TEMPERATURE OF THE UPPER ATMOSPHERE,  USING
C
C     THE ANALITICAL STATIC VERSION (LAFONTAINE -HUGHES) AND
C         -          -
C     THE DYNAMIC ATMOSPHERIC MODEL (JACCHIA 1977).
C         -       -           --
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
C       ASMADE
C       DIVARA
C       GEOACA
C       SEALAT
C       SEMIAN
C
C REFERENCES:
C
C       [1]      LAFONTAINE, J.; HUGHES, P. AN ANALYTIC VER-
C                SION OF  JACCHIA'S 1977  MODEL  ATMOSPHERE.
C                "CELESTIAL MECHANICS" 29(3-26) 1983.
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA       APR/87            V. 1.0
C            BENTO SILVA DE MATOS
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
      TYFR   = YEFR - IDINT(YEFR)
      GEAC   = SF(3)
      RLAT   = SA(2)
      RLON   = SA(1) - GSTI
      HEIG   = SA(3)
      ALTU   = HEIG/1.D03
      SUDC   = SU(2)
  
      CALL ASMADE(HEIG,SF,TE,AL,WMOL,RHOD)
  
      TEXO   = TE(1)
  
      CALL DIVARA(TEXO,SA,SU,WMOL,AL)
  
      DO I   = 1, 6
       AD(I) = AL(I)
      ENDDO
  
      TQUT   = TH(6)
  
      CALL GEOACA(TQUT,GEAC,RLAT,RLON,ALTU,AL)
  
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
  
  
      SUBROUTINE ASMADE(ALTU,SF,TE,AL,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  USES  THE  ANALYTIC  STATIC  MODEL  OF
C                                -         -       -
C     LAFONTAINE-HUGHES TO COMPUTE THE  ATMOSPHERIC  DENSITY
C                                       -            --
C     FOR HEIGHTS FROM 100 TO 2000 KM.
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
C       AMOWEI
C       TEMLO
C
C REFERENCES:
C
C       [1]      LAFONTAINE, J.; HUGHES, P. AN ANALYTIC VER-
C                SION OF  JACCHIA'S 1977  MODEL  ATMOSPHERE.
C                "CELESTIAL MECHANICS" 29(3-26) 1983.
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA        APR/87            V 1.0
C            BENTO SILVA DE MATOS
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
  
      CALL AMOWEI(THAF,HEIG,AL,WMOL,RHOD)
  
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
  
  
      SUBROUTINE DIVARA(TEXO,SA,SU,WMOL,AL)
C
C------
C
C PURPOSE:
C
C     THE SUBROUTINE  CALCULATES THE  DIURNAL VARIATIONS  IN
C                                     --      ---
C     THE ATMOSPHERE DENSITY, AS PROPOSED BY JACCHIA,  USING
C
C     THE ANALYTIC STATIC MODEL BY LAFONTAINE-HUGHES.
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
C       AMOWEI
C
C REFERENCES:
C
C       [1]      LAFONTAINE, J.; HUGHES, P. AN ANALYTIC VER-
C                SION OF  JACCHIA'S 1977  MODEL  ATMOSPHERE.
C                "CELESTIAL MECHANICS" 29(3-26) 1983.
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA       APR/87             V 1.0
C            BENTO SILVA DE MATOS
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
       CALL AMOWEI(TEMP,ALTU,AN,WEIG,RHOD)
       AL(I)  = AN(I)
      ENDDO
  
      DO   I = 1, 6
       AC(I)  = AN(I)
      ENDDO
  
      RETURN
      END
  
  
      SUBROUTINE GEOACA(TQUT,GEAC,RLAT,RLON,ALTU,DN)
C
C------
C
C PURPOSE:
C
C     THE  SUBROUTINE  GEOACA  OBTAINS  THE VARIATION OF THE
C
C     ATMOSPHERE  NUMBER-DENSITY,  DUE  TO  THE  GEOMAGNETIC
C                                                ---
C     ACTIVITY, USING THE ANALYTIC STATIC MODEL.
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
C       AMOWEI
C
C REFERENCES:
C
C       [1]      LAFONTAINE, J.; HUGHES, P. AN ANALYTIC VER-
C                SION OF  JACCHIA'S 1977  MODEL  ATMOSPHERE.
C                "CELESTIAL MECHANICS" 29(3-26) 1983.
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA       APR/87             V 1.0
C            BENTO SILVA DE MATOS
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
  
      CALL AMOWEI(TEMP,ALTU,DN,WMOL,RHOD)
  
      DELT   = 0.01D0*DELT
      DEZH   = 5000.D0*DLOG(DSQRT(1.D0 + DELT*DELT) + DELT)
      DELE   = 5.2D-04*AUXI*COL2*COL2
  
      DO   I = 1, 6
       DN(I) = DN(I) + GM(I)*DEZH + DELE
      ENDDO
  
      RETURN
      END
  
  
      SUBROUTINE AMOWEI(TEXO,ALTU,AN,WMOL,RHOD)
C
C------
C
C PURPOSE:
C
C     THIS SUBROUTINE  USES THE ANALITIC  LAFONTAINE-HUGHE'S
C                               -
C     VERSION OF THE JACCHIA'S 1977  MODEL  TO  COMPUTE  THE
C                                    --
C     MEAN MOLECULAR WEIGHT OF THE ATMOSPHERE IN  THE  LOCAL
C                    ---
C     REGION.
C
C INPUTS:
C
C       TEXO     MEAN EXOSPHERIC TEMPERATURE, IN KELVIN.
C       ALTU     GEOMETRICAL  HEIGHT   WHERE  THE  MOLECULAR
C                WEIGHT WILL BE CALCULATED, IN KM (90-2000).
C
C OUTPUTS:
C
C       AN(1)    LOGARITHM BASE 10 OF THE HE NUMBER-DENSITY.
C       AN(2)    LOGARITHM BASE 10 OF THE O2 NUMBER-DENSITY.
C       AN(3)    LOGARITHM BASE 10 OF THE N2 NUMBER-DENSITY.
C       AN(4)    LOGARITHM BASE 10 OF THE AR NUMBER-DENSITY.
C       AN(5)    LOGARITHM BASE 10 OF THE  O NUMBER-DENSITY.
C       AN(6)    LOGARITHM BASE 10 OF THE  H NUMBER-DENSITY.
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
C REFERENCES:
C
C       [1]      LAFONTAINE, J.; HUGHES, P. AN ANALYTIC VER-
C                SION OF  JACCHIA'S 1977  MODEL  ATMOSPHERE.
C                "CELESTIAL MECHANICS" 29(3-26) 1983.
C
C       [2]      JACCHIA, L. G. "THERMOSPHERIC  TEMPERATURE,
C                DENSITY AND COMPOSITION: NEW MODELS."  CAM-
C                BRIDGE, MA, SAO 1977. (SAO  SPECIAL  REPORT
C                375).
C
C AUTHOR:    VALDEMIR CARRARA       APR/87             V 1.0
C            BENTO SILVA DE MATOS
C
      IMPLICIT REAL*8 (A-H,O-Z)
  
      DIMENSION AN(6)
      DIMENSION CC(6),QI(6),AP(6),WM(6)
      DIMENSION CE(9),AA(9),BE(9),BN(6),CN(6)
C
C------
C  
      DATA PI   /3.14159265359D0/
      DATA CONV /2.302585092994D0/
      DATA AVOG /6.02217D+26/
  
      DATA CC   /28.575551D0,-0.472715D0,-0.118328D0,
     1           0.054405D0,0.009942D0,-0.005744D0/
      DATA QI   /0.78110D0,0.20955D0,0.009343D0,.5242D-5,
     1           0.D0,0.D0/
      DATA AP   /0.D0,0.D0,0.D0,-0.38D0,-0.25D0,0.D0/
      DATA WM   /28.0134,31.9988,39.948,4.0026,1.00797,
     1           15.9994/
      DATA RA   /6356.766D0/
      DATA GA   /9.806650D0/
      DATA CTGP /8.314320D0/
      DATA ZX   /125.0D0   /
      DATA ZO   /90.0D0    /
      DATA ZH   /100.0D0   /
      DATA TO   /188.0D0   /
      DATA RHOI /3.43D-06  /
      DATA EMEI /28.960D0  /
      DATA EPSI /1.D-06    /
  
      AKXI   = 2.D0*(RA + ZX)/(ZX - ZO)
      VH     = AKXI*(ZH - ZO)/(RA + ZH) - 1.D0
      WO     = (1.D0 - VH)/(1.D0 + VH)
      W1     = WO + 1.D0
  
      CE(1)  = CC(1) + WO*(CC(2) + WO*(CC(3) + WO*CC(4)))
      CE(2)  = W1*(CC(2) + WO*(2.D0*CC(3) + 3.D0*WO*CC(4)))
      CE(3)  = W1*W1*(CC(3) + 3.D0*WO*CC(4))
      CE(4)  = W1*W1*W1*CC(4)
      CE(5)  = 0.D0
      CE(6)  = 0.D0
      CE(7)  = 0.D0
      CE(8)  = 0.D0
      CE(9)  = 0.D0
  
      ARGU   = 0.0045D0*(TEXO - TO)
      AUXI   = ARGU/DSQRT(1.D0 + ARGU*ARGU)
      TX     = TO + 110.5D0*DATANH(AUXI)
      TAUO   = (TX + TO)/(TX - TO)
      TAU1   = 2.D0*TO*TX/(TO - TX)
      GDEX   = 1.9D0*(TX - TO)/(ZX - ZO)
      GDXE   = 0.475D0*TO*(RA + ZX)/TX/(RA + ZO)
  
      AA(5)  = 0.06205282D0*DLOG(TEXO + 213.9884D0) -
     1         0.6286968D0
      AA(6)  = 0.06555110D0*DLOG(TEXO - 329.6454D0) -
     1         0.1520990D0
      AA(1)  = AA(5) - GDXE - TAUO
      AA(2)  = AA(6) - GDXE + 1.5D0
      AA(3)  = GDXE - 2.D0*AA(5)
      AA(4)  = GDXE - 2.D0*AA(6) - 0.5D0
      AA(7)  = 0.D0
      AA(8)  = 0.D0
      AA(9)  = 0.D0
  
      Z      = ZH
      IF(ALTU.LT.ZH) Z      = ALTU
  
      VE     = AKXI*(Z - ZO)/(RA + Z) - 1.D0
      TAUC   = AA(1) + VE*(AA(2) + VE*(AA(3) + VE*(AA(4) +
     1         VE*(AA(5) + VE*AA(6)))))
      TCHA   = TAU1/TAUC
      GAVE   = GA*RA*RA/AKXI/(RA + ZO)
  
      DO   I = 1, 9
       BE(I)  = 0.D0
       DO  J = 1, I
        BE(I)= BE(I) + CE(J)*AA(I+1-J)
       ENDDO
      ENDDO
  
      BEAU   = 0.D0
      EMCH   = 0.D0
      VHEI   = 1.D0
      VEEI   = 1.D0
  
      DO   I = 1, 9
       VEEI  = VEEI*VE
       BEAU  = BEAU + BE(I)*(VEEI + 2.D0*(DBLE(I) -
     1         2.D0*DBLE(INT(I/2))) - 1.D0)/DBLE(I)
       EMCH  = EMCH + CE(I)*VHEI
       VHEI  = VEEI
      ENDDO
  
      BETA   = EMCH/TCHA
      ARGU   = DEXP(-GAVE*BEAU/CTGP/TAU1)
      RHOD   = BETA*ARGU*TO*RHOI/EMEI
      AUXI   = AVOG*RHOD
      DENU   = AUXI/EMCH
      ARGU   = AUXI/EMEI
  
      BN(1)  = QI(1)*ARGU
      BN(3)  = QI(3)*ARGU
      BN(4)  = QI(4)*ARGU
      CN(2)  = DENU*(EMCH*(QI(2) + 1.D0)/EMEI - 1.D0)
      CN(6)  = 2.D0*DENU*(1.D0 - EMCH/EMEI)
      AUXI   = DLOG(CN(2))/CONV -
     1         0.07D0*(1.D0 + DTANH(0.18D0*(Z - 111.D0)))
      BN(2)  = 10.D0**AUXI
      AUXI   = DLOG(CN(6))/CONV -
     1         0.24D0*DEXP(-.009D0*(Z-97.7D0)*(Z-97.7D0))
      BN(6)  = 10.D0**AUXI
      BN(5)  = 1.D0
  
      IF(ALTU.LE.ZH) GOTO 800
  
      Z      = ZX
      IF(ALTU.LT.ZX) Z      = ALTU
  
      VE     = AKXI*(Z - ZO)/(RA + Z) - 1.D0
      TCHA   = AA(1) + VE*(AA(2) + VE*(AA(3) + VE*(AA(4) +
     1         VE*(AA(5) + VE*AA(6)))))
      THAU   = AA(1) + VH*(AA(2) + VH*(AA(3) + VH*(AA(4) +
     1         VH*(AA(5) + VH*AA(6)))))
      THAU   = TAU1/THAU
      TCHA   = TAU1/TCHA
      USAR   = THAU/TCHA
      BN(2)  = CN(2)
      AUXI   = 0.D0
      VEEI   = VE
      VHEI   = VH
  
      DO  I  = 1, 6
       AUXI  = AUXI + AA(I)*(VEEI - VHEI)/DBLE(I)
       VEEI  = VEEI*VE
       VHEI  = VHEI*VH
      ENDDO
  
      DO   I = 1, 4
       ARGU  = -AUXI*WM(I)*GAVE/CTGP/TAU1
       USAR  = USAR**(1.D0 + AP(I))
       BN(I) = BN(I)*USAR*DEXP(ARGU)
      ENDDO
  
      ARGU   = -AUXI*WM(6)*GAVE/CTGP/TAU1
  
      CN(6)  = CN(6)*THAU/TCHA*DEXP(ARGU)
      CN(2)  = BN(2)
      AUXI   = DLOG(CN(2))/CONV -
     1         0.07D0*(1.D0 + DTANH(0.18D0*(Z - 111.D0)))
      BN(2)  = 10.D0**AUXI
      AUXI   = DLOG(CN(6))/CONV -
     1         0.24D0*DEXP(-.009D0*(Z-97.7D0)*(Z-97.7D0))
      BN(6)  = 10.D0**AUXI
  
      IF(ALTU.LE.ZX) GOTO 800
  
      DCAP   = DEXP(8.042617D0)*TEXO**(-0.4197668D0) +
     1         22.58421D0 -
     2         0.05719352D0*TEXO*DEXP(-TEXO/1187.417D0)
      GAMF   = (TEXO - TX)/GDEX/(RA + ZX)
      ALFA   = GAMF*DCAP/(GAMF*DCAP - 1.D0)
      ASUT   = DCAP/ALFA/(TX - TEXO)
      BSUT   = -2.D0*ASUT*TEXO - DCAP
      CSUT   = TEXO*(TEXO*ASUT + DCAP)
      ARGU   = ((ALFA/EPSI*(TX-TEXO) +
     1         1.D0)/(1.D0-ALFA))**(1.D0/DCAP)
      XINF   = DLOG(ARGU)
      ZINF   = (RA*XINF + ZX)/(1.D0 - XINF)
      HAGA   = ALTU
  
      IF(HAGA.GE.ZINF) HAGA   = ZINF
  
      XOFZ   = (HAGA - ZX)/(RA + HAGA)
      AUXI   = DEXP(-DCAP*XOFZ)
      ARGU   = ALFA*(TX - TEXO)*AUXI
      TCHA   = ARGU/(AUXI  + ALFA - 1.D0) + TEXO
      GAMA   = GA*RA*RA/CTGP/(RA + ZX)/CSUT
      BN(2)  = CN(2)
  
      ARGU   = (ASUT*TX + BSUT)*TX + CSUT
      ABBA   = ((ASUT*TCHA + BSUT)*TCHA + CSUT)/ARGU
      DTDO   = TCHA - TEXO
      DTD1   = TEXO - TX
      DTD2   = ALFA/(1.D0 - ALFA)
      AAAA   = DTD1*DTD2/DTDO + 1.D0/(1.D0 - ALFA)
  
      DO   I = 1, 4
       GAAU  = GAMA*WM(I)
       AUXI  = BN(I)*(TX/TCHA)**(1.D0 + AP(I) + GAAU)
       BN(I) = AUXI*(AAAA**(GAAU*BSUT/2.D0/DCAP))*
     1         (ABBA**(GAAU/2.D0))
      ENDDO
  
      CN(2)  = BN(2)
      AUXI   = DLOG(CN(2))/CONV -
     1         0.07D0*(1.D0 + DTANH(0.18D0*(HAGA - 111.D0)))
      BN(2)  = 10.D0**AUXI
      AUXI   = CN(6)*AAAA**(GAMA*WM(6)*BSUT/2.D0/DCAP)
      AUXI   = AUXI*ABBA**(GAMA*WM(6)/2.D0)
      CN(6)  = AUXI*(TX/TCHA)**(1.D0 + AP(6) + GAMA*WM(6))
      AUXI   = DLOG(CN(6))/CONV -  0.24D0*
     1         DEXP(-.009D0*(HAGA - 97.7D0)*(HAGA - 97.7D0))
      BN(6)  = 10.D0**AUXI
  
      HDEN   = 10.D0**(TEXO**(-.25D0)*28.9D0 + 5.94D0)
      X500   = (500.D0 - ZX)/(RA + 500.D0)
  
      AUXI   = DEXP(-X500*DCAP)
      T500   = ALFA*(TX - TEXO)*AUXI/(AUXI+ALFA-1.D0) + TEXO
      AUXI   = (TCHA*ASUT + BSUT)*TCHA + CSUT
      AUXI   = AUXI/((ASUT*T500 + BSUT)*T500 + CSUT)
      ARGU   = GAMA*WM(5)*BSUT/2.D0*(XOFZ - X500)
      AUXI   = HDEN*AUXI**(GAMA*WM(5)/2.D0)*DEXP(ARGU)
      BN(5)  = AUXI*(T500/TCHA)**(1.D0 + AP(5) + GAMA*WM(5))
  
      IF(HAGA.NE.ALTU) THEN

       X1    = 1.D0/(1.D0 + ZINF/RA)
       X2    = 1.D0/(1.D0 + ALTU/RA)
       COEF  = GA*RA/CTGP/TEXO*(X2 - X1)
       BN(2) = CN(2)
       BN(6) = CN(6)
  
       DO  I = 1, 6
        BN(I)= BN(I)*DEXP(COEF*WM(I))
        IF(BN(I).LE.0.D0) BN(I) = 1.D-32
       ENDDO
  
       AUXI  = DLOG(BN(2))/CONV -
     1         0.07D0*(1.D0 + DTANH(0.18D0*(ALTU - 111.D0)))
       BN(2) = 10.D0**AUXI
       AUXI  = DLOG(BN(6))/CONV -  0.24D0*
     1         DEXP(-.009D0*(ALTU - 97.7D0)*(ALTU - 97.7D0))
       BN(6) = 10.D0**AUXI
      ENDIF
  
      DO   I = 1, 6
       IF(BN(I).LT.1.D-30) BN(I)  = 1.D-30
      ENDDO
  
  800 CONTINUE
      WMOL   = 0.D0
      ARGU   = 0.D0
  
      DO   I = 1, 6
       AUXI  = BN(I)
       ARGU  = ARGU + AUXI
       WMOL  = WMOL + AUXI*WM(I)
      ENDDO
  
      RHOD   = WMOL/AVOG
      WMOL   = WMOL/ARGU
  
      AN(1)  = DLOG(BN(4))/CONV
      AN(2)  = DLOG(BN(2))/CONV
      AN(3)  = DLOG(BN(1))/CONV
      AN(4)  = DLOG(BN(3))/CONV
      AN(5)  = DLOG(BN(6))/CONV
      AN(6)  = DLOG(BN(5))/CONV
  
      RETURN
      END
  
