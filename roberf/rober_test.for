      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION SA(3), SU(3), AN(6), SF(3), TE(2)
	DIMENSION T(5), WM(5), RD(5), RL(5)
	DATA T / 500, 700, 1000, 1500, 2400 /
	DATA RL / -60, -30, 0, 30, 60 /
	
      OPEN(8, FILE = "robe_prof.dat")

	WRITE(6, 1000)
	WRITE(6, 1010)
	WRITE(6, 1020)
	WRITE(8, 1000)
	WRITE(8, 1010)
	WRITE(8, 1020)
      DO IJ = 90, 2000, 10
	  ALTU = IJ
        DO II = 1, 5
          TEXO = T(II)
          CALL RMOWEI (TEXO,ALTU,AN,WMOL,RHOD)   
          RD(II)  = LOG10(RHOD)
	    WM(II)  = WMOL
        ENDDO
	  WRITE (6, 1030) INT(ALTU), WM(1), RD(1), WM(2), RD(2), WM(3),
	*         RD(3), WM(4), RD(4), WM(5), RD(5)
	  WRITE (8, 1030) INT(ALTU), WM(1), RD(1), WM(2), RD(2), WM(3),
	*         RD(3), WM(4), RD(4), WM(5), RD(5)
      ENDDO

	CLOSE(8)
      OPEN(8, FILE = "robe_stat.dat")
      SF(1)   = 100
      SF(2)   = 100
	WRITE(6, 1040) SF(1), SF(2)
	WRITE(6, 1050)
	WRITE(8, 1040) SF(1), SF(2)
	WRITE(8, 1050)
      DO IJ = 90, 2000, 10
	  ALTU = IJ
        CALL RSMADE (ALTU*1000., SF, TE, AN, WMOL, RHOD)   
	  WRITE (6, 1060) INT(ALTU), WMOL, LOG10(RHOD), TE(1), TE(2),
	*  (AN(I), I=1,6)
	  WRITE (8, 1060) INT(ALTU), WMOL, LOG10(RHOD), TE(1), TE(2),
	*  (AN(I), I=1,6)
      ENDDO
	
	CLOSE(8)
      OPEN(8, FILE = "robe_dynm.dat")
      SF(1)   = 100.
      SF(2)   = 100.
	SF(3)   = 5.
	SA(1)   = 0.
	SA(3)   = 600000.
	SU(1)   = 0.
	SU(2)   = 0.
	RJUD    = 22476
	DFRA    = 0.
	GSTI    = 293.49/180*3.141592

	WRITE(6, 1070) SF(1), SF(2), SF(3), SA(1), GSTI, SA(3)/1000.,
	* RJUD, DFRA
	WRITE(6, 1080) (INT(RL(I)), I = 1, 5)
	WRITE(8, 1070) SF(1), SF(2), SF(3), SA(1), GSTI, SA(3)/1000.,
	* RJUD, DFRA
	WRITE(8, 1080) (INT(RL(I)), I = 1, 5)
      DO IJ = -180, 180, 5
	  SA(1) = IJ*3.141592/180.
	  DO II = 1, 5
	    SA(2) = RL(II)*3.141592/180.
          CALL RSDAMO (SA, SU, SF, RJUD, DAFR, GSTI, TE, AN, WMOL, RHOD)
	    WM(II) = WMOL
	    RD(II) = LOG10(RHOD)
	  ENDDO   
      WRITE (6, 1090) IJ, WM(1), RD(1), WM(2), RD(2), WM(3),
	*         RD(3), WM(4), RD(4), WM(5), RD(5)
      WRITE (8, 1090) IJ, WM(1), RD(1), WM(2), RD(2), WM(3),
	*         RD(3), WM(4), RD(4), WM(5), RD(5)
      ENDDO

	CLOSE(8)

 1000 FORMAT(1X, 32(1H*), " Static Profile ", 32(1H*), /,
     * 29X, "Exospheric Temperature", /)
 1010 FORMAT("  alt           500K           700K          1000K
     *  1500K          2400K")
 1020 FORMAT(" (km) mmolms logdsty mmolms logdsty mmolms logdsty mmolms
     *logdsty mmolms logdsty")
 1030 FORMAT(1X, I4, 5(1X, F6.2, F8.3))
 1040 FORMAT(//, 1X, 26(1H*), " Static Model - Solar Cycle ", 26(1H*),  
     */, 27X, "F10.7 = ", F6.2, " Fbar = ", F6.2, /)
 1050 FORMAT("  Alt MMolMs LogDsty ExTemp LcTemp NMolHe NMolO2 NMolN2 NM
     *olAr NMmolO  NMolH")
 1060 FORMAT(1X, I4, 1X, F6.2, F8.3, 2F7.2, 6F7.3)
 1070 FORMAT(//, 1X, 32(1H*), " Dynamic Model ", 33(1H*),  
     */, 1X, " F10.7 = ", F6.2, ", Fbar = ", F6.2, ", Kp = ", F3.1, 
     *", Right ascension (deg) = ", F6.2,	/,
     *", Greenwich Sideral Time (rad) = ", F8.3,
     * " Altitude (km) = ", F6.2, /,
     *", Modified Julian Date = ", F6.0, ", UTC (s) = ", F9.3, /)
 1080 FORMAT(1X, 31X, " Latitude (deg) ", /,
     * 1X, "SunRA", 4X, I4, 4(11X, I4), /,
     * 1X, " deg mmolms logdsty mmolms logdsty mmolms logdsty mmolms log
     *dsty mmolms logdsty")
 1090 FORMAT(1X, I4, 5(1X, F6.2, F8.3))

      END 
