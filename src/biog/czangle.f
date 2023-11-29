        SUBROUTINE  CZANGLE( JDATE, JTIME, NX, NY, LAT, LON, COSZEN,
     &                       ZFLAG )

C***********************************************************************
C  subroutine body starts at line 101
C
C  DESCRIPTION:
C       Computes cosine of zenith angle for routine HRBIO()
C
C  PRECONDITIONS REQUIRED:
C       JDATE:JTIME represented in GMT
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION  HISTORY:
C       Prototype 12/95 by Carlie J Coats, Jr., adapted from UAM-BEIS
C       subroutines SOLAR() and ZANGLE() for SMOKE-BEIS2:  produces
C       COS( ZENITH )
C
C       Revised 8/96 by SL and CJC: algorithm change to match UAM BEIS2
C       algorithm
C
C       Version 11/99: by Jeff Vukovich taken from v4.2 SMOKE prototype
C
C       Adapted 10/2023 to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***********************************************************************

        USE M3UTILIO

        IMPLICIT NONE

        !!...........   ARGUMENTS and their descriptions:

        INTEGER, INTENT (IN)  :: JDATE   !  current simulation date (YYYYDDD)
        INTEGER, INTENT (IN)  :: JTIME   !  current simulation time (HHMMSS)
        INTEGER, INTENT (IN)  :: NX  !  no. columnse
        INTEGER, INTENT (IN)  :: NY  !  no. rows

        REAL,    INTENT (IN)  :: LAT  ( NX, NY )  !  lat (deg) -90 <= LAT <= 90
        REAL,    INTENT (IN)  :: LON  ( NX, NY )  !  lon (deg) -180 <= LON <= 180

        REAL,    INTENT (OUT) :: COSZEN ( NX, NY) !  cos of zenith angle
        LOGICAL, INTENT (OUT) :: ZFLAG            !  iff sun above horizon somewhere


        !!.......   PARAMETERS:

        REAL, PARAMETER :: PI     =   3.14159265358979324
        REAL, PARAMETER :: PI180  = PI / 180.0
        REAL, PARAMETER :: AA     =   0.15
        REAL, PARAMETER :: BB     =   3.885
        REAL, PARAMETER :: CC     = - 1.253
        REAL, PARAMETER :: SIGA   = 279.9348
        REAL, PARAMETER :: D60    = 1.0 /   60.0
        REAL, PARAMETER :: D3600  = 1.0 / 3600.0
        REAL, PARAMETER :: D15    = 1.0 /   15.0
        REAL, PARAMETER :: D24    = 1.0 /   24.0
        REAL, PARAMETER :: ROTDAY = 360.0 / 365.242       ! fraction of a complete rotation per day
        REAL, PARAMETER :: SDEC   =   0.39784984          !  SIN (23^26'37.8") the declination angle

        !!.........   Scratch Local variables

        INTEGER    IOS, R, C
        REAL       SLA, GMT,  TK, DAD, DF,
     &             DESIN, DECOS, DESIN2, DECOS2, SIG, DECSIN, DECCOS,
     &             EQT, TST, HRANGL

        !!.......   SAVED local variables

        LOGICAL, SAVE ::   FIRSTIME = .TRUE.

        REAL, ALLOCATABLE, SAVE :: SINLAT ( :, : )
        REAL, ALLOCATABLE, SAVE :: COSLAT ( :, : )

        CHARACTER(16) :: PROGNAME = 'CZANGLE'   !  program name

C***********************************************************************
C................  Begin body of program  .........................

C........ compute sine of lat and lon first time through

      IF ( FIRSTIME ) THEN

          FIRSTIME = .FALSE.

          ALLOCATE( SINLAT ( NX, NY  ), COSLAT ( NX, NY ), STAT=IOS )
          IF ( IOS .NE. 0 ) THEN
              CALL M3EXIT( 'CZANGLE',0,0, 'Memory allocation error', 2 )
          END IF

          DO  R = 1, NY
          DO  C = 1, NX
              SLA = PI180 * LAT( C,R )
              SINLAT( C,R ) = SIN( SLA )
              COSLAT( C,R ) = COS( SLA )
          END DO
          END DO

      END IF    !  if firstime


C.......   Convert time to hours and add time-zone offset

        GMT    =  FLOAT( JTIME / 10000 ) +                        !  hr part
     &              D60 * ( FLOAT( MOD( JTIME / 100 , 100 ) )     !  min part
     &                    + D60 * FLOAT( MOD( JTIME, 100 ) ) )    !  sec part
        DAD    =  GMT * D24 + MOD( JDATE , 1000 )
        DF     =  ROTDAY * PI180 * DAD      !  The terrestrial-rotation angle

        DESIN  =  SIN( DF )          !  SINE   of this angle
        DECOS  =  COS( DF )          !  COSINE of this angle

        DESIN2 =  SIN( DF + DF )     !  SINE   of twice the angle
        DECOS2 =  COS( DF + DF )     !  COSINE of twice the angle

        SIG  =  DF +
     &          PI180 * ( SIGA +
     &                    1.914827 * DESIN  - 0.079525 * DECOS +
     &                    0.019938 * DESIN2 - 0.00162  * DECOS2 )

        !!  The sine of the declination

        DECSIN = SDEC * SIN( SIG )
        DECCOS = SQRT( 1.0 - DECSIN*DECSIN )

        !!  The equation of time adjustment

        EQT  =  0.123470 * DESIN  - 0.004289 * DECOS
     &        + 0.153809 * DESIN2 + 0.060783 * DECOS2

        ZFLAG = .FALSE.

        DO R = 1, NY
        DO C = 1, NX

            TK     =  GMT + LON( C,R ) * D15      !  Distance in hours from LON=0
            TST    =  TK - EQT                    !  true solar time
            HRANGL =  PI180 * 15.0 * ABS(TST - 12.0) !  hour angle

              !!  Compute the sine of the solar elevation

            COSZEN( C,R ) = DECSIN * SINLAT( C,R ) +
     &                      DECCOS * COSLAT( C,R ) * COS( HRANGL )
            ZFLAG = ( ZFLAG .OR. ( COSZEN( C,R ) .GT. -0.01 ) )  ! "sun > horizon"

        END DO
        END DO


        RETURN

      END SUBROUTINE  CZANGLE
