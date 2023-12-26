
SUBROUTINE WRTEMPROF( ODEV1, ODEV2, MYEAR, COUNTY, PMONTH,  &
                      PDTEMP, PPTEMP, THOUR, MAXTEMP,       &
                      MINTEMP, TEMPBIN, MINNORH )

    !***********************************************************************
    !  subroutine body starts at line 78
    !
    !  DESCRIPTION:
    !       Averages hourly meteorology data based on number of sources
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  none
    !
    !  REVISION  HISTORY:
    !       Created  2/2010 by B.H. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
    ! All Rights Reserved
    !
    ! Carolina Environmental Program
    ! University of North Carolina at Chapel Hill
    ! 137 E. Franklin St., CB# 6116
    ! Chapel Hill, NC 27599-6116
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the derived meteorology data for emission factors
    USE MODMET, ONLY:   RHTBIN, NRHTBIN

    USE MODMBSET, ONLY: NREFC, MCREFIDX

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    !.......  Subroutine arguments
    INTEGER     , INTENT(IN)  ::  ODEV1               ! MOVES RPP output file
    INTEGER     , INTENT(IN)  ::  ODEV2               ! MOVES RPD/RPV output file
    INTEGER     , INTENT(IN)  ::  MYEAR               ! modeling year
    CHARACTER(*), INTENT(IN)  ::  COUNTY              ! refCounty
    INTEGER     , INTENT(IN)  ::  PMONTH              ! fuel month
    INTEGER     , INTENT(IN)  ::  PDTEMP              ! RPP/RPV temperature increment
    INTEGER     , INTENT(IN)  ::  PPTEMP              ! RPP temperature increment
    REAL        , INTENT(IN)  ::  THOUR( 24 )         ! avg temp 24-hr profiles for refcounty/fuelmonth
    REAL        , INTENT(IN)  ::  MAXTEMP             ! max temp for refcounty/fuelmonth
    REAL        , INTENT(IN)  ::  MINTEMP             ! min temp for refcounty/fuelmonth
    REAL        , INTENT(IN)  ::  TEMPBIN             ! temp buffer
    INTEGER     , INTENT(IN)  ::  MINNORH             ! min no of data points for avg RH by tempbin

    !.......  Local array
    REAL   :: TKPRO( 24 ) = 0.0     !  24hr temp. profile
    REAL   :: TMPRO( 24 ) = 0.0     !  tmp 24hr temp. profile

    !.......   Other local variables
    INTEGER I, IT, J, K, L, LT, N, NF, NR, NT, NTB, MX, MN, R, S, T              ! counters and indices

    INTEGER IOS                           ! I/O status
    INTEGER IMAXT, IMINT                  ! current integer min/max temp.
    INTEGER MAXT , MINT                   ! tmp remainder of max/min

    REAL    TKMAX                         ! tmp max temperatures
    REAL    TKMIN                         ! tmp min temperatures
    REAL    TKDIF, DT, ST                 ! tmp DIFF of min/max temperatures
    REAL    TKMED                         ! tmp median temperatures
    REAL    RHAVG                         ! tmp relative humidity

    LOGICAL FIRSTIME

    CHARACTER(32)    TPROID                 ! temporal resolution header
    CHARACTER(300)   MESG                   ! message buffer
    CHARACTER(16), PARAMETER :: PROGNAME = 'WRTEMPROF'     ! program name

    !***********************************************************************
    !   begin body of subroutine WRTEMPROF

    !.......  Write out last ref. county min/max temp and avg RH
    IMAXT = INT( MAXTEMP )
    IMINT = INT( MINTEMP )

    TKMAX  = MAXVAL( THOUR )
    TKMIN  = MINVAL( THOUR )
    TKDIF  = ABS( TKMAX - TKMIN )
    TKMED  = 0.5 * ( TKMAX + TKMIN )

    TMPRO = 0.0
    DO T = 1, 24
        TMPRO( T ) = ( THOUR( T ) - TKMED ) / TKDIF
    END DO

    !.......  Calculate temp profiles based on a combination of temp bins
    MAXT = IMAXT + ( PPTEMP - ABS( MOD( IMAXT,PPTEMP ) ) )
    MINT = IMINT - ABS( MOD( IMINT,PPTEMP ) )

    IF( MINT <= 0 .AND. MINTEMP < 0  ) THEN
        MINT = IMINT - ( PPTEMP - ABS( MOD( IMINT,PPTEMP ) ) )
    END IF

    IT = 0
    !.......  Determine temperature bins based on PPTEMP
    DO MX = MAXT,MINT,-PPTEMP

        DO MN = MINT,MX,PPTEMP

            IT = IT + 1
            DT = MX - MN
            ST = MX + MN
            WRITE( TPROID,94040 ) 'M', MYEAR, PMONTH, IT

            TKPRO = 0.0
            !.......  Create 24 hr temp profile per temp bin
            DO T = 1,24
                TKPRO( T ) = TMPRO( T ) * DT + ST/2.0
            END DO

            !.......  Output for SMOKE ready input file
            WRITE( ODEV1,94050 ) COUNTY, PMONTH,&
                TRIM(TPROID), MINTEMP, MAXTEMP, (TKPRO( T ), T=1,24)

        END DO

    END DO

    !.......  Output avg RH by temperature bins for RPD/RPV modes
    !.......  Fill empty temperature bins for RPD/RPV modes
    NR = FINDC( COUNTY,NREFC, MCREFIDX( :,1 ) )
    NF = PMONTH
    NT  = 0
    NTB = 0
    DO T = -150, 200-PDTEMP, PDTEMP
        NT = NT + 1
        IF( NRHTBIN( NR,NF,NT ) < MINNORH ) THEN
            RHTBIN ( NR,NF,NT ) = 0.0
            NRHTBIN( NR,NF,NT ) = 0
        END IF
    END DO

    NT = 0
    FIRSTIME = .TRUE.
    DO T = -150, 200-PDTEMP, PDTEMP
        NT = NT + 1
        IF( FIRSTIME .AND. NRHTBIN( NR,NF,NT ) == 0 ) THEN
            NTB = NT
            FIRSTIME = .FALSE.
            CYCLE
        ELSE IF( RHTBIN( NR,NF,NT ) > 0.0 ) THEN
            RHTBIN ( NR,NF,NTB:NT ) = RHTBIN ( NR,NF,NT )
            NRHTBIN( NR,NF,NTB:NT ) = NRHTBIN( NR,NF,NT )
            NTB = NT + 1
            FIRSTIME = .TRUE.
            CYCLE
        ELSE IF( T == 200-PDTEMP ) THEN
            RHTBIN ( NR,NF,NTB:NT ) = RHTBIN ( NR,NF,NTB-1 )
            NRHTBIN( NR,NF,NTB:NT ) = NRHTBIN( NR,NF,NTB-1 )
        END IF
    END DO

    !.......  Write out refcounty min/max temp and avg RH by temperature bin for RPD/RPV
    !.......  Calculate max/min temp bins based on RPD_TEMP_INCREMENT
    MAXT = IMAXT + ( PDTEMP - ABS( MOD( IMAXT,PDTEMP ) ) )
    MINT = IMINT - ABS( MOD( IMINT,PDTEMP ) )

    IF( MINT <= 0 .AND. MINTEMP < 0  ) THEN
        MINT = IMINT - ( PDTEMP - ABS( MOD( IMINT,PDTEMP ) ) )
    END IF

    NT = 0
    DO T = -150, 200-PDTEMP, PDTEMP
        NT = NT + 1
        IF( MINT <= T .AND. T <= MAXT ) THEN

            RHAVG = RHTBIN(NR,NF,NT) / NRHTBIN(NR,NF,NT)

            WRITE( ODEV2,94060 ) COUNTY, PMONTH, RHAVG, MINTEMP,&
                   MAXTEMP, T

        END IF
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94040 FORMAT( A, I4.4, I2.2, I4.4 )

94050 FORMAT( A,',', I5,',', 3X, A, 26(',',F10.2) )

94060 FORMAT( A,',',I5,',',F12.6,',',F10.2,',',F10.2,',',I5 )

END SUBROUTINE WRTEMPROF
