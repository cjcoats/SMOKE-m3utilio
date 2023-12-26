
SUBROUTINE HOURMET( NSRC, CNTYSRC, TA, QV, PRES, JDATE, JTIME,  &
                    DAYBEGT, LDAYSAV, PDTEMP, HFLAG )

    !***********************************************************************
    !  subroutine body starts at line 92
    !
    !  DESCRIPTION:
    !       Creates summed hourly temperatures by county. Checks that temperatures
    !       are within requested minimum and maximum. Keeps track of total number
    !       of sources for averaging later.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
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
    !****************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    USE MODMBSET, ONLY: NREFC, MCREFSORT, MCREFIDX,             &
                        NREFF, FMREFSORT, NFUELC, FMREFLIST

    !.......   This module is the derived meteorology data for emission factors
    USE MODMET, ONLY: TKHOUR, NTKHOUR, RHTBIN, NRHTBIN, NFUEL,  &
                      FUELIDX, MINTSRC, MAXTSRC, MAXTDAY, MINTDAY

    !.......   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NSRGFIPS, SRGFIPS, NCELLS, FIPCELL

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER,      INTENT    (IN) :: NSRC                      ! no. sources
    CHARACTER(*), INTENT    (IN) :: CNTYSRC( NSRC )           ! no. counties
    REAL   ,      INTENT    (IN) :: TA( * )                   ! gridded temp data
    REAL   ,      INTENT    (IN) :: QV( * )                   ! gridded mixing ratio data
    REAL   ,      INTENT    (IN) :: PRES( * )                 ! gridded pressure data
    INTEGER,      INTENT    (IN) :: JDATE                     ! YYYYDDD
    INTEGER,      INTENT    (IN) :: JTIME                     ! HHMMSS
    INTEGER,      INTENT    (IN) :: DAYBEGT ( NSRC )          ! begin. time for day
    LOGICAL,      INTENT    (IN) :: LDAYSAV ( NSRC )          ! true: use daylight time
    INTEGER,      INTENT    (IN) :: PDTEMP                    ! RPP temperature increment
    LOGICAL,      INTENT    (IN) :: HFLAG                     ! true: use specific humidity (no RH)

    !.......   EXTERNAL FUNCTIONS
    REAL   , EXTERNAL :: CALCRELHUM
    INTEGER, EXTERNAL :: FINDCFIRST

    !.......   Other local variables
    INTEGER     C, K, L, LL, N, S, I, J, NR, NF, NT, T           ! counters and indices
    INTEGER     IOS             ! I/O status
    INTEGER     MONTH,DAY       ! processing month and date
    INTEGER     TIMESLOT        ! array location
    INTEGER     CURMONTH, NMON

    INTEGER, SAVE :: MXWARN     ! maximum number of warnings
    INTEGER, SAVE :: NWARN      ! total number of warnings printed

    REAL        MINTMP          ! min temperature value in Farenheight
    REAL        MAXTMP          ! max temperature value in Farenheight
    REAL        TEMPVAL         ! temperature value in Farenheight
    REAL        TEMPTMP         ! tmp temperature value in Farenheight
    REAL        RHVAL           ! RH value

    LOGICAL       :: DAYLIT  = .FALSE.      ! true: date is daylight savings
    LOGICAL, SAVE :: INITIAL = .TRUE.       ! true: first time

    CHARACTER(FIPLEN3) REFCOUNTY     ! current ref. county

    CHARACTER(300)     BUFFER        ! formatted source info for messages
    CHARACTER(300)     MESG          ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'HOURMET'     ! program name

    !***********************************************************************
    !   begin body of subroutine HOURTEMP

    !.......  For the first time, initialize all entries to zero
    IF( INITIAL ) THEN
        TKHOUR = 0.      ! array

        !.......  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
        END IF
        NWARN = 0

        INITIAL = .FALSE.
    END IF

    !.......  If last day of month, process monthly averages
    CALL DAYMON( JDATE, MONTH, DAY )

    !.......  Loop through sources
    DO S = 1, NSRC

        !.......  Apply ungridding matrix from a (possible) subgrid to data on base
        !           grid.  If no subgrid, then XOFF and YOFF will be 1 and no problem.
        LL = FINDC( CNTYSRC( S ), NSRGFIPS, SRGFIPS )

        IF( LL < 1 ) CYCLE

        IF( NCELLS( LL ) > 0 ) THEN
            RHVAL = 0.0
            TEMPVAL = 0.0
        ELSE
            RHVAL   = BADVAL3
            TEMPVAL = BADVAL3
        END IF

        N = 0
        DO I = 1, NCELLS( LL )

            !.......  Count no of cell used in county-level averaging
            N = N + 1

            !.......  Get column and row from subgrid
            C = FIPCELL( I,LL )

            !.......  Convert K to F
            TEMPTMP = 1.8 * TA( C ) - 459.67

            !.......  Store min/max by source for fuelmonth output
            MAXTSRC( S ) = MAX( TEMPTMP, MAXTSRC( S ) )
            MINTSRC( S ) = MIN( TEMPTMP, MINTSRC( S ) )

            !.......  Store min/max by cell for daily output for SMOKE RPP processing
            MAXTDAY( C ) = MAX( TEMPTMP, MAXTDAY( C ) )
            MINTDAY( C ) = MIN( TEMPTMP, MINTDAY( C ) )

            !.......  Calculate RH using Temp. Pressure and mixing ratio values
            IF( HFLAG ) THEN
                RHVAL = QV( C ) / ( QV( C ) + 1 )                 ! Specific Humidity
            ELSE
                RHVAL = CALCRELHUM( TA( C ), PRES( C ), QV( C ) )             ! relative humidity
            END IF

            !.......  Store RH into temperature bins
            NT = 0
            DO T = -150, 200-PDTEMP, PDTEMP
                NT = NT + 1
                MINTMP = REAL( T ) - ( REAL( PDTEMP ) / 2.0 )
                MAXTMP = MINTMP + REAL( PDTEMP )

                IF ( MINTMP < TEMPTMP .AND. TEMPTMP <= MAXTMP ) THEN

                    REFCOUNTY = MCREFSORT( S,2 )
                    NR = FINDC( REFCOUNTY,NREFC, MCREFIDX( :,1 ) )

                    L = FINDCFIRST( REFCOUNTY, NREFF, FMREFSORT( :,1 ) )
                    K = FINDCFIRST( REFCOUNTY, NFUELC,FMREFLIST( :,1 ) )
                    NMON = STR2INT( FMREFLIST( K, 2 ) )               ! no month of ref county

                    !.......  Loop over months per ref. county
                    DO J = L, L + NMON - 1
                        CURMONTH  = STR2INT( FMREFSORT( J,3 ) )               ! processing  current month per ref. cou
                        IF( CURMONTH == MONTH ) NF = STR2INT( FMREFSORT( J,2 ) )              ! processing fuelmonth/co
                    END DO

                    RHTBIN ( NR,NF,NT ) = RHTBIN ( NR,NF,NT )  + RHVAL
                    NRHTBIN( NR,NF,NT ) = NRHTBIN( NR,NF,NT ) + 1

                END IF

            END DO

            TEMPVAL = TEMPVAL + TEMPTMP

        END DO

        TEMPVAL = TEMPVAL / N            ! averaged Temp by source

        !.......  Calculate time slot in output array for this time step
        !               Appropriate 24 hour time will be day starting time (12 AM in local
        !               time zone ) subtracted from met data time (in GMT)
        TIMESLOT = 1 + ( JTIME - DAYBEGT( S ) ) / 10000

        !.......  Restore daylight saving time if necessay
        DAYLIT = ISDSTIME( JDATE )
        IF( DAYLIT .AND. LDAYSAV( S ) ) THEN
            TIMESLOT = TIMESLOT - 1              ! substract 1hr on DST date
        END IF

        !.......  If timeslot is less than zero, add 24; if better data comes
        !               along, the old data will get overwritten (helps in case of
        !               one running one day)
        IF( TIMESLOT <= 0 ) THEN
            TIMESLOT = TIMESLOT + 24
        END IF

        IF( TEMPVAL > AMISS3 )THEN

            !.......  Store values in hourly arrays for mothly SMOKE-ready output
            TKHOUR( S,TIMESLOT )  = TKHOUR( S,TIMESLOT ) + TEMPVAL
            NTKHOUR( S,TIMESLOT ) = NTKHOUR( S,TIMESLOT ) + 1

        END IF

    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

94020 FORMAT( A, 4( 1X, F8.2, 1X, A ) )

END SUBROUTINE HOURMET
