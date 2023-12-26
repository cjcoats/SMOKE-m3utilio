
SUBROUTINE INITSTCY

!***********************************************************************
!  subroutine INITSTCY body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to initialize the necessary fields
!      for performing state and county totals.  The first call sets up the
!      indices from each source to each county.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 8/99 by M. Houyoux
!
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
!****************************************************************************

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: NASRC,          NMSRC,  NPSRC,          &! no. of sources by category
    &                    AIFIP,          MIFIP,  PIFIP,          &! country/state/county codes
    &                    AFLAG,  BFLAG,  MFLAG,  PFLAG,  XFLAG,  &! source type flags
    &                    LREPSTA,                                &! output state total emissions flag

    &                    AEBSTA, BEBSTA, MEBSTA, PEBSTA, TEBSTA, &! state total speciated emissions
    &                    AEUSTA,         MEUSTA, PEUSTA, TEUSTA, &! state total mult-control emissions
    &                    AERSTA,         MERSTA, PERSTA, TERSTA, &! state total reac-control emissions
    &                    AECSTA,         MECSTA, PECSTA, TECSTA, &! state total all-control emissions

    &                    AEBCNY, BEBCNY, MEBCNY, PEBCNY, TEBCNY, &! county total speciated emissions
    &                    AEUCNY,         MEUCNY, PEUCNY, TEUCNY, &! county total mult-control emissions
    &                    AERCNY,         MERCNY, PERCNY, TERCNY, &! county total reac-control emissions
    &                    AECCNY,         MECCNY, PECCNY, TECCNY  ! county total all-control emissions

!.........  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: AICNY, MICNY, PICNY, NCOUNTY, CNTYCOD

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(2)    CRLF
    LOGICAL         ENVYN
    INTEGER         FINDC

    EXTERNAL   CRLF, ENVYN, FINDC

!...........   Other local variables

    INTEGER          IOS      ! i/o status
    INTEGER          J        ! counter

    LOGICAL, SAVE :: FIRSTIME = .TRUE. ! true: first time routine called

    CHARACTER(300)   MESG     ! message buffer

    CHARACTER(16) :: PROGNAME = 'INITSTCY' ! program name

!***********************************************************************
!   begin body of subroutine INITSTCY

!.........  Read surrogates (if needed) and state/county names
    IF( FIRSTIME ) THEN

!.............  Allocate memory for indices from Co/st/cy codes to counties
        ALLOCATE( AICNY( NASRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AICNY', PROGNAME )
        ALLOCATE( MICNY( NMSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MICNY', PROGNAME )
        ALLOCATE( PICNY( NPSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PICNY', PROGNAME )

!.............  Create indices to counties from Co/st/cy codes for each source
!               category
        IF( AFLAG ) THEN
            CALL SET_COUNTY_INDEX( NASRC, AIFIP, AICNY )
        END IF

        IF( MFLAG ) THEN
            CALL SET_COUNTY_INDEX( NMSRC, MIFIP, MICNY )
        END IF

        IF( PFLAG ) THEN
            CALL SET_COUNTY_INDEX( NPSRC, PIFIP, PICNY )
        END IF

        FIRSTIME = .FALSE.

    END IF

!.........  Initialize totals to zero...
!.........  State totals...
    IF( LREPSTA ) THEN

        IF( AFLAG ) THEN
            AEBSTA = 0.   ! state total inven or speciated emissions
            IF ( ALLOCATED( AEUSTA ) ) THEN
                AEUSTA = 0.   ! state total multipltv-controlled emissions
            END IF
            IF ( ALLOCATED( AERSTA ) ) THEN
                AERSTA = 0.   ! state total reactivity-controlled emissions
            END IF
            IF ( ALLOCATED( AECSTA ) ) THEN
                AECSTA = 0.   ! state total all-controlled emissions
            END IF
        END IF

        IF( BFLAG ) THEN
            BEBSTA = 0.
        END IF

        IF( MFLAG ) THEN
            MEBSTA = 0.
            IF ( ALLOCATED( MEUSTA ) ) THEN
                MEUSTA = 0.
            END IF
            IF ( ALLOCATED( MERSTA ) ) THEN
                MERSTA = 0.
            END IF
            IF ( ALLOCATED( MECSTA ) ) THEN
                MECSTA = 0.
            END IF
        END IF

        IF( PFLAG ) THEN
            PEBSTA = 0.
            IF ( ALLOCATED( PEUSTA ) ) THEN
                PEUSTA = 0.
            END IF
            IF ( ALLOCATED( PERSTA ) ) THEN
                PERSTA = 0.
            END IF
            IF ( ALLOCATED( PECSTA ) ) THEN
                PECSTA = 0.
            END IF
        END IF

        IF( XFLAG ) THEN
            TEBSTA = 0.
            IF ( ALLOCATED( TEUSTA ) ) THEN
                TEUSTA = 0.
            END IF
            IF ( ALLOCATED( TERSTA ) ) THEN
                TERSTA = 0.
            END IF
            IF ( ALLOCATED( TECSTA ) ) THEN
                TECSTA = 0.
            END IF
        END IF

    END IF

!.........  County totals...
    IF( AFLAG ) THEN
        AEBCNY = 0.   ! county total inven or speciated emissions
        IF ( ALLOCATED( AEUCNY ) ) THEN
            AEUCNY = 0.   ! county total multiplicative-controlled emissions
        END IF
        IF ( ALLOCATED( AERCNY ) ) THEN
            AERCNY = 0.   ! county total reactivity-controlled emissions
        END IF
        IF ( ALLOCATED( AECCNY ) ) THEN
            AECCNY = 0.   ! county total all-controlled emissions
        END IF
    END IF

    IF( BFLAG ) THEN
        BEBCNY = 0.
    END IF

    IF( MFLAG ) THEN
        MEBCNY = 0.
        IF ( ALLOCATED( MEUCNY ) ) THEN
            MEUCNY = 0.
        END IF
        IF ( ALLOCATED( MERCNY ) ) THEN
            MERCNY = 0.
        END IF
        IF ( ALLOCATED( MECCNY ) ) THEN
            MECCNY = 0.
        END IF
    END IF

    IF( PFLAG ) THEN
        PEBCNY = 0.
        IF ( ALLOCATED( PEUCNY ) ) THEN
            PEUCNY = 0.
        END IF
        IF ( ALLOCATED( PERCNY ) ) THEN
            PERCNY = 0.
        END IF
        IF ( ALLOCATED( PECCNY ) ) THEN
            PECCNY = 0.
        END IF
    END IF

    IF( XFLAG ) THEN
        TEBCNY = 0.
        IF ( ALLOCATED( TEUCNY ) ) THEN
            TEUCNY = 0.
        END IF
        IF ( ALLOCATED( TERCNY ) ) THEN
            TERCNY = 0.
        END IF
        IF ( ALLOCATED( TECCNY ) ) THEN
            TECCNY = 0.
        END IF
    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS   ******************************

CONTAINS

!.............  This subroutine creates the county to src FIPs index
    SUBROUTINE SET_COUNTY_INDEX( NSRC, CIFIP, ICNY )

!.............  Subprogram arguments
        INTEGER,            INTENT (IN) :: NSRC
        CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )
        INTEGER,            INTENT(OUT) :: ICNY( NSRC )

!.............  Local variables
        INTEGER     J, S     ! counters and indices

        CHARACTER(FIPLEN3) CFIP      ! tmp cy/st/co code
        CHARACTER(FIPLEN3) PFIP      ! previous cy/st/co code

!----------------------------------------------------------------------------

        PFIP = ' '
        DO S = 1, NSRC

            CFIP = CIFIP( S )

            IF( CFIP .NE. PFIP ) THEN

                J = MAX( FINDC( CFIP, NCOUNTY, CNTYCOD ), 0 )
                PFIP = CFIP

            END IF

            ICNY( S ) = J

        END DO

    END SUBROUTINE SET_COUNTY_INDEX

END SUBROUTINE INITSTCY
