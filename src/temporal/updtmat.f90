
SUBROUTINE UPDTMAT( NSRC, IGRP, NPOL, JDATE, TZONE, VIDX, HIDX,    &
                    MONTH, DAYOW, DAYOM )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !       This routine updates the temporal matrix for sources that have source-
    !       specific hourly profiles.  Since the hourly profile fractions are read
    !       in one variable and hour at a time, this routine is called for each
    !       hour and each pollutant/activity.
    !
    !  PRECONDITIONS REQUIRED:
    !       Source-specific hourly profiles read for VIDX and HIDX.  MONTH and
    !       DAYOW arrays populated
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  none
    !
    !  REVISION  HISTORY:
    !       Copied from UPDTMAT.F 2.6 by M Houyoux 1/99
    !
    !       Version 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
    !
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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

    !.......  MODULES for public variables
    !.......  MODSOURC contains the inventory arrays
    !.......  MODTMPRL contains the temporal profile tables
    !.......  MODDAYHR contains data for day- and hour-specific data

    USE MODSOURC, ONLY: TPFLAG, TZONES

    USE MODTMPRL, ONLY: MONFAC,  DOMFAC,  WEKFAC,  XWKFAC, IPOL2D,    &
                        MTHPROF, DOMPROF, WEKPROF

    USE MODDAYHR, ONLY: INDXD, INDXH, NHRSRC, NDYSRC, EMACH, LDSPOA, TMAT

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS:

    INTEGER, INTENT(IN ) :: NSRC                        ! number of sources
    INTEGER, INTENT(IN ) :: IGRP                        ! pollutant group
    INTEGER, INTENT(IN ) :: NPOL                        ! number of pollutants
    INTEGER, INTENT(IN ) :: JDATE                       ! date YYYYDDD
    INTEGER, INTENT(IN ) :: TZONE                       ! time zone (5 for Eastern)
    INTEGER, INTENT(IN ) :: VIDX                        ! pol/act index
    INTEGER, INTENT(IN ) :: HIDX                        ! hour index
    INTEGER, INTENT(IN ) :: MONTH( 24, -23:23 )           ! source time zone's month-of-year 1...12
    INTEGER, INTENT(IN ) :: DAYOW( 24, -23:23 )           ! source time zone's day-of-week   1...7
    INTEGER, INTENT(IN ) :: DAYOM( 24, -23:23 )           ! source time zone's day-of-month  1...31

    !.......   Other Local variables:

    INTEGER         I, J, S, V      ! counters and indices

    INTEGER         DAY, DOM            !  day for source and hour pointer
    INTEGER         MON                 !  month for source and hour pointer
    INTEGER         IMON, IWEK, IDOM

    REAL            FAC                 !  partial matrix factor
    REAL            YRFAC               !  year to day factor

    CHARACTER(16), PARAMETER :: PROGNAME = 'UPDTMAT'      ! program name

    !***********************************************************************
    !   begin body of subroutine  UPDTMAT

    !.......  Compute correct year-to-day conversion factor:

    YRFAC = YR2DAY( JDATE / 1000 )
    V     = IPOL2D( VIDX,IGRP )

    !.......  Compute TMAT for current group of pollutants

    DO I = 1, NHRSRC

        S = INDXH( I )

        IMON = MTHPROF( S,V )
        IWEK = WEKPROF( S,V )
        IDOM = DOMPROF( S,V )

        MON = MONTH( HIDX, TZONES( S ) )
        DAY = DAYOW( HIDX, TZONES( S ) )
        DOM = DAYOM( HIDX, TZONES( S ) )

        !.......  Skip if source index is 0 (past last source for current hour),
        !        or if profile value is missing

        IF( S .EQ. 0 .OR. EMACH( I ) .LE. AMISS3 ) CYCLE

        !.......  Use day-specific data (no adjustments for month or weekday)

        J = FIND1( S, NDYSRC, INDXD )
        IF ( LDSPOA( VIDX ) .AND. J .GT. 0 ) THEN

            TMAT( S,VIDX,HIDX ) = EMACH( I )

        !.......  Adjust for year-normal data

        ELSE IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN

            IF ( IMON .GT. 0 ) THEN
                FAC = MONFAC( MON,IMON ) * DOMFAC( DOM,MON,IDOM )
            ELSE
                FAC = MONFAC( MON,IMON ) * WEKFAC( DAY,IWEK )
            END IF
            TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

        !.......  Adjust for week-normal data assuming whole week normalizer

        ELSE IF ( MOD( TPFLAG( S ), WTPRFAC ) .EQ. 0 ) THEN

            IF ( IMON .GT. 0 ) THEN
                FAC = YRFAC * DOMFAC( DOM,MON,IDOM )
            ELSE
                FAC = YRFAC * WEKFAC( DAY,IWEK )
            END IF
            TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )


        !.......  Adjust for week-normal data assuming week-days normalizer

        ELSE IF ( MOD( TPFLAG( S ), WDTPFAC ) .EQ. 0 ) THEN

            IF ( IMON .GT. 0 ) THEN
                FAC = YRFAC * DOMFAC( DOM,MON,IDOM )
            ELSE
                FAC = YRFAC * XWKFAC( DAY,IWEK )
            END IF
            TMAT( S,VIDX,HIDX ) = FAC * EMACH( I )

        ELSE

            TMAT( S,VIDX,HIDX ) = YRFAC * EMACH( I )

        END IF

    END DO              ! end loop hour-specific sources

    RETURN

END SUBROUTINE UPDTMAT

