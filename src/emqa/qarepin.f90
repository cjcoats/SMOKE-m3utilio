
SUBROUTINE QAREPIN( RCNT, IOS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      The QAREPIN routine compares the settings for the current report and
    !      gives errors and warnings if problems are found.
    !
    !  PRECONDITIONS REQUIRED:
    !    REPCONFIG file is opened
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 3/2002 by M Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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
    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: RPT_, SDATE, STIME, EDATE, ETIME, NSTEPS,   &
                        TSTEP, RPTNSTEP, ALLRPT, AFLAG

    !.......  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: LSTCYPOP, STCYPOPYR

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: BYEAR

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN)  :: RCNT           ! report number
    INTEGER, INTENT (OUT) :: IOS            ! error status

    !.......   Other local variables
    INTEGER          I
    INTEGER          JDATE                  ! Julian date
    INTEGER          JTIME                  ! time

    LOGICAL          EFLAG                  !  true: error found

    CHARACTER(300)         MESG             !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'QAREPIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine QAREPIN

    !.......  Initialize output error status
    EFLAG   = .FALSE.
    IOS = 0

    !.......  Other checks that could be added to this routine:
    !        N: Check units (and that multiple ones aren't specified unless this
    !        N:     is supported)

    !.......  Checks when population normalization is used...
    IF( RPT_%NORMPOP ) THEN

        !.......  Ensure that report is by region of some sort
        IF ( .NOT. RPT_%BYCNRY .AND.    &
             .NOT. RPT_%BYSTAT .AND.    &
             .NOT. RPT_%BYCNTY       ) THEN

            WRITE( MESG,94010 ) BLANK5 //                               &
                   'WARNING: Invalid request for population '//         &
                   'normalization without requesting' // CRLF()//       &
                   BLANK16 // 'by country, by state, or by county. '    &
                   // 'Dropping population normalization.'
            CALL M3MSG2( MESG )

            RPT_%NORMPOP = .FALSE.
            ALLRPT( RCNT )%NORMPOP = .FALSE.

        END IF

        !.......  Ensure population is present in COSTCY file.
        IF( .NOT. LSTCYPOP ) THEN
            WRITE( MESG,94010 ) BLANK5 //                           &
                   'WARNING: Invalid request for population '//     &
                   'when population is not present' // CRLF()//     &
                   BLANK16 // 'in COSTCY file. '                    &
                   // 'Dropping population normalization.'

            CALL M3MSG2( MESG )

            RPT_%NORMPOP = .FALSE.
            ALLRPT( RCNT )%NORMPOP = .FALSE.

        !.......  Compare population year to inventory year.
        ELSE IF( STCYPOPYR .NE. BYEAR ) THEN
            WRITE( MESG,94010 ) BLANK5 //                           &
                   'NOTE: Population year ', STCYPOPYR,             &
                   'is inconsistent with inventory base year',      &
                   BYEAR
            CALL M3MSG2( MESG )

        END IF

    END IF

    !.......  Set ending date and time and number of time steps for report
    !.......  When using hourly inputs but reporting daily totals
    IF( RPT_%USEHOUR .OR. AFLAG .AND. RPT_%BYHOUR ) THEN
        JDATE = SDATE
        JTIME = STIME

        !.......  Find ending time
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

        !.......  Compare data end time with output end time
        I = SECSDIFF( EDATE, ETIME, EDATE, RPT_%OUTTIME )

        !.......  If reporting time is after data ending time, reset the no.
        !         of time steps so that the reporting ends on the previous day
        IF( I .GT. 0 ) THEN
                ! Workaround - NEXTIME does not properly subtract 24 hours so have
                !   to use two steps

            IF( .NOT. RPT_%BYHOUR ) THEN
                CALL NEXTIME( EDATE, ETIME, -23 * TSTEP )
                CALL NEXTIME( EDATE, ETIME, -1 * TSTEP )
                ETIME = RPT_%OUTTIME
            END IF

            I =  SECSDIFF( SDATE, STIME, EDATE, ETIME )
            RPTNSTEP = I / 3600 + 1

        !.......  If reporting time is before data ending time, reset the
        !         number of time steps so that the reporting ends on the
        !         reporting hour
        !.......  Also set reporting time steps for reporting time matches
        !         ending time.
        ELSE IF( I .LE. 0 ) THEN
            IF( .NOT. RPT_%BYHOUR ) RPTNSTEP = NSTEPS + I / 3600

        END IF

        !.......  Print message if time steps have changed
        IF( I .NE. 0 .AND. .NOT. RPT_%BYHOUR ) THEN
            WRITE( MESG,94010 ) BLANK5 //                               &
                   'WARNING: Resetting number of time steps for ' //    &
                   'report to ', RPTNSTEP, CRLF() // BLANK16 //         &
                   'to make output hour consistent with reporting time.'
            CALL M3MSG2( MESG )
        END IF

    !.......  When not using hourly inputs
    ELSE
        RPTNSTEP = 1
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE QAREPIN

