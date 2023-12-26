
LOGICAL FUNCTION CHKEMEPI( JDATE, JTIME, TIMESTEP, JRUNLEN )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine checks the episode parameters to make sure they
    !      are consistent with what SMOKE expects.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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
    !***************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER         JDATE           ! Julian date (DDDYYYY) (in)
    INTEGER         JTIME           ! time (HHMMSS) (in)
    INTEGER         TIMESTEP        ! duration of time step (HHMMSS) (in)
    INTEGER         JRUNLEN         ! episode duration (HHMMSS) (in)

    !.......   Other local variables
    REAL            REMNDR      ! REMNDRder for checking time units

    CHARACTER(300)  MESG       ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'CHKEMEPI'     ! program name

    !***********************************************************************
    !   begin body of function CHKEMEPI

    CHKEMEPI = .TRUE.

    IF( JDATE .LT. 1970001 .OR. JDATE .GT. 2100365 ) THEN

        CHKEMEPI = .FALSE.
        MESG = 'Date is out of acceptable modeling range:' //   &
               CRLF() // BLANK16 // '1970001 to 2100365.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

    ENDIF

    IF( JTIME .LT. 0 .OR. JTIME .GT. 235959 ) THEN

        CHKEMEPI = .FALSE.
        MESG = 'Time is out of acceptable modeling range:' //   &
               CRLF() // BLANK16 // '0 to 235959.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

    ENDIF

    REMNDR = MOD( JTIME, 10000 )

    IF( REMNDR .GT. 0 ) THEN

        CHKEMEPI = .FALSE.
        MESG = 'Time may only be in units of hours.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

    ENDIF

    IF( TIMESTEP .NE. 10000 ) THEN

        CHKEMEPI = .FALSE.
        WRITE( MESG,94010 )                 &
           'Timestep may only be one hour, but it is', TIMESTEP, '(HHMMSS).'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

    ENDIF

    REMNDR = MOD( JRUNLEN, 10000 )

    IF( REMNDR .GT. 0 ) THEN

        CHKEMEPI = .FALSE.
        MESG = 'Run duration may only be in units of hours.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

    ENDIF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************
    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END FUNCTION CHKEMEPI

