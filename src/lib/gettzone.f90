
INTEGER FUNCTION GETTZONE( CFIP )

    !***********************************************************************
    !  function body starts at line 76
    !
    !  DESCRIPTION:
    !     Returns the time zone of a FIPs code based on the input tables
    !
    !  PRECONDITIONS REQUIRED:
    !     FIP set to state and county FIPS code of interest
    !     Memory allocated for arrays in argument list and arrays populated
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !     Subroutines: I/O API subroutines
    !     Functions: I/O API functions
    !
    !  REVISION  HISTORY:
    !       Created 10/98 by M. Houyoux
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
    !****************************************************************************

    USE M3UTILIO

    !.....  MODULES for public variables
    !.....  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, CNTYTZON

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   ARGUMENTS and their descriptions:

    CHARACTER(FIPLEN3), INTENT (IN) :: CFIP     !  input FIPS code to assign zone

    !.......   LOCAL VARIABLES their descriptions:

    INTEGER          K                ! indice
    INTEGER          IOS              ! i/o status
    INTEGER, SAVE :: TZONE0           ! default time zone

    LOGICAL, SAVE :: FIRSTIME = .TRUE.

    CHARACTER        BUFFER           ! ASCII LINE from X-ref file
    CHARACTER(300)   MESG             ! Message buffer

    CHARACTER(16), PARAMETER :: PNAME = 'GETTZONE'     ! Program name

    !***********************************************************************
    !   begin body of subroutine GETTZONE

    !.....  For the first time the routine is called, retrieve the default
    !           time zone from the environment (or use the ultimate default of 5)
    IF( FIRSTIME ) THEN

        MESG = 'Default time zone for sources'
        TZONE0 = ENVINT( 'SMK_DEFAULT_TZONE', MESG, 5, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SMK_DEFAULT_TZONE"', 2 )
        END IF
        FIRSTIME = .FALSE.

    END IF

    !.....  Search for FIPS code in county-specific table
    K = FINDC( CFIP, NCOUNTY, CNTYCOD )

    IF ( K .GT. 0 ) THEN
        GETTZONE = CNTYTZON( K )

    !.....  Apply default
    ELSE
        GETTZONE = TZONE0

        WRITE( MESG,94040 )                                         &
                  'WARNING: Applying default time zone of', TZONE0, &
                  'to country, state, and county code ' // CFIP
        CALL M3MESG( MESG )

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************C.......   Internal buffering formats.....94xxx

    !.......   Internal buffering formats.... 94xxx

94040 FORMAT( A, 1X, I2.2, 1X, A )

END FUNCTION GETTZONE
