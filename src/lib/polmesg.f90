
SUBROUTINE POLMESG( NLIST, NAMES )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine writes out a message stating that the pollutants or
    !      emission types in the argument list are being processed.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutines
    !
    !  REVISION  HISTORY:
    !       Created 3/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !************************************************************************
    !
    ! Project Title: EDSS Tools Library
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

    !......   INCLUDES
    INCLUDE 'IOCNST3.h90'       !  emissions constant parameters

    !......  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NLIST              !  no. of pols or emis types
    CHARACTER(*), INTENT (IN) :: NAMES( NLIST )     !  pollutant names

    !......   Other local variables
    INTEGER       I, J, L0, L1, L2
    INTEGER       LCNT                  ! length count

    CHARACTER(5120) :: MESG             !  message buffer

    CHARACTER(20), PARAMETER :: SPACE = ' '

    !***********************************************************************
    !   begin body of subroutine POLMESG

    !......  Set up initial message.
    MESG = 'Processing data for:'

    IF( NAMES( 1 ) .NE. ' ' ) THEN
        MESG = TRIM( MESG ) // '"' // TRIM( NAMES( 1 ) ) // '"'
    END IF

    !......  Initialize length of initial message
    LCNT = L0 + L1 + 4
    DO I = 2, NLIST

        IF( NAMES( I ) .EQ. ' ' ) CYCLE

        LCNT = LCNT + L2 + 4

        IF( LCNT .GT. 80 ) THEN
            LCNT = L0 + L2 + 4
            MESG = TRIM( MESG ) // ',' // CRLF() // BLANK5 //   &
                   SPACE // '"' // TRIM( NAMES( I ) ) // '"'
        ELSE
            MESG = TRIM( MESG ) // ', "'// TRIM( NAMES( I ) )// '"'
        END IF

    END DO

    CALL M3MSG2( MESG )

    RETURN

END SUBROUTINE POLMESG
