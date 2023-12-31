
        SUBROUTINE POLMESG( NLIST, NAMES )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine writes out a message stating that the pollutants or
C      emission types in the argument list are being processed.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C      Subroutines: I/O API subroutines
C
C  REVISION  HISTORY:
C       Created 3/99 by M. Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C************************************************************************
C
C Project Title: EDSS Tools Library
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
C***************************************************************************
        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDES
        INCLUDE 'IOCNST3.EXT'   !  emissions constant parameters

C.........  SUBROUTINE ARGUMENTS
        INTEGER     , INTENT (IN) :: NLIST          !  no. of pols or emis types
        CHARACTER(*), INTENT (IN) :: NAMES( NLIST ) !  pollutant names

C...........   Other local variables
        INTEGER       I, J, L0, L1, L2
        INTEGER       LCNT              ! length count

        CHARACTER(5120) :: MESG         !  message buffer
        CHARACTER(20)   :: SPACE = ' '

        CHARACTER(16) :: PROGNAME = 'POLMESG' !  program name

C***********************************************************************
C   begin body of subroutine POLMESG

C.........  Set up initial message.
        MESG = 'Processing data for:'

        IF( NAMES( 1 ) .NE. ' ' ) THEN
            MESG = TRIM( MESG ) // '"' // TRIM( NAMES( 1 ) ) // '"'
        END IF

C.........  Initialize length of initial message
        LCNT = L0 + L1 + 4
        DO I = 2, NLIST

            IF( NAMES( I ) .EQ. ' ' ) CYCLE

            LCNT = LCNT + L2 + 4

            IF( LCNT .GT. 80 ) THEN
                LCNT = L0 + L2 + 4
                MESG = TRIM( MESG ) // ',' // CRLF() // BLANK5 //
     &                 SPACE // '"' // TRIM( NAMES( I ) ) // '"'
            ELSE
                MESG = MESG( 1:L1 )// ', "'// TRIM( NAMES( I ) )// '"'
            END IF

        END DO

        CALL M3MSG2( MESG )

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93030   FORMAT( I8.8 )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

        END SUBROUTINE POLMESG
