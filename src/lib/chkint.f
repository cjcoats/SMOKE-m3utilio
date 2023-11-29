
        LOGICAL FUNCTION CHKINT( STRING )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns true if the string appears to be an integer
C      and false otherwise
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C**************************************************************************
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

        IMPLICIT NONE


C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT(IN   ) :: STRING   ! input character string

C...........   Other local variables
        INTEGER         K, L  ! counters and indices

        LOGICAL         SPACFLAG  ! true if already encountered a space in string
        LOGICAL         NEGVFLAG  ! true if already encountered a '-' in string

        CHARACTER       CBUF
        CHARACTER(256)  BUFFER

C***********************************************************************
C   begin body of function CHKINT

        CHKINT = .TRUE.

        BUFFER = ADJUSTL( STRING )
        L = LEN_TRIM( BUFFER )

        SPACFLAG = .FALSE.
        NEGVFLAG = .FALSE.
        DO K = 1, L

            CBUF = BUFFER( K:K )

            IF( CBUF .GT. '9' .OR.
     &        ( CBUF .LT. '0' .AND.
     &          CBUF .NE. ' ' .AND.
     &          CBUF .NE. '-'       ) .OR.
     &        ( SPACFLAG .AND. CBUF .EQ. ' ' ) .OR.
     &        ( NEGVFLAG .AND. CBUF .EQ. '-' )           ) THEN

                CHKINT = .FALSE.
                RETURN

            ENDIF

            IF( CBUF .EQ. ' ' ) SPACFLAG = .TRUE.
            IF( CBUF .EQ. '-' ) NEGVFLAG = .TRUE.

        ENDDO


        RETURN

        END FUNCTION CHKINT

