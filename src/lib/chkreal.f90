
LOGICAL FUNCTION CHKREAL( STRING )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns true if the string appears to be a real
    !      and false otherwise
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

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN OUT) :: STRING       ! input character string

    !.......   Other local variables
    INTEGER         K, L      ! counters and indices

    LOGICAL         EXPOFLAG      ! true if exponential
    LOGICAL         SPACFLAG      ! true if already encountered a space in string
    LOGICAL         PERDFLAG      ! true if already encountered a period in string
    LOGICAL         NEGVFLAG      ! true if already encountered a '-' in string

    CHARACTER       CBUF
    CHARACTER(256)  BUFFER

    !***********************************************************************
    !   begin body of function CHKREAL

    CHKREAL = .TRUE.

    BUFFER = ADJUSTL( STRING )
    L = LEN_TRIM( BUFFER )

    CALL UPCASE( BUFFER )

    EXPOFLAG = .FALSE.
    SPACFLAG = .FALSE.
    PERDFLAG = .FALSE.
    NEGVFLAG = .FALSE.
    DO K = 1, L

        CBUF = BUFFER( K:K )

        IF( ( CBUF .GT. '9' .AND. CBUF .NE. 'E' ) .OR.  &
            ( CBUF .LT. '0' .AND.   &
              CBUF .NE. ' ' .AND.   &
              CBUF .NE. '.' .AND.   &
              CBUF .NE. '-' .AND. CBUF .NE. '+' ) .OR.  &
            ( ( SPACFLAG .OR. PERDFLAG .OR. NEGVFLAG ) .AND.&
                CBUF .EQ. ' ' ) .OR.                    &
            ( NEGVFLAG .AND. CBUF .EQ. '-' ) .OR.       &
            ( PERDFLAG .AND. CBUF .EQ. '.' ) .OR.       &
            ( EXPOFLAG .AND. CBUF .EQ. 'E' )      ) THEN

            CHKREAL = .FALSE.
            RETURN

        END IF

        IF( CBUF .EQ. ' ' ) SPACFLAG = .TRUE.
        IF( CBUF .EQ. '.' ) PERDFLAG = .TRUE.
        IF( CBUF .EQ. '-' ) NEGVFLAG = .TRUE.
        IF( CBUF .EQ. 'E' ) THEN
            EXPOFLAG = .TRUE.
            NEGVFLAG = .FALSE.      ! Reinitialize after exponent
        END IF
    END DO

    IF( BUFFER .EQ. '.' ) STRING = '0'

    RETURN

END FUNCTION CHKREAL

