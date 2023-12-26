
LOGICAL FUNCTION CHKINT( STRING )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns true if the string appears to be an integer
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
    CHARACTER(*), INTENT(IN   ) :: STRING       ! input character string

    !.......   Other local variables
    INTEGER         K, L      ! counters and indices

    LOGICAL         SPACFLAG      ! true if already encountered a space in string
    LOGICAL         NEGVFLAG      ! true if already encountered a '-' in string

    CHARACTER       CBUF
    CHARACTER(256)  BUFFER

    !***********************************************************************
    !   begin body of function CHKINT

    CHKINT = .TRUE.

    BUFFER = ADJUSTL( STRING )
    L = LEN_TRIM( BUFFER )

    SPACFLAG = .FALSE.
    NEGVFLAG = .FALSE.
    DO K = 1, L

        CBUF = BUFFER( K:K )

        IF( CBUF .GT. '9'                                             .OR. &
            ( CBUF .LT. '0' .AND. CBUF .NE. ' ' .AND. CBUF .NE. '-' ) .OR. &
            ( SPACFLAG .AND. CBUF .EQ. ' ' )                          .OR. &
            ( NEGVFLAG .AND. CBUF .EQ. '-' )) THEN

            CHKINT = .FALSE.
            RETURN

        ENDIF

        IF( CBUF .EQ. ' ' ) SPACFLAG = .TRUE.
        IF( CBUF .EQ. '-' ) NEGVFLAG = .TRUE.

    ENDDO


    RETURN

END FUNCTION CHKINT

