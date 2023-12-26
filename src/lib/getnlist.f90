
INTEGER FUNCTION GETNLIST( ILENGTH, STRING )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function counts the number of free-formatted strings in a list
    !      of string that may or may not have quotes.  This is used to help
    !      when a string is available for reading a series of values, but
    !      no indication is available for the number of entries to read.  It
    !      accounts for comments appended with a "    !" and skips blank lines.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 1/99
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

    USE M3UTILIO

    IMPLICIT NONE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER       ILENGTH         !  length of string
    CHARACTER(*)  STRING          !  description of source category

    !.......   Local parameters
    INTEGER,       PARAMETER :: NDELIM = 4
    CHARACTER,     PARAMETER :: DELIMLST( NDELIM ) = (/ ',', ' ', ';', '	' /)
    CHARACTER,     PARAMETER :: DOUBLEQ = '"'
    CHARACTER,     PARAMETER :: SINGLEQ = "'"
    CHARACTER,     PARAMETER :: PERIOD  = '.'
    CHARACTER(16), PARAMETER :: PROGNAME = 'GETNLIST'     ! program name

    !.......   Array of 1-char strings for processing
    CHARACTER     ARRSTR( 5120 )      ! 256 * 20

    !.......  Arrays for sorting non-delimiters on a per-machine basis
    INTEGER            NDINDX  ( NDELIM )
    CHARACTER, SAVE :: DELIMSRT( NDELIM )

    !.......   Other local variables
    INTEGER         I, J, L, L1, L2      !  counters and indices
    INTEGER         IXP                  !  index to non-delimeters
    INTEGER      :: NCNT                 !  count of fields

    LOGICAL      :: ALPHA                !  true when within alpha-numeric
    LOGICAL      :: DELIM                !  true when within or past delimiter
    LOGICAL      :: FIRSTIME = .TRUE.    !  true first time routine is called
    LOGICAL      :: PREVDELIM = .TRUE.     !  true when last char was a delim
    LOGICAL      :: NUMBER               !  true when within number in string
    LOGICAL      :: QUOTED               !  true when within quotes in string
    LOGICAL      :: THISNMBR             !  true when current iteration is numbr

    CHARACTER       CBUF                 !  temporary buffer
    CHARACTER       QUOTVAL              !  value of starting quote

    CHARACTER(300)  MESG                 ! message buffer

    !***********************************************************************
    !   begin body of function GETNLIST

    !.....  The first time the routine is called, sort the list of delimiters
    IF( FIRSTIME ) THEN
        DO I = 1, NDELIM
            NDINDX( I ) = I
        END DO

        CALL SORTIC( NDELIM, NDINDX, DELIMLST )

        DO I = 1, NDELIM
            J = NDINDX( I )
            DELIMSRT( I ) = DELIMLST( J )
        END DO

        FIRSTIME = .FALSE.

    END IF

    L2 = LEN_TRIM( STRING )

    !.....  Check for comments, and use to set the end of the line
    L = INDEX( STRING( 1:L2 ), '    !' )

    IF( L .LE. 0 ) THEN
        L = L2
    ELSE
        L = L - 1
    END IF

    !.....  Skip blank lines
    IF( L .EQ. 0 ) THEN
        GETNLIST = 0
        RETURN
    END IF

    !.....  Initialize count and flags
    NCNT    = 0
    ALPHA   = .FALSE.
    DELIM   = .TRUE.
    NUMBER  = .FALSE.
    QUOTED  = .FALSE.

    !.....  Process LINE 1-character at a time
    DO I = 1, L

        CBUF = STRING( I:I )

        !.....  Look for character in delimiters
        IXP = FINDC( CBUF, NDELIM, DELIMSRT )

        !.....  Evaluate the current character for number or not
        THISNMBR = ( CBUF .GE. '0' .AND. CBUF .LE. '9' )

        !.....  Waiting for next field...
        IF( DELIM ) THEN

            NUMBER = THISNMBR
            ALPHA  = ( .NOT. NUMBER .AND. IXP .LE. 0 )

            IF( CBUF .EQ. SINGLEQ ) THEN
                QUOTED  = .TRUE.
                DELIM   = .FALSE.
                QUOTVAL = SINGLEQ
                PREVDELIM = .FALSE.
                L1     = I + 1
                NCNT    = NCNT + 1

            ELSE IF( CBUF .EQ. DOUBLEQ ) THEN
                QUOTED  = .TRUE.
                DELIM   = .FALSE.
                QUOTVAL = DOUBLEQ
                PREVDELIM = .FALSE.
                L1      = I + 1
                NCNT    = NCNT + 1

            ELSE IF( ALPHA ) THEN
                DELIM = .FALSE.
                PREVDELIM = .FALSE.
                L1    = I
                NCNT  = NCNT + 1

            ELSE IF( NUMBER ) THEN
                DELIM  = .FALSE.
                PREVDELIM = .FALSE.
                L1     = I
                NCNT   = NCNT + 1

            !.......  If another delimeter, then another field, but last
            !                 field was blank UNLESS delim is a space
            ELSE IF( CBUF .NE. DELIMLST( 2 ) ) THEN

                IF( PREVDELIM ) THEN
                    NCNT = NCNT + 1
                ELSE
                    PREVDELIM = .TRUE.
                END IF

            END IF              ! Else its a space delimiter

            !.....  In a quoted field, skip everything unless it is an end quote
        ELSE IF( QUOTED ) THEN

            IF( CBUF .EQ. QUOTVAL ) THEN
                QUOTED  = .FALSE.
                DELIM   = .TRUE.
                PREVDELIM = .FALSE.
                L2      = I - 1

            END IF

        !.....  If start of field was a number, but adjacent character is not
        !               a delimiter, then turn field into an alpha
        ELSE IF( NUMBER .AND. .NOT. THISNMBR .AND. IXP .LE. 0 ) THEN
            ALPHA  = .TRUE.
            NUMBER = .FALSE.

        !.....  If start of field was a number or alpha, and this is a
        !               delimiter, then end of number has been reached
        ELSE IF( IXP .GT. 0 ) THEN
            ALPHA = .FALSE.
            NUMBER = .FALSE.
            DELIM  = .TRUE.
            PREVDELIM = .TRUE.
            L2     = I - 1

        END IF

    END DO

    GETNLIST = NCNT

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END FUNCTION GETNLIST
