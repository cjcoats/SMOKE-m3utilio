
SUBROUTINE PARSLINE( LINE, N, SEGMENT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine separates a "list-formatted" line of strings in which
    !      the segments may or may not have quotes.  Although fortran requires
    !      the quotes for true list-formatting, this subroutine can be used when
    !      the quotes are only present to enclose a character (such as space, comma,
    !      or semi-colon) that would otherwise be a delimiter.  If an "    !" is
    !      encountered, everything after it is treated as a comment.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 3/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !****************************************************************************/
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

    !......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: LINE             ! character string to parse
    INTEGER     , INTENT (IN) :: N                ! maximum array length
    CHARACTER(*), INTENT(OUT) :: SEGMENT( N )     ! parsed string

    !......   Local parameters
    INTEGER,   PARAMETER :: NDELIM = 4
    CHARACTER, PARAMETER :: DELIMLST( NDELIM ) = (/ ',', ' ', ';', '	' /)
    CHARACTER, PARAMETER :: DOUBLEQ = '"'
    CHARACTER, PARAMETER :: SINGLEQ = "'"
    CHARACTER, PARAMETER :: PERIOD  = '.'

    CHARACTER(16), PARAMETER :: PROGNAME = 'PARSLINE'     ! program name

    !......   Array of 1-char strings for processing
    CHARACTER   ARRSTR( 5120 )      ! 256 * 20

    !......  Arrays for sorting non-delimiters on a per-machine basis
    INTEGER            NDINDX  ( NDELIM )
    CHARACTER, SAVE :: DELIMSRT( NDELIM )

    !......   Other local variables
    INTEGER         I, J, K, L, L1, L2      !  counters and indices
    INTEGER         IXP                  !  index to non-delimeters
    INTEGER      :: NCNT                 !  count of fields

    LOGICAL      :: ALPHA                !  true when within alpha-numeric
    LOGICAL      :: DELIM                !  true when within or past delimiter
    LOGICAL, SAVE:: FIRSTIME = .TRUE.    !  true first time routine is called
    LOGICAL      :: PREVDELIM = .TRUE.   !  true when last char was a delim
    LOGICAL      :: NOSPACDLIM = .FALSE. !  true when encountered a non-space delimiter since the previous field
    LOGICAL      :: NUMBER               !  true when within number in string
    LOGICAL      :: QUOTED               !  true when within quotes in string
    LOGICAL      :: THISNMBR             !  true when current iteration is numbr

    CHARACTER    :: QUOTVAL = ' '        !  value of starting quote
    CHARACTER       CBUF                 !  temporary buffer

    CHARACTER(300)  MESG                 ! message buffer

    !***********************************************************************
    !   begin body of subroutine PARSLINE

    !......  The first time the routine is called, sort the list of delimiters
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

    L2 = LEN_TRIM( LINE )

    !......  Check for comments, and use to set the end of the line
    L = INDEX( LINE( 1:L2 ), '!' )

    IF( L .LE. 0 ) THEN
        L = L2
    ELSE
        L = L - 1
    END IF

    !......  Skip blank lines
    IF( L .EQ. 0 ) RETURN

    !......  Initialize count, flags, and segments (npte, initializing in
    !           the variable definitions is insufficient)
    NCNT    = 0
    SEGMENT = ' '     ! array
    ALPHA   = .FALSE.
    DELIM   = .TRUE.
    NUMBER  = .FALSE.
    QUOTED  = .FALSE.

    !......  Process LINE 1-character at a time
    DO I = 1, L

        CBUF = LINE( I:I )

        !......  Look for character in delimiters
        IXP = FINDC( CBUF, NDELIM, DELIMSRT )

        !......  Evaluate the current character for number or not
        THISNMBR = ( CBUF .GE. '0' .AND. CBUF .LE. '9' )

        !......  Waiting for next field...
        IF( DELIM ) THEN

            NUMBER = THISNMBR
            ALPHA  = ( .NOT. NUMBER .AND. IXP .LE. 0 )

            IF( CBUF .EQ. SINGLEQ ) THEN
                QUOTED  = .TRUE.
                DELIM   = .FALSE.
                QUOTVAL = SINGLEQ
                PREVDELIM = .FALSE.
                NOSPACDLIM = .FALSE.
                L1     = I + 1
                NCNT    = NCNT + 1
                K = 1

            ELSE IF( CBUF .EQ. DOUBLEQ ) THEN
                QUOTED  = .TRUE.
                DELIM   = .FALSE.
                QUOTVAL = DOUBLEQ
                PREVDELIM = .FALSE.
                NOSPACDLIM = .FALSE.
                L1      = I + 1
                NCNT    = NCNT + 1
                K = 1

            ELSE IF( ALPHA ) THEN
                DELIM = .FALSE.
                PREVDELIM = .FALSE.
                NOSPACDLIM = .FALSE.
                L1    = I
                NCNT  = NCNT + 1
                K = 1

            ELSE IF( NUMBER ) THEN
                DELIM  = .FALSE.
                PREVDELIM = .FALSE.
                NOSPACDLIM = .FALSE.
                L1     = I
                NCNT   = NCNT + 1
                K = 1

            !......  If hit another non-space delimiter without having
            !        hit another field with contents, then iterate the
            !        field count to create a blank space.
            ELSE IF( CBUF .NE. DELIMLST( 2 ) .AND.&
                     NOSPACDLIM                   ) THEN
                NCNT = NCNT + 1
                PREVDELIM = .TRUE.
                K = 1

            !......  If first col of line is a delimiter, fill a blank for the first col
            ELSE IF( I == 1 .AND. IXP .GT. 0 .AND.&
                     CBUF .NE. DELIMLST( 2 )     ) THEN
                ALPHA = .FALSE.
                NUMBER = .FALSE.
                DELIM  = .TRUE.
                PREVDELIM = .TRUE.
                NOSPACDLIM = .TRUE.

                NCNT = NCNT + 1
                SEGMENT( NCNT ) = ' '

            END IF              ! Else its a space delimiter

            !......  In a quoted field, skip everything unless it is an end quote
        ELSE IF( QUOTED ) THEN

            IF( CBUF .EQ. QUOTVAL ) THEN
                QUOTED  = .FALSE.
                PREVDELIM = .FALSE.
                K = 2
            END IF

            !......  If start of field was a number, but adjacent character is not
            !        a delimiter, then turn field into an alpha
        ELSE IF( NUMBER .AND. .NOT. THISNMBR .AND. IXP .LE. 0 ) THEN
            ALPHA  = .TRUE.
            NUMBER = .FALSE.
            K = 1

            !......  If start of field was a number or alpha, and this char is a
            !        delimiter, then end of the field has been reached
        ELSE IF( IXP .GT. 0 ) THEN
            ALPHA = .FALSE.
            NUMBER = .FALSE.
            DELIM  = .TRUE.
            PREVDELIM = .TRUE.
            IF( CBUF .NE. DELIMLST( 2 ) ) NOSPACDLIM = .TRUE.
            L2     = I - K

            CALL STORE_SEGMENT

        END IF

    END DO

    !......  Store final segment
    IF( CBUF .EQ. QUOTVAL ) L = L - 1
    L2 = L

    IF( IXP .LE. 0 ) CALL STORE_SEGMENT

    RETURN

    !******************  FORMAT  STATEMENTS   ************************************

    !......   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !......  This subprogram stores the segment from the input string
    SUBROUTINE STORE_SEGMENT

        IF( NCNT .LE. N ) THEN

            SEGMENT( NCNT ) = ADJUSTL( LINE( L1:L2 ) )

        ELSE

            MESG = 'ERROR: Overflow prevented while parsing line '
            CALL M3MSG2( MESG )
            MESG = 'First 200 characters of line contents are:'
            CALL M3MSG2( MESG )
            MESG = LINE( 1:200 )
            CALL M3MSG2( MESG )

            MESG = 'Formatting problem.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END SUBROUTINE STORE_SEGMENT

END SUBROUTINE PARSLINE
