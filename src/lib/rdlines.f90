
SUBROUTINE RDLINES( FDEV, DESCRIPT, NLINES, LINES )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the lines of an ASCII file to an array of strings
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
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

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT(IN   ) :: FDEV                !  file unit number
    CHARACTER(*), INTENT(IN   ) :: DESCRIPT            !  file description
    INTEGER     , INTENT(IN   ) :: NLINES              !  number of lines in file
    CHARACTER(*), INTENT(  OUT) :: LINES( NLINES )     !  ASCII lines in file

    !.......   Other local variables
    INTEGER         IOS         !  i/o status
    INTEGER         IREC        !  line counter
    INTEGER         L, LSAV, L2       !  length indices
    INTEGER         N           !  record counter

    CHARACTER(300)  LINE        !  line buffer
    CHARACTER(300)  MESG        !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDLINES'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDLINES

    L = LEN( LINES( 1 ) )

    IREC = 0
    LSAV = 0
    N    = 0
11  CONTINUE

    READ( FDEV, '(A)', END=22, IOSTAT=IOS ) LINE
    IREC = IREC + 1

    !.......  Skip blank lines
    IF( LINE .EQ. ' ' ) GO TO 11

    L2 = LEN_TRIM( LINE )

    IF( IOS .GT. 0 ) THEN
        WRITE( MESG, 94010 ) 'Error', IOS,  'reading ' //   &
               TRIM( DESCRIPT ) // ' file at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !....... Keep track of width that is larger than that allocated
    ELSE IF( L2 .GT. L ) THEN
        IF( L2 .GT. LSAV ) LSAV = L2
        GO TO 11

    !....... Skip blank lines
    ELSE IF( L2 .EQ. 0 ) THEN
        GO TO 11

    END IF

    N = N + 1
    IF( N .LE. NLINES ) THEN
        LINES( N ) = LINE
    END IF

    GO TO  11

22  CONTINUE            !  exit from loop reading FDEV

    IF( N .GT. NLINES ) THEN
        WRITE( MESG,94010 ) 'WARNING: ' // TRIM( DESCRIPT ) //      &
               CRLF() // BLANK10 // 'file only read for first ',    &
               NLINES, ' lines of ', N, ' total lines.'
        CALL M3MSG2( MESG )
    END IF

    IF( LSAV .GT. 0 ) THEN
        WRITE( MESG,94010 ) 'ERROR: ' // TRIM( DESCRIPT ) //        &
               CRLF() // BLANK10 // 'file line width is ',          &
               LSAV, ' but allocated string length is ', L
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Rewind file
    REWIND( FDEV )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDLINES
