
CHARACTER(*) FUNCTION GETCFDSC( FILEINFO, KEY, REQUIRED )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !     Retreives a character string from the FDESC array
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
    !
    !  REVISION  HISTORY:
    !       prototype 1/99 by M Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
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
    !****************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......   Include files
    INCLUDE 'IOCNST3.h90'       !  I/O API constant parameters

    !.......   ARGUMENTS and their descriptions:

    CHARACTER(*), INTENT (IN) :: FILEINFO( * )     ! Array of file information
    CHARACTER(*), INTENT (IN) :: KEY               ! Key to find in FILEINFO
    LOGICAL     , INTENT (IN) :: REQUIRED          ! true: key must be found

    !.......   LOCAL VARIABLES their descriptions:

    INTEGER       L1, L2           ! length of strings
    INTEGER       I                ! loop index
    INTEGER       K                ! description string position of key
    INTEGER       LENGTH           ! length of function

    CHARACTER(300) CVAL            ! temporary output buffer
    CHARACTER(300) BUFFER          ! key buffer
    CHARACTER(300) MESG            ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETCFDSC'        ! Program name

    !***********************************************************************
    !   begin body of subroutine GETCFDSC

    !.....  Ensure left-justified keyword
    BUFFER = ADJUSTL( KEY )

    L1 = LEN_TRIM( BUFFER )

    !.....  Get length of function
    LENGTH = LEN( GETCFDSC )

    DO I = 1, MXDESC3

        K = INDEX( FILEINFO( I ), BUFFER( 1:L1 ) )

        IF( K .GT. 0 ) THEN

            L2 = MAX( K+L1, LEN_TRIM( FILEINFO( I ) ) )
            CVAL = FILEINFO( I )( K+L1:L2 )

            L2 = LEN_TRIM( CVAL )

            IF( L2 .GT. LENGTH ) THEN
                WRITE( MESG,94010 )                                         &
                       'INTERNAL WARNING: Length of string used ' //        &
                       'to call function "' // TRIM( PROGNAME ) // '"' //   &
                       CRLF() // BLANK16 //                                 &
                       'is not long enough. Needs to be at least',          &
                       LENGTH, '. Trimming FDESC3D entry.'
                CALL M3MSG2( MESG )

            END IF

            GETCFDSC = ADJUSTL( CVAL( 1:LENGTH ) )
            RETURN

        END IF

    END DO

    !.....  If we get here, then key was not found in FDESC, so if it was
    !           required, then abort.

    IF( REQUIRED ) THEN
        MESG = 'FDESC3D packet "' // KEY( 1:L1 ) // '" was not found in file'
        CALL M3EXIT( MESG, 0, 0, PROGNAME, 2 )

    ELSE
        GETCFDSC = ' '
        RETURN

    END IF

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.....94xxx

94010 FORMAT( 10( A, :, I6, :, 2X ) )

END FUNCTION GETCFDSC
