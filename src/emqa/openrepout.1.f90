
SUBROUTINE OPENREPOUT( FILNAM, FDEV, MDEV )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!    The OPENREPOUT routine open output file FILNAM, which can be a physical
!    or logical file name.  Returns FDEV to calling program.
!
!  PRECONDITIONS REQUIRED:
!    FILNAM is defined as a physical or logical file name.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Revised 7/2003 by A. Holland
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!***********************************************************************
    USE M3UTILIO

!...........   MODULES for public variables

!.........  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: RPT_

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FILNAM   ! physical or logical file name
    INTEGER     , INTENT(OUT) :: FDEV( RPT_%NUMFILES ) ! unit number of file
    INTEGER     , INTENT(OUT) :: MDEV                  ! unit number of file for src mapping

!...........   Local variables
    INTEGER         I, L, L2        ! counters and indices
    INTEGER         IOS             ! i/o status

    CHARACTER(3)    FILENO          ! tmp file number buffer
    CHARACTER(10)   FMT             ! tmp format buffer
    CHARACTER(16)   VARBUF          ! logical file name buffer
    CHARACTER(300)  MESG            ! message buffer
    CHARACTER(300)  PNAME           ! physical file name on ENVSTR test
    CHARACTER(300)  NNAME           ! new physical file name buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENREPOUT' ! program name

!***********************************************************************
!   begin body of subroutine OPENREPOUT

    IOS = -1
!.........  If file name is less than 16 characters, check if file name is a
!           defined environment variable
    IF( L .LE. 16 ) THEN
        MESG = 'Check of file name for logical status'
        VARBUF = FILNAM
        CALL ENVSTR( VARBUF, MESG, ' ', PNAME, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble '//FILNAM, 2 )
        END IF
    ELSE
        MESG = 'Logical file name "'//TRIM(FILNAM)//&
        &       '" is exceeding max 16 characters'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Open output file(s)
    DO I = 1, RPT_%NUMFILES

!.............  If it is a defined environment variable, open as logical file name
        IF( IOS .EQ. 0 ) THEN

!.................  If muliple files are being output then open them
            IF( RPT_%RPTMODE .EQ. 1 ) THEN

!.....................  Create new file name
                IF( I .GE. 10 ) THEN
                    FMT = '( A, I2 )'
                ELSE IF( I .GE. 100 ) THEN
                    FMT = '( A, I3 )'
                ELSE
                    FMT = '( A, I1 )'
                END IF

                WRITE( NNAME, FMT ) TRIM( PNAME )//'_', I

!.....................  Set logical file name to new file name
                IF( .NOT. SETENVVAR( FILNAM, NNAME ) ) THEN
                    MESG = 'Could not set logical file '//&
                    &       'name for file ' // TRIM( NNAME )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

            FDEV( I ) = GETEFILE( FILNAM, .FALSE.,&
            &                          .TRUE., PROGNAME )

!.............  If it is not a defined environment variable, open as physical filename
        ELSE

            IF( RPT_%RPTMODE .EQ. 1 ) THEN

                WRITE( FILENO, '( I3 )' ) I
                FILENO = ADJUSTL( FILENO )
                NNAME = FILNAM // '_' // FILENO

            END IF

            FDEV( I ) = JUNIT()
            OPEN( FDEV( I ),ERR=1006,FILE=NNAME,&
            &          STATUS='UNKNOWN',RECL=2500 )

        END IF

    END DO

!.........  Open source mapping ancillary output file
    IF( RPT_%SRCMAP ) THEN

        NNAME = TRIM( PNAME )//'_src_crosswalk.txt'
        IF( .NOT. SETENVVAR( 'SRC_CROSSWALK', NNAME ) ) THEN
            MESG = 'Could not set logical file '//&
            &       'name for file ' // TRIM( NNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MDEV = GETEFILE( 'SRC_CROSSWALK', .FALSE.,&
        &                      .TRUE., PROGNAME )

    END IF

    RETURN

!.........  Error opening raw input file
1006 WRITE( MESG,94010 ) ' Could not open file:' //&
    &       CRLF() // BLANK5 // FILNAM( 1:LEN_TRIM( FILNAM ) )
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE OPENREPOUT


