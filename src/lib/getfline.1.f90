
INTEGER FUNCTION GETFLINE( IDEV, DESCRIPT )

!***********************************************************************
!  function body starts at line
!
!  DESCRIPTION:
!     Counts the number of lines in an ASCII file
!
!  PRECONDITIONS REQUIRED:
!     File opened and unit number provided
!
!  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
!
!  REVISION  HISTORY:
!       prototype 10/98 by M Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   ARGUMENTS and their descriptions:

    INTEGER       IDEV         ! Unit number for ASCII file
    CHARACTER(*)  DESCRIPT     ! Description of file

!...........   LOCAL VARIABLES their descriptions:

    INTEGER       ICNT         ! line counter
    INTEGER       IOS          ! i/o status

    CHARACTER      BUFFER      ! ASCII LINE from X-ref file
    CHARACTER(256) MESG        ! Message buffer

    CHARACTER(16) :: PROGNAME = 'GETFLINE' ! program name

!***********************************************************************
!   begin body of subroutine GETFLINE
    REWIND( IDEV )

    ICNT = 0

!.........  Loop through lines of file, counting the lines
11  CONTINUE

    READ( IDEV, 93000, IOSTAT=IOS, END=22 ) BUFFER

    ICNT = ICNT + 1

    IF ( IOS .GT. 0 ) THEN
        WRITE( MESG,94010 )&
        &      'I/O error', IOS,&
        &      'scanning ' // TRIM( DESCRIPT ) //&
        &      ' file at line', ICNT
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    GO TO 11

22  CONTINUE

    GETFLINE = ICNT

    REWIND( IDEV )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I6, :, 2X ) )

END FUNCTION GETFLINE
