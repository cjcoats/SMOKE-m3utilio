
SUBROUTINE RDMACTDSC( FDEV )

!**************************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine allocates memory for and reads the MACT descriptions.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 7/2005 by B. Baek
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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
!**************************************************************************
    USE M3UTILIO

!...........   Modules for public variables
!...........   This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: MACTDESC, NINVMACT, INVMACT

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV          ! file unit number

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

!.........  Local parameters
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDMACTDSC'    !  program name
    INTEGER      , PARAMETER :: MXSEG = 5       ! max # of potential line segments

!.........  Other arrays
    CHARACTER( SDSLEN3 ) SEGMENT( MXSEG ) ! Segment of parsed lines

!...........   Local variables
    INTEGER         J, L, N               ! indices and counters

    INTEGER         ENDLEN                ! end length for reading descriptn
    INTEGER         IOS                   ! i/o status
    INTEGER      :: IREC = 0              ! record number
    INTEGER      :: NLINES = 0            ! number of lines in input file

    LOGICAL      :: EFLAG = .FALSE.       ! true: error found

    CHARACTER(256)  LINE                  ! Read buffer for a line
    CHARACTER(300)  MESG                  ! Message buffer

    CHARACTER( MACLEN3 ) TMACT            ! tmp MACT

!***********************************************************************
!   Begin body of subroutine RDMACTDSC

    REWIND( FDEV )  ! In case of multiple calls

!.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'MACT Descriptions' )

!.........  Allocate memory for the MACT descriptions and initialize
    ALLOCATE( MACTDESC( NINVMACT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MACTDESC', PROGNAME )
    MACTDESC = 'Description unavailable'          ! array

!.........  Read the MACT descriptions, and store with MACT
    ENDLEN = MACLEN3 + SDSLEN3

    DO N = 1, NLINES

        READ ( FDEV, '( A )', END=998, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)&
            &      'I/O error', IOS, 'reading MACT '//&
            &      'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Left adjust line
        LINE = ADJUSTL( LINE )

!.............  Get MACT line
        CALL PARSLINE( LINE, 2, SEGMENT )
        TMACT = SEGMENT( 1 )( 1:MACLEN3 )

        CALL PADZERO( TMACT )

!.............  Find MACT in inventory list, and if it's in the inventory,
!               store the description.
        J = FINDC( TMACT, NINVMACT, INVMACT )

        IF ( J .GT. 0 ) THEN
            MACTDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
        END IF

    END DO

!.........  Successful completion
    RETURN

!.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR: Unexpected end of MACT description file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDMACTDSC
