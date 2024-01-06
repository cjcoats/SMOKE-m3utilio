
SUBROUTINE RDNAICSDSC( FDEV )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for and reads the NAICS descriptions.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/2005 by B. Baek
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
    !**************************************************************************
    USE M3UTILIO

    !.........   Modules for public variables
    !.........   This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NAICSDESC, NINVNAICS, INVNAICS

    IMPLICIT NONE

    !.........   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.........   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV              ! file unit number

    !.........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !.........  Local parameters
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDNAICSDSC'        !  program name
    INTEGER      , PARAMETER :: MXSEG = 5           ! max # of potential line segments

    !.........  Other arrays
    CHARACTER( SDSLEN3 ) SEGMENT( MXSEG )     ! Segment of parsed lines

    !.........   Local variables
    INTEGER         J, L, N                   ! indices and counters

    INTEGER         ENDLEN                    ! end length for reading descriptn
    INTEGER         IOS                       ! i/o status
    INTEGER         IREC                      ! record number
    INTEGER         NLINES                    ! number of lines in input file
    LOGICAL         EFLAG                     ! true: error found

    CHARACTER(256)  LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer

    CHARACTER( NAILEN3 ) TNAICS               ! tmp NAICS

    !***********************************************************************
    !   Begin body of subroutine RDNAICSDSC

    REWIND( FDEV )      ! In case of multiple calls
    NLINES = 0
    IREC   = 0
    EFLAG  = .FALSE.


    !.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'NAICS Descriptions' )

    !.........  Allocate memory for the NAICS descriptions and initialize
    ALLOCATE( NAICSDESC( NINVNAICS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NAICSDESC', PROGNAME )
    NAICSDESC = 'Description unavailable'              ! array

    !.........  Read the NAICS descriptions, and store with NAICS
    ENDLEN = NAILEN3 + SDSLEN3
    DO N = 1, NLINES

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                         &
                  'I/O error', IOS, 'reading NAICS '//  &
                  'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.........  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.........  Left adjust line
        LINE = ADJUSTL( LINE )

        !.........  Get NAICS line
        CALL PARSLINE( LINE, 2, SEGMENT )
        TNAICS = SEGMENT( 1 )( 1:NAILEN3 )
        CALL PADZERO( TNAICS )

        !.........  Find NAICSS in inventory list, and if it's in the inventory,
        !           store the description.
        J = FINDC( TNAICS, NINVNAICS, INVNAICS )
        IF ( J .GT. 0 ) THEN
            NAICSDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
        END IF

    END DO

    !.........  Successful completion
    RETURN

    !.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR: Unexpected end of NAICS description file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !.........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDNAICSDSC
