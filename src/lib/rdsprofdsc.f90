
SUBROUTINE RDSPROFDSC( FDEV )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for and reads the GSPRO descriptions.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 4/2004 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO".f90" free source format, and
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
    !.........  This module contains the speciation-profiles matrix, among other things.
    USE MODSOURC, ONLY : NGSPRO, GSPROID, GSPRDESC

    IMPLICIT NONE

    !...........   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !...........   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV              ! file unit number

    !.........  Local parameters
    INTEGER,       PARAMETER :: MXSEG = 5           ! max # of potential line segments
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSPROFDSC'        !  program name

    !.........  Other arrays
    CHARACTER( SDSLEN3 ) SEGMENT( MXSEG )     ! Segment of parsed lines

    !...........   Local variables
    INTEGER         J, L, N, NS              ! indices and counters

    INTEGER         IOS                       ! i/o status
    INTEGER      :: IREC = 0                  ! record number
    INTEGER      :: NLINES = 0                ! number of lines in input file

    LOGICAL      :: EFLAG = .FALSE.           ! true: error found

    CHARACTER(256)  LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer

    CHARACTER(SPNLEN3) CSPF                   ! tmp GSPRO

    !***********************************************************************
    !   Begin body of subroutine RDSPROFDSC

    REWIND( FDEV )      ! In case of multiple calls

    !.........  Get the number of lines in the file
    N = GETFLINE( FDEV, 'GSPRO Speciation Profile Descriptions' )

    !.........  Allocate memory for the GSPRO descriptions and initialize
    ALLOCATE( GSPRDESC( NGSPRO ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GSPRDESC', PROGNAME )
    GSPRDESC = 'Description unavailable'              ! array

    !.........  Read the GSPRO descriptions, and store with GSPRO
    DO N = 1, N

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                         &
                  'I/O error', IOS, 'reading GSPRO '//  &
                  'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.............  If it is not header, Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.............  Left adjust line
        LINE = ADJUSTL( LINE )

        !.............  Get GSPRODESC line
        CALL PARSLINE( LINE, 2, SEGMENT )
        CSPF = ADJUSTR( SEGMENT( 1 )( 1:SPNLEN3 ) )

        !.............  Look for matched GSPRO
        J = INDEX1( CSPF, NGSPRO, GSPROID )
        IF ( J > 0 ) THEN
            GSPRDESC( J ) = ADJUSTL( SEGMENT( 2 ) )
        END IF

    END DO

    !.........  Successful completion
    RETURN

    !.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR: Unexpected end of GSPRO description file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )


    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSPROFDSC
