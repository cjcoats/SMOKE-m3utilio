
SUBROUTINE RDORSDSC( FDEV )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for and reads the ORIS FIPS
    !      codes and plant names
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 11/2001 by M. Houyoux
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
    !.........  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: ORISFIP, ORISLST, ORISDSC, NORIS

    IMPLICIT NONE

    !.........   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters
    LOGICAL, EXTERNAL :: USEEXPGEO

    !.........   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV              ! file unit number

    !.........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE

    !.........   Local paramaters
    INTEGER      , PARAMETER :: NFIELD   = 11           ! no. fields in input file
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDORSDSC'   !  program name

    !.........   Local arrays
    CHARACTER(60) SEGMENT( NFIELD )            ! line parsing array

    !.........   Local variables
    INTEGER         I                         ! indices and counters

    INTEGER         ENDLEN                    ! end length for reading descriptn
    INTEGER         IOS                       ! i/o status
    INTEGER      :: IREC                      ! record number
    INTEGER      :: NLINES                    ! number of lines in input file

    LOGICAL      :: EFLAG                     ! true: error found

    CHARACTER(256)  LINE                      ! Read buffer for a line
    CHARACTER(256)  MESG                      ! Message buffer

    CHARACTER(ORSLEN3 ) CORS              ! tmp ORIS ID
    CHARACTER(DSCLEN3 ) PDSC              ! tmp plant description

    !***********************************************************************
    !   Begin body of subroutine RDORSDSC

    EFLAG  = .FALSE.
    IREC   = 0
    NLINES = 0
    !.........  Write status message
    MESG = 'Reading ORIS descriptions file...'
    CALL M3MSG2( MESG )

    !.........  Get the number of lines in the holidays file
    NLINES = GETFLINE( FDEV, 'ORIS Descriptions' )

    !.........  Allocate memory for the SCC descriptions and initialize
    ALLOCATE( ORISFIP( NLINES ),        &
              ORISLST( NLINES ),        &
              ORISDSC( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ORISDSC', PROGNAME )

    !.........  Read the SCC descriptions, and store with SCC
    I = 0
    DO IREC = 1, NLINES

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                         &
                  'I/O error', IOS, 'reading ORIS '//   &
                  'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.........  Left adjust line
        LINE = ADJUSTL( LINE )

        !.........  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.........  Get SCC line
        CALL PARSLINE( LINE, NFIELD, SEGMENT )

        !.........  Set tmp ASCII fields into arrays of the correct length
        CORS = ADJUSTL( SEGMENT( 1 ) )
        CORS = ADJUSTR( CORS )

        PDSC = ADJUSTL( SEGMENT( 5 ) )

        !.........  Check for integer field for FIPS code
        IF ( .NOT. USEEXPGEO() .AND.        &
             .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: FIPS code not an integer at line', IREC
            CALL M3MSG2( MESG )
            CYCLE
        END IF

        !.........  Store entry
        I = I + 1
        ORISFIP( I ) = SEGMENT( 2 )
        ORISLST( I ) = CORS
        ORISDSC( I ) = PDSC

    END DO

    NORIS = I

    IF ( EFLAG ) THEN

        MESG = 'Problem reading ORIS descriptions.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.........  Successful completion
    RETURN

    !.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR: Unexpected end of SCC description file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )


    !.........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030 FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

END SUBROUTINE RDORSDSC
