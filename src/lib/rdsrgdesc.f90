
SUBROUTINE RDSRGDESC( FDEV )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for surrogate code, region,
    !      description and file list from SRGDESC file.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 8/2005 by B. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO and related changes,
    !       ".f90" free source format, Bug-fixes for SORTI2, READ...END=
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

    !...........   This module contains the gridding surrogates desciption files
    USE MODSURG, ONLY: NSRGS, IDXSRGA, SRGLIST, SRGFCOD, SRGFREG,   &
                       SRGFDES, SRGFNAM, NTSRGDSC

    IMPLICIT NONE

    !...........   Subroutine arguments
    INTEGER, INTENT ( IN )  :: FDEV           ! file unit number

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !...........   Local variables
    INTEGER         I, J, K, N, NT               ! indices and counters

    INTEGER         MDEV                      ! for surrogate files
    INTEGER         ENDLEN                    ! end length for reading descriptn
    INTEGER         IOS                       ! i/o status
    INTEGER         NSRGALL                   ! No. entries in surrgoates file
    INTEGER      :: IREC = 0                  ! record number
    INTEGER      :: NLINES = 0                ! number of lines in input SRGDESC file
    INTEGER         SSC                       ! temporary spatial surrogate code
    INTEGER         LSSC                      ! previous temporary SSC

    CHARACTER(256)  LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer
    CHARACTER(60)   SEGMENT( 4 )               ! line parsing array

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSRGDESC'        !  program name

    !***********************************************************************
    !   Begin body of subroutine RDSRGDESC

    REWIND( FDEV )      ! In case of multiple calls

    !.........  Get the number of lines in the surrogate description file
    NLINES = GETFLINE( FDEV, 'SRGDESC Descriptions' )

    !.........  Allocate memory for the SRGDESC descriptions and initialize
    ALLOCATE( IDXSRGA( NLINES ),        &
              SRGFREG( NLINES ),        &
              SRGFCOD( NLINES ),        &
              SRGFDES( NLINES ),        &
              SRGFNAM( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IDXSRGA...SRGFNAM', PROGNAME )

    !.........  Read surrogate files in SRGDESC file and store
    N = 0
    DO I = 1, NLINES

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                             &
                  'I/O error', IOS, 'reading SRGDESC '//    &
                  'surrogate description file at line', I
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.............  Left adjust line
        LINE = ADJUSTL( LINE )

        !.............  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.............  Get line
        CALL PARSLINE( LINE, 4, SEGMENT )

        N = N + 1              ! actual line #s after skipping blank and comments

        IDXSRGA( N ) = N
        SRGFREG( N ) = ADJUSTL( SEGMENT( 1 ) )
        SRGFCOD( N ) = STR2INT( SEGMENT( 2 ) )
        SRGFDES( N ) = ADJUSTL( SEGMENT( 3 ) )
        SRGFNAM( N ) = ADJUSTL( SEGMENT( 4 ) )

    END DO

    NT = N
    NTSRGDSC = NT

    !.........  Sort surrogates by county code  cell  surrogate code
    CALL SORTI( NTSRGDSC, IDXSRGA, SRGFCOD )

    !.........  Count the surrogate codes and no of surrogate files
    LSSC = -1
    NSRGS = 0

    DO I = 1, NT

        J   = IDXSRGA( I )
        SSC = SRGFCOD( J )

        IF( SSC .NE. LSSC ) THEN
            NSRGS = NSRGS + 1
            LSSC  = SSC
        END IF

    END DO

    !.........  Allocate memory for derived surrogates tables
    ALLOCATE( SRGLIST( NSRGS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRGLIST', PROGNAME )

    !.........  Store the surrgage codes list and names of surrogate file
    LSSC = -1
    NSRGS = 0
    DO I = 1, NT

        J   = IDXSRGA( I )
        SSC = SRGFCOD( J )

        IF( SSC .NE. LSSC ) THEN
            NSRGS = NSRGS + 1
            SRGLIST( NSRGS ) = SSC
            LSSC = SSC
        END IF

    END DO

    !.........  Successful completion
    RETURN

    !.........  Unexpected end of file
998 MESG = 'INTERNAL ERROR : Unexpected end of surrogate file '
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    !******************  FORMAT  STATEMENTS   ******************************

93000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSRGDESC
