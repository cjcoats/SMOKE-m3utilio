
SUBROUTINE RDDATES( PFLAG, FDEV, NTPERIOD )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for and reads the time periods.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 8/2005 by B. Baek
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

    !.......   Modules for public variables
    !.......   This module contains the lists of unique source characteristics

    USE MODTMPRL, ONLY: STDATE, STTIME, RUNLEN, ITDATE

    IMPLICIT NONE

    !.......   Subroutine arguments

    LOGICAL, INTENT(IN   ) :: PFLAG
    INTEGER, INTENT(IN   ) :: FDEV           ! file unit number
    INTEGER, INTENT(  OUT) :: NTPERIOD       ! No of time periods

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !.......   Local variables
    INTEGER         I, N                      ! indices and counters
    INTEGER         TSTEP

    INTEGER         ENDLEN                    ! end length for reading descriptn
    INTEGER         IOS                       ! i/o status
    INTEGER         IREC                      ! record number
    INTEGER         NLINES                    ! number of lines in input file

    CHARACTER(256)  LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer
    CHARACTER(60)   SEGMENT( 3 )              ! line parsing array

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDDATES'        !  program name

    !***********************************************************************
    !   Begin body of subroutine RDDATES

    IF ( .NOT. PFLAG ) THEN
        NTPERIOD = 1
        ALLOCATE( ITDATE( 1 ),      &
                  STDATE( 1 ),      &
                  STTIME( 1 ),      &
                  RUNLEN( 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ITDATE,STDATE,STTIME,RUNLEN', PROGNAME )
        RETURN
    END IF


    REWIND( FDEV )      ! In case of multiple calls

    !.......  Get the number of lines in the PROCDATES file

    NLINES = GETFLINE( FDEV, 'PROCDATES Descriptions' )

    !.......  Allocate memory for the PROCDATES descriptions and initialize
    ALLOCATE( ITDATE( NLINES ),     &
              STDATE( NLINES ),     &
              STTIME( NLINES ),     &
              RUNLEN( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ITDATE,STDATE,STTIME,RUNLEN', PROGNAME )

    !.......  Read the PROCDATE file and store STDATE, ETDATE, and RUNLEN

    N    = 0
    IREC = 0
    DO I = 1, NLINES

        READ ( FDEV, 93000, END=998, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                             &
                  'I/O error', IOS, 'reading PROCDATES '//  &
                  'description file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        LINE = ADJUSTL( LINE )

        !.......  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.......  Get line
        CALL PARSLINE( LINE, 3, SEGMENT )

        N = N + 1              ! actual line #s after skipping blank and comments

        STDATE( N ) = STR2INT( SEGMENT( 1 ) )
        STTIME( N ) = STR2INT( SEGMENT( 2 ) )
        RUNLEN( N ) = STR2INT( SEGMENT( 3 ) )

    END DO

    NTPERIOD = N

    !.......  Successful completion
    RETURN

    !.......  Unexpected end of file
998 MESG = 'INTERNAL ERROR : Unexpected end of PROCDATES desc. file'

    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

93000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )


END SUBROUTINE RDDATES
