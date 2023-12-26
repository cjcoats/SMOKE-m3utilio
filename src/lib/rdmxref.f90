
SUBROUTINE RDMXREF( MDEV, NCTY, GRIDCTY )

    !***********************************************************************
    !  subroutine body starts at line 107
    !
    !  DESCRIPTION:
    !       Reads the MCREF file, checks that each inventory county is assigned
    !       a reference county, ignores counties outside of the grid, and sorts the data
    !
    !  PRECONDITIONS REQUIRED:
    !       MDEV must be opened
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  none
    !
    !  REVISION  HISTORY:
    !      10/01: Created by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
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

    !.........  MODULES for public variables
    !.........  This module is used for MOBILE6 setup information
    USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX

    IMPLICIT NONE

    !.........   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.........   SUBROUTINE ARGUMENTS
    INTEGER,        INTENT ( IN ) :: MDEV                 ! MCREF file unit no.
    INTEGER,        INTENT ( IN ) :: NCTY                 ! no. counties
    CHARACTER( * ), INTENT ( IN ) :: GRIDCTY( NCTY )      ! counties in grid

    !.........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE

    !.........   Local allocatable arrays
    INTEGER, ALLOCATABLE :: MCREFRAW( :,: )      ! raw MCREF data
    INTEGER, ALLOCATABLE :: IDX     ( : )        ! index into MCREF data

    !.........   Local arrays
    CHARACTER(3)  SEGMENT( 6 )              ! parsed input line

    !.........   Other local variables
    INTEGER I, J, N, K                    ! counters and indices


    INTEGER :: CO = 0                     ! tmp country
    INTEGER :: ST = 0                     ! tmp state
    INTEGER :: CT = 0                     ! tmp county
    INTEGER :: IOS  = 0                   ! I/O status
    INTEGER :: IREC = 0                   ! record counter
    INTEGER :: NREF = 0                   ! number of ref. counties
    INTEGER :: NLINES = 0                 ! number of lines

    INTEGER REFCOUNTY                     ! ref. county FIPS code
    INTEGER INVCOUNTY                     ! inv. county FIPS code
    INTEGER PRCOUNTY                      ! previous ref. county
    INTEGER PICOUNTY                      ! previous inv. county

    LOGICAL     EFLAG                ! true: error found

    CHARACTER(FIPLEN3) REFCNTY, INVCNTY, PRVCNTY
    CHARACTER(100)     LINE         !  line buffer
    CHARACTER(300)     MESG         !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDMXREF'       ! program name

    !***********************************************************************
    !   begin body of subroutine RDMXREF

    !.........  Get the number of lines in the file
    EFLAG  = .FALSE.
    NLINES = GETFLINE( MDEV, 'County cross-reference file' )

    !.........  Allocate arrays that use NLINES as dimension (unsorted arrays and index)
    ALLOCATE( MCREFRAW ( NLINES,2 ),            &
                   IDX ( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MCREFRAW,IDX', PROGNAME )

    !.........  Initialize arrays
    IDX = 0
    MCREFRAW = 0.
    PICOUNTY = 0
    N = 0

    DO I = 1, NLINES

        IDX( I ) = I

    !.........  Read line
        READ( MDEV, 93000, END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS /= 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )        &
                'I/O error', IOS, 'reading county cross-reference file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

    !.........  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        CALL PARSLINE( LINE, 6, SEGMENT )

    !.........  Check the format of input
        IF( .NOT. CHKINT( SEGMENT(1) ) ) EFLAG = .TRUE.
        IF( .NOT. CHKINT( SEGMENT(2) ) ) EFLAG = .TRUE.
        IF( .NOT. CHKINT( SEGMENT(3) ) ) EFLAG = .TRUE.
        IF( .NOT. CHKINT( SEGMENT(4) ) ) EFLAG = .TRUE.
        IF( .NOT. CHKINT( SEGMENT(5) ) ) EFLAG = .TRUE.
        IF( .NOT. CHKINT( SEGMENT(6) ) ) EFLAG = .TRUE.
        IF( SEGMENT(5)==' ' .OR. SEGMENT(6)==' ' ) EFLAG=.TRUE.

        IF( EFLAG ) THEN
            WRITE( MESG,94010 ) 'ERROR: Bad inventory county ' //&
                       'FIPS code at line', IREC, 'of county ' //&
                       'cross-reference file.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

    !.........  Convert inventory county to integer
        CO = STR2INT( SEGMENT( 1 ) )
        ST = STR2INT( SEGMENT( 2 ) )
        CT = STR2INT( SEGMENT( 3 ) )
        INVCOUNTY = CO*100000 + ST*1000 + CT

        CO = STR2INT( SEGMENT( 4 ) )
        ST = STR2INT( SEGMENT( 5 ) )
        CT = STR2INT( SEGMENT( 6 ) )
        REFCOUNTY = CO*100000 + ST*1000 + CT

    !.........  Store values in unsorted array
        MCREFRAW( I,1 ) = INVCOUNTY
        MCREFRAW( I,2 ) = REFCOUNTY

    !.........  Skip any entries equal to zero due to blank lines
        IF( REFCOUNTY == 0 .OR. INVCOUNTY == 0 ) CYCLE

    !.........  Check if current inventory county is duplicate (match previous)
        IF( INVCOUNTY /= PICOUNTY ) THEN

    !.........  Check that current county is inside the grid (and in the inventory)
            WRITE( INVCNTY,'(I12.12)' ) INVCOUNTY
            K = FINDC( INVCNTY, NCTY, GRIDCTY )

            IF( K > 0 ) N = N + 1

        END IF

        PICOUNTY = INVCOUNTY

    END DO      ! done reading MCREF file

    NINVC = N

    !.........  Allocate arrays that use NINVC as dimension (sorted arrays and index)
    ALLOCATE( MCREFSORT ( NINVC,2 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MCREFSORT', PROGNAME )
    MCREFSORT = ' '

    !.........  Close MCREF file
    CLOSE( MDEV )

    !.........  Abort if error found while reading cross-reference file
    IF( EFLAG ) THEN
        MESG = 'Problem reading cross-reference file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Sort MCREF index array by reference county
    CALL SORTI2( NLINES, IDX, MCREFRAW(:,2), MCREFRAW(:,1) )

    !.........  Check for duplicate entries and counties outside grid,
    !           then store sorted MCREF array
    N = 0
    PRVCNTY = ' '

    DO I = 1, NLINES

        J = IDX( I )

        WRITE( REFCNTY,'(I12.12)' )  MCREFRAW( J,2 )
        WRITE( INVCNTY,'(I12.12)' )  MCREFRAW( J,1 )

        !.........  Skip any entries equal to zero due to blank lines
        IF( STR2INT( REFCNTY ) == 0 ) CYCLE
        IF( STR2INT( INVCNTY ) == 0 ) CYCLE

        !.........  Check if current inventory county is duplicate (match previous)
        IF( INVCNTY == PRVCNTY ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Duplicate found in county ' //       &
                   'cross-reference file : ' // INVCNTY
            CALL M3MESG( MESG )

        ELSE

            !.........  Check that current county is inside the grid (and in the inventory)
            K = FINDC( INVCNTY, NCTY, GRIDCTY )

            IF( K < 0 ) THEN
                MESG = 'WARNING: Ignoring county cross-reference ' //       &
                       'for inventory county ' // INVCNTY // ' , ' //       &
                       CRLF() // BLANK10 // ' since it is not inside the grid'
                CALL M3MESG( MESG )

            ELSE
                N = N + 1
                MCREFSORT( N,1 ) = INVCNTY
                MCREFSORT( N,2 ) = REFCNTY
            END IF

        END IF

        PRVCNTY = INVCNTY

    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem(s) found in county cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Count total number of reference counties
    PRVCNTY = ' '

    DO I = 1, NINVC
        REFCNTY = MCREFSORT( I,2 )

        IF( REFCNTY /= PRVCNTY ) THEN
            NREF = NREF + 1
        END IF

        PRVCNTY = REFCNTY

    END DO

    NREFC = NREF

    !.........  Create reference county index array
    ALLOCATE( MCREFIDX ( NREFC,2 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MCREFIDX', PROGNAME )
    MCREFIDX = ' '

    N = 0
    PRVCNTY = ' '

    DO I = 1, NINVC
        REFCNTY = MCREFSORT( I,2 )

        IF( REFCNTY /= PRVCNTY ) THEN
            N = N + 1
            MCREFIDX( N,1 ) = REFCNTY
            WRITE( MCREFIDX( N,2 ),'(I8)' ) I
        END IF

        PRVCNTY = REFCNTY
    END DO

    !.........  Deallocate local memory
    DEALLOCATE( MCREFRAW, IDX )

    RETURN

999 MESG = 'End of file'
    MESG = 'End of file reached unexpectedly. ' //          &
           'Check format of MCREF' // CRLF() // BLANK5 //   &
           'input file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !.........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )
94020 FORMAT( 3( A, 1X ), I8, 1X, A, 1X )

END SUBROUTINE RDMXREF
