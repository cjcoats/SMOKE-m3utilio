
SUBROUTINE RDHDAYS( FDEV, SDATE, EDATE )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory for and reads the holiday dates.
    !      The current file structure permit holidays to be set globally for the
    !      entire inventory.  Future versions will permit setting holidays by
    !      region codes.  The data file contains the month and day of the holiday.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !      Created 9/2000 by M. Houyoux
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
    !.......   This module contains the temporal allocation information
    USE MODTMPRL, ONLY: NHOLIDAY, HOLREGN, HOLJDATE, HOLALTDY

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !.......   Subroutine arguments
    INTEGER, INTENT (IN) :: FDEV               ! file unit number
    INTEGER, INTENT (IN) :: SDATE              ! epsiode start Julian date
    INTEGER, INTENT (IN) :: EDATE              ! episode end Julian date

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local allocatable arrays
    INTEGER, ALLOCATABLE :: INDEXA  ( : )      ! sorting index
    INTEGER, ALLOCATABLE :: HOLREGNA( : )      ! unsorted holiday region
    INTEGER, ALLOCATABLE :: HDJDATEA( : )      ! unsorted holiday Julian date
    INTEGER, ALLOCATABLE :: HDALTDYA( : )      ! unsorted holiday Julian date

    !.......   Local arrays
    CHARACTER(9)    CAPDAYS( 7 )              ! names of days of week in caps
    CHARACTER(20)   SEGMENT( 5 )              ! parsed input line

    !.......   Local variables
    INTEGER         I, J, N                   ! indices and counters

    INTEGER         IOS                       ! i/o status
    INTEGER         IREC                      ! record number
    INTEGER         DAY                       ! tmp day of week
    INTEGER         DD                        ! tmp day of month
    INTEGER         JDATE                     ! Julian date
    INTEGER         NLINES                    ! no. lines in input file
    INTEGER         MM                        ! tmp month
    INTEGER         REGN                      ! tmp region code
    INTEGER         YY                        ! tmp year

    LOGICAL      :: EFLAG = .FALSE.           ! true: error found

    CHARACTER(80)   LINE                      ! Read buffer for a line
    CHARACTER(300)  MESG                      ! Message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDHDAYS'        !  program name

    !***********************************************************************
    !   Begin body of subroutine RDHDAYS

    !.......  Get the number of lines in the holidays file
    NLINES = GETFLINE( FDEV, 'Holidays file' )

    !.......  Allocate memory for the unsorted holidays data
    ALLOCATE(   INDEXA( NLINES ),       &
              HOLREGNA( NLINES ),       &
              HDJDATEA( NLINES ),       &
              HDALTDYA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'HDALTDYA', PROGNAME )

    !.......  Create array of the days of the week in capital letters
    DO I = 1, 7
        CAPDAYS( I ) = DAYS( I )
        CALL UPCASE( CAPDAYS( I ) )
    END DO

    !.......  Read the unsorted, unfiltered holidays data
    I = 0
    IREC = 0
    DO N = 1, NLINES

        READ ( FDEV, 93000, END=99, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010) 'I/O error', IOS, 'reading '//&
                                'holidays file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.......  Parse the line into its 2 segments
        CALL PARSLINE( LINE, 5, SEGMENT )

        !.......  Convert first column to region code, if possible
        IF( CHKINT( SEGMENT( 1 ) ) ) THEN
            REGN = STR2INT( SEGMENT( 1 ) )

        !.......  Temporary message for unsupported region codes
            IF( REGN .NE. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Non-zero region codes in holidays '//    &
                       'file are not yet supported.'
                CALL M3MESG( MESG )
            END IF

        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Bad region code at line',   &
                   IREC, 'of holidays file.'
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.......  Convert second column to Julian date from Gregorian date...
        !.......  Convert Month
        IF( CHKINT( SEGMENT( 2 ) ) ) THEN
            MM = STR2INT( SEGMENT( 2 ) )

        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Bad month at line',     &
                   IREC, 'of holidays file.'
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.......  Convert day
        IF( CHKINT( SEGMENT( 3 ) ) ) THEN
            DD = STR2INT( SEGMENT( 3 ) )

        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Bad day at line',       &
                   IREC, 'of holidays file.'
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.......  Convert year
        IF( CHKINT( SEGMENT( 4 ) ) ) THEN
            YY = STR2INT( SEGMENT( 4 ) )

            !.......  Convert to 4-digit year, if needed
            IF( YY .LT. 100 ) THEN
                YY = YEAR4( YY )

            ELSEIF( YY .LT. 1970 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad year at line',      &
                       IREC, 'of holidays file.'
                CALL M3MESG( MESG )
                CYCLE

            END IF

        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Bad year at line',      &
                       IREC, 'of holidays file.'
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.......  Determine Julian date
        JDATE = YY * 1000 + JULIAN( YY, MM, DD )

        !.......  If date is outside episode range, ignore
        !            IF( JDATE .LT. SDATE .OR. JDATE .GT. EDATE ) CYCLE

        !.......  Determine alternative day of the week
        CALL UPCASE( SEGMENT( 5 ) )
        DAY = INDEX1( SEGMENT( 5 ), 7, CAPDAYS )

        IF( DAY .LE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Bad alternative day of ',&
                   'week at line', IREC, 'of holidays file.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Store holiday records that match episode

        I = I + 1
        INDEXA  ( I ) = I
        HOLREGNA( I ) = REGN
        HDJDATEA( I ) = JDATE
        HDALTDYA( I ) = DAY

    END DO

    NHOLIDAY = I

    !.......  Abort if error found in holidays file
    IF( EFLAG ) THEN
        MESG = 'Problem reading holidays file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Allocate memory for the sorted holidays data
    ALLOCATE( HOLREGN( NHOLIDAY ),      &
             HOLJDATE( NHOLIDAY ),      &
             HOLALTDY( NHOLIDAY ), STAT=IOS )
    CALL CHECKMEM( IOS, 'HOLALTDY', PROGNAME )

    !.......  Sort the index of the dates
    CALL SORTI2( NHOLIDAY, INDEXA, HOLREGNA, HDJDATEA )

    !.......  Store the sorted holidays
    DO I = 1, NHOLIDAY

        J = INDEXA( I )
        HOLREGN ( I ) = HOLREGNA( J )
        HOLJDATE( I ) = HDJDATEA( J )
        HOLALTDY( I ) = HDALTDYA( J )

    END DO

    !.......  Deallocate local memory
    DEALLOCATE( INDEXA, HOLREGNA, HDJDATEA, HDALTDYA )

    !.......  Successful completion
    RETURN

    !.......  Unexpected end of file
99  MESG = 'INTERNAL ERROR: Unexpected end of holidays file'
    CALL M3MSG2( MESG )

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )


    !.......   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( A, 1X, I8, 1X, A, 1X, F10.6, 1X, A )

94030 FORMAT( A, 1X, I6.6, A, 100( ' SSC(', I2.2, '):', F10.6, : ) )

END SUBROUTINE RDHDAYS
