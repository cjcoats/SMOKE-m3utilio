
SUBROUTINE RDMRCLIST( FDEV )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!       Reads the MRCLIST file (list of emission factor files for each
!       reference county). Checks that each reference county has a
!       file and that all files can be opened. Sorts files by reference
!       county.
!
!  PRECONDITIONS REQUIRED:
!       FDEV must be opened
!       MCREF file must be read already
!
!  SUBROUTINES AND FUNCTIONS CALLED:  none
!
!  REVISION  HISTORY:
!       04/10: Created by C. Seppanen
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

!.........  MODULES for public variables
!.........  This module is used for reference county information
    USE MODMBSET, ONLY: NREFC, MCREFIDX

!.........  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: MRCLIST, MVFILDIR

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV             ! MRCLIST file unit no.

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE

!...........   Local allocatable arrays
    INTEGER           , ALLOCATABLE :: IDX( : )     ! sorting index
    INTEGER           , ALLOCATABLE :: REFFIPA( : ) ! unsorted ref. county FIPs
    CHARACTER(FIPLEN3), ALLOCATABLE :: REFFIP( : )  ! sorted ref. county FIPs
    CHARACTER(100)    , ALLOCATABLE :: FILESA( : )  ! unsorted files names
    CHARACTER(100)    , ALLOCATABLE :: FILES( : )   ! sorted file names

    INTEGER, ALLOCATABLE :: MONTHA( : )      ! unsorted month numbers
    INTEGER, ALLOCATABLE :: MONTH( : )       ! sorted month numbers

!...........   Local arrays
    CHARACTER(100)  SEGMENT( 3 )          ! parsed input line

!...........   Other local variables
    INTEGER         I, J, N     ! counters and indexes
    INTEGER         IOS         ! error status
    INTEGER      :: IREC = 0    ! record counter
    INTEGER         NLINES      ! number of lines
    INTEGER         PFIP        ! previous ref. county FIP
    INTEGER         PMONTH      ! previous fuel month
    INTEGER         TMONTH      ! tmp. fuel month
    INTEGER      :: TDEV = 0    ! tmp. file unit
    INTEGER         TIDX        ! tmp. index

    LOGICAL      :: EFLAG = .FALSE.   ! true: error found
    LOGICAL      :: FOUND = .FALSE.   ! true: found data for reference county

    CHARACTER(150)     LINE     ! line buffer
    CHARACTER(200)     FILENAME ! tmp. filename
    CHARACTER(300)     MESG     ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDMRCLIST'   ! program name

!***********************************************************************
!   begin body of subroutine RDMRCLIST

!.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'Reference county factors file list' )

    ALLOCATE( REFFIPA( NLINES ),&
    &           MONTHA( NLINES ),&
    &           FILESA( NLINES ),&
    &           REFFIP( NLINES ),&
    &            MONTH( NLINES ),&
    &            FILES( NLINES ),&
    &              IDX( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'REFFIPA...IDX', PROGNAME )

    REFFIPA = 0
    MONTHA = 0
    FILESA = ' '
    IDX = 0

    DO I = 1, NLINES

        IDX( I ) = I

!.............  Read line
        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'I/O error', IOS,&
            &  'reading reference county factors file list at line',&
            &  IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Parse the line into 3 segments
        CALL PARSLINE( LINE, 3, SEGMENT )

!.............  Convert reference county to integer
        IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Bad reference county ' //&
            &  'FIPS code at line', IREC, 'of county factors ' //&
            &  'file list.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        REFFIPA( I ) = STR2INT( SEGMENT( 1 ) )

!.............  Convert month to integer
        IF( .NOT. CHKINT( SEGMENT( 2 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Bad month number ' //&
            &  'at line', IREC, 'of county factors file list.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        TMONTH = STR2INT( SEGMENT( 2 ) )

        IF( TMONTH .LT. 1 .OR. TMONTH .GT. 12 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Invalid fuel month ' //&
            &  'at line', IREC, 'of county factors file list.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        MONTHA( I ) = TMONTH
        FILESA( I ) = SEGMENT( 3 )

!.............  Check that file can be opened
        FILENAME = TRIM( MVFILDIR ) // TRIM( FILESA( I ) )
        OPEN( TDEV, FILE=FILENAME, STATUS='OLD', IOSTAT=IOS )
        IF( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Could not open file ' //&
            &  TRIM( FILESA( I ) ) // ' at line', IREC, 'of county ' //&
            &  'factors file list.'
            CALL M3MESG( MESG )
        END IF
        CLOSE( TDEV )

    END DO

    CLOSE( FDEV )

!.........  Check for errors while reading file
    IF( EFLAG ) THEN
        MESG = 'Problem reading reference county factors list'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Sort list by reference county and fuel month
    CALL SORTI2( NLINES, IDX, REFFIPA, MONTHA )

    REFFIP = ' '
    MONTH = 0
    FILES = ' '

    PFIP = -9
    PMONTH = -9
    N = 0
    DO I = 1, NLINES
        J = IDX( I )

        IF( REFFIPA( J ) == PFIP .AND. MONTHA( J ) == PMONTH ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: Duplicate entries in ' //&
            &  'reference county factors list for ' // CRLF() //&
            &  BLANK10 // 'reference county', PFIP, 'and fuel ' //&
            &  'month', PMONTH
            CALL M3MESG( MESG )
            CYCLE
        END IF

        N = N + 1
        WRITE( REFFIP( N ),'(I12.12)' ) REFFIPA( J )
        MONTH( N ) = MONTHA( J )
        FILES( N ) = FILESA( J )

        PFIP = REFFIPA( J )
        PMONTH = MONTHA( J )

    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem found in reference county factors list.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    ALLOCATE( MRCLIST( NREFC,12 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MRCLIST', PROGNAME )
    MRCLIST = ' '

!.........  Check that each reference county has a factors file
    TIDX = 1
    DO I = 1, NREFC

!.............  Loop through sorted lines in MRCLIST file - the
!               FIPS codes are sorted in the same order as the
!               reference counties in MCREFIDX
        FOUND = .FALSE.
        DO
            IF( REFFIP( TIDX ) == MCREFIDX( I,1 ) ) THEN
                FOUND = .TRUE.
                MRCLIST( I, MONTH( TIDX ) ) = FILES( TIDX )
            ELSE
                IF( FOUND ) EXIT
            END IF

            TIDX = TIDX + 1
            IF( TIDX .GT. NLINES ) EXIT
        END DO

        IF( .NOT. FOUND ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 ) 'ERROR: No factor file found ' //&
            &  'for reference county', MCREFIDX( I,1 ), 'in ' //&
            &  'reference county factors list.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem found in reference county factors list.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    DEALLOCATE( REFFIPA, MONTHA, FILESA, IDX, REFFIP, MONTH, FILES )

    RETURN

999 MESG = 'End of file'
    MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of MRCLIST' // CRLF() // BLANK5 //&
    &       'input file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDMRCLIST
