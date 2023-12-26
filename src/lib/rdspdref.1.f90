
SUBROUTINE RDSPDREF( FDEV )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!     Reads the SPDREF file that contains cross-reference information
!     for assigning speed profiles to mobile sources
!
!  PRECONDITIONS REQUIRED:
!     File unit FDEV already is opened.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 12/02 by C. Seppanen (based off rdgref.f and rdxclude.f)
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
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
!***************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module is for cross reference tables
    USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, ISPDCDA

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV   ! SPDREF file unit no.

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE

!...........   Local parameters
    INTEGER    , PARAMETER :: MXCOL = 3

!...........   Array of input fields
    CHARACTER(SCCLEN3)  SEGMENT( MXCOL )

!...........   Other local variables
    INTEGER         I, N    !  counters and indices

    INTEGER         IOS     !  i/o status
    INTEGER         IREC    !  record counter
    INTEGER         ISPD    !  tmp spd profile code
    INTEGER         NLINES  !  number of lines
    INTEGER         NXREF   !  number of valid x-ref entries

    LOGICAL      :: EFLAG = .FALSE.   !  true: error found

    CHARACTER(10)      FIPFMT   !  formt to write co/st/cy to string
    CHARACTER(128)     LINE     !  line buffer
    CHARACTER(256)     MESG     !  message buffer
    CHARACTER(FIPLEN3) CFIP     !  buffer for CFIPS code
    CHARACTER(SCCLEN3) TSCC     !  temporary SCC

    CHARACTER(16) :: PROGNAME = 'RDSPDREF' ! program name

!***********************************************************************
!   begin body of subroutine RDSPDREF

!.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'Speeds cross-reference file' )

!.............  Set up formats
    WRITE( FIPFMT, '("(I",I2.2,".",I2.2,")")' ) FIPLEN3, FIPLEN3

!.........  Allocate memory for unsorted data used in all source categories
    ALLOCATE( ISPDCDA( NLINES ),&
    &           INDXTA( NLINES ),&
    &           CSCCTA( NLINES ),&
    &           CSRCTA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )
    ISPDCDA = 0  ! array
    INDXTA = 0   ! array
    CSCCTA = ' ' ! array
    CSRCTA = ' ' ! array

!.........  Set up constants for loop.

!.........  Second pass through file: read lines and store unsorted data for
!           the source category of interest
    IREC = 0
    N    = 0

    DO I = 1, NLINES

        READ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

!.............  Check for I/O errors
        IF( IOS > 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'I/O error', IOS,&
            &    'reading speed cross-reference file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Check for end of file
        IF( IOS < 0 ) THEN
            MESG = 'End of file reached unexpectedly. ' //&
            &       'Check format of speed' // CRLF() // BLANK5 //&
            &       'cross-reference file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Parse line into segments
        CALL PARSLINE( LINE, MXCOL, SEGMENT )

!.............  Make sure that the co/st/cy code is an integer
        IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Country/state/county ' //&
            &       'code is not an integer at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Make sure that SCC is full length
        IF( LEN_TRIM( SEGMENT( 2 ) ) /= SCCLEN3 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: SCC value is not ',&
            &       SCCLEN3, 'digits at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Make sure that the speed profile number is an integer
        IF( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //&
            &       'is not an integer at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Check that profile number is not too long
        IF( LEN_TRIM( ADJUSTL( SEGMENT( 3 ) ) ) > SPDLEN3 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //&
            &    'is longer than', SPDLEN3, 'digits at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Convert speed profile code to an integer
        ISPD = STR2INT( SEGMENT( 3 ) )

!.............  Check that profile number is greater than zero
        IF( ISPD <= 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Speed profile code ' //&
            &       'is less than or equal to zero at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Skip rest of loop if there has been an error
        IF( EFLAG ) CYCLE

!.............  Create co/st/cy code string with leading zeros
        WRITE( CFIP,FIPFMT ) STR2INT( SEGMENT( 1 ) )

!.............  Save SCC in string
        TSCC = SEGMENT( 2 )

        N = N + 1
        IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

!.............  Store fields
        INDXTA ( N ) = N
        ISPDCDA( N ) = ISPD
        CSCCTA ( N ) = TSCC
        CSRCTA ( N ) = CFIP // TSCC

    END DO   ! End of loop reading SPDREF file

!.........  Reset number of cross-reference entries in case some were dropped
    NXREF = N

!.........  Write errors for problems with input
    IF( NXREF .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No valid entries in SPDREF file!'
        CALL M3MSG2( MESG )

    ELSEIF( NXREF .GT. NLINES ) THEN
        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //&
        &       'storing SPDREF file was', NLINES,&
        &       CRLF() // BLANK10 // 'but actually needed', NXREF
        CALL M3MSG2( MESG )

    ENDIF

!.........  Check for errors reading SPDREF file and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading speed cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    CALL M3MSG2( 'Processing speed cross-reference file...' )

!.........  Sort cross-reference table
    CALL SORTIC( NXREF, INDXTA, CSRCTA )

!.........  Group cross-reference data into tables for different groups
    CALL XREFTBL( 'SPEED', NXREF )

!.........  Deallocate cross-reference sorting arrays
    DEALLOCATE( ISPDCDA, INDXTA, CSCCTA, CSRCTA )

!.........  Rewind file
    REWIND( FDEV )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSPDREF
