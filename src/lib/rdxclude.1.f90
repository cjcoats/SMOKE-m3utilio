
SUBROUTINE RDXCLUDE( FDEV )

!***********************************************************************
!  subroutine body starts at line 98
!
!  DESCRIPTION:
!     Reads the NHAPEXCLUDE file that contains a list of country/
!     state/county FIPS codes and SCCs that are to be excluded from
!     calculation of the NONHAPVOC pollutants when combining criteria
!     and toxics inventories.
!
!  PRECONDITIONS REQUIRED:
!     File unit FDEV already is opened.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 11/02 by M. Houyoux
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
    USE MODXREF, ONLY: INDXTA, CFIPTA, CSRCTA, CSCCTA, NHAP_EXCL

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV   ! NHAPEXCLUDE file unit no.

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: USEEXPGEO

!...........   Local parameters
    INTEGER    , PARAMETER :: MXCOL = 8

!...........   Array of input fields
    CHARACTER(CHRLEN3)  SEGMENT( MXCOL )

!...........   Array of point source plant characeristics
    CHARACTER(CHRLEN3) CHARS( 5 )

!...........   Other local variables
    INTEGER         I, N, L0, L1, L2   !  counters and indices

    INTEGER         IOS     !  i/o status
    INTEGER         IREC    !  record counter
    INTEGER         NLINES  !  number of lines
    INTEGER         NXREF   !  number of valid x-ref entries

    LOGICAL      :: HDRFLAG = .FALSE.
    LOGICAL      :: EFLAG  = .FALSE.      !  true: error found

    CHARACTER(128)     LINE     !  line buffer
    CHARACTER(256)     MESG     !  message buffer
    CHARACTER(FIPLEN3) CFIP     !  buffer for FIPS code
    CHARACTER(PLTLEN3) PLT      !  temporary plant ID
    CHARACTER(SCCLEN3) TSCC     !  temporary SCC
    CHARACTER(ALLLEN3) CSRCALL  !  buffer for source char, incl pol

    CHARACTER(16) :: PROGNAME = 'RDXCLUDE' ! program name

!***********************************************************************
!   begin body of subroutine RDXCLUDE

!.........  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'non-HAP inclusion/exclusions file' )

!.........  Allocate memory for unsorted data used in all source categories
    ALLOCATE( INDXTA( NLINES ),&
    &          CFIPTA( NLINES ),&
    &          CSCCTA( NLINES ),&
    &          CSRCTA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INDXTA...CSRCTA', PROGNAME )
    INDXTA = 0   ! array
    CFIPTA = ' ' ! array
    CSCCTA = ' ' ! array
    CSRCTA = ' ' ! array

!.........  Set up constants for loop.

!.........  Second pass through file: read lines and store unsorted data for
!           the source category of interest
    NHAP_EXCL = .TRUE.     ! true: default non-integrate (exclusion)
    IREC   = 0
    N      = 0
    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'I/O error', IOS,&
            &    'reading non-HAP exclusions/inclusions file at line'&
            &    , IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Search for header to define either exclusions or inclusion
        L0 = INDEX( LINE, '/EXCLUDE/' )
        L1 = INDEX( LINE, '/INCLUDE/' )
        L2 = INDEX( LINE, '/END/' )

        IF( L0 > 0 ) CYCLE         ! process non-HAP exclusions

        IF( L1 > 0 ) THEN
            NHAP_EXCL = .FALSE.    ! process non-HAP inclusions
            CYCLE
        ELSE IF( L2 > 0 ) THEN
            HDRFLAG = .TRUE.
            CYCLE
        END IF

!.............  Skip blank lines
        IF( HDRFLAG ) CYCLE        ! Skip after /END/
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Depending on source category, transfer line to temporary
!               fields.  In cases where blanks are allowed, do not use
!               STR2INT to prevent warning messages.
        CALL PARSLINE( LINE, MXCOL, SEGMENT )

        CFIP = SEGMENT( 1 )
        TSCC = SEGMENT( 2 )

!.............  Smart interpretation of SCC
        CALL FLTRNEG( TSCC )     ! Filter 0 and -9 to blank
        CALL PADZERO( TSCC )     ! Pad with zeros

!.............  Make sure that the co/st/cy code is an integer
        IF( .NOT. USEEXPGEO() .AND.&
        &    .NOT. CHKINT( CFIP ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Country/state/county ' //&
            &       'code is not an integer at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  If this record is in error, go to next iteration
        IF( EFLAG ) CYCLE

!.............  Adjust co/st/cy code
        CFIP = ADJUSTR( CFIP )
        CALL PADZERO( CFIP )

        N = N + 1
        IF( N .GT. NLINES ) CYCLE  ! Ensure no overflow

!.............  Store case-indpendent fields
        INDXTA ( N ) = N
        CFIPTA ( N ) = CFIP
        CSCCTA ( N ) = TSCC

!.............  Store sorting criteria for source.
!.............  NOTE - if point sources are added, make sure that
!               TSCC is justified correctly.

!.............  For point sources
        IF ( CATEGORY == 'POINT' ) THEN

            PLT      = SEGMENT( 3 )
            CHARS(1) = SEGMENT( 4 )
            CHARS(2) = SEGMENT( 5 )
            CHARS(3) = SEGMENT( 6 )
            CHARS(4) = SEGMENT( 7 )
            CHARS(5) = SEGMENT( 8 )

!.................  Store sorting criteria as right-justified in fields
            CALL BLDCSRC( CFIP, PLT, CHARS(1),&
            &              CHARS(2), CHARS(3), CHARS(4),&
            &              CHARS(5), POLBLNK3, CSRCALL   )

            CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC

!.............  For area and mobile sources
        ELSE

            CSRCTA( N ) = CFIP // TSCC

        END IF

    END DO      ! End of loop on I for reading in NHAPEXCLUDE file

!.........  Reset number of cross-reference entries in case some were dropped
    NXREF = N

!.........  Write errors for problems with input
    IF( NXREF .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No valid non-HAP inclusions or exclusions entries!'
        CALL M3MSG2( MESG )

    ELSEIF( NXREF .GT. NLINES ) THEN
        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //&
        &    'storing non-HAP inclusions/exclusions file was',NLINES,&
        &    CRLF() // BLANK10 // 'but actually needed', NXREF
        CALL M3MSG2( MESG )

    ENDIF

!.......  Check for errors reading XREF file, and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading non-HAP exclusions/inclusions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    IF( NHAP_EXCL ) THEN
        MESG = 'Processing non-HAP exclusions file...'
    ELSE
        MESG = 'Processing non-HAP inclusions file...'
    ENDIF

    CALL M3MSG2( MESG )

    CALL SORTIC( NXREF, INDXTA, CSRCTA )

!.........  Group cross-reference data into tables for different groups
    CALL XREFTBL( 'NONHAP', NXREF )

!.........  Deallocate cross-reference sorting arrays
    DEALLOCATE( INDXTA, CFIPTA, CSCCTA, CSRCTA )

!.........  Rewind file
    REWIND( FDEV )

    RETURN

!.........  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of non-HAP' // CRLF() // BLANK5 //&
    &       'inclusions/exclusions file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDXCLUDE
