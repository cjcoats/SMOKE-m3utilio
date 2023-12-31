
SUBROUTINE RDTREF( FDEV, FFORMAT )

    !***********************************************************************
    !  subroutine body starts at line 150
    !
    !  DESCRIPTION:
    !     Reads the temporal cross-reference file for any source category.  It
    !     allocates memory (locally) for reading the unsorted x-refs. It sorts the
    !     x-refs for processing. It allocates memory for the appropriate x-ref
    !     tables and populates the tables (passed via modules).
    !
    !  PRECONDITIONS REQUIRED:
    !     File unit FDEV already is opened... MORE
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 1/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !****************************************************************************/
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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

    !.......  MODULES for public variables
    !.......   This module is for cross reference tables
    USE MODXREF, ONLY: CSRCTA, ISPTA, CSCCTA, INDXTA, MPRNA, WPRNA, DPRNA

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NIPPA, EANAM

    USE MODTMPRL, ONLY: METPROFLAG, METPROTYPE

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: FDEV                  ! x-ref file unit no.
    CHARACTER(*), INTENT(OUT) :: FFORMAT               ! format of input file

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETNLIST
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local parameters
    INTEGER, PARAMETER :: AREATYP  = 1
    INTEGER, PARAMETER :: MOBILTYP = 2
    INTEGER, PARAMETER :: POINTTYP = 3
    INTEGER, PARAMETER :: MXTCOL   = 15

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDTREF'     ! program name
    CHARACTER(6) , PARAMETER :: LOCCATS( 3 ) = (/ 'AREA  ', 'MOBILE', 'POINT ' /)

    !.......   Sorted pollutant/emission type names
    INTEGER               INDXP  ( NIPPA )     !  sorting index for pol/etyp
    CHARACTER(NAMLEN3) :: SRTINAM( NIPPA )     !  sorted pol/act names

    !.......   Array of point source plant characeristics
    CHARACTER(CHRLEN3) CHARS( 5 )

    !.......   Array for reading temporal x-ref fields
    CHARACTER(20)      SEGMENT( MXTCOL )

    !.......   Other local variables
    INTEGER         I, J, J1, J2, J3, K, L, N        !  counters and indices

    INTEGER         COD         !  temporary pollutant/emission type code
    INTEGER         FIP         !  temporary FIPS code
    INTEGER         IDIU        !  temporary diurnal profile code
    INTEGER      :: IDUM = 0    !  tmp dummy integer
    INTEGER         IMON        !  temporary monthly profile code
    INTEGER         IOS         !  i/o status
    INTEGER         IWEK        !  temporary weekly profile code
    INTEGER         IREC        !  record counter
    INTEGER      :: JS = 0      !  position of SCC in source chars in x-ref file
    INTEGER         JSPC        !  tmp index to master pollutant/etype list
    INTEGER         LINTYPE     !  temporary source category code
    INTEGER         LPCK        !  length of point definition packet
    INTEGER      :: NCP = 0     !  input point source header parm
    INTEGER         NFIELD      !  tmp number of fields in LINE
    INTEGER         NLINES      !  number of lines
    INTEGER         NREF        !  number of x-ref entries before filtering
    INTEGER         NXREF       !  number of valid x-ref entries
    INTEGER         RDT         !  temporary road class code
    INTEGER         TMON        !  temporary met-based monthly profile code
    INTEGER         TDAY        !  temporary met-based daily profile code
    INTEGER         THRS        !  temporary met-based hourly profile code
    INTEGER         VTYPE       !  temporary vehicle type number

    LOGICAL      :: EFLAG = .FALSE.       !  true: error occurred
    LOGICAL      :: PFLAG = .FALSE.       !  true: pol/act-spec entries skipped
    LOGICAL      :: SKIPREC = .FALSE.     !  true: skip this x-ref entry

    CHARACTER          SCC1         !  1st character of SCC
    CHARACTER(5)       CPOS         !  tmp sorted position of pol/act
    CHARACTER(300)     LINE         !  line buffer
    CHARACTER(300)     MESG         !  message buffer

    CHARACTER(SICLEN3) :: CDUM = '0'     !  dummy character field for SIC
    CHARACTER(MACLEN3) :: CDUM2= '0'     !  dummy charecter field for MACT
    CHARACTER(LNKLEN3) CLNK         !  temporary link code
    CHARACTER(ALLLEN3) CSRCALL      !  buffer for source char, incl pol/act
    CHARACTER(FIPLEN3) CFIP         !  buffer for CFIPS code
    CHARACTER(CNYLEN3) CNY          !  tmp county string
    CHARACTER(FIPLEN3) FIPZERO      !  buffer for zero FIPS code
    CHARACTER(SCCLEN3) TSCC         !  temporary SCC
    CHARACTER(SCCLEN3) SCCZERO      !  buffer for zero SCC
    CHARACTER(PLTLEN3) PLT          !  tmp plant ID
    CHARACTER(NAMLEN3) CPOA         !  temporary pollutant/emission type
    CHARACTER(RWTLEN3) CRWT         !  roadway type no.
    CHARACTER(VIDLEN3) CVID         !  vehicle type ID no.

    !***********************************************************************
    !   begin body of subroutine RDTREF

    !.......  Ensure that the CATEGORY is valid
    I = INDEX1( CATEGORY, 3, LOCCATS )

    IF( I .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: category "' // TRIM( CATEGORY ) //  &
               '" is not valid in routine ' // PROGNAME
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    ENDIF

    !.......  Set up zero strings for FIPS code of zero and SCC code of zero
    FIPZERO = REPEAT( '0', FIPLEN3 )
    SCCZERO = REPEAT( '0', SCCLEN3 )

    !.......  Sort the actual list of pollutant/emission type names and store it
    DO I = 1, NIPPA
        INDXP( I ) = I
    END DO

    CALL SORTIC( NIPPA, INDXP, EANAM )

    DO I = 1, NIPPA
        J = INDXP( I )
        SRTINAM( I ) = EANAM( J )
    END DO

    !.......  Write status message
    MESG = 'Reading temporal cross-reference file...'
    CALL M3MSG2( MESG )

    !.......  Set up constants for loop...

    !.......  Length of point definition packet, plus one
    LPCK = LEN_TRIM( PDEFPCKT ) + 1

    !.......  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'Temporal cross reference file' )

    !.......  First pass through file.  Determine format and count the number
    !         of lines matching SCC list and pol/act list.  Do this so that
    !         we know how much memory to allocate for the unsorted, unprocessed
    !         arrays.
    IREC = 0
    NREF = 0
    NCP  = 6            ! ORL and IDA default (4+2)
    JS   = 6            ! ORL and IDA default (4+2)
    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )             &
                'I/O error', IOS,           &
                'reading TEMPORAL XREF file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        L = LEN_TRIM( LINE )
        NFIELD = GETNLIST( L,LINE )

        !.......  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) THEN
            CYCLE

        !.......  Auto check for LIST or EPS2.0 formated temporal x-ref...
        !.......  If header is found, read point source header information
        ELSE IF( INDEX( LINE, PDEFPCKT ) .GT. 0 ) THEN

            IF( L .GT. LPCK ) THEN
                READ( LINE( LPCK:L ), * ) NCP, JS

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Incomplete point '//    &
                       'source definition packet at line', IREC
                CALL M3MSG2( MESG )

            END IF

            !.......  Adjust for FIPS code and Plant ID, which are always there
            NCP = NCP + 2
            IF( JS .GT. 0 ) JS = JS + 2

            CYCLE

        !.......  The source-formatted will have only 3 columns
        ELSE IF( NFIELD .LE. 3 ) THEN

            FFORMAT = 'SOURCE'
            NREF = NLINES         ! equals the number of sources
            EXIT                  ! Exit from read loop

        ELSE

            !.......  Compare point source definition from header to inventory
            IF( CATEGORY .EQ. 'POINT' ) CALL CHKPTDEF( NCP, JS )

            FFORMAT = 'STANDARD'

            CALL PARSLINE( LINE, MXTCOL, SEGMENT )

            !.......  Make sure SCC is set to SCCZERO if it is missing
            TSCC = ADJUSTL( SEGMENT( 1 ) )
            CALL FLTRNEG( TSCC )
            CALL PADZERO( TSCC )

            CPOA = SEGMENT( 5 )               ! pollutant/emission type name
            CFIP = SEGMENT( 6 )               ! country/state/county code

            !.......  Post-process x-ref information to scan for '-9', pad
            !         with zeros, compare SCC version master list, compare
            !         SIC version to master list, and compare pol/act name
            !         with master list.
            CALL FLTRXREF( CFIP, CDUM, TSCC, CPOA, CDUM2,&
                           IDUM, IDUM, JSPC, PFLAG, SKIPREC )

            !.......  Skip lines that are not valid for this inven and src cat
            IF( SKIPREC ) CYCLE

            NREF = NREF + 1

        END IF         ! End format of temporal x-ref file

    END DO             ! End first pass through file

    REWIND( FDEV )

    !.......  Allocate memory for unsorted data used in all source categories and
    !         input formats
    ALLOCATE( MPRNA( NREF ),            &
              WPRNA( NREF )             &
              DPRNA( NREF ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MPRNA...DPRNA', PROGNAME )

    !.......  Special section for list-formatted cross-reference file (used only
    !         with EMS-95 inputs).
    IF( FFORMAT .EQ. 'SOURCE' ) THEN

        DO I = 1, NLINES

            READ( FDEV, *, IOSTAT=IOS ) IMON, IWEK, IDIU

            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )             &
                    'I/O error', IOS,           &
                    'reading TEMPORAL XREF file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            MPRNA( I ) = IMON
            WPRNA( I ) = IWEK
            DPRNA( I ) = IDIU

        END DO      !  End of loop for reading list-formatted xref file

    END IF

    !.......  Check for errors from reading either format
    IF( EFLAG ) THEN
        MESG = 'Problem reading temporal cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Report format for XREF file
    ELSE
        MESG = 'NOTE: File read in as ' // FFORMAT // ' format.'
        CALL M3MSG2( MESG )

    END IF

    !.......  Leave routine for list-formatted file because no grouping needed
    IF( FFORMAT .EQ. 'SOURCE' ) RETURN

    !.......  FFORMAT = 'SOURCE' DOES NOT APPLY AFTER THIS POINT

    !.......  Allocate memory for unsorted data used in all source categories
    ALLOCATE( INDXTA( NREF ),               &
               ISPTA( NREF ),               &
              CSCCTA( NREF ),               &
              CSRCTA( NREF ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSRCTA', PROGNAME )

    !.......  Initialize character strings
    CHARS   = ' '      ! array
    SEGMENT = ' '      ! array

    !.......  Second pass through file: read lines and store unsorted data for
    !         the source category of interest
    IREC = 0
    N    = 0
    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )         &
                'I/O error', IOS,       &
                'reading temporal cross-reference file at line',IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        J = INDEX( LINE, PDEFPCKT )         ! can be in middle of file
        L = LEN_TRIM( LINE )

        !.......  Skip point source header information and filter the line again
        !         for storing info
        IF( J .LE. 0 ) THEN

            !.......  Separate line from file into fields
            CALL PARSLINE( LINE, MXTCOL, SEGMENT )

            !.......  Make sure SCC is set to SCCZERO if it is missing
            TSCC = ADJUSTL( SEGMENT( 1 ) )
            CALL FLTRNEG( TSCC )
            CALL PADZERO( TSCC )

            CPOA = SEGMENT( 5 )               ! pollutant/emission type name
            CFIP = SEGMENT( 6 )               ! country/state/county code
            PLT  = SEGMENT ( 7 )

            !.......  Skip all point entries for nonpoint sectors
            IF ( CATEGORY /= 'POINT' .AND. PLT /= ' ' ) CYCLE

            !.......  Post-process x-ref information to scan for '-9', pad
            !         with zeros, compare SCC version master list, compare
            !         SIC version to master list, and compare pol/act name
            !         with master list.
            CALL FLTRXREF( CFIP, CDUM, TSCC, CPOA, CDUM2,       &
                           IDUM, IDUM, JSPC, PFLAG, SKIPREC )

            !.......  Skip lines that are not valid for this inven and src cat
            IF( SKIPREC ) CYCLE

            !.......  Give warning for entries that are plant-specific but not
            !         county-specific, since these are not supported (this can
            !         sometimes happen with Canadian records)
            CNY = CFIP( STALEN3+1:FIPLEN3 )
            IF ( CNY == '000' .AND. PLT /= ' ' ) THEN
                WRITE( MESG,94010 ) 'WARNING: temporal cross-reference'//   &
                   ' entry has plant specified but no '// CRLF() //         &
                   BLANK10 // 'county specified at line',IREC,'. '//        &
                   'Program will treat as a state or state-SCC entry.'
                CALL M3MESG( MESG )
            END IF

            !.......  Check for integers for temporal profile numbers
            IF( .NOT. CHKINT( SEGMENT( 2 ) ) .OR.       &
                .NOT. CHKINT( SEGMENT( 3 ) ) .OR.       &
                .NOT. CHKINT( SEGMENT( 4 ) )      ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )'ERROR: temporal profile '//     &
                       'code(s) at line', IREC, 'of temporal'//     &
                       CRLF() // BLANK10 // 'cross-reference file are non-integer.'
                CALL M3MESG( MESG )
                CYCLE

            !.......  Convert temporal profile numbers
            ELSE
                IMON = STR2INT( SEGMENT( 2 ) )
                IWEK = STR2INT( SEGMENT( 3 ) )
                IDIU = STR2INT( SEGMENT( 4 ) )
                TMON = STR2INT( SEGMENT( 13 ) )                ! met-based monthly profile numbers
                TDAY = STR2INT( SEGMENT( 14 ) )                ! met-based daily profile numbers
                THRS = STR2INT( SEGMENT( 15 ) )                ! met-based hourly profile numbers

                !.......  Check for bad cross-reference code
                IF( IMON .LE. 0 .OR. IWEK .LE. 0 .OR.               &
                    IDIU .LE. 0 ) THEN
                    WRITE( MESG, 94010 )                            &
                     'WARNING: Skipping bad cross-reference code '//&
                     'at line. Values are:' // CRLF() // BLANK16 // &
                     'Monthly=', IMON, 'Weekly=', IWEK, 'Diurnal=', &
                     IDIU
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

            END IF

            !.......  Write pol/act position to a character string
            WRITE( CPOS, '(I5)' ) JSPC

            N = N + 1
            IF( N .GT. NREF ) CYCLE              ! Ensure no overflow

            CSRCALL = ' '
            !.......  Store sorting criteria as right-justified in fields
            !.......  For mobile sources, retrieve link from character field
            !         and extract road class and vehicle type from SCC
            !.......  For point sources, retrieve plant + characteristics
            SELECT CASE( CATEGORY )

              CASE( 'AREA' )

                CALL BLDCSRC( CFIP, CHRBLNK3, CHRBLNK3, CHRBLNK3,   &
                              CHRBLNK3, CHRBLNK3, CHRBLNK3,         &
                              POLBLNK3, CSRCALL )

              CASE( 'MOBILE' )

                ! M Houyoux note: TSCC has been put in here instead of road type
                !     and link has been removed.  These were breaking the county-SCC specific
                !     assignments by setting CNFIP in xreftbl.f to be non-blank and not the SCC.
                !     However, this change breaks link-specific profile assignments, which
                !     are not likely to be used anyway.  I suggest that we just remove
                !     link-specific assignments from the documentation for Spcmat.

                CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,   &
                              CHRBLNK3, CHRBLNK3, CHRBLNK3,     &
                              POLBLNK3, CSRCALL )

              CASE( 'POINT' )

                !.......  Store string source characteristics
                CHARS( 1:5 ) = SEGMENT( 8:12 )

                CALL BLDCSRC( CFIP, PLT, CHARS(1),              &
                              CHARS(2), CHARS(3), CHARS(4),     &
                              CHARS(5), POLBLNK3, CSRCALL )

            END SELECT

            !.......  Store fields
            CSCCTA( N ) = TSCC
            CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC // CPOS

            !.......  Store case-indpendent fields
            INDXTA( N ) = N
            ISPTA ( N ) = JSPC             ! Save index to original EANAM or zero
            MPRNA ( N ) = IMON
            WPRNA ( N ) = IWEK
            DPRNA ( N ) = IDIU

            !.......  Reset org profID to met-based profile ID
            IF( TMON > 0 ) THEN
                MPRNA( N ) = 99999
                METPROFLAG = .TRUE.
                METPROTYPE = 'MONTHLY'
            END IF

            IF( TDAY > 0 ) THEN
                WPRNA( N ) = 99999
                METPROFLAG = .TRUE.
                METPROTYPE = 'DAILY'
            END IF

            IF( THRS > 0 ) THEN
                DPRNA( N ) = 99999
                METPROFLAG = .TRUE.
                METPROTYPE = 'HOURLY'
            END IF

        END IF      !  This line matches source category of interest

    END DO          ! End of loop on I for reading in temporal x-ref file

    !.......  Set actual number of cross-reference entries
    NXREF = N

    !.......  Check if no cross-reference entries were found
    IF( NXREF .EQ. 0 ) THEN
        MESG = 'ERROR: No valid temporal cross-reference entries were found'
        CALL M3MSG2( MESG )
        EFLAG = .TRUE.
    END IF

    !.......  Check for errors reading cross-reference file, and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading temporal cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    CALL M3MSG2( 'Processing temporal cross-reference file...' )

    !.......  Sort temporal cross-reference entries. Since CPOS was used in
    !         building CSRCTA, and CPOS will equal "0" when the x-ref entry is
    !         not pol/act-specific, the non-pol/act-specific entries will
    !         always appear first.  This is necessary for the table-generating
    !         subroutines.

    CALL SORTIC( NXREF, INDXTA, CSRCTA )

    CALL XREFTBL( 'TEMPORAL', NXREF )

    !.......  Deallocate temporary unsorted arrays
    DEALLOCATE( INDXTA, ISPTA, CSCCTA, CSRCTA, MPRNA, WPRNA, DPRNA )

    !.......  Rewind file
    REWIND( FDEV )

    RETURN

    !.......  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //              &
           'Check format of temporal' // CRLF() // BLANK5 //    &
           'cross reference file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDTREF
