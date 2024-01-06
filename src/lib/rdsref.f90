
SUBROUTINE RDSREF( FDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     Reads the speciation cross-reference file for any source category.  It
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
    !
    !       Revised 8/2016 by Carlie Coats, UNC IE, to treat XREF fractional
    !       speciation profiles

    !       Version 11/2023 by CJC:  USE M3UTILIO".f90" free source format, and
    !       related changes
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
    USE MODXREF, ONLY:  INDXTA, CSRCTA, CSCCTA, CMACTA, CISICA, CSPRNA, ISPTA, XDUPCHK

    !.........  This module contains the information about the source category
    USE MODINFO,  ONLY: CATEGORY, NIPPA, EANAM, LSCCEND

    USE MODLISTS, ONLY: NINVIFIP, INVCFIP

    USE MODMBSET, ONLY: NREFC, MCREFIDX, NINVC, MCREFSORT

    !.........  This module contains the information about the source category
    IMPLICIT NONE

    !...........   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV       ! cross-reference file unit no.

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

    !...........   Local parameters
    INTEGER,            PARAMETER :: MXCOL   = 15         !    !  add col "M" for pt-src sref-fracs
    CHARACTER(FIPLEN3), PARAMETER :: FIPZERO = REPEAT('0', FIPLEN3 )     !  buffer for zero Cy/St/Co code
    CHARACTER(SCCLEN3), PARAMETER :: SCCZERO = REPEAT('0', SCCLEN3 )     !  buffer for zero SCC
    CHARACTER(SICLEN3), PARAMETER :: SICZERO = REPEAT('0', SICLEN3 )     !  buffer for zero SIC
    CHARACTER(6),       PARAMETER :: LOCCATS( 3 ) = (/ 'AREA  ', 'MOBILE', 'POINT ' /)
    CHARACTER(16),      PARAMETER :: PROGNAME = 'RDSREF'     ! program name

    !...........   Sorted pollutant/emission type names
    INTEGER               INDXP  ( NIPPA )     !  sort index for pols/etypes
    CHARACTER(NAMLEN3) :: SRTINAM( NIPPA )     !  sorted pol/etype names

    !...........   Array of point source plant characeristics
    CHARACTER(CHRLEN3) CHARS( 5 )

    !...........   Array for parsing list-formatted inputs
    CHARACTER(50)          SEGMENT( MXCOL )

    !...........   Other local variables
    INTEGER         I, J, J1, J2, K, L, M, NN, N        !  counters and indices

    INTEGER         IDUM        !  dummy integer
    INTEGER         IOS         !  i/o status
    INTEGER         IREC        !  record counter
    INTEGER         JS          !  position of SCC in source chars in x-ref file
    INTEGER         JSPC        !  tmp index to master pollutant/etype list
    INTEGER         LPCK        !  length of point definition packet
    INTEGER         NCP         !  input point source header parm
    INTEGER         XDEV        !  unit no. of for cross-refrence county file (MCXREF)
    INTEGER         NCNTY       !  no of inv counties per ref county
    INTEGER         NLINES      !  number of lines
    INTEGER         NXREF       !  number of valid x-ref entries
    INTEGER         NCOMBO, NFRACS          !  numbers of "COMBO" and fractional-profile entries

    LOGICAL      :: EFLAG = .FALSE.       !  true: error found
    LOGICAL      :: PFLAG = .FALSE.       !  true: tmp pol-spec rec skipped
    LOGICAL      :: REFLAG = .FALSE.      !  true: use of inv-ref county mapping
    LOGICAL      :: SKIPPOL = .FALSE.     !  true: pol-spec rec skipped in x-ref
    LOGICAL         SKIPREC               !  true: record skipped in x-ref file
    LOGICAL         DUPCHECK              !  true: duplicate-check in XREFTBL

    REAL            FRAC

    CHARACTER(5)       CPOS         !  temporary pol code or position

    CHARACTER(300)     LINE         !  line buffer
    CHARACTER(300)     MESG         !  message buffer
    CHARACTER(NAMLEN3) CPOA         !  temporary pollutant/emis type name
    CHARACTER(MACLEN3) CMCT         !  temporory MACT code
    CHARACTER(SICLEN3) CSIC         !  temporary SIC
    CHARACTER(ALLLEN3) CSRCALL      !  buffer for source char, incl pol
    CHARACTER(FIPLEN3) CFIP         !  buffer for CFIPS code
    CHARACTER(SCCLEN3) TSCC         !  temporary SCC
    CHARACTER(PLTLEN3) PLT          !  tmp plant ID
    CHARACTER(SPNLEN3) SPCODE       !  tmp speciation profile code

    REAL, ALLOCATABLE :: SFRACA( : )

    !***********************************************************************
    !   begin body of subroutine RDSREF

    !.........  Ensure that the CATEGORY is valid
    I = INDEX1( CATEGORY, 3, LOCCATS )

    IF( I .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: category "' // TRIM( CATEGORY ) // '" is not valid '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Sort the actual list of pollutant/emis type names and store it
    DO I = 1, NIPPA
        INDXP( I ) = I
    END DO

    CALL SORTIC( NIPPA, INDXP, EANAM )

    DO I = 1, NIPPA
        J = INDXP( I )
        SRTINAM( I ) = EANAM( J )
    END DO

    !.........  Write status message
    MESG = 'Reading speciation cross-reference file...'
    CALL M3MSG2( MESG )

    !.........  Get the number of lines in the file
    NN = GETFLINE( FDEV, 'Speciation cross reference file' )

    !.........  Use reference to inventory county mapping file
    MESG   = 'Use the county cross-reference mapping file'
    REFLAG = ENVYN( 'USE_REF_COUNTY_MAP_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_REF_COUNTY_MAP_YN"', 2 )
    END IF
    IF( REFLAG ) THEN

        !.............  Read and store ref-inv counties mapping
        XDEV = PROMPTFFILE( 'Enter logical name for county x-reference file',   &
                            .TRUE., .TRUE., 'MCXREF', PROGNAME )
        CALL RDMXREF( XDEV, NINVIFIP, INVCFIP )

        IREC   = 0
        NLINES = NN
        DO I = 1, NN

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                WRITE( MESG,94010 ) 'I/O error', IOS,               &
                    'reading speciation x-ref file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            !.................  Skip blank lines or comments
            IF( BLKORCMT( LINE ) ) CYCLE

            !.................  Read point source header information
            J = INDEX( LINE, PDEFPCKT )             ! can be in middle of file
            IF( J .GT. 0 ) CYCLE

            CALL PARSLINE( LINE, MXCOL, SEGMENT )
            CFIP   = SEGMENT( 4 )
            IF( LEN_TRIM( CFIP ) > 0 ) CALL  PADZERO( CFIP )

            !.................  Increase no of entries by no of inv per ref county
            K = FINDC( CFIP, NREFC, MCREFIDX( :,1 ) )
            IF( K > 0 ) THEN
                L = 0
                DO M = 1, NINVC
                    IF( CFIP == MCREFSORT( M,2 ) ) L = L + 1
                END DO
                NLINES = NLINES + L
            ELSE IF( LEN_TRIM( CFIP ) > 0 ) THEN
                WRITE( MESG,94010 ) 'WARNING: ' // TRIM(CFIP) //        &
                    ' county is not a reference counties at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

        END DO

    ELSE

        NLINES = NN       ! default NLINES from GSREF input file

    END IF

    !.........  Allocate memory for unsorted data used in all source categories
    ALLOCATE( CSPRNA( NLINES ),     &
               ISPTA( NLINES ),     &
              CSCCTA( NLINES ),     &
              CSRCTA( NLINES ),     &
              CMACTA( NLINES ),     &
              CISICA( NLINES ),     &
              INDXTA( NLINES ),     &
              SFRACA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSPRNA:SPRFACA', PROGNAME )
    SFRACA(:) = -9999.0

    !.........  Set up constants for loop.
    !.........  Length of point definition packet, plus one
    LPCK = LEN_TRIM( PDEFPCKT ) + 1

    !.........  Put file read pointer at top of file
    REWIND( FDEV )

    !.........  Initialize character strings
    CHARS   = ' '       ! array
    SEGMENT = ' '       ! array

    !.........  Read lines and store unsorted data for the source category of
    !           interest
    NCOMBO = 0
    NFRACS = 0
    IREC   = 0
    N      = 0
    NCP    = 6            ! ORL and IDA default (4+2)
    JS     = 6            ! ORL and IDA default (4+2)
    DO I = 1, NN

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,               &
                'reading speciation x-ref file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.............  Skip blank lines or comments
        IF( BLKORCMT( LINE ) ) CYCLE

        J = INDEX( LINE, PDEFPCKT )     ! can be in middle of file
        L = LEN_TRIM( LINE )

            !.............  Read point source header information
        IF( J .GT. 0 ) THEN

            IF( L .GT. LPCK ) THEN
                READ( LINE( LPCK:L ), * ) NCP, JS

            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Incomplete point '//    &
                       'source definition packet at line', IREC
                CALL M3MSG2( MESG )

            END IF

            !...............  Adjust for FIPS code and Plant ID, which are always there
            NCP = NCP + 2
            IF( JS .GT. 0 ) JS = JS + 2

            CYCLE

            !...........  If not a header line, then it's a regular line.  The records
            !             that don't apply to this source category or to the current
            !             inventory will be filtered out by FLTRXREF
        ELSE

            !...............  Compare point source definition from header to inventory
            IF( CATEGORY .EQ. 'POINT' ) CALL CHKPTDEF( NCP, JS )

            CALL PARSLINE( LINE, MXCOL, SEGMENT )

            TSCC   = SEGMENT( 1 )
            SPCODE = SEGMENT( 2 )
            CPOA   = SEGMENT( 3 )
            CFIP   = SEGMENT( 4 )
            CMCT   = SEGMENT( 5 )
            CSIC   = SEGMENT( 6 )
            PLT    = SEGMENT( 7 )
            CHARS( 1:5 ) = SEGMENT( 8:12 )

            !...............  Increase no of entries by no of inv per ref county
            K = 0
            NCNTY = 0                ! default value for no ref-inv county mapping
            IF( REFLAG ) THEN
                IF( LEN_TRIM( CFIP ) > 0 ) CALL PADZERO( CFIP )
                K = FINDC( CFIP, NREFC, MCREFIDX( :,1 ) )
                IF( K > 0 ) THEN
                    L = 0
                    DO M = 1, NINVC
                        IF( CFIP == MCREFSORT( M,2 ) ) L = L + 1
                    END DO
                    NCNTY = NCNTY + L              ! count no of inv per ref
                END IF
            END IF

            IF( NCNTY == 0 ) NCNTY = 1               ! default value = 1

            DO M = 1, NCNTY

                !.................  Update ref to inv county
                IF( K > 0 )  THEN
                    J1 = STR2INT( MCREFIDX( K,2 ) ) + M - 1
                    CFIP = MCREFSORT( J1,1 )                   ! mapped inventory county
                END IF

                !.........  Adjust these for proper sorting and matching with profiles file.
                SPCODE = ADJUSTR( SPCODE )
                CPOA   = ADJUSTL( CPOA   )

                !.........  Skip all point entries for nonpoint sectors
                IF ( CATEGORY /= 'POINT' .AND. PLT /= ' ' ) CYCLE

                !.........  Post-process x-ref information to scan for '-9', pad
                !           with zeros, compare SCC version master list, and compare
                !           pollutant/emission type name with master list.
                CALL FLTRXREF( CFIP, CSIC, TSCC, CPOA, CMCT,        &
                               IDUM, IDUM, JSPC, PFLAG, SKIPREC )

                SKIPPOL = ( SKIPPOL .OR. PFLAG )

                !.......... Filter the case where the pollutant code is not present
                IF( CPOA .EQ. ' ' ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG, 94010 )                                &
                           'ERROR: Skipping cross-reference entry ' //  &
                           'at line', IREC,                             &
                           'because of missing pollutant.'
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

                IF( SKIPREC ) CYCLE                  ! Skip this record

                !.........  Write pollutant position to character string
                WRITE( CPOS, '(I5.5)' ) JSPC

                !.........  Check for bad cross-reference code
                IF( SPCODE .EQ. ' ' ) THEN
                    WRITE( MESG, 94010 )                                    &
                      'WARNING: Skipping blank profile code in cross-'//    &
                      'reference file at line ', IREC
                    CALL M3MESG( MESG )
                    CYCLE

                END IF

                !.........  If SIC is defined, make sure SCC is not and fill SCC
                !           with SIC value and special identifier
                IF( CSIC /= SICZERO .AND. TSCC /= SCCZERO ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Both SCC and SIC ' //     &
                           'values are given at line', I, '.'//CRLF() //    &
                           BLANK10 // 'Only the SCC will be used ' //       &
                           'for this cross-reference entry.'
                    CALL M3MSG2( MESG )
                    CSIC = SICZERO

                END IF

                !.................  Increment count of valid x-ref entries and check it
                N = N + 1
                IF( N .GT. NLINES ) CYCLE                  ! Ensure no overflow

                !.................  Store case-specific fields from cross reference
                CSRCALL = ' '
                SELECT CASE( CATEGORY )

                  CASE( 'AREA' )

                    CALL BLDCSRC( CFIP, TSCC, CHRBLNK3,             &
                                  CHRBLNK3, CHRBLNK3, CHRBLNK3,     &
                                  CHRBLNK3, POLBLNK3, CSRCALL   )

                    CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // CMCT // CSIC // CPOS

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

                    CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // CMCT // CSIC // CPOS

                  CASE( 'POINT' )

                    !.....................  Store sorting criteria as right-justified in fields
                    CALL BLDCSRC( CFIP, PLT, CHARS(1),              &
                                  CHARS(2), CHARS(3), CHARS(4),     &
                                  CHARS(5), POLBLNK3, CSRCALL   )

                    CSRCTA( N ) = CSRCALL( 1:SRCLEN3 ) // TSCC // CMCT // CSIC // CPOS

                END SELECT

                    !.................  Store case-indpendent fields from cross-reference
                INDXTA( N ) = N
                ISPTA ( N ) = JSPC                        ! Save index to EANAM or zero
                CSCCTA( N ) = TSCC
                CMACTA( N ) = CMCT
                CISICA( N ) = CSIC
                CSPRNA( N ) = SPCODE

                FRAC        = STR2REAL( SEGMENT( 13 ) )
                IF ( FRAC .GE. 0.0 ) THEN
                    NFRACS      = NFRACS + 1
                    SFRACA( N ) = FRAC
                END IF

                IF ( ADJUSTL( SPCODE ) .EQ. 'COMBO' ) NCOMBO = NCOMBO + 1

            END DO      ! End of loop on no of mapped inv-ref counties

        END IF        !  This line matches source category of interest

    END DO          ! End of loop on I for reading in speciation x-ref file

    !.........  Reset number of cross-reference entries in case some were dropped
    NXREF = N

    !.........  Write warning message for pollutants in cross-reference that are
    !           not in master list
    IF( SKIPPOL ) THEN
        MESG = 'Pollutant-specific entries in the speciation ' //       &
               'cross-reference file have ' // CRLF() // BLANK10 //     &
               'been skipped.'
        CALL M3WARN( PROGNAME, 0, 0, MESG )
    END IF

    WRITE( MESG, '( A, I10 )' ) 'Number of inital xrefs:', NXREF
    CALL M3MSG2( MESG )

    WRITE( MESG, '( A, I10 )' ) 'Number of COMBO xrefs:', NCOMBO
    CALL M3MSG2( MESG )

    WRITE( MESG, '( A, I10 )' ) 'Number of FRACS xrefs:', NFRACS
    CALL M3MSG2( MESG )

    IF( NXREF .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No valid speciation cross-reference entries'
        CALL M3MSG2( MESG )

    ELSE IF( NXREF .GT. NLINES ) THEN
        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //     &
               'storing speciation cross-reference was', NLINES,    &
               CRLF() // BLANK10 // 'but actually needed', NXREF
        CALL M3MSG2( MESG )

    ELSE IF( NCOMBO.GT. 0 .AND. NFRACS .GT. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Both "COMBO" and fractional-profile entries in XREF'
        CALL M3MSG2( MESG )

    END IF

    !.......  Check for errors reading XREF file, and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading speciation cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MESG = 'Processing speciation cross-reference file...'
    CALL M3MSG2( MESG )

    !.........  Sort speciation cross-reference entries. Since CPOS was used in
    !           building CSRCTA, and CPOS will equal "0" when the x-ref entry is
    !           not pollutant/emistype-specific, the these entries will
    !           always appear first.  This is necessary for the table-generating
    !           subroutines.
    CALL SORTIC( NXREF, INDXTA, CSRCTA )

    IF ( NFRACS .GT. 0 ) THEN
        XDUPCHK = .FALSE.           !    !  don't need dup-checks in XREFTBL()
        CALL XFRACTBL( NLINES, NXREF, INDXTA, CSRCTA, CSPRNA, SFRACA )
        WRITE( MESG, '( A, I10 )' ) 'Number of final xrefs:', NXREF
        CALL M3MSG2( MESG )
    END IF

    CALL XREFTBL( 'SPECIATION', NXREF )

    !.........  Deallocate other temporary unsorted arrays
    DEALLOCATE( CSCCTA, ISPTA, CMACTA, CISICA, CSRCTA, CSPRNA, INDXTA, SFRACA )

    !.........  Rewind file
    REWIND( FDEV )

    RETURN

    !.........  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //              &
           CRLF() // BLANK5 // 'Check format of speciation cross-reference file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSREF

