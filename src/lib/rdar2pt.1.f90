
SUBROUTINE RDAR2PT( FDEV, CDEV, LDEV )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!     Reads the area-to-point data in ASCII format and populates
!     the MODAR2PT arrays.  Sets up the cross-references and
!     calls the cross-reference assignments routine, which in
!     turn populates the by-source arrays that contains the
!     information needed for assigning the point source locations
!     to area sources.
!
!  PRECONDITIONS REQUIRED:
!     File unit FDEV already is opened
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     Subroutines: Models-3 subroutines
!     Functions: Models-3 functions
!
!  REVISION  HISTORY:
!       Created 11/02 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************
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
    USE MODXREF, ONLY: INDXTA, CFIPTA, CSCCTA, CSRCTA, IARPTA

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: NINVSCC, INVSCC, SCCDESC  ! Note - needed for reporting only

!.........  This module contains the arrays for the area-to-point x-form
    USE MODAR2PT, ONLY: MXROWA2P, NTABLA2P, NAR2PT, AR2PTABL,&
    &                    NA2PSCC, A2PSCC, AR2PT

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER , INTENT (IN) :: FDEV   ! area-to-point factors unit no.
    INTEGER , INTENT (IN) :: CDEV   ! SCC descriptions unit no.
    INTEGER , INTENT (IN) :: LDEV   ! log file unit no.

!...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: USEEXPGEO

!.............  Parameters
    INTEGER, PARAMETER :: NFIELDS = 10  ! no. input fields
    INTEGER, PARAMETER :: FBEG( NFIELDS ) =&
    &                    ( / 1 , 8, 11, 34, 39, 63,&
    &                        86, 97, 114, 126 / )
    INTEGER, PARAMETER :: FEND( NFIELDS ) =&
    &                    ( / 6 , 9, 32, 37, 61, 84,&
    &                        93, 109, 124, 139 / )

!...........   Local allocatable arrays...
!...........   For processing:
    INTEGER      , ALLOCATABLE :: IDXA2P  ( :,: )   ! sorting index
    INTEGER      , ALLOCATABLE :: NUMSCC  ( : )     ! no. SCCs per table
    TYPE( AR2PT ), ALLOCATABLE :: UNSRTA2P( :,: )   ! unsorted tables
    CHARACTER(FIPLEN3), ALLOCATABLE :: LOCFIP  ( : )   ! tmp FIPS codes
    CHARACTER(SCCLEN3), ALLOCATABLE :: AR2PTSCC( :,: ) ! SCCs per table

!...........   For reporting:
    INTEGER      , ALLOCATABLE :: SCCIDX  ( : )       ! sorting index
    INTEGER      , ALLOCATABLE :: SCCSECTN( : )       ! section no. for SCC
    CHARACTER(SCCLEN3), ALLOCATABLE :: INVSCCA( : ) ! unsorted SCCs


!...........   Local arrays
    CHARACTER(32) SEGMENT( NFIELDS )

!...........   Other local variables
    INTEGER         I, J, K, L, N     ! counters and indices

    INTEGER      :: CNT = 0        ! tmp record counter
    INTEGER      :: FIPCNT = 0     ! counter for no. entries per SCC/FIPs
    INTEGER         IOS            ! i/o status
    INTEGER      :: MXSCC = 0      ! max no. of SCCs per table
    INTEGER      :: NXREF = 0      ! number of cross-ref entries

    REAL         :: SUMTEST = 0.   ! value to check that factors sum to 1.

    LOGICAL      :: EFLAG  = .FALSE.  ! true: error found
    LOGICAL      :: WFLAG             ! true: convert lat-lons to Western hemisphere

    CHARACTER(256)        MESG    !  message buffer
    CHARACTER(512)     :: LINE    !  input line
    CHARACTER(FIPLEN3) :: CFIP    !  tmp char co/st/cy code
    CHARACTER(FIPLEN3) :: LFIP    !  previous co/st/cy FIPS code
    CHARACTER(SCCLEN3) :: TSCC    !  tmp SCC code

    CHARACTER(16) :: PROGNAME = 'RDAR2PT' ! program name

!***********************************************************************
!   begin body of subroutine RDAR2PT

!.........  Rewind area-to-point file
    IF ( FDEV .NE. 0 ) THEN
        REWIND( FDEV )
    ELSE
        MESG = 'Area-to-point factors file is not opened!'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Loop through the input file to determine how many header
!           lines there are and the maximum number of entries in
!           a single table
    CALL READ_ARTOPT_FILE( 'COUNT', MXSCC )

!.........  Allocate memory for unsorted input tables
!.........  Allocate memory for sorted input tables

    ALLOCATE(  IDXA2P( MXROWA2P, NTABLA2P ),&
    &         UNSRTA2P( MXROWA2P, NTABLA2P ),&
    &         AR2PTSCC( MXSCC,    NTABLA2P ),&
    &           NUMSCC( NTABLA2P ),&
    &           LOCFIP( MXROWA2P ),&
    &           NAR2PT( NTABLA2P ),&
    &         AR2PTABL( MXROWA2P, NTABLA2P ), STAT=IOS )
    CALL CHECKMEM( IOS, 'AR2PTABL', PROGNAME )

    IDXA2P = 0                ! array
    UNSRTA2P%FIP   = ' '      ! array
    UNSRTA2P%LAT   = BADVAL3  ! array
    UNSRTA2P%LON   = BADVAL3  ! array
    UNSRTA2P%ALLOC = 1.       ! array
    UNSRTA2P%NAME  = ' '      ! array
    AR2PTSCC = ' '            ! array
    NUMSCC = 0                ! array

!.........  Allocate memory for sorted input tables
    NAR2PT         = 0        ! array
    AR2PTABL%FIP   = ' '      ! array
    AR2PTABL%LAT   = BADVAL3  ! array
    AR2PTABL%LON   = BADVAL3  ! array
    AR2PTABL%ALLOC = 1.       ! array
    AR2PTABL%NAME  = ' '      ! array

!.........  Check if lat-lons should be converted to western hemisphere
    MESG = 'Western hemisphere flag'
    WFLAG = ENVYN( 'WEST_HSPHERE', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "WEST_HSPHERE"', 2 )
    END IF

!.........  Read and store contents of the file (both the SCC
!           entries as well as the tables).
!.........  In this call, populate INXA2P, UNSRTA2P, AR2PTSCC,
!           NUMSCC, and NAR2PT
    CALL READ_ARTOPT_FILE( 'STORE', MXSCC )

!.........  Sort index for finding sorted order for each table.
!.........  Compute expected number of cross-reference entries
    NXREF = 0
    DO N = 1, NTABLA2P

!.............  Transfer FIPS codes to local array
        LOCFIP = UNSRTA2P(:,N)%FIP    ! array

!.............  Sort for current table
        CALL SORTIC( NAR2PT(N), IDXA2P(1,N), LOCFIP )

!.............  Compute maximum expected x-ref entries
        LFIP = ' '
        DO I = 1, NAR2PT( N )
            J = IDXA2P( I,N )
            CFIP = UNSRTA2P( J,N )%FIP
            IF( LFIP .NE. CFIP ) THEN
                NXREF = NXREF + NUMSCC( N )
            END IF
            LFIP = CFIP
        END DO          ! end loop through rows in table

    END DO              ! end loop through tables

!........  Allocate memory for cross-reference arrays
    ALLOCATE( INDXTA( NXREF ),&
    &          CFIPTA( NXREF ),&
    &          CSCCTA( NXREF ),&
    &          CSRCTA( NXREF ),&
    &          IARPTA( NXREF,3 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IARPTA', PROGNAME )
    INDXTA = 0   ! array
    CFIPTA = ' ' ! array
    CSCCTA = ' ' ! array
    CSRCTA = ' ' ! array
    IARPTA = 0   ! array

!.........  Store sorted input tables.
!.........  Build arrays for giving to cross-reference routines for
!           grouping the cross-references.  Although, for now, there
!           will be only the FIPS//SCC group.
    CNT  = 0
    DO N = 1, NTABLA2P

!...........  Store sorted input tables
        DO I = 1, NAR2PT( N )
            J = IDXA2P( I,N )
            AR2PTABL( I,N ) = UNSRTA2P( J,N )
        END DO

!...........  Loop through SCCs for this section of the file
!             and create cross-referencing arrays
        DO K = 1, NUMSCC( N )

!...............  Reset number of entries per FIPS code
            LFIP = ' '

!...............  Loop through rows for this section of the file
            DO I = 1, NAR2PT( N )

                CFIP = AR2PTABL( I,N )%FIP
                TSCC = AR2PTSCC( K,N )

!..................  If this row is a new FIPS code
                IF( CFIP .NE. LFIP ) THEN
                    CNT = CNT + 1  ! count of cross-reference entries

!.......................  If count of entries is w/i dimensioned array
                    IF( CNT .LE. NXREF ) THEN
                        INDXTA( CNT )   = CNT
                        CFIPTA( CNT )   = CFIP
                        CSCCTA( CNT )   = TSCC
                        IARPTA( CNT,1 ) = N
                        IARPTA( CNT,2 ) = I
                        CSRCTA( CNT )   = CFIP // TSCC
                    END IF

!.......................  Reset number of entries per FIPS code
                    FIPCNT = 0

                END IF

!....................  Increment the number of entries per FIPS code
                FIPCNT = FIPCNT + 1

!....................  Store the number of entries for this FIPS/SCC
                IF( CNT .LE. NXREF ) IARPTA( CNT,3 ) = FIPCNT

!....................  Store previous FIPS code for next iteration
                LFIP = CFIP

            END DO   ! end loop through rows of current table
        END DO       ! end loop through SCCs of current table

    END DO           ! end loop through tables

!.........  Check if count exceeded maximum expected
    IF( CNT .GT. NXREF ) THEN
        WRITE( MESG,94010 )
    ELSE
        NXREF = CNT
    END IF

!.........  Sort area-to-point factors cross-reference
    CALL SORTIC( NXREF, INDXTA, CSRCTA )

!.........  Check that fractions for each FIPS/SCC combo sums to 1.
    DO I = 1, NXREF

        J    = INDXTA( I )
        CFIP = CFIPTA( J )
        TSCC = CSCCTA( J )

!............  Loop through entries for current FIPS/SCC and compute
!              sum of allocation factors
        N   = IARPTA( J,1 )
        CNT = IARPTA( J,2 ) - 1
        SUMTEST = 0.
        DO K = 1, IARPTA( J,3 )
            CNT = CNT + 1
            SUMTEST = SUMTEST + AR2PTABL( CNT,N )%ALLOC
        END DO

!............  If sum of allocation factors is outside allowable
!              range, then write error
        IF( SUMTEST .GT. 1.001 .OR.&
        &    SUMTEST .LT. 0.999      ) THEN

            EFLAG = .TRUE.
            WRITE( MESG, 94020 ) 'ERROR: Sum of factors for '//&
            &       'Co/St/Cy= '// CFIP//'and SCC= '// TSCC//&
            &       CRLF()// BLANK10 // 'is not equal to 1. '//&
            &       'Value is ', SUMTEST, '.'
            CALL M3MSG2( MESG )

        END IF

    END DO

!.........  Abort if errors found
    IF( EFLAG ) THEN
        MESG = 'Area-to-point factors file has bad data'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Call cross-referencing routine, which will also populate
!           the ARPT09 array (full FIPs and full SCC matches only)
    CALL XREFTBL( 'AR2PT', NXREF )

!.........  Report to the log file the SCCs, SCC descriptions, and
!           the section of the input file set for each SCC.
!.........  To create this report, we will create a fake INVSCC array
!           for the MODLISTS module, which will allow us to read
!           and assign the SCC descriptions for this report...

!.........  Sum SCC count
    NINVSCC = 0
    DO N = 1, NTABLA2P
        NINVSCC = NINVSCC + NUMSCC( N )
    END DO

    NA2PSCC = NINVSCC

!.........  Allocate temporary "inventory" SCC list from MODLISTS
!.........  Allocate unsorted arrays

    ALLOCATE( INVSCC( NINVSCC ),&
    &          A2PSCC( NA2PSCC ),&
    &          SCCIDX( NINVSCC ),&
    &         INVSCCA( NINVSCC ),&
    &        SCCSECTN( NINVSCC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SCCSECTN', PROGNAME )

!.........  Store unsorted SCCs
    J = 0
    DO N = 1, NTABLA2P
        DO K = 1, NUMSCC( N )
            J = J + 1
            SCCIDX  ( J ) = J
            INVSCCA ( J ) = AR2PTSCC( K,N )
            SCCSECTN( J ) = N
        END DO
    END DO

!.........  Sort SCC list, needed for reading SCC descriptions
    CALL SORTIC( NINVSCC, SCCIDX, INVSCCA )
    DO J = 1, NINVSCC
        INVSCC( J ) = INVSCCA( SCCIDX( J ) )
    END DO

    A2PSCC = INVSCC

!.........  Retrieve SCC descriptions
    CALL RDSCCDSC( CDEV )

!.........  Write header of report
    MESG = 'NOTE: The area-to-point factors file ' //&
    &       'included the following SCCs and '// CRLF()//&
    &       BLANK10 // 'section numbers'
    CALL M3MESG( MESG )
    WRITE( LDEV, 94010 ) ' '

    MESG = 'SCC                  Section   SCC Description'
    WRITE( LDEV, 94010 ) BLANK10 // TRIM( MESG )
    WRITE( LDEV, 94010 ) BLANK10 // REPEAT( '-', 70 )

!.........  Loop through SCCs and list them, section no.'s, and descriptions
    DO J = 1, NINVSCC
        K = SCCIDX( J )
        WRITE( LDEV,94675 ) INVSCC( J ), SCCSECTN( K ),&
        &                    TRIM( SCCDESC( J ) )
    END DO
    WRITE( LDEV, 94010 ) ' '

!.........  Deallocate borrowed MODLISTS arrays
    DEALLOCATE( INVSCC, SCCDESC )

!.........  Deallocate local memory
    DEALLOCATE( IDXA2P, NUMSCC, UNSRTA2P, AR2PTSCC )
    DEALLOCATE( SCCIDX, INVSCCA, SCCSECTN )

!.........  Deallocate cross-reference sorting arrays
    DEALLOCATE( INDXTA, CFIPTA, CSCCTA, CSRCTA, IARPTA )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT(  A, 1X, F10.7, 1X, A )

94675 FORMAT( 10X, A20, 4X, I2.2, 5X, A )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram counts entries in or reads
!               the contents of the area-to-point factors file
    SUBROUTINE READ_ARTOPT_FILE( STATUS, MXSCC )

!...............   EXTERNAL FUNCTIONS and their descriptions:
        LOGICAL, EXTERNAL :: CHKINT
        LOGICAL, EXTERNAL :: CHKREAL
        LOGICAL, EXTERNAL :: BLKORCMT
        INTEGER, EXTERNAL :: GETFLINE
        INTEGER, EXTERNAL :: GETNLIST

!.............  Subprogram arguments
        CHARACTER(*), INTENT(IN   ) :: STATUS  ! call status: COUNT|STORE
        INTEGER     , INTENT(INOUT) :: MXSCC   ! count the max SCC per section

!.............  Local variables
        INTEGER          I, L, L1, L2, K, N, NS   ! counters and indices

        INTEGER          CNT        ! table record counter
        INTEGER          CNY        ! tmp county code
        INTEGER          COU        ! tmp country code
        INTEGER          IOS        ! i/o status
        INTEGER          IREC       ! record number
        INTEGER       :: LINSCC = 0 ! number of SCCs on current header line
        INTEGER, SAVE :: NLINES     ! no. lines in input file
        INTEGER       :: NSCC = 0   ! SCC counter
        INTEGER          NTBL       ! no. tables
        INTEGER          STA        ! tmp state code

        REAL             ALLOC      ! tmp allocation factor
        REAL             LAT        ! tmp latitude
        REAL             LON        ! tmp longitude

        CHARACTER(FIPLEN3) CFIP     ! tmp co/st/cy FIPS code
        CHARACTER(SCCLEN3) SEGSCC   ! segment SCC

        LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time routine called
        LOGICAL       :: PREVPKT  = .FALSE. ! true: previous line was a packet
        LOGICAL       :: THISPKT  = .FALSE. ! true: this line is a packet

!......................................................................

!.............  If the first time the internal subprogram is called
        IF( FIRSTIME ) THEN

!................  Determine the number of lines in the input file
            NLINES = GETFLINE( FDEV, "ARTOPNT file" )

            FIRSTIME = .FALSE.

        END IF

!............  Loop through lines of input file and
        IREC    = 0
        NTBL    = 0
        PREVPKT = .FALSE.
        DO I = 1, NLINES

            READ( FDEV, 93000, END=2002, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                WRITE( MESG, 94010 )&
                &   'ERROR: System I/O error', IOS, 'reading ' //&
                &   'area-to-point factors file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!.................  Skip comment and blank lines
            IF( BLKORCMT ( LINE ) ) CYCLE

!................  If line is a header line, add number of SCCs for
!                  the current table
            THISPKT = .FALSE.
            L = INDEX( LINE, '/LOCATIONS/' )
            IF( L .GT. 0 ) THEN

!...................  If previous line was not another packet
                IF( .NOT. PREVPKT ) THEN

!.......................  Update the number of tables
                    NTBL = NTBL + 1

!.......................  Initialize the SCC count
                    L1 = 12               ! based on LOCATIONS packet
                    L2 = LEN_TRIM( LINE )
                    LINSCC = GETNLIST( L2 - L1 + 1, LINE( L1:L2 ) )
                    NSCC = LINSCC
                    NS   = 0

!.......................  Initialise count of entries per table
                    CNT = 0

!...................  Otherwise add to the SCC count
                ELSE
                    L1 = 12               ! based on LOCATIONS packet
                    L2 = LEN_TRIM( LINE )
                    LINSCC = GETNLIST( L2 - L1 + 1, LINE( L1:L2 ) )
                    NSCC = NSCC + LINSCC

                END IF

!....................  Store maximum SCC value
                MXSCC = MAX( MXSCC, NSCC )

!....................  Give error if no SCCs are provided in the packet
                IF( LINSCC .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: No SCCs provided '//&
                    &       'in /LOCATIONS/ packet at line', IREC
                    CALL M3MSG2( MESG )
                END IF

!....................  Set packet status
                THISPKT = .TRUE.

            END IF

!...............  Make sure that header has appeared at least once
            IF( NTBL .EQ. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: No header line found '//&
                &       'before line', IREC
                CALL M3MSG2( MESG )
                CYCLE
            END IF

!...............  Increment table counter and keep track of maximum
            IF( .NOT. THISPKT ) THEN
                CNT = CNT + 1
                MXROWA2P = MAX( MXROWA2P, CNT )
            END IF

!................  End read if just counting
            IF( STATUS .EQ. 'STORE' ) THEN

!....................  If header line, store SCCs and count
                IF( THISPKT ) THEN
                    NUMSCC( NTBL ) = NSCC

                    CALL PARSLINE( LINE( L1:L2 ), LINSCC, SEGMENT )

                    DO N = 1, LINSCC
                        NS = NS + 1

                        L = LEN_TRIM( SEGMENT( N ) )
                        IF( L .GT. SCCLEN3 ) THEN
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 ) 'ERROR: SCC "'//&
                            &  SEGMENT( 1:L ) // '" at line', IREC,&
                            &  'has', L, 'characters, which '//&
                            &  'exceeds the maximum of', SCCLEN3
                            CALL M3MSG2( MESG )
                            CYCLE
                        END IF

                        SEGSCC = TRIM( SEGMENT( N ) )
                        CALL PADZERO( SEGSCC )
                        AR2PTSCC( NS,NTBL ) = SEGSCC

                    END DO

!....................  Otherwise, parse and store table
                ELSE

                    DO N = 1, NFIELDS
                        SEGMENT(N)= ADJUSTL( LINE(FBEG(N):FEND(N)) )
                    END DO

!........................  Check that FIPS code is an integer
                    IF( .NOT. USEEXPGEO() .AND.&
                    &    .NOT. CHKINT( SEGMENT( 1 ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Country/state'//&
                        &       '/county code is not an integer '//&
                        &       'at line', IREC
                        CALL M3MSG2( MESG )
                    END IF

!........................  Check that reals are reals
                    IF( .NOT. CHKREAL( SEGMENT( 8 ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Longitude is '//&
                        &  'not a floating point value at line', IREC
                        CALL M3MSG2( MESG )
                    END IF
                    IF( .NOT. CHKREAL( SEGMENT( 9 ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Latitude is '//&
                        &  'not a floating point value at line', IREC
                        CALL M3MSG2( MESG )
                    END IF
                    IF( .NOT. CHKREAL( SEGMENT( 10 ) ) ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'ERROR: Allocation '//&
                        &  'factor is not a floating point value '//&
                        &  'at line', IREC
                        CALL M3MSG2( MESG )
                    END IF

                    IF( EFLAG ) CYCLE

!........................  Store total count of packet for this table
                    NAR2PT( NTBL ) = CNT

!........................  Store contents of this row
                    CFIP  = SEGMENT( 1 )
                    CALL PADZERO( CFIP )
                    LON   = STR2REAL( SEGMENT( 8 ) )
                    IF( WFLAG .AND. LON > 0 ) LON = -LON  ! convert to western hemisphere
                    LAT   = STR2REAL( SEGMENT( 9 ) )
                    ALLOC = STR2REAL( SEGMENT( 10 ) )

                    IDXA2P        ( CNT,NTBL ) = CNT
                    UNSRTA2P( CNT,NTBL )%FIP   = CFIP
                    UNSRTA2P( CNT,NTBL )%LAT   = LAT
                    UNSRTA2P( CNT,NTBL )%LON   = LON
                    UNSRTA2P( CNT,NTBL )%ALLOC = ALLOC
                    UNSRTA2P( CNT,NTBL )%NAME  = TRIM( SEGMENT(5) )
                END IF
            END IF

!...............  Store packet status for next iteration
            PREVPKT = THISPKT

        END DO   ! end of loop through lines

        NTABLA2P = NTBL

!.............  Abort if errors found
        IF( EFLAG ) THEN
            MESG = 'Problem reading area-to-point factors file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Rewind file for next call
        REWIND( FDEV )

        RETURN

!.............  Abort with errors...
2002    WRITE( MESG,94010 ) 'ERROR: Unexpected end of area-to-'//&
        &       'point factors file at line', IREC
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!************   SUBPROGRAM FORMAT  STATEMENTS   *************************

!...........   Formatted file I/O formats............ 93xxx

93000   FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE READ_ARTOPT_FILE

!----------------------------------------------------------------------

END SUBROUTINE RDAR2PT
