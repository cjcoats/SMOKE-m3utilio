
SUBROUTINE RDEMSPD( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES,    &
                    GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME,    &
                    EDATE, ETIME, EASTAT, SPSTAT )

    !***************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the day- or hour-specific emissions in
    !      EMS-95 format. It appends the records to the global storage from the
    !      MODDAYHR.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutine
    !
    !  REVISION  HISTORY:
    !       Created 12/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
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

    !.......  MODULES for public variables
    !.......  This module is the inventory arrays
    USE MODSOURC, ONLY: CIFIP, CSOURC, INTGRFLAG, CINTGR

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP, NINVTBL, ITFACA, ITNAMA, &
                        ITKEEPA, SORTCAS, SCASIDX, NUNIQCAS,        &
                        UCASNPOL, UNIQCAS, UCASIDX, UCASNKEP,       &
                        INVDVTS, MXIDAT, INVDNAM

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS

    !.......  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, IDXSRC, SPDIDA,      &
                        CODEA, EMISVA, DYTOTA, CIDXA

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV               ! file unit no.
    INTEGER, INTENT (IN) :: TZONE              ! output time zone
    INTEGER, INTENT (IN) :: TSTEP              ! time step HHMMSS
    INTEGER, INTENT (IN) :: MXPDSRC            ! max. day- or hr-specific source
    LOGICAL, INTENT (IN) :: GETSIZES           ! true: get no. time steps  pols
    LOGICAL, INTENT (IN) :: GETCOUNT           ! true: get max no. srcs per time
    LOGICAL, INTENT (IN) :: FIRSTCALL          ! true: first call of a loop
    LOGICAL, INTENT (IN) :: DAYFLAG            ! true: day-, false: hour-spec
    INTEGER,INTENT(INOUT):: SDATE              ! Julian starting date in TZONE
    INTEGER,INTENT(INOUT):: STIME              ! start time of data in TZONE
    INTEGER, INTENT(OUT) :: EDATE              ! Julian ending date in TZONE
    INTEGER, INTENT(OUT) :: ETIME              ! ending time of data in TZONE
    INTEGER, INTENT(OUT) :: EASTAT( NIPPA )     ! true: pol/act appears in data
    INTEGER, INTENT(OUT) :: SPSTAT( MXSPDAT )     ! true: special in data

    !.......   Local list of bad sources to prevent duplicate writing of error
    !              messages
    CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

    !.......   Local list of FIPS start/end positions to facilitate
    !              faster lookups
    INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
    INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )

    !.......   Local list of arrays for warning handling
    LOGICAL, ALLOCATABLE, SAVE :: WARNKEEP( : )     ! true: write warning for Keep = N
    LOGICAL, ALLOCATABLE, SAVE :: WARNMULT( : )     ! true: write warning for Multiple pollutants from a single pollutant in Inventory Table

    !.......   Temporary read arrays
    REAL            TDAT( 24 )           ! temporary data values

    !.......   Other local variables
    INTEGER          H, HS, I, J, L, LL, L1, L2, NV, S, T        ! counters and indices
    INTEGER          ES, NS, SS        ! end src, tmp no. src, start sourc

    INTEGER          CIDX                 ! tmp data index
    INTEGER          COD                  ! data index
    INTEGER          DAY                  ! tmp day of month
    INTEGER          FIP                  ! tmp co/st/cy code
    INTEGER, SAVE :: ICC = 0              ! tmp country code from header
    INTEGER          IOS                  ! i/o status
    INTEGER          IREC                 ! record counter
    INTEGER          JDATE                ! tmp Julian date
    INTEGER          JTIME                ! tmp HHMMSS time
    INTEGER          LFIP                 ! previous st/co FIPS code
    INTEGER, SAVE :: LOOPNO = 0           ! no. of loops
    INTEGER, SAVE :: MAXPTR               ! maximum time step reference pointer
    INTEGER, SAVE :: MINPTR               ! minimum time step reference pointer
    INTEGER          MONTH                ! tmp month number
    INTEGER, SAVE :: MXWARN       	      ! max no. warnings
    INTEGER, SAVE :: NBADSRC = 0          ! no. bad sources
    INTEGER, SAVE :: NFIELD = 0           ! number of data fields
    INTEGER, SAVE :: NFM1   = 0           ! number of data fields minus 1
    INTEGER       :: NPOA   = 0           ! unused header number of pol/act
    INTEGER, SAVE :: NSTEPS = 0           ! number of time steps
    INTEGER, SAVE :: NWARN( 5 )           ! warnings counter
    INTEGER          PTR                  ! tmp time step pointer
    INTEGER       :: RDATE = 1980001      ! reference date: Jan 1, 1980
    INTEGER       :: RTIME = 0            ! reference time
    INTEGER, SAVE :: S1 = 0               ! saved 1st position of extra field
    INTEGER, SAVE :: S2 = 0               ! saved 2nd position of extra field
    INTEGER, SAVE :: SDATESAV = 0         ! saved start date
    INTEGER, SAVE :: STIMESAV = 0         ! saved start time
    INTEGER, SAVE :: TDIVIDE  = 1         ! time step divisor
    INTEGER          WD                   ! tmp field width
    INTEGER          YEAR                 ! 4-digit year
    INTEGER       :: YR4 = 0              ! unused header year
    INTEGER          DZONE                ! time shift (ZONE-TZONE)
    INTEGER          ZONE                 ! source time zones

    REAL             CONVFAC              ! tmp conversion factor from Inventory Table
    REAL             TOTAL                ! tmp daily total of hourly file

    LOGICAL, SAVE :: DFLAG = .FALSE.      ! true: dates set by data
    LOGICAL       :: EFLAG = .FALSE.      ! TRUE iff ERROR
    LOGICAL       :: WARNOUT = .FALSE.    ! true: then output warnings
    LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called
    LOGICAL, SAVE :: LFLAG  = .FALSE.     ! true: output daily/hourly inv in local time
    LOGICAL, SAVE :: SFLAG  = .FALSE.     ! true: use daily total from hourly
    LOGICAL, SAVE :: TFLAG  = .FALSE.     ! true: use SCCs for matching with inv
    LOGICAL, SAVE :: NFLAG  = .FALSE.     ! true: Special pollutant name field used on that line
    LOGICAL, SAVE :: WIDE_FORMAT = .FALSE.     ! true: hourly data values given with 12 columns instead of 7

    CHARACTER(100) :: BUFFER = ' '        ! src description buffer
    CHARACTER(512) :: LINE   = ' '        ! line buffer
    CHARACTER(512) :: MESG   = ' '        ! message buffer

    CHARACTER(FIPLEN3) CFIP          ! tmp co/st/cy code
    CHARACTER(CASLEN3) CDAT          ! tmp Inventory data (input) name
    CHARACTER(NAMLEN3) CNAM          ! tmp SMOKE name
    CHARACTER(NAMLEN3) PNAM          ! tmp SMOKE name
    CHARACTER(CHRLEN3) CHAR4         ! tmp characteristic 4
    CHARACTER(PLTLEN3) FCID          ! tmp facility ID
    CHARACTER(CHRLEN3) SKID          ! tmp stack ID
    CHARACTER(CHRLEN3) DVID          ! tmp device ID
    CHARACTER(CHRLEN3) PRID          ! tmp process ID
    CHARACTER(SCCLEN3) TSCC          ! tmp source category code
    CHARACTER(ALLLEN3) CSRC          ! tmp source string

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDEMSPD'     !  program name

    !***********************************************************************
    !   begin body of program RDEMSPD

    !.......  First time routine called
    IF( FIRSTIME ) THEN

        !.......  Get environment variable using an hourly file as a daily file
        !.......  NOTE - the hourly file will have been assigned as a daily
        !         file when it was opened.
        MESG = 'Use daily totals only from hourly data file'
        SFLAG = ENVYN( 'HOURLY_TO_DAILY', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "HOURLY_TO_DAILY"', 2 )
        END IF

        !.......  No time zone shift for AERMOD support
        MESG = 'Outputs local time daily and/or hourly inventories (No time shift)'
        LFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTPUT_LOCAL_TIME"', 2 )
        END IF

        !.......  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET , ' ', 100, I )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
        END IF

        !.......  Determine if wideformat is being used for for hourly file
        IF( .NOT. DAYFLAG ) THEN
            MESG = 'Use 12-column format for hourly values'
            WIDE_FORMAT = ENVYN( 'HOURLY_WIDE_FMT', MESG, .FALSE., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "HOURLY_WIDE_FMT"', 2 )
            END IF
        END IF

        !.......  Give note if file is being read as a daily file
        IF( DAYFLAG .AND. SFLAG ) THEN
            MESG = 'NOTE: Daily data only being used from an hourly emissions file'
            CALL M3MSG2( MESG )

        !.......  Otherwise, ignore setting because it is an hourly file
        ELSE IF( SFLAG ) THEN
            SFLAG = .FALSE.
            MESG = 'NOTE: Ignoring HOURLY_TO_DAILY setting for reading hourly emissions data'
            CALL M3MSG2( MESG )
        END IF

        !.......  Allocate memory for bad source storage
        ALLOCATE( BADSRC( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BADSRC', PROGNAME )

        !.......  Create unique list of FIPS codes and other things
        CALL GENUSLST

        !.......  Build helper arrays for making searching faster
        ALLOCATE( STARTSRC( NINVIFIP ),    &
                    ENDSRC( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ENDSRC', PROGNAME )
        STARTSRC = 0
        ENDSRC   = 0
        S = 0
        DO I = 1, NINVIFIP
            DO
                S = S + 1
                IF ( S .GT. NSRC ) EXIT
                IF( CIFIP( S ) .EQ. INVCFIP( I ) ) THEN
                    IF( STARTSRC( I ) .EQ. 0 ) STARTSRC( I ) = S
                    ENDSRC( I ) = S
                ELSE
                    S = S - 1
                    EXIT
                END IF
            END DO
        END DO

        !.......  Initialize warnings counter
        NWARN = 0          ! array

        FIRSTIME = .FALSE.

    END IF

    !.......  For the first call in a loop of files, initialize variables
    IF( FIRSTCALL ) THEN
        MINPTR  = 99999999
        MAXPTR  = 0

        !.......  Set time step divisor
        TDIVIDE = 3600 * TSTEP / 10000

        !.......  Set the number of fields, depending on day- or hour-specific
        IF( DAYFLAG ) THEN
            NFIELD = 1
            S1 = 103           ! special pollutant field start position
            S2 = 118           ! special pollutant field end position
        ELSE IF( WIDE_FORMAT ) THEN          ! wide format hourly
            NFIELD  = 24
            S1 = 381
            S2 = 396
        ELSE                      ! standard format hourly
            NFIELD  = 24
            S1 = 261
            S2 = 276
        END IF
        NFM1   = NFIELD - 1

        !.......  If dates have been set by the data, set the number of steps
        !               steps
        IF( DFLAG ) THEN
            NSTEPS = 1+ SECSDIFF( SDATE,STIME,EDATE,ETIME )/ TDIVIDE
            SDATESAV = SDATE
            STIMESAV = STIME
        END IF

        !.......  Set switch for printing errors only the first loop through all
        !         of the input files.  The second time through is indicated
        !         for the second time that FIRSTCALL is true.
        !.......  Reset loop counter if call is to get dimensions only (because
        !         this means it is the first call or daily or hourly)
        IF( GETSIZES ) LOOPNO = 0
        LOOPNO = LOOPNO + 1
        WARNOUT = ( LOOPNO .EQ. 1 )

        !.......  Deallocate warning arrays
        IF( ALLOCATED( WARNKEEP ) ) DEALLOCATE( WARNKEEP, WARNMULT )
        ALLOCATE( WARNKEEP( NUNIQCAS ),    &
                  WARNMULT( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WARNMULT', PROGNAME )
        WARNKEEP = .TRUE.
        WARNMULT = .TRUE.

    END IF

    !.......  Loop through file and read it. In the first section, determine
    !         the minimum and maximum date. Use a reference date to do this. In
    !         the second section, determine the number of records per time
    !         step. In the third section, read and store the data.  When storing
    !         data, time step index is computed from the start date/time instead
    !         of the reference date/time so that the indexing will work properly.
    LFIP = 0
    IREC = 0
    TDAT = 0       !  array
    DO             !  Head of period-specific file read loop

        !.......  Read first line of file
        READ( FDEV, 93000, END=299 ) LINE
        IREC = IREC + 1

        L = LEN_TRIM( LINE )

        !.......  Skip blank lines
        IF( L .EQ. 0 ) CYCLE

        !.......  Determine whether special pollutant field is used
        NFLAG = .FALSE.
        IF ( L .GT. S1 ) NFLAG = .TRUE.

        !.......  Scan for header lines and check to ensure all are set
        !               properly
        CALL GETHDR( 1, .FALSE., .FALSE., .FALSE.,    &
                     LINE, ICC, YR4, NPOA, IOS )

        !.......  Interpret error status
        IF( IOS .EQ. 4 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: DATA header entry should not be used'//    &
                   'for EMS-95 day- or hour-specific files.'
            CALL M3MSG2( MESG )

        ELSE IF( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.

        END IF

        !.......  If a header line was encountered, go to next line
        IF( IOS .GE. 0 ) CYCLE

        !.......  Determine if file is day- or hour-specific by the length of the
        !         lines. Make sure day- and hour-specific data are not in the
        !         same file.
        !.......  If the file is hourly but the only the daily is to be read, then
        !         behave as if it is a daily file.
        IF( DAYFLAG .AND.    &
               ( L .GT. S2  .AND. .NOT. SFLAG ) .OR.      &        ! line longer that allowed
               ( L .LE. 101 .AND.       SFLAG )      ) THEN        ! line shorter than allowed
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: bad format or daily ' //    &
                   'data found in day-specific file at line', IREC
            CALL M3MESG( MESG )
            CYCLE

        !.......  If daily info being read from hourly file, make sure last
        !               column is present
        ELSE IF( DAYFLAG .AND. SFLAG .AND. L .LT. 240 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: bad format in hour-' //    &
                   'specific file being used for daily '//    &
                   CRLF() // BLANK10 // 'data at line', IREC
            CALL M3MESG( MESG )
            CYCLE

        !.......  Check to make sure that the last hourly field is present.
        !          The daily total field does not have to be there.
        ELSE IF( .NOT. DAYFLAG .AND. L .LT. 234 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: bad format or hourly ' //    &
                   'data found in hour-specific file at line', IREC
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.......  Set Julian day from MMDDYY8 SAS format
        MONTH = STR2INT( LINE( 62:63 ) )
        DAY   = STR2INT( LINE( 65:66 ) )
        YEAR  = YEAR4( STR2INT( LINE( 68:69 ) ) )

        JDATE = 1000 * YEAR + JULIAN( YEAR, MONTH, DAY )
        JTIME = 0

        !.......  Search for time zone name from file in master list
        CALL UPCASE( LINE( 70:72 ) )
        I = INDEX1( LINE( 70:72 ), MXTZONE, TZONNAM )

        !.......  If time zone name is not found, thenoutput error
        IF( I .LE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )    &
                  'ERROR: Unrecognized time zone "'// LINE(70:72)//    &
                  '" at line', IREC, 'in file'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Set time zone number
        IF( .NOT. LFLAG ) THEN

            !.......  If daily emissions are not in the output time zone, 
            !         print warning
            IF( WARNOUT .AND. DAYFLAG .AND. ZONE .NE. TZONE .AND.    &
                NWARN( 1 ) .LE. MXWARN ) THEN
                WRITE( MESG,94010 )    &
                      'WARNING: Time zone ', ZONE, 'in day-specific ' //    &
                      'file at line', IREC, CRLF() // BLANK10 //    &
                      'does not match output time zone', TZONE
                CALL M3MESG( MESG )
                NWARN( 1 ) = NWARN( 1 ) + 1
            END IF

        !.......  Reset time shift to 0 to correctly compute local time zone
        ELSE
            DZONE = 0

        END IF

        !.......  Convert date and time to output time zone.
        CALL NEXTIME( JDATE, JTIME, DZONE * 10000 )

        !.......  Determine time step pointer based on reference time
        PTR = SECSDIFF( RDATE, RTIME, JDATE, JTIME ) / TDIVIDE + 1

        !.......  Store minimum time step number as compared to reference
        IF( PTR .LT. MINPTR ) MINPTR = PTR

        !.......  Store maximum time step number as compared to reference
        IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

        !.......  Set key for searching sources
        FIP  = 100000 * ICC + 1000 * STR2INT( LINE( 1:2 ) ) + STR2INT( LINE( 3:5 ) )
        WRITE( CFIP,94020 ) FIP

        FCID = ADJUSTL( LINE( 6:20 ) )
        SKID = ADJUSTL( LINE( 21:32 ) )           ! point ID in IDA
        DVID = ADJUSTL( LINE( 33:44 ) )           ! stack ID in IDA
        PRID = ADJUSTL( LINE( 45:56 ) )           ! segment in IDA

        TSCC = ' '
        !.......  If FIPS code is not the same as last time, then
        !         look it up and get indidies
        IF( FIP .NE. LFIP ) THEN
            J = FINDC( CFIP, NINVIFIP, INVCFIP )
            IF( J .LE. 0 ) THEN
                WRITE( MESG,94010 ) 'INTERNAL ERROR: Could not '//    &
                       'find FIPS code', FIP, 'in internal list.'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            SS = STARTSRC( J )
            ES = ENDSRC( J )
            NS = ES - SS + 1
            LFIP = FIP

        END IF

        !.......  If SCCs are needed for matching...
        IF ( TFLAG ) THEN
            IF ( DAYFLAG ) THEN
                TSCC = ADJUSTL( LINE( 92:101) )
            ELSE IF ( WIDE_FORMAT ) THEN
                TSCC = ADJUSTL( LINE( 373:380 ) )
            ELSE         ! hourly standard format
                TSCC = ADJUSTL( LINE( 250:259 ) )
            END IF
            IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )
            CHAR4 = TSCC

            !.......  Build source characteristics field for searching inventory
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,    &
                          TSCC, CHRBLNK3, POLBLNK3, CSRC )

            !.......  Search for this record in sources
            J = FINDC( CSRC, NS, CSOURC( SS ) )

        !.......  If SCCs are not being used for matching (at least not yet)...
        ELSE

            !.......  Build source characteristics field for searching inventory
            CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,    &
                          TSCC, CHRBLNK3, POLBLNK3, CSRC )

            !.......  Search for this record in sources
            J = FINDC( CSRC, NS, CSOURC( SS ) )

            !.......  If source is not found for day-specific processing, see
            !         if reading the SCC in helps (needed for IDA format)
            IF( J .LE. 0 ) THEN

                IF ( DAYFLAG ) THEN
                    TSCC = ADJUSTL( LINE( 92:101) )
                ELSE IF ( WIDE_FORMAT ) THEN
                    TSCC = ADJUSTL( LINE( 373:380 ) )
                ELSE
                    TSCC = ADJUSTL( LINE( 250:259 ) )
                END IF
                IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )
                CHAR4 = TSCC

                !.......  Build source characteristics field for searching inventory
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,    &
                              TSCC, CHRBLNK3, POLBLNK3, CSRC )

                !.......  Search for this record in sources
                J = FINDC( CSRC, NS, CSOURC( SS ) )
                IF ( J .GT. 0 ) TFLAG = .TRUE.

            END IF

        END IF

        !.......  Store source in list of bad sources
        !.......  Print warning about sources not found in the inventory
        IF( J .LE. 0 ) THEN

            !.......  Search for source in list of bad sources
            J = INDEX1( CSRC, NBADSRC, BADSRC )

            !.......  If source is not found, give a message.  Don't need the
            !                   WARNOUT controller because this section only gets
            !                   invoked once.
            IF( J .LE. 0 ) THEN

                NBADSRC = NBADSRC + 1
                BADSRC( NBADSRC ) = CSRC

                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                IF( NWARN( 3 ) .LE. MXWARN ) THEN
                    MESG = 'WARNING: Period-specific record does '//    &
                           'not match inventory sources: '//    &
                           CRLF() // BLANK10 // BUFFER( 1:L2 )
                    CALL M3MESG( MESG )
                    NWARN( 3 ) = NWARN( 3 ) + 1
                END IF

            END IF

            CYCLE                       !  to head of read loop

        !.......  Otherwise, update master list of sources in the inventory
        ELSE
            S = SS - 1 + J                 ! calculate source number

        END IF

        !.......  Check pollutant code and set index I
        IF( NFLAG ) THEN
            CDAT = LINE( S1:S2 )

        !.......  Otherwise, use original pollutant field
        ELSE
            CDAT = LINE( 57:61 )
        END IF

        !.......  Left justify and convert pollutant name to upper case
        CDAT = ADJUSTL( CDAT )
        CALL UPCASE( CDAT )

        !.......  Look up pollutant name in unique sorted array of
        !               Inventory pollutant names
        CIDX  = FINDC( CDAT, NUNIQCAS, UNIQCAS )

        !.......  Check to see if data name is in inventory list
        COD  = INDEX1( CDAT, NIPPA, EANAM )

        !.......  If pollutant name is not in Inventory Table list
        IF ( CIDX .LE. 0 ) THEN

            !.......  Check to see if data name is in list of special names
            CIDX= INDEX1( CDAT, MXSPDAT, SPDATNAM )

            !.......  Store status of special data and flag code with
            !         special integer so can ID these records later.
            IF( CIDX .GT. 0 ) THEN
                SPSTAT( CIDX ) = CIDX
                COD = CODFLAG3 + CIDX

            ELSE IF ( CIDX .LE. 0 ) THEN
                IF( WARNOUT .AND. NWARN( 2 ) .LE. MXWARN ) THEN
                    WRITE( MESG,94010 )    &
                     'WARNING: Skipping pollutant "'// TRIM(CDAT)//    &
                     '" at line', IREC, '- not in Inventory Table'
                    CALL M3MESG( MESG )
                    NWARN( 2 ) = NWARN( 2 ) + 1
                END IF
                CYCLE              !  to head of loop

            END IF

        !.......  Otherwise, pollutant is in list of Inventory Data Names
        ELSE

            !.......  Write warning if pollutant is not kept.  Write only
            !                   one time.
            IF( UCASNKEP(CIDX) .LE. 0 .AND. WARNKEEP(CIDX) ) THEN
                WARNKEEP( CIDX ) = .FALSE.
                IF( GETSIZES ) THEN
                    WRITE( MESG,94010 )    &
                      'WARNING: Skipping all lines for pollutant "'//    &
                      TRIM( CDAT )// '" because pollutant is not '//    &
                      'kept by Inventory Table.'
                    CALL M3MESG( MESG )
                END IF
                CYCLE
            END IF

            !......  Get Inventory Data SMOKE name from Inventory Table arrays/indices
            CNAM = ITNAMA( SCASIDX( UCASIDX( CIDX ) ) )

            !......  Look up SMOKE name in list of annual EI pollutants
            COD  = INDEX1( CNAM, NIPPA, EANAM )

            !......  Check to ensure that it handles NOI and NONHAP pollutants
            !        while combining VOC + HAPs
            IF( INTGRFLAG ) THEN

            !......  Preventing processing precomputed NONHAP[VOC|TOG]
                IF( INDEX( CNAM,'NONHAP' ) > 0 ) THEN
                    MESG = 'ERROR: Can NOT process precomputed '// TRIM(CNAM)//    &
                        ' when SMK_PROCESS_HAPS was set to process anuual inventory'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG , 2 )
                END IF

                NV = INDEX1( CNAM, MXIDAT, INVDNAM )

                IF( CINTGR( S ) == 'N' .AND. INVDVTS( NV ) /= 'N' ) THEN
                    PNAM = TRIM( CNAM ) // '_NOI'
                    COD = INDEX1( PNAM, NIPPA, EANAM )

                ELSE IF( CINTGR( S ) == 'Y' ) THEN
                    L = INDEX( CNAM, ETJOIN )
                    LL= LEN_TRIM( CNAM )
                    PNAM = CNAM
                    IF( L > 0 ) PNAM = CNAM( L+2:LL )
                    IF( PNAM == 'VOC' .OR. PNAM == 'TOG' ) THEN
                        IF( L > 0 ) THEN
                            PNAM = CNAM(1:L+1) // 'NONHAP' //    &
                                   CNAM(L+2:LL)
                        ELSE
                            PNAM = 'NONHAP' // TRIM( CNAM )
                        END IF
                        COD = INDEX1( PNAM, NIPPA, EANAM )
                    END IF

                END IF

            END IF

            !......  Check to ensure that the SMOKE intermediate name
            !        set by the Inventory Table is actually in the annual
            !        inventory.  If not, write warning message and cycle.
            IF( COD .LE. 0 ) THEN
                IF( WARNOUT .AND. NWARN( 5 ) .LE. MXWARN ) THEN
                    WRITE( MESG,94010 )    &
                      'WARNING: Skipping pollutant "'// TRIM(CNAM)//    &
                      '" at line', IREC, '- not in annual inventory.'
                    CALL M3MESG( MESG )
                    NWARN( 5 ) = NWARN( 5 ) + 1
                END IF
                CYCLE

            !......  If it's found, then record that this pollutant was found
            ELSE
                EASTAT( COD ) = CIDX
            END IF

        END IF          ! if cidx le 0 or not

        !.......  If only getting dates and pollutant information, go
        !         to next loop iteration
        IF( GETSIZES ) CYCLE

        !.......  Determine time step pointer based on actual start time
        PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/ TDIVIDE + 1

        !.......  Skip record if it is out of range of output file
        !.......  NOTE - this is only useful if reading only part of data
        IF( PTR .LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE

        !.......  Count estimated record count per time step
        DO T = PTR, MIN( PTR + 23, NSTEPS )
            MXPDPT( T ) = MXPDPT( T ) + 1
        END DO

        !.......  If only counting records per time step, go to next loop iteration
        IF( GETCOUNT ) CYCLE

        !.......  Set column locations for reading data file...
        !.......  Day-specific from an hourly file - read totals
        IF( DAYFLAG .AND. SFLAG ) THEN
            L1 = 233
            L2 = 240
            WD = 8
        !.......  Day-specific from a day-specific file
        ELSE IF( DAYFLAG ) THEN
            L1 = 55
            L2 = 72
            WD = 18
        !.......  Hourly file
        ELSE
            L1 = 66
            L2 = 72
            IF ( WIDE_FORMAT ) L2 = 77
            WD = 7
        END IF

        !.......  Check and set emissions values
        DO J = 1, NFIELD

            L1 = L1 + WD
            L2 = L2 + WD

            IF ( WIDE_FORMAT ) WD = 12

            TDAT( J )  = STR2REAL( LINE( L1:L2 ) )
            IF ( TDAT( J ) .LT. 0.0 )  THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Bad line', IREC,    &
                       ': data value "' // LINE( L1:L2 ) // '"'
                CALL M3MESG( MESG )
                CYCLE          ! to head of read loop
            END IF

        END DO

        !.......  If daily data, set all TDATs with daily value
        IF( DAYFLAG ) TDAT( 2:24 ) = TDAT( 1 )          ! array

        !.......  If available, set total value from hourly file
        TOTAL = 0.
        IF( SFLAG .OR. .NOT. DAYFLAG ) THEN
            L1 = L2 + 1
            IF( WIDE_FORMAT ) THEN
                L2 = L2 + WD
            ELSE
                L2 = L2 + 8         !  note: original hourly format has 8 columns for daily total (but 7 for hourly values)
            END IF

            IF( LINE( L1:L2 ) .NE. ' ' ) THEN
                TOTAL = STR2REAL( LINE( L1:L2 ) )
                IF( TOTAL .LT. 0.0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Bad line', IREC,    &
                      ': total value "' // LINE(L1:L2) // '"'
                    CALL M3MESG( MESG )
                    CYCLE          ! to head of read loop
                END IF
            END IF
        END IF

        !.......  Store source ID
        LPDSRC( S ) = .TRUE.

        !.......  Set conversion factor from Inventory Table. Default is
        !         1., which is also what is used in all but a handful of
        !         special toxics cases.
        CONVFAC = ITFACA( SCASIDX( UCASIDX( CIDX ) ) )

        !.......  Record needed data for this source and time step
        H = 0
        DO T = PTR, MIN( PTR + 23, NSTEPS )

            H = H + 1
            NPDPT( T ) = NPDPT( T ) + 1

            HS = NPDPT( T )

            IF( HS .LE. MXPDSRC ) THEN

                IDXSRC( HS,T ) = HS
                SPDIDA( HS,T ) = S
                CIDXA ( HS,T ) = CIDX
                CODEA ( HS,T ) = COD
                EMISVA( HS,T ) = CONVFAC * TDAT( H )          ! Store data in emissions
                DYTOTA( HS,T ) = CONVFAC * TOTAL

            END IF

        END DO

    END DO

299 CONTINUE       ! Exit from read loop

    !.......  Abort if error found while reading file
    IF( EFLAG ) THEN
        MESG = 'Problem processing day- or hour-specific data'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Update output starting date/time and ending date/time
    DFLAG = .TRUE.
    SDATE = RDATE
    STIME = RTIME
    DO I = 1, MINPTR - 1
        CALL NEXTIME( SDATE, STIME, TSTEP )
    END DO

    EDATE = RDATE
    ETIME = RTIME
    DO I = 1, MAXPTR - 1
        CALL NEXTIME( EDATE, ETIME, TSTEP )
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( I6.6 )

END SUBROUTINE RDEMSPD
