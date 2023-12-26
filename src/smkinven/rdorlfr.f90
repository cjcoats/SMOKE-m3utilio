
SUBROUTINE RDORLFR( FDEV, TZONE, TSTEP, MXPDSRC, GETSIZES,      &
                    GETCOUNT, FIRSTCALL, DAYFLAG, SDATE, STIME, &
                    EDATE, ETIME, EASTAT, SPSTAT )

    !***************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the day- emissions in ORL FIREEMIS format.
    !      It appends the records to the global storage from the MODDAYHR.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutine
    !
    !  REVISION  HISTORY:
    !       Created 03/06 by B.H. Baek (based on RDEMSPD.F)
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
    USE MODSOURC, ONLY: CIFIP, CSOURC, HEATCONTENT, INTGRFLAG, CINTGR

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP, UCASNKEP, NUNIQCAS,      &
                        UNIQCAS, NINVTBL, ITNAMA, ITCASA, FIREFLAG, &
                        UCASIDX, SCASIDX, INVDVTS, MXIDAT, INVDNAM

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPPA, NSRC, EANAM, NCHARS, NMAP, MAPNAM, MAPFIL

    !.......  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: MXPDPT, LPDSRC, NPDPT, NPDPTP, IDXSRC,      &
                        SPDIDA, CODEA, CIDXA, EMISVA, DYTOTA

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......  SUBROUTINE ARGUMENTS
    INTEGER, INTENT( IN ):: FDEV               ! file unit no.
    INTEGER, INTENT( IN ):: TZONE              ! output time zone
    INTEGER, INTENT( IN ):: TSTEP              ! time step HHMMSS
    INTEGER, INTENT( IN ):: MXPDSRC            ! max. day- or hr-specific source
    LOGICAL, INTENT( IN ):: GETSIZES           ! true: get no. time steps  pols
    LOGICAL, INTENT( IN ):: GETCOUNT           ! true: get max no. srcs per time
    LOGICAL, INTENT( IN ):: FIRSTCALL          ! true: first call of a loop
    LOGICAL, INTENT( IN ):: DAYFLAG            ! true: day-specific wildfire data
    INTEGER,INTENT(INOUT):: SDATE              ! Julian starting date in TZONE
    INTEGER,INTENT(INOUT):: STIME              ! start time of data in TZONE
    INTEGER, INTENT(OUT) :: EDATE              ! Julian ending date in TZONE
    INTEGER, INTENT(OUT) :: ETIME              ! ending time of data in TZONE
    INTEGER, INTENT(OUT) :: EASTAT( NIPPA )     ! true: pol/act appears in data
    INTEGER, INTENT(OUT) :: SPSTAT( MXSPDAT )     ! true: special in data

    !.......  EXTERNAL FUNCTIONS

    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETTZONE

    !.......   SUBROUTINE PARAMETERS
    INTEGER      , PARAMETER :: NSEG     =    9         ! number of fields for ORL FIREDATA input format
    REAL         , PARAMETER :: TON2LB   = 2000.        ! pounds per short ton
    CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDORLFR'     !  program name

    !.......   Local list of bad sources to prevent duplicate writing of error
    !              messages
    CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: BADSRC( : )

    !.......   Local list of whether or not a warning was written for a pollutant or not
    LOGICAL, ALLOCATABLE, SAVE :: LCODWARN ( : )

    !.......  Local list for erroneous pollutant names in the file
    INTEGER, SAVE :: NBADPOLS
    CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: BADPOLS( : )

    !.......   Local list of FIPS start/end positions to facilitate
    !              faster lookups
    INTEGER, ALLOCATABLE, SAVE :: STARTSRC( : )
    INTEGER, ALLOCATABLE, SAVE :: ENDSRC( : )

    !.......   Local list of arrays for warning handling
    LOGICAL, ALLOCATABLE, SAVE :: WARNKEEP( : )     ! true: write warning for Keep = N
    LOGICAL, ALLOCATABLE, SAVE :: WARNMULT( : )     ! true: write warning for Multiple pollutants from a single pollutant

    !.......   Temporary read arrays
    CHARACTER(40)      SEGMENT( NSEG )     ! segments of line

    !.......   Local arrays
    REAL              , ALLOCATABLE, SAVE :: DTACBRN( : )        ! storing acre burned value (acre/day) for computing HFLUX
    REAL              , ALLOCATABLE, SAVE :: DTFUELD( : )        ! storing fuel loading value (tons/acre) for computing HFLUX

    INTEGER           , ALLOCATABLE, SAVE :: NSRCPDDAT( :,: )        ! counting number of sources per day/pollutant
    INTEGER           , ALLOCATABLE, SAVE :: IDXSD    ( : )          ! sorting index for CSRCDAYA
    CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAYA ( : )          ! unsorted source/day array
    CHARACTER(ALLLEN3), ALLOCATABLE, SAVE :: CSRCDAY  ( : )          ! sorted source/day array

    !.......   Other local variables
    INTEGER          H, HS, I, II, J, K, L, LL, N, NV, S, T, V1, V2        ! counters and indices
    INTEGER          L0, L1, L2, L3, L4, L5
    INTEGER          ES, NS, SS           ! end src, tmp no. src, start sourc

    INTEGER          D, SD
    INTEGER       :: N1 = 0
    INTEGER       :: N2 = 0

    INTEGER          CIDX                 ! CAS data index
    INTEGER          COD                  ! data index
    INTEGER          DAY                  ! tmp day of month
    INTEGER          FIP                  ! tmp co/st/cy code
    INTEGER, SAVE :: ICC = 0              ! tmp country code from header
    INTEGER          IOS                  ! i/o status
    INTEGER          IREC                 ! record counter
    INTEGER          JDATE                ! tmp Julian date
    INTEGER          JD                   ! Julian day number 1...365,366
    INTEGER          JTIME                ! tmp HHMMSS time
    INTEGER          ESTIME               ! tmp HHMMSS episode start time
    INTEGER          EETIME               ! tmp HHMMSS episode end time
    INTEGER       :: LFIP = 0             ! previous st/co FIPS code
    INTEGER, SAVE :: LOOPNO = 0           ! no. of loops
    INTEGER, SAVE :: MAXPTR               ! maximum time step reference pointer
    INTEGER, SAVE :: MINPTR               ! minimum time step reference pointer
    INTEGER          MONTH                ! tmp month number
    INTEGER, SAVE :: MXWARN               !  maximum number of warnings
    INTEGER, SAVE :: NWARN( 5 )           ! warnings counter
    INTEGER, SAVE :: NBADSRC = 0          ! no. bad sources
    INTEGER, SAVE :: NACRBND = 0          ! no. of acres burned var
    INTEGER, SAVE :: NFUELD  = 0          ! no. of fuel loading var
    INTEGER       :: NPOA    = 0          ! unused header number of pol/act
    INTEGER, SAVE :: NSRCDAY = 0          ! no. of src/day combos for computed vars
    INTEGER, SAVE :: NSTEPS  = 0          ! number of time steps
    INTEGER          PTR                  ! tmp time step pointer
    INTEGER       :: RDATE = 1980001      ! reference date: Jan 1, 1980
    INTEGER       :: RTIME = 0            ! reference time
    INTEGER, SAVE :: SDATESAV = 0         ! saved start date
    INTEGER, SAVE :: STIMESAV = 0         ! saved start time
    INTEGER, SAVE :: TDIVIDE  = 1         ! time step divisor
    INTEGER          WD                   ! tmp field width
    INTEGER          YEAR                 ! 4-digit year
    INTEGER       :: YR4 = 0              ! unused header year
    INTEGER          DZONE                ! time shift (ZONE-TZONE)
    INTEGER          ZONE                 ! source time zones

    REAL             TDAT                 ! temporary data values

    LOGICAL, SAVE :: TFLAG  = .FALSE.     ! true: use SCCs for matching with inv
    LOGICAL, SAVE :: DFLAG  = .FALSE.     ! true: dates set by data
    LOGICAL       :: EFLAG  = .FALSE.     ! TRUE iff ERROR
    LOGICAL       :: WARNOUT= .FALSE.     ! true: then output warnings
    LOGICAL, SAVE :: LFLAG = .FALSE.      ! true: output daily/hourly inv in local time
    LOGICAL, SAVE :: PRCHFX = .FALSE.     ! true: skip adding HFLUX due to precomputed heat flux
    LOGICAL       :: HFXFLAG= .FALSE.     ! true: adding HFLUX into a list
    LOGICAL       :: BNHRFLAG=.FALSE.     ! true: adding BEGHOUR into a list
    LOGICAL       :: ENHRFLAG=.FALSE.     ! true: adding ENDHOUR into a list
    LOGICAL, SAVE :: FIRSTCOUNT = .TRUE.    ! true: until after first time routine is called with GETCOUNT=TRUE
    LOGICAL, SAVE :: FIRSTIME = .TRUE.    ! true: first time routine called

    CHARACTER(256) :: BUFFER = ' '        ! src description buffer
    CHARACTER(300) :: LINE   = ' '        ! line buffer
    CHARACTER(300) :: MESG   = ' '        ! message buffer

    !.......  Temporary local character variables
    CHARACTER(FIPLEN3) CFIP          ! tmp co/st/cy code
    CHARACTER(CASLEN3) CDAT          ! tmp data name (*16)
    CHARACTER(NAMLEN3) CNAM,PNAM     ! tmp SMOKE name
    CHARACTER(NAMLEN3) CTMP          ! tmp data name (*16)
    CHARACTER(PLTLEN3) FCID          ! tmp facility ID (*15)
    CHARACTER(CHRLEN3) SKID          ! tmp stack ID (*15) = LocID
    CHARACTER(CHRLEN3) DVID          ! dummy device ID
    CHARACTER(CHRLEN3) PRID          ! dummy process ID
    CHARACTER(SCCLEN3) TSCC          ! tmp source category code (*10)
    CHARACTER(ALLLEN3) CSRC          ! tmp source string
    CHARACTER(ALLLEN3) CSRCD         ! tmp source/date string
    CHARACTER(ALLLEN3) TSRC          ! tmp source string
    CHARACTER( 8 )     DATE          ! tmp date string
    CHARACTER(NAMLEN3) PNAME         ! logical file name for data files

    !***********************************************************************
    !   begin body of program RDORLFR

    !.......  First time routine called
    IF( FIRSTIME ) THEN

        !.......  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
        END IF

        !.......  No time zone shift for AERMOD support
        MESG = 'Outputs local time daily and/or hourly inventories (No time shift)'
        LFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTPUT_LOCAL_TIME"', 2 )
        END IF

        !.......  Allocate memory for storing counting number of sources per day/pollutant
        !.......  Allocate memory for bad pollutant issues
        ALLOCATE( NSRCPDDAT( 366,NIPPA ),   &
                     BADSRC( NSRC ),        &
                   LCODWARN( NINVTBL ),     &
                    BADPOLS( NINVTBL ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NSRCPDDAT...BADPOLS', PROGNAME )

        NSRCPDDAT = 0          ! array
        LCODWARN = .FALSE.
        BADPOLS = ' '

        !.......  Create unique list of FIPS codes and other things
        CALL GENUSLST

        !.......  Build helper arrays for making searching faster
        ALLOCATE( STARTSRC( NINVIFIP ),    &
                    ENDSRC( NINVIFIP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ENDSRC', PROGNAME )
        STARTSRC = 0
        ENDSRC = 0
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

        !......  Open I/O API inventory HEATCONTENT file and store to use for
        !        computing HFLUX in PDAY intermediate output file.
        MESG = 'Reading HEATCONTENT data from inventory file...'
        CALL M3MSG2( MESG )

        !.......  Open I/O API inventory HEATCONTENT file and store
        ALLOCATE( HEATCONTENT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HEATCONENT', PROGNAME )

        CALL RDMAPPOL( NSRC, 1, 1, 'HEATCONTENT', HEATCONTENT )

        FIRSTIME = .FALSE.
        FIREFLAG = .TRUE.


        !.......  Initialize warnings counter
        NWARN = 0          ! array

    END IF      ! End first time subroutine is called

    !.......  For the first call in a loop of files, initialize variables
    IF( FIRSTCALL ) THEN
        MINPTR  = 99999999
        MAXPTR  = 0

        !.......  Set time step divisor
        TDIVIDE = 3600 * TSTEP / 10000

        !.......  If dates have been set by the data, set the number of steps
        !               steps
        IF( DFLAG ) THEN
            NSTEPS = 1+ SECSDIFF( SDATE,STIME,EDATE,ETIME )/TDIVIDE
            SDATESAV = SDATE
            STIMESAV = STIME
        END IF

        !.......  Set switch for printing errors only the first loop through all
        !         of the input files.  The second time through is indicated
        !         for the second time that FIRSTCALL is true.
        !.......  Reset loop counter if call is to get dimensions only (because
        !         this means it is the first call or daily or hourly)
        IF( GETSIZES ) LOOPNO = 0
        LOOPNO  = LOOPNO + 1
        WARNOUT = ( LOOPNO .EQ. 1 )

        !.......  Deallocate warning arrays
        IF( ALLOCATED( WARNKEEP ) ) DEALLOCATE( WARNKEEP, WARNMULT )
        ALLOCATE( WARNKEEP( NUNIQCAS ),    &
                  WARNMULT( NUNIQCAS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'WARNMULT', PROGNAME )
        WARNKEEP = .TRUE.
        WARNMULT = .TRUE.

    END IF

    !.......  For the second pass of this routine, allocate memory needed for calculating
    !           values from day-specific data on the fly.
    IF ( GETCOUNT .AND. FIRSTCOUNT ) THEN

        !.......  Determine how much memory is needed for allocating arrays
        !         This should be the maximum number of
        !         source/days of fuel load, and acres burned.
        DO D = 1, 366
            DO I = 1, NIPPA
                IF( EANAM(I)=='FUEL_LOAD' )   N1= N1+ NSRCPDDAT(D,I)
                IF( EANAM(I)=='ACRESBURNED' ) N2= N2+ NSRCPDDAT(D,I)
            END DO
        END DO
        N = MAX( N1, N2 )

        ALLOCATE(    IDXSD( N ),    &
                  CSRCDAYA( N ),    &
                   CSRCDAY( N ),    &
                   DTACBRN( N ),    &
                   DTFUELD( N ), STAT=IOS )      ! To store fuel load
        CALL CHECKMEM( IOS, 'DTFUELD', PROGNAME )

        IDXSD    = 0
        CSRCDAYA = ' '
        CSRCDAY  = ' '
        DTACBRN  = BADVAL3
        DTFUELD  = BADVAL3

        FIRSTCOUNT = .FALSE.

    END IF

    !.......  For the third pass of this routine, create sorted CSRCDAY
    !           routine.  This are used for caculating new values (e.g., PMC) from day-specific
    !           data on the fly.
    IF ( .NOT. GETSIZES .AND. .NOT. GETCOUNT ) THEN

        CALL SORTIC( NSRCDAY, IDXSD, CSRCDAYA )

        DO SD = 1, NSRCDAY
            K = IDXSD( SD )
            CSRCDAY( SD ) = CSRCDAYA( K )
        END DO

    END IF

    !.......  Loop through file and read it. In the first section, determine
    !           the minimum and maximum date. Use a reference date to do this. In
    !           the second section, determine the number of records per time
    !           step. In the third section, read and store the data.  When storing
    !           data, time step index is computed from the start date/time instead
    !           of the reference date/time so that the indexing will work properly.
    LFIP = 0
    IREC = 0
    TDAT = 0

    DO             !  Head of period-specific file read loop

        !.......  Read first line of file
        READ( FDEV, 93000, END=299 ) LINE
        IREC = IREC + 1

        L = LEN_TRIM( LINE )

        !.......  Skip blank lines
        IF( L == 0 ) CYCLE

        !.......  Scan for header lines and check to ensure all are set
        !               properly
        CALL GETHDR( 1, .TRUE., .FALSE., .FALSE., LINE, ICC, YR4, NPOA, IOS )

        !.......  Interpret error status
        IF( IOS .GT. 0 ) EFLAG = .TRUE.

        !.......  If a header line was encountered, go to next line
        IF( IOS .GE. 0 ) CYCLE

        !.......  Get lines
        CALL PARSLINE( LINE, NSEG, SEGMENT )

        !.......  Use the file format definition to parse the line into
        !               the various data fields
        CFIP = REPEAT( '0', FIPLEN3 )
        WRITE( CFIP( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) ICC          ! country code of FIPS
        CFIP( FIPEXPLEN3+2:FIPEXPLEN3+6 ) = ADJUSTR( SEGMENT( 1 )( 1:5 ) )      ! state/county code
        FIP    = STR2INT( CFIP )              ! FIP codes
        FCID   = SEGMENT( 2 )       ! fire ID
        SKID   = SEGMENT( 3 )       ! location ID
        TSCC   = SEGMENT( 4 )       ! SCC
        DVID   = ' '                ! dummy device id
        PRID   = ' '                ! dummy process id
        CDAT   = SEGMENT( 5 )       ! Pollutants(FUEL_LOAD, ACRESBURNED,,,)
        DATE   = SEGMENT( 6 )       ! Date of episode
        ESTIME = STR2INT( SEGMENT( 8 ) ) * 10000     ! episode start time
        EETIME = STR2INT( SEGMENT( 9 ) ) * 10000     ! episode end time

        IF( TSCC .NE. ' ' ) CALL PADZERO( TSCC )

        !.......  Conver CAS number to pollutant names if available
        I = INDEX1( CDAT, NINVTBL, ITCASA )
        IF( I > 0 ) CDAT = ITNAMA( I )

        !....... Check fire beginning and ending time format and print warning if necessary
        IF( EETIME > 230000 .OR. EETIME < 0 ) THEN
            MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //    &
                   ' Date: ' // DATE  // ' :: Fire ending' //    &
                   ' time not in the acceptable range of 0 to 23'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
            CYCLE
        END IF

        IF( ESTIME > 230000 .OR. ESTIME < 0  ) THEN
            MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //    &
                   ' Date: ' // DATE  // ' :: Fire starting' //    &
                   ' time not in the acceptable range of 0 to 23'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
            CYCLE
        END IF

        IF( ESTIME > EETIME ) THEN
            MESG = 'ERROR: Region: '// CFIP // ' SCC: ' // TSCC //    &
                   ' Date: ' // DATE  // ' :: Fire ending time' //    &
                   ' can not be earlier then a begining time'
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
            CYCLE
        END IF

        !.......  Check date format
        IF( DATE( 3:3 ) /= '/' .OR.    &
            DATE( 6:6 ) /= '/'      ) THEN
            MESG = 'ERROR: Incorrect date format ( MM/DD/YY ) :' //    &
                   ' Region: '// CFIP // ' SCC: ' // TSCC //    &
                   ' Date: ' // DATE
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
            CYCLE
        END IF

        !.......  Check and Set emissions values
        TDAT = STR2REAL( SEGMENT( 7 ) )             ! Day-specific total emission

        IF ( TDAT .LT. 0.0 )  THEN
            EFLAG = .TRUE.
            WRITE( MESG,94030 ) 'ERROR: Bad data value "',    &
                TDAT, '" of ' // TRIM( CDAT ) //    &
                ' for region ' // CFIP //    &
                ' and SCC ' // TSCC // ' on date '// DATE
            CALL M3MESG( MESG )
            CYCLE          ! to head of read loop
        END IF

        !.......  Counting the number of times precomputed HFLUX
        !               values appear in the input file
        IF( GETSIZES .AND. CDAT == 'HFLUX' ) PRCHFX = .TRUE.

        !.......  Counting and adding HFLUX, BEGHOUR, and ENDHOUR
        !               building a list of source characteristics and store
        IF( GETCOUNT .AND. .NOT. PRCHFX ) THEN

            IF( HFXFLAG  ) CDAT = 'HFLUX'
            IF( BNHRFLAG ) CDAT = 'BEGHOUR'
            IF( ENHRFLAG ) CDAT = 'ENDHOUR'

            HFXFLAG = .FALSE.
            IF( CDAT == 'ACRESBURNED' ) THEN
                NACRBND = NACRBND + 1
                HFXFLAG = .TRUE.        ! indicating adding HFLUX
                BACKSPACE( FDEV )
            END IF

            BNHRFLAG = .FALSE.
            IF( CDAT == 'HFLUX' ) THEN
                BNHRFLAG = .TRUE.        ! indicating adding BEGHOUR
                BACKSPACE( FDEV )
            END IF

            ENHRFLAG = .FALSE.
            IF( CDAT == 'BEGHOUR' ) THEN
                ENHRFLAG = .TRUE.        ! indicating adding ENDHOUR
                BACKSPACE( FDEV )
            END IF

            IF( CDAT == 'FUEL_LOAD' ) NFUELD = NFUELD + 1

        END IF

        !.......  Set Julian day from MMDDYY8 SAS format
        MONTH = STR2INT( DATE( 1:2 ) )
        DAY   = STR2INT( DATE( 4:5 ) )
        YEAR  = YEAR4( STR2INT( DATE( 7:8 ) ) )

        JD = JULIAN( YEAR, MONTH, DAY )
        JDATE = 1000 * YEAR + JD
        JTIME = 0


        !.......  Local time shift flag (county-specific)
        IF( .NOT. LFLAG ) THEN

            !.......  Set time zone number
            ZONE = GETTZONE( CFIP )

            !.......  If daily emissions are not in the output time zone, print
            !                   warning
            IF( GETCOUNT ) THEN
                IF( WARNOUT .AND. ( ZONE .NE. TZONE ) .AND.    &
                   ( NWARN( 1 ) .LE. MXWARN )               ) THEN
                    WRITE( MESG,94010 ) 'WARNING: Time zone ', ZONE,    &
                      'in day-specific file at line of pollutant ' //    &
                      TRIM( CDAT ) // ' on ' // TRIM( DATE ) //    &
                      ' does not match output time zone', TZONE
                    CALL M3MESG( MESG )
                    NWARN( 1 ) = NWARN( 1 ) + 1
                END IF
            END IF

            DZONE = ZONE - TZONE

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

        !.......  Store maximum time step number as compared to rference
        IF( PTR + 23 .GT. MAXPTR ) MAXPTR = PTR + 23

        !.......  If FIPS code is not the same as last time, then
        !               look it up and get indidies
        IF( FIP .NE. LFIP ) THEN
            J = FINDC( CFIP, NINVIFIP, INVCFIP )
            IF( J .LE. 0 ) THEN
                WRITE( MESG,93000 ) 'INTERNAL ERROR: Could not '//    &
                       'find FIPS code ' // CFIP // 'in internal list.'
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
            !                   if reading the SCC in helps (needed for IDA format)
            IF( J .LE. 0 ) THEN

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

        !.......  Look up pollutant name in unique sorted array of
        !               Inventory pollutant names
        CDAT  = SEGMENT( 5 )
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
            !                   special integer so can ID these records later.
            IF( CIDX .GT. 0 ) THEN
                SPSTAT( CIDX ) = CIDX
                COD = CODFLAG3 + CIDX

            !......  If not in list of special names, check to see
            !                  if it's a SMOKE pollutant name (intermediate name)
            ELSE IF ( CIDX .LE. 0 ) THEN
                IF( WARNOUT .AND. NWARN( 2 ) .LE. MXWARN ) THEN
                    WRITE( MESG,94010 )    &
                     'WARNING: Skipping pollutant "'// TRIM(CDAT)//    &
                     '" at line', IREC, '- not in Inventory Table'
                    CALL M3MESG( MESG )
                    NWARN( 2 ) = NWARN( 2 ) + 1
                END IF
                CYCLE          !  to head of loop

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
            ELSE IF ( UCASNKEP(CIDX) .GT. 1 .AND.    &
                      WARNMULT(CIDX)              ) THEN
                WARNMULT( CIDX ) = .FALSE.
                IF( GETSIZES ) THEN
                    WRITE( MESG,94010 )                                &
                      'WARNING: Skipping all lines for pollutant "'//  &
                      TRIM( CDAT )// '" because Inventory Table '//    &
                      'splits it into',UCASNKEP(CIDX),'pollutants.'//  &
                      CRLF()//BLANK10//'The SMOKE code needs to '//    &
                      'be enhanced to support this approach for '//    &
                      'day- and hour-specific data.'
                    CALL M3MESG( MESG )
                END IF
                CYCLE
            END IF

            !......  Get Inventory Data SMOKE name from Inventory Table arrays/indices
            CNAM = ITNAMA( SCASIDX( UCASIDX( CIDX ) ) )

            !......  Look up SMOKE name in list of annual EI pollutants
            COD = INDEX1( CNAM, NIPPA, EANAM )

            !......  Check to ensure that it handles NOI and NONHAP pollutants
            !                  while combining VOC + HAPs
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
                            PNAM = CNAM(1:L+1) // 'NONHAP' // CNAM(L+2:LL)
                        ELSE
                            PNAM = 'NONHAP' // TRIM( CNAM )
                        END IF
                        COD = INDEX1( PNAM, NIPPA, EANAM )
                    END IF

                END IF

            END IF

            !......  Check to ensure that the SMOKE intermediate name
            !                  set by the Inventory Table is actually in the annual
            !                  inventory.  If not, write warning message and cycle.
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

        END IF      ! if cidx le 0 or not

        !.......  Count the number of sources per day  pollutant/variable
        !.......  This will give us how many source/date combos there are for
        !               any variables, including HFLUX
        NSRCPDDAT( JD, COD ) = NSRCPDDAT( JD, COD ) + 1

        !.......  If only getting dates and pollutant information, go
        !               to next loop iteration
        IF( GETSIZES ) CYCLE

        !.......  Determine time step pointer based on actual start time
        PTR = SECSDIFF( SDATESAV,STIMESAV,JDATE,JTIME )/TDIVIDE + 1

        !.......  Skip record if it is out of range of output file
        !.......  NOTE - this is only useful if reading only part of
        IF( PTR .LT. 1 .OR. PTR .GT. NSTEPS ) CYCLE

        !.......  Count estimated record count per time step
        DO T = PTR, MIN( PTR + 23, NSTEPS )
            MXPDPT( T ) = MXPDPT( T ) + 1
        END DO

        !.......  Store variable values.  Only need to do this on the the second
        !               pass.  Need to do this before the third pass through the data because
        !               that is when the calculation is made.
        IF( GETCOUNT .AND. .NOT. PRCHFX ) THEN                  ! No precomputed formula/heat flux
            IF( ( HFXFLAG .OR. CDAT == 'FUEL_LOAD' ) ) THEN          ! Acres burned value or fuel load value

                !.......  Figure out which source/day this is for storing in correct source/day
                !.......  This code does *not* assume that the data have been sorted first.
                CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID,    &
                              TSCC, DATE, POLBLNK3, CSRCD )

                !.......  Build unsorted arrays of source/days and emissions for calculating formula
                SD = 0
                IF( NSRCDAY > 0 ) THEN
                    SD= INDEX1( CSRCD, NSRCDAY, CSRCDAYA )
                END IF

                IF( SD <= 0 ) THEN
                    NSRCDAY = NSRCDAY + 1
                    SD = NSRCDAY
                    CSRCDAYA( SD ) = CSRCD
                    IDXSD   ( SD ) = SD
                END IF

                IF( HFXFLAG ) DTACBRN( SD ) = TDAT                        ! storing acres burned
                IF( CDAT == 'FUEL_LOAD' ) DTFUELD( SD ) = TDAT         ! storing fuel load

            END IF
        END IF            ! Second pass only

        !.......  If only counting records per time step, go to next loop
        !               iteration
        IF( GETCOUNT ) CYCLE

        !.......  Store source ID
        LPDSRC( S ) = .TRUE.

        !.......  Computing HFLUX, BEGHOUR, ENDHOUR (as a default)
        IF( CDAT == 'HFLUX' .AND. .NOT. PRCHFX ) THEN

            CALL BLDCSRC( CFIP, FCID, SKID, DVID, PRID, TSCC, DATE,    &
                          POLBLNK3, CSRCD )

            !.......  Lookup source/date string in master list to get position
            SD = FINDC( CSRCD, NSRCDAY, CSRCDAY )
            K  = IDXSD( SD )

            IF( DTFUELD( K ) < AMISS3 ) THEN
                LL = LEN_TRIM( CSRCD )
                CALL FMTCSRC( CSRCD, 6, BUFFER, L2 )

                MESG = 'ERROR: Missing value of '//    &
                       'fuel load for source:'//    &
                       CRLF() // BLANK10 // BUFFER( 1:L2 ) //    &
                       ' on date ' // CSRCD( LL-7: LL )
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.

            ELSE
            !.......  Compute Heat Flux value
                TDAT = DTACBRN( K ) * DTFUELD( K ) *     &            ! computing HFLUX (BTU/day)
                       HEATCONTENT( S ) * TON2LB                 ! HEATCONTENT(BTU/lb)=8000

            END IF

        END IF

        IF( CDAT == 'BEGHOUR' ) TDAT = REAL( ESTIME )      ! storing BEGHOUR
        IF( CDAT == 'ENDHOUR' ) TDAT = REAL( EETIME )      ! storing ENDHOUR

        !.......  Record needed data for this source and time step
        H = 0
        DO T = PTR, MIN( PTR + 23, NSTEPS )
            H = H + 1
            NPDPT( T ) = NPDPT( T ) + 1
            NPDPTP( T,COD ) = NPDPTP( T,COD ) + 1

            HS = NPDPT( T )

            IF( HS .LE. MXPDSRC ) THEN

                IDXSRC( HS,T ) = HS
                SPDIDA( HS,T ) = S
                CODEA ( HS,T ) = COD
                CIDXA ( HS,T ) = CIDX
                EMISVA( HS,T ) = TDAT      ! Store data in emissions
                DYTOTA( HS,T ) = TDAT

            END IF
        END DO

    END DO         ! Main read loop of day-specific data

299 CONTINUE       ! Exit from read loop

    !.......  Warning messages for HFLUX
    IF( GETCOUNT ) THEN

        IF( PRCHFX ) THEN
            MESG = 'WARNING: Skipping internal heat flux '//    &
                   'computation due to the existence of '//    &
                   'precomputed HFLUX in PTDAY file'
            CALL M3MSG2( MESG )
        END IF

        IF( NACRBND .NE. NFUELD .AND. .NOT. PRCHFX ) THEN
            MESG = 'ERROR: No of ACRESBURNED and FUEL_LOAD are' //    &
                   ' not matched for heat flux computation.'
            CALL M3MSG2( MESG )
        END IF

        IF( NACRBND < 1 .AND. .NOT. PRCHFX ) THEN
            MESG = 'FATAL ERROR: No ACRESBURNED data are available'    &
                   // ' for internal heat flux computation.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( NFUELD < 1 .AND. .NOT. PRCHFX ) THEN
            MESG = 'FATAL ERROR: No FUEL_LOAD data are available'    &
                   // ' for internal heat flux computation.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

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
94030 FORMAT( 10( A, :, F15.0, :, 1X ) )

94020 FORMAT( I6.6 )

END SUBROUTINE RDORLFR
