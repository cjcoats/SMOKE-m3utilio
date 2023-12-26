
PROGRAM MET4MOVES

    !***********************************************************************
    !  program body starts at line 205
    !
    !  DESCRIPTION:
    !       Creates county-based 24-hour temperature profiles and relative
    !       humidity based on gridded meteorology data.
    !       Requires hourly barometric pressure values to compute RH.
    !       using barometric pressure and dew point. Tries to
    !       account for missing meteorology data as much as possible.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  REVISION  HISTORY:
    !       Created  4/10: by B.H. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
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

    !.......   MODULES for public variables
    USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX,          &
                        DAILY, WEEKLY, MONTHLY, EPISLEN,            &
                        NREFF, FMREFSORT, NFUELC, FMREFLIST

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GRDNM, NGRID, NCOLS, NROWS, COORD, GDTYP,    &
                       P_ALP, P_BET, P_GAM, XCENT, YCENT, OFFLAG

    USE MODSURG, ONLY: NSRGREC, IDXSRGA, SFIPSA, SCELLA, SSRGIDA,   &
                       IDXSRGB, NSRGFIPS, NSRGS, SRGFIPS, SRGLIST,  &
                       FIPCELL, SRGFRAC, SRGCSUM, NCELLS, SFRACA,   &
                       SRGFMT, NTSRGDSC, SRGFNAM, SRGFCOD, SRGFREG

    !.......  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: LSTCYPOP, STCYPOPYR, USEDAYLT,               &
                       NCOUNTRY, NSTATE,   NCOUNTY,                 &
                       CTRYCOD,  STATCOD,  CNTYCOD,                 &
                       CTRYNAM,  STATNAM,  CNTYNAM,                 &
                       CTRYPOPL, STATPOPL, CNTYPOPL,                &
                       CNTYTZON, CNTYTZNM

    !.......   This module is the derived meteorology data for emission factors
    USE MODMET, ONLY: MINTEMP, MAXTEMP, RHTBIN, NRHTBIN,            &
                      TKHOUR, NTKHOUR, MAXTSRC, MINTSRC,            &
                      MAXTFUEL, MINTFUEL, TKFUEL, FUELIDX,          &
                      NFUEL, MAXTDAY, MINTDAY, FUELCNTY

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables

    INCLUDE 'CONST3.EXT'        !  physical and mathematical constants (i/O API)

    !.......  EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: DSCM3GRD
    INTEGER, EXTERNAL :: GETIFDSC
    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: FINDCFIRST

    !.......   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'MET4MOVES'      !  program name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag
    INTEGER      , PARAMETER :: MXVAR    = 20

    !.......   LOCAL VARIABLES and their descriptions:
    CHARACTER(20)  SEGMENT( 20 )              ! parsed input line

    !.......   Gridded meteorology data (dim: NGRID)
    REAL   , ALLOCATABLE :: TA( : )       !  one layer of temperature
    REAL   , ALLOCATABLE :: QV( : )       !  water vapor mixing ratio
    REAL   , ALLOCATABLE :: PRES( : )     !  pressure
    REAL                    TKREFHR( 24 )     !  ref county temp profile
    REAL   , ALLOCATABLE :: MAXTCELL( : )     !  daily gridded max temp
    REAL   , ALLOCATABLE :: MINTCELL( : )     !  daily gridded min temp

    !.......   Ungridding Matrix
    INTEGER, ALLOCATABLE :: UMAT( : )     ! Contiguous ungridding matrix

    !.......  Allocatable per-source arrays
    INTEGER, ALLOCATABLE :: DAYBEGT  ( : )     ! daily start time HHMMSS
    INTEGER, ALLOCATABLE :: DAYENDT  ( : )     ! daily end time HHMMSS
    LOGICAL, ALLOCATABLE :: LDAYSAV  ( : )     ! true: src uses daylight time
    INTEGER, ALLOCATABLE :: TZONES   ( : )     ! county-specific time zones
    INTEGER                 SRGIDS  ( MXVAR )     ! listing of surrogates

    !.......  Allocatable arrays for met data
    INTEGER, ALLOCATABLE :: METDAYS( : )     ! dimension: nsteps in episode,
        ! value indicates which met data file covers that hour
    CHARACTER(256),     ALLOCATABLE :: METLIST( : )     ! listing of met file names
    CHARACTER(FIPLEN3), ALLOCATABLE :: CNTYSRC( : )     ! FIPS code and averaging values

    !.......   File units and logical names:
    INTEGER      ADEV      ! unit number for county fuelmonth file
    INTEGER      CDEV      ! unit number for co/st/cy file
    INTEGER      DDEV      ! unit number for daily group file
    INTEGER      EDEV      ! unit number for episode group file
    INTEGER      GDEV      ! tmp unit number for individual surrogate file
    INTEGER      IDEV      ! tmp unit number if ENAME is map file
    INTEGER      LDEV      ! unit number for log file
    INTEGER      MDEV      ! unit number for monthly group file
    INTEGER      PDEV      ! unit number for speeds summary file (SPDSUM)
    INTEGER      QDEV      ! unit number for surrogate(s) file
    INTEGER      SDEV      ! unit number for surrogate(s) file
    INTEGER      TDEV      ! unit number for meteorology list file
    INTEGER      WDEV      ! unit number for weekly group file
    INTEGER      ODEV      ! unit number for tmp surrogate file
    INTEGER      ODEV1     ! unit number for MOVES-ready output file
    INTEGER      ODEV2     ! unit number for RH MOVES-ready output file
    INTEGER      XDEV      ! unit number for mobile x-ref file


    CHARACTER(16) ONAME       ! logical name for daily gridded min/max temp output for SMOKE Movesmrg
    CHARACTER(16) DNAME       ! logical name for daily output ungridded hourly temps
    CHARACTER(16) INAME       ! tmp name for inven file of unknown fmt
    CHARACTER(16) METNAME     ! logical name for meteorology files
    CHARACTER(16) MNAME       ! logical name for monthly output hourly temps
    CHARACTER(16) PNAME       ! logical name for episode output hourly temps
    CHARACTER(16) UNAME       ! logical name for ungridding-matrix input file
    CHARACTER(16) WNAME       ! logical name for weekly output hourly temps

    !.......   Other local variables:
    INTEGER    C, I, NC, J, K, L, N, NF, NR, NS, S, T, TT, V      ! Counters and pointers

    INTEGER    EPI_SDATE          ! episode start date from E.V. (YYYYDDD)
    INTEGER    EPI_STIME          ! episode start time from E.V. (HHMMSS)
    INTEGER    EPI_RUNLEN         ! episode duration from E.V. (HHMMSS)
    INTEGER    EPI_NSTEPS         ! episode number of time steps
    INTEGER    EPI_EDATE          ! episode ending date based on ERUNLEN
    INTEGER    EPI_ETIME          ! episode ending time based on ERUNLEN

    INTEGER    CURCNTY         ! tmp current processing county
    INTEGER    DAY             ! tmp day of week number
    INTEGER    EDATE           ! ending input date counter (YYYYDDD) in GMT
    INTEGER    ETIME           ! ending input time counter (HHMMSS)  in GMT
    INTEGER    FILENUM         ! file number of current meteorology file
    INTEGER    IOS             ! temporary I/O status
    INTEGER    HOURIDX         ! current hour of the day
    INTEGER    JDATE           ! input date counter (YYYYDDD) in GMT
    INTEGER    JTIME           ! input time counter (HHMMSS)  in GMT
    INTEGER    LDATE,PDATE     ! date from previous loop iteration
    INTEGER    ODATE           ! output date for counties
    INTEGER    OTIME           ! output time
    INTEGER    SYEAR           ! year of start date
    INTEGER    SMONTH          ! month of start date
    INTEGER    SDAY            ! date of month (start)
    INTEGER    EYEAR           ! year of end date
    INTEGER    EMONTH          ! month of end date
    INTEGER    EDAY            ! date of month (end)
    INTEGER    MONTH           ! tmp month
    INTEGER    CURMONTH        ! current month
    INTEGER    PRVCMONTH       ! previous month
    INTEGER    FUELMONTH       ! current fuelmonth
    INTEGER    PRVFMONTH       ! previous fuelmonth
    INTEGER    METNGRID        ! no. grid cells in met data
    INTEGER    MINNORH         ! min no of RH datapoints for averaging RH by tempbin
    INTEGER    NLINES          ! no. lines in met list file
    INTEGER    NSRC            ! no. source (=counties)
    INTEGER    NVARS           ! no. surrogates
    INTEGER    NMON            ! no of fuelmonth per refcounty
    INTEGER    NFMON           ! tmp no of fuelmonths per refcounty
    INTEGER    NSTEPS          ! number of time steps to process temperature data
    INTEGER    MXTBIN          ! number of possible max temperature bins
    INTEGER    POS             ! position in time step loop
    INTEGER    PDTEMP          ! temp increment for rateperdistance lookup table
    INTEGER    PVTEMP          ! temp increment for ratepervehicle lookup table
    INTEGER    PPTEMP          ! temp increment for rateperprofile lookup table
    INTEGER    RDATE           ! date to read met file
    INTEGER    SDATE           ! output start date
    INTEGER    STIME           ! output start time
    INTEGER    TDATE           ! temporary julian date
    INTEGER    TTIME           ! temporary time
    INTEGER    SRGNCOLS        ! surrogate no of cols
    INTEGER    SRGNROWS        ! surrogate no of rows
    INTEGER    TMPMNTH         ! temporary month
    INTEGER    TSPREAD         ! time spread: difference between TZMAX and TZMIN
    INTEGER    TZONE           ! zone to determine output days
    INTEGER    TZMIN           ! minimum time zone in inventory
    INTEGER    TZMAX           ! maximum time zone in inventory
    INTEGER    WDATE           ! output date for weekly counties

    REAL       TEMPBIN                 ! temperature buffer bin
    REAL       MAXTREF, MINTREF        ! ref county min/max temperatures

    LOGICAL :: EFLAG    = .FALSE.      !  true: error found
    LOGICAL :: HFLAG    = .FALSE.      !  true: output specific humidity instead of default RH
    LOGICAL :: FIRSTIME = .FALSE.      !  true: first time
    LOGICAL :: GRID_ERR = .FALSE.      !  true: error found in grid settings
    LOGICAL :: FILEOPEN = .FALSE.      !  true: met file is open
    LOGICAL :: FND_DATA = .FALSE.      !  true: found met data for this hour
    LOGICAL :: ALT_DATA = .FALSE.      !  true: using alternate data for this hour

    CHARACTER(FIPLEN3) CFIP            ! tmp inventory county
    CHARACTER(FIPLEN3) REFCOUNTY       ! tmp ref. county
    CHARACTER(FIPLEN3) INVCOUNTY       ! tmp inv. county

    CHARACTER(16)      COORUNIT        !  coordinate system projection units
    CHARACTER(80)      GDESC           !  grid description
    CHARACTER(256)     CFNAME          !  fuelmonth input file
    CHARACTER(16)      SRG_CNTRY       !  surrogate country
    CHARACTER(256)     LINE            !  line buffer
    CHARACTER(256)     FULLNAME        !  full file name
    CHARACTER(512)     CMCXREF         !  tmp physical file name of MCXREF
    CHARACTER(512)     CMFMREF         !  tmp physical file name of MFMREF
    CHARACTER(512)     METFILE         !  tmp physical file name
    CHARACTER(512)     PREVFILE        !  previous physical file name
    CHARACTER(NAMLEN3) TVARNAME        !  temperature variable name
    CHARACTER(NAMLEN3) PRESNAME        !  pressure variable name
    CHARACTER(NAMLEN3) MIXNAME         !  mixing ratio name
    CHARACTER(200)     TEMPDIR         !  directory for output files
    CHARACTER(300)     MESG            !  message buffer

    !***********************************************************************
    !   begin body of program MET4MOVES

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Get file name for country, state, and county file, with time zones
    CDEV = PROMPTFFILE(&
           'Enter logical name for COUNTRY, STATE, AND ' //&
           'COUNTY file', .TRUE., .TRUE., 'COSTCY', PROGNAME )

    !.......  Open input surrogates file
    QDEV = PROMPTFFILE(&
           'Enter logical name for Surrogate Description file',&
           .TRUE., .TRUE., 'SRGDESC', PROGNAME )

    !.......  Open met list file
    TDEV= PROMPTFFILE(&
           'Enter logical name for METEOROLOGY LIST file'&
           ,.TRUE., .TRUE., 'METLIST', PROGNAME )

    !.......  Open mobile county x-ref file to determine representative counties
    XDEV = PROMPTFFILE(&
             'Enter logical name for MCXREF cross-reference file',&
             .TRUE., .TRUE., 'MCXREF', PROGNAME )
    CALL ENVSTR( 'MCXREF', MESG, ' ', CMCXREF, IOS )

    !.......  Obtain episode settings from the environment...

    !.......  Get temperature increments for rateperdistance lookup table
    MESG = 'Temperature increment for rateperdistance'
    PDTEMP = ENVINT( 'PD_TEMP_INCREMENT', MESG, 5, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "PD_TEMP_INCREMENT"', 2 )
    END IF

    !.......  Get temperature increments for ratepervehicle lookup table
    MESG = 'Temperature increment for ratepervehicle'
    PVTEMP = ENVINT( 'PV_TEMP_INCREMENT', MESG, 5, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "PV_TEMP_INCREMENT"', 2 )
    END IF

    IF( PDTEMP /= PVTEMP ) THEN
        MESG = 'WARNING: Current version does not support different'&
            // ' PD_TEMP_INCREMENT and PV_TEMP_INCREMENT setting'//&
            CRLF() // BLANK10 // ':: Reset RV_TEMP_INCREMENT to RD_TEMP_INCREMENT'
        CALL M3MESG( MESG )
        PVTEMP = PDTEMP        ! reset RVTEMP to PDTEMP
    END IF

    !.......  Get temperature increments for ratepervehicle lookup table
    MESG = 'Temperature increment for rateperprofile'
    PPTEMP = ENVINT( 'PP_TEMP_INCREMENT', MESG, 10, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "PP_TEMP_INCREMENT"', 2 )
    END IF

    !.......  Get episode starting date and time and ending date
    MESG = 'Temperature buffer bin'
    TEMPBIN = ENVREAL( 'TEMP_BUFFER_BIN', MESG, 10.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TEMP_BUFFER_BIN"', 2 )
    END IF

    !.......  Define minimum no of data point for calculating avg RH by tempbin
    MESG = 'Minimum no of data points for averaging RH by temperature bin'
    MINNORH = ENVINT( 'MIN_NO_RH_BY_TEMPBIN', MESG, 1, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_NO_RH_BY_TEMPBIN"', 2 )
    END IF

    !.......  Define type of humidity (RH or Specific Humidity)
    MESG = 'Use Specific Humidity ouput or not.'
    HFLAG = ENVYN( 'SPECIFIC_HUMIDITY_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SPECIFIC_HUMIDITY_YN"', 2 )
    END IF

    !.......  Get episode starting date and time and ending date
    MESG = 'Episode start date (YYYYDDD)'
    EPI_SDATE = ENVINT( 'STDATE', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "STDATE"', 2 )
    END IF

    MESG = 'Episode start time (HHMMSS)'
    EPI_STIME = ENVINT( 'STTIME', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "STTIME"', 2 )
    END IF

    MESG = 'Episode end date (YYYYDDD)'
    EPI_EDATE = ENVINT( 'ENDATE', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ENDATE"', 2 )
    END IF

    MESG = 'Episode end time (HHMMSS)'
    EPI_ETIME = ENVINT( 'ENDTIME', MESG, 230000, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ENDTIME"', 2 )
    END IF

    !.......  Find the total number of time steps
    EPI_NSTEPS = 1 + SECSDIFF( EPI_SDATE,EPI_STIME,EPI_EDATE,EPI_ETIME ) / 3600

    !.......  Get the time zone for output of the emissions
    TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_TMPDIR"', 2 )
    END IF

    !.......  Get the name of the temperature variable
    MESG = 'Temperature variable name'
    CALL ENVSTR( 'TVARNAME', MESG, 'TEMP2', TVARNAME, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TVARNAME"', 2 )
    END IF

    !.......  Set default names for additional variables
    PRESNAME = 'PRES'
    MIXNAME  = 'QV'

    !.......  Read the surrogates header and initialize the grid description
    !.......  Open/read SRGDESC file
    CALL RDSRGDESC( QDEV )

    !.......  Open output file for temporary combined surrogate file
    FULLNAME = 'TMP_COMBINED_SRG.txt'
    IF( .NOT. SETENVVAR( 'TMP_SRG_FILE', FULLNAME )) THEN
        MESG = 'Could not set logical name of file ' // TRIM( FULLNAME )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    ODEV = GETEFILE( 'TMP_SRG_FILE', .FALSE., .TRUE., PROGNAME )

    !.......  Get name of surrogate IDs to use
    MESG = 'Country of surrogate files'
    CALL ENVSTR( 'SRG_COUNTRY', MESG, 'USA', SRG_CNTRY, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SRG_COUNTRY"', 2 )
    END IF

    MESG = 'Spatial surrogate IDs'
    IF( .NOT. INTLIST( 'SRG_LIST', MESG, 20, NVARS, SRGIDS ) ) THEN
        MESG = 'Could not read list of surrogate IDs'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Prompt for and open I/O API output file(s)...
    DO I = 1, NVARS

        L = 0
        DO K = 1,NTSRGDSC
            CALL UPCASE( SRGFREG( K ) )
            IF( SRG_CNTRY    == SRGFREG( K ) ) THEN
                IF( SRGIDS( I )  == SRGFCOD( K ) ) L = K
            END IF
        END DO

        IF( L < 1 ) CYCLE

        CALL GETENV( 'SRGPRO_PATH', TEMPDIR )
        WRITE( FULLNAME, '(3A)' ) TRIM( TEMPDIR ),'/',TRIM( SRGFNAM( L ) )

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( 'TMP_SRG', FULLNAME )) THEN
            MESG = 'Could not set logical name of file ' // TRIM( FULLNAME )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  Get the number of lines in the surrogate description file desription file
        GDEV = GETEFILE( 'TMP_SRG',.TRUE., .TRUE., PROGNAME )
        NLINES = GETFLINE( GDEV, 'Surrogate file' )
        DO J = 1, NLINES
            READ(  GDEV,'(A)', END=999, IOSTAT=IOS ) LINE
            IF( LINE( 1:5 ) =='#GRID' ) THEN
                WRITE( ODEV,'(A)' ) TRIM( LINE )
            ELSE
                IF( BLKORCMT( LINE ) ) CYCLE
                CALL PARSLINE( LINE, 20, SEGMENT )
                IF( SRGIDS( I ) == STR2INT( SEGMENT( 1 ) ) ) THEN
                    WRITE( ODEV,'(A)' ) TRIM( LINE )
                END IF
            END IF
999     END DO

        CLOSE( GDEV )

    END DO

    CLOSE( ODEV )

    SDEV = GETEFILE( 'TMP_SRG_FILE',.TRUE., .TRUE., PROGNAME )

    !.......  Also get the format of the surrogates file.
    CALL RDSRGHDR( .FALSE., SDEV, SRGFMT )
    SRGNCOLS = NCOLS
    SRGNROWS = NROWS

    MESG = 'NOTE: Input surrogates are ' // TRIM( SRGFMT ) // ' format.'
    CALL M3MSG2( MESG )

    !.......  Obtain settings from the environment...
    !.......  Get grid name from the environment and read grid parameters
    IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D,    &
         COORUNIT, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,      &
         YCENT3D, XORIG3D, YORIG3D, XCELL3D,                &
         YCELL3D, NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )
    IF ( EFLAG ) THEN
        MESG = 'Problem with gridded input data.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Allocate memory for and read the gridding surrogates file,
    !           extracting data for a subgrid, if necessary
    CALL M3MSG2( 'Processing gridding surrogate(s)...' )

    CALL RDSRG( .FALSE., SDEV, SRGFMT, SRGNROWS, SRGNCOLS )

    MESG = 'NOTE: A list of surrogates used in the modeling'
    CALL M3MESG( MESG )

    DO I = 1,NSRGS
        WRITE( MESG,94010 ) 'Surrogate ID ::', SRGLIST( I )
        CALL M3MESG( MESG )
    END DO

    !.......  Read region codes file
    CALL RDSTCY( CDEV, NSRGFIPS, SRGFIPS )

    !.......  Read county x-ref file to determine representative counties
    CALL RDMXREF( XDEV, NSRGFIPS, SRGFIPS  )

    !.......  open county-specific fuelmonth input files in file
    MESG = 'Enter logical name for a reference county-' //&
           'specific FuelMonth input file'
    ADEV = PROMPTFFILE( MESG,.TRUE.,.TRUE.,'MFMREF', PROGNAME )
    CALL ENVSTR( 'MFMREF', MESG, ' ', CMFMREF, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MFMREF"', 2 )
    END IF

    CALL RDFMREF( ADEV )

    !.......  define a number of target sources (=inventory counties)
    NSRC = NINVC

    !.......  Allocate arrays for county time zone
    ALLOCATE( TZONES( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TZONES', PROGNAME )
    TZONES = 0

    !.......  Assign time zone to inventory counties
    DO I = 1, NSRC
        CFIP = MCREFSORT( I,1 )
        J = FINDC( CFIP, NCOUNTY, CNTYCOD )
        IF( J < 1 ) THEN
            WRITE( MESG,94010 ) 'ERROR: Could not find time zone '//&
                 'for county :' // CFIP // ' from COSTCY file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ELSE
            TZONES( I ) = CNTYTZON( J )
        END IF
    END DO

    !.......  Define the minimum and maximum time zones in the inventory
    TZMIN = MINVAL( TZONES )
    TZMAX = MAXVAL( TZONES )

    !.......  Adjust TZMIN for possibility of daylight savings
    TZMIN = MAX( TZMIN - 1, 0 )

    !.......  Calculate time spread based on min and max time zone
    TSPREAD = TZMAX - TZMIN

    !.......  Calculate required starting date and time based on episode settings
    SDATE = EPI_SDATE

    !.......  For MOVES model, earliest time required will be 12 a.m. in time zone closest to GMT
    STIME = TZMIN - TZONE              ! starting time in output time zone

    !.......  Make sure the starting time is between 0 and 23
    IF( STIME < 0 ) THEN
        STIME = STIME + 24
    ELSE
        STIME = MOD( STIME, 24 )
    END IF
    STIME = STIME*10000

    !.......  Calculate required ending date and time based on episode settings
    EDATE = EPI_EDATE
    ETIME = EPI_ETIME

    !.......  For MOVES model, latest time required will be 11 p.m. in time zone farthest from GMT
    ETIME = TZMAX - TZONE - 1            ! ending time in output time zone

    !.......  Make sure the ending time is between 0 and 23
    IF( ETIME < 0 ) THEN
        ETIME = ETIME + 24
    ELSE
        ETIME = MOD( ETIME, 24 )
    END IF
    ETIME = ETIME*10000

    !.......  If the episode ending time is later than calculated end time,
    !           set the ending date forward one day
    IF( EPI_ETIME > ETIME ) THEN
        CALL NEXTIME( EDATE, ETIME, 24*10000 )
    END IF

    !.......  Convert start and end dates and times back to GMT
    CALL NEXTIME( SDATE, STIME, TZONE*10000 )
    CALL NEXTIME( EDATE, ETIME, TZONE*10000 )

    !.......  Find the total number of time steps
    NSTEPS = 1 + SECSDIFF( SDATE, STIME, EDATE, ETIME ) / 3600

    !.......  Configure the month(s) of modeling period
    CALL DAYMON( SDATE, SMONTH, SDAY )
    CALL DAYMON( EDATE, EMONTH, EDAY )
    SYEAR = INT( SDATE/1000 )
    EYEAR = INT( EDATE/1000 )

    IF( SYEAR < EYEAR ) EMONTH = EMONTH + 12

    !.......  Get number of lines in met list file
    NLINES = GETFLINE( TDEV, 'METLIST file' )

    !.......  Allocate memory for storing the met file information
    ALLOCATE( METLIST( NLINES ),&
              METDAYS( NSTEPS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'METDAYS', PROGNAME )

    METLIST = ' '        ! array
    METDAYS = 0

    !.......  Store lines of METLIST file
    CALL RDLINES( TDEV, 'METLIST file', NLINES, METLIST )

    MESG = 'Checking meteorology files...'
    CALL M3MSG2( MESG )

    METNAME = 'METFILE'

    !.......  Loop through all meteorology files
    DO N = 1, NLINES

        !.......  Close previous file if needed
        IF( FILEOPEN ) THEN
            IF( .NOT. CLOSE3( METNAME ) ) THEN
                MESG = 'Could not close meteorology file ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                FILEOPEN = .FALSE.
            END IF
        END IF

        !.......  Get physical file name for current iteration
        METFILE = METLIST( N )

        !.......  Set logical file name
        IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not set logical file name ' // 'for ' // TRIM( METFILE )
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Try to open file
        IF( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not open meteorology file ' // TRIM( METFILE )
            CALL M3MESG( MESG )
            CYCLE
        ELSE
            FILEOPEN = .TRUE.
        END IF

        !.......  Try to get description from file
        IF( .NOT. DESC3( METNAME ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not get description of meteorology file ' // TRIM( METFILE )
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Check that requested variables are in the file
        J = INDEX1( TVARNAME, NVARS3D, VNAME3D )
        IF( J <= 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not find ' // TRIM( TVARNAME ) //&
                   ' in file ' // TRIM( METFILE )
            CALL M3MESG( MESG )
        END IF

        J = INDEX1( PRESNAME, NVARS3D, VNAME3D )
        IF( J <= 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not find ' // TRIM( PRESNAME ) //&
                   ' in file ' // TRIM( METFILE )
            CALL M3MESG( MESG )
        END IF

        J = INDEX1( MIXNAME, NVARS3D, VNAME3D )
        IF( J <= 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not find ' // TRIM( MIXNAME ) //&
                   ' in file ' // TRIM( METFILE )
            CALL M3MESG( MESG )
        END IF

        IF( EFLAG ) CYCLE

        !.......  Initialize (or check) reference grid with meteorology data
        CALL CHKGRID( METNAME, 'GRID', 0, GRID_ERR )
        METNGRID = NGRID

        !....... If the dimensions were in error, print message and cycle
        IF( GRID_ERR ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Grid in meteorology file ' // TRIM( METFILE ) //     &
                   ' is inconsistent with previous files.'
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Build METDAYS array linking met files to a specific date and time

        !....... Find starting position in METDAYS array
        POS = 1 + SECSDIFF( SDATE, STIME, SDATE3D, STIME3D ) / 3600

        !.......  Make sure the met data is within the episode
        IF( POS + MXREC3D - 1 <= 0 .OR. POS > NSTEPS ) CYCLE

        !.......  Fill in array for number of steps in current met file
        DO L = POS, POS + MXREC3D - 1
            IF( L <= 0 ) CYCLE
            IF( L > NSTEPS ) EXIT
            METDAYS( L ) = N
        END DO

    END DO

    !.......  Exit if there was a problem with the meteorology files
    IF( EFLAG ) THEN
        MESG = 'Problem checking meteorology files'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Check for missing meteorology data
    MESG = 'Checking for missing meteorology data...'
    CALL M3MSG2( MESG )

    !.......  Check that all days are covered
    T = 1

    !.......  Loop over all days in episode
    DO

        !.......  Make sure we're within episode bounds
        IF( T > NSTEPS ) EXIT

        !.......  Check all hours in current day, accounting for time zone spread
        DO I = 0, 23 + TSPREAD
            K = T + I

            !.......  Double check episode bounds
            IF( K > NSTEPS ) EXIT

            !.......  If no met data for current step, try to find data
            IF( METDAYS( K ) == 0 ) THEN

                !.......  No data available, exit with error
                MESG = 'ERROR: Missing meteorology data during '&
                       // 'requested episode.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        END DO           ! end loop over hours in day

        !.......  Skip to start of next day
        T = T + 24
    END DO       ! end loop over days in episode

    !.......  Get the info of counties
    ALLOCATE( CNTYSRC( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CNTYSRC', PROGNAME )
    CNTYSRC = ' '

    !.......  Assign inventory counties to new array
    CNTYSRC( : ) = MCREFSORT( :, 1 )

    !.......  Open and write out the header
    !           information for SMOKE and MOVES-ready output files
    !        CALL HDRMISS3

    FDESC3D = ''

    FDESC3D( 1 ) = 'Met4moves gridded min/max temperatures output file for SMOKE Movesrmg'
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // CVSW
    FDESC3D( 4 ) = '/GRDNAME/ ' // GDNAM3D // COORD
    FDESC3D( 5 ) = '#MCXREF : ' // TRIM( CMCXREF )
    FDESC3D( 6 ) = '#MFMREF : ' // TRIM( CMFMREF )
    WRITE( FDESC3D( 7 ), 94010 ) '/START DATE/ ', SDATE
    WRITE( FDESC3D( 8 ), 94010 ) '/END DATE/ ', EDATE

    !.......  Set header values that cannot be default
    SDATE3D = SDATE
    STIME3D = 00000
    TSTEP3D = 240000
    NVARS3D = 2

    J = 1
    VNAME3D( J ) = 'MINTEMP'
    UNITS3D( J ) = 'Fahrenheit'
    VDESC3D( J ) = 'Daily minimum temperature by cell'
    VTYPE3D( J ) = M3REAL

    J = 2
    VNAME3D( J ) = 'MAXTEMP'
    UNITS3D( J ) = 'Fahrenheit'
    VDESC3D( J ) = 'Daily maximum temperature by cell'
    VTYPE3D( J ) = M3REAL

    !.......  Open new file
    ONAME = 'SMOKE_OUTFILE'
    IF( .NOT. OPEN3( ONAME, FSUNKN3, PROGNAME ) ) THEN
        MESG = 'Could not create new output file ' // TRIM( ONAME )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Open output file for MOVES model
    ODEV1 = PROMPTFFILE(&
         'Enter logical name for MOVES-ready output files',&
         .FALSE., .TRUE., 'MOVES_OUTFILE', PROGNAME )
    WRITE( ODEV1,'(A)' )'#DESC MOVES-ready input file for RPD and RPV modes'
    WRITE( ODEV1,'(A)' )'#GRDNAME: ' // GDNAM3D // COORD
    WRITE( ODEV1,'(A)' )'#MCXREF : ' // TRIM( CMCXREF )
    WRITE( ODEV1,'(A)' )'#MFMREF : ' // TRIM( CMFMREF )
    WRITE( ODEV1,94010 )'#MODELING PERIOD : ', EPI_SDATE,'-',EPI_EDATE
    WRITE( ODEV1,94020 )'#TEMP_BUFFER_BIN : ', TEMPBIN
    WRITE( ODEV1,'(A,I5)' ) '#PP_TEMP_INCREMENT ' , PPTEMP
    WRITE( ODEV1,'(A)' )'RefCounty,FuelMonth,Temperature'//&
         'ProfileID,RefMinT,RefMaxT,Temp1,Temp2,,,,,,,,,,,,Temp24'

    !.......  Open output file
    ODEV2 = PROMPTFFILE(&
         'Enter logical name for MOVES-ready output files',&
         .FALSE., .TRUE., 'MOVES_RH_OUTFILE', PROGNAME )
    WRITE( ODEV2,'(A)' )'#DESC MOVES-ready Temperature-bin-specific'&
                     // ' averaged RH output for RPD and RPV modes'
    WRITE( ODEV2,'(A)' )'#GRDNAME: ' // GDNAM3D // COORD
    WRITE( ODEV2,'(A)' )'#MCXREF : ' // TRIM( CMCXREF )
    WRITE( ODEV2,'(A)' )'#MFMREF : ' // TRIM( CMFMREF )
    WRITE( ODEV2,94010 )'#MODELING PERIOD : ', EPI_SDATE,'-',EPI_EDATE
    WRITE( ODEV2,94020 )'#TEMP_BUFFER_BIN : ', TEMPBIN
    WRITE( ODEV2,'(A,I5)' ) '#PD_TEMP_INCREMENT ' , PDTEMP
    WRITE( ODEV2,'(A,I5)' ) '#PV_TEMP_INCREMENT ' , PVTEMP
    WRITE( ODEV2,'(A)' )'RefCounty,FuelMonth,avgRH,'&
                      //'min_temp,max_temp,tempBin'

    !.......  Estimate possible max no of temp bins
    !           Lowest ambient temperature ever measured (-128F)
    !           Highest ambient temperature ever measured (138F)
    MXTBIN = INT( 350 / PDTEMP )         ! Temp range from -150F to 200F
    NFUEL = 12
    RHTBIN  = 0.0
    NRHTBIN = 0

    !.......  Allocate met variable, dates/daylight saving,
    !.......  source met, meteorology profiles arrays
    ALLOCATE( TA( METNGRID ),&
              QV( METNGRID ),&
            PRES( METNGRID ),&
         MAXTDAY( METNGRID ),&
         MINTDAY( METNGRID ),&
        MAXTCELL( METNGRID ),&
        MINTCELL( METNGRID ),&
         DAYBEGT( NSRC ),&
         DAYENDT( NSRC ),&
         LDAYSAV( NSRC ),&
         MAXTSRC( NSRC ),&
         MINTSRC( NSRC ),&
          TKHOUR( NSRC,24 ),&
         NTKHOUR( NSRC,24 ),&
        FUELCNTY( NREFC,NFUEL ),&
          TKFUEL( NREFC,NFUEL,24 ),&
        MAXTFUEL( NREFC,NFUEL ),&
        MINTFUEL( NREFC,NFUEL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TA...MINTFUEL', PROGNAME )

    MAXTDAY  = BADVAL3
    MINTDAY  = -1.0*BADVAL3
    MAXTCELL = BADVAL3
    MINTCELL = -1.0*BADVAL3

    !.......  Create array of which sources are affected by daylight savings
    CALL GETDYSAV( NSRC, CNTYSRC, LDAYSAV )

    !.......  Source met variable arrays
    MAXTSRC =  BADVAL3
    MINTSRC = -BADVAL3

    !.......  Allocate memory for storing meteorology profiles
    TKHOUR  = 0.0
    NTKHOUR = 0
    FUELIDX  = 0
    FUELCNTY = 0
    TKFUEL   = 0.0
    TKREFHR  = 0.0
    MAXTFUEL = BADVAL3
    MINTFUEL = -1*BADVAL3


    !.......  Process meteorology data...
    MESG = 'Processing meteorology data using variables ' //&
           TRIM( TVARNAME ) // ', ' // TRIM( MIXNAME ) //&
           ', ' // TRIM( PRESNAME ) // '...'
    CALL M3MSG2( MESG )

    !.......  Loop through days/hours of meteorology files
    ODATE = SDATE
    JDATE = SDATE
    JTIME = STIME
    PDATE = -9
    OTIME = -9
    LDATE = -9

    !.......  Process max/min temperatures and avg RH
    !           fuelmonth

    !.......  loop over hours
    DO T = 1, NSTEPS

        !.......  When new day...
        IF ( JDATE /= LDATE ) THEN

        !.......  Write message for day of week and date
            DAY = WKDAY( JDATE )
            MESG = 'Processing ' // DAYS( DAY ) // ' ' // MMDDYY( JDATE )
            CALL M3MSG2( MESG )

        !.......  Set start and end hours of day for all sources
            CALL SETSRCDY( NSRC, JDATE, TZONES, LDAYSAV, .TRUE., DAYBEGT, DAYENDT )
        END IF

        !.......  Determine input file for this hour
        POS = T

        !.......  Get file number for current iteration
        FILENUM = METDAYS( POS )
        IF( FILENUM <= 0 ) THEN
            ALT_DATA = .TRUE.
        ELSE
            ALT_DATA = .FALSE.
        END IF

        !.......  Skip file opening when not doing day averaging and using alternate data
        IF( .NOT. ALT_DATA ) THEN

            !.......  Get file name
            METFILE = METLIST( ABS( FILENUM ) )

            !.......  Close previous file if needed
            IF( METFILE .NE. PREVFILE ) THEN
                IF( FILEOPEN ) THEN
                    IF( .NOT. CLOSE3( METNAME ) ) THEN
                        MESG = 'Could not close meteorology file ' // TRIM( PREVFILE )
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    ELSE
                        FILEOPEN = .FALSE.
                    END IF
                END IF

                PREVFILE = METFILE

            END IF

            !.......  Set logical file name
            IF( .NOT. SETENVVAR( METNAME, METFILE ) ) THEN
                MESG = 'Could not set logical file name for ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            !.......  Open the meteorology data file
            IF ( .NOT. OPEN3( METNAME, FSREAD3, PROGNAME ) ) THEN
                MESG = 'Could not open meteorology file ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            ELSE
                FILEOPEN = .TRUE.
            END IF

            !.......  Read current meteorology file
            IF ( .NOT. READ3( METNAME, TVARNAME, 1, JDATE, JTIME, TA ) ) THEN
                MESG = 'Could not read ' // TRIM( TVARNAME ) //&
                       ' from ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
            IF ( .NOT. READ3( METNAME, MIXNAME, 1, JDATE, JTIME, QV ) ) THEN
                MESG = 'Could not read ' // TRIM( MIXNAME ) //&
                       ' from ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            IF ( .NOT. READ3( METNAME, PRESNAME, 1, JDATE, JTIME, PRES ) ) THEN
                MESG = 'Could not read ' // TRIM( PRESNAME ) //&
                       ' from ' // TRIM( METFILE )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            !.......  Create hourly meteorology arrays by source
            CALL HOURMET( NSRC, CNTYSRC, TA, QV, PRES, JDATE, JTIME,&
                          DAYBEGT, LDAYSAV, PDTEMP, HFLAG )
        ELSE

            IF( JTIME == 230000 ) THEN
                MESG = 'WARNING: Missing meteorology file on '//&
                        MMDDYY( JDATE )
                CALL M3MSG2( MESG )
            END IF

        END IF          ! check for using alternate data or day averaging

        !.......  Make sure we've waited long enough to catch all time zones
        IF( POS > TSPREAD ) THEN

            !.......  If last day of month, process monthly averages
            CALL DAYMON( ODATE+1, TMPMNTH, DAY )
            CALL DAYMON( ODATE, MONTH, DAY )

            !.......  Processing daily SMOKE-ready output
            DO S = 1, NSRC

                !.......  compute local time
                HOURIDX = 1 + ( JTIME - DAYBEGT( S ) ) / 10000
                IF( HOURIDX <= 0 ) HOURIDX = HOURIDX + 24

                !.......  output when local time is 24hr. if not, skip
                IF( HOURIDX /= 24 ) CYCLE

                !.......  retreive inv/ref counties
                INVCOUNTY = MCREFSORT( S,1 )
                REFCOUNTY = MCREFSORT( S,2 )

                NS = FINDC( INVCOUNTY, NSRGFIPS, SRGFIPS )

                IF( NS < 1 ) CYCLE

                IF( .NOT. FIRSTIME ) THEN
                    PDATE = ODATE
                    FIRSTIME = .TRUE.
                END IF

                !.......  Output gridded min/max temp by cell
                IF( ODATE /= PDATE ) THEN

                    IF( .NOT. WRITE3( ONAME, 'MINTEMP', PDATE, 0, MINTCELL ) ) THEN
                        MESG = 'Could not write MINTEMP from ' // TRIM( ONAME )
                        CALL M3EXIT( PROGNAME, ODATE, 0, MESG, 2 )
                    END IF
                    IF( .NOT. WRITE3( ONAME, 'MAXTEMP', PDATE, 0, MAXTCELL ) ) THEN
                        MESG = 'Could not write MAXTEMP from ' // TRIM( ONAME )
                        CALL M3EXIT( PROGNAME, ODATE, 0, MESG, 2 )
                    END IF

                    !.......  initializing daily min/max temps
                    MINTCELL = -1.0*BADVAL3
                    MAXTCELL = BADVAL3

                END IF

                !.......  Update daily min/max temps by cell
                DO NC = 1, NCELLS( NS )

                    C = FIPCELL( NC,NS )

                    MAXTCELL( C ) = MAX( MAXTDAY( C ), MAXTCELL( C ) )
                    MINTCELL( C ) = MIN( MINTDAY( C ), MINTCELL( C ) )

                    !.......  initializing daily min/max temps
                    MINTDAY( C ) = -1.0*BADVAL3
                    MAXTDAY( C ) = BADVAL3

                END DO

                PDATE = ODATE

            END DO

            !.......  Estimate fuelmonth averaged monthly ref county temp and RH
            IF( TMPMNTH /= MONTH .OR. T > NSTEPS - 23 ) THEN
                CALL AVG_REF_COUNTY_TEMP( JDATE, JTIME )
            END IF

        END IF

        IF( POS > TSPREAD ) CALL NEXTIME( ODATE, OTIME, 10000 )

    !.......  Increment loop time
        LDATE = JDATE
        CALL NEXTIME( JDATE, JTIME, 10000 )

    END DO       !  End loop on hours of temperature files

    !.......  Output last date gridded min/max temp by cell
    IF( .NOT. WRITE3( ONAME, 'MINTEMP', ODATE, 0, MINTCELL ) ) THEN
        MESG = 'Could not write MINTEMP from ' // TRIM( ONAME )
        CALL M3EXIT( PROGNAME, ODATE, 0, MESG, 2 )
    END IF
    IF( .NOT. WRITE3( ONAME, 'MAXTEMP', ODATE, 0, MAXTCELL ) ) THEN
        MESG = 'Could not write MAXTEMP from ' // TRIM( ONAME )
        CALL M3EXIT( PROGNAME, ODATE, 0, MESG, 2 )
    END IF

    !.......  Compute ref county min/max temperatures per fuel month
    DO NR = 1, NREFC

        REFCOUNTY = MCREFIDX( NR,1 )

        !.......  Choose month-specific fulemonth county
        L = FINDCFIRST( REFCOUNTY, NREFF,  FMREFSORT( :,1 ) )
        K = FINDCFIRST( REFCOUNTY, NFUELC, FMREFLIST( :,1 ) )

        IF( L < 0 .OR. K < 0 ) THEN
            MESG = 'ERROR: Fuel month input file MUST contain '&
                   // 'reference county: ' //  REFCOUNTY
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        NMON = STR2INT( FMREFLIST( K, 2 ) )          ! no month available per ref county

        !.......  Loop over months per ref. county
        N = 0
        PRVCMONTH = 0       ! prv tmp processing month
        PRVFMONTH = 0       ! prv tmp fuelmonth

        DO J = L, L + NMON - 1

            FUELMONTH = STR2INT( FMREFSORT( J,2 ) )        ! processing fuelmonth/county
            CURMONTH  = STR2INT( FMREFSORT( J,3 ) )       ! processing current month per ref. county

            IF( CURMONTH <  SMONTH ) CYCLE
            IF( CURMONTH >  EMONTH ) CYCLE

            !.......  Store fuelmonth specific values into arrays
            IF( FUELMONTH /= PRVFMONTH ) THEN
                IF( N > 0 ) THEN

                    TKREFHR = TKREFHR / N
                    MAXTEMP = MAXTEMP + TEMPBIN
                    MINTEMP = MINTEMP - TEMPBIN

                    CALL WRTEMPROF( ODEV1, ODEV2, SYEAR,&
                         REFCOUNTY, PRVFMONTH, PDTEMP, PPTEMP,&
                         TKREFHR, MAXTEMP, MINTEMP, TEMPBIN, MINNORH )

                END IF

                !.......  initialize local variables
                N     = 0
                TKREFHR = 0.0
                MAXTEMP = BADVAL3
                MINTEMP = -1*BADVAL3

            END IF

            N = N + 1
            NF = FUELMONTH
            MAXTEMP = MAX( MAXTEMP, MAXTFUEL( NR,NF ) )
            MINTEMP = MIN( MINTEMP, MINTFUEL( NR,NF ) )

            DO TT = 1,24
                TKREFHR( TT ) = TKREFHR( TT ) + TKFUEL( NR,NF,TT )&
                                / FUELCNTY( NR,NF )
            END DO

            PRVCMONTH = CURMONTH
            PRVFMONTH = FUELMONTH

        END DO

        !.......  Store last fuelmonth specific values into arrays
        IF( N == 0 ) CYCLE
        TKREFHR = TKREFHR / N
        MAXTEMP = MAXTEMP + TEMPBIN
        MINTEMP = MINTEMP - TEMPBIN

        CALL WRTEMPROF( ODEV1, ODEV2, SYEAR, REFCOUNTY,&
             PRVFMONTH, PDTEMP, PPTEMP, TKREFHR, MAXTEMP,&
             MINTEMP, TEMPBIN, MINNORH )

    END DO       ! end of loop of reference couties

    !.......  Exit if there was a problem with the meteorology files
    IF( EFLAG ) THEN
        MESG = 'Problem processing meteorology files'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !....... End program successfully

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93010 FORMAT( I8, 1X, F13.5, 1X, F13.5, 1X, I8, 1X, A )

    !.......   Internal buffering fosrmats...... 94xxx

94000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( A, F8.2  )

94030 FORMAT( A, I5 )

94050 FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,&
            A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

94060 FORMAT( I6.6, I5, 3X, I5, I10, 3F10.2 )

94070 FORMAT( A, F5.1, A )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram estimates ref. county level
    !               min/max Temperatures over fuelmonth

    SUBROUTINE AVG_REF_COUNTY_TEMP( JDATE, JTIME )

        !.......  local argument
        INTEGER, INTENT( IN ) :: JDATE            ! met data date
        INTEGER, INTENT( IN ) :: JTIME            ! met data time

        !----------------------------------------------------------------------
        !.......  Loop over sources
        DO S = 1, NSRC

            !.......  compute local time
            HOURIDX = 1 + (JTIME - DAYBEGT( S ) ) / 10000
            IF( HOURIDX <= 0 ) HOURIDX = HOURIDX + 24

            !.......  output when local time is 24hr. if not, skip
            IF( HOURIDX /= 24 ) CYCLE

            !.......  retrieve inv/ref counties
            INVCOUNTY = MCREFSORT( S,1 )
            REFCOUNTY = MCREFSORT( S,2 )

            !.......  averaging Temp per inventory county
            DO TT = 1,24

                !.......  Skip sources with no days; this can happen when the
                !                       gridding surrogates do not contain data for all counties
                TKHOUR( S,TT ) = TKHOUR( S,TT ) / NTKHOUR( S,TT )

            END DO

            MAXTEMP = MAXTSRC( S )
            MINTEMP = MINTSRC( S )

            L = FINDCFIRST( REFCOUNTY, NREFF, FMREFSORT( :,1 ) )
            K = FINDCFIRST( REFCOUNTY, NFUELC,FMREFLIST( :,1 ) )
            NMON = STR2INT( FMREFLIST( K, 2 ) )               ! no month of ref county

            !.......  Loop over months per ref. county
            FUELMONTH = 0
            DO J = L, L + NMON - 1
                CURMONTH  = STR2INT( FMREFSORT( J,3 ) )                ! processing  current month per ref. county
                IF( CURMONTH == MONTH ) FUELMONTH = STR2INT( FMREFSORT( J,2 ) )              ! processing fuelmonth/county
            END DO

            IF( FUELMONTH == 0 ) THEN
                WRITE( MESG,94010 ) 'ERROR: Could not find fuel month "',&
                    MONTH, ' for reference county ' // REFCOUNTY
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

            !.......  Calculation monthly fuel month max/min temp and avg RH
            !                   per ref. county
            NR = FINDC( REFCOUNTY,NREFC, MCREFIDX( :,1 ) )
            NF = FUELMONTH

            FUELCNTY( NR,NF ) = FUELCNTY( NR,NF ) + 1
            TKFUEL( NR,NF,: ) = TKFUEL( NR,NF,: ) + TKHOUR( S,: )
            MAXTFUEL( NR,NF ) = MAX( MAXTFUEL( NR,NF ), MAXTEMP )
            MINTFUEL( NR,NF ) = MIN( MINTFUEL( NR,NF ), MINTEMP )

            !.......  reinitializing local arrays for next month averaging
            TKHOUR ( S,: ) = 0.0
            NTKHOUR( S,: ) = 0
            MAXTSRC( S ) = BADVAL3
            MINTSRC( S ) = -1.0*BADVAL3

        END DO

        !******************  FORMAT  STATEMENTS   ******************************
        !.......   Internal buffering formats...... 94xxx
94010   FORMAT( 10( A, :, I8, :, 1X ) )

94060   FORMAT( I6.6, I5, 3X, I5, I10, 3F10.2 )

    END SUBROUTINE AVG_REF_COUNTY_TEMP

END PROGRAM MET4MOVES
