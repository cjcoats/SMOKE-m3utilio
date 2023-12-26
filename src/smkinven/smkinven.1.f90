
PROGRAM SMKINVEN

!***********************************************************************
!  program body starts at line 171
!
!  DESCRIPTION:
!    The smkinven program reads the any source inventory in one of four
!    formats: EMS-95, EPS, IDA, and SMOKE list format.  It permits a flexible
!    definition of a point source, which depends on the inventory input
!    formats.  It allows any number of inventory pollutants, within the limit
!    of I/O API (this works out to only 15 pollutants per file).
!
!  PRECONDITIONS REQUIRED:
!    Set environment variables:
!      PROMPTFLAG:    If N, default inputs are used
!      RAW_DUP_CHECK: If Y, duplicate sources disallowed
!      VELOC_RECALC:  If Y, recalculates velocity based on flow
!      WEST_HSPHERE:  If N, does not reset positive longitudes to negative
!    Input files:
!      PTINV: ASCII point sources inventory
!      PSTK: Replacement stack parameters file
!      ZONES: Time zones files
!      SIPOLS: Master list of pollutant codes and names (in output order)
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!    Subroutines: I/O API subroutines, INITEM, CHECKMEM, RDTZONE, RDSIPOLS,
!       RDPTINV, FMTCSRC, FIXSTK, WRPTSCC, WRPTREF, OPENPNTS, WPNTSCHR,
!       WPNTSPOL
!    Functions: I/O API functions, GETFLINE, GETTZONE
!
!  REVISION  HISTORY:
!    started 10/98 by M Houyoux as rawpoint.f from emspoint.F 4.3
!    smkinven changes started 4/98
!    toxics changes 11/2002  A. Holland
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: CIFIP, CISIC, CSRCTYP, TZONES, CSCC, IDIU, IWEK,&
    &                    CINTGR, CEXTORL
!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MXIDAT, INVSTAT, INVDNAM, FIREFLAG, INVDCOD,&
    &                    FF10FLAG, MEDSFLAG, NCDFLAG, APIFLAG, FIREFF10
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NIPOL, NIACT, NIPPA, EIIDX, INV_MON,&
    &                   EINAM, AVIDX, ACTVTY, EANAM, NSRC

!.........  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: DAYINVFLAG, HRLINVFLAG

    IMPLICIT NONE

!...........   INCLUDES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETTZONE
    LOGICAL, EXTERNAL :: USEEXPGEO

!...........  LOCAL PARAMETERS and their descriptions:

    CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!.........  LOCAL VARIABLES and their descriptions:

!.........  Day-specific and hour-specific variable indices
    INTEGER     ::  DEAIDX( MXVARS3 ) = 0
    INTEGER     ::  DSPIDX( MXSPDAT ) = 0
    INTEGER     ::  HEAIDX( MXVARS3 ) = 0
    INTEGER     ::  HSPIDX( MXSPDAT ) = 0

!.........  Array that contains the names of the inventory variables needed for
!           this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

!.........  File units and logical/physical names

    INTEGER    :: ADEV = 0  !  unit no. for REPINVEN file
    INTEGER    :: CDEV = 0  !  unit no. for SCCs description
    INTEGER    :: DDEV = 0  !  unit no. for day-specific input file
    INTEGER    :: HDEV = 0  !  unit no. for hour-specific input file
    INTEGER    :: IDEV = 0  !  unit no. for inventory file (various formats)
    INTEGER    :: LDEV = 0  !  unit no. for log file
    INTEGER    :: MDEV = 0  !  unit no. for map inventory file
    INTEGER    :: ODEV = 0  !  unit number for ORIS description
    INTEGER    :: PDEV = 0  !  unit number for inventory data table
    INTEGER    :: RDEV = 0  !  unit no. for def stack pars or mobile codes
    INTEGER    :: UDEV = 0  !  unit no. for non-HAP inclusions/exclusions file
    INTEGER    :: SDEV = 0  !  unit no. for ASCII output inventory file
    INTEGER    :: TDEV = 0  !  unit no. for MEDS daily/hourly inventory file
    INTEGER    :: YDEV = 0  !  unit no. for area-to-point factors file
    INTEGER    :: ZDEV = 0  !  unit no. for time zone file

    CHARACTER(NAMLEN3) :: ANAME = ' '! inven ASCII output logical name
    CHARACTER(NAMLEN3) :: DNAME = ' '! day-specific input logical name
    CHARACTER(NAMLEN3) :: ENAME = ' '! inven I/O API output logical name
    CHARACTER(NAMLEN3) :: GNAME = ' '! gridded I/O API input logical
    CHARACTER(NAMLEN3) :: HNAME = ' '! hour-specific input logical name
    CHARACTER(NAMLEN3) :: INAME = ' '! inven input logical name
    CHARACTER(NAMLEN3) :: TNAME = ' '! tmp day- or hour-specific input logical name

!...........   Other local variables

    INTEGER         S, I, J, J1, J2, J3, K, L, L2, V !  counters and indices

    INTEGER      :: DNSTEP = 0 ! day-specific data time step number
    INTEGER      :: DSDATE = 0 ! day-specific data start date
    INTEGER      :: DSTIME = 0 ! day-specific data start time

    INTEGER      :: HNSTEP = 0 ! day-specific data time step number
    INTEGER      :: HSDATE = 0 ! day-specific data start date
    INTEGER      :: HSTIME = 0 ! day-specific data start time
    INTEGER         INSTEP     ! expected input time step HHMMSS
    INTEGER         IOS        ! I/O status
    INTEGER         MAXK       ! test for maximum value of K in output loop
    INTEGER      :: MXSRCDY= 0 ! max no. day-specific sources
    INTEGER      :: MXSRCHR= 0 ! max no. hour-specific sources
    INTEGER      :: NDAT = 0   ! tmp no. actual pols & activities
    INTEGER      :: NINVARR = 0! no. inventory variables to read
    INTEGER      :: NRAWBP = 0 ! number of sources with pollutants
    INTEGER      :: NRAWSRCS= 0! number of unique sources
    INTEGER      :: NVARDY = 0 ! no. day-specific variables
    INTEGER      :: NVSPDY = 0 ! no. day-specific special variables
    INTEGER      :: NVARHR = 0 ! no. hour-specific variables
    INTEGER      :: NVSPHR = 0 ! no. hour-specific special variables
    INTEGER         OUTSTEP    ! output time step HHMMSS for day/hour data
    INTEGER         TZONE      ! output time zone for day- & hour-specific

    LOGICAL         A2PFLAG          ! true: using area-to-point processing
    LOGICAL      :: CFLAG = .FALSE.  ! true: CEM processing
    LOGICAL         IFLAG            ! true: average inventory inputs used
    LOGICAL         NONPOINT         ! true: importing nonpoint inventory
    LOGICAL         ORLFLG           ! true: ORL format inventory
    LOGICAL         STKFLG           ! true: check stack parameters

    CHARACTER(5)          TYPNAM      !  'day' or 'hour' for import
    CHARACTER(256)        MESG        !  message buffer
    CHARACTER(NAMLEN3) :: GRDNM = ' ' !  I/O API input file grid name
    CHARACTER(PHYLEN3) :: VARPATH = './' ! path for pol/act files
    CHARACTER(FIPLEN3)    CFIP        ! temporary FIPS code
    CHARACTER(FIPLEN3)    PCFIP       ! previous FIPS code

    CHARACTER(16) :: PROGNAME = 'SMKINVEN'   !  program name

!***********************************************************************
!   begin body of program SMKINVEN

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Output time zone
    TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTZONE"', 2 )
    END IF

!.........  Get names of input files
    CALL OPENINVIN( CATEGORY, IDEV, DDEV, HDEV, RDEV, SDEV,&
    &                PDEV, ZDEV, CDEV, ODEV, UDEV, YDEV,&
    &                ENAME, INAME, DNAME, HNAME )

!.........  Set controller flags depending on unit numbers
    IFLAG   = ( IDEV .NE. 0 )
    A2PFLAG = ( YDEV .NE. 0 )

!.........  Set IOAPI or NetCDF gridded input file name, if available
    IF( APIFLAG .OR. NCDFLAG ) GNAME = ENAME

    MESG = 'Setting up to read inventory data...'
    CALL M3MSG2( MESG )

!.........  Read country, state, and county file for time zones
    IF( USEEXPGEO() ) THEN
        CALL RDGEOCODES( 1, I )
    ELSE IF( ZDEV .GT. 0 ) THEN
        CALL RDSTCY( ZDEV, 1, I )   !  "I" used as a dummy
    END IF

!.........  Read, sort, and store inventory data table file
    CALL RDCODNAM( PDEV )

!.........  Read/store informatino for MEDS inv processing
    IF( MEDSFLAG ) CALL RDMEDSINFO  ! read GAI_LOOKUP_TABLE (col/row to lat/lon) for MEDS format inv

!.........  Define the month of inventory processing
    MESG = 'Define the processing inventory month'
    INV_MON = ENVINT( 'SMKINVEN_MONTH', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMKINVEN_MONTH"', 2 )
    END IF

!.........  Process for ASCII average day or annual inventory
    IF( IFLAG ) THEN

!.............  Read the source information from the raw inventory files,
!               store in unsorted order, and determine source IDs
!.............  The arrays that are populated by this subroutine call
!               are contained in the module MODSOURC
        IF( INV_MON == 0 ) THEN
            MESG = 'Processing Annual inventory....'
            CALL M3MSG2( MESG )
        ELSE
            MESG = 'Processing ' // TRIM( MON_NAME( INV_MON ) ) //&
            &       ' inventory'
            CALL M3MSG2( MESG )
        END IF

        CALL M3MSG2( 'Reading inventory sources...' )

        CALL RDINVSRCS( IDEV, INAME, NRAWBP, NRAWSRCS, ORLFLG )

!.............  Process source information and store in sorted order
        CALL M3MSG2( 'Processing inventory sources...' )

        CALL PROCINVSRCS( NRAWSRCS )

!.............  Read the data from the raw inventory files and store in
!               sorted order
        CALL M3MSG2( 'Reading inventory data...' )

        CALL RDINVDATA( IDEV, INAME, NRAWBP, NONPOINT )

!.............  Check if ORL inventory and reset NHAPEXCLUDE if needed
        IF( .NOT. ORLFLG .AND. UDEV > 0 ) THEN
            MESG = 'NOTE: Ignoring SMK_PROCESS_HAPS setting ' //&
            &       'since an ORL or FF10 inventory is not being processed'
            CALL M3MSG2( MESG )
            UDEV = 0
        END IF

!.............  Process inventory records and store in sorted order
        CALL M3MSG2( 'Processing inventory data...' )

        CALL PROCINVEN( NRAWBP, UDEV, YDEV, CDEV, LDEV )

!.............  Integrate criteria and toxic pollutants
        IF( ORLFLG ) THEN
            CALL SETNONHAP( NRAWBP )
        END IF

!.............  Determine memory needed for actual pollutants list and actual
!               activities list and allocate them. Invstat has been updated
!               to be +/- 2 depending on whether the pollutant or activity was
!               present in the inventory.
        NIPOL = 0
        NIACT = 0
        DO I = 1, MXIDAT
            IF( INVSTAT( I ) .GT.  1 ) NIPOL = NIPOL + 1
            IF( INVSTAT( I ) .LT. -1 ) NIACT = NIACT + 1
        ENDDO

        NIPPA = NIPOL + NIACT

        ALLOCATE( EIIDX( NIPOL ),&
        &          EINAM( NIPOL ),&
        &          AVIDX( NIACT ),&
        &         ACTVTY( NIACT ),&
        &          EANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EIIDX...EANAM', PROGNAME )

!.............  Create list of actual pollutants and activities and indexes to
!               the master list. The order in EINAM and ACTVTY will be the
!               output order. The indexes are for accessing INVDCOD, if needed.
!.............  These are for opening output file and processing output data
        J1 = 0
        J2 = 0
        J3 = 0
        DO I = 1, MXIDAT

            IF( INVSTAT( I ) .GT. 0 ) THEN
                J1 = J1 + 1
                J3 = J3 + 1
                EIIDX( J3 ) = I
                EINAM( J3 ) = INVDNAM( I )
                EANAM( J1 ) = INVDNAM( I )
            END IF

            IF( INVSTAT( I ) .LT. 0 ) THEN
                J1 = J1 + 1
                J2 = J2 + 1
                AVIDX ( J2 ) = I
                ACTVTY( J2 ) = INVDNAM( I )
                EANAM ( J1 ) = INVDNAM( I )
            END IF

        END DO

!.............   Fix stack parameters for point sources
!.............   Some of these arguments are variables that are defined in the
!                module MODSOURC
        IF( CATEGORY .EQ. 'POINT' ) THEN
            MESG = 'Check stack parameter values'
            STKFLG = ENVYN( 'CHECK_STACKS_YN', MESG, .TRUE., IOS )
            IF( FIREFLAG .OR. FIREFF10 ) STKFLG = .FALSE.    ! skip check stack para when wildfire
            IF( STKFLG ) CALL FIXSTK( RDEV, NSRC )
        END IF

!.............  Set time zones based on country/state/county code. Note that a
!               few counties in the Western U.S. are divided by a time zone, so
!               this is not perfectly accurate for all counties.
        ALLOCATE( TZONES( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TZONES', PROGNAME )

        PCFIP = ''
        DO S = 1, NSRC
            CFIP = CIFIP( S )

            IF( CFIP /= PCFIP ) THEN
                TZONES( S ) = GETTZONE( CFIP )
                PCFIP = CFIP
            ELSE
                TZONES( S ) = TZONES( S - 1 )
            END IF
        END DO

!.............  Write out primary inventory files. Do this before the day- or
!               hour-specific processing so that if there is a problem, the
!               lengthy inventory import does not need to be redone...

!.............  Write out SCC file
        CALL WRCHRSCC( CSCC )

    END IF  ! For ASCII annual/ave-day inputs

!.........  Input gridded I/O API inventory data
    IF( APIFLAG .OR. NCDFLAG ) THEN

        IF( APIFLAG ) CALL RDGRDAPI( GNAME, GRDNM )
        IF( NCDFLAG ) CALL RDGRDNCF( INAME, ENAME )

!.............  initialize arrays for later
        ALLOCATE( CISIC ( NSRC ),&
        &         CSRCTYP( NSRC ),&
        &          CINTGR( NSRC ),&
        &         CEXTORL( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CISIC...CEXTORL', PROGNAME )

        CISIC   = ' '       ! array
        CSRCTYP = ' '       ! array
        CINTGR  = ' '       ! array
        CEXTORL = ' '       ! array

        IFLAG = .TRUE.

    END IF  ! For gridded I/O API NetCDF inventory

!.........  Output SMOKE inventory files
    IF( IFLAG ) THEN

!.............  Generate message to use just before writing out inventory files
!.............  Open output I/O API and ASCII files

        CALL OPENINVOUT( A2PFLAG, GRDNM, ENAME, ANAME, MDEV, SDEV,&
        &                 ADEV, VARPATH )

        MESG = 'Writing SMOKE ' // TRIM( CATEGORY ) //&
        &       ' SOURCE INVENTORY file...'

        CALL M3MSG2( MESG )

!.............  Write source characteristics to inventory files (I/O API and
!               ASCII)
        CALL WRINVCHR( ENAME, SDEV, A2PFLAG, NONPOINT )

!.............  Deallocate sorted inventory info arrays, except CSOURC
        CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .FALSE., 1, 1, 1 )

!.............  Write out average inventory data values
!.............  Compute inventory data values, if needed
        CALL WRINVEMIS( MDEV, VARPATH )

!.............  Deallocate sorted inventory info arrays
        CALL SRCMEM( CATEGORY, 'SORTED', .FALSE., .TRUE., 1, 1, 1 )

!.........  If the inventory is not being created, then read necessary
!           information from existing inventory files, which will be used
!           for day- and hour-specific data import.
    ELSE

!.............  Store source-category-specific header information,
!           including the inventory pollutants in the file (if any).  Note that
!           the I/O API head info is passed by include file and the
!           results are stored in module MODINFO.
        CALL GETSINFO( ENAME )

!.............  Since not reading inventory files yet, but need to know
!               if this is a fire inventory, scan inventory data names
!               to see if HFLUX is present. If so, FIREFLAG = .true.
        IF ( INDEX1( 'HFLUX',NIPPA,EANAM ) .GT. 0 ) FIREFLAG= .TRUE.

        NINVARR = 5
        IVARNAMS( 1 ) = 'CIFIP'   ! In case CEM input
        IVARNAMS( 2 ) = 'CSOURC'  ! In case non-CEM input
        IVARNAMS( 3 ) = 'CSCC'    ! In case CEM input (for reporting)
        IVARNAMS( 4 ) = 'CPDESC'  ! In case CEM input
        IVARNAMS( 5 ) = 'CINTGR'  ! for VOC + HAPs integration

        IF ( .NOT. FIREFLAG ) THEN
            NINVARR = 7
            IVARNAMS( 6 ) = 'CORIS'   ! In case CEM input
            IVARNAMS( 7 ) = 'CBLRID'  ! In case CEM input
        END IF

        CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR,&
        &               IVARNAMS )

    END IF !   End processing of average annual import or not

!.........  Read in daily or hourly MEDS emission values and output to a SMOKE inter output file
    IF( MEDSFLAG ) THEN
!.............  Read and output day-specific data
        IF( DAYINVFLAG ) THEN
            TYPNAM = 'day'
            TDEV   = DDEV
            CALL GENMEDSOUT( TDEV, TNAME, TZONE, TYPNAM )

        ELSE IF( HRLINVFLAG ) THEN
            TYPNAM = 'hour'
            TDEV   = HDEV
            CALL GENMEDSOUT( TDEV, TNAME, TZONE, TYPNAM )

        END IF

    END IF

!.........  Read in daily emission values and output to a SMOKE file
    IF( DAYINVFLAG .AND. .NOT. MEDSFLAG ) THEN

        INSTEP  = 240000
        OUTSTEP = 10000
        TYPNAM  = 'day'

!.............  Preprocess day-specific file(s) to determine memory needs.
!               Also determine maximum and minimum dates for output file.
        CALL GETPDINFO( DDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, DNAME,&
        &                DSDATE, DSTIME, DNSTEP, NVARDY, NVSPDY,&
        &                MXSRCDY, DEAIDX, DSPIDX, CFLAG )

!.............  Read and output day-specific data
        CALL GENPDOUT( DDEV, CDEV, ODEV, ADEV, TZONE, DSDATE,DSTIME,&
        &               DNSTEP, INSTEP, OUTSTEP, NVARDY, NVSPDY,&
        &               MXSRCDY, TYPNAM, DNAME, DEAIDX, DSPIDX, CFLAG )

    END IF

!.........  Read in hourly emission values and output to a SMOKE file
    IF( HRLINVFLAG .AND. .NOT. MEDSFLAG ) THEN

        INSTEP  = 10000
        OUTSTEP = 10000
        TYPNAM  = 'hour'

!.............  Preprocess hour-specific file(s) to determine memory needs.
!               Also determine maximum and minimum dates for output file.
        CALL GETPDINFO( HDEV, TZONE, INSTEP, OUTSTEP, TYPNAM, HNAME,&
        &                HSDATE, HSTIME, HNSTEP, NVARHR, NVSPHR,&
        &                MXSRCHR, HEAIDX, HSPIDX, CFLAG )

!.............  Read and output hour-specific data
        CALL GENPDOUT( HDEV, CDEV, ODEV, ADEV,TZONE, HSDATE,HSTIME,&
        &               HNSTEP, INSTEP, OUTSTEP, NVARHR, NVSPHR,&
        &               MXSRCHR, TYPNAM, HNAME, HEAIDX, HSPIDX, CFLAG )

    END IF

!.............  Write inventory report file
!        CALL M3MSG2( 'Writing inventory report file...' )
!        IF( .NOT. MEDSFLAG ) CALL WREPINVEN( ADEV, CDEV )

!.........  End program successfully
    CALL M3EXIT( PROGNAME, 0, 0, 'Completion', 0 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Informational (LOG) message formats... 92xxx

92000 FORMAT( 5X, A )

92010 FORMAT( 5X, A, :, I10 )


!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!93041   FORMAT( I5, X, I5, X, I3, X, I8, X, I5.5, 3( X, I3 ) )

93060 FORMAT( 10( A, :, E10.3, :, 1X ) )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94020 FORMAT( A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,&
    &        A, I6, 1X, A, I6, 1X, A, :, I6 )

94040 FORMAT( A, I2.2 )

94060 FORMAT( 10( A, :, E10.3, :, 1X ) )

94080 FORMAT( '************  ', A, I7, ' ,  ' , A, I12 )

END PROGRAM SMKINVEN
