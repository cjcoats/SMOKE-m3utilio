
PROGRAM TEMPORAL

!***********************************************************************
!
!  DESCRIPTION:
!    This program computes the hourly emissions data from inventoryg emissions
!    and/or activity and emission factor data. It can read average-inventory,
!    day-specific and hour-specific emissions and activity data.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       copied by: M. Houyoux 01/99
!       origin: tmppoint.F 4.3
!
!       Created by Marc Houyoux
!
!       Revised 07/2014 by C.Coats for  new GENTPRO CSV profiles and cross-references
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
!*************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: TZONES, TPFLAG, FLTRDAYL

!.........  This module contains the temporal profile tables
    USE MODTMPRL, ONLY: NMON, NWEK, NHRL, STDATE, STTIME, RUNLEN,&
    &                    ITDATE, METPROFLAG, METPROTYPE, IPOL2D, LTFLAG

!.........  This module contains emission factor tables and related
    USE MODEMFAC, ONLY: EMTPOL, NEPOL, TEMPEF, USETIME

!.........  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: DYPNAM, DYPDSC, NDYPOA, NDYSRC,&
    &                    HRPNAM, HRPDSC, NHRPOA, NHRSRC,&
    &                    LDSPOA, LHSPOA, LHPROF,&
    &                    INDXD, EMACD, INDXH, EMACH,&
    &                    EMAC, EMACV, EMIST, EMFAC, TMAT

!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP, MXIDAT, INVDNAM, INVDVTS

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, BYEAR, NIPPA, EANAM, NSRC,&
    &                   INVPIDX, NIPOL, EAREAD, EINAM, ACTVTY

!.........  This module is used for MOBILE setup information
    USE MODMBSET, ONLY: DAILY, WEEKLY, MONTHLY, EPISLEN

    IMPLICIT NONE

!.........  INCLUDES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!..........  EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: RDTPROF
    LOGICAL, EXTERNAL :: USEEXPGEO

!.........  LOCAL PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'TEMPORAL' ! program name
    CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$'  ! CVS revision tag

!.........  Array that contains the names of the inventory variables needed for
!           this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

!.........  Day-specific, hour-specific data, and elevated sources data.
!.........  These need only to allow enough dimensions for one read per
!           pollutant per time step.

    INTEGER                 NPELV        ! optional elevated source-count
    INTEGER, ALLOCATABLE :: INDXE( : )   ! SMOKE source IDs
    REAL   , ALLOCATABLE :: EMISE( : )   ! elevated source emissions

!.........  Names of pollutants and activities associated with output variables
    CHARACTER(NAMLEN3), ALLOCATABLE:: ALLIN( : )

!.........  Reshaped input variables and output variables
    INTEGER         NGRP                ! no. of pol/emis-types groups
    INTEGER         NGSZ                ! no. of pols/emis-types per group
    CHARACTER(NAMLEN3), ALLOCATABLE:: ALLIN2D( :,: )
    CHARACTER(NAMLEN3), ALLOCATABLE:: EANAM2D( :,: )
    CHARACTER(NAMLEN3), ALLOCATABLE:: EAREAD2D( : )

!...........   Logical names and unit numbers

    INTEGER      :: CDEV = 0!  unit number for region codes file
    INTEGER      :: HDEV = 0!  unit number for holidays file
    INTEGER         KDEV    !  unit number for time periods file
    INTEGER         LDEV    !  unit number for log file
    INTEGER         PDEV    !  unit number for supplemental tmprl file
    INTEGER         SDEV    !  unit number for ASCII inventory file

    CHARACTER(16) :: ANAME = ' '    !  logical name for ASCII inven input
    CHARACTER(16) :: DNAME = 'NONE' !  day-specific  input file, or "NONE"
    CHARACTER(16) :: ENAME = ' '    !  logical name for I/O API inven input
    CHARACTER(16) :: FNAME = ' '    !  emission factors file
    CHARACTER(16) :: HNAME = 'NONE' !  hour-specific input file, or "NONE"
    CHARACTER(16) :: PNAME = 'NONE' !  hour-specific profile input file, or "NONE"
    CHARACTER(16) :: TNAME = ' '    !  timestepped (low-level) output file

!...........   Other local variables

    INTEGER         I, II, IS, J, K, L, L1, L2, N, S, T

    INTEGER         IOS, IOS1, IOS2, IOS3, IOS4 ! i/o status
    INTEGER         IOS6, IOS7, IOS8, IOS9      ! i/o status
    INTEGER         AVERTYPE            ! time period averaging type
    INTEGER         DYSTPOS, DYENDPOS   ! start and end position in file name string
    INTEGER         EARLYDATE           ! earliest starting date based on time zones
    INTEGER         EARLYTIME           ! earliest starting time based on time zones
    INTEGER         EARLST              ! earliest starting date within entire episode periods
    INTEGER         EDATE, ETIME        ! ending Julian date and time
    INTEGER         EFSDATE, EFEDATE    ! start and end date of current ef file
    INTEGER         ENLEN               ! length of ENAME string
    INTEGER         ENDPOS              ! ending position in ef day array
    INTEGER         FIRSTPOS            ! temporary position in file name string
    INTEGER         FDATE, FTIME        ! emission factor date and time
    INTEGER         HYPPOS              ! position of hyphen in file name string
    INTEGER         JDATE, JTIME        ! Julian date and time
    INTEGER         NDAYS               ! no. days in episode
    INTEGER         JYEAR, JMNTH, JDAYS ! tmp year, month, and date
    INTEGER         LATEDATE            ! latest ending date based on time zones
    INTEGER         LATETIME            ! latest ending time
    INTEGER         LATEST              ! latest starting date within entire episode periods
    INTEGER         NINVARR             ! no. inventory variables to read
    INTEGER         NLINES              ! no. lines in ef list file
    INTEGER         NMATX               ! size of ungridding matrix
    INTEGER         NMAJOR              ! no. major sources
    INTEGER         NPING               ! no. ping sources
    INTEGER         NSTEPS              ! number of output time steps
    INTEGER         NTPERIOD            ! No of time periods
    INTEGER      :: PYEAR = 0           ! projected year
    INTEGER         SDATE, STIME        ! starting Julian date and time
    INTEGER         STPOS               ! starting position in ef day array
    INTEGER         TSTEP               ! output time step
    INTEGER         TZONE               ! output-file time zone
    INTEGER         TZMIN               ! minimum time zone in inventory
    INTEGER         TZMAX               ! maximum time zone in inventory
    INTEGER      :: TDMAX = 0           ! maximum episode days
    INTEGER         LDATE               ! date used in previous subroutine call

    REAL            RTMP                ! tmp float

    LOGICAL      :: DAYLIT    = .FALSE. ! true: TZONES are in daylight time
    LOGICAL         DFLAG               !  true: day-specific  file available
    LOGICAL      :: EFLAG = .FALSE.     !  error-flag
    LOGICAL      :: EFLAG2= .FALSE.     !  error-flag (2)
    LOGICAL         ENDFLAG             !  true: couldn't find file end date
    LOGICAL      :: FNDOUTPUT = .FALSE. ! true: found output hydrocarbon
    LOGICAL         HFLAG               !  true: hour-specific file available
    LOGICAL         NFLAG               !  true: use all uniform temporal profiles
    LOGICAL         PFLAG               !  true: episode time periods needed
    LOGICAL         WFLAG               !  true: write QA on current time step

    CHARACTER(8)         TREFFMT   ! tmprl x-ref format (SOURCE|STANDARD)
    CHARACTER(8)         SCDATE    ! tmprl date
    CHARACTER(14)        DTBUF     ! buffer for MMDDYY
    CHARACTER(3)         INTBUF    ! buffer for integer
    CHARACTER(20)        MODELNAM  ! emission factor model name
    CHARACTER(256)       CURFNM    ! current emission factor file name
    CHARACTER(16)        CURLNM    ! current ef logical file name
    CHARACTER(NAMLEN3)   VOLNAM    ! volatile pollutant name
    CHARACTER(300)       MESG      ! buffer for M3EXIT() messages
    CHARACTER(NAMLEN3)   CBUF      ! pollutant name temporary buffer
    CHARACTER(NAMLEN3)   EBUF      ! pollutant name temporary buffer
    CHARACTER(20)        SEARCHSTR ! string used in search
    CHARACTER(MXDLEN3)   TEMPLINE  ! line from file description


!***********************************************************************
!   begin body of program TEMPORAL

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

!.........  Obtain settings from the environment...

!.........  Get the time zone for output of the emissions
    TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTZONE"', 2 )
    END IF

!.........  Output hourly emission in Local Time.
    MESG = 'Outputs hourly emissions in local time  (No time shift)'
    LTFLAG = ENVYN( 'OUTPUT_LOCAL_TIME', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTPUT_LOCAL_TIME"', 2 )
    END IF

!.........  Get environment variable that overrides temporal profiles and
!               uses only uniform profiles.
    NFLAG = ENVYN( 'UNIFORM_TPROF_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UNIFORM_TPROF_YN"', 2 )
    END IF

!.........  Set source category based on environment variable setting
    CALL GETCTGRY

!.........  Get the name of the emission factor model to use for one run
!        IF ( CATEGORY .EQ. 'MOBILE' ) THEN
!            MESG = 'Emission factor model'
!            CALL ENVSTR( 'SMK_EF_MODEL', MESG, '', MODELNAM, IOS)
!        ELSE
!            MODELNAM = ' '
!        END IF

!.........  Get inventory file names given source category
    CALL GETINAME( CATEGORY, ENAME, ANAME )

!.........  Prompt for and open input files
!.........  Also, store source-category specific information in the MODINFO
!           module.
    CALL OPENTMPIN( NFLAG, PFLAG, ENAME, ANAME, DNAME,&
    &                HNAME, SDEV, CDEV, HDEV, KDEV, PYEAR )

!.........  Determine status of some files for program control purposes
    DFLAG = ( DNAME .NE. 'NONE' )  ! Day-specific emissions
    HFLAG = ( HNAME .NE. 'NONE' )  ! Hour-specific emissions

!.........  Get length of inventory file name
    ENLEN = LEN_TRIM( ENAME )

!.........  Set inventory variables to read for all source categories
    IVARNAMS( 1 ) = 'CIFIP'
    IVARNAMS( 2 ) = 'TZONES'
    IVARNAMS( 3 ) = 'TPFLAG'
    IVARNAMS( 4 ) = 'CSCC'
    IVARNAMS( 5 ) = 'CSOURC'

!.........  Set inventory variables to read for specific source categories
    IF( CATEGORY .EQ. 'AREA' ) THEN
        NINVARR = 5

    ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
        NINVARR = 9
        IVARNAMS( 6 ) = 'IRCLAS'
        IVARNAMS( 7 ) = 'IVTYPE'
        IVARNAMS( 8 ) = 'CLINK'
        IVARNAMS( 9 ) = 'CVTYPE'

    ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
        NINVARR = 5

    END IF

!.........  Allocate memory for and read in required inventory characteristics
    CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

!.........  Reset TPFLAG if average day emissions are being used since
!           we don't want to apply the monthly adjustment factors in this case.
    IF ( INVPIDX .EQ. 1 ) THEN
        DO S = 1, NSRC
            IF ( MOD( TPFLAG( S ), MTPRFAC ) .EQ. 0 ) THEN
                TPFLAG( S ) = TPFLAG( S ) / MTPRFAC
            END IF
        END DO
    END IF

!.........  Build unique lists of SCCs per SIC from the inventory arrays
    CALL GENUSLST

!.........  Define the minimum and maximum time zones in the inventory
    TZMIN = MINVAL( TZONES )
    TZMAX = MAXVAL( TZONES )

!.........  Adjust TZMIN and TZMAX for possibility of daylight savings
    TZMIN = TZMIN - 1
    TZMAX = TZMAX + 1

!.........  Read special files...
!.........  Read region codes file
    IF( USEEXPGEO() ) THEN
        CALL RDGEOCODES( NINVIFIP, INVCFIP )
    ELSE
        CALL RDSTCY( CDEV, NINVIFIP, INVCFIP )
    END IF

!.........  Populate filter for sources that use daylight time
    CALL SETDAYLT

!.........  Output hourly emissions in local time
    IF( LTFLAG ) THEN
        FLTRDAYL = 0         ! reset daily light saving value to zero
        TZONE  = 0           ! reset output time zone to zero
        TZONES = 0           ! reset county-specific time zone to zero
        TZMIN  = 0
        TZMAX  = 0
    END IF

!.........  Read holidays file
    CALL RDHDAYS( HDEV, EARLST, LATEST )

!.........  Determine all of the variables to be output based on the
!           activities and input pollutants.
    CALL TMNAMUNT

!.........  Reset the number of all output variables as the number of pollutants
!           and emission types, instead of the number of pollutants and
!           activities
    NIPPA = NIPOL

!.........  Allocate memory for I/O pol names, activities, & emission types
!.........  Will be resetting EANAM to include the emission types instead
!           of the activities
    DEALLOCATE( EANAM )
    ALLOCATE( EANAM( NIPPA ),&
    &          ALLIN( NIPPA ),&
    &       EAREAD2D( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EAREAD2D', PROGNAME )
    EAREAD2D = ' '

!.........  Create 1-d arrays of I/O pol names
!.........  If pollutant is created by one or more activities, then give a
!           warning.
    N = 0
    DO I = 1, NIPOL

!.............  Look for pollutant in list of pollutants created by activities
!               First, make sure that EMTPOL has been allocated; it won't be
!               when not using emission factors
        IF( ALLOCATED( EMTPOL ) ) THEN
            J = INDEX1( EINAM( I ), NEPOL, EMTPOL )
        ELSE
            J = 0
        END IF

!.............  If pollutant created by activity, skip from this list, unless
!               pollutant is also part of the inventory pollutants
        IF( J .GT. 0 ) THEN
            MESG = 'WARNING: Pollutant "' // TRIM( EINAM ( I ) )//&
            &       '" is explicitly in the inventory and' //&
            &       CRLF() // BLANK10 // 'it is also generated by '&
            &       // 'activity data.'
            CALL M3MSG2( MESG )
        END IF

        N = N + 1
        ALLIN( N ) = EAREAD( I )
        EANAM( N ) = EINAM ( I )

    END DO

!.........  Reset number of pollutants and emission types based on those used

    NIPPA = N

!.........  Read episode time period lists from PROCDATES.txt
!.........  or set NTPERIOD=1 and allocate RUNLEN(:)
!.........  Get episode settings from the Models-3 environment variables
!.........  when $GE_DAT/procdates.txt is not available for episode time periods

    CALL RDDATES( PFLAG, KDEV, NTPERIOD )

    IF ( .NOT.PFLAG ) THEN
        SDATE  = 0
        STIME  = 0
        NSTEPS = 1
        CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )

        TSTEP = 10000           ! Only 1-hour time steps supported
        JYEAR = SDATE / 1000
        CALL DAYMON( SDATE, JMNTH, JDAYS )
        STDATE(1) = JDAYS + 100*( JMNTH + 100*JYEAR )
        STTIME(1) = STIME
        RUNLEN(1) = NSTEPS
        ITDATE(1) = SDATE
    END IF

    DO II = 1, NTPERIOD

        NDAYS = 0

!.............  Get episode settings from episode time periods file
        IF( PFLAG ) THEN

            WRITE( SCDATE,'(I8)' ) STDATE( II )
            JYEAR = STR2INT( SCDATE( 1:4 ) )
            JMNTH = STR2INT( SCDATE( 5:6 ) )
            JDAYS = STR2INT( SCDATE( 7:  ) )
            SDATE = JULIAN( JYEAR, JMNTH, JDAYS )
            ITDATE( II ) = JYEAR*1000 + SDATE

            SDATE = ITDATE( II )
            STIME = STTIME( II )
            TSTEP  = 10000  ! Only 1-hour time steps supported
            NSTEPS = RUNLEN ( II ) / TSTEP
        END IF

!............  Determine number of days in episode
        IF( II == 1 ) THEN
            EARLST = SDATE
            LATEST = SDATE
        END IF

!............  Earliest day is start time in maximum time zone
        EARLYDATE = SDATE
        EARLYTIME = STIME
        CALL NEXTIME( EARLYDATE, EARLYTIME,&
        &             -( TZMAX - TZONE )*10000 )

!............  Latest day is end time in minimum time zone
!............  Calculate the ending date and time
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

        LATEDATE = EDATE
        LATETIME = ETIME
        CALL NEXTIME( LATEDATE, LATETIME,&
        &             -( TZMIN - TZONE )*10000 )

        NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
        NDAYS = NDAYS + 1
        IF( NDAYS > TDMAX ) TDMAX = NDAYS

!..............  Check earliest and latest dates during episode
        IF( EARLYDATE < EARLST ) EARLST = EARLYDATE
        IF( LATEDATE  > LATEST ) LATEST = LATEDATE

    END DO ! end of entire time episode periodes loops


!.........  Read temporal-profile cross-reference file and profiles, and
!.........  put into tables
!.........  Only read entries for pollutants that are in the inventory.
!.........  Only read if not using uniform temporal profiles
!.........  Assign temporal profiles to sources.

    CALL PROCTPRO( NFLAG, METPROFLAG, PNAME )

    IF ( .NOT.NFLAG ) CALL ASGNTPRO()

!.........  Check requested episode against available emission factors

!.........  For day-specific data input...
    IF( DFLAG ) THEN

!.............  Get header description of day-specific input file
        IF( .NOT. DESC3( DNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0,&
            &             'Could not get description of file "'&
            &             // DNAME( 1:LEN_TRIM( DNAME ) ) // '"', 2 )
        END IF

!.............  Allocate memory for pollutant pointer
        ALLOCATE( DYPNAM( NVARS3D ),&
        &          DYPDSC( NVARS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DYPDSC', PROGNAME )
        DYPNAM = ' '  ! array
        DYPDSC = ' '  ! array

!C.............  Set day-specific file dates, check dates, and report problems
        CALL PDSETUP( DNAME, EARLST, STIME, LATEST, ETIME, TZONE,&
        &              NIPPA, EANAM, NDYPOA, NDYSRC, EFLAG, DYPNAM,&
        &              DYPDSC )

    ENDIF

!.........  Allocate memory for reading day-specific emissions data
!.........  NDYSRC is initialized to zero in case DFLAG is false
    ALLOCATE( INDXD( NDYSRC ),&
    &          EMACD( NDYSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMACD', PROGNAME )

!.........  For hour-specific data input...
    IF( HFLAG ) THEN

!............. Get header description of hour-specific input file
        IF( .NOT. DESC3( HNAME ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0,&
            &             'Could not get description of file "'&
            &             // HNAME( 1:LEN_TRIM( HNAME ) ) // '"', 2 )
        ENDIF

!.............  Allocate memory for pollutant pointer
        ALLOCATE( HRPNAM( NVARS3D ),&
        &          HRPDSC( NVARS3D ), STAT=IOS )
        CALL CHECKMEM( IOS, 'HRPDSC', PROGNAME )
        HRPNAM = ' '  ! array
        HRPDSC = ' '  ! array

!.............  Set day-specific file dates, check dates, and report problems
        CALL PDSETUP( HNAME, EARLST, STIME, LATEST, ETIME, TZONE,&
        &              NIPPA, EANAM, NHRPOA, NHRSRC, EFLAG2, HRPNAM,&
        &              HRPDSC )

    ENDIF

!.........  Allocate memory for reading hour-specific emissions data
!.........  NHRSRC is initialized to 0 in case HFLAG is false
    ALLOCATE( INDXH( NHRSRC ),&
    &          EMACH( NHRSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMACH', PROGNAME )

    IF( EFLAG .OR. EFLAG2 ) THEN
        MESG = 'Problem with day- or hour-specific inputs'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

!......  Get episode settings from episode time periods file

    DO II = 1, NTPERIOD

        LDATE = -1     ! initializing previous date

!.........  Initialize EAREAD2D every processing date
        DO I = 1, NIPPA
            EAREAD2D( I ) = EAREAD( I )
        END DO

!.........  It is important that all major arrays must be allocated by this
!           point because the next memory allocation step is going to pick a
!           data structure that will fit within the limits of the host.

!.........  Allocate memory, but allow flexibility in memory allocation
!           for second dimension.
!.........  The second dimension (the number of pollutants and emission types)
!            can be different depending on the memory available.
!.........  To determine the approproate size, first attempt to allocate memory
!           for all pollutants & emission types to start, and if this fails,
!           divide pollutants into even groups and try again.

        NGSZ = NIPPA            ! No. of pollutant & emis types in each group
        NGRP = 1                ! Number of groups

        DO

            ALLOCATE( TMAT ( NSRC, NGSZ, 24 ),&
            &          EMAC ( NSRC, NGSZ )    ,&
            &          EMACV( NSRC, NGSZ )    ,&
            &          EMIST( NSRC, NGSZ )    ,&
            &          EMFAC( NSRC, NGSZ )    , STAT=IOS9 )

            IF( IOS9 .GT. 0 ) THEN

                IF( NGSZ .EQ. 1 ) THEN
                    J = 8 * NSRC * 31    ! Assume 8-byte reals
                    WRITE( MESG,94010 )&
                    &  'Insufficient memory to run program.' //&
                    &  CRLF() // BLANK5 // 'Could not allocate ' //&
                    &  'pollutant-dependent block of', J, 'bytes.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NGRP = NGRP + 1
                NGSZ = ( NIPPA + NGRP - 1 ) / NGRP

                IF( ALLOCATED( TMAT  ) ) DEALLOCATE( TMAT )
                IF( ALLOCATED( EMAC  ) ) DEALLOCATE( EMAC )
                IF( ALLOCATED( EMACV ) ) DEALLOCATE( EMACV )
                IF( ALLOCATED( EMIST ) ) DEALLOCATE( EMIST )
                IF( ALLOCATED( EMFAC ) ) DEALLOCATE( EMFAC )

            ELSE
                EXIT

            END IF

        END DO

        WRITE( MESG, '( 2( A, I4, 2X ) )' )&
        &   'Processing', NGRP, 'pollutant groups of size', NGSZ
        CALL M3MESG( MESG )

!.........  Allocate a few small arrays based on the size of the groups
!.........  NOTE that this has a small potential for a problem if these little
!           arrays exceed the total memory limit.
        IF( ALLOCATED( ALLIN2D ) ) DEALLOCATE( ALLIN2D )
        IF( ALLOCATED( EANAM2D ) ) DEALLOCATE( EANAM2D )
        IF( ALLOCATED( IPOL2D  ) ) DEALLOCATE( IPOL2D  )
        IF( ALLOCATED( LDSPOA  ) ) DEALLOCATE( LDSPOA  )
        IF( ALLOCATED( LHSPOA  ) ) DEALLOCATE( LHSPOA  )
        IF( ALLOCATED( LHPROF  ) ) DEALLOCATE( LHPROF  )

        ALLOCATE( ALLIN2D( NGSZ, NGRP ),&
        &          EANAM2D( NGSZ, NGRP ),&
        &           IPOL2D( NGSZ, NGRP ),&
        &           LDSPOA( NGSZ ),&
        &           LHSPOA( NGSZ ),&
        &           LHPROF( NGSZ ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ALLIN2D...LHPROF', PROGNAME )

!.........  Create 2-d arrays of I/O pol names, activities, & emission types

        I = 0
        DO N = 1, NGRP
            DO J = 1, NGSZ
                I = I + 1
                IF ( I .LE. NIPPA ) THEN
                    IPOL2D ( J,N ) = I
                    EANAM2D( J,N ) = EANAM( I )
                    ALLIN2D( J,N ) = ALLIN( I )
                ELSE
                    IPOL2D ( J,N ) = 0
                    EANAM2D( J,N ) = ' '
                    ALLIN2D( J,N ) = ' '
                END IF
            END DO
        END DO

!......  Determine number of days in episode

        IF( PFLAG ) THEN
            SDATE  = ITDATE( II )
            STIME  = STTIME( II )
            TSTEP  = 10000  ! Only 1-hour time steps supported
            NSTEPS = RUNLEN ( II ) / TSTEP
        END IF

!.........  Earliest day is start time in maximum time zone
        EARLYDATE = SDATE
        EARLYTIME = STIME
        CALL NEXTIME( EARLYDATE, EARLYTIME,&
        &             -( TZMAX - TZONE )*10000 )

!.........  Latest day is end time in minimum time zone
        EDATE = SDATE
        ETIME = STIME
        CALL NEXTIME( EDATE, ETIME, NSTEPS * 10000 )

        LATEDATE = EDATE
        LATETIME = ETIME
        CALL NEXTIME( LATEDATE, LATETIME,&
        &             -( TZMIN - TZONE )*10000 )

        NDAYS = SECSDIFF( EARLYDATE, 0, LATEDATE, 0 ) / ( 24*3600 )
        NDAYS = NDAYS + 1

!.....  Allocate memory for emission factor arrays
        IF( ALLOCATED( TEMPEF ) ) DEALLOCATE ( TEMPEF )
        ALLOCATE( TEMPEF( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TEMPEF', PROGNAME )

        TEMPEF = 0.

!.........  Compare base year with episode and warn if not consistent
        IF( SDATE / 1000 .NE. BYEAR ) THEN

            WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR,&
            &       'is inconsistent with year ' // CRLF() // BLANK10 //&
            &       'of episode start date', SDATE/1000
            CALL M3MSG2( MESG )

        ENDIF

!.........  Set up and open I/O API output file(s) ...
        CALL OPENTMP( II, ENAME, SDATE, STIME, TSTEP, NSTEPS, TZONE,&
        &              NPELV, TNAME, PDEV, PFLAG )

!.........  Give a note if running for a projected year
        IF( PYEAR .GT. 0 ) THEN

            WRITE( MESG,94010 ) 'NOTE: Emissions based on projected '//&
            &       'year', PYEAR
            CALL M3MSG2( MESG )

        END IF

!.........  Loop through pollutant/emission-type groups
        DO N = 1, NGRP

!.............  If this is the last group, NGSZ may be larger than actual
!               number of variables, so reset based on total number
            IF( N == NGRP ) THEN
                NGSZ = NIPPA - ( NGRP - 1 )*NGSZ
            END IF

!.............  Skip group if the first pollutant in group is blank (this
!               shouldn't happen, but it is happening, and it's easier to
!               make this fix).
            IF ( EANAM2D( 1,N ) .EQ. ' ' ) CYCLE

!.............  Write message stating the pols/emission-types being processed
            CALL POLMESG( NGSZ, EANAM2D( 1,N ) )

!.............  Set up logical arrays that indicate which pollutants/activities
!               are day-specific and which are hour-specific.
!.............  Also set flag for which hour-specific pollutants/activities
!               are actually diurnal profiles instead of emissions
            LDSPOA = .FALSE.   ! array
            DO I = 1, NDYPOA
                J = INDEX1( DYPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                IF ( J /= 0 ) LDSPOA( J ) = .TRUE.
            END DO

            LHSPOA = .FALSE.   ! array
            LHPROF = .FALSE.   ! array
            DO I = 1, NHRPOA
                J = INDEX1( HRPNAM( I ), NGSZ,  EANAM2D( 1,N ) )
                IF ( J /= 0 ) LHSPOA( J ) = .TRUE.

                CALL UPCASE( HRPDSC( I ) )
                K = INDEX( HRPDSC( I ), 'PROFILE' )
                IF( K .GT. 0 ) LHPROF( J ) = .TRUE.
            END DO

!.............  Initialize emissions, activities, and other arrays for this
!.............  pollutant/emission-type group

            TMAT  = 0.
            EMAC  = 0.
            EMACV = 0.
            EMIST = 0.

!.............  Read in pollutant emissions or activities from inventory for
!               current group

            DO I = 1, NGSZ

                EBUF = EANAM2D( I,N )
                CBUF = ALLIN2D( I,N )
                L1   = LEN_TRIM( CBUF )

!.................  Skip blanks that can occur when NGRP > 1
                IF ( CBUF .EQ. ' ' ) CYCLE

!.................  Read the emissions data in either map format
!                   or the old format.
                CALL RDMAPPOL( NSRC, 1, 1, CBUF, EMAC( 1,I ) )

!...............  If there are any missing values in the data, give an
!                 error to avoid problems in genhemis routine
                RTMP = MINVAL( EMAC( 1:NSRC,I ) )
                IF( RTMP .LT. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Missing or negative emission '//&
                    &       'value(s) in inventory for "' //&
                    &       CBUF( 1:L1 ) // '".'
                    CALL M3MSG2( MESG )
                END IF

!...............  If pollutant name is average day-based, remove the
!                 prefix from the input pollutant name
                K = INDEX1( CBUF, NIPPA, EAREAD )
                J = INDEX( CBUF, AVEDAYRT )
                IF( J .GT. 0 ) THEN
                    CBUF = CBUF( CPRTLEN3+1:L1 )
                    ALLIN2D( I,N ) = CBUF
                    EAREAD2D ( K ) = CBUF
                END IF

            END DO

!.............  Abort if error found
            IF( EFLAG ) THEN
                MESG = 'Problem with input data.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!.............  For each time step and pollutant or emission type in current
!               group, generate hourly emissions, and write layer-1 emissions
!               file (or all data).
            JDATE = SDATE
            JTIME = STIME

            DO T = 1, NSTEPS        !  Loop through time steps for current pollutant group

                WRITE( MESG, '( A, I4, 2X, A, I9.7, A, I6.6 )' )&
                &      'Processing pol-group', N, 'for', JDATE, ':', JTIME
                CALL M3MESG( MESG )

!.................  Adjust sources' time zones to account for daylight time...
!                   Subtract 1 if date is daylight time and TZONES is not already
!                   converted.  Add 1 if date is standard and TZONES has been converted.
!                   FLTRDAYL is a source-array of 0s and 1s to permit sources
!                   to not get daylight time conversion.

                IF( ISDSTIME( JDATE ) .AND. .NOT. DAYLIT ) THEN

                    DAYLIT = .TRUE.
                    TZONES = TZONES - 1 * FLTRDAYL   ! arrays

                ELSE IF( .NOT. ISDSTIME( JDATE ) .AND. DAYLIT ) THEN

                    DAYLIT = .FALSE.
                    TZONES = TZONES + 1 * FLTRDAYL   ! arrays

                END IF

!.................  Generate and write hourly emissions for current hour

                CALL GENHEMIS( N, NGRP, NGSZ, JDATE, JTIME, TZONE, DNAME, HNAME,&
                &        PNAME, ALLIN2D( 1,N ), EANAM2D( 1,N ), EAREAD2D,&
                &        LDATE )

                DO I = 1, NGSZ      !  Loop through pollutants/emission-types in this group
                !  Skip blanks that can occur when NGRP > 1
                    CBUF = EANAM2D( I,N )
                    IF ( CBUF .EQ. ' ' ) CYCLE

                    IF( .NOT. WRITESET( TNAME, CBUF, ALLFILES, JDATE,&
                    &                    JTIME, EMIST( 1,I ) )    ) THEN
                        MESG = 'Could not write "' // TRIM( CBUF ) //&
                        &       '" to file "' // TRIM( TNAME ) // '"'
                        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

                    END IF

                END DO  ! End loop on pollutants/emission-types I in this group

!.................  Advance the output date/time by one time step

                CALL NEXTIME( JDATE, JTIME, TSTEP )

            END DO      ! End loop on time steps T

        END DO         ! End loop on pollutant groups N

!.............  Write supplemental temporal profiles file

        CALL WRTSUP( PDEV, NSRC, NIPPA, EANAM )

        IF( .NOT. CLOSESET( TNAME ) ) THEN
            MESG = 'Could not close file "' // TRIM( TNAME ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        DEALLOCATE( TMAT, EMAC, EMACV, EMIST, EMFAC )

    END DO            ! End loop on time period II

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94050 FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,&
    &        A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

END PROGRAM TEMPORAL

