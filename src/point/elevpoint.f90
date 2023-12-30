
PROGRAM ELEVPOINT

    !***********************************************************************
    !  program body starts at line
    !
    !  DESCRIPTION:
    !       Identifies sources as elevated (major or plume-in-grid) or minor.
    !       Major sources will get plume rise and minor sources will not, however,
    !       the major/minor distinction is not required because SMOKE will compute
    !       layer fractions for all sources efficiently when needed.  If desired,
    !       the program can use an analytical computation (the PLUMRIS routine)
    !       along with a cutoff height to determine the major sources.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Copied from elevpoint.F 4.2 by M Houyoux
    !       12/5/2007: 1) Added code to select fire sources so the stack parameter related
    !                  calculations can be skipped for in-line pt source processing CMAQv4.7
    !
    !                  2) disabled checking of stack paramters for analytical plume rise because
    !                  default values are used when stack diameter is zero for example
    !
    !                  3) added new option for select sources for in-line plume rise in CMAQ
    !                      SMK_ELEV_METHOD=2
    !                   George Pouliot
    !
    !                  4) Allow the passing through of ACRES from PDAY file into stack
    !                      groups file 3/18/08
    !
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !************************************************************************
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
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: XLOCA, YLOCA, STKDM, STKHT, STKTK, STKVE,       &
                        CSOURC, CIFIP, CPDESC, CSCC, CNAICS

    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: LMAJOR, LPING, LCUTOFF, GROUPID, GINDEX,         &
                       NGRPVAR, NEVPVAR, NELVCRIT, MXELVCHK,            &
                       NPNGCRIT, MXPNGCHK,  NGRPCRIT, MXGRPCHK,         &
                       GRPVALS, GRPTYPES, NEVPEMV, NGROUP,              &
                       ELVVALS, ELVCHRS, ELVTYPES, PNGVALS, PNGCHRS,    &
                       PNGTYPES, GRPGID, GRPCNT, GRPLAT, GRPLON,        &
                       GRPDM, GRPHT, GRPTK, GRPVE, GRPFL, GRPGIDA,      &
                       GRPIDX, GRPCOL, GRPROW, GRPXL, GRPYL, RISE,      &
                       GRPFIP, GRPLMAJOR, GRPLPING,                     &
                       MXEMIS, MXRANK, EVPEMIDX, SRCXL, SRCYL,          &
                       DAY_ACRES, FFLAG, DAY_INDEX, ACRES, GRPACRES

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, CATLEN, NSRC, MXCHRS,     &
                       SC_BEGP, SC_ENDP, NCHARS, EINAM, JSTACK

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GDTYP, GRDNM, P_ALP, P_BET, P_GAM,       &
                       XCENT, YCENT, NCOLS, NROWS, NGRID

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'CONST3.EXT'        !  physical and mathematical constants from I/O API

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: DSCM3GRD
    LOGICAL, EXTERNAL :: EVALCRIT
    REAL   , EXTERNAL :: PLUMRIS

    !.......  LOCAL PARAMETERS and their descriptions:
    CHARACTER(16), PARAMETER :: PROGNAME = 'ELEVPOINT'       !  program name
    CHARACTER(50), PARAMETER ::     CVSW = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    INTEGER, PARAMETER :: LAYPOINT_APPROACH   = 0
    INTEGER, PARAMETER :: NOPING_APPROACH     = 0
    INTEGER, PARAMETER :: PELVCONFIG_APPROACH = 2


    !.......   Indicator for which public inventory arrays need to be read
    INTEGER,            PARAMETER :: NINVARR = 12
    CHARACTER(NAMLEN3), PARAMETER :: IVARNAMS( NINVARR ) =  &
                                   (/ 'CIFIP          '     &
                                     , 'TZONES         '    &
                                     , 'XLOCA          '    &
                                     , 'YLOCA          '    &
                                     , 'STKHT          '    &
                                     , 'STKDM          '    &
                                     , 'STKTK          '    &
                                     , 'STKVE          '    &
                                     , 'CSOURC         '    &
                                     , 'CNAICS         '    &
                                     , 'CPDESC         '    &
                                     , 'CSCC           ' /)

    !.......   Descriptions of plume-in-grid and elevated source approaches
    INTEGER, PARAMETER :: NELEVMTHD = 3
    INTEGER, PARAMETER :: NPINGMTHD = 2
    CHARACTER(60), PARAMETER :: ELEVMTHD( 0:NELEVMTHD-1 ) = (/              &
        'Allow Laypoint to determine elevated sources                ',     &
        'Use PELVCONFIG file to determine elevated sources           ',     &
        'Use PELVCONFIG file to determine srcs for in-line plume rise'      /)

    CHARACTER(60), PARAMETER :: PINGMTHD( 0:NPINGMTHD-1 ) = (/              &
        'No PinG sources                                             ',     &
        'Use PELVCONFIG file to determine PinG sources               '      /)

    !.......   Allocateable arrays for using GENPTCEL routine to get grid-cell
    !              numbers based on current projection
    INTEGER, ALLOCATABLE :: INDX ( : )      ! sorting index (unused)
    INTEGER, ALLOCATABLE :: GN   ( : )      ! cell numbers
    INTEGER, ALLOCATABLE :: SN   ( : )      ! stack group pos in list (unused)
    INTEGER, ALLOCATABLE :: NX   ( : )      ! no. stack groups per cell (unused)

    !.......   Allocatable arrays for temporary storage of inventory group info
    INTEGER, ALLOCATABLE :: LOCGID ( : )      ! group number
    INTEGER, ALLOCATABLE :: LOCCNT ( : )      ! number sources in group
    LOGICAL, ALLOCATABLE :: LOCSTAT( : )      ! true: group reset during renumbering

    !.......   Temporary by-source arrays
    INTEGER, ALLOCATABLE :: SRCGROUP( : )     ! group number in sorted src order

    !.......   Other allocatable arrays
    REAL   , ALLOCATABLE :: VALS( : )      ! tmp test values
    REAL   , ALLOCATABLE :: RANK( : )      ! tmp ranked value
    LOGICAL, ALLOCATABLE :: EVSTAT( :,:,: )     ! tmp status of elev checks
    LOGICAL, ALLOCATABLE :: PGSTAT( :,:,: )     ! tmp status of PiNG checks
    LOGICAL, ALLOCATABLE :: SMOLDER (:)
    CHARACTER(PLTLEN3), ALLOCATABLE :: CHRS( : )     ! tmp test character strings

    !.......   File units and logical/physical names
    INTEGER         CDEV        !  elevated source configuration file
    INTEGER         IDEV        !  tmp unit number if ENAME is map file
    INTEGER         LDEV        !  log-device
    INTEGER         PDEV        !  for output major/mepse src ID file
    INTEGER         RDEV        !  ASCII output report
    INTEGER         SDEV        !  ASCII part of inventory unit no.

    CHARACTER(16)   ANAME       !  logical name for ASCII inventory input file
    CHARACTER(16)   ENAME       !  logical name for i/o api inventory input file
    CHARACTER(16)   INAME       !  tmp name for inven file of unknown fmt
    CHARACTER(16)   MNAME       !  plume-in-grid srcs stack groups output file

    !.......   Other local variables
    INTEGER         G, I, J, K, DS, S, L, L2, N, V, I1,J1,T        ! indices and counters

    INTEGER         COL               ! tmp column number
    INTEGER      :: ELEVTYPE = 0      ! code for elevated source approach
    INTEGER         IGRP              ! tmp group ID
    INTEGER         IOS               ! i/o status
    INTEGER         IOSCUT            ! i/o status for cutoff E.V.
    INTEGER         IREC              ! record counter
    INTEGER      :: NEXCLD = 0        ! no. stack groups exlcuded from the grid
    INTEGER         NINVGRP           ! no. inventory groups
    INTEGER      :: NMAJOR = 0        ! no. major sources
    INTEGER      :: NMJRGRP = 0       ! no. major sources incl groups
    INTEGER      :: NPING  = 0        ! no. plume-in-grid sources
    INTEGER      :: NPINGGRP = 0      ! no. plume-in-grid sources incl groups
    INTEGER         NSLINES           ! no. lines in stack splits file
    INTEGER      :: NSTEPS = 24       ! no. time steps
    INTEGER         OUTG              ! group number for output report
    INTEGER         MS                ! tmp src ID for major sources
    INTEGER         PEGRP             ! grp no. for elev/ping from prev iteratn
    INTEGER         PGRP              ! group no. from previous iteration
    INTEGER      :: PINGTYPE = 0      ! code for PinG source approach
    INTEGER         PLTEND            ! end position for plant string
    INTEGER         PS                ! tmp src ID for plume in grid sources
    INTEGER         ROW               ! tmp row number
    INTEGER      :: SDATE = 0         ! Julian start date
    INTEGER      :: STIME = 0         ! start time
    INTEGER      :: TSTEP = 10000     ! time step HHMMSS
    INTEGER      :: TSTEP_T           ! unsued timestep from environment
    INTEGER         TZONE             ! output time zone
    INTEGER      :: DAY_NSRC
    INTEGER      :: JDATE, JTIME

    REAL            DM, DMVAL         ! tmp inside stack diameter [m]
    REAL            FL                ! tmp stack exit flow rate [m^3/s]
    REAL            HT                ! tmp inside stack diameter [m]
    REAL            LAT               ! tmp latitude [degrees]
    REAL            LON               ! tmp longitude [degrees]
    REAL            MINDM             ! min stack group diam
    REAL            MINFL             ! min stack group flow
    REAL            MINHT             ! min stack group height
    REAL            MINTK             ! min stack group temperature
    REAL            MINVE             ! min stack group velocity
    REAL            TK                ! tmp stack exit temperature [K]
    REAL            VE, VEVAL         ! tmp stack exit velocity diameter [m/s]

    LOGICAL :: EFLAG    = .FALSE.     ! true: error detected
    LOGICAL :: SFLAG    = .FALSE.     ! true: store group info
    LOGICAL    VFLAG                  ! true: use variable grid
    LOGICAL :: RFLAG    = .TRUE.      ! true: skip checking fake source
    LOGICAL :: MFLAG    = .TRUE.      ! true: fill fake source
    LOGICAL :: LFLAG    = .TRUE.      ! true: write out lat/lon info

    CHARACTER(FIPLEN3) CFIP           ! tmp country/st/county code
    CHARACTER(SCCLEN3) SCC
    CHARACTER(80)   GDESC         !  grid description
    CHARACTER(512)  BUFFER
    CHARACTER(256)  MESG
    CHARACTER(16) DAYNAME       !  daily inventory file name


    CHARACTER(NAMLEN3) COORD3D      !  coordinate system name
    CHARACTER(NAMLEN3) COORUN3D     !  coordinate system units
    CHARACTER(ALLLEN3) CSRC         !  buffer for source char, incl pol/act
    CHARACTER(PLTLEN3) PLT          !  tmp plant code
    CHARACTER(CHRLEN3) STK          !  tmp stack code

    !***********************************************************************
    !   begin body of program ELEVPOINT

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Get environment variables that control this program
    MESG = 'Approach for create plume-in-grid outputs'
    PINGTYPE = ENVINT( 'SMK_PING_METHOD', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_PING_METHOD"', 2 )
    END IF

    MESG = 'Approach for defining major/minor sources'
    ELEVTYPE = ENVINT( 'SMK_ELEV_METHOD', MESG, 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_ELEV_METHOD"', 2 )
    END IF

    !.......  Define whether write out a fake source when there is no group
    MESG = 'Define whether write out a fake source when ' //    &
           'there is no group to output'
    MFLAG = ENVYN( 'ELEV_WRITE_FAKE_SRC', MESG, .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ELEV_WRITE_FAKE_SRC"', 2 )
    END IF

    !.......  Define whether write out lat/lon info for the elevated sources
    MESG = 'Define whether write out lat/lon info for the elevated sources or not'
    LFLAG = ENVYN( 'ELEV_WRITE_LATLON', MESG, .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ELEV_WRITE_LATLON"', 2 )
    END IF

    IF( .NOT. LFLAG .AND. ( PINGTYPE == PELVCONFIG_APPROACH - 1 ) ) THEN
        MESG = 'ERROR: Cannot set ELEV_WRITE_LATLON to N when ' //    &
               'processing plume-in-grid (PinG) sources'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    !.......  End program for invalid settings
    IF ( ELEVTYPE .GT. PELVCONFIG_APPROACH .OR.    &
         PINGTYPE .GT. PELVCONFIG_APPROACH-1      ) THEN

        IF ( ELEVTYPE .GT. PELVCONFIG_APPROACH ) THEN
            MESG = 'WARNING: SMK_ELEV_METHOD value is invalid'//                            &
               CRLF() // BLANK10 // 'Valid values are: '//                                  &
               CRLF() // BLANK10 // '0 = Allow Laypoint to determine elevated sources'//    &
               CRLF() // BLANK10 // '1 = Use PELVCONFIG to determine elevated sources.'//   &
               CRLF() // BLANK10 // '2 = Use PELVCONFIG to determine sources for in-line.'//    &
               CRLF() // BLANK10 // 'Setting to 0.'

            ELEVTYPE = 0
            CALL M3MSG2( MESG )
        END IF

        IF ( PINGTYPE .GT. PELVCONFIG_APPROACH-1 ) THEN
            MESG = 'WARNING: SMK_PING_METHOD value is invalid'//        &
               CRLF() // BLANK10 // 'Valid values are: '//              &
               CRLF() // BLANK10 // '0 = No plume-in-grid sources'//    &
               CRLF() // BLANK10 // '1 = Use PELVCONFIG to determine plume-in-grid sources.' //    &
               CRLF() // BLANK10 // 'Setting to 0.'
            PINGTYPE = 0
            CALL M3MSG2( MESG )
        END IF

    END IF

    !.......  End program run if program inputs indicate it does not need to be
    !           run
    IF( ELEVTYPE .LE. LAYPOINT_APPROACH .AND.    &
        PINGTYPE .LE. NOPING_APPROACH         ) THEN
        MESG = 'ERROR: Neither plume-in-grid nor elevated sources will be identified '//        &
               CRLF()// BLANK10// 'based on SMK_ELEV_METHOD and SMK_PING_METHOD settings. ' //  &
               CRLF()// BLANK10//     'Elevpoint is not needed. Ending...'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    END IF

    !.......  Check if using a variable grid
    VFLAG = ENVYN( 'USE_VARIABLE_GRID',    &
                   'Use variable grid definition', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

    !.......  Set source category based on environment variable setting
    CALL GETCTGRY

    !.......  Get inventory file names given source category
    CALL GETINAME( CATEGORY, ENAME, ANAME )

    !.......  Make sure only run for point sources
    IF( CATEGORY .NE. 'POINT' ) THEN
        MESG = 'ERROR: ' // TRIM( PROGNAME ) //    &
               ' is not valid for ' // CATEGORY( 1:CATLEN ) //    &
               ' sources'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Write out notes about elevated and PinG approachs
    IF ( ELEVTYPE .GT. 0 ) THEN
        MESG = 'NOTE: Elevated source approach is...' //    &
                CRLF() // BLANK10 // ELEVMTHD( ELEVTYPE  )
        CALL M3MSG2( MESG )
    END IF

    MESG = 'NOTE: Plume-in-grid (PinG) approach is...'//    &
           CRLF() // BLANK10 // PINGMTHD( PINGTYPE )
    CALL M3MSG2( MESG )

    !.......   Get file name; open input point source and output
    !.......   elevated points files

    !.......  Prompt for and open inventory file
    INAME = ENAME
    MESG = 'Enter logical name for the MAP INVENTORY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.......  Open and read map file
    CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    !.......  Get elevated source configuration file, if needed
    IF( PINGTYPE .EQ.  (PELVCONFIG_APPROACH - 1) .OR.    &
        ELEVTYPE .EQ.  PELVCONFIG_APPROACH .OR.    &
        ELEVTYPE .EQ. (PELVCONFIG_APPROACH -1)      ) THEN

        CDEV = PROMPTFFILE(    &
           'Enter logical name for the ELEVATED SOURCE CONFIGURATION file',    &
           .TRUE., .TRUE., CRL // 'ELVCONFIG', PROGNAME )

    END IF

    !.......  Open ASCII report file
    RDEV = PROMPTFFILE(    &
          'Enter name for ELEVATED SELECTION REPORT file',    &
          .FALSE., .TRUE., 'REP' // CRL // 'ELV', PROGNAME )

    !.......  Store source-category-specific header information,
    !           including the inventory pollutants in the file (if any).  Note that
    !           the I/O API header info is passed by include file and the
    !           results are stored in module MODINFO.
    CALL GETSINFO( ENAME )

    !.......  Allocate memory for and read in required inventory characteristics
    CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

    !.......  If at least one stack parameters is missing, then we have a fire inventory
    DO J = 1, NSRC
        IF (STKHT(J) .NE. BADVAL3) FFLAG = .FALSE.
    END DO

    !.......  Allocate memory for source status arrays and group numbers
    ALLOCATE( LMAJOR( NSRC ),    &
               LPING( NSRC ),    &
             GROUPID( NSRC ),    &
              GINDEX( NSRC ),    &
            SRCGROUP( NSRC ),    &
               SRCXL( NSRC ),    &
               SRCYL( NSRC ),    &
              SMOLDER( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LMAJOR...SMOLDER', PROGNAME )
    SMOLDER(1:NSRC) = .FALSE.            ! array

    !.......  Initialize source status and group number arrays
    LMAJOR  = .FALSE.       ! array
    LPING   = .FALSE.       ! array
    GROUPID = 0             ! array
    GINDEX  = 0             ! array

    !.......  Get grid description for converting the stack coordinates
    !           to grid cells for the STACK_GROUPS file.
    IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D, &
                  P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,      &
                  XORIG3D, YORIG3D, XCELL3D, YCELL3D,               &
                  NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE
        GDTYP = GDTYP3D
        GRDNM = GDNAM3D
        P_ALP = P_ALP3D
        P_BET = P_BET3D
        P_GAM = P_GAM3D
        XCENT = XCENT3D
        YCENT = YCENT3D
        NCOLS = NCOLS3D
        NROWS = NROWS3D
        NGRID = NCOLS * NROWS
    END IF

    !.......  Convert source x,y locations to coordinates of the projected grid
    SRCXL = XLOCA
    SRCYL = YLOCA

    CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM, XCENT, YCENT, SRCXL, SRCYL )

    !.......  Allocate memory so that we can use the GENPTCEL
    ALLOCATE( NX( NGRID ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NX', PROGNAME )

    !.......  If needed, read config file

    IF( (ELEVTYPE .EQ.  PELVCONFIG_APPROACH   ) .OR.    &
        (PINGTYPE .EQ. (PELVCONFIG_APPROACH-1)) .OR.    &
        (ELEVTYPE .EQ. (PELVCONFIG_APPROACH-1))     ) THEN

        CALL RPELVCFG( CDEV )

    END IF

    !.......  Allocate memory for temporary arrays for use in Evalcrit
    I = MAX( NGRPVAR, NEVPVAR )
    ALLOCATE( VALS( I ),    &
              CHRS( I ),    &
              RANK( I ),    &
              EVSTAT( NELVCRIT, MXELVCHK, NEVPVAR ),    &
              PGSTAT( NPNGCRIT, MXPNGCHK, NEVPVAR ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VALS...PGSTAT', PROGNAME )
    VALS = 0.
    CHRS = ' '
    RANK = 0

    !.......  Allocate memory for analytical plume rise array, if it's needed
    IF( LCUTOFF ) THEN
        ALLOCATE( RISE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'RISE', PROGNAME )
    END IF

    !.......  Determine valid stack groups based on inventory. The stacks at the
    !           same plant with the same stack parameters can be assigned to the
    !           same group.
    !.......  There will be more groups possibly because of
    !           major and PinG sources that aren't in the inventory groups will
    !           get assigned their own group.
    !.......  This routine populates arrays in MODELEV that contain the group-
    !           specific information.  It allocates some of the group arrays.
    !.......  NGRPCRIT may be zero if no grouping criteria have been set, but
    !           the routine will still set groups for facility stacks that match
    !           exactly

    IF (.NOT. FFLAG ) THEN
        CALL ASGNGRPS( NGRPVAR, NGRPCRIT, MXGRPCHK,    &
                       GRPVALS, GRPTYPES, NINVGRP   )
    ELSE
        NINVGRP = NSRC
        DO J = 1, NSRC
            GINDEX(J) = J
        ENDDO

        ALLOCATE( GRPCNT( NINVGRP ),    &
                  GRPGID( NINVGRP ),    &
                  GRPLAT( NINVGRP ),    &
                  GRPLON( NINVGRP ),    &
                   GRPDM( NINVGRP ),    &
                   GRPHT( NINVGRP ),    &
                   GRPTK( NINVGRP ),    &
                   GRPVE( NINVGRP ),    &
                   GRPFL( NINVGRP ),    &
                  GRPFIP( NINVGRP ),    &
                GRPACRES( NINVGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPCNT...GRPACRES', PROGNAME )

        GRPACRES = 0.0
        GRPGID = 0
        GRPCNT = 0

    ENDIF

    !.......  If emissions are needed as a criteria, determine maximum daily
    !           emissions over the time period being processed for each source
    !           and inventory group.
    !.......  Also get date and time information a different way if using
    !           an emissions file on input
    IF( NEVPEMV .GT. 0 ) THEN

        !.......  Get episode information for setting date and time of STACK_PING file
        MESG = 'NOTE: Getting date/time information to constrain '//    &
               'time period for emissions input file'
        CALL M3MSG2( MESG )

        CALL GETM3EPI( -9, SDATE, STIME, TSTEP_T, NSTEPS )

        !.......  Create maximum daily emissions by stack group.
        !         The stack groups have already been set, and now the emissions for those
        !         groups must be computed to assign MEPSEs and MPSs.
        CALL MXGRPEMIS( NINVGRP, TSTEP, SDATE, STIME, NSTEPS )

    ELSE
        MESG = 'NOTE: Getting date/time information only for use in STACK_PING file'
        CALL M3MSG2( MESG )

        CALL GETM3EPI( -9, SDATE, STIME, TSTEP_T, -9 )

    END IF      ! End of whether emissions are needed as a criteria


    IF( FFLAG )THEN
        DAYNAME = PROMPTMFILE( 'Enter logical name for DAY-SPECIFIC file',    &
                               FSREAD3, CRL // 'DAY', PROGNAME )

        !.......  Check to see if appropriate variable list exists
        CALL RETRIEVE_IOAPI_HEADER( DAYNAME )

        I1 = INDEX1( 'ACRESBURNED', NVARS3D, VNAME3D )
        J1 = INDEX1( 'AREA', NVARS3D, VNAME3D )

        IF( I1 <= 0 .AND. J1 <= 0  ) THEN
            MESG = 'ERROR: Cannot find acres burned ' //    &
                   'variable "ACRESBURNED" or "AREA" in daily ' //    &
                    CRLF() // BLANK10 // 'inventory file '    &
                    // TRIM( DAYNAME )
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        DAY_NSRC = NROWS3D
        WRITE( MESG,94010 )'NOTE: Number of Sources in Daily File',    &
                           DAY_NSRC
        CALL M3MSG2( MESG )

        IF( EFLAG ) THEN
            MESG = 'Problem with hourly fire data inputs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


        ALLOCATE( DAY_ACRES( DAY_NSRC ),    &
                  DAY_INDEX( DAY_NSRC ),    &
                      ACRES( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'DAY_ACRES...ACRES', PROGNAME )
        ACRES = 0.0

        JDATE = SDATE
        JTIME = STIME
        DO T = 1, NSTEPS

            CALL SAFE_READ3( DAYNAME, 'ACRESBURNED', ALLAYS3, JDATE, JTIME, DAY_ACRES )       ! Wildfire inventory format

            IF ( .NOT. READ3( DAYNAME, 'INDXD', ALLAYS3, JDATE, JTIME, DAY_INDEX ) ) THEN

                MESG = 'Could not read "INDXD" from file "'// TRIM( DAYNAME ) // '".'
                CALL M3EXIT( PROGNAME, SDATE, STIME, MESG, 2 )

            END IF

            DO DS = 1, DAY_NSRC
                S = DAY_INDEX( DS )
                IF( S > 0 ) ACRES( S ) = DAY_ACRES( DS )
            ENDDO

            CALL NEXTIME( JDATE, JTIME, 10000 )

        ENDDO

    ENDIF

    !.......  Loop through sources to determine elevated and PinG sources.  If
    !           source is in a stack group, use group settings to compare to
    !           the elevated and/or PinG criteria.
    !.......  Also, update stack parameters and emissions permanently for
    !           program duration if source is in an inventory group
    NGROUP = 0
    PGRP   = -9
    PEGRP  = -9
    PLTEND = CH1POS3 - 1
    DO J = 1, NSRC

        S    = GINDEX ( J )
        IGRP = GROUPID( S )

        !.......  Exclude sources that are outside of the grid
        IF( .NOT. INGRID( SRCXL( S ), SRCYL( S ), NCOLS, NROWS, COL, ROW ) ) CYCLE

        !.......  For sources in an inventory group...
        IF ( IGRP .GT. 0 ) THEN

            !.......  Update stack parameters, if needed
            XLOCA( S ) = GRPLON( IGRP )
            YLOCA( S ) = GRPLAT( IGRP )
            STKDM( S ) = GRPDM ( IGRP )
            STKHT( S ) = GRPHT ( IGRP )
            STKTK( S ) = GRPTK ( IGRP )
            STKVE( S ) = GRPVE ( IGRP )
            CIFIP( S ) = GRPFIP( IGRP )

        END IF

        !.......  Store reordered group IDs
        SRCGROUP( S ) = IGRP

        !.......  Set temporary values for the current source
        CFIP = CIFIP  ( S )
        CSRC = CSOURC ( S )
        PLT  = CSRC   ( PLTPOS3:PLTEND )
        HT   = STKHT  ( S )
        DM   = STKDM  ( S )
        TK   = STKTK  ( S )
        VE   = STKVE  ( S )
        SCC  = CSCC   ( S )
        IF (SCC(9:9) .eq. 'S') SMOLDER(S) = .TRUE.

        !.......  Set values to be compared in selection formulas for this
        !         source, depending on whether source is in a group or not...
        !.......  Source arrays have already been updated with group info
        !         (include emissions TOTAL for group).
        VALS( HT_IDX ) = HT
        VALS( DM_IDX ) = DM
        VALS( TK_IDX ) = TK
        VALS( VE_IDX ) = VE
        DMVAL = DM
        VEVAL = VE
        IF( DM == BADVAL3 ) DMVAL = 0.0
        IF( VE == BADVAL3 ) VEVAL = 0.0
        VALS( FL_IDX ) = 0.25 * PI * DMVAL * DMVAL * VEVAL
        VALS( SRC_IDX )= S
            !VALS( FIP_IDX )= CFIP
        CHRS( PLT_IDX )= ADJUSTL( PLT )

        !.......  If cutoff approach is used, compute and store plume rise
        IF( LCUTOFF ) THEN

            !.......  Calculate estimated plume rise
            RISE( S ) = PLUMRIS( HT, TK, VE, DM )
            VALS( RISE_IDX ) = RISE( S )


        !.......  Otherwise, set value of rise to zero
        ELSE
            VALS( RISE_IDX ) = 0.

        END IF                    ! end cutoff approach or not

        !.......  Add pollutant value to VALS and set RANK for pollutants
        IF( NEVPEMV .GT. 0 ) THEN
            N = MAX( HT_IDX,DM_IDX,TK_IDX,VE_IDX,FL_IDX,RISE_IDX,    &
                     SRC_IDX, FIP_IDX, PLT_IDX  )         ! in case of code alteration
            DO K = 1, NEVPEMV
                N = N + 1
                VALS( N ) = MXEMIS( S,K )
                RANK( N ) = REAL( MXRANK( S,K ) )
            END DO
        END IF

        !.......  If PELVCONFIG used for elevated sources, check if source matches
        !               criteria given
        IF(  (ELEVTYPE .EQ. PELVCONFIG_APPROACH-1) ) THEN

        !.......  See if source matches criteria for elevated sources
            EVSTAT = .FALSE.          ! array
            IF ( EVALCRIT( NEVPVAR, NELVCRIT, MXELVCHK, VALS, VALS,     &
                           RANK, CHRS, ELVVALS, ELVCHRS, ELVTYPES,      &
                           EVSTAT ) ) THEN
                IF ( IGRP .NE. PGRP ) NMJRGRP = NMJRGRP + 1
                NMAJOR = NMAJOR + 1
                LMAJOR( S ) = .TRUE.

            END IF

        END IF                ! End elevated sources approach

        IF( ELEVTYPE .EQ. PELVCONFIG_APPROACH ) THEN

    !.......  See if source matches criteria for elevated sources
            EVSTAT = .FALSE.      ! array

            IF ( FFLAG .AND. SMOLDER( S ) ) THEN
                LMAJOR( S ) = .FALSE.

            ELSE
                IF ( EVALCRIT( NEVPVAR, NELVCRIT, MXELVCHK,VALS,VALS,   &
                            RANK, CHRS, ELVVALS, ELVCHRS, ELVTYPES,     &
                            EVSTAT ) ) THEN
                    NMAJOR = NMAJOR + 1
                    IF ( IGRP .NE. PGRP ) NMJRGRP = NMJRGRP + 1
                    LMAJOR( S ) = .TRUE.
                END IF
            END IF
        END IF                    ! End elevated sources approach


        !.......  If PELVCONFIG used for PinG sources, check if source matches
        !               criteria given
        IF( PINGTYPE .EQ. ( PELVCONFIG_APPROACH-1 ) ) THEN

            !.......  See if source matches criteria for PinG sources
            PGSTAT = .FALSE.              ! array

            IF ( FFLAG .AND. SMOLDER( S ) ) THEN
                LPING(S) = .FALSE.

            ELSE
                IF ( EVALCRIT( NEVPVAR, NPNGCRIT, MXPNGCHK,VALS,VALS,    &
                            RANK, CHRS, PNGVALS, PNGVALS, PNGTYPES,    &
                            PGSTAT )  ) THEN
                    NPING = NPING + 1
                    IF ( IGRP .NE. PGRP ) NPINGGRP = NPINGGRP + 1
                    LPING( S ) = .TRUE.

                END IF

            END IF
        END IF         ! End whether PinG approach is to use PELVCONFIG or not

        !.......  If source is a major source or a PinG source, but it's not in
        !               a group, increase the total maximum group count.
        IF( ( LMAJOR( S ) .OR.    &
              LPING ( S )      ) ) THEN

            IF ( IGRP .EQ. 0 .OR. IGRP .NE. PEGRP )  NGROUP = NGROUP + 1

            PEGRP = IGRP

        END IF

        PGRP = IGRP

    END DO             ! End loop over sources

    !.......  Assign a fake source when there is no group for elevated and/or ping
    IF( MFLAG .AND. NGROUP < 1 ) THEN
        RFLAG = .FALSE.
        NGROUP = 1
        DO J = 1, NSRC
            S    = GINDEX ( J )
            IGRP = GROUPID( S )

    !.......  Select a source that are inside of the grid
            IF( INGRID( SRCXL( S ), SRCYL( S ), NCOLS, NROWS, COL, ROW ) ) THEN

                IF( (ELEVTYPE .EQ. PELVCONFIG_APPROACH) .OR.    &
                    (ELEVTYPE .EQ. PELVCONFIG_APPROACH-1) ) THEN
                    NMAJOR = 1
                    LMAJOR( S ) = .TRUE.
                END IF

                IF( (PINGTYPE .EQ. PELVCONFIG_APPROACH-1) ) THEN
                    NPING = 1
                    LPING ( S ) = .TRUE.
                END IF

                EXIT       ! exit once fill a fake source

            END IF
        END DO
    END IF

    !.......  Now reset group arrays for major and PinG sources only. Groups
    !           are not used by SMOKE for other point sources.

    !.......  Allocate memory for and save inventory groups in local arrays
    ALLOCATE( LOCGID( NINVGRP ),    &
              LOCCNT( NINVGRP ),    &
             LOCSTAT( NINVGRP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LOCGID...LOCSTAT', PROGNAME )

    LOCGID = GRPGID
    LOCCNT = GRPCNT
    LOCSTAT = .FALSE.

    !.......  Deallocate inventory groups so that these can be allocated
    IF ( ALLOCATED( GRPLAT ) ) THEN

        DEALLOCATE( GRPLAT, GRPLON, GRPDM, GRPHT, GRPTK,    &
                    GRPVE, GRPFL, GRPCNT, GRPFIP )
        IF (FFLAG) DEALLOCATE (GRPACRES)

    END IF

    !.......  Reallocate group arrays based on inventory groups, major, and
    !           PinG settings.
    !.......  The group sorting index is here in case we need to add back
    !           in the reading of PSPLIT and PGROUP files, which might be
    !           unsorted.  The WPINGSTK routine uses this index
    ALLOCATE( GRPGIDA( NGROUP ),    &
               GRPIDX( NGROUP ),    &
               GRPLAT( NGROUP ),    &
               GRPLON( NGROUP ),    &
                GRPDM( NGROUP ),    &
                GRPHT( NGROUP ),    &
                GRPTK( NGROUP ),    &
                GRPVE( NGROUP ),    &
                GRPFL( NGROUP ),    &
               GRPCNT( NGROUP ),    &
               GRPCOL( NGROUP ),    &
               GRPROW( NGROUP ),    &
                GRPXL( NGROUP ),    &
                GRPYL( NGROUP ),    &
               GRPFIP( NGROUP ),    &
            GRPLMAJOR( NGROUP ),    &
             GRPLPING( NGROUP ),    &
                 INDX( NGROUP ),    &
                   GN( NGROUP ),    &
                   SN( NGROUP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GRPGIDA...SN', PROGNAME )
    IF (FFLAG) THEN
        ALLOCATE( GRPACRES(NGROUP), STAT=IOS )
        CALL CHECKMEM( IOS, 'GRPACRES', PROGNAME )
    ENDIF

    GRPGIDA = 0
    GRPIDX  = 0
    GRPLAT  = BADVAL3
    GRPLON  = BADVAL3
    GRPDM   = BADVAL3
    GRPHT   = BADVAL3
    GRPTK   = BADVAL3
    GRPVE   = BADVAL3
    GRPFL   = BADVAL3
    GRPFIP  = ' '
    GRPCNT  = 0
    GRPCOL  = 0
    GRPROW  = 0
    GRPXL   = BADVAL3
    GRPYL   = BADVAL3
    GN      = 0
    SN      = 0
    GRPLMAJOR = 0
    GRPLPING  = 0
    IF (FFLAG) GRPACRES = BADVAL3

    !.......  Loop over sources to fill in group settings with new group numbers
    !           and to populate group arrays for major and PinG sources.
    !           Groups that aren't major or PinG sources will be dropped.
    !.......  Make sure to keep sources together that share a group ID from
    !           the inventory grouping.
    !.......  Source arrays have already been updated with group info.
    !.......  OUTG is used because G is not necessarily going to stay in order
    !           with the LOCGID construct.
    G = 0
    CHRS = ' '          ! array
    DO S = 1, NSRC

        !.......  If major or PinG
        IF ( LMAJOR( S ) .OR. LPING( S ) ) THEN

            !.......  Retrieve old group number
            IGRP = SRCGROUP( S )

            !.......  If source is in an inventory group...
            IF ( IGRP .GT. 0 ) THEN

                !.......  If  it has already been assigned a new number,
                !                      retrieve number for the report.
                IF ( LOCSTAT( IGRP ) ) THEN

                    GROUPID( S ) = LOCGID( IGRP )
                    OUTG = LOCGID( IGRP )
                    SFLAG = .FALSE.                    ! controller for later in loop

                !.......  Otherwise, group information has not yet been stored and
                !                       it needs to be
                ELSE

                    !...........  If source is in an inventory group...
                    G = G + 1
                    OUTG = G
                    LOCGID ( IGRP ) = G                          ! store for next iteration
                    LOCSTAT( IGRP ) = .TRUE.
                    SFLAG = .TRUE.                    ! controller for later in loop

                    IF ( G .LE. NGROUP ) THEN
                        GROUPID( S ) = G
                        GRPIDX ( G ) = LOCGID( IGRP )
                        GRPGIDA( G ) = LOCGID( IGRP )
                        GRPCNT ( G ) = LOCCNT( IGRP )
                    END IF

                END IF

            !.......  If source not in an inventory group, it needs a grp no.
            ELSE
                G = G + 1
                OUTG = G
                SFLAG = .TRUE.                ! controller for later in loop

                IF ( G .LE. NGROUP ) THEN
                    GROUPID( S ) = G
                    GRPIDX ( G ) = G
                    GRPGIDA( G ) = G
                    GRPCNT ( G ) = 1
                END IF

            END IF

            !.......  Ensure no internal error overflow
            IF ( G .GT. NGROUP ) CYCLE

            !.......  Store the rest of the group settings in output arrays
            IF( SFLAG ) THEN
                GRPLAT( G ) = YLOCA ( S )
                GRPLON( G ) = XLOCA ( S )
                GRPDM ( G ) = STKDM ( S )
                GRPHT ( G ) = STKHT ( S )
                GRPTK ( G ) = STKTK ( S )
                GRPVE ( G ) = STKVE ( S )
                DM   = STKDM  ( S )
                VE   = STKVE  ( S )
                DMVAL = DM
                VEVAL = VE
                IF( GRPDM( G ) == BADVAL3 ) DMVAL = 0.0
                IF( GRPVE( G ) == BADVAL3 ) VEVAL = 0.0
                GRPFL ( G ) = 0.25 * PI * DMVAL * DMVAL * VEVAL
                GRPFIP( G ) = CIFIP ( S )
                IF (FFLAG) GRPACRES( G ) = ACRES( S)
                IF (LMAJOR(S)) GRPLMAJOR( G ) = 1
                IF (LPING(S)) GRPLPING ( G ) = 1
            END IF

            !.......  Write out report information...

            !.......  Get setup for another call to EVALCRIT to get STATUS
            VALS = 0.                       ! array
            VALS( HT_IDX ) = GRPHT ( OUTG )
            VALS( DM_IDX ) = GRPDM ( OUTG )
            VALS( TK_IDX ) = GRPTK ( OUTG )
            VALS( VE_IDX ) = GRPVE ( OUTG )
            VALS( FL_IDX ) = GRPFL ( OUTG )
            IF( LCUTOFF ) VALS( RISE_IDX ) = RISE( S )
            VALS( SRC_IDX )= S
                        !VALS( FIP_IDX )= CIFIP( S )

            PLT = CSOURC( S )( PLTPOS3:PLTEND )
            CHRS( PLT_IDX )= ADJUSTL( PLT )

            !.......  Add pollutant value to VALS and set RANK for pollutants
            IF( NEVPEMV .GT. 0 ) THEN
                N = MAX(HT_IDX,DM_IDX,TK_IDX,VE_IDX,FL_IDX,RISE_IDX,    &
                        SRC_IDX, FIP_IDX, PLT_IDX  )             ! in case of code alteration
                DO K = 1, NEVPEMV
                    N = N + 1
                    VALS( N ) = MXEMIS( S,K )
                    RANK( N ) = REAL( MXRANK( S,K ) )
                END DO
            END IF

            !.......  If source is PinG, write out for PinG
            IF ( RFLAG ) THEN
                IF ( LPING( S ) ) THEN

                    !....... Evaluate PinG criteria again to get PGSTAT for writing;
                    !                      if valid, then write report fields
                    IF ( EVALCRIT( NEVPVAR, NPNGCRIT, MXPNGCHK, VALS,    &
                                   VALS, RANK, CHRS, PNGVALS, PNGCHRS,    &
                                   PNGTYPES, PGSTAT ) ) THEN

                        CALL WRITE_REPORT( RDEV, S, OUTG, NEVPVAR,    &
                             NPNGCRIT, MXPNGCHK, 'P', VALS, RANK, CHRS,    &
                             PNGVALS, PNGCHRS, PNGTYPES, PGSTAT )

                    !.......  Otherwise, internal error
                    ELSE
                        EFLAG = .TRUE.

                        WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //    &
                               'Second evaluation of PinG source ', S,    &
                               'inconsistent with first evaluation.'
                        CALL M3MESG( MESG )

                    END IF

                !....... Evaluate elevated criteria again to get PGSTAT for
                !                      writing; if valid, then write report fields
                ELSE

                    IF ( EVALCRIT( NEVPVAR, NELVCRIT, MXELVCHK, VALS,    &
                                   VALS, RANK, CHRS, ELVVALS, ELVCHRS,    &
                                   ELVTYPES, EVSTAT )  ) THEN

                        !...........  Add source to report
                        CALL WRITE_REPORT( RDEV, S, OUTG, NEVPVAR,    &
                             NELVCRIT, MXELVCHK, 'E', VALS, RANK, CHRS,    &
                             ELVVALS, ELVCHRS, ELVTYPES, EVSTAT )

                    !.......  Otherwise, internal error
                    ELSE
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 ) 'INTERNAL ERROR: ' //    &
                               'Second evaluation of elevated source ',    &
                               S, 'inconsistent with first evaluation.'
                        CALL M3MESG( MESG )

                    END IF


                END IF
            END IF

        END IF          ! End if major or PinG sources

    END DO          ! End loop on sources

    !.......  Ensure that all is well with memory allocation
    IF ( G .NE. NGROUP ) THEN
        WRITE( MESG,94010 ) 'INTERNAL ERROR: Expected number of groups was',    &
                            NGROUP, 'but actual number was',  G
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  If all values are zero, give error
    IF ( NMAJOR .EQ. 0 .AND. NPING  .EQ. 0 ) THEN
        MESG = 'No groups, major sources, or plume-in-grid sources.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Process the stack group coordinates for the current grid...

    !.......  Convert x,y location to coordinates of the projected grid
    GRPXL = GRPLON
    GRPYL = GRPLAT
    CALL CONVRTXY( NGROUP, GDTYP, GRDNM, P_ALP, P_BET, P_GAM, XCENT, YCENT, GRPXL, GRPYL )

    !.......  Determine grid cells for these coordinate locations
    IF( VFLAG ) THEN
        CALL GENPTVCEL( NGROUP, NGRID, GRPXL, GRPYL, NEXCLD, NX, INDX, GN, SN )
    ELSE

        CALL GENPTCEL( NGROUP, NGRID, GRPXL, GRPYL, NEXCLD, NX, INDX, GN, SN )
    END IF

    !.......  Convert grid cells to row and columns numbers
    DO I = 1, NGROUP

        ROW = 0
        COL = 0
        N   = GN( I )

        IF( N .GT. 0 ) THEN
            ROW = N / NCOLS              ! note: integer math
            IF( MOD( N, NCOLS ) .GT. 0. ) ROW = ROW + 1
            COL = N - ( ROW-1 ) * NCOLS
        END IF

        GRPROW( I ) = ROW
        GRPCOL( I ) = COL

    END DO

    !.......  Abort if an error occurred
    IF( EFLAG ) THEN
        MESG = 'Problem selecting major/plume-in-grid sources'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Write status of processing
    IF( NGROUP .GT. 0 ) THEN
        WRITE( MESG,94010 ) 'Number of stack groups:', NGROUP
        CALL M3MSG2( MESG )
    END IF

    IF( NMAJOR .GT. 0 ) THEN
        WRITE( MESG,94010 ) 'Number of major sources:', NMAJOR
        CALL M3MSG2( MESG )
    END IF

    IF( NMJRGRP .GT. 0 ) THEN
        WRITE(MESG,94010) 'Number of major source groups:', NMJRGRP
        CALL M3MSG2( MESG )
    END IF

    IF( NPING .GT. 0 ) THEN
        WRITE( MESG,94010 ) 'Number of plume-in-grid sources:',NPING
        CALL M3MSG2( MESG )
    END IF

    IF( NPINGGRP .GT. 0 ) THEN
        WRITE( MESG,94010 ) 'Number of plume-in-grid source groups:', NPINGGRP
        CALL M3MSG2( MESG )
    END IF

    !.......  Open output files
    CALL OPENEOUT( NGROUP, SDATE, STIME, ENAME, VFLAG, LFLAG,    &
                   PDEV, MNAME )

    !.......  Write ASCII file
    MESG = 'Writing ELEVATED POINT SOURCE output file...'
    CALL M3MSG2( MESG )

    COORUN3D = 'METERS '
    IF ( GDTYP3D .EQ. LATGRD3 ) THEN
        COORD3D = 'LAT-LON '
        COORUN3D = 'DEGREES '
    ELSE IF ( GDTYP3D .EQ. LAMGRD3 ) THEN
        COORD3D = 'LAMBERT '
    ELSE IF ( GDTYP3D .EQ. MERGRD3 ) THEN
        COORD3D = 'MERCATOR '
    ELSE IF ( GDTYP3D .EQ. EQMGRD3 ) THEN
        COORD3D = 'MERCATOR '
    ELSE IF ( GDTYP3D .EQ. STEGRD3 ) THEN
        COORD3D = 'STEREOGRAPHIC '
    ELSE IF ( GDTYP3D .EQ. UTMGRD3 ) THEN
        COORD3D = 'UTM '
    ELSE IF ( GDTYP3D .EQ. POLGRD3 ) THEN
        COORD3D = 'POLAR '
    ELSE
        MESG = 'Current projection code is not supported'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    WRITE( PDEV, 94500 ) GDNAM3D, XORIG3D, YORIG3D, XCELL3D,    &
                  YCELL3D, NCOLS, NROWS, NTHIK3D, COORD3D,      &
                  COORUN3D, P_ALP, P_BET, P_GAM, XCENT, YCENT

    DO S = 1, NSRC

        IF( LMAJOR( S ) .OR. LPING( S ) ) THEN

            MS   = 0
            PS   = 0
            IF( LMAJOR( S ) ) MS = S
            IF( LPING ( S ) ) THEN
                MS = 0
                PS = S
            END IF
            IGRP = GROUPID( S )

            IF( LPING( S ) .AND. ( NEVPEMV .GT. 0 ) ) THEN
                CSRC = CSOURC( S )
                PLT = CSRC( PTBEGL3( 2 ):PTENDL3( 2 ) )
                STK = CSRC( PTBEGL3( JSTACK ):PTENDL3( JSTACK ) )
                WRITE( PDEV, 93630 ) MS, PS, IGRP, CIFIP( S ),    &
                    PLT, STK, MXEMIS( S,1 )
            ELSE
                WRITE( PDEV, 93620 ) MS, PS, IGRP
            END IF

        END IF

    END DO

    !.......  Sort and write plume-in-grid output file for Models-3 processing
    !           or for elevated source identidied by cutoff method
    IF( NGROUP .GT. 0 ) THEN

    !.......  Make sure that stack parameters are set for all groups
    !.......  This is a simplistic way of doing this for now,  later
    !               add call to FIXSTK routine
        MINDM = MINVAL( GRPDM( 1:NGROUP ) )
        MINHT = MINVAL( GRPHT( 1:NGROUP ) )
        MINTK = MINVAL( GRPTK( 1:NGROUP ) )
        MINVE = MINVAL( GRPVE( 1:NGROUP ) )
        MINFL = MINVAL( GRPFL( 1:NGROUP ) )

        IF ( .NOT. FFLAG ) THEN
            IF((MIN( MINDM, MINHT, MINTK, MINVE, MINFL ) .LT. 0.) ) THEN

                MESG = 'Bad stack group or stack split file. ' //       &
                       'Unable to assign stack ' // CRLF()//BLANK10//   &
                       'parameters to all stack groups. Could be '//    &
                       'a source matching problem.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF
        END IF

    !.......  Give warning if any plume-in-grid stack groups are outside the
    !               grid
        IF( NEXCLD .GT. 0 ) THEN
            WRITE( MESG,94010 ) 'WARNING: ', NEXCLD,     &
                   'stack groups are outside of grid "' // TRIM( GDNAM3D ) // '"'
            CALL M3MSG2( MESG )
        END IF

        MESG='Writing ELEVATED/PING STACK PARAMETERS output file...'
        CALL M3MSG2( MESG )

        CALL WPINGSTK( MNAME, SDATE, STIME, LFLAG )

    END IF

    !.......  Normal completion of program
    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93500 FORMAT( I6, A1, 21X, I5, F9.0, F9.0, 3X, F8.0, F7.0, F7.0,    &
            F7.0, F10.0 )

93550 FORMAT( 6X, I6, A1, A1, I1, I2, I3, A15, A15, A11, 7X,    &
            F9.0, F9.0, F8.0, F7.0, F7.0, F7.0, F10.0 )

93620 FORMAT( 3(I8,1X) )

93630 FORMAT( 3(I8,1X), A, 1X, 2(A20,1X), F10.3 )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

94020 FORMAT( A, 1X, F8.2, 1X, A )

94030 FORMAT( 'H[m]:', 1X, F6.2, 1X, 'D[m]:'  , 1X, F4.2, 1X,    &
            'T[K]:', 1X, F7.1, 1X, 'V[m/s]:', 1X, F10.1 )

94300 FORMAT( A, I2.2, A, I2.2, A )

94500 FORMAT( '#GRID ', A, 4F15.4, 3I5, 2(2X,A), 5F10.4 )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram writes one report record based
    !               on interpretation of the STATUS argument using the
    !               contents of the other array arguments.
    SUBROUTINE WRITE_REPORT( FDEV, S, IGRP, NV, NORS, MXAND,    &
                             LABEL, VALS, RANK, CHRS, COMPARE,  &
                             COMCHRS, TYPES, STATUS )

        !.......  Subprogram arguments
        INTEGER     ,  INTENT(IN):: FDEV            ! report file unit number
        INTEGER     ,  INTENT(IN):: S               ! source number
        INTEGER     ,  INTENT(IN):: IGRP            ! group number
        INTEGER     ,  INTENT(IN):: NV              ! Number of values
        INTEGER     ,  INTENT(IN):: NORS            ! Number of OR conditions
        INTEGER     ,  INTENT(IN):: MXAND           ! Max no.  ANDs for single data val
        CHARACTER(*),  INTENT(IN):: LABEL           ! E=elevated,  P=PinG
        REAL        ,  INTENT(IN):: VALS   ( NV )               ! Data values
        REAL        ,  INTENT(IN):: RANK   ( NV )               ! Ranking order
        CHARACTER(*),  INTENT(IN):: CHRS   ( NV )               ! String values
        REAL        ,  INTENT(IN):: COMPARE( NORS,  MXAND,  NV )         ! Formula values
        CHARACTER(*),  INTENT(IN):: COMCHRS( NORS,  MXAND,  NV )         ! Condition
        CHARACTER(*),  INTENT(IN):: TYPES  ( NORS,  MXAND,  NV )         ! Condition
        LOGICAL     ,  INTENT(IN):: STATUS ( NORS,  MXAND,  NV )         ! true: condition met

        INTEGER,  PARAMETER :: NHEADER  = 17
        CHARACTER(15),  PARAMETER :: HEADERS( NHEADER ) =   &
                                (/ 'Source ID      ',       &
                                    'Region         ',      &
                                    'Plant          ',      &
                                    'Char 1         ',      &
                                    'Char 2         ',      &
                                    'Char 3         ',      &
                                    'Char 4         ',      &
                                    'Char 5         ',      &
                                    'NAICS          ',      &
                                    'Plt Name       ',      &
                                    'Elevstat       ',      &
                                    'Group          ',      &
                                    'Stk Ht         ',      &
                                    'Stk Dm         ',      &
                                    'Stk Tmp        ',      &
                                    'Stk Vel        ',      &
                                    'Stk Flw        '       /)

        !.......  Subprogram local allocatable arrays
        LOGICAL,             ALLOCATABLE,  SAVE ::    LF( : )          ! true: output source chars
        CHARACTER(NAMLEN3),  ALLOCATABLE,  SAVE :: VNAME( : )          ! var names

        !.......  Subprogram local static arrays
        CHARACTER(32) CHARS ( MXCHRS )         !  source fields for output

        !.......  Local subprogram variables
        INTEGER      K,  L,  L1,  L2,  M,  N          ! indices and counters
        INTEGER      MX                         ! max of MXPNGCHK and MXELVCHK
        INTEGER      NC                         ! local no. of src chars to output
        INTEGER   :: NM = 0                     ! local no. of max vals before emis

        LOGICAL      DFLAG                      ! true: OR was true
        LOGICAL   :: FIRSTIME = .TRUE.          ! true: first time subprogram called

        CHARACTER(1024) BUFFER
        CHARACTER(256)  FMTBUF

        !----------------------------------------------------------------------

        !.......  If firstime routine called
        IF ( FIRSTIME ) THEN

            !.......  Allocate local memory
            ALLOCATE( LF( MXCHRS ),     &
               VNAME( NV ),  STAT=IOS )
            CALL CHECKMEM( IOS,  'VNAME',  PROGNAME )

            !.......  Initialize output status of source characteristics
            LF = .FALSE.                ! array
            LF( 1:NCHARS ) = .TRUE.

            !.......  Initialize variable names for report
            VNAME( HT_IDX )   = 'HT'
            VNAME( DM_IDX )   = 'DM'
            VNAME( TK_IDX )   = 'TK'
            VNAME( VE_IDX )   = 'VE'
            VNAME( FL_IDX )   = 'FL'
            VNAME( RISE_IDX ) = 'RISE'
            VNAME( SRC_IDX )  = 'SOURCE'
            VNAME( FIP_IDX )  = 'FIPS'
            VNAME( PLT_IDX )  = 'PLANT'

            NM = MAX( HT_IDX, DM_IDX, TK_IDX, VE_IDX, FL_IDX, RISE_IDX,     &
                      SRC_IDX,  FIP_IDX,  PLT_IDX  )             ! in case of code alteration
            N = NM
            DO K = 1,  NEVPEMV
                N = N + 1
                VNAME( N ) = EINAM( EVPEMIDX( K ) )
            END DO

            !.......  Build header (this is sloppy job for now)
            N = 1
            M = 0
            BUFFER = HEADERS( 1 )
            DO WHILE( N < NHEADER )

                N = N + 1

                IF ( N .GE. 4 .AND. N .LE. 8 ) THEN
                    M = M + 1
                    IF ( M .GT. NCHARS-2 ) CYCLE
                END IF

                BUFFER = TRIM( BUFFER ) //'; ' // HEADERS( N )

            END DO

            !.......  Add plume rise onto header
            IF( LCUTOFF ) THEN
                BUFFER = TRIM( BUFFER ) // '; Rise '
            END IF

            !.......  Add pollutants onto header
            DO K = 1,  NEVPEMV
                BUFFER = TRIM( BUFFER ) // '; Group ' // EINAM( EVPEMIDX( K ) )
            END DO

            !.......  Add results onto header
            MX = MAX( MXELVCHK,  MXPNGCHK )
            DO K = 1,  MX
                WRITE( BUFFER,  '(4(A, I1))' ) TRIM( BUFFER ) //    &
                       '; Var',  K,  '; Type',  K,  '; Test',  K,   &
                       '; Val',  K
                IF( K .LT. MX ) THEN
                    WRITE( BUFFER,  '(A)' ) TRIM( BUFFER ) // '; AND'
                END IF

            END DO

            !.......  Write out header
            WRITE( FDEV,  '(A)' ) TRIM( BUFFER )

            FIRSTIME = .FALSE.

        END IF

        !.......  Subdivide source description
        CALL PARSCSRC( CSOURC( S ),  MXCHRS,  SC_BEGP,  SC_ENDP, LF,  NC,  CHARS )

        !.......  Write source information format and then use format
        WRITE( FMTBUF,  94790 ) FIPLEN3,  PLTLEN3, ( CHRLEN3,  N=1, NCHARS-2 ),  NAILEN3,  DSCLEN3
        FMTBUF = TRIM( FMTBUF ) // ')'

        WRITE( BUFFER,  FMTBUF ) S,  ( CHARS( N ),  N = 1,  NCHARS ), CNAICS( S ),  CPDESC( S )

        !.......  Add label,  group number,  stack parameters,  and emissions
        WRITE( BUFFER,  94791 ) TRIM( BUFFER ),  LABEL,  IGRP,      &
               VALS( HT_IDX ),  VALS( DM_IDX ),  VALS( TK_IDX ),    &
               VALS( VE_IDX ),  VALS( FL_IDX )

        !.......  If needed,  add plume rise value
        IF ( LCUTOFF ) THEN
            WRITE( BUFFER,  94792 ) TRIM( BUFFER ),  RISE( S )
        END IF

        !.......  Add emissions
        IF ( NEVPEMV .GT. 0 ) THEN
            WRITE( BUFFER,  94792 ) TRIM( BUFFER ), ( VALS( K ),  K = NM+1,  NV )
        END IF

        !.......  Write characteristics that caused matching
        DFLAG = .FALSE.
        DO L = 1,  NORS
            DO M = 1,  MXAND
                DO N = 1,  NV

                    !...........  Check if the status was used to include source
                    !                           or not
                    IF ( STATUS( L, M, N ) ) THEN

                        !...........  Exit after this OR
                        DFLAG = .TRUE.

                        !...........  Add to report for this OR and AND (if any)
                        IF ( TYPES( L, M, N ) .EQ. 'TOP' ) THEN
                            WRITE( BUFFER,  94793 ) TRIM( BUFFER ),     &
                                   VNAME( N ),  ' RANK;      =;', INT( RANK( N ) )

                        !...........  For integer values stored as reals
                        ELSE IF ( N .EQ. SRC_IDX .OR.    &
                                  N .EQ. FIP_IDX      ) THEN
                            WRITE( BUFFER,  94796 ) TRIM( BUFFER ),     &
                                   VNAME( N ),  TYPES( L, M, N ), INT( COMPARE( L, M, N ) )

                        !...........  Use "IS" type as way to I.D. string criteria
                        ELSE IF ( TYPES( L, M, N ) .EQ. 'IS' ) THEN
                            WRITE( BUFFER,  94795 ) TRIM( BUFFER ),     &
                                   VNAME( N ),  TYPES( L, M, N ), COMCHRS( L, M, N )

                        !...........  For all reals
                        ELSE IF ( TYPES( L, M, N ) .NE. ' ' ) THEN
                            WRITE( BUFFER,  94794 ) TRIM( BUFFER ),     &
                                   VNAME( N ),  TYPES( L, M, N ), COMPARE( L, M, N )

                        END IF

                    END IF

                END DO

                !.......  If this OR is valid and more than one AND,  add AND
                !         to output buffer
                IF( DFLAG .AND. M+1 .LE. MXAND ) THEN

                    DO N = 1,  NV
                        IF ( STATUS( L, M+1, N ) ) THEN

                            BUFFER = TRIM( BUFFER ) // ' AND;'
                            EXIT                  ! end loop

                        END IF
                    END DO

                END IF

            END DO

            !.......  Exit from "OR" loop if one of the OR criteria has been met
            IF( DFLAG ) EXIT

        END DO

        !.......  Write buffer to report file
        WRITE( FDEV,  '(A)' ) TRIM( BUFFER )

        !---------------------  FORMAT  STATEMENTS  -------------------------

94790   FORMAT( '(I7, ";", A',  I2.2,  ', ";"',  10(', A',  I2.2, ', ";"') )

94791   FORMAT( A,  1X,  A1,  '; ',  I6,  ';',  5( F10.2,  ';' ) )

94792   FORMAT( A,  1X,  20( F10.2,  ';' ) )

94793   FORMAT( A,  1X,  A16,  ';',  A,  I10,  ';' )

94794   FORMAT( A,  1X,  A16,  ';     ;',  1X,  A6,  '; ',  F10.2,  ';' )

94795   FORMAT( A,  1X,  A16,  ';     ;',  1X,  A6,  '; ',  A15,  ';' )

94796   FORMAT( A, 1X, A16, ';     ;', 1X, A6, '; ', I10, ';' )

    END SUBROUTINE WRITE_REPORT


    !----------------------------------------------------------------------


    SUBROUTINE RETRIEVE_IOAPI_HEADER( FILNAM )

        !.......  Subprogram arguments
        CHARACTER(*) FILNAM

        IF ( .NOT. DESC3( FILNAM ) ) THEN

            MESG = 'Could not get description of file "' //    &
                   FILNAM( 1:LEN_TRIM( FILNAM ) ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END SUBROUTINE RETRIEVE_IOAPI_HEADER


    !----------------------------------------------------------------------


    SUBROUTINE SAFE_READ3( FILNAM, VARNAM, LAYER, JDATE, JTIME, XBUF     )

        !.......  Subprogram arguments
        CHARACTER(*), INTENT(IN   ) :: FILNAM            ! logical file name
        CHARACTER(*), INTENT(IN   ) :: VARNAM        ! variable name
        INTEGER     , INTENT(IN   ) :: LAYER         ! layer number (or ALLAYS3)
        INTEGER     , INTENT(IN   ) :: JDATE         ! Julian date
        INTEGER     , INTENT(IN   ) :: JTIME         ! time
        REAL        , INTENT(  OUT) :: XBUF( * )     ! read buffer

        IF ( .NOT. READ3( FILNAM, VARNAM, LAYER, JDATE, JTIME, XBUF ) ) THEN

            IF( VARNAM == 'TEMP2' .OR. VARNAM == 'TEMP1P5' ) THEN
                MESG = 'Please reset PLUME_GTEMP_NAME to match ' //    &
                   'to a variable name from file '// TRIM( FILNAM )
                CALL M3MSG2( MESG )
            END IF

            MESG = 'Could not read "' // TRIM( VARNAM ) //    &
                   '" from file "' // TRIM( FILNAM ) // '".'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )


        END IF

    END SUBROUTINE SAFE_READ3

END PROGRAM ELEVPOINT
