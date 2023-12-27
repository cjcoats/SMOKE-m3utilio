
SUBROUTINE TMPBEIS360( CVSW )

    !***********************************************************************
    !  program body starts at line  187
    !
    !  DESCRIPTION:
    !        Computes hourly time stepped gridded biogenic emissions using
    !        normalized gridded emissions from Normbeis3 (3.12) and postprocessed MM5
    !        meteorology.
    !
    !  PRECONDITIONS REQUIRED:
    !        Postprocessed MM5 meteorology that contains temperature,
    !        solar radiation, and pressure data.
    !        Normalized gridded emissions B3GRD from Normbeis3 (3.12)
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !        HRBIO, PREBMET
    !
    !  REVISION  HISTORY:
    !       3/01: Prototype by Jeff Vukovich
    !        Tested only on 36km Lambert domain
    !        Summer/winter switch file option not tested
    !       8/04: Updated for BEIS v3.12
    !       7/14:  Updated with J. Bash's Temperature code
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !        System
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

    !......  Modules for public variables
    !......  This module contains the speciation profile tables
    USE MODSPRO, ONLY: MXSPEC, SPCNAMES

    !......  This module contains BEIS3 specific arrays
    USE MODBEIS3, ONLY: AVGEMIS, AVGLAI, NOEMIS

    IMPLICIT NONE

    !......  INCLUDES:
    INCLUDE 'EMCNST3.h90'         !
    INCLUDE 'B3V14DIMS3.h90'      ! biogenic-related constants

    !......  EXTERNAL FUNCTIONS and their descriptions:
    CHARACTER(50), EXTERNAL :: GETCFDSC
    INTEGER      , EXTERNAL :: GETFLINE
    CHARACTER(16), EXTERNAL :: VERCHAR

    !......  ARGUMENTS and their descriptions
    CHARACTER(50), INTENT(IN) :: CVSW        ! CVS release tag

    !......  PARAMETERs
    CHARACTER(16), PARAMETER :: PROGNAME = 'TMPBEIS312'       !  program name

    CHARACTER(80), PARAMETER :: PARMENU( 1 ) = (/ 'Use MM5 generated radiation' /)

    CHARACTER(16), PARAMETER :: UNITSMENU( 2 ) = (/ 'moles/hr ', 'moles/s  '  /)

    !......  Latitude and longitude for zenith angle calculation
    REAL, ALLOCATABLE :: LAT  ( :, : )        !  grid lat (deg) -90 <= LAT <= 90
    REAL, ALLOCATABLE :: LON  ( :, : )        !  grid lon (deg) -180 <= LON <= 180

    !......  Gridded meteorology data
    REAL, ALLOCATABLE :: TASFC ( :, : )         ! level-1 air  temperature (K)
    REAL, ALLOCATABLE :: COSZEN( :, : )         ! cosine of zenith angle
    REAL, ALLOCATABLE :: PRES  ( :, : )         ! surface pressure
    REAL, ALLOCATABLE :: ISLTYP( :, : )         ! soil type
    REAL, ALLOCATABLE :: SOILM ( :, : )         ! soil moisture
    REAL, ALLOCATABLE :: SOILT ( :, : )         ! soil temperature
    REAL, ALLOCATABLE :: RN    ( :, : )         ! nonconvective rainfall
    REAL, ALLOCATABLE :: RC    ( :, : )         ! convective rainfall
    REAL, ALLOCATABLE :: RAINFALL( :, :, : )    ! rainfall for 24 hours
    REAL, ALLOCATABLE :: RADYNI( :, : )         ! inverse of aerodynamic resistance
    REAL, ALLOCATABLE ::   RATM( :, : )         ! aerodynamic resistance
    REAL, ALLOCATABLE ::  TEMPG( :, : )         !ground temperature
    REAL, ALLOCATABLE ::    LAI( :, : )         ! LAI from WRF
    REAL, ALLOCATABLE ::     Q2( :, : )         ! two meter water vapor mixing ratio
    REAL, ALLOCATABLE ::  USTAR( :, : )         ! surface friction velocity
    REAL, ALLOCATABLE :: RSTOMI( :, : )         ! inverse of stomatal resistance
    REAL, ALLOCATABLE ::  RSTOM( :, : )         ! stomatal resistance
    REAL, ALLOCATABLE ::  RGRND( :, : )         ! solar radiation reaching ground

    INTEGER, ALLOCATABLE ::    SWITCH( :, : )         ! Seasonal switch
    INTEGER, ALLOCATABLE ::    PTYPE ( :, : )         ! NO emissions 'pulse type'
    INTEGER, ALLOCATABLE :: PULSEDATE( :, : )      ! date when NO emission pulse begins
    INTEGER, ALLOCATABLE :: PULSETIME( :, : )      ! time when NO emission pulse begins

    !......  Gridded normalized emissions
    REAL, ALLOCATABLE :: SEMIS( :, :, : )          ! temporary emissions
    REAL, ALLOCATABLE :: SLAI ( :, :, : )          ! temporary LAI
    REAL, ALLOCATABLE :: NONAGNO  ( :, : )         ! non agriculture NO emis
    REAL, ALLOCATABLE :: NGROWAGNO( :, : )         ! non growing season ag NO emis
    REAL, ALLOCATABLE :: GROWAGNO ( :, : )         ! growing season ag NO emis

    !......  Mole and mass factors
    REAL, ALLOCATABLE :: MLFAC( :, : )             ! mole factors
    REAL, ALLOCATABLE :: MSFAC( :, : )             ! mass factors

    !......  BEIS3 internal, output species
    REAL, ALLOCATABLE :: EMPOL( :, :, : )          ! emissions of biogenic categories
    REAL, ALLOCATABLE :: EMISL( :, :, : )          ! emissions in moles/hour
    REAL, ALLOCATABLE :: EMISS( :, :, : )          ! emissions in tons/hour

    CHARACTER(16), ALLOCATABLE :: EMSPC( : )       ! names of emitting species

    !......  Logical names and unit numbers
    INTEGER         RDEV        !  unit number for speciation profiles file

    CHARACTER(16)   LNAME       !  logical name for emissions output (moles)
    CHARACTER(16)   SNAME       !  logical name for emissions output (mass)
    CHARACTER(16)   NNAME       !  logical name for normalized-emissions input
    CHARACTER(16)   GNAME       !  logical name for GRID_CRO_2D
    CHARACTER(16)   BNAME       !  logical name for frost switch input
    CHARACTER(16)   M3NAME      !  logical name for MET_FILE1
    CHARACTER(16)   M2NAME      !  logical name for MET_FILE2
    CHARACTER(16)   PNAME       !  logical name for file with pressure variable
    CHARACTER(16)   SOILINP     !  logical name for input NO soil information
    CHARACTER(16)   SOILOUT     !  logical name for output NO soil information

    !......  Other variables and their descriptions:
    INTEGER         B, M        !  counters for biogenic, model species
    INTEGER         BSTEPS      !  no. of hourly time steps for output
    INTEGER         HR          !  current simulation hour
    INTEGER         I, J, K, L, N, C, R      !  loop counters and subscripts
    INTEGER         INDEX
    INTEGER         IOS         !  temporay IO status
    INTEGER         JDATE       !  current simulation date (YYYYDDD)
    INTEGER         JRUNLEN     !  run length from environment
    INTEGER         JTIME       !  current simulation time (HHMMSS)
    INTEGER         LDATE       !  previous simulation date
    INTEGER         MDATE       !  met file 1 start date
    INTEGER         MSPCS       !  no. of emitting species
    INTEGER         MTIME       !  met file 1 start time
    INTEGER         MXSTEPS     !  maximum number of time steps
    INTEGER         NCOLS       !  no. of grid columns
    INTEGER         NDATE       !  date of SOILINP file
    INTEGER         NGRID       !  no. of grid cells
    INTEGER         NLINES      !  no. of lines in GSPRO speciation profiles file
    INTEGER         NROWS       !  no. of grid rows
    INTEGER         NSTEPS      !  duration of met file
    INTEGER         NTIME       !  time of SOILINP file
    INTEGER         PARTYPE     !  method number to calculate PAR
    INTEGER         RDATE       !  met file 2 start date
    INTEGER         RTIME       !  met file 2 start time
    INTEGER         TSTEP       !  time step set by environment
    INTEGER         TZONE       !  output-file time zone ; not used in program
    INTEGER         UNITTYPE    ! define output units
    INTEGER ::      RHOURS = 24      ! no. of rainfall hours
    INTEGER ::      RVARS       ! number of valid rainfall hours (if less than one day run)
    INTEGER ::      RINDEX      ! pointer to RAINFALL BUFFER
    LOGICAL         EFLAG       !  error flag
    LOGICAL ::      SAMEFILE = .TRUE.       ! radiation/cld and tmpr data in same file
    LOGICAL ::      SWITCH_FILE = .TRUE.      ! use frost switch file
    LOGICAL ::      ASSUME_SUMMER = .TRUE.     ! use summer normalized emissions
    LOGICAL         GETATN
    LOGICAL         PX_VERSION         ! true: using PX version of MCIP
    LOGICAL         INITIAL_RUN, INITIAL_HOUR
    LOGICAL         VFLAG             ! true: use variable grid
    LOGICAL ::      FDFLAG= .FALSE.        ! true: output file meta desc

    CHARACTER(8)       CHKNAME        ! check for virtual mode
    CHARACTER(50)      CLOUDSHM       ! cloud scheme name
    CHARACTER(5)       CTZONE         ! time zone
    CHARACTER(256)     EQNAME         ! equivalent filename
    CHARACTER(50)      LUSE           ! land use description
    CHARACTER(300)     MESG           ! message buffer
    CHARACTER(50)      METSCEN        ! met scenario name
    CHARACTER(128)     METADESC       ! output meta description
    CHARACTER(16)      MRGFDESC       ! name for met data description
    CHARACTER(16)      PRESNAM        ! surface pressure variable name
    CHARACTER(16)      RCNAM          ! convective rainfall variable name
    CHARACTER(16)      RNNAM          ! nonconvective rainfall variable name
    CHARACTER(16)      SMNAM          ! soil moisture variable name
    CHARACTER(SPNLEN3) SPPRO          ! speciation profile to use
    CHARACTER(16)      STMPNAM        ! soil temperature variable name
    CHARACTER(16)      STNAM          ! soil type variable name
    CHARACTER(16)      TMPRNAM        ! temperature variable name
    CHARACTER(16)      RADYNINAM        ! inverse aero resistance variable name
    CHARACTER(16)      TEMPGNAM        ! ground temp variable name
    CHARACTER(16)      RGRNDNAM        ! radiation reaching ground variable name
    CHARACTER(16)      LAINAM        ! LAI (from WRF) variable name
    CHARACTER(16)      Q2NAM        ! 2m mixing ratio variable name
    CHARACTER(16)      RSTOMINAM        ! inverse stomatal resistance variable name
    CHARACTER(16)      USTARNAM         ! surface friction variable name
    CHARACTER(16)      VTMP           ! temporary variable name

    !***********************************************************************
    !   begin body of subroutine TMPBEIS312

    !......  Evaluate the environment variables...

    !......  Check if processing variable grid data
    MESG = 'Use variable grid definition'
    VFLAG = ENVYN( 'USE_VARIABLE_GRID', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

    !......  Get the time zone for output of the emissions
    TZONE = ENVINT( 'OUTZONE', 'Output time zone', 0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "OUTZONE"', 2 )
    END IF

    !......  Write time zone to character string
    WRITE( CTZONE,94000 ) TZONE

    !......  Get the speciation profile to use
    MESG = 'Speciation profile to use for biogenics'
    CALL ENVSTR( 'BIOG_SPRO', MESG, 'B3V10', SPPRO, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BIOG_SPRO"', 2 )
    END IF

    !......  Are the radiation/cloud data in the same file as temperature data
    MESG = 'Radiation/cloud in same file as temperature data?'
    SAMEFILE = ENVYN ( 'BIOMET_SAME', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BIOMET_SAME"', 2 )
    END IF

    !......  Determine output units; only moles/hr and moles/s supported
    UNITTYPE = ENVINT( 'OUT_UNITS', 'BEIS3 output units ', 1, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "OUT_UNITS"', 2 )
    END IF

    IF( UNITTYPE /= 1 .AND. UNITTYPE /= 2 ) THEN
        MESG = 'ERROR: Unsupported output units setting; valid ' //    &
               'settings are 1 (moles/hr) and 2 (moles/s)'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MESG = 'BEIS3 output will be in ' // UNITSMENU( UNITTYPE )
    CALL M3MESG( MESG )

    !......  Open speciation profiles file
    RDEV = PROMPTFFILE( 'Enter logical name for SPECIATION PROFILES file',    &
                        .TRUE., .TRUE., 'GSPRO', PROGNAME )

    !......  Scan speciation profiles file to get all of the pollutant-species
    !        combinations that are valid for the pollutants in the inventory.
    !......  The species names are sorted in ABC order for each pollutant, and
    !        and the pollutants are in the same order as BIOTYPES.
    !......  Also retrieve the maximum number of species per pollutant and
    !        maximum number of profile entries per pollutant.

    CALL DSCSPROF( RDEV, NSEF, BIOTYPES )


    NLINES = GETFLINE( RDEV, 'Speciation profile file' )

    ALLOCATE( EMSPC( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMSPC', PROGNAME )

    MSPCS = 0


    !......  Find emitting species names
    DO I = 1, NSEF
        DO J = 1, MXSPEC
            IF ( SPCNAMES( J, I ) .EQ. ' ' ) CYCLE

            K = INDEX1 ( SPCNAMES( J, I ), MSPCS, EMSPC )

            IF ( K .EQ. 0 ) THEN
                MSPCS = MSPCS + 1
                EMSPC( MSPCS ) = SPCNAMES( J, I )
            END IF
        END DO
    END DO

    !......  Allocate memory for storing mole- and mass-based factors
    ALLOCATE( MLFAC ( MSPCS, NSEF ),    &
              MSFAC ( MSPCS, NSEF ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MLFAC,MSFAC', PROGNAME )

    !......  Read speciation profiles file
    MESG = 'Reading biogenic speciation profile...'
    CALL M3MSG2( MESG )

    CALL RDBPRO( RDEV, SPPRO, NSEF, BIOTYPES, MSPCS, EMSPC,    &
                 MLFAC, MSFAC )

    !......  Get the method to calculate PAR
    PARTYPE = ENVINT( 'BG_CLOUD_TYPE',    &
                      'How PAR will be calculated', 1, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BG_CLOUD_TYPE"', 2 )
    ELSE IF( PARTYPE /= 1 ) THEN
        MESG = 'ERROR: Unsupported PAR calculation method; ' //    &
               'valid method is 1 (Use MM5 generated radiation)'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MESG =  'PAR calculation will ' // PARMENU( PARTYPE )
    CALL M3MESG( MESG )

    !......  Check to see if frost date switch file to be used
    MESG = 'Using a frost date switch file?'
    SWITCH_FILE = ENVYN ( 'BIOSW_YN', MESG, .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BIOSW_YN"', 2 )
    END IF

    !......  Get normalized emissions file, BGRD
    NNAME = PROMPTMFILE(    &
            'Enter name for NORMALIZED EMISSIONS input file',    &
            FSREAD3, 'B3GRD', PROGNAME )

    !......  Read description of normalized emissions file
    IF ( .NOT. DESC3( NNAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( NNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Initialize grid definition
    CALL CHKGRID( NNAME, 'GRID', 0, EFLAG )

    NCOLS = NCOLS3D
    NROWS = NROWS3D
    NGRID = NCOLS3D * NROWS3D
    LUSE  = GETCFDSC( FDESC3D, '/LANDUSE/', .FALSE. )

    !......  Open and check bioseason file if using
    IF ( SWITCH_FILE ) THEN

        BNAME = PROMPTMFILE( 'Enter name for season switch input file',    &
                             FSREAD3, 'BIOSEASON', PROGNAME )

        !......  Read description of switch file
        IF ( .NOT. DESC3( BNAME ) ) THEN

            MESG = 'Could not get description of file "' // TRIM( BNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

        !......  Check grid definition
        CALL CHKGRID( BNAME, 'GRID' , 0 , EFLAG )
        EFLAG = .FALSE.

        IF ( EFLAG ) THEN
            MESG = 'Grid in file "' // TRIM( BNAME ) // '" does not match previously set grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( SWITCH( NCOLS, NROWS), STAT=IOS )
        CALL CHECKMEM( IOS, 'SWITCH', PROGNAME )
        SWITCH = 0       ! array

    ELSE

        MESG = 'Use summer normalized emissions?'
        ASSUME_SUMMER = ENVYN ( 'SUMMER_YN', MESG, .TRUE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "SUMMER_YN"', 2 )
        END IF

    END IF

    !......  Open temperature file
    M3NAME = PROMPTMFILE( 'Enter name for gridded temperature input file',    &
                          FSREAD3, 'MET_FILE1', PROGNAME )

    !......  Get description of temperature file
    IF ( .NOT. DESC3( M3NAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( M3NAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Check that grid description matches BGRD file
    CALL CHKGRID( M3NAME, 'GRID' , 0 , EFLAG )
    EFLAG = .FALSE.
    IF ( EFLAG ) THEN
        MESG = 'Grid in file "'//TRIM( M3NAME )//'" does not match previously set grid.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Check for temperature variable
    MESG = 'Variable name for temperature'
    CALL ENVSTR( 'TMPR_VAR', MESG, 'TA', TMPRNAM, IOS )
    CALL CHECK_VARIABLE( TMPRNAM, M3NAME )

    !......  Check if using PX version of MCIP
    PX_VERSION = ENVYN( 'PX_VERSION', 'MCIP is PX version?', .FALSE., IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "PX_VERSION"', 2 )
    END IF

    IF (PX_VERSION) THEN

        !......  Check for soil moisture variable
        MESG = 'Variable name for soil moisture'
        CALL ENVSTR( 'SOILM_VAR', MESG, 'SOIM1', SMNAM, IOS )
        CALL CHECK_VARIABLE( SMNAM, M3NAME )

        !......  Check for soil temperature variable
        MESG = 'Variable name for soil temperature'
        CALL ENVSTR( 'SOILT_VAR', MESG, 'SOIT1', STMPNAM, IOS )
        CALL CHECK_VARIABLE( STMPNAM, M3NAME )

        !......  Check for soil type variable
        MESG = 'Variable name for soil type'
        CALL ENVSTR( 'ISLTYP_VAR', MESG, 'SLTYP', STNAM, IOS )
        CALL CHECK_VARIABLE( STNAM, M3NAME )

        ALLOCATE( SOILM( NCOLS, NROWS ),    &
                  SOILT( NCOLS, NROWS ),    &
                 ISLTYP( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SOILM...ISLTYP', PROGNAME )

    END IF

    !......  Check for nonconvective rainfall variable
    MESG = 'Variable name for nonconvective rainfall'
    CALL ENVSTR( 'RN_VAR', MESG, 'RN', RNNAM, IOS )
    CALL CHECK_VARIABLE( RNNAM, M3NAME )

    !......  Check for convective rainfall variable
    MESG = 'Variable name for convective rainfall'
    CALL ENVSTR( 'RC_VAR', MESG, 'RC', RCNAM, IOS )
    CALL CHECK_VARIABLE( RCNAM, M3NAME )

    !......  Check for pressure variable; could be in either MET_FILE1 or MET_FILE2
    MESG = 'Variable name for surface pressure'
    CALL ENVSTR( 'PRES_VAR', MESG, 'PRES', PRESNAM, IOS )

    J = INDEX1( PRESNAM, NVARS3D, VNAME3D )
    IF( J <= 0 ) THEN
        IF( SAMEFILE ) THEN
            MESG = 'Could not find "'//TRIM( PRESNAM )//'" in file "'//TRIM( M3NAME )//'"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    ELSE
        PNAME = M3NAME
    END IF

    !......  Get met and cloud scheme descriptions from M3NAME if they exist
    METSCEN  = GETCFDSC( FDESC3D, '/MET SCENARIO/', .FALSE. )
    CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

    !......  Open second met file if needed
    IF ( .NOT. SAMEFILE ) THEN

        M2NAME = PROMPTMFILE( 'Enter name for gridded radiation/cloud input file',    &
                              FSREAD3, 'MET_FILE2', PROGNAME )

        !......  Get description of radiation/cloud file
        IF ( .NOT. DESC3( M2NAME ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( M2NAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !......  Check that grid description matches BGRD file
        CALL CHKGRID( M2NAME, 'GRID' , 0 , EFLAG )
        EFLAG = .FALSE.
        IF ( EFLAG ) THEN
            MESG = 'Grid in file "'//TRIM( M2NAME )//'" does not match previously set grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CLOUDSHM = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .FALSE. )

        !......  Check for pressure variable
        IF ( PNAME .NE. M3NAME ) THEN
            CALL CHECK_VARIABLE( PRESNAM, M2NAME )
            PNAME = M2NAME
        END IF

    !......  If not using separate files, set M2NAME equal to first file
    ELSE

        M2NAME = M3NAME

    END IF

    !......  Check for radiation variable
    MESG = 'Variable name for radiation'
    CALL ENVSTR( 'RGRND_VAR', MESG, 'RGRND', RGRNDNAM, IOS )
    CALL CHECK_VARIABLE( RGRNDNAM, M2NAME )


    !......  Check for inverse aerodynamic resistance variable
    MESG = 'Variable name for inverse aerodynamic resistance'
    CALL ENVSTR( 'RADYNI_VAR', MESG, 'RADYNI', RADYNINAM, IOS )
    CALL CHECK_VARIABLE( RADYNINAM, M3NAME )

    !......  Check for ground temp variable
    MESG = 'Variable name for ground temperature'
    CALL ENVSTR( 'TEMPG_VAR', MESG, 'TEMPG', TEMPGNAM, IOS )
    CALL CHECK_VARIABLE( TEMPGNAM, M3NAME )

    !......  Check for LAI variable
    MESG = 'Variable name for LAI '
    CALL ENVSTR( 'LAI_VAR', MESG, 'LAI', LAINAM, IOS )
    CALL CHECK_VARIABLE( LAINAM, M3NAME )

    !......  Check for Q2 variable
    MESG = 'Variable name for Q2 '
    CALL ENVSTR( 'Q2_VAR', MESG, 'Q2', Q2NAM, IOS )
    CALL CHECK_VARIABLE( Q2NAM, M3NAME )

    !......  Check for RSTOMI variable
    MESG = 'Variable name for RSTOMI '
    CALL ENVSTR( 'RSTOMI_VAR', MESG, 'RSTOMI', RSTOMINAM, IOS )
    CALL CHECK_VARIABLE( RSTOMINAM, M3NAME )

    !......  Check for USTAR variable
    MESG = 'Variable name for USTAR '
    CALL ENVSTR( 'USTAR_VAR', MESG, 'USTAR', USTARNAM, IOS )
    CALL CHECK_VARIABLE( USTARNAM, M3NAME )

    !......  Set the EV name for met data description
    MESG = 'Setting for the environment variable name for file ' //    &
             'meta description for output file'
    CALL ENVSTR( 'MRG_FILEDESC', MESG, ' ', MRGFDESC, IOS  )

    IF( IOS >= 0 ) THEN
        FDFLAG = .TRUE.
        MESG = 'Use this file meta description for output file'
        CALL ENVSTR( MRGFDESC, MESG, ' ', METADESC, IOS )
        IF( IOS < 0 ) FDFLAG = .FALSE.
    END IF

    !......  Get default time characteristic for output file:
    !        If we're going to prompt, then set the defaults based on met
    !        otherwise, use environment variables to set defaults
    JDATE  = SDATE3D
    JTIME  = STIME3D

    !......  Check if first meteorology file is virtual
    CALL NAMEVAL( M3NAME, EQNAME )
    CHKNAME = EQNAME( 1:8 )
    CALL UPCASE( CHKNAME )

    IF ( CHKNAME .EQ. 'VIRTUAL ' ) THEN
        JRUNLEN = ENVINT( 'G_RUNLEN', 'Duration(HHMMSS)', 490000,    &
                          IOS )
        NSTEPS  = JRUNLEN / 10000

        WRITE( MESG,94010 ) 'WARNING: Assuming execution in '//    &
               'virtual mode; duration set to ', NSTEPS, 'hours.'

        CALL M3MSG2( MESG )
    ELSE
        NSTEPS = MXREC3D
    END IF
    TSTEP = 10000

    CALL GETM3EPI( TZONE, JDATE, JTIME, TSTEP, NSTEPS )

    !......  Build description for, and create/open output file
    !        (all but variables-table in description is borrowed from M3NAME)
    SDATE3D = JDATE
    STIME3D = JTIME
    BSTEPS  = NSTEPS
    MXREC3D = BSTEPS
    NVARS3D = MSPCS
    NLAYS3D = 1
    TSTEP3D = 10000       ! only 1-hour timestep supported

    DO M = 1, MSPCS
        VNAME3D( M ) = EMSPC( M )
        UNITS3D( M ) = UNITSMENU( UNITTYPE )
        VDESC3D( M ) = 'biogenic emissions of the indicated species'
        VTYPE3D( M ) = M3REAL
    END DO

    FDESC3D = ' '       ! array

    FDESC3D( 1 ) = 'Gridded biogenic emissions from SMOKE-BEIS3'
    IF( FDFLAG ) FDESC3D( 1 ) = METADESC
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
    FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
    FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
    FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
    FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM

    IF( VFLAG ) THEN
        FDESC3D( 8 ) = '/VARIABLE GRID/ ' // GDNAM3D
    END IF

    !......  Open first output file (moles/time period)
    LNAME = PROMPTMFILE( 'Enter name for B3GTS output file - moles',    &
                         FSUNKN3, 'B3GTS_L', PROGNAME )

    !......  Open second output file (tons/hour)
    DO M = 1, MSPCS
        UNITS3D( M ) = 'tons/hr'
    END DO

    SNAME = PROMPTMFILE( 'Enter name for B3GTS output file - mass',    &
                         FSUNKN3, 'B3GTS_S', PROGNAME )

    !......  Build name table for variables in normalized emissions file
    !......  Allocate memory for data and read
    ALLOCATE( AVGEMIS( NCOLS, NROWS, NSEF, NSEASONS ),    &
               NOEMIS( NCOLS, NROWS, NNO ),    &
               AVGLAI( NCOLS, NROWS, NLAI, NSEASONS ),    &
                  LAT( NCOLS, NROWS ),    &
                  LON( NCOLS, NROWS ),    &
               COSZEN( NCOLS, NROWS ),    &
                TASFC( NCOLS, NROWS ),    &
                   RN( NCOLS, NROWS ),    &
                   RC( NCOLS, NROWS ),    &
                TEMPG( NCOLS, NROWS ),    &
               RADYNI( NCOLS, NROWS ),    &
                 RATM( NCOLS, NROWS ),    &
                  LAI( NCOLS, NROWS ),    &
                   Q2( NCOLS, NROWS ),    &
               RSTOMI( NCOLS, NROWS ),    &
                RSTOM( NCOLS, NROWS ),    &
                USTAR( NCOLS, NROWS ),    &
                RGRND( NCOLS, NROWS ),    &
             RAINFALL( NCOLS, NROWS, RHOURS ),    &
                PTYPE( NCOLS, NROWS ),    &
            PULSEDATE( NCOLS, NROWS ),    &
            PULSETIME( NCOLS, NROWS ),    &
                 PRES( NCOLS, NROWS ),    &
                EMPOL( NCOLS, NROWS, NSEF ),    &
                EMISL( NCOLS, NROWS, MSPCS ),    &
                EMISS( NCOLS, NROWS, MSPCS ),    &
                SEMIS( NCOLS, NROWS, NSEF ),    &
              NONAGNO( NCOLS, NROWS ),    &
            NGROWAGNO( NCOLS, NROWS ),    &
             GROWAGNO( NCOLS, NROWS ),    &
                 SLAI( NCOLS, NROWS, NLAI ), STAT=IOS )
    CALL CHECKMEM( IOS, 'AVGEMIS...SLAI', PROGNAME )

    AVGEMIS = 0.0     ! array
    NOEMIS  = 0.0     ! array
    AVGLAI  = 0.0     ! array
    !......  Build name table for variables in normalized emissions file
    !......  Open 2-D grid parameters file to get LAT and LON
    GNAME = PROMPTMFILE( 'Enter name for 2D GRID PARAMETERS input file',    &
                         FSREAD3, 'GRID_CRO_2D', PROGNAME )

    IF( .NOT. DESC3( GNAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !......  Check grid description against BGRD File
    CALL CHKGRID( GNAME, 'GRID' , 0 , EFLAG )
    EFLAG = .FALSE.
    IF( EFLAG ) THEN
        MESG = 'Grid in file "'//TRIM(GNAME)//    &
               '" does not match previously set grid.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF


    IF( .NOT. READ3( GNAME, 'LAT', 1,0,0, LAT ) ) THEN
        MESG = 'Could not read LAT from "' //    &
                TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( .NOT. READ3( GNAME, 'LON', 1,0,0, LON ) ) THEN
        MESG = 'Could not read LON from "' //    &
                TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF



    !......  Read the various categories of normalized emissions
    DO M = 1, NSEASONS

        DO B = 1, NSEF-2
            VTMP = 'AVG_' // TRIM( BIOTYPES( B ) ) // SEASON( M )

            IF( .NOT. READ3( NNAME, VTMP, 1,0,0, AVGEMIS(1,1,B,M) ) ) THEN
                MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( NNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO

        DO B = NSEF, NSEF
            VTMP = 'AVG_' // TRIM( BIOTYPES( B ) ) // SEASON( M )

            IF( .NOT. READ3( NNAME, VTMP, 1,0,0, AVGEMIS(1,1,B,M) ) ) THEN
                MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( NNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO


        DO N = 1, NLAI
            VTMP = 'LAI_' // TRIM( LAITYPES( N ) ) // SEASON( M )

            IF ( .NOT. READ3( NNAME, VTMP, 1,0,0, AVGLAI(1,1,N,M) ) ) THEN
                MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( NNAME ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END DO

    END DO       ! end loop over seasons

    VTMP = 'AVG_NOAG_GROW'

    IF ( .NOT. READ3( NNAME, VTMP, 1,0,0,    &
                      NOEMIS( 1, 1, 1 ) ) ) THEN
        MESG = 'Could not read "' // TRIM(VTMP) //    &
               '" from "' // TRIM( NNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    VTMP = 'AVG_NOAG_NONGROW'

    IF ( .NOT. READ3( NNAME, VTMP, 1,0,0,    &
                      NOEMIS( 1, 1, 2 ) ) ) THEN
        MESG = 'Could not read "' // TRIM(VTMP) //    &
               '" from "' // TRIM( NNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    VTMP = 'AVG_NONONAG'

    IF ( .NOT. READ3( NNAME, VTMP, 1,0,0,    &
                      NOEMIS( 1, 1, 3 ) ) ) THEN
        MESG = 'Could not read "' // TRIM(VTMP) //    &
               '" from "' // TRIM( NNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MESG = 'Initial run?'
    INITIAL_RUN = ENVYN( 'INITIAL_RUN', MESG, .FALSE., IOS)

    !......  If initial run, initialize some variables, otherwise get them from
    IF( INITIAL_RUN ) THEN
        PULSEDATE = 0       ! array
        PULSETIME = 0       ! array
        PTYPE     = 0       ! array
        RVARS     = 0
        RINDEX    = 1
    ELSE

        !......  Open saved NO emissions file
        SOILINP = PROMPTMFILE( 'Enter name for NO EMISSIONS SAVE file',    &
                               FSREAD3, 'SOILINP', PROGNAME )

        !......  Get description of NO emissions file
        IF( .NOT. DESC3( SOILINP ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( SOILINP ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        RVARS = NVARS3D-3
        RINDEX = NTHIK3D
        NDATE = SDATE3D
        NTIME = STIME3D

        !......  Check that file's start date and time are consistent
        IF( NDATE /= JDATE ) THEN
            WRITE( MESG,94010 ) 'Cannot use SOILINP file; ' //    &
                'found date ', NDATE, ' expecting ', JDATE
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( NTIME /= JTIME ) THEN
            WRITE( MESG,94010 ) 'Cannot use SOILINP file; ' //    &
                'found time ', NTIME, ' expecting ', JTIME
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !......  Read data from
        VTMP = 'PTYPE'
        IF( .NOT. READ3( SOILINP, VTMP, 1, JDATE, JTIME, PTYPE ) ) THEN
            MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( SOILINP ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        VTMP = 'PULSEDATE'
        IF( .NOT. READ3( SOILINP, VTMP,1, JDATE, JTIME, PULSEDATE ) ) THEN
            MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( SOILINP ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        VTMP = 'PULSETIME'
        IF( .NOT. READ3( SOILINP, VTMP, 1, JDATE, JTIME, PULSETIME ) ) THEN
            MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( SOILINP ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        DO I = 1, RVARS

            !......  Build variable name
            WRITE( VTMP, '(A8,I2.2)' ) 'RAINFALL', I

            IF( .NOT. READ3( SOILINP, VTMP, 1,JDATE,JTIME, RAINFALL( 1,1,I ) ) ) THEN
                MESG = 'Could not read "' // TRIM(VTMP) // '" from "' // TRIM( SOILINP ) // '"'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF
        END DO

    END IF


    !......  Set up gridded met file(s) dates and times for specific time zone
    CALL PREBMET( M3NAME, M2NAME, SAMEFILE, TZONE, TSTEP3D,    &
                  JDATE, JTIME, NSTEPS, MDATE, MTIME, RDATE, RTIME )

    IF( M2NAME .EQ. M3NAME ) THEN
        RDATE = MDATE
        RTIME = MTIME
    ENDIF

    !...... Initialize normalized emissons to be used
    IF( ASSUME_SUMMER ) THEN

        SEMIS = AVGEMIS( 1:NCOLS, 1:NROWS, 1:NSEF,   NSUMMER )
        SLAI  = AVGLAI ( 1:NCOLS, 1:NROWS, 1:NLAI,   NSUMMER )

    ELSE

        SEMIS = AVGEMIS( 1:NCOLS, 1:NROWS, 1:NSEF,   NWINTER )
        SLAI  = AVGLAI ( 1:NCOLS, 1:NROWS, 1:NLAI,   NWINTER )

    END IF

    GROWAGNO  = NOEMIS( 1:NCOLS, 1:NROWS, 1 )
    NGROWAGNO = NOEMIS( 1:NCOLS, 1:NROWS, 2 )
    NONAGNO   = NOEMIS( 1:NCOLS, 1:NROWS, 3 )

    !......  Loop thru the number of time steps (hourly)
    LDATE = 0

    DO HR = 1, BSTEPS

        EMISL = 0       !  array
        EMISS = 0       !  array
        EMPOL = 0       !  array

        IF( JDATE .NE. LDATE ) THEN

            CALL WRDAYMSG( JDATE, MESG )

            !......  If new date, read season switch
            IF ( SWITCH_FILE ) THEN
                MESG = 'Reading gridded season switch data...'
                CALL M3MSG2( MESG )

                IF( .NOT. READ3( BNAME, 'SEASON', 1, JDATE, 0, SWITCH ) ) THEN
                    MESG = 'Could not read SEASON from "' // TRIM( BNAME ) // '"'
                    CALL M3EXIT( PROGNAME, JDATE, 0, MESG, 2 )
                END IF

                MESG = 'Applying gridded season switch data...'
                CALL M3MSG2( MESG )

                DO J = 1, NROWS
                DO I = 1, NCOLS

                    !......  If switch equal to 0 use winter normalized emissions
                    IF( SWITCH( I,J ) == 0 ) THEN
                        SEMIS( I, J, 1:NSEF ) = AVGEMIS( I, J, 1:NSEF, NWINTER )
                        SLAI ( I, J, 1:NLAI ) = AVGLAI ( I, J, 1:NLAI, NWINTER )
                    ELSE
                        SEMIS( I, J, 1:NSEF ) = AVGEMIS( I, J, 1:NSEF, NSUMMER )
                        SLAI ( I, J, 1:NLAI ) = AVGLAI ( I, J, 1:NLAI, NSUMMER )
                    END IF
                END DO                 ! loop over rows
                END DO     ! loop over columns

            END IF      ! if using switch file

        END IF      ! if new day

        WRITE( MESG,94030 ) HHMMSS( JTIME )
        CALL M3MSG2( MESG )


        !......  Compute zenith angle
        CALL CZANGLE( JDATE, JTIME, NCOLS, NROWS, LAT, LON,    &
                      COSZEN, GETATN )


        !......  Read temperature data
        IF( .NOT. READ3( M3NAME, TMPRNAM, 1, MDATE, MTIME, TASFC ) ) THEN
            MESG = 'Could not read "' // TMPRNAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
        END IF

        IF( PX_VERSION ) THEN

            !......  Read soil moisture data
            IF( .NOT. READ3( M3NAME, SMNAM, 1, MDATE, MTIME, SOILM ) ) THEN
                MESG = 'Could not read "' // SMNAM // '" from "' // TRIM( M3NAME ) // '"'
                CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
            END IF

            !......  Read soil temperature data
            IF( .NOT. READ3( M3NAME, STMPNAM, 1, MDATE, MTIME, SOILT ) ) THEN
                MESG = 'Could not read "' // STMPNAM // '" from "' // TRIM( M3NAME ) // '"'
                CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
            END IF

            !......  Read soil type data
            IF( .NOT. READ3( M3NAME, STNAM, 1, MDATE, MTIME, ISLTYP ) ) THEN
                MESG = 'Could not read "' // STNAM // '" from "' // TRIM( M3NAME ) // '"'
                CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
            END IF

        END IF

        !......  Read rainfall data
        IF( .NOT. READ3( M3NAME, RNNAM, 1, MDATE, MTIME, RN ) ) THEN
            MESG = 'Could not read "' //  RNNAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, RCNAM, 1, MDATE, MTIME, RC ) ) THEN
            MESG = 'Could not read "' //  RCNAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
        END IF

        !......  Read surface pressure data
        IF( .NOT. READ3( PNAME, PRESNAM, 1, MDATE, MTIME, PRES ) ) THEN
            MESG = 'Could not read "' // PRESNAM // '" from "' // TRIM( PNAME ) // '"'
            CALL M3EXIT( PROGNAME, MDATE, MTIME, MESG, 2 )
        END IF

        !......  Convert pressure from Pa to millibars
        !        PRES = PRES * 0.010           ! array   (CONVERSION NOT NEEDED IN VERSION 3.50)

        !......  Read radiation data
        IF( .NOT. READ3( M3NAME, RGRNDNAM, 1, RDATE, RTIME, RGRND ) ) THEN
            MESG = 'Could not read "' // RGRNDNAM //'" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, TEMPGNAM, 1, RDATE, RTIME, TEMPG ) ) THEN
            MESG = 'Could not read "' // TEMPGNAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, RADYNINAM, 1,RDATE, RTIME, RADYNI ) ) THEN
            MESG = 'Could not read "' // RADYNINAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, LAINAM, 1, RDATE, RTIME, LAI ) ) THEN
            MESG = 'Could not read "' // LAINAM // '" from "' // TRIM( M3NAME ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, Q2NAM, 1, RDATE, RTIME, Q2 ) ) THEN
            MESG = 'Could not read "' // Q2NAM // '" from "' // TRIM( Q2NAM ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, RSTOMINAM, 1, RDATE, RTIME, RSTOMI ) ) THEN
            MESG = 'Could not read "' // RSTOMINAM // '" from "' // TRIM( RSTOMINAM ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        IF( .NOT. READ3( M3NAME, USTARNAM, 1, RDATE, RTIME, USTAR ) ) THEN
            MESG = 'Could not read "' // USTARNAM // '" from "' // TRIM( USTARNAM ) // '"'
            CALL M3EXIT( PROGNAME, RDATE, RTIME, MESG, 2 )
        END IF

        ! ARRAY ARITHMETIC


        RSTOM( : , : ) = 1.0/RSTOMI( : , : )
        RATM ( : , : ) = 1.0/RADYNI( : , : )

        IF(( INITIAL_RUN ) .OR. (RVARS <= RHOURS)) THEN
        !        IF( (HR+RVARS) <= (RHOURS - 1) ) THEN

            IF( (HR+RVARS) <= (RHOURS ) ) THEN
                INITIAL_HOUR = .TRUE.
            ELSE
                INITIAL_HOUR = .FALSE.
            ENDIF
        END IF

        !......  Calculate hourly rainfall totals
        !        INDEX = MOD( HR-1+RVARS+RINDEX, RHOURS ) + 1           ! inital run RVARS = 0, subsequent runs RVARS
        RAINFALL( 1:NCOLS, 1:NROWS, RINDEX ) = RN + RC


        RN = 0.
        DO I = 1, MIN(RHOURS,RVARS+HR)           ! only sum rainfall if we have at least 24 hrs
            RN = RN + RAINFALL( 1:NCOLS, 1:NROWS, I )
        END DO



        IF ((HR+RVARS > RHOURS)) THEN
        !	    IF ((HR+RVARS >= RHOURS)) THEN
            INITIAL_HOUR = .FALSE.
            INITIAL_RUN = .FALSE.
        ENDIF


        !......  Calculate non-speciated emissions
        CALL HRBEIS360( MDATE, MTIME, NCOLS, NROWS, MSPCS,         &
                        PX_VERSION, INITIAL_HOUR, COSZEN, SEMIS,   &
                        GROWAGNO, NGROWAGNO, NONAGNO, SLAI, TASFC, &
                        SOILM, SOILT, ISLTYP, RN, PRES,            &
                        RGRND,LAI,USTAR,RSTOM,RATM,Q2,TEMPG,       &
                        PTYPE, PULSEDATE, PULSETIME, EMPOL )

        !...... Speciate emissions
        DO L = 1, MSPCS
        DO K = 1, NSEF
        DO J = 1, NROWS
        DO I = 1, NCOLS
            EMISL( I, J, L  ) = EMISL( I ,J, L ) + EMPOL( I, J, K ) * MLFAC( L, K )
            EMISS( I, J, L  ) = EMISS( I, J, L ) + EMPOL( I, J, K ) * MSFAC( L, K )
        END DO
        END DO
        END DO
        END DO

        !......  Convert to moles/second if necessary
        IF ( UNITTYPE .EQ. 2 ) THEN
            EMISL = EMISL * HR2SEC            ! array
        END IF

        !......  Write out speciated emissions
        IF( .NOT. WRITE3( LNAME, 'ALL', JDATE, JTIME, EMISL ) ) THEN
            MESG = 'Could not write to output file "' // TRIM( LNAME ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        IF( .NOT. WRITE3( SNAME, 'ALL', JDATE, JTIME, EMISS ) ) THEN
            MESG = 'Could not write to output file "' // TRIM( SNAME ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF

        !......  Increment time step
        LDATE = JDATE

        IF( HR /= BSTEPS ) THEN
            CALL NEXTIME( JDATE, JTIME, 10000 )
        END IF

        CALL NEXTIME( MDATE, MTIME, 10000 )

        IF ( M2NAME .EQ. M3NAME ) THEN
            RDATE = MDATE
            RTIME = MTIME
        ELSE
            CALL NEXTIME( RDATE, RTIME, 10000 )
        END IF

        RINDEX = MOD(RINDEX,RHOURS) + 1

    END DO     ! loop over hours

    !......  Create saved NO emissions file

    !......  Build description for, and create/open output file
    NCOLS3D = NCOLS
    NROWS3D = NROWS
    NLAYS3D = 1
    SDATE3D = JDATE
    STIME3D = JTIME
    MXREC3D = 1

    RVARS = BSTEPS + RVARS -1       ! don't save last hour of rainfall
    IF (RVARS .GE. 24) THEN
        RVARS = 24
    ENDIF
    NTHIK3D= RINDEX-1


    NVARS3D = 3 + RVARS
    TSTEP3D = 10000

    VNAME3D = ' '
    VNAME3D( 1 ) = 'PTYPE'
    VNAME3D( 2 ) = 'PULSEDATE'
    VNAME3D( 3 ) = 'PULSETIME'

    DO I = 1, RVARS
        WRITE( VTMP, '(A8,I2.2)' ) 'RAINFALL', I
        VNAME3D( I+3 ) = VTMP
    END DO

    UNITS3D = ' '
    UNITS3D( 1 ) = 'NUMBER'
    UNITS3D( 2 ) = 'M3DATE'
    UNITS3D( 3 ) = 'M3TIME'
    UNITS3D( 4:4+RVARS-1 ) = 'CM'

    VDESC3D( 1 ) = 'NO EMISSION PULSE TYPE'
    VDESC3D( 2 ) = 'MODELS-3 STARTING DATE FOR NO EMISSION PULSE'
    VDESC3D( 3 ) = 'MODELS-3 STARTING TIME FOR NO EMISSION PULSE'
    VDESC3D( 4:4+RVARS-1 ) = 'TOTAL RAINFALL for 6 HOURS'

    VTYPE3D = 0
    VTYPE3D( 1 ) = M3INT
    VTYPE3D( 2 ) = M3INT
    VTYPE3D( 3 ) = M3INT
    VTYPE3D( 4:4+RVARS-1 ) = M3REAL

    FDESC3D = ' '
    FDESC3D( 1 ) = 'Gridded NO emission storage'
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
    FDESC3D( 4 ) = '/TZONE/ '   // CTZONE
    FDESC3D( 5 ) = '/LANDUSE/ ' // LUSE
    FDESC3D( 6 ) = '/MET SCENARIO/ ' // METSCEN
    FDESC3D( 7 ) = '/CLOUD SCHEME/ ' // CLOUDSHM


    !......  Open emissions file
    SOILOUT = PROMPTMFILE( 'Enter name for NO EMISSIONS SOIL file',    &
                           FSUNKN3, 'SOILOUT', PROGNAME )

    !......  Write emissions file
    VTMP = 'PTYPE'
    IF( .NOT. WRITE3( SOILOUT, VTMP, JDATE, JTIME, PTYPE ) ) THEN
        MESG = 'Could not write "' // TRIM(VTMP) // '" to "' // TRIM( SOILOUT ) // '"'
        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
    END IF

    VTMP = 'PULSEDATE'
    IF( .NOT. WRITE3( SOILOUT, VTMP, JDATE, JTIME, PULSEDATE ) ) THEN
        MESG = 'Could not write "' // TRIM(VTMP) // '" to "' // TRIM( SOILOUT ) // '"'
        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
    END IF

    VTMP = 'PULSETIME'
    IF( .NOT. WRITE3( SOILOUT, VTMP, JDATE, JTIME, PULSETIME ) ) THEN
        MESG = 'Could not write "' // TRIM(VTMP) // '" to "' // TRIM( SOILOUT ) // '"'
        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
    END IF

    DO I = 1, RVARS
        WRITE( VTMP, '(A8,I2.2)' ) 'RAINFALL', I

        IF( .NOT. WRITE3( SOILOUT, VTMP, JDATE, JTIME, RAINFALL( 1,1,I ) ) ) THEN
            MESG = 'Could not write "' // TRIM(VTMP) // '" to "' // TRIM( SOILOUT ) // '"'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
        END IF
    END DO

    !......   End of subroutine
    RETURN


    !******************  FORMAT  STATEMENTS   ******************************

    !......   Informational (LOG) message formats... 92xxx

92000 FORMAT ( 5X , A )

    !......   Internal buffering formats...... 94xxx

94000 FORMAT( I2.2 )
94010 FORMAT( 10( A, :, I8, :, 1X ) )

94030 FORMAT( 8X, 'at time ', A8 )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !......  This internal subprogram checks if the given variable
    !        is available in the I/O API file
    SUBROUTINE CHECK_VARIABLE( VNAME, FNAME )

        !......  Subprogram arguments
        CHARACTER(*), INTENT(IN   ) ::  VNAME           ! variable name
        CHARACTER(*), INTENT(IN   ) :: FNAME          ! file name

        J = INDEX1( VNAME, NVARS3D, VNAME3D )

        IF( J <= 0 ) THEN
            MESG = 'Could not find var "' // TRIM( VNAME ) // '" in "' // TRIM( FNAME ) // '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

    END SUBROUTINE CHECK_VARIABLE

    !----------------------------------------------------------------------

END SUBROUTINE TMPBEIS360
