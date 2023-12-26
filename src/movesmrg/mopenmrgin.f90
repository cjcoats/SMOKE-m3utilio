
SUBROUTINE OPENMRGIN

    !***********************************************************************
    !  subroutine OPENMRGIN body starts at line
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to open all of the necessary
    !      files for the merge routine and set the episode information
    !      for the calling program.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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
    !****************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains the major data structure and control flags
    USE MODMERGE, ONLY:                                     &
            MENAME, MSDEV, CFDEV,                           &
            NMSRC, MPRJFLAG, MFLAG_BD, MTNAME, MSDATE,      &
            MNIPPA, MEANAM, MGNAME, MNGMAT,                 &
            PDEV, CDEV, TZONE, SDATE,                       &
            STIME, TSTEP, NSTEPS, EDATE, ETIME, BYEAR, PYEAR,&
            VARFLAG, SRCGRPFLAG, SGDEV, SUBSECFLAG

    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, RPPFLAG, RPHFLAG, ONIFLAG, RPSFLAG,  &
            TVARNAME, METNAME, XDEV, MDEV, FDEV, CFFLAG, SPDISTFLAG,            &
            SPDPROFLAG, MSNAME_L, MSNAME_S, MNSMATV_L, MNSMATV_S,               &
            MSVDESC_L, MSVDESC_S, MSVUNIT_L, MSVUNIT_S, ETABLEFLAG

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL, NIACT, NSRC, CATEGORY, ACTVTY

    !.......  This module contains the inventory arrays
    USE MODSOURC, ONLY: SPEED, VPOP, CIFIP

    !.......   This module contains emission factor information
    USE MODEMFAC, ONLY: MXETYPE, EMTNAM

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, VGTYP, VGTOP, VGLVS

    !.......  This module is required for the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......  INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......  EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL      , EXTERNAL :: DSCM3GRD
    CHARACTER(50), EXTERNAL :: GETCFDSC
    INTEGER      , EXTERNAL :: GETIFDSC

    !.......   LOCAL VARIABLES and their descriptions:

    !.......  Array that contains the names of the inventory variables to read
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

    !.......  Other local variables

    INTEGER         I, J, K, M, N, V         ! counters and indices

    INTEGER         IDEV              ! tmp unit number if ENAME is map file
    INTEGER         TDEV              ! unit number for MEPROC file
    INTEGER         SPDEV             ! unit number for SPDPRO file
    INTEGER         SDDEV             ! unit number for SPDIST file
    INTEGER         IOS               ! tmp I/O status
    INTEGER         ISECS             ! tmp duration in seconds
    INTEGER         NPACT             ! no. variables per activity
    INTEGER         NPPOL             ! no. variables per pollutant
    INTEGER         NDIM              ! tmp dimensioning variable
    INTEGER         NVAR              ! tmp no. variables
    INTEGER         NINVARR           ! number inventory variables to read

    LOGICAL      :: EFLAG = .FALSE.      ! true: error in routine
    LOGICAL      :: IFLAG = .FALSE.      ! true: episode settings have been init
    LOGICAL      :: OFLAG = .FALSE.      ! true: met info has been init
    LOGICAL      :: YFLAG = .FALSE.      ! true: year/projection info been init
    LOGICAL      :: ZFLAG = .FALSE.      ! true: time zone has been init

    CHARACTER(16)   DUMNAME          ! tmp file name
    CHARACTER(16)   INAME            ! tmp name for inven file of unknown fmt
    CHARACTER(50)   METSCENR         ! met scenario name
    CHARACTER(50)   METCLOUD         ! met cloud scheme name
    CHARACTER(50)   METTMP           ! temporary buffer for met info
    CHARACTER(80)   GDESC            ! grid description
    CHARACTER(256)  MESG             ! message buffer
    CHARACTER(NAMLEN3) COORD3D       ! coordinate system name
    CHARACTER(NAMLEN3) COORUN3D      ! coordinate system projection units
    CHARACTER(NAMLEN3) PROJTYPE      ! projection type
    CHARACTER(NAMLEN3) OUTGRDNM      ! output grid name

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENMRGIN'     ! program name

    !***********************************************************************
    !   begin body of subroutine OPENMRGIN

    !.......  Initialize gridded information with grid description file
    IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,     &
                        P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,    &
                        XORIG3D, YORIG3D, XCELL3D, YCELL3D,             &
                        NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF
    OUTGRDNM = GDNAM3D

    !.......  Check or initialize the grid; do not allow subgrids
    !           when using a variable grid
    IF( VARFLAG ) THEN
        CALL CHKGRID( 'general', 'GRIDDESC', 0, EFLAG )
    ELSE
        CALL CHKGRID( 'general', 'GRIDDESC', 1, EFLAG )
    END IF

    !.......  Get inventory file names given source category
    CALL GETINAME( 'MOBILE', MENAME, DUMNAME )

    !.......  Prompt for and open inventory file
    MESG= 'Enter logical name for the MAP MOBILE INVENTORY file'
    INAME = MENAME
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.......  Read map-formatted inventory file
    CALL RDINVMAP( INAME, IDEV, MENAME, DUMNAME, MSDEV )

    !.......  Store source-category-specific header information,
    !           including the inventory pollutants in the file (if any).  Note that
    !           the I/O API header info is passed by include file and the
    !           results are stored in module MODINFO.
    CALL GETSINFO( MENAME )

    !.......  Ensure that there is at least one activity in the inventory
    !           file, or else this program does not need to be run
    IF( NIACT == 0 ) THEN
        MESG = 'No activities are found in the inventory file!  Program cannot be used.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Set inventory variables to read
    IVARNAMS( 1 ) = 'CIFIP'
    IVARNAMS( 2 ) = 'CSCC'
    IVARNAMS( 3 ) = 'TZONES'
    NINVARR = 3

    !.......  Allocate memory for and read required inventory characteristics
    CALL RDINVCHR( CATEGORY, MENAME, MSDEV, NSRC, NINVARR, IVARNAMS )

    !.......  Read speed and vehicle population data from the inventory
    IF( RPDFLAG ) THEN
        M = INDEX1( 'SPEED', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include speed data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( SPEED( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPEED', PROGNAME )
        CALL RDMAPPOL( NSRC, 1, 1, 'SPEED', SPEED )

        !.......  Make sure inventory has VMT as activity (won't be using this
        !               data but it needs to be there to make emission processes work)
        M = INDEX1( 'VMT', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include VMT data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    !.......  Read hotelling data from the inventory
    IF( RPHFLAG ) THEN
        M = INDEX1( 'HOTELLING', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include hotelling data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    IF( RPSFLAG ) THEN
        M = INDEX1( 'STARTS', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include engine start data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    IF( ONIFLAG ) THEN
        M = INDEX1( 'IDLING', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include idling data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    IF( RPVFLAG .OR. RPPFLAG ) THEN
        M = INDEX1( 'VPOP', NMAP, MAPNAM )
        IF( M <= 0 ) THEN
            MESG = 'Mobile inventory does not include vehicle ' //&
                   'population data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        ALLOCATE( VPOP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VPOP', PROGNAME )
        CALL RDMAPPOL( NSRC, 1, 1, 'VPOP', VPOP )
    END IF

    !.......  Build unique lists of SCCs and country/state/county codes
    !           from the inventory arrays
    CALL GENUSLST

    !.......  Get number of sources from MODINFO and store in MODMERGE variable
    NMSRC = NSRC

    !.......  Determine the year and projection status of the inventory
    CALL CHECK_INVYEAR( MENAME, MPRJFLAG, FDESC3D )

    IF( RPDFLAG .OR. RPHFLAG .OR. RPSFLAG .OR. ONIFLAG  ) THEN

        !.......  Open all temporal files for either by-day or standard
        !               processing.
        !.......  Compare headers to make sure files are consistent.
        CALL OPEN_TMP_FILES( 'MOBILE', MFLAG_BD, MTNAME, MSDATE)

        !.......  Determine the year and projection status of the hourly
        CALL CHECK_INVYEAR( MTNAME( 1 ), MPRJFLAG, FDESC3D )

    END IF

    !.......  Open gridding matrix, compare number of sources, and
    !           compare grid information
    MGNAME = PROMPTMFILE( 'Enter logical name for the MOBILE GRIDDING MATRIX',&
                          FSREAD3, 'MGMAT', PROGNAME )

    IF( .NOT. DESC3( MGNAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( MGNAME ) // '" '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( VARFLAG ) THEN
        DUMNAME = GETCFDSC( FDESC3D, '/VARIABLE GRID/', .TRUE. )
    END IF

    CALL CHKSRCNO( 'mobile', 'MGMAT', NTHIK3D, NMSRC, EFLAG )

    !.......  Check the grid definition; do not allow subgrids if using
    !           a variable grid
    IF( VARFLAG ) THEN
        CALL CHKGRID( 'mobile', 'GMAT', 0, EFLAG )
    ELSE
        CALL CHKGRID( 'mobile', 'GMAT', 1, EFLAG )
    END IF

    MNGMAT = NCOLS3D

    !.......  Open mole-based speciation matrix, compare number of sources, and store
    !           speciation variable descriptions.
    MSNAME_L = PROMPTSET( 'Enter logical name for the MOLE-BASED SPECIATION MATRIX',&
                          FSREAD3, 'MSMAT_L', PROGNAME )

    IF ( .NOT. DESCSET( MSNAME_L, ALLFILES ) ) THEN
        MESG = 'Could not get description of file set "' // TRIM( MSNAME_L ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    CALL CHKSRCNO( 'mobile', 'MSMAT_L', NROWS3D, NMSRC, EFLAG)
    MNSMATV_L = NVARSET
    ALLOCATE( MSVDESC_L( MNSMATV_L ),   &
              MSVUNIT_L( MNSMATV_L ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MSVDESC_L,MSVUNIT_L', PROGNAME )
    CALL STORE_VDESCS( 1, 1, MNSMATV_L, .TRUE., MSVDESC_L )
    CALL STORE_VUNITS( 1, 1, MNSMATV_L, .TRUE., MSVUNIT_L )

    !.......  Open mass-based speciation matrix, compare number of sources, and store
    !           speciation variable descriptions.
    MSNAME_S = PROMPTSET( 'Enter logical name for the MASS-BASED SPECIATION MATRIX',&
                          FSREAD3, 'MSMAT_S', PROGNAME )

    IF ( .NOT. DESCSET( MSNAME_S, ALLFILES ) ) THEN
        MESG = 'Could not get description of file set "' // TRIM( MSNAME_S ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    CALL CHKSRCNO( 'mobile', 'MSMAT_S', NROWS3D, NMSRC, EFLAG)
    MNSMATV_S = NVARSET
    ALLOCATE( MSVDESC_S( MNSMATV_S ),   &
              MSVUNIT_S( MNSMATV_S ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MNSMATV_S,MSVUNIT_S', PROGNAME )
    CALL STORE_VDESCS( 1, 1, MNSMATV_S, .TRUE., MSVDESC_S )
    CALL STORE_VUNITS( 1, 1, MNSMATV_S, .TRUE., MSVUNIT_S )

    !.......  Check that variables in mole and mass speciation matrices match
    IF( MNSMATV_L .NE. MNSMATV_S ) THEN
        WRITE( MESG,94010 )                                     &
           'ERROR: Mole-based speciation matrix contains ',     &
           MNSMATV_L, 'variables but mass-based speciation ' // &
           'matrix contains ', MNSMATV_S, 'variables.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    DO I = 1, MNSMATV_L
        IF( MSVDESC_L( I ) .NE. MSVDESC_S( I ) ) THEN
            MESG = 'ERROR: Variable descriptions are not ' //   &
              'consistent between mole- and mass-based ' //     &
              'speciation matrices.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END DO

    !.......  Open meteorology file
    IF( RPDFLAG .OR. RPVFLAG .OR. RPHFLAG .OR. RPSFLAG .OR. ONIFLAG ) THEN
        IF( .NOT. ETABLEFLAG ) THEN
            METNAME = PROMPTMFILE( 'Enter logical name for the METCRO2D meteorology file',&
                                   FSREAD3, 'MET_CRO_2D', PROGNAME )

            IF( .NOT. DESC3( METNAME ) ) THEN
                MESG = 'Could not get description of file "' //&
                    METNAME( 1:LEN_TRIM( METNAME ) ) // '" '
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CALL CHKGRID( 'mobile', 'GRID', 0, EFLAG )

            !.......  Check modeling period
            CALL UPDATE_TIME_INFO( METNAME, .FALSE. )

            !.......  Make sure met file contains requested temperature variable
            J = INDEX1( TVARNAME, NVARS3D, VNAME3D )
            IF( J <= 0 ) THEN
                MESG = 'ERROR: Could not find "' // TRIM( TVARNAME ) //&
                     '" in file "' // TRIM( METNAME )
                CALL M3MESG( MESG )
            END IF

        END IF

    ELSE
        METNAME = PROMPTMFILE( 'Enter logical name for the METMOVES meteorology file',&
                               FSREAD3, 'METMOVES', PROGNAME )

        IF( .NOT. DESC3( METNAME ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( METNAME ) // '" '
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL CHKGRID( 'mobile', 'GRID', 0, EFLAG )
    END IF

    !.......  Get file name for inventory pollutants codes/names
    MESG = 'Enter logical name for INVENTORY DATA TABLE file'
    PDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INVTABLE', PROGNAME )

    !.......  Get country, state, and county names no matter what, because it is
    !           needed to allocate memory for the state and county totals, even
    !           when they aren't going to be output
    MESG = 'Enter logical name for COUNTRY, STATE, AND COUNTY file'
    CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'COSTCY', PROGNAME )

    !.......  Open source groups file if needed
    IF( SRCGRPFLAG ) THEN
        MESG = 'Enter logical name for SOURCE GROUPS file'
        SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SOURCE_GROUPS', PROGNAME )
    END IF

    !.......  Open sub-sector source groups file if needed
    IF( SUBSECFLAG ) THEN
        MESG = 'Enter logical name for SUB-SECTOR SOURCE GROUPS file'
        SGDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SUB_SEC_SOURCES', PROGNAME )
    END IF

    !.......  Get emission processes file name
    TDEV = PROMPTFFILE( 'Enter logical name for EMISSION PROCESSES file',&
                        .TRUE., .TRUE., 'MEPROC', PROGNAME )

    CALL RDEPROC( TDEV )

    !.......  Store process/pollutants combinations for correct activity
    IF( RPDFLAG ) THEN
        M = INDEX1( 'VMT', NIACT, ACTVTY )
    END IF

    IF( RPHFLAG ) THEN
        M = INDEX1( 'HOTELLING', NIACT, ACTVTY )
    END IF

    IF( RPSFLAG ) THEN
        M = INDEX1( 'STARTS', NIACT, ACTVTY )
    END IF

    IF( ONIFLAG ) THEN
        M = INDEX1( 'IDLING', NIACT, ACTVTY )
    END IF

    IF( RPVFLAG .OR. RPPFLAG ) THEN
        M = INDEX1( 'VPOP', NIACT, ACTVTY )
    END IF

    IF( M <= 0 ) THEN
        MESG = 'INTERNAL ERROR: Could not find expected activity'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MNIPPA = 0
    DO I = 1, MXETYPE
        IF( EMTNAM( I,M ) .NE. ' ' ) THEN
            MNIPPA = MNIPPA + 1
        END IF
    END DO

    ALLOCATE( MEANAM( MNIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MEANAM', PROGNAME )

    DO I = 1, MNIPPA
        MEANAM( I ) = EMTNAM( I,M )
    END DO

    !.......  Get county cross-reference file
    XDEV = PROMPTFFILE( 'Enter logical name for MCXREF cross-reference file',&
                        .TRUE., .TRUE., 'MCXREF', PROGNAME )

    !.......  Get county fuel month file
    MDEV = PROMPTFFILE( 'Enter logical name for fuel month reference file',&
                        .TRUE., .TRUE., 'MFMREF', PROGNAME )

    !.......  Get reference county emission factors file list
    FDEV = PROMPTFFILE( 'Enter logical name for reference county file list',&
                        .TRUE., .TRUE., 'MRCLIST', PROGNAME )

    !.......  Open and read hourly speed data
    IF( RPDFLAG .AND. SPDPROFLAG ) THEN
        SPDEV = PROMPTFFILE( 'Enter logical name for speed profiles file',&
                             .TRUE., .TRUE., 'SPDPRO', PROGNAME )
        CALL RDSPDPRO( SPDEV )
    END IF

    !.......  Open and read avereage speed distribution data
    IF( RPDFLAG .AND. SPDISTFLAG ) THEN
        SDDEV = PROMPTFFILE( 'Enter logical name for average speed distribution file',&
                             .TRUE., .TRUE., 'SPDIST', PROGNAME )
        CALL RDSPDIST( SDDEV )
    END IF

    !.......  Get control factor file
    IF( CFFLAG ) THEN
        CFDEV = PROMPTFFILE( 'Enter logical name for control factor file',&
                             .TRUE., .TRUE., 'CFPRO', PROGNAME )
    END IF

    !.......  If there were any errors inputing files or while comparing
    !           with one another, then abort
    IF( EFLAG ) THEN
        MESG = 'Problems opening input files. See ERROR(S) above.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  If we are using temporalized emissions, then update date/time and
    !           duration using environment variable settings, then prompt.
    CALL GETM3EPI( TZONE, SDATE, STIME, TSTEP, NSTEPS )
    TSTEP = 10000       ! only 1-hour time steps supported
    EDATE = SDATE
    ETIME = STIME
    CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP )

    !.......  Compare base year with episode and warn if not consistent
    IF( BYEAR .NE. 0 .AND. SDATE / 1000 .NE. BYEAR ) THEN

        WRITE( MESG,94010 ) 'WARNING: Inventory base year ', BYEAR,&
               'is inconsistent with year ' // CRLF() // BLANK10 //&
               'of episode start date', SDATE/1000
        CALL M3MSG2( MESG )

    ENDIF

    !.......  Give a note if running for a projected year
    IF( PYEAR .NE. BYEAR ) THEN

        WRITE( MESG,94010 ) 'NOTE: Emissions based on projected year', PYEAR
        CALL M3MSG2( MESG )

    END IF

    !.......  Reset output grid name in case meteorology files use different name
    GRDNM = OUTGRDNM

    !.......  Write message stating grid name and description
    MESG = 'NOTE: Output grid "' // TRIM( GRDNM ) //&
           '" set; described as' // CRLF() // BLANK10 // GDESC
    CALL M3MSG2( MESG )

    RETURN

94010 FORMAT( 10( A, :, I8, :, 1X ) )


    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This subprogram updates the time (episode) information
    !               and compares to the existing information, if it has been
    !               previously set.
    SUBROUTINE UPDATE_TIME_INFO( FILNAM, CHKTZONE )

        !.......  Subprogram arguments
        CHARACTER(*) FILNAM
        LOGICAL      CHKTZONE

        !.......  Local variables
        INTEGER ISECS           ! number of seconds different between dates/times
        INTEGER ED              ! tmp ending date
        INTEGER ET              ! tmp ending time
        INTEGER LOCZONE         ! tmp time zone

        !----------------------------------------------------------------------

        !.......  If time information has already been initialized...
        IF( IFLAG ) THEN
            ISECS = SECSDIFF( SDATE, STIME, SDATE3D, STIME3D )

            IF( ISECS .GT. 0 ) THEN          ! SDATE3D/STIME3D are later
                SDATE = SDATE3D
                STIME = STIME3D
            END IF

            ED = SDATE3D
            ET = STIME3D
            CALL NEXTIME( ED, ET, ( MXREC3D-1 ) * TSTEP3D )

            ISECS = SECSDIFF( EDATE, ETIME, ED, ET )

            IF( ISECS .LT. 0 ) THEN      ! ED/ET are earlier
                EDATE = ED
                ETIME = ET
            END IF

            NSTEPS = 1+ SECSDIFF( SDATE, STIME, EDATE, ETIME )/ 3600

            IF( NSTEPS .LE. 0 ) THEN
                MESG = 'Because of file ' // FILNAM //&
                       ', dates and times do not overlap at all    !'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        !.......  If time information needs to be initialized...
        ELSE
            SDATE  = SDATE3D
            STIME  = STIME3D
            NSTEPS = MXREC3D

            EDATE  = SDATE
            ETIME  = STIME
            CALL NEXTIME( EDATE, ETIME, ( NSTEPS-1 ) * TSTEP3D )

            IFLAG = .TRUE.

        END IF

        !.......  Make sure that time step is one hour
        IF( TSTEP3D .NE. 10000 ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Time step is not one hour in ' // FILNAM // ' file        !'
            CALL M3MSG2( MESG )

        END IF

        !.......  Retrieve and compare time zone
        IF( CHKTZONE ) THEN

            LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

            IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                   'Time zone ', LOCZONE, 'in ' // FILNAM //&
                   ' hourly emissions file is not consistent ' //&
                   'with initialized value of', TZONE
                CALL M3MSG2( MESG )

            ELSE IF( .NOT. ZFLAG ) THEN
                ZFLAG = .TRUE.
                TZONE = LOCZONE

                MESG = 'NOTE: Time zone initialized using ' //&
                       FILNAM // ' hourly emissions file.'

                CALL M3MSG2( MESG )
            END IF
        END IF

        !------------------  FORMAT  STATEMENTS   -----------------------------

        !.......   Internal buffering formats.......94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE UPDATE_TIME_INFO

        !----------------------------------------------------------------------
        !----------------------------------------------------------------------
        !.......  This subprogram opens the temporal emissions files. If their
        !               are multiple files, it compares the files to make sure that they
        !               are consistent with each other.  The number of sources
        !               are compared to the master number of sources.
    SUBROUTINE OPEN_TMP_FILES( LOCCAT, LBDSTAT, FNAME, SDATE )

        !.......  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: LOCCAT
        LOGICAL     , INTENT (IN) :: LBDSTAT
        CHARACTER(*), INTENT(OUT) :: FNAME( 7 )
        INTEGER     , INTENT(OUT) :: SDATE( 7 )

        !.......  Local parameters
        CHARACTER(3), PARAMETER :: SUFFIX( 7 ) =    &
            (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)

        !.......  Local allocatable arrays
        CHARACTER(NAMLEN3), ALLOCATABLE :: LOCVNAM ( : )
        CHARACTER(NAMLEN3), ALLOCATABLE :: LOCVUNIT( : )

        !.......  Local arrays
        INTEGER        IDX( 7 )             ! index for per-file arrays

        !.......  Local variables
        INTEGER        D, L, N              ! counters and indices

        INTEGER        LOCZONE           ! tmp time zone
        INTEGER        LOCNVAR           ! tmp local number of variables in file
        INTEGER        NFILE             ! no. hourly emission files

        LOGICAL     :: NFLAG = .FALSE.          ! true: no. vars inconsistent
        LOGICAL     :: VFLAG = .FALSE.          ! true: var names inconsistent
        LOGICAL     :: UFLAG = .FALSE.          ! true: var units inconsistent

        CHARACTER      CRL              ! 1-letter src category indicator
        CHARACTER(16)  TMPNAM           ! temporary logical file name
        CHARACTER(300) MESG             ! message buffer

        !----------------------------------------------------------------------

        IF( LOCCAT .EQ. 'MOBILE' ) CRL = 'M'

        !.......  Set the number of files and open the files...
        !.......  For by-day processing...
        IF( LBDSTAT ) THEN
            NFILE = 7

            DO D = 1, NFILE

                MESG = 'Enter logical name for the ' // SUFFIX( D )&
                       // ' ' // LOCCAT // ' HOURLY EMISSIONS file'
                TMPNAM = CRL // 'TMP_' // SUFFIX( D )

                FNAME( D ) = PROMPTSET( MESG,FSREAD3,&
                                          TMPNAM,PROGNAME )
                IDX( D ) = D
            END DO

        !.......  For standard processing...
        ELSE
            NFILE = 1

            MESG = 'Enter logical name for the ' // LOCCAT //&
                   ' HOURLY EMISSIONS file'
            TMPNAM = CRL // 'TMP'

            FNAME = PROMPTSET( MESG,FSREAD3,TMPNAM,PROGNAME )         ! array
            IDX( NFILE ) = 1

        END IF

        !.......  Loop through each file and ensure they are consistent
        DO D = 1, NFILE

            TMPNAM = FNAME( IDX( D ) )

            !.......  Get header and compare source number and time range
            IF ( .NOT. DESCSET( TMPNAM, ALLFILES ) ) THEN
                MESG = 'Could not get description of file set "' //&
                       TMPNAM( 1:LEN_TRIM( TMPNAM ) ) // '"'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ENDIF

            !.......  Store the starting date
            SDATE( IDX( D ) ) = SDATE3D

            !.......  Check the number of sources
            CALL CHKSRCNO( 'mobile', TMPNAM, NROWS3D,&
                           NMSRC, EFLAG )

            !.......  For standard processing, compare time info to master
            IF( .NOT. LBDSTAT .AND. D .EQ. 1 ) THEN
                CALL UPDATE_TIME_INFO( TMPNAM, .TRUE. )
            END IF

            !.......  For by-day files, make sure that the file starts at hour 0
            IF( LBDSTAT .AND. STIME3D .NE. 0 ) THEN
                EFLAG = .TRUE.
                L = LEN_TRIM( TMPNAM )
                WRITE( MESG,94010 ) 'ERROR: Start time of', STIME3D,&
                       'in file "'// TMPNAM( 1:L ) //&
                       '" is invalid.' // CRLF() // BLANK10 //&
                       'Only start time of 000000 is valid for' //&
                       'processing by day.'
                CALL M3MSG2( MESG )

            END IF

            !.......  Make sure that the file has at least 24 hours
            IF( LBDSTAT .AND. MXREC3D .LT. 24 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Number of hours', MXREC3D,      &
                       'in file "'// TRIM( TMPNAM )// '" is invalid.' //    &
                       CRLF() // BLANK10 //&
                       'Minimum number of 24 hours is needed for processing by day.'
                CALL M3MSG2( MESG )

            END IF

            LOCZONE = GETIFDSC( FDESC3D, '/TZONE/', .TRUE. )

            IF( ZFLAG .AND. LOCZONE .NE. TZONE ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                   'Time zone ', LOCZONE, 'in ' // TMPNAM //&
                   ' hourly emissions file is not consistent ' //&
                   'with initialized value of', TZONE
                CALL M3MSG2( MESG )

            ELSE IF( .NOT. ZFLAG ) THEN
                ZFLAG = .TRUE.
                TZONE = LOCZONE

                MESG = 'NOTE: Time zone initialized using ' //&
                       TMPNAM // ' hourly emissions file.'

                CALL M3MSG2( MESG )
            END IF

            !.......  For first file, store the pollutant names and units for
            !                   making comparisons with other files.
            IF( D .EQ. 1 ) THEN

                LOCNVAR = NVARS3D
                ALLOCATE( LOCVNAM( LOCNVAR ),   &
                         LOCVUNIT( LOCNVAR ), STAT=IOS )
                CALL CHECKMEM( IOS, 'LOCVUNIT', PROGNAME )

                LOCVNAM ( 1:LOCNVAR ) = VNAMESET( 1:LOCNVAR )
                LOCVUNIT( 1:LOCNVAR ) = VUNITSET( 1:LOCNVAR )

            !.......  Compare the pollutant names and units
            ELSE

                !.......  Check to make sure the number is consistent first
                IF( NVARSET .NE. LOCNVAR ) NFLAG = .TRUE.

                !.......  Make sure no overflows
                N = MIN( NVARSET, LOCNVAR )

                !.......  compare variable names and units among files
                DO V = 1, N
                    IF( LOCVNAM( V ) .NE. VNAMESET( V ) ) THEN
                        VFLAG = .TRUE.
                    END IF

                    IF( LOCVUNIT( V ) .NE. VUNITSET( V ) ) THEN
                        UFLAG = .TRUE.
                    END IF
                END DO

            END IF

        END DO

        !.......  Write message and set error if any inconsistencies
        IF( NFLAG ) THEN
        ! bbh               EFLAG = .TRUE.          ! removed to prevent false errer of odd nubmer
        !                                     of species for tmp files
            MESG = 'WARNING: ' // LOCCAT // ' source hourly ' //&
                   'emission files have inconsistent ' //&
                   CRLF() // BLANK10 // 'number of variables.'
            CALL M3MSG2( MESG )
        END IF

        IF( VFLAG ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //&
                   'emission files have inconsistent ' //&
                   CRLF() // BLANK10 // 'variable names.'
            CALL M3MSG2( MESG )
        END IF

        IF( UFLAG ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: ' // LOCCAT // ' source hourly ' //&
                   'emission files have inconsistent ' //&
                   CRLF() // BLANK10 // 'variable units.'
            CALL M3MSG2( MESG )
        END IF

        !.......  Deallocate local memory
        DEALLOCATE( LOCVNAM, LOCVUNIT )

        RETURN

        !------------------  FORMAT  STATEMENTS   -----------------------------

        !.......   Internal buffering formats.......94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE OPEN_TMP_FILES

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram initializes and checks the inventory year
    !               of the emissions and the projection status
    SUBROUTINE CHECK_INVYEAR( FNAME, PRJFLAG, IODESC )

        !.......  Subprogram arguments
        CHARACTER(*), INTENT (IN)     :: FNAME
        LOGICAL     , INTENT (IN OUT) :: PRJFLAG
        CHARACTER(*), INTENT (IN)     :: IODESC( * )

        !.......  Local variables
        INTEGER           L
        INTEGER           YY              ! tmp year
        LOGICAL           STRICT          ! flag for strict checks or not
        CHARACTER(20)     BUFFER          ! program name buffer
        INTEGER,  SAVE :: FLEN            ! name length of savnam
        CHARACTER(NAMLEN3), SAVE :: SAVNAM          ! name of file used to init

        !----------------------------------------------------------------------

        STRICT = .TRUE.

        !.......  First determine whether to abort when projected year does not
        !               match.  This is used for reactivity matrices, which will
        !               always have a projection year, even if the inventory isn't
        !               projected.
        IF( .NOT. PRJFLAG ) THEN
            BUFFER = GETCFDSC( FDESC3D, '/FROM/', .FALSE. )
            IF( BUFFER .EQ. 'OPENRMAT' ) STRICT = .FALSE.
        END IF

        !.......  If time information has already been initialized...
        IF( YFLAG ) THEN

            YY = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )
            IF( YY .LE. 0 ) THEN

                YY = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. )
                IF( YY .NE. BYEAR ) THEN
                    WRITE( MESG,94010 )&
                          'Base year of ' // FNAME // ' file:', YY,&
                          CRLF() // BLANK10 //&
                          ', does not equal emissions year of ' //&
                          SAVNAM( 1:FLEN ) // ' file:', BYEAR

                    !.......  If there is projection, abort
                    IF ( PRJFLAG ) THEN
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                    !.......  Otherwise, make it a warning
                    ELSE
                        L = LEN_TRIM( MESG )
                        MESG = 'WARNING: ' // MESG( 1:L )
                        CALL M3MSG2( MESG )
                    END IF

                END IF

            ELSE IF ( STRICT            .AND.&
                      YY     .GT. 0     .AND.&
                      YY     .NE. PYEAR      ) THEN

                WRITE( MESG,94010 )&
                      'Projected year of ' // FNAME // ' file:', YY,&
                      CRLF() // BLANK10 //&
                      ', does not equal emissions year of ' //&
                      SAVNAM( 1:FLEN ) // ' file:', PYEAR
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

        !.......  If year information needs to be initialized...
        ELSE

            BYEAR = GETIFDSC( IODESC, '/BASE YEAR/', .FALSE. )
            PYEAR = GETIFDSC( IODESC, '/PROJECTED YEAR/', .FALSE. )

            IF( PYEAR .GT. 0 ) THEN
                PRJFLAG = .TRUE.
            ELSE
                PYEAR = BYEAR
            END IF

            SAVNAM = FNAME
            FLEN   = LEN_TRIM( SAVNAM )
            YFLAG  = .TRUE.

        END IF

        !------------------  FORMAT  STATEMENTS   -----------------------------

        !.......   Internal buffering formats.......94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE CHECK_INVYEAR

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram updates the met information and compares to
    !               the existing information, if it has been previously set.
    SUBROUTINE CHECK_MET_INFO( CATDESC )

        !.......  Subprogram arguments
        CHARACTER(*) CATDESC          ! category descriptions

        !.......  Local variables
        INTEGER       L, L1, L2          ! length of strings
        CHARACTER(30) FILDESC            ! description of input file

        !----------------------------------------------------------------------

        !.......  Set tmp rows, columns, and total cells depending on file type
        IF( CATDESC .EQ. 'biogenics' ) THEN
            FILDESC = 'gridded emissions file'

        ELSEIF( CATDESC .EQ. 'mobile' ) THEN
            FILDESC = 'hourly emissions file'

        ELSEIF( CATDESC .EQ. 'point' ) THEN
            FILDESC = 'layer fractions file'

        ELSE
            MESG= 'INTERNAL ERROR: Category description "' //&
                  CATDESC// '" not known in call to CHECK_MET_INFO        !'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

        !.......  If met information has already been initialized, then compare
        !               existing to this file.
        IF( OFLAG ) THEN

            METTMP = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
            IF ( METTMP .NE. METSCENR ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Meteorology scenario name "' // TRIM( METTMP ) //    &
                      '" in ' // CATDESC //&
                       TRIM( FILDESC ) // ' is inconsistent with '//&
                       'initialized value "'// TRIM( METSCENR ) // '"'
                CALL M3MSG2( MESG )

            END IF

            METTMP = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )
            IF ( METTMP .NE. METCLOUD ) THEN

                EFLAG = .TRUE.
                MESG = 'ERROR: Meteorology cloud scheme "' // TRIM( METTMP ) // &
                       '" in ' // CATDESC //&
                       TRIM( FILDESC ) // ' is inconsistent with '//&
                       'initialized value "'// TRIM( METCLOUD ) // '"'
                CALL M3MSG2( MESG )

            END IF

        !.......  Initialize meteorology information
        ELSE

            OFLAG    = .TRUE.
            METSCENR = GETCFDSC( FDESC3D, '/MET SCENARIO/', .TRUE. )
            METCLOUD = GETCFDSC( FDESC3D, '/CLOUD SCHEME/', .TRUE. )

            MESG = 'NOTE: Meteorology description initialized '//&
                   'using '// CATDESC// ' '// TRIM( FILDESC )// '.'
            CALL M3MSG2( MESG )

        ENDIF

    END SUBROUTINE CHECK_MET_INFO

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram stores I/O API NetCDF variable descriptions into
    !               a local array based on indices in subprogram call.
    SUBROUTINE STORE_VDESCS( ISTART,INCRMT,NDESC,LFSET,DESCS )

        INCLUDE 'SETDECL.h90'           !  FileSetAPI variables and functions

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: ISTART           ! starting position in VDESCS of names
        INTEGER     , INTENT (IN) :: INCRMT           ! increment of VDESCS for names
        INTEGER     , INTENT (IN) :: NDESC            ! number of descriptions
        LOGICAL     , INTENT (IN) :: LFSET            ! number of descriptions
        CHARACTER(*), INTENT(OUT) :: DESCS( NDESC )        ! stored variable descriptions

        !.......  Local variables
        INTEGER  I, J

        !----------------------------------------------------------------------

        DESCS = ' '

        J = ISTART
        DO I = 1, NDESC

            IF( LFSET ) THEN      ! From FileSetAPI
                DESCS( I ) = TRIM( VDESCSET( J ) )
            ELSE                  ! From standard I/O API
                DESCS( I ) = TRIM( VDESC3D( J ) )
            END IF

            J = J + INCRMT

        END DO

    END SUBROUTINE STORE_VDESCS

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------
    !.......  This subprogram stores I/O API NetCDF variable units into
    !               a local array based on indices in subprogram call.
    SUBROUTINE STORE_VUNITS( ISTART,INCRMT,NUNIT,LFSET,UNITS )

        INCLUDE 'SETDECL.h90'           !  FileSetAPI variables and functions

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: ISTART                ! starting position in VDESCS of names
        INTEGER     , INTENT (IN) :: INCRMT                ! increment of VDESCS for names
        INTEGER     , INTENT (IN) :: NUNIT                 ! number of units
        LOGICAL     , INTENT (IN) :: LFSET                 ! number of descriptions
        CHARACTER(*), INTENT(OUT) :: UNITS( NUNIT )        ! stored variable units

        !.......  Local variables
        INTEGER  I, J, L

        !----------------------------------------------------------------------

        UNITS = ' '

        J = ISTART
        DO I = 1, NUNIT

            IF( LFSET ) THEN      ! From FileSetAPI
                UNITS( I ) = TRIM( VUNITSET( J ) )
            ELSE                  ! From standard I/O API
                UNITS( I ) = TRIM( UNITS3D( J ) )
            END IF

            J = J + INCRMT

        END DO

    END SUBROUTINE STORE_VUNITS

END SUBROUTINE OPENMRGIN
