
PROGRAM GRDMAT

    !***********************************************************************
    !  program body starts at line 148
    !
    !  DESCRIPTION:
    !     Creates the gridding matrix for any source category and creates the
    !     "ungridding" matrix for mobile sources.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
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

    !...........   MODULES for public variables
    !...........   This module is the source inventory arrays
    USE MODSOURC, ONLY: XLOCA, YLOCA, XLOC1, YLOC1, XLOC2, YLOC2,       &
                        CELLID, CIFIP, CSCC

    !...........   This module contains the cross-reference tables
    USE MODXREF, ONLY: ASRGID, SRGFIPIDX

    !.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, COORD, GDTYP,               &
                       P_ALP, P_BET, P_GAM, XCENT, YCENT, NGRID,        &
                       OFFLAG

    !.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NIPPA, EANAM

    !.........  This module is required by the FileSetAPI
    USE MODFILESET

    !...........   This module contains the gridding surrogates desciption files
    USE MODSURG, ONLY: NSRGS, SRGLIST, NTSRGDSC, SRGFNAM, NSRGFIPS,     &
                       SRGFIPS, SRGFCOD, SRGFMT, SRGNCOLS, SRGNROWS

    USE MODGRDLIB

    IMPLICIT NONE

    !...........   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL           , EXTERNAL :: DSCM3GRD
    CHARACTER(MXDLEN3), EXTERNAL :: GETCFDSC
    CHARACTER(16)     , EXTERNAL :: VERCHAR
    INTEGER           , EXTERNAL :: GETFLINE
    LOGICAL           , EXTERNAL :: BLKORCMT

    !...........   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'GRDMAT'       !  program name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !...........   LOCAL VARIABLES and their descriptions:
    !...........   Gridding Matrix

    INTEGER, ALLOCATABLE :: GMAT( : )              ! Contiguous gridding matrix

    !...........   Ungridding Matrix

    INTEGER, ALLOCATABLE :: UMAT( : )              ! Contiguous ungridding matrix

    !.........  Array that contains the names of the inventory variables needed
    !           for this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

    !...........   File units and logical/physical names
    INTEGER         IDEV        ! tmp unit number if ENAME is map file
    INTEGER         LDEV        !  log-device
    INTEGER         KDEV        !  for link defs file
    INTEGER         GDEV        !  for surrogate coeff file
    INTEGER      :: MDEV = 0    !  mobile sources codes file
    INTEGER      :: RDEV = 0    !  report file for surrogates usage
    INTEGER      :: SDEV = 0    !  ASCII part of inventory unit no.
    INTEGER         XDEV        !  for surrogate xref  file
    INTEGER         FDEV        !  for surrogate description file

    CHARACTER(NAMLEN3)   ANAME       !  logical name for ASCII inventory input file
    CHARACTER(NAMLEN3)   ENAME       !  logical name for i/o api inventory input file
    CHARACTER(NAMLEN3)   INAME       !  tmp name for inven file of unknown fmt
    CHARACTER(NAMLEN3)   GNAME       !  logical name for grid matrix output file
    CHARACTER(NAMLEN3)   UNAME       !  logical name for ungrid matrix output file

    !...........   Other local variables

    INTEGER         I, II, IT, J, L1, L2, NT, N, NN, K, KK, S     !  indices and counters.

    INTEGER         COL           ! tmp column
    INTEGER         CMAX          ! max number srcs per cell
    INTEGER         CMIN          ! min number srcs per cell
    INTEGER         ENLEN         ! length of the emissions inven name
    INTEGER         LSSC          ! previous tmp surrogate codes
    INTEGER         IOS           ! i/o status
    INTEGER         IREC          ! Record counter
    INTEGER         NK            ! Number of gridding coefficients
    INTEGER         NKU           ! Number of ungridding coefficients
    INTEGER         NINVARR       ! no. of inventory characteristics
    INTEGER      :: NMATX = 0     ! no cell-source intersections
    INTEGER      :: NMATXU= 0     ! no county-source intrsctns for all sources
    INTEGER      :: NLINES = 0    ! number of lines in input file
    INTEGER         NLIST         ! tmp surrogate list
    INTEGER         NTLINES       ! max no. of line buffers
    INTEGER         MXCCL         ! max no cells per county or link
    INTEGER         MXCSRC        ! max no cells per source
    INTEGER         MXSCEL        ! max no sources per cell
    INTEGER         ROW           ! tmp row
    INTEGER         SSC           ! surrogates code
    INTEGER         LSRGNCOLS     ! tmp surrogate colum
    INTEGER         LSRGNROWS     ! tmp surrogate row
    INTEGER         DEFSRGID      !  default surrogate ID
    INTEGER         ISDEF         !  default surrogate ID code index

    REAL            CAVG       ! average number sources per cell
    REAL*8          XX, YY

    LOGICAL      :: A2PFLAG = .FALSE.      ! true: inv has ar-to-pt locations
    LOGICAL      :: AFLAG   = .FALSE.      ! true: use grid adjustments file
    LOGICAL      :: DFLAG   = .FALSE.      ! true: use link defs file
    LOGICAL      :: EFLAG   = .FALSE.      ! true: error found
    LOGICAL      :: SRGFLAG = .FALSE.      ! true: surrogates are needed
    LOGICAL      :: NFLAG   = .FALSE.      ! true: process gridded I/O API inv file
    LOGICAL      :: UFLAG   = .FALSE.      ! true: create ungridding matrix
    LOGICAL      :: VFLAG   = .FALSE.      ! true: use variable grid
    LOGICAL      :: FSGFLAG = .FALSE.      ! true: use default fallback surrogate

    !...........   Local parameters
    CHARACTER(SCCLEN3)       PSCC                 !  tmp and previous SCC
    CHARACTER(NAMLEN3)       COORUNIT             !  coordinate system projection units
    CHARACTER(NAMLEN3)    :: INVGRDNM  = ' '      !  inventory grid name
    CHARACTER(NAMLEN3)    :: SRGGRDNM  = ' '      !  surrogates file grid name
    CHARACTER(NAMLEN3)    :: LSRGGRDNM = ' '      !  surrogates file grid name
    CHARACTER(NAMLEN3)       GRDSRG               !  reading type of surrogate files
    CHARACTER(MXDLEN3)       GDESC                !  grid description
    CHARACTER(196)           NAMBUF               !  surrogate file name buffer
    CHARACTER(256)           NAMBUFT              !  tmp surrogate file name buffer
    CHARACTER(256)           TSRGFNAM             !  tmp surrogate file name buffer
    CHARACTER(300)           MESG                 !  message buffer
    CHARACTER(80)            LINE                 ! Read buffer for a line

    CHARACTER(MXDLEN3)  IFDESC2, IFDESC3     !  fields 2  3 from PNTS FDESC

    !***********************************************************************
    !   begin body of program GRDMAT

    LDEV = INIT3()

    !.........  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.........  Get environment variables that control this program
    AFLAG = ENVYN( 'GRDMAT_ADJUST',&
                   'Use grid adjustments file or not',&
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "GRDMAT_ADJUST"', 2 )
    END IF

    DFLAG = ENVYN( 'GRDMAT_LINKDEFS',&
                   'Use link definitions file or not',&
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "GRDMAT_LINKDEFS"', 2 )
    END IF

    VFLAG = ENVYN( 'USE_VARIABLE_GRID',&
                   'Use variable grid definition',&
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

    NFLAG = ENVYN( 'IMPORT_GRDNETCDF_YN',&
                   'Process native NetCDF gridded inventory file',&
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "IMPORT_GRDNETCDF_YN"', 2 )
    END IF

    !.........  Temporary section for disallowing optional files
    IF( AFLAG ) THEN
        MESG = 'NOTE: Grid adjustments file is not supported yet'
        CALL M3MSG2( MESG )
        AFLAG = .FALSE.
    END IF

    IF( DFLAG ) THEN
        MESG = 'NOTE: Link definitions file is not supported yet'
        CALL M3MSG2( MESG )
        DFLAG = .FALSE.
    END IF

    !.........  Set source category based on environment variable setting
    CALL GETCTGRY

    !.........  Get inventory file names given source category
    CALL GETINAME( CATEGORY, ENAME, ANAME )

    !.........  Prompt for and open inventory file
    INAME = ENAME
    MESG = 'Enter logical name for the MAP INVENTORY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.........  Open and read map file
    CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    !        IF( AFLAG )
    !       ADEV = PROMPTFFILE(
    !                'Enter logical name for ADJUSTMENT FACTORS file',
    !                .TRUE., .TRUE., CRLF // 'ADJUST', PROGNAME )

    !.............  Store source specific information based on header
    CALL GETSINFO( ENAME )

    IFDESC2  = GETCFDSC( FDESC3D, '/FROM/',     .TRUE. )
    IFDESC3  = GETCFDSC( FDESC3D, '/VERSION/',  .TRUE. )
    INVGRDNM = GETCFDSC( FDESC3D, '/GRIDNAME/', .FALSE. )

    !.........  Set inventory variables to read for all source categories
    IVARNAMS( 1 ) = 'CIFIP'

    SELECT CASE ( CATEGORY )

      CASE ( 'AREA' )
        NINVARR = 4
        IVARNAMS( 2 ) = 'CSCC'
        IVARNAMS( 3 ) = 'CELLID'
        IVARNAMS( 4 ) = 'CSOURC'

    !............  Check to see if point locations are in the AREA
    !              file (i.e. we have area-to-point sources)
        K = INDEX1( 'XLOCA', NVARSET, VNAMESET )
        IF ( K .GT. 0 ) THEN

    !.................  Make sure we're not using a variable grid
            IF( VFLAG ) THEN
                MESG = 'Cannot use area-to-point sources ' //&
                       'with a variable grid.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            A2PFLAG = .TRUE.
            NINVARR = 6
            IVARNAMS( 5 ) = 'XLOCA'
            IVARNAMS( 6 ) = 'YLOCA'
        END IF

      CASE ( 'MOBILE' )
        NINVARR = 10
        IVARNAMS( 2 ) = 'IRCLAS'
        IVARNAMS( 3 ) = 'IVTYPE'
        IVARNAMS( 4 ) = 'CSOURC'
        IVARNAMS( 5 ) = 'CSCC'
        IVARNAMS( 6 ) = 'CLINK'
        IVARNAMS( 7 ) = 'XLOC1'
        IVARNAMS( 8 ) = 'YLOC1'
        IVARNAMS( 9 ) = 'XLOC2'
        IVARNAMS( 10 ) = 'YLOC2'

      CASE ( 'POINT' )
        NINVARR = 3
        IVARNAMS( 2 ) = 'XLOCA'
        IVARNAMS( 3 ) = 'YLOCA'

    END SELECT

    !.........  Allocate memory for and read in required inventory characteristics
    CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

    !.........  Define source-category-specific settings
    SELECT CASE( CATEGORY )

      CASE( 'AREA' )

    !.............  Determine if surrogates are needed by checking whether all cell
    !               values are defined or not
        DO S = 1, NSRC
            IF( CELLID( S ) .LE. 0 ) SRGFLAG = .TRUE.
            IF( SRGFLAG ) EXIT
        END DO

      CASE( 'MOBILE' )

    !.............  Check for link-based sources when using a variable grid
        IF( MAXVAL( XLOC1 ) .GT. AMISS3 .AND. VFLAG ) THEN
            MESG = 'Cannot use link-based data with a variable grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    !.............  If all sources are link-based, don't need surrogates (needs to be tested)
    !            IF( MINVAL( XLOC1 ) .GT. AMISS3 ) THEN
    !                SRGFLAG = .FALSE.
    !            ELSE
        SRGFLAG = .TRUE.
    !            END IF

      CASE( 'POINT' )

    END SELECT

    !.........  Get gridding surrogates and x-ref file, if needed
    IF( SRGFLAG ) THEN

        XDEV = PROMPTFFILE( 'Enter logical name for GRIDDING SURROGATE XREF file',&
                            .TRUE., .TRUE., CRL // 'GREF', PROGNAME )
    END IF

    !.........  Get mobile-specific files
    IF( CATEGORY .EQ. 'MOBILE' ) THEN

        IF( DFLAG ) KDEV = PROMPTFFILE( 'Enter logical name for LINK DEFINITIONS file',&
                                        .TRUE., .TRUE., CRL // 'GLNK', PROGNAME )

    END IF      ! End of mobile file opening

    !.........  If surrogates are in use, read the gridding cross-reference,
    !           and the header of the surrogates file to get the grid on which
    !           the surrogates are available.
    IF( SRGFLAG ) THEN

        !.............  Build unique lists of SCCs and country/state/county codes
        !               from the inventory arrays
        CALL GENUSLST

        !..............  For mobile sources, read the mobile codes
        CALL M3MSG2( 'Reading gridding cross-reference file...' )

        !.............  Read the gridding cross-reference
        CALL RDGREF( XDEV )
        CALL M3MSG2( 'Reading gridding surrogates header...' )

        !.............  Read surrogate description file
        MESG = 'Enter logical name for surrogate description file'
        FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'SRGDESC', PROGNAME )

        CALL M3MSG2( 'Reading gridding surrogate description file... ' )
        CALL RDSRGDESC( FDEV )

        !..............  Read the link definition file
        !            CALL RDLNKDEF( )

        !.............  Determine default surrogate number from by environment variable
        !.............  Default surrogate code 100 is population

        FSGFLAG  = ENVYN( 'SMK_USE_FALLBACK',       &
                          'Using default fallback surrogate', .FALSE., IOS )

        DEFSRGID = ENVINT( 'SMK_DEFAULT_SRGID', 'Default surrogate', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_DEFAULT_SRGID"', 2 )
        END IF

        ISDEF = FIND1( DEFSRGID, NSRGS, SRGLIST )

        IF( ISDEF .LT. 1 ) THEN
            WRITE( MESG, 93000 ) 'ERROR: Default surrogate code' //     &
                ' is not available in your surrogate description '//    &
                 CRLF()//BLANK10//': User MUST choose a population'//   &
                ' surrogate code (i.e.,100) as default for a proper'//  &
                ' processing'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !..............  Allocate memory for indices to surrogates tables for each source
        ALLOCATE( ASRGID( NSRC ),       &
               SRGFIPIDX( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ASRGID,SRGFIPIDX', PROGNAME )

        !.............  Assigns the index of the surrogate to each source (stored
        !           in SRGIDPOS passed through MODXREF)
        CALL ASGNSURG

        !............  Sort surrogates by county code  cell  surrogate code
        CALL SORTI1( NSRC, SRGFIPIDX, ASRGID )

        !.............  Creating surrogate list of assigned surrogate codes
        !.............  Count surrogate codes in surrogates file
        LSSC  = -1
        NSRGS = 0
        DO I = 1, NSRC

            J   = SRGFIPIDX( I )
            SSC = ASRGID( J )

            IF( SSC .NE. LSSC ) THEN
                NSRGS = NSRGS + 1
                LSSC = SSC
            END IF

        END DO

        DEALLOCATE( SRGLIST )

        !.............  Allocate memory for creating assigned surrogates tables
        ALLOCATE( SRGLIST( NSRGS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SRGLIST', PROGNAME )

        !.............  Initialize arrays
        SRGLIST = 0          ! array

        !.............  Store the surrogate ID list
        LSSC  = -1
        NN = 0
        DO I = 1, NSRC

            J   = SRGFIPIDX( I )
            SSC = ASRGID( J )

            IF( SSC .NE. LSSC ) THEN
                NN = NN + 1
                SRGLIST( NN ) = SSC
                LSSC = SSC
            END IF

        END DO

        !.............  Re-create SRGLIST including default surrogate
        !.............  if it doesn't exit in a new assigned SRGLIST
        ISDEF = FIND1( DEFSRGID, NSRGS, SRGLIST )

        IF( ISDEF < 1 ) THEN
            NSRGS = NSRGS + 1

            !.................  Allocate memory for creating assigned surrogates tables
            DEALLOCATE ( SRGLIST )
            ALLOCATE( SRGLIST( NSRGS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRGLIST', PROGNAME )

            !.................  Initialize arrays
            SRGLIST = 0              ! array

            !.................  Store the surrogate ID list
            LSSC  = -1
            NN = 0
            DO I = 1, NSRC

                J   = SRGFIPIDX( I )
                SSC = ASRGID( J )

                IF( SSC .NE. LSSC ) THEN

                    IF( LSSC < DEFSRGID .AND. DEFSRGID < SSC ) THEN
                        NN = NN + 1
                        SRGLIST( NN ) = DEFSRGID
                        LSSC = DEFSRGID
                    END IF

                    NN = NN + 1
                    SRGLIST( NN ) = SSC
                    LSSC = SSC

                END IF

            END DO

            !.................  When default surrogate is greater than any other srgs
            IF( DEFSRGID > SRGLIST( NN ) ) THEN
                NN = NN + 1
                SRGLIST( NN ) = DEFSRGID
            END IF

        END IF

        ISDEF = FIND1( DEFSRGID, NSRGS, SRGLIST )

        !.............  Read the surrogates header and initialize the grid description
        !.............  Also, obtain the format of the file.
        !.............  Save the name of the input grid
        !.............  Ensure the grid name from the surrogate description files are consistent
        LSRGNCOLS = -1
        LSRGNROWS = -1
        LSRGGRDNM = ''
        KK = 0

        DO I = 1, NTSRGDSC      ! Open all surrogate files using the same srg code

            !.................  Prompt for and open surrogate file(s)...
            CALL GETENV( 'SRGPRO_PATH', NAMBUF )
            WRITE( NAMBUFT, '( 2A )' ) TRIM( NAMBUF ), SRGFNAM( I )

            IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                CALL OPEN_SRGFILE

                CALL RDSRGHDR( VFLAG, GDEV, SRGFMT )
                SRGGRDNM = GRDNM
                SRGNCOLS = NCOLS
                SRGNROWS = NROWS

                IF( SRGNCOLS .NE. LSRGNCOLS .OR.&
                    SRGNROWS .NE. LSRGNROWS .OR.&
                    SRGGRDNM .NE. LSRGGRDNM ) THEN

                    LSRGGRDNM = SRGGRDNM
                    LSRGNCOLS = SRGNCOLS
                    LSRGNROWS = SRGNROWS
                    KK = KK + 1

                END IF

                IF( KK .GT. 1 ) THEN
                    MESG = '#GRID Description from surrogate files'&
                           // ' are not consistent'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                CLOSE( GDEV )

            END IF

            TSRGFNAM = NAMBUFT        ! store a previous surrogate file name buffer

        END DO

    END IF

    !.........  Get grid name from the environment and read grid parameters
    IF ( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD, GDTYP3D,    &
         COORUNIT, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D,      &
         YCENT3D, XORIG3D, YORIG3D, XCELL3D,                &
         YCELL3D, NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Check or initialize the output grid settings (depends on
    !           if a surrogates file is being used). For variable grids,
    !           do not allow subgrids.
    IF( VFLAG ) THEN
        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 0, EFLAG )

        IF( EFLAG ) THEN
            MESG = 'Problem with variable grid input data'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    ELSE
        CALL CHKGRID( GDNAM3D, 'GRIDDESC', 1, EFLAG )

        IF ( EFLAG ) THEN
            MESG = 'Problem with gridded input data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    !.........  Ensure that the output grid is consistent with the input grid
    IF ( INVGRDNM .NE. ' ' .AND.    &
         INVGRDNM .NE. GRDNM     ) THEN

        MESG = 'Output grid "' // TRIM( GRDNM ) //                  &
               '" is inconsistent with inventory input grid "'//    &
               TRIM( INVGRDNM )// '".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.........  Write message stating grid name and description
    MESG = 'NOTE: Output grid "' // TRIM( GRDNM ) //    &
           '" set; described as' // CRLF() // BLANK10 // GDESC
    CALL M3MSG2( MESG )

    !.........  If output grid is different from surrogates, write message
    IF ( OFFLAG ) THEN
        MESG = 'NOTE: gridding surrogates extracted for output '//&
               'grid from grid "' // TRIM( SRGGRDNM ) // '"'
        CALL M3MSG2( MESG )
    END IF

    !.........  If area or mobile inventory has point source locations,
    !           convert point source coordinates from lat-lon to output grid
    IF( A2PFLAG ) THEN
        CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,&
                       XCENT, YCENT, XLOCA, YLOCA )
    END IF

    !.........  Depending on source category, convert coordinates, determine size
    !           of gridding matrix, and allocate gridding matrix.
    SELECT CASE( CATEGORY )

      CASE( 'AREA' )
        !.............  Determine sizes for allocating area gridding matrix
        IF( SRGFLAG ) THEN
            CALL SIZGMAT( CATEGORY, NSRC, VFLAG, DEFSRGID, FSGFLAG,&
                          MXSCEL, MXCSRC, MXCCL, NMATX, NMATXU)

        ELSE          ! processing pregridded IOAPI file

            MXSCEL = 1
            MXCSRC = 1
            MXCCL  = 1
            NMATX  = NSRC
            NMATXU = NSRC

        ENDIF

        !.............  Allocate memory for mobile source gridding matrix
        IF( .NOT. NFLAG ) THEN
            ALLOCATE( GMAT( NGRID + 2*NMATX ), STAT=IOS )
            CALL CHECKMEM( IOS, 'GMAT', PROGNAME )
        END IF

      CASE( 'MOBILE' )

        !.............  Convert mobile source coordinates from lat-lon to output grid
        CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,&
                       XCENT, YCENT, XLOC1, YLOC1 )
        CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,&
                       XCENT, YCENT, XLOC2, YLOC2 )

        !.............  Determine sizes for allocating mobile gridding matrix
        CALL SIZGMAT( CATEGORY, NSRC, VFLAG, DEFSRGID, FSGFLAG,&
                      MXSCEL, MXCSRC, MXCCL, NMATX, NMATXU)

        !.............  Allocate memory for mobile source gridding matrix
        ALLOCATE( GMAT( NGRID + 2*NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

        !.............  Allocate memory for mobile source ungridding matrix
        IF( .NOT. UFLAG ) NMATXU = 1
        ALLOCATE( UMAT( NSRC + 2*NMATXU ), STAT=IOS )
        CALL CHECKMEM( IOS, 'UMAT', PROGNAME )

      CASE( 'POINT' )

        !.............  Convert point source coordinates from lat-lon to output grid
        CALL CONVRTXY( NSRC, GDTYP, GRDNM, P_ALP, P_BET, P_GAM,&
                       XCENT, YCENT, XLOCA, YLOCA )

        !.............  Set the number of source-cell intersections
        !.............  Use double-precision INGRID to minimize round-off errors
        DO S = 1, NSRC
            IF( INGRID( XLOCA( S ), YLOCA( S ), NCOLS, NROWS, COL, ROW  ) ) THEN
                NMATX = NMATX + 1
            END IF
        END DO

        !.............  Allocate memory for point source gridding matrix
        ALLOCATE( GMAT( NGRID + NMATX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

    END SELECT

    !.........  Abort of there are no source-cell intersections
    IF( NMATX .EQ. 0 ) THEN
        MESG = 'No source-cell intersections found.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Get file names; open output gridding matrix (and ungridding matrix
    !           for mobile) using grid characteristics from DSCM3GRD() above
    !.........  Also open report file
    CALL OPENGMAT( NMATX, NMATXU, UFLAG, IFDESC2, IFDESC3, VFLAG,&
                   NFLAG, GNAME, UNAME, RDEV )

    CALL M3MSG2( 'Generating gridding matrix...' )

    !.........  Generate gridding matrix for given source category, and write it
    !           out.  It is necessary to write it out while in the subroutine,
    !           because of the type transformation from real to integer that
    !           is done so the sparse i/o api format can be used.

    SELECT CASE( CATEGORY )

      CASE( 'AREA' )

        IF( SRGFLAG ) THEN
            CALL GENAGMAT( GNAME, RDEV, MXSCEL, NSRC, NMATX, VFLAG,     &
                       DEFSRGID, FSGFLAG, GMAT( 1 ),GMAT( NGRID+1 ),    &
                       GMAT( NGRID+NMATX+1 ), NK, CMAX, CMIN )

        ELSE IF( NFLAG ) THEN        ! processing pregridded Raw NetCDF file
            CALL GENLGMAT( GNAME, RDEV, NK, CMAX, CMIN )

        ELSE         ! processing pregridded IOAPI file
            CALL GENGGMAT( GNAME, RDEV, MXSCEL, NSRC, NMATX,            &
                          GMAT( 1 ),GMAT( NGRID+1 ),                    &
                          GMAT( NGRID+NMATX+1 ), NK, CMAX, CMIN )

        ENDIF

      CASE( 'MOBILE' )

        CALL GENMGMAT( ENAME, GNAME, UNAME, RDEV, MXSCEL, MXCSRC,       &
                       MXCCL, NSRC, NMATX, NMATXU, UFLAG, VFLAG,        &
                       DEFSRGID, FSGFLAG, GMAT( 1 ),GMAT( NGRID+1 ),    &
                       GMAT( NGRID+NMATX+1 ), UMAT( 1 ),                &
                       UMAT( NSRC+1 ), UMAT( NSRC+NMATXU+1 ), NK,       &
                       CMAX, CMIN, NKU )

      CASE( 'POINT' )

        CALL GENPGMAT( GNAME, NSRC, NGRID, XLOCA, YLOCA, VFLAG,         &
                       GMAT( 1 ), GMAT( NGRID+1 ), NK, CMAX, CMIN )

    END SELECT

    !.........  Report statistics for gridding matrix

    CAVG = FLOAT( NK ) / FLOAT( NGRID )
    CALL M3MSG2( 'GRIDDING-MATRIX statistics:' )

    WRITE( MESG,94010 ) 'Total number of coefficients    :',        NK,     &
           CRLF() // BLANK5 // 'Max  number of sources per cell :', CMAX,   &
           CRLF() // BLANK5 // 'Min  number sources per cell > 0:', CMIN

    L1 = LEN_TRIM( MESG )
    WRITE( MESG,94020 ) MESG( 1:L1 ) // CRLF() // BLANK5 // &
           'Mean number of sources per cell :', CAVG

    CALL M3MSG2( MESG )

    !.........  Report statistics for ungridding matrix
    IF( UFLAG ) THEN

        CAVG = FLOAT( NKU ) / FLOAT( NSRC )
        CALL M3MSG2( 'UNGRIDDING-MATRIX statistics:' )

        WRITE( MESG,94010 ) 'Total number of coefficients    :', NKU

        WRITE( MESG, 94020 ) TRIM( MESG ) // CRLF() // BLANK5 //    &
          'Mean number of cells per source :', CAVG

        CALL M3MSG2( MESG )

    END IF

    !.........  End of program

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Informational (LOG) message formats... 92xxx

92000 FORMAT( 5X, A )

92010 FORMAT( 5X, A, :, I12 )

92020 FORMAT( 5X, A, :, F17.4 )


    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )
93001 FORMAT( 2A )

93010 FORMAT( A16 )


    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

94020 FORMAT( A, :, F10.2 )

    !******************  INTERNAL SUBPROGRAMS  *****************************
CONTAINS

    !.........  This internal subprogram opens individual surrogate file

    SUBROUTINE OPEN_SRGFILE
        !----------------------------------------------------------------------

        !.........  Set logical file name
        IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUFT )) THEN
            MESG = 'Could not set logical file name of ' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.........  Get the number of lines in the surrogate description file desription file
        GDEV = GETEFILE( 'SRGPRO_PATH',.TRUE., .TRUE., PROGNAME )

        IF( GDEV .LT. 0 ) THEN
            MESG = 'Could not open input surrogate file' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        REWIND( GDEV )

        NLINES = GETFLINE( GDEV, 'Reading surrogate files' )

        IF( .NOT. SETENVVAR( 'SRGPRO_PATH', NAMBUF )) THEN
            MESG = 'Could not set logical file name of ' // TRIM( NAMBUF )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

    !------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

    END SUBROUTINE OPEN_SRGFILE

END PROGRAM GRDMAT
