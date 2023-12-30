
PROGRAM BELD3TO2

    !***********************************************************************
    !
    !  DESCRIPTION: Converts IOAPI netCDF landuse files created from SMOKE
    !               tool into BELD2 gridded ASCII format.  The SMOKE tool
    !               takes BELD3 1km resolution data and creates landuse
    !               files for user-defined domain for use with BEIS3.  BELD3
    !               databases have over 230 landuse types.  This database
    !               cannot be used for SMOKE-BEIS2 modeling since this model
    !               uses BELD2 format data that has 126 landuse types.  This
    !               program maps the 230 BELD3 types to the 126 BELD2 landuse
    !               types and outputs an ASCII format file that SMOKE-BEIS2
    !               can use.  There are 4 BELD2 categories that are output :
    !                     1: Urban Forest
    !                     2: Rural Forest
    !                     3: Agricultural
    !                     4: Other
    !
    !               SMOKE-BEIS2 treats Urban and Rural Forest the same.  So,
    !               the program outputs all forest land use types for a grid
    !               cell into the first category (Urban Forest).  The BELD3
    !               to BELD2 xref file includes the BELD2 category for each
    !               BELD3 landuse types.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !              02/01 : prototype created by J. Vukovich
    !                      tested on Lambert projection
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
    !***********************************************************************
    USE M3UTILIO

    IMPLICIT NONEf

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'         !

    !.......   PARAMETERS and their descriptions:

    !.......   LOCAL PARAMETERS

    CHARACTER(16), PARAMETER :: PROGNAME = 'BELD3TO2'        !  program name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag
    REAL,          PARAMETER :: TOHECT   = 0.0001      ! factor to convert to hectares
    CHARACTER,     PARAMETER :: AQUT     = "'"
    INTEGER,       PARAMETER :: IZERO    = 0

    !.......  EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(16),EXTERNAL :: VERCHAR
    INTEGER      ,EXTERNAL :: GETFLINE

    !.......   LOCAL VARIABLES and their descriptions:

    INTEGER         B, C, R, I, J, K, L, M, N     ! loop counters and subscripts
    INTEGER         IOS         !  I/O status result
    INTEGER         JCAT        !  temporary BELD2 category
    INTEGER         FDEV        !  unit number for BELD3 to BELD2 xref file
    INTEGER         LDEV        !  unit number for log file:
    INTEGER         XDEV        !  BELD2 output file

    CHARACTER(256)  MESG        !  message buffer for M3EXIT()
    CHARACTER(16)   LUNAMA( MXVARS3 )       ! Landuse file A var. names
    CHARACTER(16)   LUNAMB( MXVARS3 )       ! Landuse file B var. names
    CHARACTER(16), ALLOCATABLE  :: VEG3( : )         ! BELD3 landuse types
    CHARACTER(16)   LUNAM       !  temporary landuse name
    CHARACTER(16)   GRDNM       !  grid name
    CHARACTER(16)   GNAMEA      !  unit number for gridded land use file A
    CHARACTER(16)   GNAMEB      !  unit number for gridded land use file B
    CHARACTER(16)   GNAMET      !  unit number for gridded land use totals file
    CHARACTER(16)   COORD       !  coordinate system used
    CHARACTER(8)    CUNIT       !  units

    CHARACTER(4), ALLOCATABLE :: VEG2( : )          ! BELD2 xref record
    CHARACTER(4), ALLOCATABLE :: B2TYPES ( : )      ! unsorted BELD2 types
    CHARACTER(4), ALLOCATABLE :: B2SORT( : )        ! sorted BELD2 types
    CHARACTER(4)                 TMPVEG             ! temporary veg type

    INTEGER         NCOLS       ! no. of grid columns
    INTEGER         NGRID       ! no. of grid cells
    INTEGER         NROWS       ! no. of grid rows
    INTEGER         NVARSA      ! no. of variables in A file
    INTEGER         NVARSB      ! no. of variables in B file
    INTEGER         NVARST      ! total of variables in A and B
    INTEGER         IFOUND      ! temporary variable to help with indexing
    INTEGER         IXREFS      ! no. of records in xref file
    INTEGER         IB2         ! no. of BELD2 landuse types to output
    INTEGER         IUSDEC, IUSBRD, IUSCON, IUSMIX
    INTEGER         JFLAG       ! flag for counting BELD2 types

    INTEGER, ALLOCATABLE :: ICAT  ( : )          ! BELD2 category flag
    INTEGER, ALLOCATABLE :: IFLAG ( : )          ! flag for aggregration
    INTEGER, ALLOCATABLE :: IB2CAT( : )          !
    INTEGER, ALLOCATABLE :: LUINDX( : )          ! landuse name index
    INTEGER, ALLOCATABLE :: NLUS( :, :, : )      ! no. BELD2 landuse per category
    INTEGER, ALLOCATABLE :: INDXB2( : )          ! index sort of BELD2 types

    REAL, ALLOCATABLE    :: LUSE ( :, :, : )     ! BELD3 landuse
    REAL, ALLOCATABLE    :: LUB2 ( :, :, : )     ! BELD2 landuse
    REAL, ALLOCATABLE    :: TCAT ( :, :, : )     ! Total landuse per BELD2 category
    REAL, ALLOCATABLE    :: TOTAL( :, : )        ! Total landuse per grid cell
    REAL, ALLOCATABLE    :: FIA ( :, : )         ! Forest inv % in each cell
    REAL  XRES, YRES, BTMP, TTMP

    LOGICAL         EFLAG                        ! error flag

    !***********************************************************************
    !   begin body of program NORMBIO

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.

    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......   Open BELD3 to BELD2 xref file and BELD2 landuse names file

    FDEV = PROMPTFFILE( 'Enter logical name for Land Use xref file',&
                        .TRUE., .TRUE., 'B3XRF', PROGNAME )

    !.......   Loop:  read BELD3 to BELD2 xref file

    CALL M3MSG2( 'Reading BELD3to2 XREF file' )

    IXREFS = GETFLINE( FDEV, 'BELD3 to BELD2 xref file ' )

    !.......  Allocate memory

    ALLOCATE( VEG2 ( IXREFS ),      &
              VEG3 ( IXREFS ),      &
              ICAT ( IXREFS ),      &
           B2TYPES ( IXREFS ),      &
            IB2CAT ( IXREFS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VEG2...IB2CAT', PROGNAME )

    !......   Read in BELD3 to BELD2 xref and keep
    !......   track of number of BELD2 types to output.

    IB2 = 0
    JFLAG = 0

    DO I = 1, IXREFS

        READ( FDEV, 93010 ) VEG2( I ), JCAT, VEG3( I )

        !.......  Assign BELD2 category for each output type

        IF ( JCAT .GT. 4 ) THEN
            ICAT( I ) = 1
        ELSE
            ICAT( I ) = JCAT
        ENDIF

        !......  Check to see if this is new BELD2 type

        DO J = 1, IB2
            IF ( VEG2( I ) .EQ. B2TYPES( J ) ) JFLAG = 1
        ENDDO

        !......  If it is a new type assign values to appropriate arrays

        IF ( JFLAG .EQ. 0 ) THEN
            IB2 = IB2 + 1
            B2TYPES( IB2 ) = VEG2( I )
            IB2CAT( IB2 ) = ICAT( I )
        ELSE
            JFLAG = 0
        ENDIF

    ENDDO

    !....... Allocate memory and initialize

    ALLOCATE( INDXB2( IB2 ),        &
               IFLAG( IB2 ),        &
              B2SORT( IB2 ), STAT=IOS )
    CALL CHECKMEM( IOS, 'B2SORT,INDXB2', PROGNAME )
    INDXB2 = 0     ! array

    !....... Setup for alphabetically sorting BELD2 types

    DO I = 1, IB2
        B2SORT ( I ) = B2TYPES( I )
        INDXB2 ( I ) = I
    ENDDO

    !.......  Sort BELD2 types alphabetically

    CALL SORTIC ( IB2, INDXB2 , B2SORT )

    !.......   Open gridded landuse files

    GNAMET = PROMPTMFILE( 'Enter logical name for GRIDDED LANDUSE totals file',&
                          FSREAD3, 'BELD3_TOT', PROGNAME )

    IF ( .NOT. DESC3( GNAMET ) ) THEN
        MESG = 'Could not get description of"' // TRIM( GNAMET ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    !.......   Initialize grid definition

    CALL CHKGRID( GNAMET, 'GRID' , 0, EFLAG )

    GNAMEA = PROMPTMFILE( 'Enter logical name for GRIDDED LANDUSE A file',&
                          FSREAD3, 'BELD3_A', PROGNAME )

    IF ( .NOT. DESC3( GNAMEA ) ) THEN
        MESG = 'Could not get description of "' // TRIM( GNAMEA ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    !.......   Check grid definition

    CALL CHKGRID( GNAMEA, 'GRID' , 0, EFLAG )

    !......  If grid definition does not match first landuse file then stop

    IF ( EFLAG ) THEN
        MESG = 'Problems opening input files. See ERROR(S) above.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF


    NVARSA = NVARS3D

    !.......   Resolution variables assumed to be in meters

    NCOLS = NCOLS3D
    NROWS = NROWS3D
    XRES  = XCELL3D
    YRES  = YCELL3D

    DO I = 1, NVARSA
        LUNAMA ( I ) = VNAME3D ( I )
    ENDDO

    GNAMEB = PROMPTMFILE( 'Enter logical name for GRIDDED LANDUSE B file',&
                          FSREAD3, 'BELD3_B', PROGNAME )

    IF ( .NOT. DESC3( GNAMEB ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( GNAMEB ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    !.......   Check grid definition

    CALL CHKGRID( GNAMEB, 'GRID' , 0, EFLAG )

    !......  If grid definition does not match first landuse file then stop

    IF ( EFLAG ) THEN
        MESG = 'Problems opening input files. See ERROR(S) above.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    NVARSB = NVARS3D

    DO I = 1, NVARSB
        LUNAMB ( I ) = VNAME3D ( I )
    ENDDO

    !....... some grid description info obtained from call to RDSRG above

    XDEV = PROMPTFFILE( 'Enter logical name for output BELD2 LANDUSE file',&
                        .FALSE., .TRUE., 'BGUSE', PROGNAME )

    NVARST = NVARSA + NVARSB

    !...... prepare info to output for header of BELD2 output file
    !...... this header will contain grid information

    CUNIT = 'meters  '
    IF ( GDTYP3D .EQ. 1 ) THEN
        COORD = 'LAT-LON '
        CUNIT = 'degrees '
    ELSE IF ( GDTYP3D .EQ. 2 ) THEN
        COORD = 'LAMBERT '
    ELSE
        MESG = 'Code not tested for this grid type'
        CALL M3WARN( PROGNAME, 0, 0, MESG )
    ENDIF

    WRITE( XDEV, 93100 ) GDNAM3D,XORIG3D, YORIG3D, XCELL3D, YCELL3D,    &
                         NCOLS3D, NROWS3D, IZERO, COORD, CUNIT,         &
                         P_ALP3D, P_BET3D,                              &
                         P_GAM3D, XCENT3D, YCENT3D

    !.......  Allocate memory for BELD3 land use and index

    ALLOCATE(   LUSE( NCOLS, NROWS, NVARST ),        &
              LUINDX( NVARST ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LUINDX', PROGNAME )

    LUINDX = -9       ! array

    !.......  Check to see if all BELD3 types in xref are available in
    !.......  xref file

    DO I = 1, NVARST

        IFOUND = 0

        IF ( I .LE. NVARSA ) THEN
            LUNAM = LUNAMA( I )
        ELSE
            K = I - NVARSA
            LUNAM = LUNAMB( K )
        ENDIF

        DO J = 1, IXREFS

            IF ( VEG3 ( J ) .EQ. LUNAM ) THEN
                IFOUND = 1
                LUINDX( I ) = J
            ENDIF

            !.......  Assign landuse type value for 4 USGS classes to be used later

            IF ( VEG3( J ) .EQ. 'USGS_decidforest' ) IUSDEC = J
            IF ( VEG3( J ) .EQ. 'USGS_evbrdleaf  ' ) IUSBRD = J
            IF ( VEG3( J ) .EQ. 'USGS_coniferfor ' ) IUSCON = J
            IF ( VEG3( J ) .EQ. 'USGS_mxforest   ' ) IUSMIX = J

        ENDDO

        IF ( IFOUND .EQ. 0 ) THEN

            MESG =   LUNAM //  ' does NOT have xref record in B3TO2 xref file'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        ENDIF

    ENDDO


    !.......  Read the gridded landuse from the landuse files

    DO M = 1, NVARST
        N = LUINDX ( M )
        IF ( N .GT. 0 ) THEN
            IF ( M .LE. NVARSA ) THEN

                MESG = ' Reading ' // LUNAMA( M )
                CALL M3MESG( MESG )
                IF ( .NOT. READ3( GNAMEA, LUNAMA( M ) , 1,0,0 , LUSE( 1,1,N ) ) ) THEN

                    MESG = 'Could not find ' // LUNAMA( M ) // ' in file ' // GNAMEA

                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF              !  if read3() failed

            ELSE

                K = M - NVARSA
                MESG = ' Reading ' // LUNAMB( K )
                CALL M3MESG( MESG )

                IF ( .NOT. READ3( GNAMEB, LUNAMB( K ) , 1,0,0 , LUSE( 1,1,N ) ) ) THEN

                    MESG = 'Could not find ' // LUNAMB( K ) // ' in file '// GNAMEB

                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                ENDIF              !  if read3() failed

            ENDIF
        ENDIF

    ENDDO

    !.......   Allocate memory for output purposes
 
    ALLOCATE(   FIA( NCOLS, NROWS ),        &
              LUB2 ( NCOLS, NROWS, IB2 ),   &
              TCAT ( NCOLS, NROWS, 4 ),     &
              NLUS ( NCOLS, NROWS, 4 ),     &
             TOTAL ( NCOLS, NROWS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LUB2...TOTAL', PROGNAME )

    !....... Read FIA percentage from landuse totals file

    MESG = ' Reading FIA from totals file '
    CALL M3MESG( MESG )
    IF ( .NOT. READ3( GNAMET, 'FOREST', 1, 0, 0, FIA ) ) THEN
        MESG = 'Could not read FOREST from file "' // TRIM( GNAMET ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Initialize arrays

    LUB2 =  0.000     ! array
    TCAT =  0.000     ! array
    NLUS =  0         ! array
    TOTAL = 0.000     ! array

    !.......   Calculate category totals, landuse type totals for
    !.......   BELD2 output file

    DO J = 1, NROWS
    DO I = 1, NCOLS

        TTMP = 0
        IFLAG = 0               !array

        !.......  If Forest inventory used greater than 0

        IF ( FIA( I, J ) .GT. 0.000 ) THEN

            !.......   Convert from percent to hectares

            DO N = 1, IXREFS

                BTMP = LUSE( I, J, N ) * 0.01 * XRES * YRES * TOHECT

                !.......  Assume these 4 USGS classes are Scrub

                IF ( N .EQ. IUSDEC .OR. N .EQ. IUSBRD .OR.&
                     N .EQ. IUSCON .OR. N .EQ. IUSMIX )  THEN

                    TMPVEG = 'Scru'
                    K = 4
                ELSE
                    TMPVEG = VEG2( N )
                    K = ICAT( N )
                ENDIF

                IF ( BTMP .GT. 0.00049 ) THEN

                    TTMP = TTMP + BTMP

                    DO M = 1, IB2
                        L = INDXB2 ( M )
                        IF ( TMPVEG .EQ. B2SORT( L ) ) THEN

                            LUB2( I, J, L ) = LUB2( I, J, L ) + BTMP
                            TCAT( I, J, K ) = TCAT( I, J, K ) + BTMP

                            IF ( IFLAG( L ) .EQ. 0 ) THEN
                                NLUS( I, J, K ) = NLUS( I, J, K ) + 1
                                IFLAG( L ) = 1
                            ENDIF

                            TOTAL( I, J ) = TOTAL( I, J ) + BTMP
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO

        !.......  If FIA less than or equal 0 (Canada) then use USGS xrefs as is

        ELSE

            DO N = 1, IXREFS
                BTMP = LUSE( I, J, N ) * 0.01 * XRES * YRES * TOHECT

                TTMP = TTMP + BTMP

                IF ( BTMP .GT. 0.00049 ) THEN
                    K = ICAT( N )
                    TMPVEG = VEG2( N )

                    DO M = 1, IB2
                        L = INDXB2 ( M )
                        IF ( TMPVEG .EQ. B2SORT( L ) ) THEN

                            LUB2( I, J, L ) = LUB2( I, J, L ) + BTMP
                            TCAT( I, J, K ) = TCAT( I, J, K ) + BTMP

                            IF ( IFLAG( L ) .EQ. 0 ) THEN
                                NLUS( I, J, K ) = NLUS( I, J, K ) + 1
                                IFLAG( L ) = 1
                            ENDIF

                            TOTAL( I, J ) = TOTAL( I, J ) + BTMP
                        ENDIF
                    ENDDO
                ENDIF
            ENDDO
        ENDIF

    ENDDO
    ENDDO


    !.......   Write BELD2 output file

    DO J = 1, NROWS
    DO I = 1, NCOLS

        WRITE( XDEV, 93030 ) I, J , TOTAL( I, J )
        DO K = 1, 4
            WRITE( XDEV, 93040 ) TCAT( I ,J, K ), NLUS( I, J, K )
            DO M = 1, IB2
                L = INDXB2( M )
                IF ( IB2CAT( L ) .EQ. K ) THEN
                    IF ( LUB2( I, J, L ) .GT. 0.00049 ) THEN

                        WRITE( XDEV, 93050 ) AQUT, B2SORT( L ), AQUT, LUB2( I, J, L )

                    ENDIF
                ENDIF
            ENDDO
        ENDDO

    ENDDO
    ENDDO

    !.......   End of program:

    CALL M3EXIT( PROGNAME, 0, 0, 'Successful completion', 0 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Informational (LOG) message formats... 92xxx

92000 FORMAT ( 5X , A )

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )
93010 FORMAT( A4, 1x, I3, 1x , A16 )
93020 FORMAT( 1x, A4 )
93030 FORMAT( 2(I5,1x), F13.3 )
93040 FORMAT( F13.3, 1x, I3 )
93050 FORMAT( 11x, A, A4, A, F13.3 )
93100 FORMAT('#GRID ',A, 2F9.0,2F7.0,3I4, 2(1x, A), 5F5.0 )

    !.......   Internal buffering formats...... 94xxx
94010 FORMAT( 10 ( A, :, I5, :, 2X ) )

END PROGRAM  BELD3TO2

