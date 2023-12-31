
SUBROUTINE NORMBEIS370( CVSW )

    !***********************************************************************
    !
    !  DESCRIPTION:  Produces normalized biogenic emissions for use with
    !        SMOKE-BEIS3 v3.12.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED: RDB3FAC
    !
    !  REVISION  HISTORY:
    !       4/00 Prototype, Jeff Vukovich
    !       1/03 changes to NO, George Pouliot
    !       4/12 update to BELD4 using comma delmited emission factor file
    !       8/2012  update for BELD4 and soil NOx for irrigated crop lands
    !       8/2020  update for BELD5 and BEISFAC for BEISv3.7
    !       3/2022  fix bugs for landuse names (MODIS_14,
    !        Potatoes and Other_Crops), Vukovich
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  
    !       bug-fix for bogus "missing: test, and related changes
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

    !......  MODULES for public variables
    !......  This module contains biogenic variables
    USE MODBEIS3, ONLY: NVEG, VEGID, AVGEMIS, AVGLAI, NOEMIS,    &
                        EMFAC, LAI, SLW, WFAC, LFBIO

    USE MODGRDLIB

    IMPLICIT NONE

    !......  INCLUDES
    INCLUDE 'B3V14DIMS3.h90'         ! BEIS3-related declarations

    !......  ARGUMENTS and their descriptions
    CHARACTER(50), INTENT(IN) :: CVSW        ! CVS release tag

    !......  EXTERNAL FUNCTIONS and their descriptions
    INTEGER      , EXTERNAL :: GETFLINE
    CHARACTER(16), EXTERNAL :: VERCHAR

    !......  PARAMETERs:
    CHARACTER(16) :: PNAME = 'NORMBEIS370'      ! Program name

    !......  LOCAL VARIABLES and their descriptions
    INTEGER         B, C, R, I, J, K, L, M, N     ! loop counters and subscripts
    INTEGER         IOS         !  I/O status result
    INTEGER         FDEV        !  unit number for emissions factor file
    INTEGER         LDEV        !  unit number for log file

    CHARACTER(16)   ENAME       !  logical name for normalized emissions output
    CHARACTER(16)   DOTNAME     !  logical name for gridded 2D dot file
    CHARACTER(16)   GRDNM       !  grid name
    CHARACTER(16), ALLOCATABLE  :: LUNAM( : )      ! land use type names
    CHARACTER(16), ALLOCATABLE  :: VGLIST( : )      ! land use type names
    CHARACTER(16)   LUNAM_TMP       ! temporary tag for land use names
    CHARACTER(16)   GNAMET      !  unit number for gridded land use totals file

    CHARACTER(256)  MESG        !  message buffer for M3EXIT()
    CHARACTER(5)    BTMP        ! temporary tag for naming output variables

    INTEGER, ALLOCATABLE      :: LUINDX( : )       ! index for land use types
    INTEGER         NCOLS       ! no. of grid columns
    INTEGER         NROWS       ! no. of grid rows
    INTEGER         NCDOT       ! no. of columns in dot file
    INTEGER         NRDOT       ! no. of rows in dot file

    INTEGER         NVARS4      ! total number of variables in input file > 120
    INTEGER         IFOUND      ! used for checking land use vs. emis facs
    INTEGER         IUSDEC      ! USGS decid forest
    INTEGER         IUSBRD      ! USGS evbrdleaf
    INTEGER         IUSCON      ! USGS coniferfor
    INTEGER         IUSMIX      ! USGS mixed forest
    INTEGER         IUSSHR      ! USGS shrubland
    INTEGER         IUSCGS      ! USGS Cropgrass
    INTEGER         IUSCWD      ! USGS Cropwoodland
    INTEGER         IUSCDY      ! USGS Drycrop
    INTEGER         IUSCIR      ! USGS Irrcrop
    INTEGER         IALFAL      ! Alfalfa
    INTEGER         IBARLE      ! Barley
    INTEGER         ICORNG       ! Corn
    INTEGER         ICORNS       ! Corn
    INTEGER         ICOTTO      ! Cotton
    INTEGER         IGRASS      ! Grass
    INTEGER         IHAY        ! Hay
    INTEGER         IMISCC      ! Misc_crop
    INTEGER         IOATS       ! Oats
    INTEGER         IPEANU      ! Peanuts
    INTEGER         IPOTAT      ! Potatoes
    INTEGER         IRICE       ! Rice
    INTEGER         IRYE        ! Rye
    INTEGER         ISORGHG      ! Sorghum
    INTEGER         ISORGHS      ! Sorghum
    INTEGER         ISOYBE      ! Soybeans
    INTEGER         IBEANSED     ! Edible Benas
    INTEGER         IBEANS      ! Beans
    INTEGER         ICANOLA     ! Canola
    INTEGER         IWHEATW      ! Wheat
    INTEGER         IWHEATS      ! Wheat
    
    INTEGER         IALFAL_IR      ! Alfalfa
    INTEGER         IBARLE_IR      ! Barley
    INTEGER         ICORNG_IR       ! Corn
    INTEGER         ICORNS_IR       ! Corn
    INTEGER         ICOTTO_IR      ! Cotton
    INTEGER         IGRASS_IR      ! Grass
    INTEGER         IHAY_IR        ! Hay
    INTEGER         IMISCC_IR      ! Misc_crop
    INTEGER         IOATS_IR       ! Oats
    INTEGER         IPEANU_IR      ! Peanuts
    INTEGER         IPOTAT_IR      ! Potatoes
    INTEGER         IRICE_IR       ! Rice
    INTEGER         IRYE_IR        ! Rye
    INTEGER         ISORGHG_IR      ! Sorghum
    INTEGER         ISORGHS_IR      ! Sorghum
    INTEGER         ISOYBE_IR      ! Soybeans
    INTEGER         IBEANSED_IR     ! Edible Benas
    INTEGER         IBEANS_IR      ! Beans
    INTEGER         ICANOLA_IR     ! Canola
    INTEGER         IWHEATS_IR      ! Wheat Spring Irrigated
    INTEGER         IWHEATW_IR      ! Wheat Winter Irrigated

    INTEGER         NLCD82      ! Cultivated Crops
    INTEGER         NLCD81      ! Pasture Hay Grass

    INTEGER         MODIS14      ! modis mosaic 1/3 (grass + mxforest + drycrop)
    INTEGER         MODIS12      ! irrigated crops
    INTEGER         LAI_SAVE_INDEX(3)
    INTEGER         NLINES

    REAL,   ALLOCATABLE ::    XVALS ( :,: )         ! x values for grid cell boundaries
    REAL,   ALLOCATABLE ::    YVALS ( :,: )         ! y values for grid cell boundaries
    REAL,   ALLOCATABLE ::     LUSE ( :,:,:  )      ! BELD4 land use data
    REAL,   ALLOCATABLE ::     SUMEM( : )           ! Summer emissions
    REAL,   ALLOCATABLE ::    SUMEMW( : )           ! Winter emissions
    REAL,   ALLOCATABLE ::      NOEM( : )           ! NO emissions
    REAL,   ALLOCATABLE ::    SUMLAI( : )           ! Summer LAIs
    REAL,   ALLOCATABLE ::   SUMLAIW( : )           ! Winter LAIs
    REAL*8, ALLOCATABLE :: PRCNT2KM2( :,: )         ! Prcnt to km**2

    REAL        VEGAREA             ! Veg. area for
    REAL        EFTMP               ! Emission factors
    REAL        XCELL               ! x cell size
    REAL        YCELL               ! y cell size

    LOGICAL     USE_SHRUB
    LOGICAL     VFLAG               ! variable grid flag
    LOGICAL     EFLAG               ! Error flag

    !***********************************************************************
    !   begin body of subroutine NORMBEIS370

    EFLAG = .FALSE.

    !......  Check for variable grid data
    VFLAG = ENVYN( 'USE_VARIABLE_GRID',    &
                   'Use variable grid definition',    &
                   .FALSE., IOS )
    IF ( IOS .NE. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

    !......  Get file name; open emission factors file
    FDEV = PROMPTFFILE(    &
             'Enter logical name for EMISSION FACTORS file',    &
             .TRUE., .TRUE., 'BEISFAC', PNAME )
    IF ( FDEV .LT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Could not open "B3FAC"', 2 )
    END IF

    !......  Open gridded landuse files
    GNAMET = PROMPTMFILE(    &
             'Enter logical name for GRIDDED LANDUSE totals file',    &
             FSREAD3, 'BELD5', PNAME )

    IF ( .NOT. DESC3( GNAMET ) ) THEN
        MESG = 'Could not get description of file "' //    &
               TRIM( GNAMET ) // '"'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF
    NVARS4 = NVARS3D

    !......  Initialize grid definition
    CALL CHKGRID( GNAMET, 'GRID' , 0, EFLAG )

    !......  Grid cell resolution assumed to be in meters
    !        Input land use will be in percentages
    !        Compute conversion factor to go from percentages
    !        to km**2
    NCOLS = NCOLS3D
    NROWS = NROWS3D
    GRDNM = GDNAM3D

    IF ( VFLAG ) THEN
        !......  Open GRIDDOT2D file
        DOTNAME = PROMPTMFILE(    &
         'Enter logical name for DOT-POINT SURFACE GRID file',    &
         FSREAD3, 'GRID_DOT_2D', PNAME )

        If ( .NOT. DESC3( DOTNAME ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( DOTNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        ENDIF

        !......  Check grid definition
        CALL CHKGRID( DOTNAME, 'DOT', 0, EFLAG )

        IF ( EFLAG ) THEN
            MESG = 'Grid in file "' // TRIM( DOTNAME ) //    &
                   '" does not match previously set grid.'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        NCDOT = NCOLS + 1
        NRDOT = NROWS + 1

        !......  Allocate memory for grid cell coordinates
        ALLOCATE( XVALS( NCDOT, NRDOT ),    &
                  YVALS( NCDOT, NRDOT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'XVALS,YVALS', PNAME )

        !......  Read grid cell coordinates
        IF( .NOT. READ3( DOTNAME, 'LON', 1, 0, 0, XVALS ) ) THEN
            MESG = 'Could not read LON from file "' //    &
                   TRIM( DOTNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. READ3( DOTNAME, 'LAT', 1, 0, 0, YVALS ) ) THEN
            MESG = 'Could not read LAT from file "' //    &
                   TRIM( DOTNAME ) // '"'
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        !......  Convert coordinates to map projection units
        CALL CONVRTXY( NCDOT, NRDOT, GDTYP3D, GRDNM,    &
                       P_ALP3D, P_BET3D, P_GAM3D,    &
                       XCENT3D, YCENT3D, XVALS, YVALS )

        !......  Calculate cell size for each cell and conversion factor
        ALLOCATE( PRCNT2KM2( NCOLS, NROWS ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRCNT2KM2', PNAME )
        DO J = 1, NROWS
            DO I = 1, NCOLS
                XCELL = ABS( XVALS( I + 1, J ) - XVALS( I, J ) )
                YCELL = ABS( YVALS( I, J + 1 ) - YVALS( I, J ) )
                PRCNT2KM2( I,J ) = XCELL * YCELL * 1E-08
            END DO
        END DO
    ELSE
        ALLOCATE( PRCNT2KM2( 1, 1 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PRCNT2KM2', PNAME )
        PRCNT2KM2( 1, 1 ) = XCELL3D * YCELL3D * 1E-08
    END IF

    !......  Store landuse variable names from first file
    NVEG = NVARS4
    ALLOCATE( LUNAM( NVEG ),    &
             VGLIST( NVEG ),    &
             VEGID ( NVEG ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LUNAM...VEGID', PNAME )

    DO I = 1, NVEG
        LUNAM ( I ) = TRIM(VNAME3D ( I ) )
        VGLIST( I ) = TRIM(VNAME3D ( I ) )
        VEGID ( I ) = VGLIST ( I )
    END DO

    !......  Set up header variables for output file B3GRD
    NROWS3D = NROWS
    NCOLS3D = NCOLS
    GDNAM3D = GRDNM
    FTYPE3D = GRDDED3

    SDATE3D = 0           !  n/a
    STIME3D = 0           !  n/a
    TSTEP3D = 0           !  time independent
    NVARS3D = ( (NSEF-1) + NLAI ) * NSEASONS + NNO        ! special treatment of NO
    NLAYS3D = 1
    NTHIK3D = 1
    VGTYP3D = IMISS3
    VGTOP3D = AMISS3

    FDESC3D = ' '       ! array

    FDESC3D( 1 ) = 'BEIS3 normalized emissions values.'
    FDESC3D( 2 ) = '/FROM/ '    // PNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
    FDESC3D( 4 ) = '/LANDUSE/ BELD5 '                     ! BELD5 only works with BEISv3.7 or later

    I = 0

    !......  Set up variable names and output units
    DO M = 1, NSEASONS
        DO B = 1, NSEF

            IF( BIOTYPES( B ) /= 'NO' ) THEN
                BTMP = BIOTYPES( B )
                I    = I + 1
                VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // SEASON( M )
                VDESC3D( I ) = 'normalized emissions'
                UNITS3D( I ) = 'gramsC/hour'
                VTYPE3D( I ) = M3REAL
            END IF
        END DO

        DO N = 1, NLAI

            BTMP = LAITYPES( N )
            I    = I + 1

            VNAME3D( I ) = 'LAI_' // TRIM( BTMP ) // SEASON( M )
            VDESC3D( I ) = 'normalized emissions'
            UNITS3D( I ) = 'index'
            VTYPE3D( I ) = M3REAL
        END DO
    END DO

    !......  Handle NO types (not dependent on season)
    DO B = 1, NSEF

        BTMP = BIOTYPES( B )

        IF( BIOTYPES( B ) == 'NO' ) THEN

            I    = I + 1
            BTMP = BIOTYPES( B )
            VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_GROW'
            VDESC3D( I ) = 'normalized emissions for NO AG_GROW'
            UNITS3D( I ) = 'gramsN/hour'
            VTYPE3D( I ) = M3REAL

            I = I + 1
            VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'AG_NONGROW'
            VDESC3D( I ) = 'normalized emissions for NO AG_NONGROW'
            UNITS3D( I ) = 'gramsN/hour'
            VTYPE3D( I ) = M3REAL

            I = I + 1
            VNAME3D( I ) = 'AVG_' // TRIM( BTMP ) // 'NONAG'
            VDESC3D( I ) = 'normalized emissions for NO NONAG'
            UNITS3D( I ) = 'gramsN/hour'
            VTYPE3D( I ) = M3REAL

        END IF
    END DO

    !......  Open output file
    ENAME = PROMPTMFILE(    &
          'Enter logical name for NORMALIZED emissions output file',    &
          FSUNKN3, 'B3GRD', PNAME )

    !......  Get length of BFAC file
    NLINES = GETFLINE( FDEV, 'Emissions factor file' )

    !......  Allocate memory for emission factor variables
    ALLOCATE(  EMFAC( NVEG, NSEF ), &
                 LAI( NVEG ),       &
                 SLW( NVEG ),       &
                WFAC( NVEG ),       &
               LFBIO( NVEG ),       &
              LUINDX( NVARS4 ),     &
                LUSE( NCOLS, NROWS, NVEG ),    &
              AVGLAI( NCOLS, NROWS, NLAI, NSEASONS ),    &
             AVGEMIS( NCOLS, NROWS, NSEF, NSEASONS ),    &
              NOEMIS( NCOLS, NROWS, NNO ),    &
               SUMEM( NSEF ),       &
                NOEM( NNO ),        &
              SUMEMW( NSEF ),       &
              SUMLAI( NLAI ),       &
             SUMLAIW( NLAI ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMFAC...SUMLAIW', PNAME )

    EMFAC  = BADVAL3
    LAI    = IMISS3
    SLW    = BADVAL3
    WFAC   = BADVAL3
    LFBIO  = BADVAL3

    !......  Read emissions factor file
    MESG = 'Reading emissions factor file...'
    CALL M3MSG2( MESG )

    WRITE( MESG,94010 )    &
        'Number of landuse types in factor file: ', NVEG
    CALL M3MSG2( MESG )

    !......  This routine reads in emission factors
    CALL RDB4FAC_CSV( PNAME, NLINES, NSEF, FDEV, NVARS4,  VGLIST,    &
                  BIOTYPES, LAI, LFBIO, WFAC, SLW, EMFAC)

    DO I = 1, NVEG
        IF (LAI(I) .EQ. IMISS3 ) THEN
            MESG = 'ERROR: MISSING LAI FOR VEG TYPE: '//VGLIST(I)
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        ENDIF

        IF (SLW(I) .LT. AMISS3 ) THEN
            MESG = 'ERROR: MISSING SLW FOR VEG TYPE: '//VGLIST(I)
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        ENDIF

        IF (WFAC(I) .LT. AMISS3 ) THEN
            MESG = 'ERROR: MISSING WINTER FAC FOR VEG TYPE: '//VGLIST(I)
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        ENDIF

        IF (LFBIO(I) .LT. AMISS3 ) THEN
            MESG = 'ERROR: MISSING LFBIO FOR VEG TYPE: '//VGLIST(I)
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        ENDIF
    ENDDO

    DO J = 1, NSEF
    DO I = 1, NVEG
        IF (EMFAC(I,J) .LT. AMISS3 ) THEN
            MESG ='ERROR: MISSING EMISSION FACTOR FOR VEG TYPE:'    &
                   // VGLIST(I)// 'AND SPECIES: '//BIOTYPES(J)
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        ENDIF
    ENDDO
    ENDDO

    IF( EFLAG ) THEN
        MESG = 'ERRORs are occurred above. Check the messages'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    LUINDX = -9       ! array

    !......  Check to see if there are emissions factors for all landuse types
    DO I = 1, NVARS4

        IFOUND = 0

        LUNAM_TMP = LUNAM( I )

        DO J = 1, NVEG

            IF (VEGID ( J ) .EQ. LUNAM_TMP ) THEN
                IFOUND = 1
                LUINDX( I ) = J
            END IF

            !C......  Store vegids for certain USGS categories for use later
            !        IF ( VEGID( J ) .EQ. 'USGS_decidforest' ) IUSDEC = J
            !        IF ( VEGID( J ) .EQ. 'USGS_evbrdleaf  ' ) IUSBRD = J
            !        IF ( VEGID( J ) .EQ. 'USGS_coniferfor ' ) IUSCON = J
            !        IF ( VEGID( J ) .EQ. 'USGS_mxforest   ' ) IUSMIX = J
            !        IF ( VEGID( J ) .EQ. 'USGS_shrubland  ' ) IUSSHR = J

            IF ( VEGID( J ) .EQ. 'NLCD_81         ' ) NLCD81 = J
            IF ( VEGID( J ) .EQ. 'NLCD_82         ' ) NLCD82 = J
            IF ( VEGID( J ) .EQ. 'MODIS_12        ' ) MODIS12 = J
            ! 3/2022 was MODIS_47 before corrected here
            IF ( VEGID( J ) .EQ. 'MODIS_14        ' ) MODIS14 = J

            IF ( VEGID( J ) .EQ. 'Alfalfa         ' ) IALFAL   = J
            IF ( VEGID( J ) .EQ. 'Alfalfa_ir      ' ) IALFAL_IR = J
            IF ( VEGID( J ) .EQ. 'Barley          ' ) IBARLE   = J
            IF ( VEGID( J ) .EQ. 'Barley_ir       ' ) IBARLE_IR = J
            IF ( VEGID( J ) .EQ. 'CornGrain       ' ) ICORNG     = J
            IF ( VEGID( J ) .EQ. 'CornGrain_ir    ' ) ICORNG_IR  = J
            IF ( VEGID( J ) .EQ. 'CornSilage      ' ) ICORNS     = J
            IF ( VEGID( J ) .EQ. 'CornSilage_ir   ' ) ICORNS_IR  = J
            IF ( VEGID( J ) .EQ. 'Cotton          ' ) ICOTTO    = J
            IF ( VEGID( J ) .EQ. 'Cotton_ir       ' ) ICOTTO_IR = J
            IF ( VEGID( J ) .EQ. 'Other_Grass     ' ) IGRASS = J
            IF ( VEGID( J ) .EQ. 'Other_Grass_ir  ' ) IGRASS_IR = J
            IF ( VEGID( J ) .EQ. 'Hay             ' ) IHAY      = J
            IF ( VEGID( J ) .EQ. 'Hay_ir          ' ) IHAY_IR   = J
            ! 3/2022 Other_crop should have been Other_Crop was Other_crop
            IF ( VEGID( J ) .EQ. 'Other_Crop      ' ) IMISCC = J
            IF ( VEGID( J ) .EQ. 'Other_Crop_ir   ' ) IMISCC_IR = J
            IF ( VEGID( J ) .EQ. 'Oats            ' ) IOATS  = J
            IF ( VEGID( J ) .EQ. 'Oats_ir         ' ) IOATS_IR  = J
            IF ( VEGID( J ) .EQ. 'Peanuts         ' ) IPEANU    = J
            IF ( VEGID( J ) .EQ. 'Peanuts_ir      ' ) IPEANU_IR = J
            ! 3/2022 Potatoes misspelled in earlier versions
            IF ( VEGID( J ) .EQ. 'Potatoes       ' ) IPOTAT    = J
            IF ( VEGID( J ) .EQ. 'Potatoes_ir    ' ) IPOTAT_IR = J
            IF ( VEGID( J ) .EQ. 'Rice            ' ) IRICE     = J
            IF ( VEGID( J ) .EQ. 'Rice_ir         ' ) IRICE_IR  = J
            IF ( VEGID( J ) .EQ. 'Rye             ' ) IRYE      = J
            IF ( VEGID( J ) .EQ. 'Rye_ir          ' ) IRYE_IR   = J
            IF ( VEGID( J ) .EQ. 'SorghumGrain    ' ) ISORGHG    = J
            IF ( VEGID( J ) .EQ. 'SorghumGrain_ir ' ) ISORGHG_IR = J
            IF ( VEGID( J ) .EQ. 'SorghumSilage   ' ) ISORGHS    = J
            IF ( VEGID( J ) .EQ. 'SorghumSilage_ir' ) ISORGHS_IR = J
            IF ( VEGID( J ) .EQ. 'Soybeans        ' ) ISOYBE = J
            IF ( VEGID( J ) .EQ. 'Soybeans_ir     ' ) ISOYBE_IR = J
            IF ( VEGID( J ) .EQ. 'Beans           ' ) IBEANS = J
            IF ( VEGID( J ) .EQ. 'Beans_ir        ' ) IBEANS_IR = J
            IF ( VEGID( J ) .EQ. 'BeansEdible     ' ) IBEANSED = J
            IF ( VEGID( J ) .EQ. 'BeansEdible_ir  ' ) IBEANSED_IR =J
            IF ( VEGID( J ) .EQ. 'Canola          ' ) ICANOLA = J
            IF ( VEGID( J ) .EQ. 'Canola_ir       ' ) ICANOLA_IR = J
            IF ( VEGID( J ) .EQ. 'Wheat_Spring    ' ) IWHEATS    = J
            IF ( VEGID( J ) .EQ. 'Wheat_Spring_ir ' ) IWHEATS_IR = J
            IF ( VEGID( J ) .EQ. 'Wheat_Winter    ' ) IWHEATW    = J
            IF ( VEGID( J ) .EQ. 'Wheat_Winter_ir ' ) IWHEATW_IR = J

        END DO
        IF( IFOUND .EQ. 0 ) THEN
            MESG =   LUNAM_TMP // ' does NOT have emissions factors in BFAC file'
            CALL M3WARN( PNAME, 0, 0, MESG )
        END IF

    END DO

    !......  Read the gridded landuse from the landuse files
    DO M = 1, NVARS4

        N = LUINDX( M )

        IF( N > 0 ) THEN

            CALL M3MESG( 'Reading ' // LUNAM( M ) )

            IF( .NOT. READ3( GNAMET, LUNAM( M ), 1,0,0, LUSE(1,1,N) ) ) THEN
                MESG = 'Could not find "' // LUNAM( M ) //    &
                       '" in file "' // TRIM( GNAMET ) // '"'
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF
        END IF

    END DO

    AVGEMIS = 0.0      !  array
    AVGLAI  = 0.0      !  array
    NOEMIS  = 0.0      !  array
    LAI_SAVE_INDEX(1:3) = 0
    !......  Calculate normalized fluxes

    DO J = 1, NROWS
    DO I = 1, NCOLS

    !......  Initialize variables
        SUMEM   = 0.0         ! array
        SUMEMW  = 0.0         ! array
        NOEM    = 0.0         ! array
        SUMLAI  = 0.0         ! array
        SUMLAIW = 0.0         ! array

        DO M = 1, NVEG

            !......  Assuming that land use is in percentages
            !        Converting to units of area (km**2)
            IF ( VFLAG ) THEN
                VEGAREA = LUSE( I, J, M ) * PRCNT2KM2( I, J )
            ELSE
                VEGAREA = LUSE( I, J, M ) * PRCNT2KM2( 1, 1 )
            END IF

            DO N = 1, NSEF

                BTMP = BIOTYPES( N )

                !......  Special handling for NO emissions
                IF ( BTMP == 'NO' ) THEN

                    IF ( VEGAREA > 0. ) THEN

                        IF ( IS_BG( M, MODIS12 , MODIS14, NLCD81,   &
                               NLCD82,IALFAL  , IBARLE,  ICORNG,    &
                               ICORNS, ICOTTO  , IGRASS,  IHAY,     &
                               IMISCC,                              &
                               IOATS   , IPEANU, IPOTAT , IRICE,    &
                               IRYE    ,ISORGHG, ISORGHS, ISOYBE,   &
                               IBEANSED, IBEANS, ICANOLA,           &
                               IWHEATW, IWHEATS,                    &
                               IALFAL_IR  , IBARLE_IR,              &
                               ICORNG_IR, ICORNS_IR,                &
                               ICOTTO_IR  , IGRASS_IR,              &
                               IHAY_IR  , IMISCC_IR,                &
                               IOATS_IR   , IPEANU_IR,              &
                               IPOTAT_IR , IRICE_IR ,               &
                               IRYE_IR    ,ISORGHG_IR,              &
                               ISORGHS_IR, ISOYBE_IR,               &
                               IBEANSED_IR, IBEANS_IR,              &
                               ICANOLA_IR, IWHEATW_IR,              &
                               IWHEATS_IR ) ) THEN


                            !...... Mosaics of small-scale cultivation 40-60% with
                            !......natural tree, shrub, or herbaceous vegetation.
                            !......(from latest MODIS doc for MODIS_14 at
                            !      https://lpdaac.usgs.gov/documents/101/MCD12_User_Guide_V6.pdf
                            !...... using 50% below

                            !......  Compute NO emissions for agriculture regions
                            !        during growing season
                            IF( IS_MAG (M,MODIS14) ) THEN
                                NOEM( 1 ) = NOEM( 1 ) + VEGAREA * EMFAC(M,N)     *0.500
                                NOEM( 2 ) = NOEM( 2 ) + VEGAREA * EMFAC(IGRASS,N)*0.500
                                NOEM( 3 ) = NOEM( 3 ) + VEGAREA * EMFAC(M,N)     *0.500
                            ELSE
                                NOEM( 1 ) = NOEM( 1 ) + VEGAREA * EMFAC(M,N)
                                !......  Compute NO emissions for agriculture regions
                                !        outside growing season
                                NOEM( 2 ) = NOEM( 2 ) + VEGAREA * EMFAC(IGRASS,N)
                            ENDIF

                        ELSE

                            !......  Compute NO emissions for Non-Agriculture regions
                            NOEM( 3 ) = NOEM( 3 ) + VEGAREA * EMFAC(M,N)

                        END IF

                    END IF

                ELSE

                    EFTMP = EMFAC( M, N )

                    IF ( VEGAREA > 0. ) THEN

                        !......  Compute summer emissions
                        SUMEM( N ) = SUMEM( N ) + VEGAREA * EFTMP

                        !......  Compute winter emissions
                        SUMEMW( N ) = SUMEMW( N ) + VEGAREA * EFTMP * WFAC( M )
                    END IF

                        !......  Compute LAI on ISOP and MBO and METH

                    IF ( BTMP == 'ISOP') THEN
                        LAI_SAVE_INDEX(1) = N
                        SUMLAI ( 1 ) = SUMLAI ( 1 ) + VEGAREA * LAI( M ) * EFTMP
                        SUMLAIW( 1 ) = SUMLAIW( 1 ) + VEGAREA * LAI( M ) * EFTMP * WFAC( M )
                    ELSE IF ( BTMP == 'MBO' ) THEN
                        LAI_SAVE_INDEX(2) = N
                        SUMLAI ( 2 ) = SUMLAI ( 2 ) + VEGAREA * LAI( M ) * EFTMP
                        SUMLAIW( 2 ) = SUMLAIW( 2 ) + VEGAREA * LAI( M ) * EFTMP * WFAC( M )
                    ELSE IF ( TRIM( BTMP ) == 'METH' ) THEN
                        LAI_SAVE_INDEX(3) = N
                        SUMLAI ( 3 ) = SUMLAI ( 3 ) + VEGAREA * LAI( M ) * EFTMP
                        SUMLAIW( 3 ) = SUMLAIW( 3 ) + VEGAREA * LAI( M ) * EFTMP * WFAC( M )

                    END IF

                END IF          ! check if NO emissions

            END DO          ! end of emis fac loop
        END DO          ! end of veg land use loop2


        DO K = 1, NLAI

            IF ( SUMLAI( K ) <= 1E-06 ) THEN
                AVGLAI( I, J, K, 1 ) = 0.
            ELSE IF ( SUMEM(  LAI_SAVE_INDEX(K) ) == 0. ) THEN
                AVGLAI(  I, J, K, 1 ) = 0.
            ELSE
                AVGLAI( I, J, K, 1 ) =  SUMLAI( K )/SUMEM( LAI_SAVE_INDEX(K) )
            ENDIF

            IF ( SUMLAIW( K ) <= 1E-06 ) THEN
                AVGLAI( I, J, K, 2 ) = 0.
            ELSE IF ( SUMEMW(  LAI_SAVE_INDEX(K) ) == 0. ) THEN
                AVGLAI(  I, J, K, 2 ) = 0.
            ELSE
                AVGLAI( I, J, K, 2) =  SUMLAIW( K ) /SUMEMW( LAI_SAVE_INDEX(K) )
            ENDIF

        END DO

        DO N = 1, NSEF

            BTMP = BIOTYPES( N )

            !......  Check for NO emissions
            IF ( BTMP == 'NO' ) THEN
                NOEMIS( I, J, 1 ) = NOEM( 1 )
                NOEMIS( I, J, 2 ) = NOEM( 2 )
                NOEMIS( I, J, 3 ) = NOEM( 3 )
            ELSE
                AVGEMIS( I, J, N, 1 ) = SUMEM( N )
                AVGEMIS( I, J, N, 2 ) = SUMEMW( N )
            END IF

        END DO          ! end loop over emission factors

    END DO          ! end loop over rows
    END DO      ! end loop over columns

    !......  Write output file
    I = 0
    DO M = 1, NSEASONS
        DO B = 1, NSEF
            BTMP = BIOTYPES( B )

            !......  Handle types other than NO
            IF ( BTMP /= 'NO' ) THEN

                I = I + 1
                IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0,0, AVGEMIS(1,1,B,M) ) ) THEN
                    MESG = 'Could not write "' // TRIM( VNAME3D( I ) ) //    &
                            '" to "' // TRIM( ENAME ) // '"'
                    CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
                END IF
            END IF

        END DO

        DO N = 1, NLAI

            I = I + 1
            IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0,0, AVGLAI(1,1,N,M) ) ) THEN
                MESG = 'Could not write "' // TRIM( VNAME3D( I ) ) //    &
                        '" to "' // TRIM( ENAME ) // '"'
                CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
            END IF

        END DO

    END DO

    I = I + 1
    IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0, NOEMIS(1,1,1) ) ) THEN
        MESG = 'Could not write "' // TRIM( VNAME3D( I ) ) //    &
                '" to "' // TRIM( ENAME ) // '"'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    I = I + 1
    IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0, NOEMIS(1,1,2) ) ) THEN
        MESG = 'Could not write "' // TRIM( VNAME3D( I ) ) //    &
                '" to "' // TRIM( ENAME ) // '"'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    I = I + 1
    IF ( .NOT. WRITE3( ENAME, VNAME3D( I ), 0, 0, NOEMIS(1,1,3) ) ) THEN
        MESG = 'Could not write "' // TRIM( VNAME3D( I ) ) //    &
                '" to "' // TRIM( ENAME ) // '"'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    !......  End of subroutine
    RETURN

    !******************  FORMAT  STATEMENTS   ******************************
    !......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I5, :, 2X ) )


CONTAINS


    !-----------------------------------------------------------------------------
    !......  This internal function checks for "one-third" agricultural areas

    LOGICAL FUNCTION IS_MAG( M,MODIS14 )

        IMPLICIT NONE

    !......  Function arguments
        INTEGER, INTENT(IN) :: M
        INTEGER, INTENT(IN) :: MODIS14
        
        IS_MAG = ( M == MODIS14 )
        RETURN

    END FUNCTION IS_MAG

    !-----------------------------------------------------------------------------
    !......  This internal function checks for agricultural areas

    LOGICAL FUNCTION IS_BG( M, MODIS12 , MODIS14, NLCD81,   &
        NLCD82, IALFAL  , IBARLE,  ICORNG, ICORNS,          &
        ICOTTO  , IGRASS,  IHAY  , IMISCC,                  &
        IOATS   , IPEANU, IPOTAT , IRICE ,                  &
        IRYE    ,ISORGHG, ISORGHS, ISOYBE,                  &
        IBEANSED, IBEANS, ICANOLA, IWHEATW, IWHEATS,        &
        IALFAL_IR  , IBARLE_IR,  ICORNG_IR, ICORNS_IR,      &
        ICOTTO_IR  , IGRASS_IR,  IHAY_IR  , IMISCC_IR,      &
        IOATS_IR   , IPEANU_IR, IPOTAT_IR , IRICE_IR ,      &
        IRYE_IR    ,ISORGHG_IR, ISORGHS_IR, ISOYBE_IR,      &
        IBEANSED_IR, IBEANS_IR, ICANOLA_IR, IWHEATW_IR,     &
        IWHEATS_IR)

        IMPLICIT NONE

        !......  Function arguments
        INTEGER, INTENT(IN) :: M
        INTEGER, INTENT(IN) :: MODIS12
        INTEGER, INTENT(IN) :: MODIS14
        INTEGER, INTENT(IN) :: NLCD81
        INTEGER, INTENT(IN) :: NLCD82

        INTEGER, INTENT(IN) ::         IALFAL          ! Alfalfa
        INTEGER, INTENT(IN) ::         IBARLE          ! Barley
        INTEGER, INTENT(IN) ::         ICORNG           ! Corn
        INTEGER, INTENT(IN) ::         ICORNS           ! Corn
        INTEGER, INTENT(IN) ::         ICOTTO          ! Cotton
        INTEGER, INTENT(IN) ::         IGRASS          ! Grass
        INTEGER, INTENT(IN) ::         IHAY            ! Hay
        INTEGER, INTENT(IN) ::         IMISCC          ! Misc_crop
        INTEGER, INTENT(IN) ::         IOATS           ! Oats
        INTEGER, INTENT(IN) ::         IPEANU          ! Peanuts
        INTEGER, INTENT(IN) ::         IPOTAT          ! Potatoes
        INTEGER, INTENT(IN) ::         IRICE           ! Rice
        INTEGER, INTENT(IN) ::         IRYE            ! Rye
        INTEGER, INTENT(IN) ::         ISORGHG          ! Sorghum
        INTEGER, INTENT(IN) ::         ISORGHS          ! Sorghum
        INTEGER, INTENT(IN) ::         ISOYBE          ! Soybeans
        INTEGER, INTENT(IN) ::         IBEANSED         ! Edible Benas
        INTEGER, INTENT(IN) ::         IBEANS          ! Beans
        INTEGER, INTENT(IN) ::         ICANOLA         ! Canola
        INTEGER, INTENT(IN) ::         IWHEATW          ! Wheat
        INTEGER, INTENT(IN) ::         IWHEATS          ! Wheat


        INTEGER, INTENT(IN) ::         IALFAL_IR          ! Alfalfa
        INTEGER, INTENT(IN) ::         IBARLE_IR          ! Barley
        INTEGER, INTENT(IN) ::         ICORNG_IR           ! Corn
        INTEGER, INTENT(IN) ::         ICORNS_IR           ! Corn
        INTEGER, INTENT(IN) ::         ICOTTO_IR          ! Cotton
        INTEGER, INTENT(IN) ::         IGRASS_IR          ! Grass
        INTEGER, INTENT(IN) ::         IHAY_IR            ! Hay
        INTEGER, INTENT(IN) ::         IMISCC_IR          ! Misc_crop
        INTEGER, INTENT(IN) ::         IOATS_IR           ! Oats
        INTEGER, INTENT(IN) ::         IPEANU_IR          ! Peanuts
        INTEGER, INTENT(IN) ::         IPOTAT_IR          ! Potatoes
        INTEGER, INTENT(IN) ::         IRICE_IR           ! Rice
        INTEGER, INTENT(IN) ::         IRYE_IR            ! Rye
        INTEGER, INTENT(IN) ::         ISORGHG_IR          ! Sorghum
        INTEGER, INTENT(IN) ::         ISORGHS_IR          ! Sorghum
        INTEGER, INTENT(IN) ::         ISOYBE_IR          ! Soybeans
        INTEGER, INTENT(IN) ::         IBEANSED_IR         ! Edible Benas
        INTEGER, INTENT(IN) ::         IBEANS_IR          ! Beans
        INTEGER, INTENT(IN) ::         ICANOLA_IR         ! Canola
        INTEGER, INTENT(IN) ::         IWHEATS_IR          ! Wheat Spring Irrigated
        INTEGER, INTENT(IN) ::         IWHEATW_IR          ! Wheat Winter Irrigated
        !-----------------------------------------------------------------------------

        IS_BG = .FALSE.

        IF( M == MODIS12 ) IS_BG = .TRUE.
        IF( M == MODIS14 ) IS_BG = .TRUE.
        IF( M == NLCD81  ) IS_BG = .TRUE.
        IF( M == NLCD82  ) IS_BG = .TRUE.
        IF( M == IALFAL  ) IS_BG = .TRUE.
        IF( M == IBARLE  ) IS_BG = .TRUE.
        IF( M == ICORNG  ) IS_BG = .TRUE.
        IF( M == ICORNS  ) IS_BG = .TRUE.
        IF( M == ICOTTO ) IS_BG = .TRUE.
        IF( M == IGRASS ) IS_BG = .TRUE.
        IF( M == IHAY   ) IS_BG = .TRUE.
        IF( M == IMISCC ) IS_BG = .TRUE.
        IF( M == IOATS  ) IS_BG = .TRUE.
        IF( M == IPEANU ) IS_BG = .TRUE.
        IF( M == IPOTAT ) IS_BG = .TRUE.
        IF( M == IRICE  ) IS_BG = .TRUE.
        IF( M == IRYE   ) IS_BG = .TRUE.
        IF( M == ISORGHG ) IS_BG = .TRUE.
        IF( M == ISORGHS ) IS_BG = .TRUE.
        IF( M == ISOYBE ) IS_BG = .TRUE.
        IF( M == IBEANS ) IS_BG = .TRUE.
        IF( M == IBEANSED ) IS_BG = .TRUE.
        IF( M == ICANOLA ) IS_BG = .TRUE.
        IF( M == IWHEATS ) IS_BG = .TRUE.
        IF( M == IWHEATW ) IS_BG = .TRUE.
        IF( M == IALFAL_IR  ) IS_BG = .TRUE.
        IF( M == IBARLE_IR  ) IS_BG = .TRUE.
        IF( M == ICORNG_IR  ) IS_BG = .TRUE.
        IF( M == ICORNS_IR  ) IS_BG = .TRUE.
        IF( M == ICOTTO_IR ) IS_BG = .TRUE.
        IF( M == IGRASS_IR ) IS_BG = .TRUE.
        IF( M == IHAY_IR   ) IS_BG = .TRUE.
        IF( M == IMISCC_IR ) IS_BG = .TRUE.
        IF( M == IOATS_IR  ) IS_BG = .TRUE.
        IF( M == IPEANU_IR ) IS_BG = .TRUE.
        IF( M == IPOTAT_IR ) IS_BG = .TRUE.
        IF( M == IRICE_IR  ) IS_BG = .TRUE.
        IF( M == IRYE_IR   ) IS_BG = .TRUE.
        IF( M == ISORGHG_IR ) IS_BG = .TRUE.
        IF( M == ISORGHS_IR ) IS_BG = .TRUE.
        IF( M == ISOYBE_IR ) IS_BG = .TRUE.
        IF( M == IBEANS_IR ) IS_BG = .TRUE.
        IF( M == IBEANSED_IR ) IS_BG = .TRUE.
        IF( M == ICANOLA_IR ) IS_BG = .TRUE.
        IF( M == IWHEATS_IR ) IS_BG = .TRUE.
        IF( M == IWHEATW_IR ) IS_BG = .TRUE.
        RETURN

    END FUNCTION IS_BG


END SUBROUTINE NORMBEIS370
