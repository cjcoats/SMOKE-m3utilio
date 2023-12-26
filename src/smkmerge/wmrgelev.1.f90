
SUBROUTINE WMRGELEV( VNAME, NSRC, NMAJOR, JDATE, JTIME )

!***********************************************************************
!  subroutine WMRGELEV body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to write out an ASCII elevated point
!      source emissions file. The default format is the UAM/CAMx elevated
!      sources file.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 4/2000 by M. Houyoux
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
!****************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: EXPLFLAG, EMLAYS, PVSDATE, PVSTIME, PVNAME,&
    &                    PENAME, PSDEV, NPSRC, NMSPC, EMNAM,&
    &                    SDATE, STIME, EDATE, ETIME, TSTEP,&
    &                    JSTACK, EVDEV, LFRAC

!...........   This module is the source inventory arrays
    USE MODSOURC, ONLY: CSOURC

!.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: NHRSRC, GRPCOL, GRPROW, GRPXX, GRPYY,&
    &                   GRPHT, GRPDM, GRPTK, GRPVE,&
    &                   NGROUP, GRPGID, ELEVSRC, LPING,&
    &                    ELEVSIDX, GROUPID, INDXH, ELEVEMIS

!.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NCOLS, NROWS, XCELL, YCELL, XORIG, YORIG,&
    &                   GRDNM, P_ALP, GDTYP, VGTYP, VGLVS

    IMPLICIT NONE

!.........  INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: VNAME   ! variable name to output
    INTEGER     , INTENT (IN) :: NSRC    ! no. sources
    INTEGER     , INTENT (IN) :: NMAJOR  ! no. elevated sources
    INTEGER     , INTENT (IN) :: JDATE   ! Julian date to output (YYYYDDD)
    INTEGER     , INTENT (IN) :: JTIME   ! time to output (HHMMSS)

!.........  Local parameters
    INTEGER, PARAMETER :: STKLEN3 = FIPLEN3 + PLTLEN3 + CHRLEN3

!.........  Variables allocated by module settings...
    INTEGER, ALLOCATABLE, SAVE :: ESRTIDX( : ) ! major srcs sorting index
    INTEGER, ALLOCATABLE, SAVE :: EIDXSRC( : ) ! Source ID per major source
    INTEGER, ALLOCATABLE, SAVE :: EIDX2  ( : ) ! another index
    INTEGER, ALLOCATABLE, SAVE :: ELAYER ( : ) ! addt'l srcs for explicit

    REAL, ALLOCATABLE :: EOUTHT ( : ) ! output stack heights [m]
    REAL, ALLOCATABLE :: EOUTDM ( : ) ! output stack diameters [m]
    REAL, ALLOCATABLE :: EOUTTK ( : ) ! output exit tmprs [K]
    REAL, ALLOCATABLE :: EOUTVE ( : ) ! output exit velocities [m/hr]
    REAL, ALLOCATABLE :: EOUTXL ( : ) ! output x-loc [units of grid]
    REAL, ALLOCATABLE :: EOUTYL ( : ) ! output y-loc [units of grid]

    CHARACTER(STKLEN3), ALLOCATABLE, SAVE :: ECSRCA  (:)! FIPS//plt//stk
    CHARACTER(STKLEN3), ALLOCATABLE       :: EOUTCSRC(:)! FIPS//plt//stk

!.........  Allocatable array for fake stack heights for explicit plums
    REAL, ALLOCATABLE :: LAYRMID( : )

!.........  Fixed size arrays
    REAL            ESUM   ( MXLAYS3 )  ! emissions data sum

!.........  UAM-format specific variables
    INTEGER         DDEVOUT             ! diffbreak file unit number
    INTEGER         ESWITCH             ! 0=no, 1=yes - print vert mthds tbl
    INTEGER         GSWITCH             ! 0=no, 1=yes - print output grid
    INTEGER         LSWITCH             ! 0=no, 1=yes - print src locs tbl
    INTEGER      :: MDEVOUT = 0         ! metscalars file unit number
    INTEGER         MSWITCH             ! 0=no, 1=yes - print methods tbl
    INTEGER         NH                  ! tmp no. explicit sources
    INTEGER         NPARAM              ! no. paramaters for control pkt
    INTEGER         NULAYS              ! no. UAM model layers
    INTEGER         NZLOWR              ! no. layers below Diffbreak height
    INTEGER         NZUPPR              ! no. layers above Diffbreak height
    INTEGER      :: PDEVOUT = 0         ! PTSRCE output unit no.
    INTEGER      :: RDEVOUT = 0         ! regiontop file unit number
    INTEGER      :: TDEVOUT = 0         ! temperatur file unit number
    INTEGER         USWITCH             ! 0=no, 1=yes - print units table
    INTEGER         VSWITCH             ! 0=no, 1=yes - print values tbl
    INTEGER      :: WDEVOUT = 0         ! wind file unit number

    REAL            HTSUR               ! height of surface layer [m]
    REAL            HTLOWR              ! min cell ht b/w sfc and diffbr [m]
    REAL            HTUPPR              ! min cell ht b/w diffbr and top [m]

    CHARACTER(10)       :: SPCNAM       ! UAM-format species name
    CHARACTER(10), SAVE :: VTYPE        ! User-spec vertical method type
    CHARACTER(44)          NOTEDEF      ! Default note
    CHARACTER(44)          UNOTE        ! UAM file note from env variable
    CHARACTER(60)       :: FNOTE        ! UAM file header note

!.........  Other local variables
    INTEGER          I, J, K, KK, L, LL, LM, LN, M, N, S

    INTEGER          COL                 ! tmp column
    INTEGER          EML                 ! tmp emissions layers
    INTEGER          FIP                 ! tmp co/st/cy code
    INTEGER          GID                 ! tmp stack group ID
    INTEGER          IBD                 ! beg 5-digit Julian date of step
    INTEGER          IBT                 ! beginning time of time step
    INTEGER          IED                 ! end 5-digit Julian date of step
    INTEGER          IET                 ! ending time of time step
    INTEGER          IOS                 ! i/o status
    INTEGER          JDATEP1             ! julian date, plus 1 time step
    INTEGER          JTIMEP1             ! time, plus 1 time step
    INTEGER          NEXPLOOP            ! no. for explicit plume loop
    INTEGER, SAVE :: NOUT                ! number output stacks
    INTEGER, SAVE :: PTIME = -9          ! previous call's time (HHMMSS)
    INTEGER          ROW                 ! tmp row

    REAL             DM                  ! tmp stack diameter [m]
    REAL             HT                  ! tmp stack height [m]
    REAL             XLOC                ! tmp x location
    REAL             XMAX                ! rightmost x location
    REAL             YLOC                ! tmp y location
    REAL             YMAX                ! topmost y location

    LOGICAL       :: EFLAG    = .FALSE.  ! true: error occurred
    LOGICAL       :: FIRSTIME = .TRUE.   ! true: first time routine called

    CHARACTER(16),SAVE :: OUTFMT       ! output format for elevated ASCII
    CHARACTER(100)     EFMT         ! output emissions foamat
    CHARACTER(200)     BUFFER       ! source chars buffer
    CHARACTER(300)     MESG         ! message buffer

    CHARACTER(FIPLEN3) CFIP         ! tmp country/state/county code
    CHARACTER(PLTLEN3) FCID         ! tmp facility ID
    CHARACTER(CHRLEN3) SKID         ! tmp stack ID
    CHARACTER(STKLEN3) ECS          ! stack elevated source chars
    CHARACTER(STKLEN3) PECS         ! tmp previous ECS
    CHARACTER(NAMLEN3) GRDENV       ! gridded output units from envrmt

    CHARACTER(16), PARAMETER :: PROGNAME = 'WMRGELEV' ! program name

!***********************************************************************
!   begin body of subroutine WMRGELEV

!.........  For the first time the subroutine is called...
    IF( FIRSTIME ) THEN

        MESG = 'Setting up to output ASCII elevated file...'
        CALL M3MSG2( MESG )

!.............  Get environment variable settings

!.............  Get type of plume rise calculation from environment
        MESG = 'Elevated ASCII output file format'
        CALL ENVSTR( 'SMK_ASCIIELEV_FMT', MESG, 'UAM', OUTFMT, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_ASCIIELEV_FMT"', 2 )
        END IF
        IF( VTYPE .EQ. 'UAM4' ) VTYPE = 'UAM'

!.............  Set size for allocating output elevated arrays
!.............  Adjustment for EXPLFLAG so that one record for each layer of
!               each explicit source can be inserted.
        I = NMAJOR
        IF ( EXPLFLAG ) I = I - NHRSRC + EMLAYS * NHRSRC

!.............  Allocate memory for local elevated sources arrays
        ALLOCATE( ESRTIDX( I ),&
        &          EIDXSRC( I ),&
        &           ECSRCA( I ),&
        &         EOUTCSRC( I ),&
        &           EOUTHT( I ),&
        &           EOUTDM( I ),&
        &           EOUTTK( I ),&
        &           EOUTVE( I ),&
        &           EOUTXL( I ),&
        &           EOUTYL( I ),&
        &            EIDX2( I ),&
        &           ELAYER( I ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ESRTIDX...ELAYER', PROGNAME )

        EIDX2 = 0   ! array
        ELAYER = 0  ! array

!.............  If needed for explicit plume rise, allocate and set local array
!               with the mid-point of layers to use as fake stack heights.
        IF ( EXPLFLAG ) THEN

            ALLOCATE( LAYRMID( EMLAYS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'LAYRMID', PROGNAME )

!.................  Method of midpoint depends on vertical structure
            SELECT CASE ( VGTYP )
              CASE ( VGHVAL3 )

                DO L = 1, EMLAYS
                    LAYRMID( L ) = VGLVS( L-1 ) + 0.5 *&
                    &             ( VGLVS( L ) - VGLVS( L-1 ) )
                END DO

              CASE DEFAULT
                WRITE( MESG,94010 ) 'Do not know ' //&
                &       'how to set layer midpoints for layer '//&
                &       'structure type', VGTYP
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SELECT

        END IF

!.............  Read stack groups file variables. Read in group ID using x-loc
!               array because SAFE_READ3 is expecting a real array to be
!               passed to it.
        I = PVSDATE
        J = PVSTIME
        CALL SAFE_READ3_I( PVNAME, 'COL' , 1, I, J, GRPCOL )
        CALL SAFE_READ3_I( PVNAME, 'ROW' , 1, I, J, GRPROW )
        CALL SAFE_READ3( PVNAME, 'XLOCA' , 1, I, J, GRPXX )
        CALL SAFE_READ3( PVNAME, 'YLOCA' , 1, I, J, GRPYY )
        CALL SAFE_READ3( PVNAME, 'STKHT' , 1, I, J, GRPHT )
        CALL SAFE_READ3( PVNAME, 'STKDM' , 1, I, J, GRPDM )
        CALL SAFE_READ3( PVNAME, 'STKTK' , 1, I, J, GRPTK )
        CALL SAFE_READ3( PVNAME, 'STKVE' , 1, I, J, GRPVE )

!.............  Read source characteristics from the point source inventory
!               file
        CALL RDINVCHR( 'POINT', PENAME, PSDEV, NPSRC, 1, 'CSOURC' )

!.............  Copy source characteristics from sources to elevated list with
!               only country/state/county, plant, and stack.  For IDA
!               inventories, the stack ID comes after the point ID, so need to
!               use position of stack ID in source definition to properly
!               build the elevated sources arrays
        DO S = 1, NSRC

            I    = ELEVSIDX( S )
            IF ( I .EQ. 0 ) CYCLE
            CFIP = CSOURC( S )( PTBEGL3( 1 ):PTENDL3( 1 ) )
            FCID = CSOURC( S )( PTBEGL3( 2 ):PTENDL3( 2 ) )
            SKID = CSOURC( S )( PTBEGL3(JSTACK):PTENDL3(JSTACK) )
            ESRTIDX( I ) = I
            EIDXSRC( I ) = S
            ECSRCA ( I ) = CFIP // FCID // SKID

        END DO

!.............  Sort elevated sources
        CALL SORTIC( NMAJOR, ESRTIDX, ECSRCA )

        XMAX = XORIG + XCELL * NCOLS
        YMAX = YORIG + YCELL * NROWS

!.............  Create indices from major-sources list to output list
!.............  Eliminate sources that are not in the domain
!.............  Add sources that have explicit plume rise and require fake
!               stacks to ensure that their emissions get put in the correct
!               model layers.
        K  = 0
        KK = 0
        M  = 0
        PECS = ' '
        DO I = 1, NMAJOR

            J   = ESRTIDX ( I )
            S   = EIDXSRC ( J )
            ECS = ECSRCA  ( J )

!.................  Find group for this record in stack groups list
            GID = GROUPID( S )
            N = FIND1( GID, NGROUP, GRPGID )

!.................  If group not found, error...
            IF( N .LE. 0 ) THEN
                EFLAG = .TRUE.
!                    CALL FMTCSRC( ECS, 3, BUFFER, L )
! note: FMTCSRC doesn't work here because MODINFO is not populated
                WRITE( MESG,94010 )&
                &       'ERROR: Group ID', GID, 'in elevated ' //&
                &       'source ASCII input not in stack ' //&
                &       'groups file for:' // CRLF() // BLANK10 //&
                &       ECS
                CALL M3MSG2( MESG )
                CYCLE
            END IF

            COL  = GRPCOL( N )
            ROW  = GRPROW( N )
            XLOC = GRPXX ( N )
            YLOC = GRPYY ( N )

!.................  Skip records that are outside the grid
            IF ( COL .EQ. 0 .OR. ROW .EQ. 0 ) THEN

!                    CALL FMTCSRC( ECS, 3, BUFFER, L )
                WRITE( MESG,94010 ) 'NOTE: Source is outside the '//&
                &       'domain boundaries - excluding from UAM ' //&
                &       CRLF() // BLANK10 // 'elevated file for:' //&
                &       CRLF()// BLANK10// ECS
                CALL M3MESG( MESG )
                CYCLE

!.................  As per UAM's PTSRCE specifications, eliminate sources
!                   where the lat/lon coordinates are on the domain boundaries
            ELSE IF ( OUTFMT .EQ. 'UAM' .AND.&
            &        ( XLOC   .EQ. XMAX  .OR.&
            &          XLOC   .EQ. XORIG .OR.&
            &          YLOC   .EQ. YMAX  .OR.&
            &          YLOC   .EQ. YORIG       ) ) THEN

!                    CALL FMTCSRC( ECS, 3, BUFFER, L )
                WRITE( MESG,94010 ) 'NOTE: Source is on the ' //&
                &       'domain boundaries - excluding from UAM ' //&
                &       CRLF() // BLANK10 // 'elevated file for:' //&
                &       CRLF()// BLANK10// ECS
                CALL M3MESG( MESG )
                CYCLE

            END IF

!...............  Check to see if source has explicit plume rise and
!                 needs fake stacks
            IF ( EXPLFLAG ) M = FIND1( S, NHRSRC, ELEVSRC )
            IF ( M .GT. 0 ) THEN
                NEXPLOOP = EMLAYS
            ELSE
                NEXPLOOP = 1
            END IF

!...............  Loop through layers for explicit plume srcs (1 for
!                 standard sources)
            DO L = 1, NEXPLOOP

!.........................  Set stack height depending on whether the source
!                           has explicit plume rise or not
                IF( M .GT. 0 ) THEN
                    HT = LAYRMID( L )
                ELSE
                    HT = GRPHT ( N )
                END IF

!.........................  Set stack diameter depending on whether the source
!                           is a PinG source (for UAM-V & CAMx), or whether
!                           it's an explicit plume rise source (REMSAD).
                DM = GRPDM ( N )
                IF( ALLOCATED( LPING ) ) THEN
                    IF( LPING( S ) ) DM = -DM
                END IF
                IF( M .GT. 0 ) DM = -DM  ! Explicit source

!..................  If FIPS/plant/stack is not the same as the
!                    previous FIP/plant/stack, then store this unique
!                    list as needed for output to UAM format header
                IF( ECS .NE. PECS ) THEN

                    KK = KK + 1
                    EOUTCSRC( KK ) = ECS

!.........................  Set output stack parameters and coordinates
                    EOUTHT ( KK ) = HT
                    EOUTDM ( KK ) = DM
                    EOUTTK ( KK ) = GRPTK ( N )
                    EOUTVE ( KK ) = GRPVE ( N ) * 3600.    ! m/s -> m/hr
                    EOUTXL ( KK ) = XLOC
                    EOUTYL ( KK ) = YLOC

                END IF  ! end new FIP/plant/stack

!..................  For PinG sources, may need to overwrite stored diameter;
!                    first source in group may not be PinG, but later sources
!                    can be
                IF( DM < 0 ) EOUTDM ( KK ) = DM

!..................  Store arrays for writing emissions
                ELAYER ( I ) = L
                EIDX2  ( I ) = KK

            END DO      ! end loop over layers or 1

            PECS = ECS

        END DO          ! end of major sources

        NOUT = KK

        IF( EFLAG ) THEN
            MESG = 'Problem processing elevated sources.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( OUTFMT .EQ. 'UAM' ) THEN

!.................  Get settings from the environment...

!.................  Get type of plume rise calculation from environment
            MESG = 'Vertical method: "PLUMERISE" or "STACKHGT"'
            CALL ENVSTR( 'UAM_VERTICAL_METHOD', MESG, 'STACKHGT',&
            &             VTYPE, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_VERTICAL_METHOD"', 2 )
            END IF

!.................  Get UAM file note
            NOTEDEF = 'UAM elev pt emis from ' //&
            &          PROGNAME
            MESG = 'UAM file note for output file'
            CALL ENVSTR( 'UAM_NOTE', MESG, NOTEDEF , UNOTE, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_NOTE"', 2 )
            END IF

!.................  Insert grid name
            L = LEN_TRIM( GRDNM )
            FNOTE  = GRDNM( 1:L ) // ' ' // ADJUSTL( UNOTE )

!.................  Vertical structure information
            EML = EMLAYS
            MESG = 'Number of emission layers'
            IF( EML .EQ. 1 ) THEN
                EML = ENVINT( 'SMK_EMLAYS', MESG, 8, IOS )
                IF ( IOS .GT. 0 ) THEN
                    CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_EMLAYS"', 2 )
                END IF
            END IF

            MESG = 'Number of UAM vertical layers'
            NULAYS = ENVINT( 'UAM_LAYERS', MESG,&
            &                 EML , IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_LAYERS"', 2 )
            END IF

            NZLOWR = ENVINT( 'UAM4_LAYBELOW',&
            &                 'Layers below diffbreak', 3, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_LAYBELOW"', 2 )
            END IF

            NZUPPR = NULAYS - NZLOWR
            NZUPPR = ENVINT( 'UAM4_LAYABOVE',&
            &                 'Layers above diffbreak', NZUPPR, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_LAYABOVE"', 2 )
            END IF

            MESG = 'Height of surface layer [meters]'
            HTSUR = ENVREAL( 'UAM4_HTSFC', MESG, 0., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTSFC"', 2 )
            END IF

            MESG= 'Min. ht. of cells b/w sfc and diffbreak [meters]'
            HTLOWR = ENVREAL( 'UAM4_HTLOWR', MESG, 20., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTLOWR"', 2 )
            END IF

            MESG= 'Min. ht. of cells b/w diffbreak and top [meters]'
            HTUPPR = ENVREAL( 'UAM4_HTUPPR', MESG, 100., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_HTUPPR"', 2 )
            END IF

!.................  Control packet settings
            MESG = 'Number of parameters for UAM control packet'
            NPARAM = ENVINT( 'UAM_NPARAM', MESG, 20, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_NPARAM"', 2 )
            END IF

            MESG = 'PTSRCE output file unit number'
            PDEVOUT = ENVINT( 'UAM_PTSRCE_OUTUNIT', MESG, 20, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PTSRCE_OUTUNIT"', 2 )
            END IF

            MESG = 'Print output grid, 0=no, 1=yes'
            GSWITCH = ENVINT( 'UAM_PRINT_OUTGRD', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_OUTGRD"', 2 )
            END IF

            MESG = 'Print units table, 0=no, 1=yes'
            USWITCH = ENVINT( 'UAM_PRINT_UNITS', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_UNITS"', 2 )
            END IF

            MESG = 'Print source locations table, 0=no, 1=yes'
            LSWITCH = ENVINT( 'UAM_PRINT_LOCATIONS', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_LOCATIONS"', 2 )
            END IF

            MESG = 'Print methods table, 0=no, 1=yes'
            MSWITCH = ENVINT( 'UAM_PRINT_METHODS', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_METHODS"', 2 )
            END IF

            MESG = 'Print source values table, 0=no, 1=yes'
            VSWITCH = ENVINT( 'UAM_PRINT_VALUES', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_VALUES"', 2 )
            END IF

            MESG = 'Print vertical methods table, 0=no, 1=yes'
            ESWITCH = ENVINT( 'UAM_PRINT_VERTMETH', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM_PRINT_VERTMETH"', 2 )
            END IF

            MESG = 'DIFFBREAK file unit number'
            DDEVOUT = ENVINT( 'UAM4_DIFFBREAK_UNIT', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_DIFFBREAK_UNIT"', 2 )
            END IF

            MESG = 'REGIONTOP file unit number'
            RDEVOUT = ENVINT( 'UAM4_REGIONTOP_UNIT', MESG, 0, IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_REGIONTOP_UNIT"', 2 )
            END IF

!.................  If PLUMERISE method has been selected, get file unit
!                   assignments
            IF( VTYPE .EQ. 'PLUMERISE' ) THEN

                MESG = 'TEMPERATUR file unit number'
                TDEVOUT= ENVINT( 'UAM4_TEMPERATUR_UNIT',&
                &                 MESG, 14, IOS )
                IF ( IOS .GT. 0 ) THEN
                    CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_TEMPERATUR_UNIT"', 2 )
                END IF

                MESG = 'METSCALARS file unit number'
                MDEVOUT= ENVINT( 'UAM4_METSCALARS_UNIT',&
                &                 MESG, 15, IOS )
                IF ( IOS .GT. 0 ) THEN
                    CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_METSCALARS_UNIT"', 2 )
                END IF

                MESG = 'WIND file unit number'
                WDEVOUT = ENVINT( 'UAM4_WIND_UNIT', MESG, 16, IOS )
                IF ( IOS .GT. 0 ) THEN
                    CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "UAM4_WIND_UNIT"', 2 )
                END IF

            END IF

!.................  Get output units from the environment
            BUFFER = 'Units for output gridded emissions'
            CALL ENVSTR( 'MRG_GRDOUT_UNIT', BUFFER,&
            &               ' ', GRDENV, IOS)
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_GRDOUT_UNIT"', 2 )
            END IF

!.................  Write header information to the file
            WRITE( EVDEV,93010 ) 'CONTROL', GRDENV
            WRITE( EVDEV,93000 ) 'PTSOURCE'
            WRITE( EVDEV,93000 ) FNOTE

!.................  Insert crrect number of sources into 1st line of ESTRING
            WRITE( EVDEV,93015 ) NMSPC, 0, NOUT, 1, NPARAM
            WRITE( EVDEV,93015 ) PDEVOUT, 0, GSWITCH
            WRITE( EVDEV,93015 ) USWITCH, LSWITCH, 0, MSWITCH,&
            &                     VSWITCH
            WRITE( EVDEV,93015 ) 1, 0, ESWITCH
            WRITE( EVDEV,93015 ) DDEVOUT, RDEVOUT, 0, TDEVOUT,&
            &                     MDEVOUT, WDEVOUT

!.................  Write species names
            DO I = 1, NMSPC
                SPCNAM = EMNAM( I )( 1:LEN( SPCNAM ) )

                IF( SPCNAM .NE. EMNAM( I ) ) THEN
                    L = LEN_TRIM( EMNAM( I ) )
                    MESG = 'WARNING: Species name "' //&
                    &       EMNAM( I )( 1:L ) //&
                    &       '" is being truncated to "' // SPCNAM //&
                    &       CRLF() // BLANK10 // '" on output to ' //&
                    &       'UAM elevated point sources file'
                    CALL M3MSG2( MESG )

                END IF

                WRITE( EVDEV, 93020, IOSTAT=IOS ) SPCNAM

            END DO

!.................  Output the beginning and ending dates and times
            IBD = REMOVE_4DIGIT_YEAR( SDATE )
            IED = REMOVE_4DIGIT_YEAR( EDATE )
            IBT = STIME / 100
            IET = ETIME / 100
            WRITE( EVDEV, 93030 ) IBD, IBT, IED, IET

            WRITE( EVDEV, 93000 ) 'END'

!.................  Region packet:
            WRITE( EVDEV, 93000 ) 'REGION'
            WRITE( EVDEV, 93040 ) 0., 0., INT( P_ALP )
            WRITE( EVDEV, 93040 ) XORIG, YORIG

            IF( GDTYP .EQ. LATGRD3 ) THEN           ! Print out for Lat-Lon
                WRITE( EVDEV, 93045 ) XCELL, YCELL
            ELSE                                    ! other coordinate sys
                WRITE( EVDEV, 93040 ) XCELL, YCELL
            END IF

            WRITE( EVDEV, 93050 ) NCOLS, NROWS, NULAYS
            WRITE( EVDEV, 93060 ) NZLOWR, NZUPPR, HTSUR,&
            &                      HTLOWR, HTUPPR
            WRITE( EVDEV, 93000 ) 'END'

!.................  Point Sources packet:

            WRITE( EVDEV, 93000 ) 'POINT SOURCES'

            DO I = 1, NOUT

                ECS  = EOUTCSRC( I )
                CFIP = ECS( 1:FIPLEN3 )
                FCID = ADJUSTL( ECS( PTBEGL3( 2 ):PTENDL3( 2 ) ) )
                SKID = ADJUSTL( ECS( PTBEGL3( 3 ):PTENDL3( 3 ) ) )

                L = LEN_TRIM( FCID )
                IF( L .GT. 10 ) THEN
                    MESG = 'WARNING: Plant name truncated from "' //&
                    &       FCID( 1:L ) // '" to "' // FCID( 1:10 )//&
                    &       '" on output to ASCII elevated file'
                    CALL M3MESG( MESG )
                END IF

                L = LEN_TRIM( SKID )
                IF( L .GT. 10 ) THEN
                    MESG = 'WARNING: Stack name truncated from "' //&
                    &       SKID( 1:L ) // '" to "' // SKID( 1:10 )//&
                    &       '" on output to ASCII elevated file'
                    CALL M3MESG( MESG )
                END IF

                FIP  = STR2INT( CFIP )
                FCID = ADJUSTR( FCID( 1:10 ) )
                SKID = ADJUSTR( SKID( 1:10 ) )

!.....................  Output source informat in a different format if
!                       the grid is lat-lon or not
                IF( GDTYP .EQ. LATGRD3 ) THEN
                    WRITE( EVDEV, 93065 )&
                    &    I, 'STD       ', EOUTXL( I ), EOUTYL( I ),&
                    &    FCID( 1:10 ), SKID( 1:10 ), FIP

                ELSE
                    WRITE( EVDEV, 93070, IOSTAT=IOS )&
                    &    I, 'STD       ', EOUTXL( I ), EOUTYL( I ),&
                    &    FCID( 1:10 ), SKID( 1:10 ), FIP

                END IF

                WRITE( EVDEV, 93080, IOSTAT=IOS )&
                &       EOUTHT ( I ), EOUTDM ( I ),&
                &       EOUTTK ( I ), EOUTVE ( I )

            END DO

            WRITE( EVDEV, 93000 ) 'END'

        END IF    ! End of output format for UAM file

        FIRSTIME = .FALSE.

    END IF        ! end if firstime routine is called

!.........  Write out date/time-specific data...
!.........  Data previously merged in MRGELEV
    IF( JTIME .NE. PTIME .AND. OUTFMT .EQ. 'UAM' ) THEN

!.............  Time Interval packet:
!.............  Set required time parameters
        JDATEP1 = JDATE
        JTIMEP1 = JTIME
        CALL NEXTIME( JDATEP1, JTIMEP1, TSTEP )

        IBT = JTIME   / 100
        IET = JTIMEP1 / 100
        IBD = REMOVE_4DIGIT_YEAR( JDATE )
        IED = REMOVE_4DIGIT_YEAR( JDATEP1 )

        WRITE( EVDEV, 93000 ) 'TIME INTERVAL'
        WRITE( EVDEV, 93050 ) IBD, IBT, IED, IET

!.............  Method packet:
!.............  Provide same output as PSTPNT

        WRITE( EVDEV, 93000 ) 'METHOD'
        WRITE( EVDEV, 93000 )&
        &      'STD       ALL       EMVALUES  0.        50000.'
        WRITE( EVDEV, 93000 ) 'END'

!.............  Vertical Method packet:
!.............  Provide same output as PSTPNT with variable plume-rise method

        WRITE( EVDEV, 93000 ) 'VERTICAL METHOD'
        WRITE( EVDEV, 93000 )&
        &       'STD       ALL       ' // VTYPE // ' 0.       10000.'
        WRITE( EVDEV, 93000 ) 'END'

    END IF    ! End UAM/CAMx elevated point sources format

!.........  But first determine how many explicit sources
!           there are for the current hour.
    DO N = 1, NHRSRC
        IF ( INDXH( N ) .EQ. 0 ) EXIT
    END DO
    NH = N - 1

!.........  Write out emissions data for UAM/CAMx elevated point sources

    IF( OUTFMT .EQ. 'UAM' ) THEN

!.............  Emissions Values packet.  Only write this if the time is
!               new, otherwise, it will write each time routine is called
!               for the next species.
        IF( JTIME .NE. PTIME ) THEN
            WRITE( EVDEV, 93000 ) 'EMISSIONS VALUES'
            WRITE( EVDEV, 93000 ) 'ALL       ALL            0.000'
        END IF

        LN = 0
        LM = 0
        LL = 0
        M  = 0
        ESUM = 0.  ! array
        PECS = ' '
        DO I = 1, NMAJOR

            J   = ESRTIDX ( I )  ! sort for UAM output index
            L   = ELAYER  ( I )  ! layer, 1 or 0 if outside grid
            N   = EIDX2   ( I )  ! UAM elev src no. or 0 if outside grd

            IF ( N .LE. 0 ) CYCLE   ! skip because not in grid

            S   = EIDXSRC( J )
            ECS = ECSRCA ( J )

            IF ( EXPLFLAG ) M = FIND1( S, NHRSRC, ELEVSRC )

! NOTE: This algorithm will not work if there are explicit elevated sources
!    N: that have duplicates.

!.................  If FIP/plant/stack is different from last time,
!                   write emissions for previous elevated output src
!                   and initialize
            IF ( PECS .NE. ECS .OR. L .NE. LL ) THEN

!.....................  Write emissions for standard output srcs
                IF( LM .LE. 0 ) THEN

!......................  Only write emissions if not zero
                    IF( ESUM( 1 ) .NE. 0. ) THEN
                        CALL GET_ESUM_FORMAT( VNAME, ESUM(1), EFMT )
                        WRITE( EVDEV, EFMT ) LN, VNAME, ESUM(1)
                    END IF
                    ESUM( 1 ) = 0.

!.....................  Write emissions for explicit source records
                ELSE

                    IF( ESUM( LL ) .NE. 0. ) THEN
                        CALL GET_ESUM_FORMAT( VNAME,ESUM(LL),EFMT )
                        WRITE( EVDEV,EFMT ) LN, VNAME, ESUM( LL )
                    END IF
                    ESUM( L ) = 0.

                END IF

!.....................  Initialize for standard FIP/plant/stack
                IF( M .LE. 0 ) THEN
                    ESUM( 1 ) = ELEVEMIS( J )

!.....................  Initialize for explicit sources
                ELSE
                    ESUM( L ) = ELEVEMIS( J ) * LFRAC( M,L )

                END IF

!.................  If FIP/plant/stack is the same as the previous
!                   iteration, sum emissions
            ELSE

!.....................  For standard sources
                IF( M .LE. 0 ) THEN

                    ESUM( 1 ) = ESUM( 1 ) + ELEVEMIS( J )

!.....................  For explicit sources
                ELSE
                    ESUM(L) = ESUM(L) + ELEVEMIS(J) * LFRAC(M,L)

                END IF

            END IF

            LN = N       ! index to grouped UAM source number
            LL = L       ! layer number or 1
            LM = M
            PECS = ECS

        END DO

!.............  Output for final record(s)...
!.............  For standard sources
        IF( M .LE. 0 ) THEN

            IF( ESUM( 1 ) .NE. 0. ) THEN
                CALL GET_ESUM_FORMAT( VNAME, ESUM(1), EFMT )
                WRITE( EVDEV, EFMT ) NOUT, VNAME, ESUM(1)
            END IF

!.............  For explicit sources
        ELSE
            IF( ESUM( L ) .NE. 0. ) THEN
                CALL GET_ESUM_FORMAT( VNAME, ESUM( L ), EFMT )
                WRITE( EVDEV,EFMT ) NOUT, VNAME, ESUM( L )
            END IF

        END IF

!.............  After last species, write packet end fields
        IF( VNAME .EQ. EMNAM( NMSPC ) ) THEN
            WRITE( EVDEV, 93000 ) 'END'
            WRITE( EVDEV, 93000 ) 'ENDTIME'
        END IF

    END IF   ! End UAM/CAMx elevated point sources format

    PTIME = JTIME

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93010 FORMAT( A7, 3X, A )

93015 FORMAT( 6I10 )

93020 FORMAT( A10 )

93030 FORMAT( 4I10 )

93040 FORMAT( 2F10.0, I10 )

93045 FORMAT( 2F10.7 )

93050 FORMAT( 4I10 )

93060 FORMAT( 2I10, 3F10.0 )

93065 FORMAT( I10, A10, F10.5, F10.5, 2A10, I10.5 )

93070 FORMAT( I10, A10, F10.0, F10.0, 2A10, I10.5 )

93080 FORMAT( F10.1, F10.2, F10.1, F10.0 )

93090 FORMAT( I10, A10, F10.3 )

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94300 FORMAT( A, I10, A )

!*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

!.............  This internal subprogram tries to read a variable from an
!               I/O API file, and aborts if not successful.
    SUBROUTINE SAFE_READ3( FILNAM, VARNAM, LAYER,&
    &                       JDATE, JTIME, XBUF     )

!.............  Subprogram arguments
        CHARACTER(*) FILNAM    ! logical file name
        CHARACTER(*) VARNAM    ! variable name
        INTEGER      LAYER     ! layer number (or ALLAYS3)
        INTEGER      JDATE     ! Julian date
        INTEGER      JTIME     ! time
        REAL         XBUF( * ) ! read buffer

!----------------------------------------------------------------------

        IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,&
        &                  JDATE, JTIME, XBUF ) ) THEN

            MESG = 'Could not read "' // TRIM( VARNAM ) //&
            &       '" from file "' // TRIM( FILNAM ) // '."'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

    END SUBROUTINE SAFE_READ3

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.............  This internal subprogram tries to read a variable from an
!               I/O API file, and aborts if not successful.
    SUBROUTINE SAFE_READ3_I( FILNAM, VARNAM, LAYER,&
    &                        JDATE, JTIME, IBUF     )

!.............  Subprogram arguments
        CHARACTER(*) FILNAM    ! logical file name
        CHARACTER(*) VARNAM    ! variable name
        INTEGER      LAYER     ! layer number (or ALLAYS3)
        INTEGER      JDATE     ! Julian date
        INTEGER      JTIME     ! time
        INTEGER      IBUF( * ) ! read buffer

!.............  Local variables
        INTEGER L1, L2

!----------------------------------------------------------------------

        IF ( .NOT. READ3( FILNAM, VARNAM, LAYER,&
        &                  JDATE, JTIME, IBUF ) ) THEN

            MESG = 'Could not read "' // TRIM( VARNAM ) //&
            &       '" from file "' // TRIM( FILNAM ) // '."'
            CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

        END IF

    END SUBROUTINE SAFE_READ3_I

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.............  This internal subprogram removes the first two digits of the
!               year from a 7-digit Julian date
    INTEGER FUNCTION REMOVE_4DIGIT_YEAR( JDATE )

!.............  Subroutine arguments
        INTEGER, INTENT (IN) :: JDATE

!.............  Local variables
        INTEGER    YRREMOVE

!----------------------------------------------------------------------

        YRREMOVE = ( JDATE / 100000 ) * 100000  ! integer math
        REMOVE_4DIGIT_YEAR = JDATE - YRREMOVE

        RETURN

    END FUNCTION REMOVE_4DIGIT_YEAR

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.............  This internal subprogram determines the format to output
!               emission values
    SUBROUTINE GET_ESUM_FORMAT( VBUF, VAL, FMT )

!.............  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: VBUF
        REAL        , INTENT (IN) :: VAL
        CHARACTER(*), INTENT(OUT) :: FMT

!----------------------------------------------------------------------

        IOS = 0

!.............  Value is too large for
        IF( VAL .GT. 999999999. ) THEN
            FMT = '( I10, A10, E10.3 )'

            WRITE( MESG,94020 ) 'WARNING: "' // TRIM( VBUF ) //&
            &  '"Emissions value of', ESUM, CRLF()// BLANK10//&
            &  '" is too large for file format, so writing ' //&
            &  'in scientific notation for source:',&
            &  CRLF() // BLANK10 // ECSRCA( J )

            CALL M3MESG( MESG )

        ELSE IF( VAL .GT. 99999999. ) THEN
            FMT = '( I10, A10, F10.0 )'

        ELSE IF( VAL .GT. 9999999. ) THEN
            FMT = '( I10, A10, F10.1 )'

        ELSE IF( VAL .GT. 999999. ) THEN
            FMT = '( I10, A10, F10.2 )'

        ELSE IF( VAL .GT. 0. .AND. VAL .LT. 1. ) THEN
            FMT = '( I10, A10, E10.4 )'

        ELSE
            FMT = '( I10, A10, F10.3 )'

        END IF

        RETURN

94020   FORMAT( 10( A, :, E12.5, :, 1X ) )

    END SUBROUTINE GET_ESUM_FORMAT

END SUBROUTINE WMRGELEV
