
SUBROUTINE RDGRDNCF( FNAME, GNAME )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads an native (raw)  NetCDF gridded inventory files and
    !      GRIDMASK input file that holds Geocode and timezone by grid cell.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 02/2016 by B.H. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: POLVAL, TZONES, CIFIP, CELLID, TPFLAG, INVYR,    &
                        NPCNT, CSCC, IPOSCOD, CSOURC

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPOL, NIPPA, NPPOL, CATEGORY, NEM, NDY,    &
                       EIIDX, EINAM, EANAM, EAUNIT, EADESC, NSRC

    USE M3UTILIO
    USE MODNCFIO

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME         ! logical name of input file
    CHARACTER(*), INTENT (IN) :: GNAME         ! logical name of input file

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETIFDSC, GETFORMT, GETFLINE

    !.......   Local allocatable arrays
    REAL,               ALLOCATABLE :: EMIS( : )      ! tmp gridded emissions data
    INTEGER,            ALLOCATABLE :: INXLST( : )     ! sort index for LST files
    INTEGER,            ALLOCATABLE :: IDXFIP( : )     ! sort index for FIPS
    INTEGER,            ALLOCATABLE :: IFIP( : )      ! tmp gridded input data
    INTEGER,            ALLOCATABLE :: ZONE( : )      ! tmp gridded time zone
    CHARACTER(NAMLEN3), ALLOCATABLE :: CPOL( : )      ! tmp gridded pol names
    CHARACTER(SCCLEN3), ALLOCATABLE :: TSCC( : )      ! tmp gridded sccs
    CHARACTER(300),     ALLOCATABLE :: SORTLST( : )
    CHARACTER(300),     ALLOCATABLE :: EDGARLST( : )
    CHARACTER(SCCLEN3+NAMLEN3), ALLOCATABLE :: SCCPOL( : )

    CHARACTER(256)  SEGMENT( 5 )

    !.......   Other local variables
    INTEGER         C, ES, I, J, K, L, N, NF, S0, S, V         !  counters and indices

    INTEGER         FLEN            !  file name length
    INTEGER         FDEV            !  unit no for ARINV
    INTEGER         IDX             !  emissions array index (ann or ave day)
    INTEGER      :: INY = 0         !  tmp inventory year
    INTEGER         IOS             !  i/o status
    INTEGER         INVFMT          !  inventory format code
    INTEGER      :: NDAYS = 0       !  days of year
    INTEGER         NSCC            !  number of SCCs
    INTEGER         NLINE           !  number of lines in list format file
    INTEGER      :: NCELL           !  tmp cell numbers
    INTEGER      :: TPF = 0         !  tmp temporal adjustments setting
    INTEGER         WKSET           !  setting for wkly profile TPFLAG component
    INTEGER         NCOLS, NROWS      !  no. variables in gridded input file

    REAL         :: CNVFAC = 0.0      !  conversion factor

    LOGICAL      :: ANNFLAG= .FALSE.     ! true: annual inventory
    LOGICAL      :: EFLAG  = .FALSE.     ! true: error occurred
    LOGICAL      :: DFLAG  = .FALSE.     ! true: weekday (not full week) nrmlizr

    CHARACTER(300)  MESG, LINE      !  message buffer

    CHARACTER( 2 )     CMON         ! monthly inv flag
    CHARACTER(CELLEN3) CCELL        ! tmp cell ID
    CHARACTER(POLLEN3) CCOD         ! character pollutant index
    CHARACTER(SCCLEN3) PSCC, SCC          ! source category code
    CHARACTER(NAMLEN3) CBUF, CVAR         ! tmp variable name

    CHARACTER(16), PARAMETER :: PROGNAME =  'RDGRDNCF'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDGRDNCF

    !.......  Set default inventory characteristics (declared in MODINFO)
    CALL INITINFO( IOGFMT )

    !.......  Retrieve raw netcdf inventory base year
    MESG = 'Define the year of NetCDF gridded inventory'
    INY = ENVINT( 'NETCDF_INV_YEAR', MESG, 0, IOS )
    IF( INY == 0 ) THEN
        MESG = 'ERROR: MUST define the year of gridded NetCDF inventory file [NETCDF_INV_YEAR]'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        NDAYS = INT( 1 / YR2DAY ( INY ) )
    END IF

    MESG = 'Enter logical name for the inventory list file'
    FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., FNAME, PROGNAME )

    !.......  Determine file format of inventory file
    INVFMT = GETFORMT( FDEV, -1 )

    !.......  If SMOKE list format, read file and check file for formats.
    !           NOTE- LSTFMT defined in EMCNST3.h90
    IF( INVFMT == LSTFMT ) THEN

        !.......  Generate message for GETFLINE and RDLINES calls
        MESG = TRIM( CATEGORY ) // ' inventory file, ' //    &
               TRIM( FNAME ) // ', in list format'

        !.......  Get number of lines of inventory files in list format
        NLINE = GETFLINE( FDEV, MESG )

        !.......  Allocate memory for storing contents of list-formatted file
        ALLOCATE( SORTLST( NLINE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SORTLST', PROGNAME )
        SORTLST = ' '         ! array

        !.......  Store lines of inventory list file
        CALL RDLINES( FDEV, FNAME, NLINE, SORTLST )

        NF = 0
        DO I = 1, NLINE
            IF( BLKORCMT( SORTLST(I) ) ) CYCLE
            NF = NF + 1
        END DO

        !.......  Allocate local arrays
        ALLOCATE( INXLST( NF ),    &
                    TSCC( NF ),    &
                    CPOL( NF ),    &
                  SCCPOL( NF ),    &
                EDGARLST( NF ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INXLST...EDGARLST', PROGNAME )
        INXLST = 0
        TSCC = ' '
        CPOL = ' '
        SCCPOL = ' '
        EDGARLST = ' '

        !.......  Store sorted input list
        NF = 0
        DO I = 1, NLINE

            IF( BLKORCMT ( SORTLST( I ) ) ) CYCLE

            READ( SORTLST( I ), * ) SCC, CBUF, CVAR, CMON, LINE
            CALL PADZERO( SCC )

            NF = NF + 1
            INXLST( NF ) = NF
            SCCPOL( NF ) = SCC // CBUF
            EDGARLST( NF ) = SORTLST( I )

        END DO

        CALL SORTIC( NF, INXLST, SCCPOL )

        !.......  count no of SCCs and Pollutants
        NSCC = 0
        NIPOL = 0
        DO I = 1, NF

            J = INXLST( I )

            READ( EDGARLST( J ), * ) SCC, CBUF, CVAR, CMON, LINE

            IF( INDEX1( CBUF, NF, CPOL ) < 1 ) THEN
                NIPOL = NIPOL + 1
                CPOL( NIPOL ) = CBUF
            END IF

            IF( INDEX1( SCC, NF, TSCC ) < 1 ) THEN
                NSCC = NSCC + 1
                TSCC( NSCC ) = SCC
            END IF

        END DO

        IF( CMON == '0' .OR. CMON == 'N' .OR. CMON == ' ' ) THEN
            ANNFLAG = .TRUE.       ! Determine inv temp resolution (annual or monthly)
        END IF

        NIPPA = NIPOL

        CLOSE( FDEV )

    END IF

    !.......  Allocate memory for variables
    !.......  NOTE - Both EINAM and EANAM are created to support WRINVEMIS
    ALLOCATE( EIIDX( NIPOL ),    &
              EINAM( NIPOL ),    &
              EANAM( NIPPA ),    &
             EAUNIT( NIPPA ),    &
             EADESC( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EIIDX...EADESC', PROGNAME )

    DO V = 1, NIPOL
        EIIDX ( V ) = V
        EINAM ( V ) = CPOL( V )
        EANAM ( V ) = CPOL( V )
        IF( ANNFLAG ) THEN
            EAUNIT( V ) = 'tons/yr'
            EADESC( V ) = 'Annual Emissions for ' // TRIM( CPOL(V) )
        ELSE
            EAUNIT( V ) = 'tons/day'
            EADESC( V ) = 'Average-day Emissions for ' // TRIM( CPOL(V) )
        END IF
    END DO

    DEALLOCATE( CPOL, TSCC, SORTLST )

    !.......  Read the header of the GRIDMASK input file

    IF( .NOT. DESC3( GNAME ) ) THEN
        MESG = 'Could not get description of file "' // TRIM( GNAME ) // '"'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Give error if there are no variables
    IF( NVARS3D .LT. 1 ) THEN
        MESG = 'No variables found in gridded input file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Give warning if layers are greater than 1
    IF( NLAYS3D .GT. 1 ) THEN
        WRITE( MESG,94010 ) 'WARNING: Only the first layer out of'    &
              , NLAYS3D, 'will be imported.'
        CALL M3MSG2( MESG )
    END IF

    !.......  Set the number of cells, variables, grid name, etc.
    NCOLS = NCOLS3D
    NROWS = NROWS3D
    NCELL = NROWS * NCOLS
    NSRC  = NCELL * NSCC

    !.......  Allocate local memory for temporary gridded data file
    ALLOCATE( IFIP( NCELL ),    &
              ZONE( NCELL ),    &
            IDXFIP( NCELL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IFIP...IDXFIP', PROGNAME )
    IFIP = 0
    ZONE = 0

    !.......  Read local variables
    IF( .NOT. READ3( GNAME, 'GEOCODE', 1, 0, 0, IFIP ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not read GEOCODE from file ' // TRIM( GNAME ) // '".'
        CALL M3MSG2( MESG )
    END IF

    !.......  Read local variables
    IF( .NOT. READ3( GNAME, 'TZONES', 1, 0, 0, ZONE ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not read TZONES from file ' // TRIM( GNAME ) // '".'
        CALL M3MSG2( MESG )
    END IF

    !.......  Sort IFIP list
    IDXFIP = 0

    DO I = 1, NCELL
        IDXFIP( I ) = I
    END DO

    CALL SORTI( NCELL, IDXFIP, IFIP )

    !.......  Allocate memory for (sorted) output inventory characteristics.
    !.......  The sorted arrays can be allocated right away because the only
    !           source characteristic for this type of data is the grid cell, and
    !           it is sorted already.
    CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .FALSE., NSRC, NSRC*NIPOL, NPPOL )

    CALL SRCMEM( CATEGORY, 'SORTED', .TRUE., .TRUE.,  NSRC, NSRC*NIPOL, NPPOL )

    !.......  Get setting for interpreting weekly temporal profiles from the
    !           environment.
    DFLAG = .FALSE.
    MESG = 'Use weekdays only to normalize weekly profiles'
    DFLAG = ENVYN( 'WKDAY_NORMALIZE', MESG, DFLAG, IOS )

    !.......  Set weekly profile interpretation flag...
    !.......  Weekday normalized
    IF( DFLAG ) THEN
        WKSET = WDTPFAC
        MESG = 'NOTE: Setting inventory to use weekday '//    &
               'normalizer for weekly profiles'
    !.......  Full-week normalized
    ELSE
        WKSET = WTPRFAC
        MESG = 'NOTE: Setting inventory to use full-week '//    &
               'normalizer for weekly profiles'
    END IF

    !.......  Write message
    CALL M3MSG2( MESG )

    !.......  Initialize emissions data
    ALLOCATE( EMIS( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMIS', PROGNAME )
    EMIS = 0.0

    POLVAL = BADVAL3        ! array

    S = 0
    S0 = 0
    ES = 0
    PSCC = ' '
    !.......  Loop over pollutant-specific raw NetCDF gridded inventory file
    DO I = 1, NF

        J = INXLST( I )

        CALL PARSLINE( EDGARLST( J ), 5, SEGMENT )

        SCC  = TRIM( SEGMENT( 1 ) )
        CBUF = TRIM( SEGMENT( 2 ) )
        CVAR = TRIM( SEGMENT( 3 ) )
        CMON = TRIM( SEGMENT( 4 ) )

        CALL PADZERO( SCC )

        V = INDEX1( CBUF, NIPOL, EANAM )

        !.......  Define inventory temporal resoltion (annual or avg-day)
        IF( CMON == '0' .OR. CMON == 'N' .OR. CMON == ' ' ) THEN           ! Annual inventory
            TPF = MTPRFAC * WKSET
            IDX = NEM
            CNVFAC = 1000. * GM2TON * NDAYS * DAY2SEC          ! kg/m2/s -> tons/m2/year
        ELSE
            TPF = WKSET
            IDX = NDY
            CNVFAC = 1000. * GM2TON * DAY2SEC          ! kg/m2/s -> tons/m2/day
        END IF

        !.......  Set output logical file name
        IF( .NOT. SETENVVAR( 'TMPFILE', SEGMENT( 5 ) ) ) THEN
            MESG = 'ERROR: Could not set logical name ' //    &
               'for native gridded NetCDF inventory file:'//    &
               CRLF()// BLANK10// TRIM( LINE )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. READNCVAR( 'TMPFILE', CVAR, NCOLS, NROWS, EMIS ) ) THEN
            EFLAG = .TRUE.
            MESG = 'Could not read "' // TRIM( CVAR ) // '" from gridded inventory file'
            CALL M3MESG( MESG )
        END IF

        IF( SCC /= PSCC ) THEN

            S0 = S + 1

            DO J = 1, NCELL

                C = IDXFIP( J )

                S = S + 1
                CELLID( S ) = C
                INVYR ( S ) = INY
                NPCNT ( S ) = NIPPA
                CSCC  ( S ) = SCC
                TPFLAG( S ) = TPF
                TZONES( S ) = ZONE( C )
                WRITE( CIFIP( S ), '( I6.6,I6.6 )' ) IFIP( C ),0
                WRITE( CCELL, 94130 ) CELLID( S )

                CALL BLDCSRC( CIFIP( S ), SCC, CCELL, CHRBLNK3,    &
                              CHRBLNK3, CHRBLNK3, CHRBLNK3, CBUF, CSOURC( S ) )
            END DO

            PSCC = SCC

        END IF

        !.......  Read data from gridded file and store in appropriate data structure
        !               for use by the rest of the programs
        C = 0
        DO L = S0, S
            C = C + 1                                   ! increment cell ID
            ES = ( L-1 ) * NIPPA + V
            IPOSCOD( ES )     = V
            POLVAL ( ES,IDX ) = CNVFAC * EMIS( IDXFIP( C ) )      ! convert kg/m2/s to tons/m2/yr-hr
        END DO

    END DO

    !.......  Abort if there was a reading error
    IF( EFLAG ) THEN
        MESG = 'Problem reading gridded inventory file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94125 FORMAT( I5 )

94130 FORMAT( I8 )

END SUBROUTINE RDGRDNCF
