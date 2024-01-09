
PROGRAM SPCMAT

    !***********************************************************************
    !  program body starts at line
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Copied from spcpmat.F 1/99 by M. Houyoux 
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90", and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !    System
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
    !************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains the speciation profiles
    USE MODSPRO, ONLY: MXSPEC, MXSPFUL, SPCNAMES, INPRF, SPECID,    &
                       MASSFACT, MOLEFACT, MOLUNITS, NSPFUL,        &
                       SPROFN, IDXSPRO, NSPECIES, IDXSSPEC,         &
                       NPOLSPRO, NSPCALL, SPCLIST

    !.......  This module contains emission factor tables and related
    USE MODEMFAC, ONLY: INPUTHC, OUTPUTHC, EMTNAM, EMTPOL, NEPOL, NETYPE

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CRL, NSRC, NIACT, NIPPA, NIPOL, EANAM, EINAM

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDVTS

    !.......  This module contains the tagging arrays
    USE MODTAG, ONLY: TAGNUM, TAGNAME, MXTAG

    !.......  This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE

    !.......  LOCAL PARAMETERS and their descriptions:

    CHARACTER(1 ), PARAMETER :: QUOTE    = "'"
    CHARACTER(16), PARAMETER :: PROGNAME = 'SPCMAT'     ! program name
    CHARACTER(50), PARAMETER ::     CVSW = '$Name SMOKEv5.0_Jun2023$'      ! CVS revision tag

    !.......   LOCAL VARIABLES and their descriptions:

    !.......  Array that contains the names of the inventory variables needed for
    !    this program
    CHARACTER(NAMLEN3) IVARNAMS( MXINVARR )

    !.......   Speciation matrices:

    REAL, ALLOCATABLE :: MASSMATX( :,: )        !  mass-speciation coefficients
    REAL, ALLOCATABLE :: MOLEMATX( :,: )        !  mole-speciation coefficients

    !.......  Local output arrays for mole- and mass-based matrices.
    REAL, ALLOCATABLE :: MASSARR( : )        !  mass-speciation for write
    REAL, ALLOCATABLE :: MOLEARR( : )        !  mole-speciation for write

    !.......  Local tagging arrays
    INTEGER, ALLOCATABLE :: SRCTAGS( : )     ! tag number by source

    !.......  Inventory pollutants actually in the inventory
    LOGICAL           , ALLOCATABLE :: SPCOUT( : )     ! true: output spcs
    LOGICAL           , ALLOCATABLE :: IDXCHK( : )     ! true: EAIDX value accounted for
    LOGICAL           , ALLOCATABLE :: LEMIS ( : )     ! true: emissions, not activity-based
    CHARACTER(NAMLEN3), ALLOCATABLE :: IINAM ( : )     ! initial pols
    CHARACTER(NAMLEN3), ALLOCATABLE :: SINAM ( : )     ! output pollutants

    !.......  Activities and pollutants
    INTEGER, ALLOCATABLE :: EAIDX( : )               ! Index from EANAM to IINAM
    CHARACTER(NAMLEN3), ALLOCATABLE :: SANAM ( : )     ! output pollutants

    !.......  Names for output variables in mass-based and mole-based files
    CHARACTER(NAMLEN3), ALLOCATABLE :: MASSONAM( :,:,: )
    CHARACTER(NAMLEN3), ALLOCATABLE :: MOLEONAM( :,:,: )

    !.......  Unit numbers and logical file names
    INTEGER      :: GDEV = 0    !  unit number for tagging file
    INTEGER         IDEV        !  tmp unit number if ENAME is map file
    INTEGER      :: KDEV = 0    !  unit no. for optional pol to pol conversion
    INTEGER         LDEV        !  unit number for log file
    INTEGER      :: MDEV = 0    !  unit number for mobile codes file
    INTEGER         RDEV        !  unit number for speciation profiles file
    INTEGER         SDEV        !  unit number for ASCII inventory file
    INTEGER         TDEV        !  unit number for ASCII emission process file
    INTEGER         UDEV        !  unit number for ASCII supplemental file
    INTEGER         VDEV        !  unit no. for inventory data table
    INTEGER         XDEV        !  unit no. for cross-reference file

    CHARACTER(16)   ANAME       !  logical name for additive control matrix
    CHARACTER(16)   ENAME       !  logical name for point-source input file
    CHARACTER(16)   INAME       !  tmp name for inven file of unknown fmt
    CHARACTER(16)   SNAME       !  logical name for mass spec matrix output file
    CHARACTER(16)   LNAME       !  logical name for mole spec matrix output file

    !.......   Other local variables
    INTEGER          I, J, K, L, L1, L2, L3, L4, LT, N, S, T, V     !  counters and indices

    INTEGER          IDX                   ! tmp index value
    INTEGER          IOS                   ! i/o status
    INTEGER          NINVARR               ! number inventory variables to input
    INTEGER          NMSPC                 ! number of model species
    INTEGER          NOPOL                 ! no. pollutants for output
    INTEGER       :: NP    = 0             ! no. all pollutants (size of IINAM)
    INTEGER          NTAG                  ! tmp number of tags
    INTEGER          PIDX                  ! previous index value
    INTEGER          TAGCNT                ! count of sources with tags for a species

    LOGICAL       :: EFLAG   = .FALSE.     !  error flag
    LOGICAL       :: KFLAG   = .FALSE.     !  if pol to pol convert file or not
    LOGICAL       :: TFLAG   = .FALSE.     !  true: tagging file needed
    LOGICAL       :: MASSOUT = .TRUE.      !  true: output mass-based matrix
    LOGICAL       :: MOLEOUT = .TRUE.      !  true: output mole-based matrix
    LOGICAL       :: DEFREPRT= .TRUE.      !  true: report default spc profiles
    LOGICAL       :: MULTIPRO= .TRUE.      !  true: multiple profs for pollutant
    LOGICAL       :: FNDOUTPUT=.FALSE.     !  true: found output hydrocarbon

    CHARACTER(300)     MESG          !  message buffer
    CHARACTER(NAMLEN3) CBUF          !  smat output name temporary buffer
    CHARACTER(NAMLEN3) ENAM          !  tmp emission types name
    CHARACTER(NAMLEN3) PNAM          !  input pol name temporary buffer
    CHARACTER(NAMLEN3) SBUF          !  species name temporary buffer
    CHARACTER(NAMLEN3) SNAM          !  speciation pol name temporary buffer


    !***********************************************************************
    !   begin body of program SPCMAT

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !    to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Get environment variables that control program behavior
    !.......  Retrieve the whether to prompt for and use pollutant conversion file
    KFLAG = ENVYN( 'POLLUTANT_CONVERSION',    &
                   'Use pollutant-to-pollutant conversion file',    &
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "POLLUTANT_CONVERSION"', 2 )
    END IF

    !.......  Retrieve the whether to prompt for and use pollutant conversion file
    TFLAG = ENVYN( 'SRC_TAGGING', 'Use source tagging file',    &
                   .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SRC_TAGGING"', 2 )
    END IF

    !.......  Set source category based on environment variable setting
    CALL GETCTGRY

    !.......  Get inventory file names given source category
    CALL GETINAME( CATEGORY, ENAME, ANAME )

    !.......   Get file names and units; open input files

    !.......  Prompt for and open inventory file
    INAME = ENAME
    MESG = 'Enter logical name for the MAP INVENTORY file'
    IDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., INAME, PROGNAME )

    !.......  Open and read map file
    CALL RDINVMAP( INAME, IDEV, ENAME, ANAME, SDEV )

    XDEV = PROMPTFFILE(    &
             'Enter logical name for SPECIATION XREF file',    &
             .TRUE., .TRUE., 'GSREF', PROGNAME )

    RDEV = PROMPTFFILE(    &
             'Enter logical name for SPECIATION PROFILES file',    &
             .TRUE., .TRUE., 'GSPRO', PROGNAME )

    IF( KFLAG )    &
        KDEV = PROMPTFFILE(    &
             'Enter logical name for POLLUTANT CONVERSION file',    &
             .TRUE., .TRUE., 'GSCNV', PROGNAME )

    IF( TFLAG )    &
        GDEV = PROMPTFFILE(    &
             'Enter logical name for TAGGING file',    &
             .TRUE., .TRUE., 'GSTAG', PROGNAME )

    VDEV = PROMPTFFILE(    &
             'Enter logical name for INVENTORY DATA TABLE file',    &
             .TRUE., .TRUE., 'INVTABLE', PROGNAME )

    !.......  Otherwise, store source-category-specific header information,
    !    including the inventory pollutants in the file (if any).
    !    Note that the I/O API head info is passed by include file
    !    and the results are stored in module MODINFO.
    CALL GETSINFO( ENAME )

    !.......   Open files that depend on inventory characteristics
    IF( NIACT .GT. 0 ) THEN
        TDEV = PROMPTFFILE(    &
               'Enter logical name for EMISSION PROCESSES file',    &
               .TRUE., .TRUE., CRL // 'EPROC', PROGNAME )

    END IF

    !.......  Set inventory variables to read for all source categories
    IVARNAMS( 1 ) = 'CSCC'
    IVARNAMS( 2 ) = 'CSOURC'

    !.......  Set inventory variables to read for specific source categories
    IF( CATEGORY .EQ. 'AREA' ) THEN
        NINVARR = 5
        IVARNAMS( 3 ) = 'CMACT'
        IVARNAMS( 4 ) = 'CISIC'
        IVARNAMS( 5 ) = 'CIFIP'

    ELSE IF( CATEGORY .EQ. 'MOBILE' ) THEN
        NINVARR = 6
        IVARNAMS( 3 ) = 'IRCLAS'
        IVARNAMS( 4 ) = 'IVTYPE'
        IVARNAMS( 5 ) = 'CVTYPE'
        IVARNAMS( 6 ) = 'CIFIP'

    ELSE IF( CATEGORY .EQ. 'POINT' ) THEN
        NINVARR = 5
        IVARNAMS( 3 ) = 'CMACT'
        IVARNAMS( 4 ) = 'CISIC'
        IVARNAMS( 5 ) = 'CIFIP'
    END IF

    !.......  Allocate memory for and read required inventory characteristics
    CALL RDINVCHR( CATEGORY, ENAME, SDEV, NSRC, NINVARR, IVARNAMS )

    !.......  Build unique lists of SCCs per SIC from the inventory arrays
    CALL GENUSLST

    !.......  When mobile codes file is being used read mobile codes file
    !    IF( MDEV .GT. 0 ) CALL RDMVINFO( MDEV )

    !.......  Read inventory table (used for NONHAP checks)
    CALL RDCODNAM( VDEV )

    !.......  Perform the steps needed for using activities and emission types
    !    instead of pollutants
    !.......  Read emission processes file.  Populate array in MODEMFAC and
    !    set NETYPE
    IF( NIACT .GT. 0 ) THEN

        CALL RDEPROC( TDEV )

        !.......  Find input hydrocarbon name
        INPUTHC = ' '
        DO I = 1, NEPOL
            SELECT CASE( EMTPOL( I ) )
              CASE( 'VOC', 'THC', 'NMHC', 'TOG', 'NMOG' )
                INPUTHC = EMTPOL( I )
                EXIT
            END SELECT
        END DO

        IF( INPUTHC == ' ' ) THEN
            MESG = 'No valid hydrocarbon pollutant specified in emission processes file'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
        OUTPUTHC = 'NONHAP' // TRIM( INPUTHC )

        FNDOUTPUT = .FALSE.
        K = 0

        !.......  Check if NONHAP values are processed
        DO I = 1, MXIDAT

            IF( INVDNAM( I ) == OUTPUTHC ) THEN
                FNDOUTPUT = .TRUE.
                CYCLE
            END IF

            !.......  If requested hydrocarbon is not TOG or VOC, skip rest of loop
            IF( INPUTHC /= 'TOG' .AND. INPUTHC /= 'VOC' ) EXIT

            IF( INVDVTS( I ) /= 'N' ) THEN

                !.......  Check that pollutant is generated by MOBILE6
                DO J = 1, NEPOL
                    IF( INVDNAM( I ) == EMTPOL( J ) ) THEN
                        IF( INVDVTS( I ) == 'V' ) THEN
                            K = K + 1
                        ELSE IF( INPUTHC == 'TOG' ) THEN
                            K = K + 1
                        END IF
                        EXIT
                    END IF
                END DO
            END IF
        END DO

        !.......  If output was not found, set name to blank
        IF( .NOT. FNDOUTPUT .OR. K == 0 ) THEN
            OUTPUTHC = ' '
        END IF

        !.......  Rename emission factors if necessary
        IF( OUTPUTHC /= ' ' ) THEN
            DO J = 1, NIACT
            DO I = 1, SIZE( EMTNAM,1 )
                IF( EMTNAM( I,J ) == INPUTHC ) THEN
                    EMTNAM( I,J ) = OUTPUTHC
                END IF
            END DO
            END DO
        END IF

    END IF

    !.......  Reset the number of all input variables as the number of pollutants
    !    and emission types, instead of the number of pollutants and
    !    activities
    NIPPA = NIPOL
    DO I = 1, NIACT
        NIPPA = NIPPA + NETYPE( I )
    END DO

    !.......  Set the number of output pollutants based on the number from the
    !    inventory file and the number coming from emission types
    NOPOL = NIPOL + NEPOL

    !.......  Allocate memory for pollutant names, emission types, and associated
    !       pollutants
    DEALLOCATE( EANAM )
    ALLOCATE( EANAM( NIPPA ),    &
              EAIDX( NIPPA ),    &
              SANAM( NIPPA ),    &
              LEMIS( NIPPA ),    &
              IINAM( NOPOL ),    &
              SINAM( NOPOL ),    &
             SPCOUT( NOPOL ),    &
             IDXCHK( NOPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IDXCHK', PROGNAME )

    !.......  Initialize arrays
    EANAM  = ' '          ! array
    EAIDX  = 1            ! array
    SANAM  = ' '          ! array
    LEMIS  = .TRUE.       ! array
    IINAM  = ' '          ! array
    SINAM  = ' '          ! array
    SPCOUT = .TRUE.       ! array
    IDXCHK = .FALSE.      ! array

    !.......  Create array of pollutant names from emission types and pollutants
    !.......  Loop through emission types and add pollutants to master list
    !       in order of appearance
    J  = 0
    NP = 0
    DO I = 1, NIACT
    DO K = 1, NETYPE( I )

        J = J + 1
        ENAM = EMTNAM( K,I )
        EANAM( J ) = ENAM

        PNAM = ENAM

        !.......  If it does not already appear in list, store pollutant name
        N = INDEX1( PNAM, NP, IINAM )

        IF( N .LE. 0 ) THEN
            NP = NP + 1
            IINAM( NP ) = PNAM
            N = NP
        END IF

        !.......  Set index back to master list (IINAM) position
        EAIDX( J ) = N

        !.......  Flag that this entry is activity-based, not raw emissions
        LEMIS( J ) = .FALSE.

    END DO
    END DO

    !.......  Now assign index numbers for pollutants and add pollutants (which
    !    may also be emission types) to EANAM
    K    = 0
    DO I = 1, NIPOL

        J = J + 1

        PNAM = EINAM( I )
        EANAM( J ) = PNAM

        !.......  First extract pollutant name from emission type, if any
        L1 = INDEX( PNAM, ETJOIN )
        L2 = LEN_TRIM( PNAM )
        IF( L1 .GT. 0 ) PNAM = PNAM( L1+2:L2 )

        !.......  Find pollutant name in master list
        N = INDEX1( PNAM, NP, IINAM )

        !.......  If not in list, then add to it
        IF( N .LE. 0 ) THEN
            NP = NP + 1
            IINAM( NP ) = PNAM
            N = NP
        END IF

        !.......  Set index back to master list (IINAM) position
        EAIDX( J ) = N

    END DO

    !.......  Confirm that counter is equal to total number of emission types
    !    and pollutants
    IF( J .NE. NIPPA ) THEN
        WRITE( MESG,94010 ) 'INTERNAL ERROR: number of emission ' //    &
               'types and pollutants ', J, CRLF() // BLANK10 //    &
               'is inconsistent with dimensioned value', NIPPA
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
    END IF

    !.......  Read the speciation cross-reference file
    CALL RDSREF( XDEV )

    !.......  Initialize output emission types/pollutants array
    SANAM = EANAM

    !.......  Read the pollutant to pollutant conversion file, if any
    !.......  resulting tables are passed via MODSPRO
    IF ( KFLAG ) THEN

        CALL RDSCONV( KDEV, NIPPA, EANAM, SANAM )

    END IF

    !.......  Create input and output pollutant names based on output
    !         emission types/pollutants names for input and output.
    J = 0
    K = 0
    LT = LEN_TRIM( ETJOIN )
    DO I = 1, NIPPA

        IDX = EAIDX( I )

        IF( .NOT. IDXCHK( IDX ) ) THEN

            L1 = INDEX( SANAM( I ), ETJOIN )
            L1 = MAX( L1, 1 )
            IF( L1 .GT. 1 ) L1 = L1 + LT

            L2 = LEN_TRIM( SANAM( I ) )
            CBUF = SANAM( I )( L1:L2 )

            !.......  Ensure that pollutant isn't already in the list, which can
            !    happen when the emission type from an activity and
            !    the pollutant are both in the inventory for different
            !    sources.
            IF( K .GT. 0 ) J = INDEX1( CBUF, K, SINAM )
            IF ( J .LE. 0 ) THEN
                K = K + 1
                SINAM( K ) = CBUF
            END IF

            IDXCHK( IDX ) = .TRUE.
        END IF

    END DO

    !.......  Reset number of output pollutants based on new count from removing
    !    the emission processes
    NOPOL = K

    !.......  Scan speciation profiles file to get all of the pollutant-species
    !    combinations that are valid for the pollutants in the inventory.
    !.......  The species names are sorted in ABC order for each pollutant, and
    !    and the pollutants are in the same order as SINAM.
    !.......  Also retrieve the maximum number of species per pollutant and
    !    maximum number of profile entries per pollutant.
    MESG = 'Scanning speciation profiles file for species...'
    CALL M3MSG2( MESG )
    CALL DSCSPROF( RDEV, NOPOL, SINAM )

    !.......  Give warning if some pollutants won't be speciated, and keep track
    !    of which ones don't get species.
    J = 0
    DO V = 1, NOPOL

        IF( SPCNAMES( 1,V ) .EQ. ' ' ) THEN
            J    = J + 1
            MESG = 'WARNING: No speciation profiles found ' //    &
                   'for pollutant "' // TRIM( IINAM( V ) ) // '"' //    &
                   CRLF( )// BLANK10// 'Pollutant ignored!'
            CALL M3MSG2( MESG )

            SPCOUT( V ) = .FALSE.

        END IF

    END DO

    !.......  Make sure at least one pollutant will be speciated
    IF( J .EQ. NIPPA ) THEN
        MESG = 'No speciation profiles for'
        L1 = LEN_TRIM( MESG )

        IF( NEPOL .EQ. 0 ) THEN
            MESG = MESG( 1:L1 ) // ' pollutants in inventory!'
        ELSE IF( NIPOL .EQ. 0 ) THEN
            MESG = MESG( 1:L1 ) // ' emission types!'
        ELSE
            MESG = MESG( 1:L1 ) // ' pollutants or emission types!'
        END IF

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Allocate memory for master 1-d array of species
    NSPCALL = MXSPEC * NOPOL      ! this is reset below
    ALLOCATE ( SPCLIST( NSPCALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPCLIST', PROGNAME )

    K = 0
    DO V = 1, NOPOL
    DO J = 1, MXSPEC
        IF ( SPCNAMES( J,V ) == ' ' ) CYCLE
        K = K + 1
        SPCLIST( K ) = SPCNAMES( J,V )
    END DO
    END DO
    NSPCALL = K

    !.......  Allocate TAGNUM and initialize for all cases (even no tagging)
    !.......  Allocate TAGNAME and initialize for all cases (even no tagging)
    !    Note that the zero position of TAGNAME will always have no tag,
    !    so that the untagged name will be preserved.
    ALLOCATE( TAGNUM(   MXSPEC, NOPOL ),    &
              TAGNAME( 0:MXTAG, MXSPEC, NOPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGNAME', PROGNAME )
    TAGNUM = 0          ! array
    TAGNAME = ' '       ! array

    !.......  Read the tagging file and build the arrays needed for tagging
    IF( TFLAG ) THEN

        CALL RDTAG( GDEV, NOPOL )

    !.......  If any tags assigned, then allocate memory for source tags.
    !    Will be initialized by species in ASGNTAG
        IF( MAXVAL( TAGNUM ) > 0 ) THEN
            ALLOCATE( SRCTAGS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SRCTAGS', PROGNAME )
        END IF

    END IF

    !.......  Determine the number of output species for all pollutants and
    !    emission types. For the case of no emission types, this value
    !    is the same as MXSPEC

    !.......  Allocate memory for speciation factors by source using the maximum
    !    number of species per pollutant, MXSPEC. Also, initialize.
    IF( MASSOUT ) THEN
        ALLOCATE( MASSMATX( NSRC, MXSPEC ),     &
                   MASSARR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MASSARR', PROGNAME )
    END IF
    IF( MOLEOUT ) THEN
        ALLOCATE( MOLEMATX( NSRC, MXSPEC ),     &
                   MOLEARR( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MOLEARR', PROGNAME )
    END IF

    !.......  Allocate memory for arrays of speciation tables and unique lists
    !    using the maximum number of profile table entires per pollutant,
    !    MXSPFUL, which is from module MODSPRO
    !.......  Allocate memory for names of output variables
    ALLOCATE( INPRF( MXSPFUL ),                  &
             SPECID( MXSPFUL ),                  &
           MOLEFACT( MXSPFUL ),                  &
           MASSFACT( MXSPFUL ),                  &
           MASSONAM( 0:MXTAG, MXSPEC, NIPPA  ),  &
           MOLEONAM( 0:MXTAG, MXSPEC, NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INPRF...MOLEONAM', PROGNAME )
    MOLEFACT = 0.0
    MASSFACT = 0.0

    !.......  Open speciation matrix file(s).  Depending on MASSOUT and MOLEOUT,
    !    the mass-based and/or mole-based files will be set up and opened.
    CALL OPENSMAT( ENAME, MASSOUT, MOLEOUT, NOPOL, MXSPEC,     &
                   MXTAG, EANAM, EAIDX, SPCNAMES, MOLUNITS,    &
                   UDEV, SNAME, LNAME, MASSONAM, MOLEONAM )


    !.......  Loop through inventory pollutants to create speciation factors for
    !    each species used for each pollutant. In some cases, a pollutant
    !    will have no species (e.g., CO), so the factors will simply be a
    !    converstion factor from tons to grams or moles.
    DO K = 1, NIPPA

        V = EAIDX( K )
        ENAM = EANAM( K )      ! Emission type name
        PNAM = IINAM( V )      ! Input pollutant name
        SNAM = SINAM( V )          ! Speciation pollutant name

        !.......  Only process for pollutants we know have at least one species,
        !    otherwise, go to end of loop
        IF( SPCOUT( V ) ) THEN

            !.......  Build message stating the pollutant/emission type being
            !    processed.
            !.......  For pollutant data...
            IF( ENAM .EQ. PNAM ) THEN
                L1 = LEN_TRIM( PNAM )
                MESG = 'Processing pollutant "'// PNAM( 1:L1 )// '"'

                IF( PNAM .NE. SNAM ) THEN
                    MESG = TRIM( MESG )// ' using pollutant "' //    &
                           TRIM( SNAM ) // '" for profiles'
                END IF

            !.......  For emission type data...
            ELSE
                MESG = 'Processing emission type "'//    &
                       TRIM( ENAM ) // '" using pollutant "' //    &
                       TRIM( SNAM ) // '" for profiles'
            END IF

            !.......  Write message
            CALL M3MSG2( MESG )

        !.......  Build message stating the pollutant being skipped
        ELSE

            !.......  For pollutant data...
            IF( ENAM .EQ. PNAM ) THEN
                L1 = LEN_TRIM( PNAM )
                MESG = 'Skipping pollutant "' // PNAM( 1:L1 ) // '"'

            !.......  For emission type data...
            ELSE
                L1 = LEN_TRIM( ENAM )
                MESG = 'Skipping emission type "'//ENAM( 1:L1 )//'"'
            END IF

            !.......  Write message stating the pollutant being skipped
            CALL M3MSG2( MESG )
            CYCLE

        END IF

        !.......  Read speciation profiles file
        MESG = BLANK5 // 'Reading speciation profiles file...'
        CALL M3MSG2( MESG )

        CALL RDSPROF( RDEV, SNAM, NMSPC )

        !.......  If current pollutant is a NONHAP* pollutant, compare
        !    its definition from the inventory table to the
        !    definition in the GSPRO file. Note that if the user
        !    provides a different INVTABLE to Spcmat than the one
        !    used to create the inventory, this error will not
        !    be detected.
        CALL CHKNONHAP( PNAM, EFLAG )

        !.......  Initilialize multiple profiles and default reporting to true
        MULTIPRO = .TRUE.
        DEFREPRT = .TRUE.

        !.......  When one profile used for all sources, so do simple processing.
        !    The one exception is when there is a pollutant-to-pollutant
        !    conversion, then we must still do the standard processing.
        IF( NPOLSPRO .EQ. 1 ) THEN

            IF( NMSPC .EQ. 1 ) THEN
                MESG = '     NOTE: "' // TRIM( PNAM ) //     &
                       '" only has a unit conversion ' //    &
                       'using profile "' // INPRF( 1 ) // '"'
            ELSE
                MESG = '     NOTE: "' // TRIM( PNAM ) //    &
                       '" is split for all sources ' //     &
                        'using profile "' // INPRF( 1 ) // '"'
            END IF

            !.......  If there is no pollutant-to-pollutant conversion, then
            !    set speciation matrices using the one profile and continue
            IF( PNAM .EQ. SNAM ) THEN

                IF( MASSOUT ) THEN
                    DO J = 1, NMSPC
                        MASSMATX( :,J ) = MASSFACT( J )
                    END DO
                END IF

                IF( MOLEOUT ) THEN
                    DO J = 1, NMSPC
                        MOLEMATX( :,J ) = MOLEFACT( J )
                    END DO
                END IF

                MULTIPRO = .FALSE.      ! no multiple profiles

            !.......  Otherwise, need to continue so that pollutant-to-pollutant
            !    conversion is done (so don't reset multipro)
            ELSE

                MESG = TRIM( MESG ) // CRLF() // BLANK10 //    &
                       'and a pollutant conversion to "' // TRIM( SNAM ) // '"'

                DEFREPRT = .FALSE.      ! no default reporting

            END IF      ! End EINAM .EQ. SINAM or not

            CALL M3MSG2( MESG )

        END IF          ! End single profile or nots

        !.......  If this pollutant has multiple profiles...
        IF( MULTIPRO ) THEN

            !.......  Abridge profiles so that there is an array of unique profiles
            CALL PROCSPRO( NMSPC, SPCNAMES( 1,V ) )

            !.......  Assign speciation profile and populate speciation matrices
            !    for all sources for this pollutant.
            CALL ASGNSPRO( MASSOUT, MOLEOUT, DEFREPRT, NSRC, UDEV,    &
                           ENAM, LEMIS( K ), MASSMATX, MOLEMATX )

            !.......  Deallocate memory for unique profiles arrays
            DEALLOCATE( SPROFN, IDXSPRO, NSPECIES, IDXSSPEC )

        END IF              ! End multi-profile processing

        !.......  Write out the speciation matrix for current pollutant

        MESG= '     Writing MASS-BASED SPECIATION MATRIX...'
        IF( MASSOUT ) CALL M3MSG2( MESG )

        MESG= '     Writing MOLE-BASED SPECIATION MATRIX...'
        IF( MOLEOUT ) CALL M3MSG2( MESG )

        !.......  Loop through species, apply tags where needed, and write
        !    out matrix arrays for all untagged and tagged species
        DO J = 1, NMSPC

            SBUF = SPCNAMES( J,V )
            NTAG = TAGNUM  ( J,V )

            !.......  If this species has any tags, then assign the tags to
            !    the sources
            IF ( NTAG > 0 ) THEN

                CALL ASGNTAG( SBUF, NSRC, NTAG, TAGNAME( 1,J,V ), SRCTAGS )

            !.......  If no tags matched, then just use matrices as-is
            ELSE IF ( NTAG == 0 ) THEN
                IF( MASSOUT ) MASSARR = MASSMATX(1:NSRC,J)              ! array
                IF( MOLEOUT ) MOLEARR = MOLEMATX(1:NSRC,J)              ! array

            END IF

            !.......  Loop through tags that matches sources for this species
            !.......  Add an extra loop (T=0) for the case of NTAG = 0, so that
            !    the same loop can be used to write the species and to
            !    write the untagged sources for a tagged species.
            DO T = 0, NTAG

                !.......  If tags are actually being used for species J...
                IF( NTAG > 0 ) THEN

                    !.......  Edit species name (note, length has already been
                    !    checked when reading tags file)
                    IF( T > 0 ) SBUF= TRIM( SBUF//TAGNAME( T,J,V ) )

                    !.......  Initialize mass and mole output arrays
                    IF( MASSOUT ) MASSARR = 0.                      ! array
                    IF( MOLEOUT ) MOLEARR = 0.                      ! array

                    !.......  Loop through sources and enter matrices values
                    !    for any sources that match tag index.
                    TAGCNT = 0
                    DO S = 1, NSRC

                        IF ( SRCTAGS( S ) == T ) THEN
                            TAGCNT = TAGCNT + 1
                            IF( MASSOUT ) MASSARR(S) = MASSMATX(S,J)
                            IF( MOLEOUT ) MOLEARR(S) = MOLEMATX(S,J)
                        END IF

                    END DO       ! end loop on sources

                    IF ( TAGCNT == 0 ) THEN
                        MESG = 'WARNING: No sources matched ' //    &
                         'tagging criteria in GSTAG for tag "' //    &
                         TRIM( TAGNAME(T,J,V))// '" and species "'//    &
                         TRIM( SBUF ) // '"'

                        CALL M3MESG( MESG )
                    END IF

                END IF

                !.......  If writing mass-based matrix, then write it
                IF( MASSOUT ) THEN
                    CBUF = MASSONAM( T,J,K )
                    IF( .NOT. WRITESET( SNAME, CBUF, ALLFILES,0,0, MASSARR )) THEN

                        EFLAG = .TRUE.

                        MESG = '     Could not write "' //              &
                          TRIM( ENAM ) // '"-to-"' // TRIM( SBUF ) //   &
                          '" speciation factor using name "' //         &
                          TRIM( CBUF ) // '" to file "' // TRIM( SNAME ) // '"'

                        CALL M3MSG2( MESG )
                        CYCLE

                    END IF

                END IF

                IF( MOLEOUT ) THEN

                    CBUF = MOLEONAM( T,J,K )
                    IF( .NOT.WRITESET( LNAME, CBUF, ALLFILES,0,0, MOLEARR )) THEN

                        EFLAG = .TRUE.

                        MESG = 'Could not write "' // TRIM( ENAM ) //   &
                           '"-to-"' // TRIM( SBUF ) //                  &
                           '" speciation factor using name "' //        &
                           TRIM( CBUF ) // '" to file "' // TRIM( LNAME ) // '"'

                        CALL M3MSG2( MESG )
                        CYCLE

                    END IF    ! if not written out

                END IF        ! End mole-based output

                !.......  If get here, then all writes were successful, so write
                !    message to log file that all is well
                CBUF = MASSONAM( T,J,K )
                MESG = BLANK10 // TRIM( ENAM ) // '-to-' //    &
                        TRIM( SBUF ) // ' written to ' //    &
                        TRIM( SNAME )// ' as variable '// TRIM(CBUF)

            END DO     ! End loop over tags (or single iteration T=0 for no tags)

        END DO         ! End write out loop for model species for current pollutant

    END DO             ! End loop through inventory pollutants

    !.......  Check error flag for problems and end
    IF( EFLAG ) THEN

        MESG = 'Problem running speciation program.See errors above.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Exit program with normal completion
    CALL M3EXIT( PROGNAME, 0, 0, 'Normal completion', 0 )


    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT ( 10 ( A, :, I10, :, 2X ) )

END PROGRAM SPCMAT

