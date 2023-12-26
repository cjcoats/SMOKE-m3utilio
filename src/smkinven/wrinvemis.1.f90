
SUBROUTINE WRINVEMIS( IDEV, DATPATH )

!***********************************************************************
!  subroutine body starts at line 135
!
!  DESCRIPTION:
!      This subroutine writes the average inventory emissions to the inventory
!      files
!
!  PRECONDITIONS REQUIRED:
!      Logical name of output file defined
!      Emission arrays populated
!      MODINFO values assigned
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutine
!
!  REVISION  HISTORY:
!       Created 12/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!*************************************************************************
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

!.........  MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: CSOURC, NPCNT, IPOSCOD, POLVAL

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVDUNT, FIREFLAG, NINVTBL,&
    &                    ITNAMA, ITCASA, FF10FLAG, MEDSFLAG

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, CATDESC, NSRC, NMAP,&
    &                   NIPOL, NIACT, NIPPA, NPPOL, EANAM, ACTVTY,&
    &                   MAPNAM, MAPFIL, EINAM, EIIDX, NPACT, AVIDX,&
    &                   NDY, NC1, NC2, NCE, NRE, NRP, VAR_FORMULA,&
    &                   NCOMP, CHKPLUS, CHKMINUS, FORMULAS, VIN_A,&
    &                   VIN_B, VNAME

!.........  This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!.........  SUBROUTINE ARGUMENTS
    INTEGER           , INTENT (IN) :: IDEV     ! unit number for map file
    CHARACTER(PHYLEN3), INTENT (IN) :: DATPATH  ! path for pol/act files

!...........   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

!.........  Inventory temporary arrays
    INTEGER, ALLOCATABLE:: IPPTR ( : ) ! position in POLVAL sparse array (for input pols/actv)
    INTEGER, ALLOCATABLE:: IPPTR2( : ) ! position in POLVAL sparse array (for computed values)
    INTEGER, ALLOCATABLE:: IPMAX ( : ) ! max IPPTR by source
    INTEGER, ALLOCATABLE:: SRCID ( : ) ! source index
    REAL   , ALLOCATABLE:: SRCPOL( :,: )  ! data-spec values by source

!...........   Names, Units, types, & descriptions for pollutant-specific
!              output variables.  NOTE - second dimension will work so long
!              as NPPOL > NPACT, which is expected to always be the case

    CHARACTER(NAMLEN3), ALLOCATABLE :: EONAMES( :,: ) ! Names
    INTEGER           , ALLOCATABLE :: EOTYPES( :,: ) ! Types (Real|Int)
    CHARACTER(NAMLEN3), ALLOCATABLE :: EOUNITS( :,: ) ! Units
    CHARACTER(MXDLEN3), ALLOCATABLE :: EODESCS( :,: ) ! Dscriptions

!...........   Other local allocatable arrays
    CHARACTER(NAMLEN3), ALLOCATABLE :: SAVEANAM( : ) ! tmp variables

!...........   Other local variables
    INTEGER         F, I, J, S, L, L2, N     ! counters and indices

    INTEGER         IOS       ! i/o status

    INTEGER         MCNT      ! count of actual mapped pol/act files
    INTEGER         MXWARN    ! maximum number of warnings of each type to write

    INTEGER         RIMISS3          ! real value of integer missing

    LOGICAL      :: EFLAG    = .FALSE. ! true: error found
    LOGICAL      :: FFLAG    = .FALSE. ! true: formula in use
    LOGICAL      :: NEGOK    = .FALSE. ! true: okay to output negative emission values
    LOGICAL         ZFLAG              ! true: write zeros to output file

    CHARACTER(80)   NAME1            ! tmp file name component
    CHARACTER(80)   NAME2            ! tmp file name component
    CHARACTER(128)  BUFFER           ! message buffer
    CHARACTER(256)  MESG             ! message buffer

    CHARACTER(PHYLEN3) :: BPATH  ! base path
    CHARACTER(PHYLEN3) :: RPATH  ! relative path

    CHARACTER(16) :: PROGNAME = 'WRINVEMIS' !  program name

!***********************************************************************
!   begin body of program WRINVEMIS

!..........  Get environment variables
    MESG = 'First name of output inventory files'
    CALL ENVSTR( 'INVNAME1', MESG, ' ', NAME1, IOS )
    IF( IOS .NE. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: INVNAME1 environment variable is not' //&
        &       'defined for output file name'
        CALL M3MSG2( MESG )
    END IF

    MESG = 'Second name of output inventory files'
    CALL ENVSTR( 'INVNAME2', MESG, ' ', NAME2, IOS )
    IF( IOS .NE. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: INVNAME2 environment variable is not' //&
        &       'defined for output file name'
        CALL M3MSG2( MESG )
    END IF

    MXWARN  = ENVINT( WARNSET  , ' ', 100, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
    END IF

    NEGOK = ENVYN( 'ALLOW_NEGATIVE',&
    &               'Allow negative output data',&
    &               .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ALLOW_NEGATIVE"', 2 )
    END IF

    IF( EFLAG ) THEN
        MESG = 'Problem with input environment variables'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  If processing fires or FF10 formats, then must write out zero values, since
!           zero annual will be filled with daliy/hourly fire/F10 inventories
    IF ( FIREFLAG .OR. FF10FLAG .OR. MEDSFLAG ) THEN
        ZFLAG = .TRUE.

!.........  Otherwise, check environment to see if user wants to output
!           zero values
    ELSE
        MESG = 'Write zero values to annual emissions inventory'
        ZFLAG = ENVYN( 'WRITE_ANN_ZERO', MESG, .FALSE., IOS )
    END IF

!.........  Compute real value of integer missing
    RIMISS3 = REAL( IMISS3 )

!.........  Set maximum number of map variables for map-formatted outputs
    NMAP = NIPOL + NIACT
    MCNT = 0              ! initialize actual pol/act file count for later

!.........  If there is one or more computed output variable, get set up
    IF( LEN_TRIM( VAR_FORMULA ) .GT. 0 ) THEN
        CALL FORMLIST
        FFLAG = .TRUE.
        NMAP = NMAP + NCOMP  ! NCOMP more variable(s) to map
    END IF

!.........  Allocate memory for map file arrays
    ALLOCATE( MAPNAM( NMAP ),&
    &          MAPFIL( NMAP ),&
    &         EONAMES( NIPOL,NPPOL ),&
    &         EOUNITS( NIPOL,NPPOL ),&
    &         EOTYPES( NIPOL,NPPOL ),&
    &         EODESCS( NIPOL,NPPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EODESCS', PROGNAME )

    MAPNAM = ' '
    MAPFIL = ' '

!.........  Get names, units, etc. of output pollutant-specific records
    CALL BLDENAMS( CATEGORY, NIPOL, NPPOL, EINAM,&
    &               EONAMES, EOUNITS, EOTYPES, EODESCS )

!.........  Set up for opening I/O API sparse pollutant output files
    CALL HDRMISS3  ! Initialize for emissions

!.........  Set number of variables and allocate file description arrays
    NVARSET = 1 + NPPOL   ! Additional 1 for the SRCID variable
    WRITE( FDESC3D( 1 ), '(A,1X,I8)' ) '/NSRC/', NSRC

    IF( ALLOCATED( VTYPESET ) )&
    &    DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )
    ALLOCATE( VTYPESET( NVARSET ),&
    &          VNAMESET( NVARSET ),&
    &          VUNITSET( NVARSET ),&
    &          VDESCSET( NVARSET ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VNAMESET...VDESCSET', PROGNAME )

    IF( ALLOCATED( VARS_PER_FILE ) ) DEALLOCATE( VARS_PER_FILE )

    VTYPESET = 0    ! array initialization
    VNAMESET = ' '  ! array initialization
    VUNITSET = ' '  ! array initialization
    VDESCSET = ' '  ! array initialization

!.........  Set up source ID variable information, which is the same for
!           all pollutant files
    VTYPESET( 1 ) = M3INT
    VNAMESET( 1 ) = 'SRCID'
    VUNITSET( 1 ) = 'n/a'
    VDESCSET( 1 ) = 'Source ID number'

!.........  Separate the pol/act file path into two parts to be
!           able to build relative file names for the map file.
    L = LEN_TRIM( DATPATH )
    DO N = L, 1, -1

        IF( DATPATH( N:N ) .EQ. '/' .OR.&
        &    DATPATH( N:N ) .EQ. '\'      ) THEN
            BPATH = DATPATH( 1:N )
            RPATH = DATPATH( N+1:L )
            EXIT
        END IF

    END DO

!.........  Allocate memory for indices IPPTR & IPMAX for pointing to position
!           in sparsely stored data array.
!.........  Allocate memory for storing and writing emissions
    ALLOCATE( IPPTR( NSRC ),&
    &         IPPTR2( NSRC ),&
    &          IPMAX( NSRC ),&
    &         SRCPOL( NSRC, NPPOL ),&
    &          SRCID( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SRCID', PROGNAME )

!.........  Loop through non-formula pollutants, store, and write to inventory file
    IF( NIPOL .GT. 0 )&
    &    CALL LOOP_FOR_OUTPUT( NIPOL, NPPOL, IDEV, EIIDX,&
    &                          EINAM, MCNT )

!.........  Loop through formula-based pollutants, compute, and write to inventory file
    IF( NCOMP .GT. 0 )&
    &    CALL OUTPUT_FORMULAS( NIPOL, NPPOL, IDEV, EIIDX,&
    &                          EINAM, MCNT )

!.........  If computed variable, update EANAM, in case needed for day-
!           and hour-specific data
    IF( FFLAG ) THEN

        ALLOCATE( SAVEANAM( NIPPA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SAVEANAM', PROGNAME )
        SAVEANAM = EANAM       ! array

        DEALLOCATE( EANAM )
        ALLOCATE( EANAM( NIPPA+NCOMP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EANAM', PROGNAME )

!.............  Save names of existing pollutants
        EANAM( 1:NIPOL ) = SAVEANAM( 1:NIPOL )

!.............  Add names of computed pollutants
        DO F = 1, NCOMP
            EANAM( NIPOL+F ) = VNAME( F )
        END DO

!.............  Save names of existing activities
        EANAM( NIPOL+NCOMP+1:NIPPA+NCOMP )=SAVEANAM( NIPOL+1:NIPPA )

!.............  Update pollutant/activity counts
        NIPOL = NIPOL + NCOMP
        NIPPA = NIPPA + NCOMP

!.............  Release memory for temporary array
        DEALLOCATE( SAVEANAM )

    END IF

!.........  Output for activities...

!.........  Deallocate arrays for variable names
    DEALLOCATE( EONAMES, EOUNITS, EOTYPES, EODESCS )

!.........  Allocate memory for temporary variable names etc.
    ALLOCATE( EONAMES( NIACT,NPACT ),&
    &          EOUNITS( NIACT,NPACT ),&
    &          EOTYPES( NIACT,NPACT ),&
    &          EODESCS( NIACT,NPACT ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EONAMES...EODESCS', PROGNAME )

!.........  Get names, units, etc. of output activity-specific records
    CALL BLDENAMS( CATEGORY, NIACT, NPACT, ACTVTY,&
    &               EONAMES, EOUNITS, EOTYPES, EODESCS )

!.........  Set up for opening I/O API sparse pollutant output files
    CALL HDRMISS3  ! Initialize for emissions

!.........  Set number of variables and allocate file description arrays
    NVARSET = 1 + NPACT
    IF ( CATEGORY .EQ. 'POINT' )  NVARSET = NVARSET + 4     !!  FUG_* variables...
    WRITE( FDESC3D( 1 ), '(A,1X,I8)' ) '/NSRC/', NSRC

    IF( ALLOCATED( VARS_PER_FILE ) ) DEALLOCATE( VARS_PER_FILE )

!.........  Loop through activity data, store, and write to inventory file
    IF( NIACT .GT. 0 )&
    &    CALL LOOP_FOR_OUTPUT( NIACT, NPACT, IDEV, AVIDX,&
    &                          ACTVTY, MCNT )
!.........  Deallocate local arrays
    IF( ALLOCATED( IPPTR ) )    DEALLOCATE( IPPTR )
    IF( ALLOCATED( IPPTR2 ) )   DEALLOCATE( IPPTR2 )
    IF( ALLOCATED( IPMAX ) )    DEALLOCATE( IPMAX )
    IF( ALLOCATED( SRCPOL ) )   DEALLOCATE( SRCPOL )
    IF( ALLOCATED( SRCID ) )    DEALLOCATE( SRCID )

    DEALLOCATE( EONAMES, EOTYPES, EOUNITS, EODESCS )

!.........  Reset the number of map records, in case some had no data
    NMAP = MCNT

!.........  Write the map inventory file
    WRITE( IDEV, '(A)' ) CATDESC
    WRITE( IDEV, '(A)' ) '/IOAPI/ ' // TRIM( NAME1 ) // '.ncf'
    WRITE( IDEV, '(A)' ) '/TEXT/ ' // TRIM( NAME2 ) // '.txt'
    WRITE( IDEV, '(A,I8)' ) '/NDAT/ ', NMAP
    WRITE( IDEV, '(A)' ) '/DATMAP/'

    DO I = 1, NMAP
!.............  Write the relative file name
        WRITE( IDEV, '(A)' ) MAPNAM( I ) // ' ' // TRIM( MAPFIL(I) )

!.............  Reset with the full path, in case needed for hour-specific
!               or day-specific processing emissions reading
        MAPFIL( I ) = TRIM( BPATH ) // TRIM( MAPFIL(I) )

    END DO

    WRITE( IDEV, '(A)' ) '/END/'

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram is for writing out the inventory
!               data, whether it is the pollutant data or the activity data.
!.............  Most variables are defined from the main subroutine
    SUBROUTINE LOOP_FOR_OUTPUT( NOUT, NPVAR, IDEV, INDX,&
    &                            NAMES, MCNT )

!.............  Subroutine arguments
        INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
        INTEGER     , INTENT (IN) :: NPVAR         ! no. vars per data
        INTEGER     , INTENT (IN) :: IDEV          ! unit no. for map file
        INTEGER     , INTENT (IN) :: INDX ( NOUT ) ! index to master list
        CHARACTER(*), INTENT (IN) :: NAMES( NOUT ) ! names of pols/act
        INTEGER , INTENT (IN OUT) :: MCNT          ! map counter

!.............  Local variables
        INTEGER    I, J, K, L, M

        INTEGER          NREC      ! number of output records for each pollutant

        CHARACTER(PHYLEN3) :: RFNAME ! relative physical file name

!----------------------------------------------------------------------

!.............  Loop through output variables (pollutants or activities)
        DO I = 1, NOUT

            NREC = 0
            CALL FILL_OUTPUT_DATA( INDX(I), NPVAR, 'ORDERED', NREC )

!.................  Give warning if no emissions data available
            IF( NREC .EQ. 0 ) THEN
                MESG = 'WARNING: no sources with data for output '//&
                &       'of data variable "' // TRIM( NAMES(I) )//'"'&
                &       //CRLF()//BLANK10//&
                &       'File will not be written.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

!.................  Set up to open output file...
!.................  Set the pollutant-specific items for the I/O API header
            NROWS3D = NREC

            M = 1
            DO J = 1, NPVAR
                M = M + 1
                VNAMESET( M ) = EONAMES( I,J )
                VTYPESET( M ) = EOTYPES( I,J )
                VUNITSET( M ) = EOUNITS( I,J )
                VDESCSET( M ) = EODESCS( I,J )
            END DO

!.................  Reset units for the primary (annual) data value
            K = INDEX1( EONAMES( I,1 ), MXIDAT, INVDNAM )
            IF( K .LE. 0 ) THEN
                MESG = 'ERROR: Cannot find pollutant "'//&
                &       TRIM( EONAMES( I,1 ) ) // '" from '//&
                &       'inventory in INVTABLE file.'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

            VUNITSET( 2 ) = INVDUNT( K )

!.................  Build name and output pollutant or activity file
            RFNAME = TRIM( RPATH )// '/'// TRIM( NAMES(I) )// '.ncf'

            CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NSRC, NPVAR,&
            &               NAMES(I), SRCID, SRCPOL, EFLAG )

!...............  Update map file counter and map file names
            MCNT = MCNT + 1
            MAPNAM( MCNT ) = NAMES( I )
            MAPFIL( MCNT ) = TRIM( RFNAME )

        END DO  ! end loop I through output variables

!.............  Exit program
        IF( EFLAG ) THEN
            MESG = 'Problem writing data.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

!----------------------  FORMAT  STATEMENTS   --------------------------

!...........   Internal buffering formats............ 94xxx

94020   FORMAT( 10( A, :, E10.2, :, 1X ) )

    END SUBROUTINE LOOP_FOR_OUTPUT

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!.............  This internal subprogram is for writing out the computed formulas'
!               inventory data.
!.............  Most variables are defined from the main subroutine
    SUBROUTINE OUTPUT_FORMULAS( NOUT, NPVAR, IDEV, INDX,&
    &                            NAMES, MCNT )

!.............  Subroutine arguments
        INTEGER     , INTENT (IN) :: NOUT          ! no. pols/act for output
        INTEGER     , INTENT (IN) :: NPVAR         ! no. vars per data
        INTEGER     , INTENT (IN) :: IDEV          ! unit no. for map file
        INTEGER     , INTENT (IN) :: INDX ( NOUT ) ! index to master list
        CHARACTER(*), INTENT (IN) :: NAMES( NOUT ) ! names of pols/act
        INTEGER , INTENT (IN OUT) :: MCNT          ! map counter

!.............  Arrays for emissions or activities output
        REAL            COMPUTED( NSRC, NPTPPOL3 ) ! computed data-spec values by src

!.............  Local variables
        INTEGER    F, J, K, L, M, N, S

        INTEGER          IDXA      ! position of first formula variable in list
        INTEGER          IDXB      ! position of 2nd formula variable in list
        INTEGER          NREC      ! number of output records for each pollutant
        INTEGER          WARNCNT_A ! number of times output warning A
        INTEGER          WARNCNT_B ! number of times output warning B

        REAL           MINANN      ! min negative annual value
        REAL           MINAVD      ! min negative seasonal/ave day value

        LOGICAL     :: FFLAGA = .FALSE.  ! true: first part of formula processed
        LOGICAL     :: FFLAGB = .FALSE.  ! true: 2nd part of formula processed

        CHARACTER*256      :: BUFFER
        CHARACTER(PHYLEN3) :: RFNAME ! relative physical file name


!----------------------------------------------------------------------

!.............  Write message saying that formula variables are being computed
        WRITE( MESG,94010 ) 'NOTE: Computing variables for',&
        &                     NCOMP, 'formulas.'
        CALL M3MSG2( MESG )

!.............  Deallocate arrays for variable names for output to SMOKE intermediate files
        DEALLOCATE( EONAMES, EOUNITS, EOTYPES, EODESCS )

!.............  Allocate memory for temporary variable names etc.
        ALLOCATE( EONAMES( NCOMP,NPVAR ),&
        &          EOUNITS( NCOMP,NPVAR ),&
        &          EOTYPES( NCOMP,NPVAR ),&
        &          EODESCS( NCOMP,NPVAR ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EEONAMES...ODESCS', PROGNAME )

!.............  Get names, units, etc. of output activity-specific records
        CALL BLDENAMS( CATEGORY, NCOMP, NPVAR, VNAME,&
        &               EONAMES, EOUNITS, EOTYPES, EODESCS )

!.............  Loop through formulas. Note that the variables have already
!               been checked in the main routine to make sure the names match if
!               the code has gotten this far, so no additional checks on the INDEX1
!               calls included here.
        DO F = 1, NCOMP

!.................  Write message about current formula
            MESG = 'NOTE: Computing variable "'//TRIM( VNAME(F) )//&
            &       '" using formula "'// TRIM( FORMULAS(F) )// '"'
            CALL M3MSG2( MESG )

!.................  Initialize values for minimum negative values
            MINANN = 0.
            MINAVD = 0.

!.................  Initialize values for warning counts (once per formula)
            WARNCNT_A = 0
            WARNCNT_B = 0

!.................  Initialize array to account for sparse storage when filling in
            COMPUTED = 0.  ! array

!.................  Find first variable in formula in list of pollutants
            IDXA = INDEX1( VIN_A( F ), NOUT, NAMES )

!.................  Populate SRCPOL array from first input variable
            NREC = 0
            CALL FILL_OUTPUT_DATA( INDX(IDXA),NPVAR,'SEARCH',NREC )

!.................  Give error if no emissions data available
            IF( NREC .EQ. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: no sources with data for output '//&
                &       'of data variable "'// TRIM(NAMES(IDXA))//'"'&
                &       //CRLF()//BLANK10//&
                &       'Variable "'//TRIM(VNAME(F))//'" will '//&
                &       'not be computed.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

!.................  Initialize full-source computational array with first variable
            DO N = 1, NREC
                S = SRCID( N )
                COMPUTED( S,1 ) = MAX( 0., SRCPOL( N,1 ) )

                IF( NPVAR .GT. 1 ) COMPUTED( S,2 ) =&
                &                   MAX( 0., SRCPOL( N,2 ) ) ! Average day value

                IF( NPVAR .GT. 2 )                                &! Other info
                &    COMPUTED( S,3:NPVAR ) = SRCPOL( N,3:NPVAR )

            END DO

!.................  Find second variable in formula in list of pollutants
            IDXB = INDEX1( VIN_B( F ), NOUT, NAMES )

!.................  Populate SRCPOL array with second input variable
            NREC = 0
            CALL FILL_OUTPUT_DATA( INDX(IDXB),NPVAR,'SEARCH',NREC )

!.................  Give error if no emissions data available
            IF( NREC .EQ. 0 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: no sources with data for output '//&
                &       'of data variable "'// TRIM(NAMES(IDXB))//'"'&
                &       //CRLF()//BLANK10//&
                &       'Variable "'//TRIM(VNAME(F))//'" will '//&
                &       'not be computed.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF

!.................  Compute variable by adding or subtracting second variable
            DO N = 1, NREC
                S = SRCID( N )

                IF( CHKPLUS( F ) ) THEN     ! Formula is an addition

                    COMPUTED( S,1 ) = COMPUTED( S,1 ) +&
                    &                  MAX( 0., SRCPOL( N,1 ) )

                    IF ( NPVAR .GT. 1 )&
                    &     COMPUTED( S,2 ) =&
                    &     COMPUTED( S,2 ) + MAX( 0.,SRCPOL( N,2 ) )

                ELSE IF( CHKMINUS( F ) ) THEN    ! Formula is a subtraction

                    COMPUTED( S,1 ) = COMPUTED( S,1 ) -&
                    &                  MAX( 0., SRCPOL( N,1 ) )

                    IF ( NPVAR .GT. 1 )&
                    &     COMPUTED( S,2 ) =&
                    &     COMPUTED( S,2 ) - MAX( 0.,SRCPOL( N,2 ) )

!.........................  Check to see if the computed value is now negative, and
!                           if so, reset to zero.
                    IF( COMPUTED( S,1 ) .LT. 0 ) THEN
                        WARNCNT_A = WARNCNT_A + 1

!.............................  If warning count is less than max and not
!                               fires.
                        IF ( WARNCNT_A .LE. MXWARN .AND.&
                        &     .NOT. FIREFLAG ) THEN

                            CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )

                            IF( NEGOK ) THEN

                                WRITE( MESG,94020 ) 'WARNING: '//&
                                &'Retaining negative value of annual "'//&
                                &TRIM(VNAME(F))//'" of', COMPUTED(S,1),&
                                &' for source:'//CRLF()//BLANK10//&
                                &BUFFER( 1:L2 )

                            ELSE

                                WRITE( MESG,94020 ) 'WARNING: '//&
                                &'Resetting negative value of annual "'//&
                                &TRIM(VNAME(F))//'" from', COMPUTED(S,1),&
                                &'to 0. for source:'//CRLF()//BLANK10//&
                                &BUFFER( 1:L2 )

                            END IF

                            CALL M3MESG( MESG )
                        END IF

                        MINANN = MIN( MINANN, COMPUTED( S,1 ) )
                        IF( .NOT. NEGOK ) COMPUTED( S,1 ) = 0.

                    END IF

!..........................  Check for negative values for average day value
                    IF( COMPUTED( S,2 ) .LT. 0 ) THEN
                        WARNCNT_B = WARNCNT_B + 1

                        IF ( WARNCNT_B .LE. MXWARN ) THEN
                            CALL FMTCSRC( CSOURC( S ), 7, BUFFER, L2 )

                            IF( NEGOK ) THEN
                                WRITE( MESG,94020 ) 'WARNING: '//&
                                & 'Retaining negative value of "'//&
                                & 'average-day "'//TRIM(VNAME(F))//&
                                & '" of',COMPUTED(S,2),'for source:'//&
                                & CRLF()//BLANK10//BUFFER(1:L2)
                            ELSE
                                WRITE( MESG,94020 ) 'WARNING: '//&
                                & 'Resetting negative value of "'//&
                                & 'average-day "'//TRIM(VNAME(F))//&
                                & '" from',COMPUTED(S,2),'to 0. for '//&
                                &'source:'//CRLF()//BLANK10//BUFFER(1:L2)
                            END IF
                            CALL M3MESG( MESG )
                        END IF

                        MINAVD = MIN( MINAVD, COMPUTED( S,2 ) )
                        IF ( .NOT. NEGOK ) COMPUTED( S,2 ) = 0.

                    END IF

                END IF   ! If formula is a plus or a minus

            END DO       ! End loop over values in second variable

!................   Copy data to output structure (no condensing output structure)
            IF ( ZFLAG ) THEN
                DO S = 1, NSRC
                    SRCID( S ) = S
                END DO

                SRCPOL( 1:NSRC,1:NPVAR ) = COMPUTED( 1:NSRC,1:NPVAR)   ! array
                NREC = NSRC

!................   Condense data to output structure by excluding zeros
            ELSE
                N = 0
                DO S = 1, NSRC

!.........................  Skip records with zero data, unless option (ZFLAG) says not to
                    IF( COMPUTED( S,1 ) == 0 .AND.&
                    &    COMPUTED( S,2 ) == 0       ) CYCLE

                    N = N + 1
                    SRCID ( N ) = S
                    SRCPOL( N,1:NPVAR ) = COMPUTED( S,1:NPVAR )

                END DO          ! End loop through sources for sparse storage
                NREC = N

            END IF

!.................  Give an error when all computed values are zero
            IF( NREC == 0 ) THEN
                MESG = 'WARNING: All computed values are zero. Nothing to output'
                CALL M3MESG( MESG )
                CYCLE
            END IF

!...............  Give warning for that includes negative annual value
            IF( MINANN .LT. 0. ) THEN
                WRITE( MESG,94020 ) 'WARNING: Largest negative '//&
                &       'value for "'// TRIM( EONAMES(F,1) ) //&
                &       '" was', MINANN
                IF( .NOT. NEGOK )&
                &    MESG = TRIM( MESG ) // ' and was reset to 0.'

                CALL M3MSG2( MESG )
            END IF

!...............  Give warning for that includes negative seasonal value
            IF( MINAVD .LT. 0. ) THEN
                WRITE( MESG,94020 ) 'WARNING: Largest negative '//&
                &       'value for "'// TRIM( EONAMES(F,2) ) //&
                &       '" was', MINAVD
                IF( .NOT. NEGOK )&
                &    MESG = TRIM( MESG ) // ' and was reset to 0.'

                CALL M3MSG2( MESG )
            END IF

!.................  Set up to open output file...
!.................  Set the pollutant-specific items for the I/O API header
            NROWS3D = NREC

            M = 1
            DO J = 1, NPVAR
                M = M + 1
                VNAMESET( M ) = EONAMES( F,J )
                VTYPESET( M ) = EOTYPES( F,J )
                VUNITSET( M ) = EOUNITS( F,J )
                VDESCSET( M ) = EODESCS( F,J )
            END DO

!.................  Reset units for the primary (annual) data value
            K = INDEX1( EONAMES( F,1 ), MXIDAT, INVDNAM )
            IF( K .LE. 0 ) THEN
                MESG = 'ERROR: Cannot find pollutant "'//&
                &       TRIM( EONAMES( F,1 ) ) // '" from '//&
                &       'inventory in INVTABLE file.'
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
                CYCLE
            END IF

            VUNITSET( 2 ) = INVDUNT( K )

!...............  Build name and output pollutant file
            RFNAME = TRIM( RPATH )// '/'// TRIM( VNAME(F) )// '.ncf'
            CALL WRINVPOL( CATEGORY, BPATH, RFNAME, NSRC, NPVAR,&
            &               VNAME(F), SRCID, SRCPOL, EFLAG )

!...............  Update map file counter and map relative file names
            MCNT = MCNT + 1
            MAPNAM( MCNT ) = VNAME( F )
            MAPFIL( MCNT ) = TRIM( RFNAME )

        END DO   !  End of formulas, loop on F

!.............  Exit if error found (all formulas couldn't be computed)
        IF ( EFLAG ) THEN
            MESG = 'Could not output all formulas. ' //&
            &       'See previous error messages.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        RETURN

!......................  FORMAT  STATEMENTS   ..........................

!...............   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

94020   FORMAT( 10( A, :, E10.2, :, 1X ) )

    END SUBROUTINE OUTPUT_FORMULAS

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!.............  This internal subprogram is for writing out the computed formulas'
!               inventory data.
!.............  Most variables are defined from the main subroutine
    SUBROUTINE FILL_OUTPUT_DATA( POLINDEX, NPVAR,&
    &                             CALLTYPE, NREC   )

!.............  Subroutine arguments
        INTEGER     , INTENT(IN) :: POLINDEX ! index of pollutant or activity to fill
        INTEGER     , INTENT(IN) :: NPVAR    ! no. variables per pol/act
        CHARACTER(*), INTENT(IN) :: CALLTYPE ! ordered or search
        INTEGER     , INTENT(OUT):: NREC     ! number of non-zero records in array for current pollutant

!.............  Local variables
        INTEGER    I, J, K, N, S

        LOGICAL       :: FFLAG             ! pollutant has been found
        LOGICAL, SAVE :: FIRSTIME = .TRUE. ! first time called

        CHARACTER*256    MESG  ! error message buffer

!----------------------------------------------------------------------

!.............  Since the emission values are already sorted in the output order of
!               the pollutants/activities, can use the IPPTR and IPMAX indices
!               to keep track of what pollutant/activity we are on for each source
!               (sources can have different numbers and types of output pollutants
!               and activities).
!.............  IPPTR2 has been added for use in "search" mode.  It is never
!               reset from its initial value, unlike IPPTR.  Need to have both
!               because IPPTR is used for the input pollutants and activities,
!               whereas IPPTR2 is used for just the computed pollutants. Since
!               the computed pollutants are inserted between the input pollutants
!               and the activities, it's necessary to have both arrays
        IF ( FIRSTIME ) THEN

            IPPTR ( 1 ) = 1
            IPPTR2( 1 ) = 1
            IPMAX( 1 ) = NPCNT( 1 )
            DO S = 2, NSRC
                IPPTR ( S ) = IPPTR( S-1 ) + NPCNT( S-1 )
                IPPTR2( S ) = IPPTR( S )
                IPMAX ( S ) = IPPTR( S )   + NPCNT( S ) - 1
            END DO
            FIRSTIME = .FALSE.

        END IF

!.............  Initialize SRCPOL pollutant/activity-specific data array
!               as missing
        SRCPOL( :,1 ) = 0.      ! array
        IF( NDY .GT. 0 ) SRCPOL( :,NDY ) = 0.      ! array
        IF( NC1 .GT. 0 ) SRCPOL( :,NC1 ) = RIMISS3 ! array
        IF( NC2 .GT. 0 ) SRCPOL( :,NC2 ) = RIMISS3 ! array
        IF( NCE .GT. 0 ) SRCPOL( :,NCE ) = 0.      ! array
        IF( NRE .GT. 0 ) SRCPOL( :,NRE ) = 100.    ! array
        IF( NRP .GT. 0 ) SRCPOL( :,NRP ) = 100.    ! array

!.................  Transfer emissions or activity data to output SRCPOL array
        K = 0
        N = 0
        DO S = 1, NSRC

            FFLAG = .FALSE.
            SELECT CASE( CALLTYPE )
              CASE ( 'ORDERED' )

!.....................  Set position in sparse array by comparing pointer to max
                K = MIN( IPPTR( S ), IPMAX( S ) )
                IF( IPOSCOD( K ) .EQ. POLINDEX ) THEN
                    IPPTR( S ) = K + 1  ! pointer for source S
                    FFLAG = .TRUE.
                END IF

              CASE ( 'SEARCH' )

                K = MIN( IPPTR2( S ), IPMAX( S ) )
                DO J = 1, NPCNT( S )

                    IF( IPOSCOD( K ) .EQ. POLINDEX ) THEN
                        FFLAG = .TRUE.
                        EXIT            ! pollutant has been found, break out of loop
                    ELSE
                        K = MIN( K + 1, IPMAX( S ) )

                    END IF

                END DO

            END SELECT

!.....................  Retrieve emissions from the pollutant array if the
!                       source pollutant ID equals the pollutant ID of the
!                       pollutant of iteration I.
            IF( FFLAG ) THEN

!......................  Skip records with zero data or missing data
                IF( .NOT. ZFLAG ) THEN
                    IF( POLVAL( K,1 ) .LE. 0 . AND.&
                    &    POLVAL( K,2 ) .LE. 0        ) CYCLE
                END IF

                N = N + 1
                SRCID( N ) = S

                DO J = 1, NPVAR     ! rearrange pollutant-specific info
                    SRCPOL( N,J ) = POLVAL( K,J )
                END DO

            END IF           ! If pol/act available for current source

        END DO  ! end of loop through sources
        NREC = N

        RETURN

    END SUBROUTINE FILL_OUTPUT_DATA

END SUBROUTINE WRINVEMIS
