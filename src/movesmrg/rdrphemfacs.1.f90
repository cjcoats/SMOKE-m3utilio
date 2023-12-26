
SUBROUTINE RDRPHEMFACS( REFIDX, MONTH )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!       Reads the emission factor for the given county and month
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:  none
!
!  REVISION  HISTORY:
!       06/14: Created by B.H. Baek
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
!***********************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module is used for reference county information
    USE MODMBSET, ONLY: MCREFIDX

!.........  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: MRCLIST, MVFILDIR, EMPOLIDX,&
    &                     NEMTEMPS, EMTEMPS, RPHEMFACS,&
    &                     NMVSPOLS, MVSPOLNAMS, NOXADJFLAG

!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: NSMATV, TSVDESC, NMSPC, NIPPA, EMNAM, EANAM

!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVSCC, INVSCC

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'MVSCNST3.EXT'  !  MOVES contants

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT(IN) :: REFIDX       ! ref. county index
    INTEGER, INTENT(IN) :: MONTH        ! current processing month

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL

!...........   Local parameters
    INTEGER, PARAMETER :: NSEG = 200    ! number of segments
    INTEGER, PARAMETER :: NNONPOL = 6   ! number of non-pollutant fields in file

!...........   Local allocatable arrays
    CHARACTER(50),  ALLOCATABLE :: SEGMENT( : )    ! parsed input line
    LOGICAL,        ALLOCATABLE :: LINVSCC( : )    ! check inv SCC availability in lookup table

!...........   Other local variables
    INTEGER     I, J, L, LJ, L1, N, P, V  ! counters and indexes
    INTEGER     IOS         ! error status
    INTEGER  :: IREC = 0    ! record counter
    INTEGER  :: TDEV = 0    ! tmp. file unit
    INTEGER     SCCIDX
    INTEGER     TMPIDX
    INTEGER     NSCC        ! no of processing SCCs
    INTEGER  :: TOGIDX = 0  ! index of TOG pollutant
    INTEGER  :: NHTOGIDX = 0  ! index of NONHAPTOG pollutant

    REAL        TMPVAL      ! temperature value
    REAL        PTMP        ! previous temperature value
    REAL        EMVAL       ! emission factor value

    LOGICAL     NOX_ADJUSTED  ! true? header for NOx adj was found
    LOGICAL     FOUND       ! true: header record was found
    LOGICAL     SKIPSCC     ! true: current SCC is not in inventory
    LOGICAL  :: EFLAG = .FALSE.

    CHARACTER(PLSLEN3)  SVBUF     ! tmp speciation name buffer
    CHARACTER(NAMLEN3)  CPOL      ! tmp pollutant buffer
    CHARACTER(NAMLEN3)  CSPC      ! tmp species buffer
    CHARACTER(SCCLEN3)  TSCC      ! current SCC
    CHARACTER(SCCLEN3)  PSCC      ! previous SCC
    CHARACTER(FIPLEN3)  CFIP      ! current ref county

    CHARACTER(10000)    LINE          ! line buffer
    CHARACTER(100)      FILENAME      ! tmp. filename
    CHARACTER(200)      FULLFILE      ! tmp. filename with path
    CHARACTER(300)      MESG          ! message buffer

    CHARACTER(16) :: PROGNAME = 'RDRPHEMFACS'    ! program name

!***********************************************************************
!   begin body of subroutine RDRPHEMFACS

!.........  Open emission factors file based on MRCLIST file
    FILENAME = TRIM( MRCLIST( REFIDX, MONTH ) )

    IF( FILENAME .EQ. ' ' ) THEN
        WRITE( MESG, 94010 ) 'ERROR: No emission factors file ' //&
        &  'for reference county ' // MCREFIDX( REFIDX,1 ) //&
        &  ' and fuel month', MONTH
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    FULLFILE = TRIM( MVFILDIR ) // FILENAME
    OPEN( TDEV, FILE=FULLFILE, STATUS='OLD', IOSTAT=IOS )
    IF( IOS .NE. 0 ) THEN
        MESG = 'ERROR: Could not open emission factors file ' //&
        &  FILENAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        MESG = 'Reading emission factors file ' //&
        &  FILENAME
        CALL M3MESG( MESG )
    END IF

!.........  Allocate memory to parse lines
    ALLOCATE( SEGMENT( NSEG ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

!.........  Read header line to get list of pollutants in file
    NOX_ADJUSTED = .TRUE.
    FOUND = .FALSE.
    IREC = 0
    NEMTEMPS = 0
    DO

        READ( TDEV, 93000, END=100, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 ) 'I/O error', IOS,&
            &  'reading emission factors file ' //&
            &  FILENAME // ' at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Check for header line
        IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) THEN
            LJ = LEN_TRIM( LINE )
            NEMTEMPS = STR2INT( LINE( 13:LJ ) )

        ELSE IF( LINE( 1:21 ) .EQ. 'HUMIDITY_ADJUSTED_NOX' ) THEN
            LJ = LEN_TRIM( LINE )
            IF( TRIM( LINE( 23:LJ ) ) == 'N' ) NOX_ADJUSTED = .FALSE.  ! Do not apply NOx humidity adj

        ELSE IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) THEN
            FOUND = .TRUE.

            SEGMENT = ' '  ! array
            CALL PARSLINE( LINE, NSEG, SEGMENT )

!.................  Count number of pollutants
            NMVSPOLS  = 0
            DO J = NNONPOL + 1, NSEG

                IF( SEGMENT( J ) .NE. ' ' ) THEN
                    NMVSPOLS = NMVSPOLS + 1
                    IF( SEGMENT( J ) == 'NONHAPTOG' ) THEN
                        NHTOGIDX = INDEX1( 'NONHAPTOG', NIPPA, EANAM )
                    ELSE IF( SEGMENT( J ) == 'TOG' ) THEN
                        TOGIDX = INDEX1( 'TOG', NIPPA, EANAM )
                    END IF
                ELSE
                    EXIT
                END IF

            END DO

            IF( ALLOCATED( MVSPOLNAMS ) ) DEALLOCATE( MVSPOLNAMS )
            ALLOCATE( MVSPOLNAMS( NMVSPOLS  ), STAT=IOS )
            CALL CHECKMEM( IOS, 'MVSPOLNAMS', PROGNAME )
            MVSPOLNAMS = ''

!.................  Store pollutant names
            DO J = 1, NMVSPOLS
                MVSPOLNAMS( J ) = SEGMENT( NNONPOL + J )
            END DO

            EXIT

        END IF

    END DO

!.........  Do not apply NOx humidity adjustment since it has been already applied
    IF( NOX_ADJUSTED .AND. NOXADJFLAG ) THEN
        MESG = 'WARNING: OVERRIDE the setting of APPLY_NOX_HUMIDITY_ADJ to N '//&
        &       'since NOx adustment has been already made by MOVES'
        CALL M3MESG( MESG )
        NOXADJFLAG = .FALSE.
    END IF

100 CONTINUE

    REWIND( TDEV )

    IF( .NOT. FOUND ) THEN
        MESG = 'ERROR: Missing header line in emission factors file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( NMVSPOLS  == 0 ) THEN
        MESG = 'ERROR: Emission factors file does not contain ' //&
        &       'any pollutants'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( TOGIDX == 0 .AND. NHTOGIDX == 0 ) THEN
        MESG = 'ERROR: Emission factors file does not contain ' //&
        &       'data for both TOG/NONHAPTOG pollutants'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Build pollutant mapping table
!.........  Find emission pollutant in list of pollutants
    DO V = 1, NIPPA

        CPOL = TRIM( EANAM( V ) )

        J = INDEX1( CPOL, NMVSPOLS , MVSPOLNAMS )
        IF( J .LE. 0 ) THEN
            MESG = 'WARNING: Emission factors file does not ' //&
            &  'contain requested inventory pollutant '//TRIM( CPOL )
            CALL M3MESG( MESG )
        ELSE
            EMPOLIDX( V ) = J

        END IF

    END DO

!.............  Find model species in list of pollutants
    DO V = 1, NMSPC

        CSPC = EMNAM( V )
        J = INDEX1( CSPC, NMVSPOLS , MVSPOLNAMS )
        IF( J .LE. 0 ) THEN
            MESG = 'WARNING: Emission factors file does not ' //&
            &  'contain requested model species ' // TRIM( CSPC )
            CALL M3MESG( MESG )
        ELSE
            EMPOLIDX( NIPPA + V ) = J
        END IF

    END DO

!.........  Error message
    IF( EFLAG ) THEN
        MESG = 'Error occurred while reading emission factors file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Allocate memory to parse lines
    DEALLOCATE( SEGMENT )
    ALLOCATE( SEGMENT( NNONPOL + NMVSPOLS  ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

!.........  Read through file to determine maximum number of temperatures

!.........  Assumptions:
!             File will contain data for all SCCs in the inventory.
!             Lines are sorted by:
!                 temperature
!                 SCC (matching sorting of INVSCC)
!             Each SCC will have data for same set of emission temperatures, and
!               pollutants.

!.........  Limitations:
!             Program doesn't know if emission factors file is missing values.

!.........  Expected columns:
! MOVESScenarioID,yearID,monthID,FIPS,SCCsmoke,smokeProcID,avgSpeedBinID,temperature,relHumidity,THC,CO ...

    IF( NEMTEMPS .EQ. 0 ) THEN

        IREC = 0
        NEMTEMPS = 0
        PTMP = -999
        DO

            READ( TDEV, 93000, END=200, IOSTAT=IOS ) LINE

            IREC = IREC + 1

            IF( IOS .NE. 0 ) THEN
                WRITE( MESG, 94010 ) 'I/O error', IOS,&
                &  'reading emission factors file ' //&
                &  FILENAME // ' at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!.............  Skip blank or comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

!.............  Skip header line
            IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) CYCLE
            IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) CYCLE
            IF( LINE( 1:21 ) .EQ. 'HUMIDITY_ADJUSTED_NOX' ) CYCLE

!.............  Parse line into segments
            CALL PARSLINE( LINE, NNONPOL + NMVSPOLS , SEGMENT )

!.............  Check that county matches requested county
            IF( .NOT. CHKINT( SEGMENT( 4 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad reference county ' //&
                &  'FIPS code at line', IREC, 'of emission factors ' //&
                &  'file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            CFIP = SEGMENT( 4 )
            CALL PADZERO( CFIP )
            IF( CFIP .NE. MCREFIDX( REFIDX,1 ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Reference county ' //&
                &  'at line', IREC, 'of emission factors file ' //&
                &  'does not match county listed in MRCLIST file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!.............  Check that fuel month matches requested month
            IF( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad fuel month ' //&
                &  'at line', IREC, 'of emission factors file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( STR2INT( SEGMENT( 3 ) ) .NE. MONTH ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Fuel month at line',&
                &  IREC, 'of emission factors file does not match ' //&
                &  'fuel month listed in MRCLIST file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

!.............  Check temperature value
            IF( .NOT. CHKREAL( SEGMENT( 6 ) ) ) THEN
                WRITE( MESG, 94010 ) 'ERROR: Bad temperature value ' //&
                &  'at line', IREC, 'of emission factors file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            TMPVAL = STR2REAL( SEGMENT( 6 ) )
            IF( TMPVAL .NE. PTMP ) THEN

!.................  Check that temperatures are sorted
                IF( TMPVAL .LT. PTMP ) THEN
                    WRITE( MESG, 94010 ) 'ERROR: Temperature value ' //&
                    &  'at line', IREC, 'of emission factors file is ' //&
                    &  'smaller than previous temperature.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

                NEMTEMPS = NEMTEMPS + 1
                PTMP = TMPVAL
            END IF

        END DO

200     CONTINUE

        REWIND( TDEV )
    END IF

!.........  Allocate memory to store emission factors, temperature values
    IF( ALLOCATED( RPHEMFACS ) )  DEALLOCATE( RPHEMFACS )
    IF( ALLOCATED( EMTEMPS ) )    DEALLOCATE( EMTEMPS )

    ALLOCATE( RPHEMFACS( NINVSCC, NEMTEMPS, NMVSPOLS  ),&
    &           EMTEMPS( NEMTEMPS ),&
    &           LINVSCC( NINVSCC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'RPHEMFACS...LINVSCC', PROGNAME )

    RPHEMFACS = 0.  ! array
    EMTEMPS = 0.  ! array
    LINVSCC = .FALSE.

!.........  Read and store emission factors
    IREC = 0
    PSCC = ' '
    SCCIDX = 0
    PTMP = -999
    TMPIDX = 0
    NSCC = 0
    DO
        READ( TDEV, 93000, END=300, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 ) 'I/O error', IOS,&
            &  'reading emission factors file ' //&
            &  FILENAME // ' at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Skip header line
        IF( LINE( 1:12 ) .EQ. 'NUM_TEMP_BIN' ) CYCLE
        IF( LINE( 1:15 ) .EQ. 'MOVESScenarioID' ) CYCLE
        IF( LINE( 1:21 ) .EQ. 'HUMIDITY_ADJUSTED_NOX' ) CYCLE

!.............  Parse line into segments
        CALL PARSLINE( LINE, NNONPOL + NMVSPOLS , SEGMENT )

!.............  Set SCC index for current line
        TSCC = TRIM( SEGMENT( 5 ) )
        CALL PADZERO( TSCC )

        IF( TSCC .NE. PSCC ) THEN
            SKIPSCC = .FALSE.

            SCCIDX = SCCIDX + 1
            IF( SCCIDX .GT. NINVSCC ) THEN
                SCCIDX = 1
            END IF

!.................  Make sure the SCC at this line is what we expect
            IF( TSCC .EQ. INVSCC( SCCIDX ) ) THEN
                LINVSCC( SCCIDX ) = .TRUE.

            ELSE
!.....................  If the SCC is in the inventory, it means that the
!                       emissions factors are out of order and that's a
!                       problem; otherwise, for SCCs in the emission
!                       factors file that aren't in the inventory, we
!                       can skip this line
                J = FINDC( TSCC, NINVSCC, INVSCC )
                IF( J .LE. 0 ) THEN

!.........................  Emission factor SCC is not in the inventory
                    SKIPSCC = .TRUE.
                    SCCIDX = SCCIDX - 1
                    IF( SCCIDX .LT. 1 ) THEN
                        SCCIDX = NINVSCC
                    END IF
                ELSE
                    WRITE( MESG, 94010 ) 'Expected SCC ' //&
                    &  TRIM( INVSCC( SCCIDX ) ) // ' in emission ' //&
                    &  'factors file at line', IREC
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

            PSCC = TSCC
        END IF

        IF( SKIPSCC ) CYCLE
        NSCC = NSCC + 1

!.............  Set temperature index for current line
        TMPVAL = STR2REAL( SEGMENT( 6 ) )
        IF( TMPVAL .NE. PTMP ) THEN
            TMPIDX = TMPIDX + 1
            IF( TMPIDX .GT. NEMTEMPS ) THEN
                WRITE( MESG, 94010 ) 'Unexpected temperature value ' //&
                &  'in emission factors file at line', IREC
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            PTMP = TMPVAL
        END IF

        EMTEMPS( TMPIDX ) = TMPVAL

!.............  Store emission factors for each pollutant
        DO P = 1, NMVSPOLS

            EMVAL = STR2REAL( SEGMENT( NNONPOL + P ) )

            IF( EMVAL < 0 ) THEN    ! reset negative EF to zero
                EMVAL = 0.0
                WRITE( MESG, 94010 ) 'WARNING: Resetting negative ' //&
                &  ' emission factor to zero at line', IREC
                CALL M3MESG( MESG )
            END IF

            RPHEMFACS( SCCIDX, TMPIDX, P ) = EMVAL

        END DO

    END DO

300 CONTINUE

!.........  Error message when inventory SCC is missing in the lookup table
    DO I = 1, NINVSCC
        IF( .NOT. LINVSCC( I ) ) THEN
            MESG = 'WARNING: Following SCC "' // INVSCC( I ) //&
            &    '" is missing in this emission factors file'
            CALL M3MESG( MESG )
        END IF
    END DO

    IF( NSCC < 1 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No inventory SCC is found in emission factors file'
        CALL M3MSG2( MESG )
    END IF

    IF( EFLAG )  CALL M3EXIT( PROGNAME, 0, 0, '', 2 )

    CLOSE( TDEV )

    DEALLOCATE( SEGMENT, LINVSCC )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDRPHEMFACS
