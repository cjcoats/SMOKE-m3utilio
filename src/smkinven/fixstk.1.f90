
SUBROUTINE  FIXSTK( FDEV, NSRC )

!***********************************************************************
!  subroutine body starts at line 159
!
!  DESCRIPTION:
!       Read and use replacement stack parameters from file PSTK to fill in
!       stack parameters which are "missing" (i.e., negative). Also use
!       ultimate defaults (set as local parameters) when no other stack
!       parameters are available.
!
!  PRECONDITIONS REQUIRED:
!       Opened file with unit FDEV
!       Memory of arrays already allocated and with NSRC
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!       Subroutines: I/O API subroutines, CHECKMEM, FMTCSRC
!       Function: I/O API functions, GETFLINE
!
!  REVISION  HISTORY:
!       prototype 12/95 by CJC
!       copied by: mhouyoux
!       origin: fixstk.F 4.3
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
!***************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC, CSCC, CIFIP, CERPTYP,&
    &                    STKHT, STKDM, STKVE, STKTK

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   ARGUMENTS and their descriptions:

    INTEGER, INTENT( IN ) :: FDEV      ! unit no. for stack parameter file PSTK
    INTEGER, INTENT( IN ) :: NSRC      ! actual number of sources

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT
    LOGICAL, EXTERNAL :: USEEXPGEO

!...........    LOCAL VARIABLES and their descriptions:
!...........   SUBROUTINE PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'FIXSTK'  ! program name
    INTEGER      , PARAMETER :: NSEG = 6        ! number of fields for ORL FIREDATA input format

!...........   Temporary read arrays
    CHARACTER(10)      SEGMENT( NSEG ) ! segments of line

!.........  Unsorted arrays from stack replacements file
    INTEGER           , ALLOCATABLE :: INDXA( : ) ! sorting index
    REAL              , ALLOCATABLE :: SHTA ( : ) ! stack height
    REAL              , ALLOCATABLE :: SDMA ( : ) ! stack diameter
    REAL              , ALLOCATABLE :: STKA ( : ) ! stack temperature
    REAL              , ALLOCATABLE :: SVEA ( : ) ! stack velocity
    CHARACTER(FPSLEN3), ALLOCATABLE :: SFSCA( : ) ! FIPS code // SCC

!.........  Tables of stack parameter updates
!.........  Ultimate defaults
    REAL            :: HT0 = -1     !  ultimate fallback height
    REAL            :: DM0 = -1     !  ultimate fallback diameter
    REAL            :: TK0 = -1     !  ultimate fallback temperature
    REAL            :: VE0 = -1     !  ultimate fallback velocity

!.........  SCC-only table
    INTEGER                            NR1       ! number
    INTEGER           , ALLOCATABLE :: ID1 ( : ) ! index to unsorted
    CHARACTER(SCCLEN3), ALLOCATABLE :: TBL1( : ) ! SCC

!.........  SCC-country/state table
    INTEGER                            NR2       ! number
    INTEGER           , ALLOCATABLE :: ID2 ( : ) ! index to unsorted
    CHARACTER(STSLEN3), ALLOCATABLE :: TBL2( : ) ! co/st // scc

!.........  SCC-FIPS code table
    INTEGER                            NR3       ! number
    INTEGER           , ALLOCATABLE :: ID3 ( : ) ! index to unsorted
    CHARACTER(FPSLEN3), ALLOCATABLE :: TBL3( : ) ! FIPS code // scc

!.........  Other local variables
    INTEGER         I, J, K, L1, L2, L3, S   !  counters and indices

    INTEGER         NCNT    !  number of entry lines
    INTEGER         IOS     !  I/O error status
    INTEGER         IREC    !  current record number
    INTEGER         LDEV    !  log file unit number
    INTEGER         NLINE   !  number of lines
    INTEGER         NPSTK   !  number of PSTK entries
    INTEGER         SID     !  temporary state ID

    INTEGER, SAVE :: NWARN = 0            ! warning count
    INTEGER, SAVE :: MXWARN               ! max no. warnings

    REAL            HT      !  temporary height
    REAL            DM      !  temporary diameter
    REAL            TK      !  temporary exit temperature
    REAL            VE      !  temporary velocity
    REAL   MINHT, MAXHT      ! min/max stack heights
    REAL   MINDM, MAXDM      ! min/max stack diameter
    REAL   MINTK, MAXTK      !  min/max stack exit temperature
    REAL   MINVE, MAXVE      !  min/max stack velocity

    LOGICAL      :: EFLAG = .FALSE.  !  error flag
    LOGICAL      :: SFLAG = .FALSE.  ! true: skip applying stack para for fugitive src
    LOGICAL         DFLAG( NSRC )    ! true if source getting default parms

    CHARACTER(300)     BUFFER    !  temporary buffer
    CHARACTER(300)     MESG      !  message buffer
    CHARACTER(300)     LINE      !  read buffer for a line
    CHARACTER(FIPLEN3) CFIP      !  tmp character-string FIP
    CHARACTER(FIPLEN3) FIPZERO   !  zero buffer for FIPS code
    CHARACTER(STALEN3) CSTA      !  tmp country/state
    CHARACTER(CNYLEN3) CCNY      !  tmp county
    CHARACTER(SCCLEN3) TSCC      !  tmp SCC
    CHARACTER(SCCLEN3) SCCZERO   !  zero buffer for SCC
    CHARACTER(STSLEN3) CSTASCC   !  tmp country/state // SCC
    CHARACTER(STSLEN3) CSTASCCZ  !  zero buffer for cntry/state // SCC
    CHARACTER(FPSLEN3) CFIPSCC   !  tmp FIPS code // SCC
    CHARACTER(FPSLEN3) CFIPSCCZ  !  zero buffer for FIPS code // SCC

!***********************************************************************
!   begin body of subroutine FIXSTK

!.........   Get LOG file unit, so can write to directly (using M3MESG would
!            add too many spaces for some messages)
    LDEV = INIT3()

!.........   Get maximum number of warnings
    MXWARN = ENVINT( WARNSET , ' ', 100, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
    END IF

    CALL M3MSG2( 'Reading default stack parameters...' )

!.........  Skip if sources are fugitive
    MESG  = 'Skip replacing stack parameters for fugitive sources'
    SFLAG = ENVYN( 'NO_STACK_REPLACE_FUGITIVE',MESG,.FALSE.,IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "NO_STACK_REPLACE_FUGITIVE"', 2 )
    END IF

!.........  Define min/max ranges of stack parameters
    MESG = 'Define minimum stack height in unit of meter'
    MINHT = ENVREAL( 'MIN_STK_HEIGHT', MESG, 0.5, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_STK_HEIGHT"', 2 )
    END IF

    MESG = 'Define maximum stack height in unit of meter'
    MAXHT = ENVREAL( 'MAX_STK_HEIGHT', MESG, 5100.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MAX_STK_HEIGHT"', 2 )
    END IF

    MESG = 'Define minimum stack diameter in unit of meter'
    MINDM = ENVREAL( 'MIN_STK_DIAMETER', MESG, 0.01, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_STK_DIAMETER"', 2 )
    END IF

    MESG = 'Define maximum stack diameter in unit of meter'
    MAXDM = ENVREAL( 'MAX_STK_DIAMETER', MESG, 100.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MAX_STK_DIAMETER"', 2 )
    END IF

    MESG = 'Define minimum stack exit temperature in unit of Kelvin'
    MINTK = ENVREAL( 'MIN_STK_TEMPERATURE', MESG, 260.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_STK_TEMPERATURE"', 2 )
    END IF

    MESG = 'Define maximum stack exit temperature in unit of Kelvin'
    MAXTK = ENVREAL( 'MAX_STK_TEMPERATURE', MESG, 2000.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MAX_STK_TEMPERATURE"', 2 )
    END IF

    MESG = 'Define minimum stack exit velocity in unit of m/sec'
    MINVE = ENVREAL( 'MIN_STK_VELOCITY', MESG, 0.0001, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_STK_VELOCITY"', 2 )
    END IF

    MESG = 'Define maximum stack exit velocity in unit of m/sec'
    MAXVE = ENVREAL( 'MAX_STK_VELOCITY', MESG, 500.0, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MAX_STK_VELOCITY"', 2 )
    END IF

!.........  Get dimensions of input file
    NLINE = GETFLINE( FDEV, 'Stack replacement file')

!.........  Allocate memory for arrays.  Since this file is not likely to be
!           large, allocate all arrays based on number of lines in the file.
    ALLOCATE( INDXA( NLINE ),&
    &           SHTA( NLINE ),&
    &           SDMA( NLINE ),&
    &           STKA( NLINE ),&
    &           SVEA( NLINE ),&
    &          SFSCA( NLINE ),&
    &           TBL1( NLINE ),&
    &            ID1( NLINE ),&
    &           TBL2( NLINE ),&
    &            ID2( NLINE ),&
    &           TBL3( NLINE ),&
    &            ID3( NLINE ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INDXA...ID3', PROGNAME )

!.........  Create zero-filled buffers
    FIPZERO  = REPEAT( '0', FIPLEN3 )
    SCCZERO  = REPEAT( '0', SCCLEN3 )
    CSTASCCZ = REPEAT( '0', STSLEN3 )
    CFIPSCCZ = REPEAT( '0', FPSLEN3 )

!.........  Read the PSTK file until hit the end of the file
!.........  For now, require the SCC to be in quotes to use list formatting
    NCNT  = 0
    IREC  = 0
    DO I = 1, NLINE       !  head of input loop

        READ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN      !  I/O error
            WRITE( MESG,94010 ) 'I/O Error', IOS,&
            &       'reading stack replacements file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( BLKORCMT( LINE ) ) CYCLE

!            READ( FDEV, *, IOSTAT=IOS ) FIP, TSCC, HT, DM, TK, VE
!.............  Get lines
        CALL PARSLINE( LINE, NSEG, SEGMENT )

        NCNT = NCNT + 1      ! actual line #s after skipping blank and comments

        CFIP = SEGMENT( 1 )
        TSCC = SEGMENT( 2 )

        CALL PADZERO( CFIP )
        CALL PADZERO( TSCC )

        INDXA( NCNT ) = NCNT
        SFSCA( NCNT ) = CFIP // TSCC
        SHTA ( NCNT ) = STR2REAL( SEGMENT( 3 ) )
        SDMA ( NCNT ) = STR2REAL( SEGMENT( 4 ) )
        STKA ( NCNT ) = STR2REAL( SEGMENT( 5 ) )
        SVEA ( NCNT ) = STR2REAL( SEGMENT( 6 ) )

    END DO

    NPSTK = NCNT

!.........  Sort PSTK data
    CALL SORTIC( NPSTK, INDXA, SFSCA )

!.........  Disaggregate PSTK data into 4 categories
    NR1   = 0
    NR2   = 0
    NR3   = 0
    L1    = STALEN3  ! without county
    L2    = L1 + 1
    L3    = FIPLEN3 + 1
    DO I = 1, NPSTK

        J    = INDXA( I )
        CFIP = SFSCA( J )( 1 :FIPLEN3 )
        CSTA = SFSCA( J )( 1 :L1 )
        CCNY = SFSCA( J )( L2:FIPLEN3 )
        TSCC = SFSCA( J )( L3:FPSLEN3 )

        IF( CFIP .EQ. FIPZERO .AND.&
        &    TSCC .EQ. SCCZERO       ) THEN  !  fallback default
            HT0 = SHTA ( J )
            DM0 = SDMA ( J )
            TK0 = STKA ( J )
            VE0 = SVEA ( J )

        ELSEIF( CFIP .EQ. FIPZERO ) THEN    !  SCC only
            NR1 = NR1 + 1
            TBL1( NR1 ) = TSCC
            ID1 ( NR1 ) = J

        ELSE IF( CCNY .EQ. '000' .AND. .NOT. USEEXPGEO() ) THEN !  state and SCC
            NR2 = NR2 + 1
            TBL2( NR2 ) = CSTA // TSCC
            ID2 ( NR2 ) = J

        ELSE                                !  FIP and SCC
            NR3 = NR3 + 1
            TBL3( NR3 ) = CFIP // TSCC
            ID3 ( NR3 ) = J

        END IF

    END DO ! End loop on NPSTK

!.........  Bound stack parameters to minima and maxima values
!.........  This is in a separate loop to permit better reporting
!.........  Watch out for negative numbers or zeroes, because these are
!.........  the missing stack parameters, which should get defaults.

    CALL M3MSG2( 'Bounding MIN and MAX stack parameters...' )

    DO S = 1, NSRC

        HT   = STKHT( S )
        DM   = STKDM( S )
        TK   = STKTK( S )
        VE   = STKVE( S )
        TSCC = CSCC ( S )

        IF ( HT .GT. MAXHT .OR.&
        &   ( HT .LT. MINHT .AND. HT .GT. 0 ) .OR.&
        &     DM .GT. MAXDM .OR.&
        &   ( DM .LT. MINDM .AND. DM .GT. 0 ) .OR.&
        &     TK .GT. MAXTK .OR.&
        &   ( TK .LT. MINTK .AND. TK .GT. 0 ) .OR.&
        &     VE .GT. MAXVE .OR.&
        &   ( VE .LT. MINVE .AND. VE .GT. 0 ) ) THEN

            NWARN = NWARN + 1
            IF( NWARN <= MXWARN ) THEN
                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2 )
                WRITE( MESG,94010 ) BUFFER( 1:L2 ) // ' SCC: '// TSCC
                CALL M3MESG( MESG )
            END IF

        END IF

        IF ( HT .GT. MAXHT ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94030 ) 'Height', HT, MAXHT
            HT = MAXHT

        ELSEIF( HT .LT. MINHT .AND. HT .GT. 0 ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94040 ) 'Height', HT, MINHT
            HT = MINHT

        END IF

        IF ( DM .GT. MAXDM ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94030 ) '  Diam', DM, MAXDM
            DM = MAXDM

        ELSEIF( DM .LT. MINDM .AND. DM .GT. 0 ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94040 ) '  Diam', DM, MINDM
            DM = MINDM

        END IF

        IF ( TK .GT. MAXTK )THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94030 ) '  Temp', TK, MAXTK
            TK = MAXTK

        ELSEIF( TK .LT. MINTK .AND. TK .GT. 0 ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94040 ) '  Temp', TK, MINTK
            TK = MINTK

        END IF

        IF ( VE .GT. MAXVE )THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94030 ) ' Veloc', VE, MAXVE
            VE = MAXVE

        ELSEIF( VE .LT. MINVE .AND. VE .GT. 0 ) THEN
            IF( NWARN <= MXWARN )&
            &    WRITE( LDEV,94040 ) ' Veloc', VE, MINVE
            VE = MINVE

        END IF

        STKHT( S ) = HT
        STKDM( S ) = DM
        STKTK( S ) = TK
        STKVE( S ) = VE

    END DO ! Loop on sources

    CALL M3MSG2( 'Fixing MISSING stack parameters...' )

!...........   Now do replacements of MISSING stack parameters:
!...........   4 passes -- ht, dm, tk, ve
!...........   Treat parameters equal to 0 as missing

    NWARN = 0
    DO S = 1, NSRC

        K = 0                ! Initialize K to test if replacements made
        DFLAG( S ) = .FALSE. ! Initialize DFLAG to test if defaults used

!.............  Skip for fugitive source (01)
        IF( SFLAG ) THEN
            IF( STR2INT( CERPTYP(S) ) == 1 ) THEN
                STKHT( S ) = 0.0
                STKDM( S ) = 0.0
                STKTK( S ) = 0.0
                STKVE( S ) = 0.0
                CYCLE
            END IF
        END IF

!.............  Set up temporary character strings
        CFIP = CIFIP( S )
        CSTA = CFIP( 1:STALEN3 )
        TSCC = CSCC( S )
        CFIPSCC  = CFIP // TSCC
        CFIPSCCZ = CFIP // SCCZERO
        CSTASCC  = CSTA // TSCC
        CSTASCCZ = CSTA // SCCZERO

        IF ( STKHT( S ) .LE. 0.0 ) THEN

            K = FINDC( CFIPSCC, NR3, TBL3 )

            IF( K .LE. 0 )&
            &    K = FINDC( CFIPSCCZ, NR3, TBL3 )

            IF ( K .GT. 0 ) THEN
                J  = ID3 ( K )
                HT = SHTA( J )
                DM = SDMA( J )
                TK = STKA( J )
                VE = SVEA( J )
            ELSE
                K = FINDC( CSTASCC, NR2, TBL2 )

                IF( K .LE. 0 )&
                &    K = FINDC( CSTASCCZ, NR2, TBL2 )

                IF ( K .GT. 0 ) THEN
                    J  = ID2 ( K )
                    HT = SHTA( J )
                    DM = SDMA( J )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE

                    K = FINDC( TSCC, NR1, TBL1 )
                    IF ( K .GT. 0 ) THEN
                        J  = ID1 ( K )
                        HT = SHTA( J )
                        DM = SDMA( J )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE
                        DFLAG( S ) = .TRUE.

                    END IF
                END IF
            END IF

            IF( .NOT. DFLAG( S ) ) THEN

                NWARN = NWARN + 1

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                WRITE( MESG,94010 )&
                &       BUFFER( 1:L2 ) // ' SCC: ' // TSCC //&
                &       CRLF() // BLANK5 //&
                &       '             Old        New'
                IF( NWARN <= MXWARN ) CALL M3MESG( MESG )

                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) 'Height', STKHT( S ), HT
                STKHT( S ) = HT

                IF ( STKDM( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM
                    STKDM( S ) = DM
                END IF

                IF ( STKTK( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                    STKTK( S ) = TK
                END IF

                IF ( STKVE( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                    STKVE( S ) = VE
                END IF
            END IF

        END IF      !  if stack height bad

        IF ( STKDM( S ) .LE. 0.0 ) THEN
            K = FINDC( CFIPSCC, NR3, TBL3 )

            IF( K .LE. 0 )&
            &    K = FINDC( CFIPSCCZ, NR3, TBL3 )

            IF ( K .GT. 0 ) THEN
                J  = ID3 ( K )
                DM = SDMA( J )
                TK = STKA( J )
                VE = SVEA( J )
            ELSE
                K = FINDC( CSTASCC, NR2, TBL2 )

                IF( K .LE. 0 )&
                &    K = FINDC( CSTASCCZ, NR2, TBL2 )

                IF ( K .GT. 0 ) THEN
                    J  = ID2 ( K )
                    DM = SDMA( J )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE
                    K = FINDC( TSCC, NR1, TBL1 )
                    IF ( K .GT. 0 ) THEN
                        J  = ID1 ( K )
                        DM = SDMA( J )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE
                        DFLAG( S ) = .TRUE.

                    END IF
                END IF
            END IF

            IF( .NOT. DFLAG( S ) ) THEN

                NWARN = NWARN + 1

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                WRITE( MESG,94010 )&
                &       BUFFER( 1:L2 ) // ' SCC: ' // TSCC //&
                &       CRLF() // BLANK5 //&
                &       '             Old        New'
                IF( NWARN <= MXWARN ) CALL M3MESG( MESG )

                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM
                STKDM( S ) = DM

                IF ( STKTK( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                    STKTK( S ) = TK
                END IF

                IF ( STKVE( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                    STKVE( S ) = VE
                END IF
            END IF

        END IF      !  if stack diameter bad

        IF ( STKTK( S ) .LE. 0.0 ) THEN
            K = FINDC( CFIPSCC, NR3, TBL3 )

            IF( K .LE. 0 )&
            &    K = FINDC( CFIPSCCZ, NR3, TBL3 )

            IF ( K .GT. 0 ) THEN
                J  = ID3 ( K )
                TK = STKA( J )
                VE = SVEA( J )
            ELSE
                K   = FINDC( CSTASCC, NR2, TBL2 )

                IF( K .LE. 0 )&
                &    K = FINDC( CSTASCCZ, NR2, TBL2 )

                IF ( K .GT. 0 ) THEN
                    J  = ID2 ( K )
                    TK = STKA( J )
                    VE = SVEA( J )
                ELSE
                    K = FINDC( TSCC, NR1, TBL1 )
                    IF ( K .GT. 0 ) THEN
                        J  = ID1 ( K )
                        TK = STKA( J )
                        VE = SVEA( J )
                    ELSE
                        DFLAG( S ) = .TRUE.

                    END IF
                END IF
            END IF

            IF( .NOT. DFLAG( S ) ) THEN

                NWARN = NWARN + 1

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                WRITE( MESG,94010 )&
                &       BUFFER( 1:L2 ) // ' SCC: ' // TSCC //&
                &       CRLF() // BLANK5 //&
                &       '             Old        New'
                IF( NWARN <= MXWARN ) CALL M3MESG( MESG )

                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK
                STKTK( S ) = TK

                IF ( STKVE( S ) .LE. 0 ) THEN
                    IF( NWARN <= MXWARN )&
                    &    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                    STKVE( S ) = VE
                END IF
            END IF

        END IF      !  if stack exhaust temperature bad

        IF ( STKVE( S ) .LE. 0.0 ) THEN
            K = FINDC( CFIPSCC, NR3, TBL3 )

            IF( K .LE. 0 )&
            &    K = FINDC( CFIPSCCZ, NR3, TBL3 )

            IF ( K .GT. 0 ) THEN
                J  = ID3 ( K )
                VE = SVEA( J )
            ELSE
                K = FINDC( CSTASCC, NR2, TBL2 )

                IF( K .LE. 0 )&
                &    K = FINDC( CSTASCCZ, NR2, TBL2 )

                IF ( K .GT. 0 ) THEN
                    J  = ID2 ( K )
                    VE = SVEA( J )
                ELSE
                    K = FINDC( TSCC, NR1, TBL1 )
                    IF ( K .GT. 0 ) THEN
                        J  = ID1 ( K )
                        VE = SVEA( J )
                    ELSE
                        DFLAG( S ) = .TRUE.

                    END IF
                END IF
            END IF

            IF( .NOT. DFLAG( S ) ) THEN

                NWARN = NWARN + 1

                CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2)
                WRITE( MESG,94010 )&
                &       BUFFER( 1:L2 ) // ' SCC: ' // TSCC //&
                &       CRLF() // BLANK5 //&
                &       '             Old        New'
                IF( NWARN <= MXWARN ) CALL M3MESG( MESG )

                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE
                STKVE( S ) = VE
            END IF

        END IF      !  if stack exhaust velocity bad

    END DO  !  end loop on sources for fixing missing stack parameters

!.........  Apply ultimate fallback parameters, and write report
!.........  This is in a separate loop to permit better reporting

    CALL M3MESG( 'Ultimate fallback stack parameters report:' )

    NWARN = 0
    DO S = 1, NSRC

        IF( DFLAG( S ) ) THEN

            CALL FMTCSRC( CSOURC(S), 7, BUFFER, L2 )

!.................  Error msg when there are no default x-reference available
            IF( HT0 < 0.0 .AND. DM0 < 0.0 .AND.&
            &    TK0 < 0.0 .AND. VE0 < 0.0       ) THEN

                WRITE( MESG,94010 )&
                &  'ERROR: No PSTK cross-reference ' //&
                &  'available (and no default) for:' //&
                &  CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            NWARN = NWARN + 1
            WRITE( MESG,94010 )&
            &       BUFFER( 1:L2 ) // ' SCC: ' // TSCC //&
            &       CRLF() // BLANK5 //&
            &       '             Old        New'
            IF( NWARN <= MXWARN ) CALL M3MESG( MESG )

            IF ( STKHT( S ) .LE. 0 ) THEN
                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) 'Height', STKHT( S ), HT0
                STKHT( S ) = HT0
            END IF

            IF ( STKDM( S ) .LE. 0 ) THEN
                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) '  Diam', STKDM( S ), DM0
                STKDM( S ) = DM0
            END IF

            IF ( STKTK( S ) .LE. 0 ) THEN
                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) '  Temp', STKTK( S ), TK0
                STKTK( S ) = TK0
            END IF

            IF ( STKVE( S ) .LE. 0 ) THEN
                IF( NWARN <= MXWARN )&
                &    WRITE( LDEV,94020 ) ' Veloc', STKVE( S ), VE0
                STKVE( S ) = VE0
            END IF

        END IF

    END DO  ! Loop through sources for applying ultimate fallbacks

    DEALLOCATE( INDXA )
    DEALLOCATE( SFSCA )
    DEALLOCATE( SHTA )
    DEALLOCATE( SDMA )
    DEALLOCATE( STKA )
    DEALLOCATE( SVEA )
    DEALLOCATE( TBL1 )
    DEALLOCATE( ID1 )
    DEALLOCATE( TBL2 )
    DEALLOCATE( ID2 )
    DEALLOCATE( TBL3 )
    DEALLOCATE( ID3 )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

93000 FORMAT( A )

94010 FORMAT( 10 ( A, :, I8, :, 1X ) )

94020 FORMAT( 7X, A, 2X, E10.3, 1X, E10.3 )

94030 FORMAT( 7X, A6, 1X, '> max.  Change from ',&
    &        E10.3, ' to ', E10.3 )

94040 FORMAT( 7X, A6, 1X, '< min.  Change from ',&
    &        E10.3, ' to ', E10.3 )

94050 FORMAT( A4, A6 )

94300 FORMAT( A, I2.2, A, I2.2, A )

END SUBROUTINE FIXSTK

