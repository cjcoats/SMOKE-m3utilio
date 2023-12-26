
SUBROUTINE RPELVCFG( FDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      Subroutine RPELVCFG reads the PELVCONFIG file that is used to set
    !      the elevated groups, the elevated sources, and the plume-in-grid
    !      (PinG) sources.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/2001 by M Houyoux
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

    !.......   MODULES for public variables
    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: PKT_IDX, ELG_IDX, PNG_IDX, ELV_IDX, &
                        SPCF_NAND, INSPCIFY, LIN_SPCIFY,    &
                        SPCF_NOR, NCRTSYBL, CRTSYBL

    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: GRPVALS, GRPTYPES,                   &
                       ELVVALS, ELVTYPES, ELVCHRS,          &
                       PNGVALS, PNGTYPES, PNGCHRS,          &
                       EVPEMIDX, EVPESTAT, EVPPSTAT,        &
                       NGRPVAR, NEVPVAR, NEVPEMV,           &
                       NGRPCRIT, NELVCRIT, NPNGCRIT,        &
                       MXGRPCHK, MXELVCHK, MXPNGCHK,        &
                       LCUTOFF, LPNGRNK, LELVRNK

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: EINAM, NIPOL

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: FDEV           ! input file unit

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local arrays
    INTEGER     TMPIDX( NIPOL )      ! tmp index for flagging pollutants
    LOGICAL     EISTAT( NIPOL )      ! true: pollutant used as criteria

    !.......   Other local variables
    INTEGER         I, V       ! indices and counters

    INTEGER         IOS                      ! i/o status
    INTEGER,SAVE :: NALLPOL = 0              ! number of pols used as criteria
    INTEGER         NLINES                   ! number of lines

    LOGICAL,SAVE :: EFLAG = .FALSE.          ! true: error found

    CHARACTER(300)  MESG                     ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'RPELVCFG'     ! program name

    !***********************************************************************
    !   begin body of subroutine RPELVCFG

    !.......  Write status message
    MESG = 'Reading elevated source configuration file...'
    CALL M3MSG2( MESG )

    !.......  Initialize the number of variables for grouping, selecting PinG
    !           sources, and selecting elevated sources.
    NGRPVAR = 5         ! The number of stack parmeters
    NEVPVAR = MAX( HT_IDX, DM_IDX, TK_IDX, VE_IDX, FL_IDX,    &
                   SRC_IDX, FIP_IDX, PLT_IDX )

    !.......  Get number of lines in file
    NLINES = GETFLINE( FDEV, 'Elevated configuration' )
    CALL CHECKMEM( IOS, 'EISTAT', PROGNAME )
    TMPIDX = 0
    EISTAT = .FALSE.

    !.......  Read through file to determine how many ORs and ANDs for each
    !           section of file
    CALL READ_PELVCONFIG( FDEV, 'COUNT', EFLAG )

    !.......  Abort if error
    IF ( EFLAG ) THEN
        MESG = 'Problem scanning elevated configuration file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    NEVPVAR = NEVPVAR + NALLPOL

    NEVPEMV = 0
    DO V = 1, NIPOL
        IF ( EISTAT( V ) ) NEVPEMV = NEVPEMV + 1
    END DO

    !.......  Limit group criteria to 1 for later usage in ASGNGRPS
    NGRPCRIT = MAX( NGRPCRIT, 1 )
    MXGRPCHK = MAX( MXGRPCHK, 1 )

    !.......  Allocate memory for criteria arrays
    ALLOCATE( GRPVALS( NGRPCRIT, MXGRPCHK, NGRPVAR ),   &
             GRPTYPES( NGRPCRIT, MXGRPCHK, NGRPVAR ),   &
              ELVVALS( NELVCRIT, MXELVCHK, NEVPVAR ),   &
              ELVCHRS( NELVCRIT, MXELVCHK, NEVPVAR ),   &
             ELVTYPES( NELVCRIT, MXELVCHK, NEVPVAR ),   &
              PNGVALS( NPNGCRIT, MXPNGCHK, NEVPVAR ),   &
              PNGCHRS( NPNGCRIT, MXPNGCHK, NEVPVAR ),   &
             PNGTYPES( NPNGCRIT, MXPNGCHK, NEVPVAR ),   &
             EVPEMIDX( NEVPEMV ),                       &
             EVPESTAT( NIPOL ),                         &
             EVPPSTAT( NIPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EVPPSTAT', PROGNAME )

    GRPVALS  = 0.
    GRPTYPES = '<'     ! bogus comparison type. By default, no groups get set.
    ELVVALS  = 0.
    ELVCHRS  = ' '
    ELVTYPES = ' '
    PNGVALS  = 0.
    PNGCHRS  = ' '
    PNGTYPES = ' '
    EVPEMIDX = 0
    EVPESTAT = .FALSE.
    EVPPSTAT = .FALSE.

    !.......  Store index from pollutants used as selection criteria to master
    I = 0
    DO V = 1, NIPOL
        IF ( EISTAT( V ) ) THEN
            I = I + 1
            EVPEMIDX( I ) = V
            TMPIDX  ( V ) = I
        END IF
    END DO

    !.......  Read through file to store criteria arrays
    CALL READ_PELVCONFIG( FDEV, 'STORE', EFLAG )

    IF ( EFLAG ) THEN
        MESG = 'Problem reading elevated configuration file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN


CONTAINS        !    !******************  INTERNAL SUBPROGRAMS  ********************

    SUBROUTINE READ_PELVCONFIG( FDEV, READTYPE, EFLAG )

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: FDEV               ! input file unit
        CHARACTER(*), INTENT (IN) :: READTYPE           ! Reading type
        LOGICAL  , INTENT(IN OUT) :: EFLAG              ! true: error found

        !.......  External functions
        LOGICAL, EXTERNAL :: BLKORCMT
        LOGICAL, EXTERNAL :: CHKREAL
        INTEGER, EXTERNAL :: GETNLIST

        !.......  Subprogram local allocatable arrays
        CHARACTER(32), ALLOCATABLE :: SEGMENT( : )

        !.......  Local variables
        INTEGER       I, I1, I2, I3, J, K, L, N, V          ! counters and indices

        INTEGER       IOS              ! i/o status
        INTEGER       IREC             ! record counter
        INTEGER       MX_IDX           ! maximum of stack parm indices
        INTEGER    :: NS = 1           ! no. segments in line
        INTEGER       RCNT             ! record count

        REAL          VAL              ! tmp value

        LOGICAL, SAVE :: FIRSTSET = .TRUE.         ! true: 1st time group crit found

        CHARACTER(600) BUFFER           ! tmp line buffer as uppercase
        CHARACTER(600) LINE             ! tmp line buffer
        CHARACTER(300) MESG             ! mesg buffer

        !----------------------------------------------------------------------

        !.......  Rewind input file
        REWIND( FDEV )

        !.......  Set local constant variables
        MX_IDX = MAX( HT_IDX, DM_IDX, TK_IDX, VE_IDX, FL_IDX,    &
                      RISE_IDX, SRC_IDX, FIP_IDX, PLT_IDX )

        !.......  Read file with different steps depending on READTYPE arguments
        DO IREC = 1, NLINES

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,               &
                  'reading elevated source configuration file '//   &
                  'at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Skip blank lines and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

            !.......  Screen for appended comments and remove them
            CALL RMCOMMNT( '##', LINE )

            !.......  Left-justify and convert line to upper case
            BUFFER = ADJUSTL( LINE )
            CALL UPCASE( BUFFER )

            !.......  Deallocate segment from previous iteration
            IF ( ALLOCATED( SEGMENT ) ) DEALLOCATE( SEGMENT )

            !.......  Allocate memory for segments and parse line into segments
            L = LEN_TRIM( BUFFER )
            NS = GETNLIST( L, BUFFER )
            ALLOCATE( SEGMENT( NS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'SEGMENT', PROGNAME )

            CALL PARSLINE( BUFFER, NS, SEGMENT )

            !.......  Interpret line of code.  Set global variables in MODREPRT.
            CALL PRCLINRC( IREC, NS, BUFFER, SEGMENT )

            !.......  Skip line if it is the start or end of a packet
            IF ( .NOT. INSPCIFY .OR. LIN_SPCIFY ) CYCLE

            !.......  Store the maximum number of ANDs, set the number of ORs (the
            !                   number of records), and update the number of variables
            SELECT CASE( PKT_IDX )
              CASE( ELG_IDX )
                NGRPCRIT = SPCF_NOR
                MXGRPCHK = MAX( MXGRPCHK, SPCF_NAND )

                !.......  Select PinG sources
              CASE( PNG_IDX )
                NPNGCRIT = SPCF_NOR
                MXPNGCHK = MAX( MXPNGCHK, SPCF_NAND )

                !.......  Search segments for a match with pollutant names
                N = 0
                DO I = 1, NS, 4
                    N = N + 1
                    I1 = ( N-1 )*4 + 1
                    I2 = ( N-1 )*4 + 2

                    J = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )

                    !...........  If pollutant found with ranking, increase pollutant
                    !                           ranking counter
                    IF ( J             .GT. 0     .AND.    &
                         SEGMENT( I2 ) .EQ. 'TOP'       ) THEN
                        IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                        EISTAT( J ) = .TRUE.

                    !...........  If pollutant found without ranking, increase
                    !                           pollutant value counter
                    ELSE IF ( J .GT. 0 ) THEN
                        IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                        EISTAT( J ) = .TRUE.

                    END IF
                END DO

                !.......  Select elevated sources
              CASE( ELV_IDX )
                NELVCRIT = SPCF_NOR
                MXELVCHK = MAX( MXELVCHK, SPCF_NAND )

                !.......  Search segments for a match with pollutant names, and
                !                       increase counter if found
                N = 0
                DO I = 1, NS, 4
                    N = N + 1
                    I1 = ( N-1 )*4 + 1
                    I2 = ( N-1 )*4 + 2

                    J = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )

                    !...........  If pollutant found with ranking, increase pollutant
                    !                           ranking counter
                    IF ( J             .GT. 0     .AND.    &
                         SEGMENT( I2 ) .EQ. 'TOP'       ) THEN
                        IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                        EISTAT( J ) = .TRUE.

                    !...........  If pollutant found without ranking, increase
                    !                           pollutant value counter
                    ELSE IF ( J .GT. 0 ) THEN
                        IF ( .NOT. EISTAT(J) ) NALLPOL = NALLPOL + 1
                        EISTAT( J ) = .TRUE.

                    END IF
                END DO

            END SELECT

            !.......  End loop if just counting the number of records
            IF ( READTYPE .EQ. 'COUNT' ) CYCLE

            !.......  Check fields
            N = 0
            DO I = 1, NS, 4
                N  = N + 1
                I1 = ( N-1 )*4 + 1
                I2 = ( N-1 )*4 + 2
                I3 = ( N-1 )*4 + 3

                IF ( I2 .GT. NS .OR. I3 .GT. NS ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 ) 'ERROR: Bad use of ' //    &
                           'delimeters at line ', IREC,    &
                           '. Line skipped.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

                !.......  Check and store first part of each AND component
                SELECT CASE( SEGMENT( I1 ) )
                  CASE( 'HEIGHT', 'HT' )
                    K = HT_IDX
                  CASE( 'DIAMETER', 'DM' )
                    K = DM_IDX
                  CASE( 'TEMPERATURE', 'TK' )
                    K = TK_IDX
                  CASE( 'VELOCITY', 'VE' )
                    K = VE_IDX
                  CASE( 'FLOW', 'FL' )
                    K = FL_IDX
                  CASE( 'RISE' )
                    K = RISE_IDX

                    !...........  Make sure RISE is not used for setting groups
                    IF( PKT_IDX .EQ. ELG_IDX ) THEN
                        WRITE( MESG,94010 )    &
                              'WARNING: Variable "RISE" '//    &
                              'cannot be used to set elevated '//    &
                              'groups at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                    !...........  Flag file as using cutoff approach
                    LCUTOFF = .TRUE.

                  CASE( 'SOURCE' )
                    K = SRC_IDX

                    !...........  Make sure SOURCE is not used for setting groups
                    IF( PKT_IDX .EQ. ELG_IDX ) THEN
                        WRITE( MESG,94010 )    &
                              'WARNING: Variable "SOURCE" '//    &
                              'cannot be used to set elevated '//    &
                              'groups at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                  CASE( 'FIPS' )
                    K = FIP_IDX

                    !...........  Make sure FIPS is not used for setting groups
                    IF( PKT_IDX .EQ. ELG_IDX ) THEN
                        WRITE( MESG,94010 )    &
                              'WARNING: Variable "FIPS" '//    &
                              'cannot be used to set elevated '//    &
                              'groups at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                  CASE( 'PLANT' )
                    K = PLT_IDX

                    !...........  Make sure PLANT is not used for setting groups
                    IF( PKT_IDX .EQ. ELG_IDX ) THEN
                        WRITE( MESG,94010 )    &
                              'WARNING: Variable "PLANT" '//    &
                              'cannot be used to set elevated '//    &
                              'groups at line', IREC
                        CALL M3MSG2( MESG )
                        CYCLE
                    END IF

                    !.......  Make sure comparison type is consistent with
                    !                       plant
                    IF( SEGMENT( I2 ) .NE. 'IS' ) THEN
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )    &
                          'ERROR: Comparison type "'// TRIM( SEGMENT( I2 ) ) //     &
                          '" at line', IREC, 'is not allowed with PLANT criteria.'
                        CALL M3MSG2( MESG )
                    END IF

                    !.......  Otherwise, determine if the field is a pollutant, and
                    !                       if so, store index accordingly
                  CASE DEFAULT

                    V = INDEX1( SEGMENT( I1 ), NIPOL, EINAM )
                    IF ( V .GT. 0 ) THEN

                        SELECT CASE( PKT_IDX )
                          CASE( ELG_IDX )
                            WRITE( MESG,94010 )    &
                              'WARNING: Pollutant "'// TRIM( EINAM( V ) ) //    &
                              '" cannot be used to set elevated groups at line', IREC
                            CALL M3MSG2( MESG )
                            CYCLE                   ! To next line of file

                          CASE( PNG_IDX )
                            EVPPSTAT( V ) = .TRUE.

                          CASE( ELV_IDX )
                            EVPESTAT( V ) = .TRUE.

                        END SELECT

                        K = MX_IDX + TMPIDX( V )

                    !............  Otherwise, error
                    ELSE
                        EFLAG = .TRUE.
                        WRITE( MESG,94010 )    &
                          'ERROR: Criteria "'// TRIM( SEGMENT( I1 ) ) //    &
                          '" at line', IREC, 'is not recognized.'
                        CALL M3MSG2( MESG )

                    END IF

                END SELECT

                !.......  Check and store second part of each AND component
                !.......  These must match what is recognized in Evalcrit routine
                J = INDEX1( SEGMENT( I2 ), NCRTSYBL, CRTSYBL )
                IF ( J .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )                                             &
                           'ERROR: Comparison type "'// TRIM( SEGMENT( I2 ) ) //    &
                           '" at line', IREC, 'is not recognized.'
                    CALL M3MSG2( MESG )

                !.......  Store type
                ELSE
                    SELECT CASE( PKT_IDX )
                      CASE( ELG_IDX )

                        !...........  Make sure TOP is not used for setting groups
                        IF( SEGMENT( I2 ) .EQ. 'TOP' ) THEN
                            WRITE( MESG,94010 )    &
                              'WARNING: Comparison type "TOP" '//    &
                              'ignored at line', IREC, 'for ' //    &
                              'stack group specification'
                            CALL M3MSG2( MESG )
                            CYCLE
                        !...........  Otherwise store it
                        ELSE

                            IF( FIRSTSET ) THEN
                                GRPTYPES = ' '
                                FIRSTSET = .FALSE.
                            END IF
                            GRPTYPES(NGRPCRIT,N,K) = SEGMENT( I2 )

                        END IF

                      CASE( PNG_IDX )
                        PNGTYPES( NPNGCRIT, N, K ) = SEGMENT( I2 )
                        IF( SEGMENT(I2) .EQ. 'TOP' ) LPNGRNK= .TRUE.
                      CASE( ELV_IDX )
                        ELVTYPES( NELVCRIT, N, K ) = SEGMENT( I2 )
                        IF( SEGMENT(I2) .EQ. 'TOP' ) LELVRNK= .TRUE.
                    END SELECT

                END IF

                !.......  Check value of each AND component...
                !.......  If plant is specified, then store string information...
                IF ( K .EQ. PLT_IDX ) THEN

                    SELECT CASE( PKT_IDX )
                      CASE( ELG_IDX )

                      CASE( PNG_IDX )
                        PNGCHRS( NPNGCRIT, N, K ) = SEGMENT( I3 )
                      CASE( ELV_IDX )
                        ELVCHRS( NELVCRIT, N, K ) = SEGMENT( I3 )
                    END SELECT

                !.......  If plant not specified, then check if integer or real
                ELSE IF ( .NOT. CHKREAL( SEGMENT( I3 ) ) ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )                                     &
                           'ERROR: Value "'// TRIM( SEGMENT( I3 ) ) //      &
                           '" at line', IREC, 'is an invalid numeric value.'
                    CALL M3MSG2( MESG )

                !.......  Store value of each AND component as a real value
                ELSE

                    VAL = STR2REAL( SEGMENT( I3 ) )

                    SELECT CASE( PKT_IDX )
                      CASE( ELG_IDX )
                        GRPVALS( NGRPCRIT, N, K ) = VAL
                      CASE( PNG_IDX )
                        PNGVALS( NPNGCRIT, N, K ) = VAL
                      CASE( ELV_IDX )
                        ELVVALS( NELVCRIT, N, K ) = VAL
                    END SELECT

                END IF

            END DO              ! End loop over ANDs on line

        END DO              ! End loop over lines of file

        !.......  Check if there are no group criteria and write warning
        IF( READTYPE .EQ. 'COUNT' .AND.    &
            NGRPCRIT .LE. 0             ) THEN
            MESG = 'WARNING: No grouping criteria in PELVCONFIG ' //    &
                   'will likely result in wasteful ' //                 &
                   CRLF() // BLANK10 // ' separation of PinG '//        &
                   'and/or elevated records that could be grouped.'
            CALL M3MSG2( MESG )
        END IF

        RETURN

        !.......  Problem(s) reading input file...
999     WRITE( MESG,94010 )         &
            'INTERNAL ERROR: Unexpected end of file at line', IREC
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        !........  FORMAT  STATEMENTS   ............

        !.......   Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE READ_PELVCONFIG

END SUBROUTINE RPELVCFG
