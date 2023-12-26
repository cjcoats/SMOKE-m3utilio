
SUBROUTINE GETMRGEV

!***********************************************************************
!  subroutine GETMRGEV body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 2/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!           Bug-fix in ENVINT calls at lines 253,260
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
    USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,&
    &                    AUFLAG, MUFLAG, PUFLAG,&
    &                    ARFLAG, MRFLAG, PRFLAG,&
    &                    AFLAG_BD, MFLAG_BD, PFLAG_BD,&
    &                    XFLAG, TFLAG, SFLAG, LFLAG,&
    &                    PINGFLAG, ELEVFLAG, EXPLFLAG, VARFLAG,&
    &                    LMETCHK, LMKTPON, LGRDOUT, LREPCNY,&
    &                    LREPSTA, LREPINV, LREPSPC, LREPCTL,&
    &                    LREPANY, LAVEDAY, INLINEFLAG, SRCGRPFLAG,&
    &                    SUBSECFLAG

!...........  This module contains the information about the source category
    USE MODINFO, ONLY:   INVPIDX

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Other local variables

    INTEGER         IOS      ! tmp I/O status
    INTEGER         INVIOS   ! i/o status for MEG_REPINV_YN
    INTEGER         SPCIOS   ! i/o status for MEG_REPSPC_YN
    INTEGER         I, J     ! counter

    CHARACTER(5)    MRGSRC   ! value of MRG_SOURCE E.V.
    CHARACTER(5)    CTLMULT  ! value of MRG_CTLMAT_MULT E.V.
    CHARACTER(5)    CTLREAC  ! value of MRG_CTLMAT_REAC E.V.
    CHARACTER(5)    TMPBYDAY ! value of MRG_BYDAY E.V.
    CHARACTER(100)  BUFFER   ! text buffer
    CHARACTER(300)  MESG     ! message buffer

    CHARACTER(16) :: PROGNAME = 'GETMRGEV' ! program name

!***********************************************************************
!   begin body of subroutine GETMRGEV

!.........  Retrieve the environment variables that can be set to contain any
!           combination of "A" "B" "M" "P".
!.........  The default value of these will be blank, so that by not setting
!           the environment variable(s), nothing will happen with the program

    BUFFER = 'Control for which source categories are in merge'
    CALL ENVSTR( 'MRG_SOURCE', BUFFER, ' ', MRGSRC, IOS  )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_SOURCE"', 2 )
    END IF

    IF( MRGSRC .EQ. ' ' ) THEN
        MESG = 'No source categories selected! You must define ' //&
        &       'the MRG_SOURCE ' // CRLF() // BLANK10 //&
        &       'environment variable using the letters A, B, ' //&
        &       'M, and P'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    BUFFER = 'Setting for which source categories get ' //&
    &         'multiplicative control matrix'
    CALL ENVSTR( 'MRG_CTLMAT_MULT', BUFFER, ' ', CTLMULT, IOS  )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_CTLMAT_MULT"', 2 )
    END IF

    BUFFER = 'Setting for which source categories get ' //&
    &         'reactivity control matrix'
    CALL ENVSTR( 'MRG_CTLMAT_REAC', BUFFER, ' ', CTLREAC, IOS  )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_CTLMAT_REAC"', 2 )
    END IF

    BUFFER = 'Setting for which source categories get ' //&
    &         'by-day processing'
    CALL ENVSTR( 'MRG_BYDAY', BUFFER, ' ', TMPBYDAY, IOS  )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_BYDAY"', 2 )
    END IF

!.........  Process the character strings to set the various control flags
    AFLAG = ( INDEX( MRGSRC, 'A' ) .GT. 0 )
    BFLAG = ( INDEX( MRGSRC, 'B' ) .GT. 0 )
    MFLAG = ( INDEX( MRGSRC, 'M' ) .GT. 0 )
    PFLAG = ( INDEX( MRGSRC, 'P' ) .GT. 0 )

    AUFLAG = ( INDEX( CTLMULT, 'A' ) .GT. 0 .AND. AFLAG )
    MUFLAG = ( INDEX( CTLMULT, 'M' ) .GT. 0 .AND. MFLAG )
    PUFLAG = ( INDEX( CTLMULT, 'P' ) .GT. 0 .AND. PFLAG )

    ARFLAG = ( INDEX( CTLREAC, 'A' ) .GT. 0 .AND. AFLAG )
    MRFLAG = ( INDEX( CTLREAC, 'M' ) .GT. 0 .AND. MFLAG )
    PRFLAG = ( INDEX( CTLREAC, 'P' ) .GT. 0 .AND. PFLAG )

!.........  Determine if there are multiple source categories
    J = 0
    IF( AFLAG ) J = J + 1
    IF( BFLAG ) J = J + 1
    IF( MFLAG ) J = J + 1
    IF( PFLAG ) J = J + 1
    XFLAG = ( J .GT. 1 )

!.........  Retrieve the on/off environment variables

    TFLAG   = ENVYN( 'MRG_TEMPORAL_YN',&
    &                 'Use hourly emissions or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_TEMPORAL_YN"', 2 )
    END IF

    SFLAG   = ENVYN( 'MRG_SPCMAT_YN',&
    &                 'Use speciation matrices or not', .FALSE., IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_SPCMAT_YN"', 2 )
    END IF

    LMETCHK = ENVYN( 'MRG_METCHK_YN', 'Check consistency ' //&
    &                 'of met headers or not', .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_METCHK_YN"', 2 )
    END IF

    LMKTPON = ENVYN( 'MRG_MARKETPEN_YN', 'Use market penetration '//&
    &                 'from reactivity matrices', .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_MARKETPEN_YN"', 2 )
    END IF

    LGRDOUT = ENVYN( 'MRG_GRDOUT_YN', 'Output gridded I/O API ' //&
    &                 'NetCDF files or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_GRDOUT_YN"', 2 )
    END IF

    LREPSTA = ENVYN( 'MRG_REPSTA_YN', 'Output state total ' //&
    &                 'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPSTA_YN"', 2 )
    END IF

    LREPCNY = ENVYN( 'MRG_REPCNY_YN', 'Output county total ' //&
    &                 'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPCNY_YN"', 2 )
    END IF

    LREPINV = ENVYN( 'MRG_REPINV_YN', 'Report inventory ' //&
    &                 'emissions or not', .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPINV_YN"', 2 )
    END IF
    INVIOS = IOS

    LREPSPC = ENVYN( 'MRG_REPSPC_YN', 'Report speciated ' //&
    &                 'emissions or not', .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPSPC_YN"', 2 )
    END IF
    SPCIOS = IOS

    LREPCTL = ENVYN( 'MRG_REPCTL_YN', 'Report controlled ' //&
    &                 'emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPCTL_YN"', 2 )
    END IF

    LREPANY = ( LREPSTA .OR. LREPCNY )

!.........  Check for variable grid
    VARFLAG = ENVYN( 'USE_VARIABLE_GRID', 'Use variable grid ' //&
    &               'definition', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

!.........  Check if source grouping should be used
    SRCGRPFLAG = ENVYN( 'SMK_SRCGROUP_OUTPUT_YN', 'Use source ' //&
    &                    'grouping', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_SRCGROUP_OUTPUT_YN"', 2 )
    END IF

    IF( SRCGRPFLAG .AND. XFLAG ) THEN
        MESG = 'Source grouping cannot be used with multiple ' //&
        &       'source categories.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Check if source sub-grouping should be used
    SUBSECFLAG = ENVYN( 'SMK_SUB_SECTOR_OUTPUT_YN', 'Use sub-sector ' //&
    &                    'source grouping', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_SUB_SECTOR_OUTPUT_YN"', 2 )
    END IF

!.........  Error if both source group apportinment and sub-sector options
    IF( SRCGRPFLAG .AND. SUBSECFLAG ) THEN
        MESG = 'ERROR: Can NOT process both SMK_SRCGROUP_OUTPUT_YN ' //&
        &       ' & SMK_SUB_SECTOR_OUTPUT_YN flags at the same time'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Point-source specific environment variables
    IF ( PFLAG ) THEN

        MESG = 'Use layer fractions or not'
        LFLAG   = ENVYN( 'MRG_LAYERS_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_LAYERS_YN"', 2 )
        END IF

        MESG = 'Create CMAQ plume-in-grid outputs or not'
        I = ENVINT( 'SMK_PING_METHOD', MESG, 1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_PING_METHOD"', 2 )
        END IF
        PINGFLAG = ( I .EQ. 1 )

        MESG = 'Create CMAQ in-line point source outputs or not'
        I = ENVINT( 'SMK_ELEV_METHOD', MESG, 1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_ELEV_METHOD"', 2 )
        END IF
        INLINEFLAG = ( I .EQ. 2 )

        MESG = 'Processing CMAQ in-line point source ouptut by '//&
        &       'setting SMK_ELEV_METHOD to 2'
        IF( INLINEFLAG ) CALL M3MSG2( MESG )

        IF( LFLAG .AND. INLINEFLAG ) THEN
            MESG = 'ERROR: MUST set MRG_LAYERS_YN to N to '//&
            &   'output in-line CMAQ outputs ' // CRLF() // BLANK10
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        MESG = 'Create ASCII elevated sources file or not'
        ELEVFLAG = ENVYN( 'SMK_ASCIIELEV_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_ASCIIELEV_YN"', 2 )
        END IF

        IF( ELEVFLAG .AND. INLINEFLAG ) THEN
            MESG = 'ERROR: ASCII elevated output cannot be switched on'//&
            &       ' with in-line CMAQ outputs ' // CRLF() // BLANK10
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF


        IF ( ELEVFLAG ) PINGFLAG = .FALSE.

        MESG = 'Indicator for including explicit plume ' //&
        &       'rise sources'
        EXPLFLAG = ENVYN( 'EXPLICIT_PLUMES_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "EXPLICIT_PLUMES_YN"', 2 )
        END IF

!.............  Must be running for UAM-style processing to use explicit...
        IF( EXPLFLAG .AND. .NOT. ELEVFLAG ) THEN
            ELEVFLAG = .TRUE.
            MESG = 'NOTE: ASCII elevated output switched on to be'//&
            &       'consitent with ' // CRLF() // BLANK10//&
            &       'EXPLICIT_PLUMES_YN = Y setting.'
            CALL M3MSG2( MESG )
        END IF

!.............  Check source group restrictions
        IF ( SRCGRPFLAG ) THEN
            IF ( LFLAG ) THEN
                MESG = 'Source grouping cannot be used with layer ' //&
                &       'fractions.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF ( .NOT. INLINEFLAG ) THEN
                MESG = 'Source grouping must be used with in-line ' //&
                &       'CMAQ outputs.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

    END IF

!.........  If temporal processing, set which source categories get by-day
!           processing
    IF( TFLAG ) THEN
        AFLAG_BD = ( INDEX( TMPBYDAY, 'A' ) .GT. 0 .AND. AFLAG )
        MFLAG_BD = ( INDEX( TMPBYDAY, 'M' ) .GT. 0 .AND. MFLAG )
        PFLAG_BD = ( INDEX( TMPBYDAY, 'P' ) .GT. 0 .AND. PFLAG )
    END IF

!.........  Retrieve variable to indicate whether to use annual or average day data
    MESG = 'Use annual or average day emissions'
    LAVEDAY = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )

!.........  Set index for extracting pollutant data
    INVPIDX = 1
    IF( ( AFLAG .OR. PFLAG ) .AND.&
    &    .NOT. TFLAG          .AND.&
    &    LAVEDAY                    ) INVPIDX = 2

!.........  Check output flags to ensure at least some output
    IF( .NOT. LGRDOUT .AND.&
    &    .NOT. LREPSTA .AND.&
    &    .NOT. LREPCNY       ) THEN
        MESG = 'No output flags are set!  You must set at least ' //&
        &       'one of the' // CRLF() // BLANK10 //&
        &       'following to "Y": '//&
        &       'MRG_GRDOUT_YN, MRG_REPSTA_YN, or MRG_REPCNY_YN.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

!.........  Make biogenics flag consitent with speciation flag (must have
!           speciation to be able to include biogenic emissions)
    IF( BFLAG .AND. XFLAG .AND. .NOT. SFLAG ) THEN

        MESG = 'Speciation control environment variable ' //&
        &       '"MRG_SPCMAT_YN" indicates no'//&
        &       CRLF()// BLANK10// 'speciation, so biogenics will '//&
        &       'not be included despite setting of environment ' //&
        &       CRLF() // BLANK10 // 'variable "MRG_SOURCE"'
        CALL M3WARN( PROGNAME, 0, 0, MESG )

        BFLAG = .FALSE.

    ELSE IF ( BFLAG ) THEN
        SFLAG = .TRUE.

    ENDIF

!.........  Do not output ASCII elevated unless speciation is being done
    IF( ELEVFLAG .AND. .NOT. SFLAG ) THEN
        MESG = 'NOTE: Turning off ASCII elevated point source ' //&
        &       'outputs because speciation is turned off.'
        CALL M3MSG2( MESG )
        ELEVFLAG = .FALSE.

    END IF

!.........  Strange to have elevated ASCII and layer merge at the same time
    IF( ELEVFLAG .AND. LFLAG ) THEN
        MESG = 'WARNING: Elevated ASCII output and 3-d ' //&
        &       'output at the same time!'
        CALL M3MSG2( MESG )

    END IF

!.........  Cannot have layer merge at the same time as explicit plume rise
    IF( LFLAG .AND. EXPLFLAG ) THEN

        LFLAG = .FALSE.
        MESG = 'WARNING: Turning off layered merge (MRG_LAYERS_YN'//&
        &       ' = Y) because' // CRLF() // BLANK10 //&
        &       'explicit plume rise being used ' //&
        &       '(EXPLICIT_PLUMES_YN = Y).'
        CALL M3MSG2( MESG )

    END IF

!.........  If requested controlled reports, but no control matrices
!           provided, then override setting and give warning.
    IF( LREPCTL .AND. ( .NOT. AUFLAG .AND. .NOT. ARFLAG .AND.&
    &                    .NOT. MUFLAG .AND. .NOT. MRFLAG .AND.&
    &                    .NOT. PUFLAG .AND. .NOT. PRFLAG      )) THEN

        MESG = 'WARNING: Turning off requested controlled '//&
        &       'emissions reports, since no control matrices '//&
        &       'provided.'
        CALL M3MESG( MESG )

    END IF

!.........  Cannot have BFLAG without TFLAG
    IF( BFLAG ) TFLAG = .TRUE.

!.........  Report that flags for reporting inventory emissions, speciated
!           emissions or not, and controlled emissions or not do not work yet
    IF( INVIOS .GE. 0 ) THEN   ! If it was set to something
        MESG = 'NOTE: MRG_REPINV_YN control is not yet functional.'
        CALL M3MSG2( MESG )
    END IF

    IF( SPCIOS .GE. 0 ) THEN
        MESG = 'NOTE: MRG_REPSPC_YN control is not yet functional.'
        CALL M3MSG2( MESG )
    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE GETMRGEV
