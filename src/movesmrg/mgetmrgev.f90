
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
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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

    !.......  MODULES for public variables
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    !.......  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,     &
                        MFLAG_BD, VARFLAG, LREPANY,     &
                        LMETCHK, LGRDOUT, LREPCNY,      &
                        LREPSTA, LREPSCC, LREPSRC,      &
                        SRCGRPFLAG, SUBSECFLAG, SMATCHK

    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: RPDFLAG, RPHFLAG, ONIFLAG, RPVFLAG, RPPFLAG,   &
                         MVFILDIR, TVARNAME, SPDPROFLAG, SPDISTFLAG,    &
                         CFFLAG, EXPCFFLAG, REFCFFLAG, TEMPBIN,         &
                         MOPTIMIZE, GRDENV, TOTENV, MTMP_OUT,           &
                         NOXADJFLAG, NOXADJEQS, RPSFLAG, ETABLEFLAG,    &
                         MXTEMP, MNTEMP, TMPINC, NTBINS

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   Other local variables

    INTEGER         IOS          ! tmp I/O status
    INTEGER         I, J, L      ! counter

    CHARACTER(5)    TMPBYDAY     ! value of MRG_BYDAY E.V.
    CHARACTER(100)  BUFFER       ! text buffer
    CHARACTER(300)  MESG         ! message buffer

    CHARACTER(16) , PARAMETER:: PROGNAME = 'GETMRGEV'     ! program name

    !***********************************************************************
    !   begin body of subroutine GETMRGEV

    !.......  Set source category based on environment variable setting
    CALL GETCTGRY

    !.......  End program if source category is not mobile sources
    IF( CATEGORY /= 'MOBILE' ) THEN
        MESG = 'Program ' // TRIM( PROGNAME ) // ' does not support ' //&
               TRIM( CATEGORY ) // ' sources.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Retrieve the environment variables that can be set to contain any
    !         combination of "A" "B" "M" "P".
    !.......  The default value of these will be blank, so that by not setting
    !         the environment variable(s), nothing will happen with the program

    BUFFER = 'Setting for which source categories get by-day processing'
    CALL ENVSTR( 'MRG_BYDAY', BUFFER, ' ', TMPBYDAY, IOS  )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_BYDAY"', 2 )
    END IF

    !.......  Set merge source category flags (impacts library routines shared
    !         between Movesmrg and Smkmerge)
    AFLAG = .FALSE.
    BFLAG = .FALSE.
    MFLAG = .TRUE.
    PFLAG = .FALSE.

    !.......  Retrieve the on/off environment variables
    SMATCHK = ENVYN( 'USE_SPCMAT_SPC_YN', 'Use Speciation Matrix '//&
                     'output file to calculate model species',&
                     .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_SPCMAT_SPC_YN"', 2 )
    END IF

    LMETCHK = ENVYN( 'MRG_METCHK_YN', 'Check consistency ' //&
                     'of met headers or not', .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_METCHK_YN"', 2 )
    END IF

    LGRDOUT = ENVYN( 'MRG_GRDOUT_YN', 'Output gridded I/O API ' //&
                     'files or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_GRDOUT_YN"', 2 )
    END IF

    LREPSTA = ENVYN( 'MRG_REPSTA_YN', 'Output state total ' //&
                     'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPSTA_YN"', 2 )
    END IF

    LREPCNY = ENVYN( 'MRG_REPCNY_YN', 'Output county total ' //&
                     'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPCNY_YN"', 2 )
    END IF

    LREPSCC = ENVYN( 'MRG_REPSCC_YN', 'Output SCC total ' //&
                     'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPSCC_YN"', 2 )
    END IF

    LREPSRC = ENVYN( 'MRG_REPSRC_YN', 'Output source total ' //&
                     'ASCII reports or not', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_REPSRC_YN"', 2 )
    END IF

    TEMPBIN = ENVREAL( 'TEMP_BUFFER_BIN', 'Buffer for ' //&
                     ' min/max temperature ', 0., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TEMP_BUFFER_BIN"', 2 )
    END IF

    LREPANY = ( LREPSTA .OR. LREPCNY .OR. LREPSCC .OR. LREPSRC )

    !.......  Retrieve variables for setting the output units for gridded and
    !         country/state/county total data
    BUFFER = 'Units for output gridded emissions'
    CALL ENVSTR( 'MRG_GRDOUT_UNIT', BUFFER, ' ', GRDENV, IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_GRDOUT_UNIT"', 2 )
    END IF

    BUFFER = 'Units for output state/county total emissions'
    CALL ENVSTR( 'MRG_TOTOUT_UNIT', BUFFER, ' ', TOTENV, IOS)
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MRG_TOTOUT_UNIT"', 2 )
    END IF

    IF( LREPSRC ) THEN
        LREPCNY = .FALSE.
        LREPSTA = .FALSE.
        LREPSCC = .FALSE.
        MESG = 'WARNING: Reset MRG_REPSCC_YN, MRG_REPCNY, and ' //&
             'MRG_REPSTA_YN to N since MRG_REPSRC_YN is set to Y '
        CALL M3MESG( MESG )
    END IF

    !.......  Check for variable grid
    VARFLAG = ENVYN( 'USE_VARIABLE_GRID', 'Use variable grid definition', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_VARIABLE_GRID"', 2 )
    END IF

    !.......  Check if source grouping should be used
    SRCGRPFLAG = ENVYN( 'SMK_SRCGROUP_OUTPUT_YN', 'Use source grouping', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_SRCGROUP_OUTPUT_YN"', 2 )
    END IF

    !.......  Check if source sub-grouping should be used
    SUBSECFLAG = ENVYN( 'SMK_SUB_SECTOR_OUTPUT_YN', 'Use sub-sector ' //&
                        'source grouping', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_SUB_SECTOR_OUTPUT_YN"', 2 )
    END IF

    !.......  Error if both source group apportinment and sub-sector options
    IF( SRCGRPFLAG .AND. SUBSECFLAG ) THEN
        MESG = 'ERROR: Can NOT process both SMK_SRCGROUP_OUTPUT_YN ' //&
               '  SMK_SUB_SECTOR_OUTPUT_YN flags at the same time'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Check for rate-per-distance, rate-per-vehicle, or rate-per-profile processing
    RPDFLAG = ENVYN( 'RPD_MODE', 'Calculate rate-per-distance emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RPD_MODE"', 2 )
    END IF

    RPHFLAG = ENVYN( 'RPH_MODE', 'Calculate rate-per-hour emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RPH_MODE"', 2 )
    END IF

    RPVFLAG = ENVYN( 'RPV_MODE', 'Calculate rate-per-vehicle emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RPV_MODE"', 2 )
    END IF

    RPPFLAG = ENVYN( 'RPP_MODE', 'Calculate rate-per-profile  emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RPP_MODE"', 2 )
    END IF

    RPSFLAG = ENVYN( 'RPS_MODE', 'Calculate rate-per-start emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RPS_MODE"', 2 )
    END IF

    ONIFLAG = ENVYN( 'ONI_MODE', 'Calculate Off-Network Idling emissions', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ONI_MODE"', 2 )
    END IF

    MOPTIMIZE = ENVYN( 'MEMORY_OPTIMIZE_YN', 'Optimize Memory usage', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MEMORY_OPTIMIZE_YN"', 2 )
    END IF

    MTMP_OUT = ENVYN( 'MTMP_OUTPUT_YN', 'Output mobile hourly emissions ', .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MTMP_OUTPUT_YN"', 2 )
    END IF

    MESG = 'Apply control factors to emissions'
    CFFLAG = ENVYN( 'USE_CONTROL_FACTORS', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_CONTROL_FACTORS"', 2 )
    END IF

    IF( CFFLAG ) THEN
        MESG = 'Use pollutant/species-specific control factor'
        EXPCFFLAG = ENVYN( 'USE_EXP_CONTROL_FAC_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_EXP_CONTROL_FAC_YN"', 2 )
        END IF

        MESG = 'Use reference county-specific control factor'
        REFCFFLAG = ENVYN( 'USE_REF_CONTROL_FAC_YN', MESG, .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_REF_CONTROL_FAC_YN"', 2 )
        END IF
    END IF

    IF( .NOT. RPDFLAG .AND.     &
        .NOT. RPSFLAG .AND.     &
        .NOT. ONIFLAG .AND.     &
        .NOT. RPHFLAG .AND.     &
        .NOT. RPVFLAG .AND.     &
        .NOT. RPVFLAG .AND.     &
        .NOT. RPPFLAG ) THEN
        MESG = 'No mode selected    !  You must set either RPD_MODE, ' //&
               'RPH_MODE, RPV_MODE, RPS_MODE, RPP_MODE or ONI_MODE to "Y".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Select default processing mode
    IF( RPDFLAG ) THEN
        RPHFLAG = .FALSE.
        RPVFLAG = .FALSE.
        RPSFLAG = .FALSE.
        RPPFLAG = .FALSE.
        ONIFLAG = .FALSE.
    ELSE IF( RPVFLAG ) THEN
        RPDFLAG = .FALSE.
        RPSFLAG = .FALSE.
        RPHFLAG = .FALSE.
        RPPFLAG = .FALSE.
        ONIFLAG = .FALSE.
    ELSE IF( RPPFLAG ) THEN
        RPDFLAG = .FALSE.
        RPSFLAG = .FALSE.
        RPHFLAG = .FALSE.
        RPVFLAG = .FALSE.
        ONIFLAG = .FALSE.
    ELSE IF( RPHFLAG ) THEN
        RPDFLAG = .FALSE.
        RPSFLAG = .FALSE.
        RPVFLAG = .FALSE.
        RPPFLAG = .FALSE.
        ONIFLAG = .FALSE.
    ELSE IF( RPSFLAG ) THEN
        RPDFLAG = .FALSE.
        RPVFLAG = .FALSE.
        RPHFLAG = .FALSE.
        RPPFLAG = .FALSE.
        ONIFLAG = .FALSE.
    END IF

    !.......  Check if hourly speeds should be used
    IF( RPDFLAG ) THEN
        SPDPROFLAG = ENVYN( 'USE_HOURLY_SPEEDS',    &
                            'Use hourly speed profile data', .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_HOURLY_SPEEDS"', 2 )
        END IF

        SPDISTFLAG = ENVYN( 'USE_AVG_SPD_DIST',    &
                            'Use average speed distribution profile data', .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_AVG_SPD_DIST"', 2 )
        END IF
    END IF

    !.......  Get if NOx adjustment should be applied
    IF( RPDFLAG .OR. RPHFLAG .OR. ONIFLAG .OR. RPSFLAG ) THEN
        NOXADJFLAG = ENVYN( 'APPLY_NOX_HUMIDITY_ADJ',    &
                            'Apply humidity adjusment to NOx emissions', .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "APPLY_NOX_HUMIDITY_ADJ"', 2 )
        END IF
    END IF

    !......  Define the type of NOx adj eqs (MOVES 3 or older version)
    IF( NOXADJFLAG ) THEN
        NOXADJEQS = .FALSE.         ! True: Use the latest MOVES4 NOx adj eqs
        NOXADJEQS = ENVYN( 'USE_MOVES4_NOX_ADJ_EQS',    &
                           'Use the MOVES4 NOx humidity correction equations',&
                           .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "USE_MOVES4_NOX_ADJ_EQS"', 2 )
        END IF
    END IF
    !.......  Output precomputed gridded hourly emissions by temp. bins
    MESG = 'Output temperature-bin precomputed gridded hourly emissions tables or not'
    ETABLEFLAG = ENVYN( 'OUTPUT_EMIS_TABLE_YN', MESG,  .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTPUT_EMIS_TABLE_YN"', 2 )
    END IF

    IF( RPPFLAG .AND. ETABLEFLAG ) THEN
        MESG = 'WARNING: OUTPUT_EMIS_TABLE_YN flag has been reset to N ' //&
               'since it is not implemented for RPP mode run yet'
        CALL M3MESG( MESG )
        ETABLEFLAG = .FALSE.
    END IF

    !.......   Disabling several options to generate the emissisons table
    IF( ETABLEFLAG ) THEN

        MESG = 'Lowest temperature for precomputed emissions table output file'
        MNTEMP = ENVREAL( 'MIN_TEMP_EMIS_TABLE', MESG, 0.0, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MIN_TEMP_EMIS_TABLE"', 2 )
        END IF

        MESG = 'Highest temperature for precomputed emissions table output file'
        MXTEMP = ENVREAL( 'MAX_TEMP_EMIS_TABLE', MESG, 120.0, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MAX_TEMP_EMIS_TABLE"', 2 )
        END IF

        MESG = 'Temperature bin increment for precomputed emissions table output file'
        TMPINC = ENVREAL( 'TEMP_INCREMENT_EMIS_TABLE', MESG, 5., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TEMP_INCREMENT_EMIS_TABLE"', 2 )
        END IF

        NTBINS = INT( ( MXTEMP - MNTEMP ) / TMPINC ) + 1      ! no of temperature bins

        IF( SRCGRPFLAG .OR. SUBSECFLAG ) THEN
            MESG = 'WARNING: SMK_SRCGROUP_OUTPUT_YN and SMK_SUB_SECTOR_OUTPUT_YN' //&
                   ' flags have been reset to N'
            CALL M3MSG2( MESG )
            SRCGRPFLAG = .FALSE.        ! reset SMK_SRCGROUP_OUTPUT_YN=N
            SUBSECFLAG = .FALSE.        ! reset SMK_SUB_SECTOR_OUTPUT_YN=N
        END IF

        IF( MTMP_OUT ) THEN
            MESG = 'WARNING: MTMP_OUTPUT_YN flag has been reset to N'
            CALL M3MSG2( MESG )
            MTMP_OUT = .FALSE.
        END IF

        IF( NOXADJFLAG ) THEN
            MESG = 'WARNING: APPLY_NOX_HUMIDITY_ADJ flag has been reset to N'
            CALL M3MSG2( MESG )
            NOXADJFLAG = .FALSE.
            NOXADJEQS  = .FALSE.
        END IF

        IF( LREPANY ) THEN
            MESG = 'WARNING: Disabling any summary reports generation'
            CALL M3MSG2( MESG )
            LREPANY = .FALSE.
            LREPSTA = .FALSE.
            LREPCNY = .FALSE.
            LREPSCC = .FALSE.
            LREPSRC = .FALSE.
        END IF

    END IF

    !.......  Get directory where MOVES output files are stored
    BUFFER = 'Location of MOVES output files'
    CALL ENVSTR( 'SMK_MVSPATH', BUFFER, '.', MVFILDIR, IOS )

    IF( IOS /= 0 ) THEN
        MESG = 'ERROR: Could not read value for SMK_MVSPATH.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  If temporal processing, set which source categories get by-day
    !         processing
    MFLAG_BD = ( INDEX( TMPBYDAY, 'M' ) .GT. 0 )

    !.......  Get name of temperature variable to read from meteorology files
    BUFFER = 'Name of temperature variable'
    CALL ENVSTR( 'TVARNAME', BUFFER, 'TEMP2', TVARNAME, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "TVARNAME"', 2 )
    END IF

    !.......  Check output flags to ensure at least some output
    IF( .NOT. LGRDOUT .AND.     &
        .NOT. LREPSTA .AND.     &
        .NOT. LREPCNY .AND.     &
        .NOT. LREPSRC .AND.     &
        .NOT. LREPSCC       ) THEN
        MESG = 'No output flags are set    !  You must set at least one of the' // &
               CRLF() // BLANK10 // 'following to "Y": '//&
               'MRG_GRDOUT_YN, MRG_REPSTA_YN, MRG_REPCNY_YN, ' //&
               'MRG_REPSCC_YN, or MRG_REPSRC_YN.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.......94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE GETMRGEV
