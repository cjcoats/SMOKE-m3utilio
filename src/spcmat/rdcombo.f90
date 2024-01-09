
SUBROUTINE RDCOMBO( CDEV, ENAM )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !    This program
    !
    !  REVISION  HISTORY:
    !    Created 12/07 by M. Houyoux
    !
    !    Version 8/2016 by C.Coats:
    !    Check whether weights sum to 1
    !    Needed status-checks 
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90", and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !    System
    ! File: @(#)$Id$
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module contains the speciation profile tables
    USE MODSPRO, ONLY: CMBNP, CMBSPCD, CMBWGHT, CMBMAX

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT    (IN) :: CDEV        ! unit number of input file
    CHARACTER(*), INTENT    (IN) :: ENAM        ! pol/emis type name of interest

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !.......  Local parameters
    INTEGER      , PARAMETER :: DEFLTTYP = 0
    INTEGER      , PARAMETER :: STATETYP = 1
    INTEGER      , PARAMETER :: CNTYTYP  = 2
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDCOMBO'     ! program name

    !.......  Local arrays
    INTEGER       :: CMBTYP( NINVIFIP )      ! type of last entry applied to county (state- or county-specific)
    INTEGER, SAVE :: ERRCNT( 6 )
    INTEGER, SAVE :: WARNCNT( 2 )
    REAL          :: CWEIGHT( CMBMAX )              ! tmp array for profile weights
    CHARACTER(NAMLEN3) :: SEGMENT( 24 )      ! Segments of parsed lines

    !.......  Other local variables
    INTEGER         F, I, N           !  counters and indices

    INTEGER         FCNT        !  tmp count of FIPs codes missing records in GSPRO_COMBO
    INTEGER         IOS         !  i/o status
    INTEGER         IREC        !  record counter
    INTEGER         NLINES      !  number of lines
    INTEGER, SAVE ::PERIOD      !  period to use to match to GSPRO_COMBO file
    INTEGER         PER         !  tmp period from input file
    INTEGER, SAVE ::MXERR                    ! max no. errors
    INTEGER, SAVE ::MXWARN                   ! max no. warnings
    INTEGER         NP          !  tmp number of profiles per combo
    INTEGER         NRECS       !  number of records read in current call

    LOGICAL          TFLAG                  ! out-of-bounds
    LOGICAL          EFLAG                  ! true: error detected
    LOGICAL, SAVE :: FIRSTIME = .TRUE.
    LOGICAL, SAVE :: FIRSTSTA = .TRUE.
    LOGICAL, SAVE :: CMBCHECK

    CHARACTER(FIPLEN3), SAVE ::  FIPZERO      ! zero Cy/St/Co code

    CHARACTER(FIPLEN3) CFIP         !  tmp buffer for state/county FIPS code
    CHARACTER(STALEN3) CSTA         !  tmp buff for state FIPS code
    CHARACTER(STALEN3) PSTA         !  tmp buff for previous state FIPS code
    CHARACTER(NAMLEN3) CPOL         !  tmp buffer for pollutant code
    CHARACTER(256)     MESG         !  message buffer
    CHARACTER(512)     LINE         !  line buffer

    !***********************************************************************
    !   begin body of subroutine RDCOMBO

    EFLAG    = .FALSE.
    !.......  Perform one-time steps
    IF ( FIRSTIME ) THEN

        MXWARN = ENVINT( WARNSET , ' ', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'Bad environment variable "' // WARNSET // '"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        MXERR  = ENVINT( ERRSET  , ' ', 100, IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'Bad environment variable "' // ERRSET // '"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        CMBCHECK = ENVYN( 'COMBO_CHKFRACS',    &
                   'FATAL for bad sum-of-weights for profile-combos?',    &
                   .TRUE., IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'Bad environment variable "COMBO_CHKSFRACS"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF

        MESG   = 'Period to read from GSPRO_COMBO file'
        PERIOD = ENVINT( 'SPCMAT_PERIOD', MESG, 1, IOS )
        IF ( IOS .GT. 0 ) THEN
            MESG = 'Bad environment variable "SPCMAT_PERIOD"'
            CALL M3EXIT( PROGNAME, 0,0, MESG, 2 )
        END IF
        !.......  Set up zero strings for FIPS code

        FIPZERO  = REPEAT( '0', FIPLEN3 )

        ALLOCATE( CMBNP( NINVIFIP ),            &
                CMBSPCD( NINVIFIP, CMBMAX ),    &
                CMBWGHT( NINVIFIP, CMBMAX ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CMBWGHT', PROGNAME )

        FIRSTIME = .FALSE.

    END IF

    !.......  Initialize public and local arrays to 0

    CMBNP   = 0       ! array
    CMBSPCD = ' '     ! array
    CMBWGHT = 0.      ! array
    CMBTYP  = -1       ! array

    !.......  Other initializations
    PSTA = '-9'

    !.......  Get the number of lines in the input file
    NLINES = GETFLINE( CDEV, 'GSPRO Combos file' )

    !.......  Loop through file and read until the end.
    IREC  = 0
    NRECS = 0
    TFLAG = .FALSE.

    DO I = 1, NLINES

        READ ( CDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,    &
                'reading profile combinations file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.......  Split out columns from line
        CALL PARSLINE( LINE, 24, SEGMENT )

        !.......  Get period, pollutant, and number of profiles from line
        CPOL = SEGMENT( 1 )
        CFIP = SEGMENT( 2 )
        CALL PADZERO( CFIP )
        PER  = STR2INT( SEGMENT( 3 ) )
        NP   = STR2INT( SEGMENT( 4 ) )

        !.......  If pollutant doesn't match, or if number
        !    of profiles is <=0, then skip
        IF ( CPOL .NE. ENAM .OR. NP <= 0 ) CYCLE

        !.......  If period is greater than 0 and doesn't match, then skip
        !.......  This intentionally allows all entries with PERIOD=0 to apply
        !    to all periods (if CPOL and NP criteria are met)
        IF ( PER > 0 .AND. PER .NE. PERIOD ) CYCLE

        !.......  Error message if a number of profile is greater than 10
        IF ( NP > 10 ) THEN
            WRITE( MESG,94010 ) 'ERROR: A number of COMBO '//    &
                   'profiles for county '//CFIP//' is greater '//    &
                   'than a max. 10 at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.......  If entry assigns profile to all FIPs codes, then implement it
        IF ( CFIP == FIPZERO ) THEN

            DO F = 1, NINVIFIP           ! Loop through inventory state/county FIPS codes

                !.......  Give error if any state- or county-specific entries already
                !    applied in the input file
                IF( CMBTYP( F ) >= STATETYP .AND.    &
                    ERRCNT(1) <=  MXERR ) THEN

                    IF( CMBTYP(F) == STATETYP ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Default' //    &
                          'record at line',IREC,'comes after '//    &
                          'state-specific record matching same '//    &
                          'source.'
                    ELSE IF ( CMBTYP(F) == CNTYTYP ) THEN
                        WRITE( MESG,94010 ) 'ERROR: Default ' //    &
                          'record at line',IREC,'comes after '//    &
                          'county-specific record matching same '//    &
                          'source.'
                    END IF

                    CALL M3MESG( MESG )
                    ERRCNT(1) = ERRCNT(1) + 1
                    EFLAG = .TRUE.
                    CYCLE

                END IF

                !.......  Give error if duplicate (non-zero) entry
                IF ( CMBTYP( F ) == DEFLTTYP .AND.    &
                     ERRCNT(2) <= MXERR        ) THEN

                    WRITE( MESG,94010 )    &
                      'ERROR: Duplicate default (FIPS=0) '//    &
                      'entry found at line',IREC,    &
                      'of GSPRO_COMBO file'
                    CALL M3MESG( MESG )
                    ERRCNT(2) = ERRCNT(2) + 1
                    EFLAG = .TRUE.
                    CYCLE

                END IF

                !.......  Record  flag for default entry (FIPS=00000)
                CMBTYP( F ) = DEFLTTYP

                !.......  Store entry for current FIPs, pollutant, period
                CMBNP( F ) = NP

                !.......  Convert fractions from strings to reals
                DO N = 1,  NP
                    CWEIGHT( N ) = STR2REAL(SEGMENT(4+N*2))
                END DO

                !.......  Check if profile fractions meet the +/- 0.001 criterion
                !    and renormalize if needed. Provide a warning if need
                !    to renormalize.
                CALL CHECK_AND_SET_FRACS

                DO N = 1, NP

                    CMBSPCD( F,N )= SEGMENT( 3+N*2 )
                    CMBWGHT( F,N )= CWEIGHT( N )

                END DO

            END DO          ! loop of inventory FIPs codes

        !.......  If state-specific entry (county set to 000), then
        !    loop through all FIPs codes and apply information.
        ELSE IF ( CFIP(STALEN3+1:FIPLEN3) == '000' ) THEN

            !.......  Extract state code from state/county FIPS code
            CSTA = CFIP(1:STALEN3)

            FIRSTSTA = .TRUE.
            DO F = 1, NINVIFIP               ! Loop through inventory state/county FIPS codes

                !.......  If inventory state code matches COMBO record state,
                !    the store information for current county
                IF( INVCFIP(F)(1:STALEN3) .EQ. CSTA ) THEN

                    !.......  Give error if any county-specific entries already
                    !    applied in the input file
                    IF( CMBTYP( F ) == CNTYTYP .AND.    &
                        ERRCNT(1) <=  MXERR ) THEN

                        WRITE( MESG,94010 ) 'ERROR: State-specific '    &
                          //'record at line',IREC,'comes after '//      &
                          'county-specific record matching same '//     &
                          'source.'
                        CALL M3MESG( MESG )
                        ERRCNT(1) = ERRCNT(1) + 1
                        EFLAG = .TRUE.
                        CYCLE

                    END IF

                    !.......  Give error if duplicate (non-zero) entry
                    IF ( CMBTYP( F ) == STATETYP .AND.    &
                         CSTA == PSTA .AND. FIRSTSTA .AND.    &
                         ERRCNT(2) <= MXERR        ) THEN

                        WRITE( MESG,94010 )    &
                          'ERROR: Duplicate entry found at line',    &
                          IREC,'for state "' // CSTA //    &
                          '000" of GSPRO_COMBO file'
                        CALL M3MESG( MESG )
                        ERRCNT(2) = ERRCNT(2) + 1
                        FIRSTSTA = .FALSE.
                        EFLAG = .TRUE.
                        CYCLE

                    END IF

                    !.......  Record  flag for state-specific entry
                    CMBTYP( F ) = STATETYP

                    !.......  Store entry for current FIPs, pollutant, period
                    CMBNP( F ) = NP

                    !.......  Convert fractions from strings to reals
                    DO N = 1,  NP
                        CWEIGHT( N ) = STR2REAL(SEGMENT(4+N*2))
                    END DO

                    !.......  Check if profile fractions meet the +/- 0.001 criterion
                    !    and renormalize if needed. Provide a warning if need
                    !    to renormalize.
                    CALL CHECK_AND_SET_FRACS

                    DO N = 1, NP

                        CMBSPCD( F,N )= SEGMENT( 3+N*2 )
                        CMBWGHT( F,N )= CWEIGHT( N )

                    END DO

                END IF


                !.......  Set previous state code for next iteration
                PSTA = CSTA

            END DO

        !.......  Otherwise, apply entry to a single county
        ELSE

            !.......  Search for state/county FIPS code in list of these from
            !    the inventory. If not found, then skip this line.
            F = FINDC( CFIP, NINVIFIP, INVCFIP )
            IF ( F .LE. 0 ) CYCLE

            !.......  Give error if duplicate (non-zero) entry
            IF ( CMBTYP( F ) == CNTYTYP .AND.    &
                 ERRCNT(3)   <= MXERR         ) THEN
                WRITE( MESG,94010 )    &
                  'ERROR: Duplicate entry found at line', IREC,    &
                  'for county "'//CFIP//'" of GSPRO_COMBO file'
                CALL M3MESG( MESG )
                ERRCNT(3) = ERRCNT(3) + 1
                EFLAG = .TRUE.
                CYCLE
            END IF

            !.......  Record  flag for county-specific entry
            CMBTYP( F ) = CNTYTYP

            !.......  Store entry for current FIPs, pollutant, period
            CMBNP( F ) = NP

            !.......  Convert fractions from strings to reals
            DO N = 1,  NP
                CWEIGHT( N ) = STR2REAL( SEGMENT( 4+N*2 ) )
            END DO

            !.......  Check if profile fractions meet the +/- 0.001 criterion
            !    and renormalize if needed. Provide a warning if need
            !    to renormalize.
            CALL CHECK_AND_SET_FRACS

            DO N = 1, NP

                CMBSPCD( F,N ) = SEGMENT( 3 + N*2 )
                CMBWGHT( F,N ) = CWEIGHT( N )

            END DO

        END IF

        NRECS = NRECS + 1

    END DO

    !.......  Check that at least one record present for the selected PERIOD
    !.......  Note: want this always to be reported, regardless of MAXERROR
    IF( NRECS == 0 ) THEN

        WRITE( MESG, 94010 ) 'ERROR: No records found for period',    &
               PERIOD, 'and pollutant '//TRIM( ENAM ) // ' in '//    &
               'GSPRO_COMBO file.'
        CALL M3MESG( MESG )
        EFLAG = .TRUE.

    !.......  Record to log how many counties received the COMBO approach
    ELSE

        FCNT = 0
        DO F = 1, NINVIFIP

        !.......  Count counties found
            IF ( CMBTYP( F ) .GT. -1 ) FCNT = FCNT + 1

        END DO

        !.......  Write error if nothing matched (there should have been
        !    at least one thing or else this routine wouldn't have
        !    been called)
        !.......  Note: want this always to be reported, regardless of MAXERROR
        IF ( FCNT == 0 ) THEN
            WRITE( MESG,94010 ) 'ERROR: No counties assigned COMBO profiles'//  &
               CRLF()//BLANK10// 'for period', PERIOD,'and pollutant "'// TRIM( ENAM ) //'"'
            EFLAG = .TRUE.

        !.......  Write message with total number of counties receiving approach
        !.......  Note: want this always to be reported, regardless of MAXERROR
        ELSE
            WRITE( MESG,94010 ) 'NOTE:',FCNT, 'counties assigned COMBO profiles'//  &
               CRLF()//BLANK10// 'for period', PERIOD,'and pollutant "'// TRIM( ENAM ) //'"'
            CALL M3MESG( MESG )
        END IF

    END IF

    !.......  If error found, abort program
    IF( EFLAG ) THEN
        MESG = 'Problem reading or applying profile combos data from GSPRO_COMBO'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE IF ( TFLAG .AND. CMBCHECK ) THEN
        MESG = 'Non-normalized profile combos in data from GSPRO_COMBO'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !.......  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //    &
           'Check format of GSPRO_COMBO file.'

    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subroutine writes the message when a default
    !    speciation profile is unavailable for a given pollutant
    SUBROUTINE CHECK_AND_SET_FRACS

        !....... Inherits the following variables:
        !    CWEIGHT
        !    NP
        !    IREC

        !.......  Local parameters
        REAL, PARAMETER :: TOLERANCE = 0.001

        !.......  Local variables
        INTEGER   N
        REAL      SUM            ! sum of weights
        REAL      SUMINV         ! inverse of SUM

        !----------------------------------------------------------------------
        SUM = 0.0
        DO N = 1,  NP
            SUM = SUM + CWEIGHT( N )
        END DO

        IF( SUM > 1.0+TOLERANCE .OR. SUM < 1.0-TOLERANCE ) THEN

            !.......  Non-normalization is an error:
            IF ( CMBCHECK ) THEN
                WRITE( MESG,94010 ) 'WARNING: GSPRO_COMBO '//    &
                   'Non-normalized fractions at line',IREC
                CALL M3MESG( MESG )
                TFLAG = .TRUE.
                RETURN
            END IF

            !.......  Renormalize fractions
            !.......  Send warning message that entry is being renormalized
            SUMINV = 1.0 / SUM
            DO N = 1, NP
                CWEIGHT( N ) = CWEIGHT( N ) * SUMINV
            END DO

            IF( SUM > 1+TOLERANCE .AND. WARNCNT(1) <= MXWARN ) THEN
                WRITE( MESG,94010 ) 'WARNING: GSPRO_COMBO '//    &
                   'fractions summed > 1.001 at line',IREC,    &
                   'and were renormalized'
                CALL M3MESG( MESG )
                WARNCNT(1) = WARNCNT(1) + 1
                TFLAG = .TRUE.

            ELSE IF ( SUM < 1-TOLERANCE .AND.    &
                      WARNCNT(2) <= MXWARN    ) THEN
                WRITE( MESG,94010 ) 'WARNING: GSPRO_COMBO '//    &
                   'fractions summed < 0.999 at line',IREC,    &
                   'and were renormalized'
                CALL M3MESG( MESG )
                WARNCNT(2) = WARNCNT(2) + 1
                TFLAG = .TRUE.
            END IF

        END IF

        !------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE CHECK_AND_SET_FRACS

END SUBROUTINE RDCOMBO
