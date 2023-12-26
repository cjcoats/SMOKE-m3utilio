
PROGRAM GCNTL4CARB

    !***********************************************************************
    !  program body starts at line
    !
    !  DESCRIPTION:
    !       This program will develop GCNTL input files for Cntlmat program in SMOKEv4.0
    !       based on California Air Resource Board (CARB) specific control adjustment
    !       excel spreadsheet
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created  Apr 2014 by B.H. Baek
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: $Id$
    !
    ! COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
    ! All Rights Reserved
    !
    ! Institute for the Environment
    ! University of North Carolina at Chapel Hill
    ! Chapel Hill, NC 27599
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !*************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......  This module contains the information about the source category

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'         ! emissions constant parameters

    !.......   PARAMETERS and their descriptions:

    CHARACTER(50), PARAMETER :: SCCSW    = '%W%'
    CHARACTER(16), PARAMETER :: PROGNAME = 'GCNTL4CARB'       !  program name

    !....... EXTERNAL FUNCTIONS and their descriptions:
    INTEGER,EXTERNAL :: GETFLINE
    LOGICAL,EXTERNAL :: BLKORCMT
    LOGICAL,EXTERNAL :: CHKINT

    !.......  local arrays
    CHARACTER( SCCLEN3 ), ALLOCATABLE :: FIPSLIST( : )          ! state codes

    !.......  Logical file names and unit numbers
    INTEGER   CDEV             ! control adjustment factor input file
    INTEGER   FDEV             ! co-ab-dis description input file
    INTEGER   ODEV             ! Output file
    INTEGER   LDEV             ! log file unit number

    !.......  Local arrays
    CHARACTER( 20 ) SEGMENT( 20 )

    !.......  Local variables
    INTEGER   I, CNT, F, J, PF, ST             ! indices and counters

    INTEGER   IOS              ! i/o status
    INTEGER   IREC             ! record counter
    INTEGER   NLINES           ! number of entries in split file
    INTEGER   NFIPS, NEIC      ! no of master COABDIS and EIC list
    INTEGER   REGID            ! regional flag ID (2, 4, 5, 6)
    INTEGER   YEAR, YR         ! modeling year
    INTEGER   FIPS, CNTY       ! tmp county code (6-digit)

    REAL      FACS             ! real control factor

    LOGICAL :: EFLAG = .FALSE.       ! true: error found

    CHARACTER( 256 )  MESG        ! temporary message array
    CHARACTER( 256 )  LINE        ! line buffer
    CHARACTER( 100 )  CBUF        ! tmp output buffer

    CHARACTER( CNYLEN3 )  ARBN, DSTR, SUBC
    CHARACTER( FIPLEN3 )  GEOCODE        ! tmp 12-digit co/ab/dis code
    CHARACTER( SCCLEN3 )  EIC            ! EIC
    CHARACTER( NAMLEN3 )  POL            ! Pollutant name

    !***********************************************************************
    !   begin body of program GCNTL4CARB

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, SCCSW, PROGNAME )

    !.......  Prompt for name of input files
    MESG = 'Enter logical name of the CO-AB-DIS description input file'
    FDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'COABDIS', PROGNAME )

    !.......  Open CARB control factors input file
    MESG = 'Enter logical name of the control adjustment input file'
    CDEV = PROMPTFFILE( MESG, .TRUE., .TRUE., 'INFILE', PROGNAME )

    !.......   Get size for co-ab-dis description file
    NLINES = GETFLINE( FDEV, 'CO-AB-DIS Definition' )

    !.......  Allocate memory for arrays
    ALLOCATE( FIPSLIST( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FIPSLIST', PROGNAME )

    !.......  Read split definitions and create list of output file
    !           numbers
    IREC = 0
    CNT = 0
    DO I = 1, NLINES

        READ( FDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010)&
                'I/O error', IOS, 'reading COABDIS definitions '//&
                'file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        CALL PARSLINE( LINE, 20, SEGMENT )

        IF( .NOT. CHKINT( SEGMENT( 1 ) ) ) CYCLE

        CNT = CNT + 1

        FIPS = 6000 + STR2INT( SEGMENT(1) )
        ARBN = ADJUSTR( TRIM( SEGMENT(2) ) )
        DSTR = ADJUSTR( TRIM( SEGMENT(3) ) )
        CALL PADZERO( ARBN )
        CALL PADZERO( DSTR )

        WRITE( FIPSLIST( CNT ),'(A3,I6.6,A3)' ) ARBN, FIPS, DSTR

    END DO

    NFIPS = CNT

    CLOSE( FDEV )

    !.......  Exit if problem with COABDIS input file
    IF( EFLAG ) THEN
        MESG = 'Problem reading COABDIS definitions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  define modeling year
    MESG = 'Episode start time (HHMMSS)'
    YEAR = ENVINT( 'CONTROL_YEAR', MESG, 0, IOS )
    IF( YEAR < 1 ) THEN
        MESG = 'MUST define control modeling year.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    MESG = 'Define state code [default: 06]'
    ST = ENVINT( 'STATE_CODE', MESG, 6, IOS )

    !.......   Open output file
    MESG = 'Enter logical name of output file'
    ODEV = PROMPTFFILE( MESG, .FALSE., .TRUE., 'OUTFILE', PROGNAME )
    WRITE( ODEV,93000 ) '/CONTROL/'

    !.......   Get size for co-ab-dis description file
    NLINES = GETFLINE( CDEV, 'Control factor input file' )

    !.......  Read EIC definitions and create list of output file
    !           numbers
    IREC = 0
    CNT = 0
    DO I = 1, NLINES

        READ( CDEV, 93000, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010) 'I/O error', IOS,       &
                'reading EIC definitions file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        CALL PARSLINE( LINE, 20, SEGMENT )

        IF( .NOT. CHKINT( SEGMENT( 4 ) ) ) CYCLE

        YR = STR2INT( SEGMENT( 12 ) )       ! 2=full co/ab/dis, 4=co-specific, 5=ab-specific, 6=dis-specific
        IF( YEAR /= YR ) CYCLE             ! skip if modeling year is not met

        CNT = CNT + 1

        EIC = TRIM( SEGMENT( 9 ) )
        POL = TRIM( SEGMENT( 10 ) )
        FACS = STR2REAL( SEGMENT( 14 ) ) * 100.0
        IF( FACS < 100.0 ) FACS = ( 100.0 - FACS )

        !.......  Define regional flag
        REGID = STR2INT( SEGMENT( 4 ) )       ! 2=full co/ab/dis, 4=co-specific, 5=ab-specific, 6=dis-specific

        GEOCODE = ''
        IF( REGID == 2 ) THEN                        ! full co/ab/dis specific entry from CONTROL_FACS input file
            FIPS = ST * 1000 + STR2INT( SEGMENT(7) )
            ARBN = ADJUSTR( TRIM( SEGMENT(6) ) )
            DSTR = ADJUSTR( TRIM( SEGMENT(5) ) )
            CALL PADZERO( ARBN )
            CALL PADZERO( DSTR )
            WRITE( GEOCODE,'(A3,I6.6,A3)' ) ARBN, FIPS, DSTR

            WRITE( ODEV,94050 ) GEOCODE, TRIM(EIC), TRIM(POL), FACS, ',100,100,,,Y'

        ELSE IF( REGID == 4 ) THEN                   ! full county-specific entry from CONTROL_FACS input file

            DO J = 1, NFIPS
                FIPS = ST * 1000 + STR2INT( SEGMENT(7) )
                CNTY = STR2INT( FIPSLIST(J)(4:9) )
                IF( FIPS == CNTY ) THEN
                    GEOCODE = FIPSLIST( J )
                ELSE
                    CYCLE
                END IF

                WRITE( ODEV,94050 ) GEOCODE, TRIM(EIC), TRIM(POL), FACS, ',100,100,,,Y'

            END DO

        ELSE IF( REGID == 5 ) THEN                   ! full county-specific entry from CONTROL_FACS input file

            DO J = 1, NFIPS
                SUBC = FIPSLIST(J)( 1:3 )
                ARBN = ADJUSTR( TRIM( SEGMENT(6) ) )
                CALL PADZERO( ARBN )
                IF( ARBN == SUBC ) THEN
                    GEOCODE = FIPSLIST( J )
                ELSE
                    CYCLE
                END IF

                WRITE( ODEV,94050 ) GEOCODE, TRIM(EIC), TRIM(POL), FACS, ',100,100,,,Y'

            END DO

        ELSE IF( REGID == 6 ) THEN                   ! full county-specific entry from CONTROL_FACS input file

            DO J = 1, NFIPS
                SUBC = FIPSLIST(J)( 10:12 )
                DSTR = ADJUSTR( TRIM( SEGMENT(5) ) )
                CALL PADZERO( DSTR )
                IF( DSTR == SUBC ) THEN
                    GEOCODE = FIPSLIST( J )
                ELSE
                    CYCLE
                END IF

                WRITE( ODEV,94050 ) GEOCODE, TRIM(EIC), TRIM(POL), FACS, ',100,100,,,Y'

            END DO

        ELSE

            WRITE( MESG, 94010) 'ERROR: Missing regional flag value at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END DO

    WRITE( ODEV,93000 ) '/END/'

    CLOSE( ODEV )

    !.......  Exit if problem with splits input file
    IF( EFLAG ) THEN
        MESG = 'Problem converting CARB Adjustment factors to GCNTL input file.'
        IOS  = 2
    ELSE
        MESG = 'Successful completion'
        IOS  =0
    END IF
    CALL M3EXIT( PROGNAME, 0, 0, MESG, IOS )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I5, :, 2X ) )
94050 FORMAT( 3(A,','),',', F10.5, A )

END PROGRAM GCNTL4CARB

