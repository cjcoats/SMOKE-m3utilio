
SUBROUTINE RDSTCY( FDEV, NDIM, INCNTYS )

    !***********************************************************************
    !  subroutine RDSTCY body starts at line 120
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to read the state and county names
    !      from the file for the specified list of county codes or for all
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" free source format,
    !        and related changes.  Bug-fix in INDEX1 call, line 456
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
    !.........  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: LSTCYPOP, STCYPOPYR, USEDAYLT,   &
                       NCOUNTRY, NSTATE,   NCOUNTY,     &
                       CTRYCOD,  STATCOD,  CNTYCOD,     &    
                       CTRYNAM,  STATNAM,  CNTYNAM,     &
                       CTRYPOPL, STATPOPL, CNTYPOPL,    &
                       CNTYTZON, CNTYTZNM

    IMPLICIT NONE

    !...........   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !...........   Subroutine arguments
    INTEGER,            INTENT (IN):: FDEV                ! county file unit no.
    INTEGER,            INTENT (IN):: NDIM                ! dim. for arrays or zero
    CHARACTER(FIPLEN3), INTENT (IN):: INCNTYS( NDIM )     ! input county codes or empty

    !...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT
    INTEGER, EXTERNAL :: GETFLINE

    !...........   Local parameters
    CHARACTER(16), PARAMETER :: CNTYPKT  = '/COUNTY/'
    CHARACTER(16), PARAMETER :: STATPKT  = '/STATE/'
    CHARACTER(16), PARAMETER :: CTRYPKT  = '/COUNTRY/'
    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSTCY'     ! program name

    !...........   Array for filtering state information
    INTEGER, ALLOCATABLE :: INSTATE( : )

    !...........   Arrays for determining sections of file
    INTEGER         ISKIP  ( 4 )      ! position of packets + total file length
    CHARACTER(16)   SECTION( 3 )      ! name of each section

    !...........   Other local variables

    INTEGER         F, I, J, K, L2, N            ! counters

    INTEGER         CNY             ! tmp county code
    INTEGER         COU             ! tmp country code
    INTEGER         COUST           ! tmp country/state code
    INTEGER         FIP             ! tmp country/state/county
    INTEGER      :: ISKIPCTR = 0    ! no. lines in FDEV to skip until countries
    INTEGER      :: ISKIPCNY = 0    ! no. lines in FDEV to skip until countys
    INTEGER      :: ISKIPSTA = 0    ! no. lines in FDEV to skip until states
    INTEGER         IOS             ! i/o status
    INTEGER         IREC            ! line number counter
    INTEGER         LCOU            ! country from previous iteration
    INTEGER         LCOUST          ! country/state from previous iteration
    INTEGER         LFIP            ! country/state/county from previous iter
    INTEGER         LSTA            ! state from previous iteration
    INTEGER         NDIMCY          ! number for allocating county table
    INTEGER         NDIMST          ! number for allocating state table
    INTEGER         NCNTY           ! no. entries in county section
    INTEGER         NCTRY           ! no. entries in country section
    INTEGER         NSTAT           ! no. entries in state section
    INTEGER         STA             ! tmp state code
    INTEGER, SAVE :: TZONE0          ! default time zone

    REAL            POP                    ! tmp population data

    LOGICAL      :: EFLAG                  ! true: error found
    LOGICAL      :: FILTER   = .FALSE.     ! true: filter county codes by input
    LOGICAL, SAVE :: FIRSTIME = .TRUE.      ! true: first time routine called
    LOGICAL      :: FOUNDCTR = .FALSE.     ! true: country packet found
    LOGICAL      :: FOUNDSTA = .FALSE.     ! true: state packet found
    LOGICAL      :: FOUNDCNY = .FALSE.     ! true: county packet found

    CHARACTER(FIPLEN3) CFIP      ! tmp character FIPS code
    CHARACTER(FIPLEN3) CSTA      ! tmp character state code
    CHARACTER(FIPLEN3) CCOU      ! tmp character country code
    CHARACTER       DLCHR        ! tmp daylight time exemptions flag
    CHARACTER(3)    TZN          ! tmp time zone
    CHARACTER(300)  LINE         ! line read buffer
    CHARACTER(300)  MESG         ! message buffer

    CHARACTER(FIPLEN3)  FBUF         ! message buffer

    !***********************************************************************
    !   begin body of subroutine RDSTCY

    IF( FIRSTIME ) THEN
        MESG = 'Default time zone for sources'
        TZONE0 = ENVINT( 'SMK_DEFAULT_TZONE', MESG, 5, IOS )
        FIRSTIME = .FALSE.
    END IF

    EFLAG = .FALSE.
    MESG  = 'Reading state and county names and time zones...'
    CALL M3MSG2( MESG )

    !.........  Loop through lines of file and determine dimension for county
    !           arrays, state arrays, and country arrays

    IREC = 0
    K = 0
    DO

        READ ( FDEV, 93000, END=101, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG, 94010 )                                    &
                  'I/O error', IOS, 'reading country, state,  '//   &
                  'county names file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        !.............  Look for header lines
        IF( LINE( 1:1 ) .EQ. CINVHDR ) THEN

            !.................  Look for population header, and if found...
            L2 = LEN_TRIM( LINE )
            IF( LINE( 2:11 ) .EQ. 'POPULATION' ) THEN

                !.....................  If space for year is blank, then error
                IF( LINE( 12:L2 ) .EQ. ' ' ) THEN

                    MESG = 'Year must be included at end of #POPULATION  header.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                !.....................  If space for year is integer, then convert and store it
                ELSE IF( CHKINT( LINE( 12:L2 ) ) ) THEN
                    LSTCYPOP = .TRUE.
                    STCYPOPYR = STR2INT( LINE( 12:L2 ) )

                    WRITE( MESG,94010 ) 'NOTE: Population data '//      &
                           'read from country/state/county file '//     &
                           'for year ', STCYPOPYR
                    CALL M3MSG2( MESG )

                !.....................  If space for year is non-integer, then error
                ELSE

                    MESG = 'Non-integer year value at end of #POPULATION  header.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

                !.................  Population header not found
            ELSE
                MESG = 'NOTE: No population data will be read ' //      &
                       'from country/state/county file.'
                CALL M3MSG2( MESG )

            END IF

            CYCLE

        END IF

        !.............  Search for packets in every line for error checking
        FOUNDCTR = ( INDEX( LINE, CTRYPKT ) .GT. 0 )
        FOUNDSTA = ( INDEX( LINE, STATPKT ) .GT. 0 )
        FOUNDCNY = ( INDEX( LINE, CNTYPKT ) .GT. 0 )

        !.............  Error if country packet is found twice
        IF ( FOUNDCTR .AND. ISKIPCTR .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Country packet found ' //   &
                   'again at line', IREC, 'but also at line',       &
                   ISKIPCTR
            CALL M3MSG2( MESG )

        !.............  Store position and order if found for the first time
        ELSE IF( FOUNDCTR ) THEN
            K = K + 1
            ISKIPCTR     = IREC
            ISKIP  ( K ) = IREC
            SECTION( K ) = CTRYPKT
        END IF

        !.............  Error if state packet is found twice
        IF ( FOUNDSTA .AND. ISKIPSTA .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: State packet found ' //     &
                   'again at line', IREC, 'but also at line', ISKIPSTA
            CALL M3MSG2( MESG )

        !.............  Store position and order if found for the first time
        ELSE IF( FOUNDSTA ) THEN
            K = K + 1
            ISKIPSTA     = IREC
            ISKIP  ( K ) = IREC
            SECTION( K ) = STATPKT
        END IF

        !.............  Error if county packet is found twice
        IF ( FOUNDCNY .AND. ISKIPCNY .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: County packet found ' //    &
                   'again at line', IREC, 'but also at line', ISKIPCNY
            CALL M3MSG2( MESG )

        !.............  Store position and order if found for the first time
        ELSE IF( FOUNDCNY ) THEN
            K = K + 1
            ISKIPCNY     = IREC
            ISKIP  ( K ) = IREC
            SECTION( K ) = CNTYPKT
        END IF

    END DO

101 CONTINUE           ! exit from read loop

    !.........  Make sure the file had all sections
    IF( K .NE. 3 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Country, state, and county sections must ' //&
               'be present in file'
        CALL M3MESG( MESG )
    END IF

    !.........  Abort if error encountered
    IF( EFLAG ) THEN
        MESG = 'Problem found in country, state, county file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.........  Determine the number of entries in each section.
    MESG = 'Country, state  county names file'
    ISKIP( 4 ) = GETFLINE( FDEV, MESG ) + 1
    DO K = 1, 3

        IF( SECTION( K ) .EQ. CTRYPKT ) THEN
            NCTRY = ISKIP( K+1 ) - ISKIP( K ) - 1

        ELSE IF( SECTION( K ) .EQ. STATPKT ) THEN
            NSTAT = ISKIP( K+1 ) - ISKIP( K ) - 1

        ELSE IF( SECTION( K ) .EQ. CNTYPKT ) THEN
            NCNTY = ISKIP( K+1 ) - ISKIP( K ) - 1

        END IF

    END DO

    !.........  Set up size for county read depending on the input dimension
    IF( NDIM .GT. 1 ) THEN
        NDIMCY = NDIM
        FILTER = .TRUE.
    ELSE
        NDIMCY = NCNTY
    END IF

    !.........  Create state table for filtering
    NDIMST = NSTAT
    IF( FILTER ) THEN

        !.............  First count the number of states
        LSTA = -9
        N    = 0
        DO I = 1, NDIM

            STA = STR2INT( INCNTYS( I ) ) / 1000
            IF( STA .NE. LSTA ) THEN
                N = N + 1
                LSTA = STA
            END IF

        END DO
        NDIMST = N

        ALLOCATE( INSTATE( NDIMST ), STAT=IOS )
        CALL CHECKMEM( IOS, 'INSTATE', PROGNAME )

        LSTA = -9
        N    = 0
        DO I = 1, NDIM

            STA = STR2INT( INCNTYS( I ) ) / 1000
            IF( STA .NE. LSTA ) THEN
                N = N + 1
                INSTATE( N ) = STA
                LSTA = STA
            END IF

        END DO

    END IF

    !.........  Allocate memory for data arrays from MODSTCY module
    ALLOCATE( CTRYCOD( NCTRY ),     &
              CTRYNAM( NCTRY ),     &
             CTRYPOPL( NCTRY ),     &
              STATCOD( NDIMST ),    &
              STATNAM( NDIMST ),    &
             STATPOPL( NDIMST ),    &
              CNTYCOD( NDIMCY ),    &
              CNTYNAM( NDIMCY ),    &
             CNTYTZON( NDIMCY ),    &
             CNTYTZNM( NDIMCY ),    &
             CNTYPOPL( NDIMCY ),    &
             USEDAYLT( NDIMCY ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTRYCOD...USEDAYLT', PROGNAME )

    !.........  Initialize
    CTRYCOD  = ' '
    STATCOD  = ' '
    CNTYCOD  = ' '
    CNTYTZON = -9        ! array
    CTRYNAM  = ' '
    STATNAM  = ' '
    CNTYNAM  = ' '
    CTRYPOPL = BADVAL3
    STATPOPL = BADVAL3
    CNTYPOPL = BADVAL3
    USEDAYLT = .FALSE.

    !.........  Read in country section of the file...
    !.........  Skip ahead in file to country section
    REWIND ( FDEV )
    CALL SKIPL( FDEV, ISKIPCTR )

    !.........  Loop through and read country data
    IREC = ISKIPCTR
    N    = 0
    LCOU = -9
    DO N = 1, NCTRY

        READ ( FDEV, 93000, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( .NOT. CHECK_READ_STATUS() ) CYCLE

        COU = STR2INT( LINE( 1:1 ) )
        IF( COU .LT. LCOU ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Unsorted country code '//   &
                   'at line', IREC
            CALL M3MESG( MESG )

        ELSE IF( COU .EQ. LCOU ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Duplicate country code '//  &
                   'at line', IREC
            CALL M3MESG( MESG )
        END IF

        CTRYCOD( N ) = REPEAT( '0', FIPLEN3 )
        CTRYCOD( N )( FIPEXPLEN3+1:FIPEXPLEN3+1 ) = LINE( 1:1 )
        CTRYNAM( N ) = ADJUSTL( LINE( 3:22 ) )

        LCOU = COU

    END DO

    !.........  Read in state section of the file...
    !.........  Skip ahead in file to state section
    REWIND( FDEV )
    CALL SKIPL( FDEV, ISKIPSTA )

    !.........  Loop through and read state data
    IREC = ISKIPSTA
    K = 0
    LCOUST = -9
    DO N = 1, NSTAT

        READ ( FDEV, 93000, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( .NOT. CHECK_READ_STATUS() ) CYCLE

        COU = MAX( STR2INT( LINE( 1:1 ) ), 0 )
        STA = MAX( STR2INT( LINE( 2:3 ) ), 0 )
        COUST = COU * 100 + STA

        !.............  Find state code in valid list
        IF( FILTER ) THEN
            J = FIND1( COUST, NDIMST, INSTATE )
            IF( J .LE. 0 ) CYCLE
        END IF

        IF( COUST .LT. LCOUST ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Unsorted state code '//     &
                   'at line', IREC
            CALL M3MESG( MESG )

        ELSE IF( COUST .EQ. LCOUST ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Duplicate state code '//    &
                   'at line', IREC
            CALL M3MESG( MESG )
        END IF

        K = K + 1
        STATCOD( K ) = REPEAT( '0', FIPLEN3 )
        WRITE( STATCOD( K )( FIPEXPLEN3+1:FIPEXPLEN3+1 ), '(I1)' ) COU
        WRITE( STATCOD( K )( FIPEXPLEN3+2:FIPEXPLEN3+3 ), '(I2.2)' ) STA
        STATNAM( K ) = ADJUSTL( LINE( 7:26 ) )

        LCOUST = COUST

    END DO

    !.........  Check if input states all have information in state codes file
    IF( K .NE. NDIMST ) THEN

        IF( FILTER ) THEN

            !.................  Loop through input states and report missing
            DO N = 1, NDIMST
                WRITE( FBUF, '(I12)' ) INSTATE( N ) * 1000
                J = INDEX1( FBUF, K, STATCOD )
                IF( J .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,'(A,1X,I3.3,A)' )                           &
                      'ERROR: Input data contains country/state code ',     &
                      INSTATE( N ), ', but it '// CRLF()// BLANK10//        &
                      'is not found in state/county codes file.'
                    CALL M3MSG2( MESG )
                END IF
            END DO

        END IF

        MESG = 'INTERNAL ERROR: Actual count of state codes in error'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.........  Read in county section of the file...
    !.........  Skip ahead in file to county section
    REWIND( FDEV )
    CALL SKIPL( FDEV, ISKIPCNY )

    !.........  Loop through and read county data
    IREC = ISKIPCNY
    K    = 0
    LFIP = -9
    DO N = 1, NCNTY

        READ ( FDEV, 93000, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( .NOT. CHECK_READ_STATUS() ) CYCLE

        COU = MAX( STR2INT( LINE( 26:26 ) ), 0 )
        STA = MAX( STR2INT( LINE( 27:28 ) ), 0 )
        CNY = MAX( STR2INT( LINE( 29:31 ) ), 0 )
        TZN = LINE( 40:42 )
        DLCHR = LINE( 43:43 )
        IF( LSTCYPOP ) POP = STR2REAL( LINE( 114:128 ) )

        POP = MAX( 0., POP )                   ! Remove missing values

        FIP = COU * 100000 + STA * 1000 + CNY
        WRITE( CFIP, '(I12.12)' ) FIP

        !.............  If input codes have been provided, find current code in the
        !               list, or skip to next iteration
        IF( FILTER ) THEN

            J = FINDC( CFIP, NDIM, INCNTYS )

            IF( J .LE. 0 ) CYCLE            ! to next iteration

        END IF

        !.............  Find the time zone in the list and retrieve the integer value
        J = INDEX1( TZN, MXTZONE, TZONNAM )

        IF( FIP .LT. LFIP ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Unsorted county code '//    &
                   'at line', IREC
            CALL M3MESG( MESG )

        ELSE IF( FIP .EQ. LFIP ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Duplicate county code '//   &
                   'at line', IREC
            CALL M3MESG( MESG )
        END IF

        !.............  Store the county-specific information
        K = K + 1

        IF( K .LE. NDIMCY ) THEN

            !.................  Store region code and name
            CNTYCOD( K ) = CFIP
            CNTYNAM( K ) = ADJUSTL( LINE( 5:24 ) )

            !.................  Store population data, if included
            IF( LSTCYPOP ) CNTYPOPL( K ) = POP

            !.................  Store the status of the daylight time exemptions
            USEDAYLT( K ) = ( DLCHR .EQ. ' ' )

            !.................  If the time zone is not defined, apply default.
            IF( J .GT. 0 ) THEN
                CNTYTZON( K ) = TZONNUM( J )

            ELSE
                IF( FIP .NE. LFIP ) THEN
                    WRITE( MESG,94010 )                                 &
                      'WARNING: Applying default time zone', TZONE0,    &
                      'to country/state/county code:', FIP
                    CALL M3MESG( MESG )
                END IF

                CNTYTZON( K ) = TZONE0

            END IF

            !.................  Store time zone name (do not apply default time zone)
            CNTYTZNM( K ) = TZN

        END IF

        LFIP = FIP

    END DO               !  end of loop through counties

    !.........  Error if the counties read are less than the expected number
    IF( NDIM .LE. 0 .AND. K .NE. NDIMCY ) THEN
        MESG = 'INTERNAL ERROR: Actual count of county codes in error'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.........  Error if the counties read are less than the number requested
    !           by the calling program.
    ELSE IF( K .LT. NDIMCY ) THEN
        MESG = 'ERROR:  Some requested counties not found in county names file:'
        CALL M3MSG2( MESG )

        DO J = 1, NDIM
            CFIP = INCNTYS( J )
            I = INDEX1( CFIP, K, CNTYCOD )

            IF( I .LE. 0 ) THEN
                WRITE( MESG,93000 ) BLANK10 // 'Code:' // CFIP
                CALL M3MSG2( MESG )
            END IF

        END DO

        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

    END IF

    !.........  Set final sizes that are stored in module
    NCOUNTRY = NCTRY
    NSTATE   = NDIMST
    NCOUNTY  = NDIMCY

    !.........  Store state and country populations
    IF ( LSTCYPOP ) THEN

        CTRYPOPL = 0.      ! array
        STATPOPL = 0.      ! array

        DO J = 1, NCOUNTY

            CFIP = CNTYCOD( J )
            CCOU = CFIP
            CCOU( FIPEXPLEN3+2:FIPLEN3 ) = REPEAT( '0', FIPLEN3-FIPEXPLEN3+2 )
            CSTA = CFIP
            CSTA( STALEN3+1:FIPLEN3 ) = REPEAT( '0', FIPLEN3-STALEN3+1 )

            !.................  Add to country population
            F = INDEX1( CCOU, NCOUNTRY, CTRYCOD )
            IF ( F .GT. 0 ) THEN
                CTRYPOPL( F ) = CTRYPOPL( F ) + CNTYPOPL( J )
            ELSE
                MESG = 'INTERNAL ERROR: Could not find country '//  &
                       CCOU // ' in list while computing country population.'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

            !.................  Add to state population
            F = INDEX1( CSTA, NSTATE, STATCOD )
            IF ( F .GT. 0 ) THEN
                STATPOPL( F ) = STATPOPL( F ) + CNTYPOPL( J )
            ELSE
                MESG = 'INTERNAL ERROR: Could not find state '//    &
                       CSTA // ' in list while computing state population.'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
            END IF

        END DO

    END IF        ! if population data are present

    !.........  Deallocate locally allocated memory
    IF( ALLOCATED( INSTATE ) ) DEALLOCATE( INSTATE )

    IF ( EFLAG ) THEN
        MESG = 'Problem reading country/state/county file. '//      &
               'See previously listed errors.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )


CONTAINS    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


    !..............  This internal subprogram checks the read status and reports
    !                and error if there was one.  It sets EFLAG and skips to
    !                next record

    LOGICAL FUNCTION CHECK_READ_STATUS()

        !.............................................................................

        CHECK_READ_STATUS = .TRUE.
        IF( IOS .NE. 0 ) THEN

            WRITE( MESG, 94010 )                                    &
                  'I/O error', IOS, 'reading country, state,  '//   &
                  'county names file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            CHECK_READ_STATUS = .FALSE.

        END IF

        !.............................................................................

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END FUNCTION CHECK_READ_STATUS

END SUBROUTINE RDSTCY
