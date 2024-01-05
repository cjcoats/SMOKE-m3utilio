
SUBROUTINE GETDYSAV( NSRC, CIFIP, LDAYSAV )

    !***********************************************************************
    !  subroutine GETDYSAV body starts at line < >
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !      Creates an array of sources that are affected by daylight savings time.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO and related changes
    !
    !***************************************************************************
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

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER           , INTENT (IN) :: NSRC                 ! no. sources
    CHARACTER(FIPLEN3), INTENT (IN) :: CIFIP( NSRC )        ! country/state/co codes
    LOGICAL           , INTENT(OUT) :: LDAYSAV( NSRC )      ! true: source gets DLS

    !.......   EXTERNAL FUNCTIONS
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local allocatable arrays
    CHARACTER(FIPLEN3), ALLOCATABLE :: EXMPTCSC( : )

    !.......   Logical names and unit numbers
    INTEGER       FDEV       ! unit no. for daylight savings exemption file

    !.......   Other local variables
    INTEGER       K, N, S           ! counters and indices

    INTEGER       IOS               ! i/o status
    INTEGER       IREC              ! record counter
    INTEGER       NEXEMPT           ! no. entries in the exemptions file
    INTEGER       NLINES            ! no. lines in input file

    LOGICAL       EFLAG              ! true: error found
    LOGICAL       FFLAG              ! true: use daylight time exemption file

    CHARACTER(FIPLEN3) CSC           ! tmp country/state/county code
    CHARACTER(FIPLEN3) STA           ! tmp country/state code
    CHARACTER(FIPLEN3) CNY           ! tmp country code
    CHARACTER(300) TRAILER           ! ending part of output log message
    CHARACTER(300) MESG              ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETDYSAV'     ! program name

    !***********************************************************************
    !   begin body of subroutine GETDYSAV

    !.....  Initialize daylight savings indicator array
    EFLAG   = .TRUE.
    LDAYSAV = .TRUE.      ! array

    !.....  Get environment variable to decide whether file prompt is
    !           needed
    MESG  = 'Daylight time exemption switch'
    FFLAG = ENVYN( 'DAYLIGHT_EXEMPT', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "DAYLIGHT_EXEMPT"', 2 )
    ELSE IF( .NOT.FFLAG ) THEN
        RETURN      !.....  Return if no file was needed
    END IF

    FDEV = PROMPTFFILE( 'Enter logical name for DAYLIGHT TIME EXEMPT file', &
                        .TRUE., .TRUE., 'DLTEXEMPT', PROGNAME )

    !.....  Get number of lines in the exemptions file
    NLINES = GETFLINE( FDEV, 'Daylight exemptions file' )

    !.....  Allocate memory for the temporary arrays for reading in file
    ALLOCATE( EXMPTCSC( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EXMPTCSC', PROGNAME )

    !.....  Loop through file and read it
    DO N = 1, NLINES

        READ( FDEV, *, END=999, IOSTAT = IOS ) CSC

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )                                 &
                'I/O error', IOS, 'reading daylight ' //        &
                'savings exemptions file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        EXMPTCSC( N ) = CSC

    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem reading daylight savings exemptions file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    NEXEMPT = N

    TRAILER = ' exempted from daylight savings time'

    !.....  Loop through sources and see if any of the countries, states, or
    !           counties do not use daylight time
    DO S = 1, NSRC

        CSC =  CIFIP( S )

        !.....  Search for county
        K = FINDC( CSC,  NEXEMPT, EXMPTCSC )

        IF( K .GT. 0 ) THEN

            LDAYSAV( S ) = .FALSE.
            MESG = 'County ' // CSC // TRAILER
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.....  Search for state
        STA = CSC
        STA( STALEN3+1:FIPLEN3 ) = REPEAT( '0', FIPLEN3-STALEN3+1 )
        K = FINDC( STA,  NEXEMPT, EXMPTCSC )

        IF( K .GT. 0 ) THEN

            LDAYSAV( S ) = .FALSE.
            MESG = 'State ' // STA // TRAILER
            CALL M3MESG( MESG )
            CYCLE

        END IF

        !.....  Search for country
        CNY = CSC
        CNY( FIPEXPLEN3+2:FIPLEN3 ) = REPEAT( '0', FIPLEN3-FIPEXPLEN3+2 )
        K = FINDC( CNY, NEXEMPT, EXMPTCSC )

        IF( K .GT. 0 ) THEN

            LDAYSAV( S ) = .FALSE.
            MESG = 'Country ' // CNY // TRAILER
            CALL M3MESG( MESG )
            CYCLE

        END IF

    END DO

    !.....  Deallocate local memory
    DEALLOCATE( EXMPTCSC )

    RETURN

    !.....  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //              &
           'Check format of daylight' // CRLF() // BLANK5 //    &
           'savings exemptions file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )

94020 FORMAT( 10 ( A, 1X, I5, 2( 1X, A, 1X, F8.2 ) ) )

END SUBROUTINE GETDYSAV
