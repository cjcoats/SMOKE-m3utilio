
SUBROUTINE RDSPDPROF( FDEV )

!***********************************************************************
!  subroutine body starts at line 117
!
!  DESCRIPTION:
!      This subroutine reads the speed profiles file, then creates
!      sorted arrays by profile number.
!
!  PRECONDITIONS REQUIRED:
!      FDEV is opened
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!**************************************************************************
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

!...........   Modules for public variables
!...........  This module is for mobile-specific data
    USE MODMOBIL, ONLY: SPDPROFS, SPDNUMS, NSPDPROF

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Subroutine arguments
    INTEGER      , INTENT  (IN) :: FDEV       ! File unit number

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

!...........   Local parameters
    INTEGER    , PARAMETER :: MXCOL = 25

!...........   Array of input fields
    CHARACTER(SPDLEN3)  SEGMENT( MXCOL )

!...........   Local, allocatable arrays
    REAL,    ALLOCATABLE :: SPDPROFA( :,: )   ! unsorted hourly speed profiles
    INTEGER, ALLOCATABLE :: SPDNUMA ( : )     ! unsorted profile numbers
    INTEGER, ALLOCATABLE :: SPDIDXA ( : )     ! sorting index

!...........   Local variables
    INTEGER         I,J,K                     ! counters and indices

    INTEGER         IOS                       ! I/O status
    INTEGER         NLINES                    ! no. lines in file
    INTEGER         PREVNUM                   ! previous profile number

    LOGICAL      :: EFLAG = .FALSE.           ! true: error found

    CHARACTER(150)  LINE                      ! Read buffer for a line
    CHARACTER(256)  MESG                      ! message buffer

    CHARACTER(16) :: PROGNAME = 'RDSPDPROF'    ! program name

!***********************************************************************
!   Begin body of subroutine RDSPDPROF

!........  Count number of lines in file
    NLINES = GETFLINE( FDEV, 'Speed profiles file')

!........  Allocate unsorted arrays
    ALLOCATE( SPDPROFA( NLINES, 24 ),&
    &           SPDNUMA( NLINES ),&
    &           SPDIDXA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPDPROFA...SPDIDXA', PROGNAME )

    K = 0

!........  Loop through file
    DO I = 1, NLINES

        READ( FDEV, 93000, IOSTAT=IOS ) LINE

!............  Check for I/O errors
        IF( IOS > 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )&
            &   'I/O error', IOS,&
            &   'reading speed profiles file at line', I
            CALL M3MESG( MESG )
            CYCLE
        END IF

!............  Check for end of file
        IF( IOS < 0 ) THEN
            MESG = 'End of file reached unexpectedly. ' //&
            &       'Check format of speed ' // CRLF() // BLANK5 //&
            &       'profiles file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!............  Skip blank lines and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!............  Parse line into segments
        CALL PARSLINE( LINE, MXCOL, SEGMENT )

!............  Check that the profile number is a positive integer
        IF( .NOT. CHKINT( SEGMENT( 1 ) ) .OR.&
        &    STR2INT( SEGMENT( 1 ) ) < 1 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Profile number is not ' //&
            &       'an integer or is less than 1 at line', I
            CALL M3MESG( MESG )
        END IF

!............  Check that hourly speeds are valid
        DO J = 2, 25
            IF( .NOT. CHKREAL( SEGMENT( J ) ) .OR.&
            &    STR2REAL( SEGMENT( J ) ) < 0. ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'ERROR: Invalid hourly ' //&
                &       'speed for hour', J-1, 'at line', I
                CALL M3MESG( MESG )
            END IF
        END DO

!............  Skip rest of loop if there has been an error
        IF( EFLAG ) CYCLE

        K = K + 1
        IF( K .GT. NLINES ) CYCLE  ! Ensure no overflow

!............  Store profile information
        SPDIDXA( K ) = K
        SPDNUMA( K ) = STR2INT( SEGMENT( 1 ) )

        DO J = 2, 25
            SPDPROFA( K,J-1 ) = STR2REAL( SEGMENT( J ) )
        END DO

    END DO  ! End loop reading SPDPROF file

    NSPDPROF = K

    IF( NSPDPROF == 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No valid entries in SPDPROF file!'
        CALL M3MSG2( MESG )

    ELSEIF( NSPDPROF > NLINES ) THEN
        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //&
        &       'storing SPDPROF file was', NLINES,&
        &       CRLF() // BLANK10 // 'but actually needed', NSPDPROF
        CALL M3MSG2( MESG )

    END IF

!........  Check for errors reading SPDPROF file and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading speed profiles file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!........  Sort speed profiles by profile number
    CALL SORTI1( NSPDPROF, SPDIDXA, SPDNUMA )

!........  Allocate memory for sorted arrays
    ALLOCATE( SPDPROFS( NSPDPROF, 24 ),&
    &           SPDNUMS( NSPDPROF ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPDNUMS', PROGNAME )

!........  Create sorted profile arrays
    PREVNUM = 0

    DO I = 1, NSPDPROF
        J = SPDIDXA( I )

!............  Check that duplicate profile numbers have not been used
        IF( SPDNUMA( J ) == PREVNUM ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: At least two speed ' //&
            &       'profiles are assigned the profile number',&
            &       PREVNUM
            CALL M3MESG( MESG )
            CYCLE
        END IF

        SPDNUMS( I ) = SPDNUMA( J )
        PREVNUM = SPDNUMA( J )

        SPDPROFS( I,: ) = SPDPROFA( J,: )
    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem with speed profiles file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!........  Deallocate local arrays
    DEALLOCATE( SPDIDXA, SPDNUMA, SPDPROFA )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSPDPROF

