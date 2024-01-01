
SUBROUTINE INITEM( LDEV, NAMEVERS, INPROGNM )

    !**********************************************************************
    !  program body starts at line 81
    !
    !  DESCRIPTION:
    !       Writes out abridged copyright information, calling program version,
    !       web address for documentation, general program and info, program-
    !       specific info, and prompts to continue running the calling program
    !
    !  PRECONDITIONS REQUIRED:
    !       Unit number LDEV defined
    !       Version name set in calling program using SCCS
    !       Calling program name defined
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !     Subroutines: Models-3 subroutines, PROGDESC
    !     Functions: Models-3 functions
    !
    !  REVISION  HISTORY:
    !       Copied from prototype version 1.2  in 10/98 by M Houyoux
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
    !****************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.....  Subroutine arguments
    INTEGER,       INTENT (INOUT) :: LDEV       ! Log file unit number
    CHARACTER(50), INTENT (IN) :: NAMEVERS      ! ASCII field e/ version number appended to end
    CHARACTER(16), INTENT (IN) :: INPROGNM      ! Calling program name

    !.....  Parameters
    INTEGER, PARAMETER :: STDOUT = 6
    INTEGER, PARAMETER :: YEAR   = 2004

    !.....  Local variables
    REAL          VERSION

    INTEGER       IOUT( 2 )      ! output unit numbers for stdout and logfile
    INTEGER       I, J, L
    INTEGER       NLOOP          ! Number of times to loop through output

    CHARACTER(50)  VERCHAR
    CHARACTER(300) LINE0, LINE1, LINE2, LINE3, LINE4, LINE5

    CHARACTER(16), PARAMETER :: PROGNAME = 'INITEM'      ! program name

    !***********************************************************************
    !   begin body of program INITEM

    LINE0 = 'SMOKE ---------------'
    WRITE( LINE1,94020 ) 'Copyright (c)', YEAR, 'Center for Environmental Modeling for Policy Development'

    LINE2 = 'All rights reserved'

    LINE3 = 'Online documentation available at:'
    LINE4 = '    https://cmascenter.org/smoke/'

    !.....  Set up program version information
    VERCHAR = ADJUSTL( NAMEVERS )

    IF( VERCHAR( 1:1 ) .EQ. '%' ) THEN
        VERCHAR = 'Dev'
    ELSE
        J = INDEX( VERCHAR, ' ' )
        L = MAX( J+1,LEN_TRIM( VERCHAR ) )
        VERCHAR = ADJUSTL( VERCHAR( J+1:L ) )
        J = INDEX( VERCHAR, '$' )
        IF( J .GT. 1 ) VERCHAR = VERCHAR( 1:J-1 )
    ENDIF

    !.....  Set up writing loop
    IOUT( 1 ) = STDOUT
    IOUT( 2 ) = LDEV
    IF( LDEV .NE. STDOUT ) THEN
        NLOOP = 2
    ELSE
        NLOOP = 1
    ENDIF

    DO I = 1, NLOOP

        LDEV = IOUT( I )

        !.....  Write copyright information

        WRITE( LDEV,92000 ) TRIM( LINE0 )
        WRITE( LDEV,92000 ) TRIM( LINE1 )
        WRITE( LDEV,92000 ) TRIM( LINE2 )
        WRITE( LDEV,92000 )

        !.....  Write program version information

        L = MAX( LEN_TRIM( VERCHAR ), 1 )
        WRITE( LDEV,92010 ) TRIM( INPROGNM ), VERCHAR( 1:L )

        !.....  Write web site information
        WRITE( LDEV,92000 ) LINE3( 1:LEN_TRIM( LINE0 ) )
        WRITE( LDEV,92000 ) LINE4( 1:LEN_TRIM( LINE1 ) )

        !.....  Write program-specific information
        CALL PROGDESC( LDEV, INPROGNM )

        !.....  Write general information for SMOKE programs

        WRITE( LDEV,92000 )                                             &
        &' ',                                                           &
        &'You will need to enter the logical names for the input and',  &
        &'output files (and to have set them prior to program start,',  &
        &'using "setenv <logicalname> <pathname>").',                   &
        &' ',                                                           &
        &'You may use END_OF-FILE (control-D) to quit the program',     &
        &'during logical-name entry. Default responses are given in',   &
        &'brackets [LIKE THIS] and can be accepted by hitting the', &
        &'<RETURN> key.',&
        &' '

    ENDDO      ! End of write loop

    IF ( .NOT. GETYN( 'Continue with program?', .TRUE. ) ) THEN
        CALL M3EXIT( INPROGNM, 0, 0, 'Ending program.', 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Informational (LOG) message formats... 92xxx

92000 FORMAT( 5X, A )

92010 FORMAT( 5X, 'Program ', A, ', Version ', A )

    !.......   Internal buffering formats.... 94xxx

94020 FORMAT( 10( A, :, I4, :, 1X )  )

END SUBROUTINE INITEM
