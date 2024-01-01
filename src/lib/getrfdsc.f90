
REAL FUNCTION GETRFDSC( FILEINFO, KEY, REQUIRED )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !     Retreives a real value from the FDESC array or returns BADVAL3 if key is
    !     not required to be found and key is not present
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
    !
    !  REVISION  HISTORY:
    !       Created 7/99 by M Houyoux
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

    !.......   ARGUMENTS and their descriptions:

    CHARACTER(*), INTENT (IN) :: FILEINFO( * )     ! Array of file information
    CHARACTER(*), INTENT (IN) :: KEY               ! Key to find in FILEINFO
    LOGICAL     , INTENT (IN) :: REQUIRED          ! true: key must be found

    !.......   LOCAL VARIABLES their descriptions:

    INTEGER       L1, L2            ! length of strings
    INTEGER       I                 ! loop index
    INTEGER       K                 ! description string position of key

    REAL          RVAL              ! temporary real value

    CHARACTER(300) BUFFER           ! Key buffer
    CHARACTER(300) MESG             ! Message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETRFDSC'        ! Program name

    !***********************************************************************
    !   begin body of function GETRFDSC

    BUFFER = ADJUSTL( KEY )

    L1 = LEN_TRIM( BUFFER )

    DO I = 1, MXDESC3

        K = INDEX( FILEINFO( I ), BUFFER( 1:L1 ) )

        IF( K .GT. 0 ) THEN

            L2 = LEN_TRIM( FILEINFO( I ) )
            RVAL = STR2REAL( FILEINFO( I )( K+L1:L2 ) )

            IF( RVAL .LE. AMISS3 ) THEN
                MESG = 'ERROR: non-real result found at FDESC '//   &
                       'entry "'// KEY( 1:L1 )// '" in I/O API file'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            ELSE
                GETRFDSC = RVAL
                RETURN

            END IF
        END IF

    END DO

    !.....  If we get here, then key was not found in FDESC, so if it was
    !           required, then abort.

    IF( REQUIRED ) THEN
        MESG = 'FDESC3D packet "' // TRIM( KEY ) //     &
               '" was not found in I/O API file!'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE
        GETRFDSC = BADVAL3
        RETURN

    END IF

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.....94xxx

94010 FORMAT( 10( A, :, I6, :, 2X ) )

END FUNCTION GETRFDSC
