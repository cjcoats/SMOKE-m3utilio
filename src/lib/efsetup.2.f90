
SUBROUTINE EFSETUP( FNAME, MODELNAM, VOLNAM )

    !***********************************************************************
    !  subroutine EFSETUP body starts at line < >
    !
    !  DESCRIPTION:
    !      Get the names of the emission factors given the emission factor
    !      model name and locally obtained environment variable settings
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION HISTORY:
    !       Created ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   MODULES for public variables
    !.....  This module contains emission factor tables and related
    USE MODEMFAC, ONLY: NEFS, EFSNAM, EFSUNT, EFSDSC

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'M6CNST3.h90'       !  Mobile6 constants

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FNAME               ! logical file or 'NONE'
    CHARACTER(*), INTENT (IN) :: MODELNAM            ! name of EF model
    CHARACTER(*), INTENT(OUT) :: VOLNAM              ! volatile pollutant name

    !.......   Local variables
    INTEGER         I, J, K, L, L2, N     ! counters and indices

    INTEGER         IOS         ! status from retrieving E.V.s
    INTEGER         LJ          ! length of emission type joiner
    INTEGER         ML          ! length of MODELNAM buffer
    INTEGER         PINDX       ! polluntant index

    LOGICAL, SAVE :: FIRSTIME = .TRUE.       ! true: first time routine called

    CHARACTER(256)  MESG             ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'EFSETUP'     ! program name

    !***********************************************************************
    !   begin body of subroutine EFSETUP

    !.....  Set flag for no input file available
    IFLAG = ( FNAME .EQ. 'NONE' )

    IF ( MODELNAM .NE. 'MOBILE6' ) THEN
        MESG = 'Model name "' // TRIM( MODELNAM ) // '" not recognized'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE IF ( FNAME .NE 'NONE' ) THEN
        RETURN
    END IF
    
    IF ( FIRSTIME ) THEN

        !.....  For new file, get environment variable for the volatile pol
        !      for mobile sources. For now, this routine only knows MOBILE6

        MESG = 'Volatile pollutant type'
        CALL ENVSTR( 'MB_HC_TYPE', MESG, 'TOG', VOLNAM, IOS)

        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MB_HC_TYPE"', 2 )
        ELSE  IF ( IOS .LT. 0 ) THEN
            MESG = 'WARNING: Default mobile volatile pollutant "'// TRIM( VOLNAM ) // '" used.'
        ELSE
            MESG = 'NOTE: Using mobile volatile pollutant "' // TRIM( VOLNAM ) // '".'
        END IF

        FIRSTIME = .FALSE.

    END IF      !!  if firstime

    !.....  Confirm that the pollutant of interest is valid for the model selected.

    PINDX = INDEX1( VOLNAM( 1:4 ), NM6VPOL, M6VPOLS )

    IF( PINDX .LE. 0 ) THEN

        MESG  = 'ERROR: Volatile pollutant type "' // TRIM( VOLNAM ) //  &
                '" is invalid for the ' // TRIM( MODELNAM ) // ' model'
        CALL M3MSG2( MESG )

        MESG = 'Problem configuring for ' // MODELNAM( 1:ML ) // ' emission factor model.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF      !!  if iflag

    !.....  Store variable names, units, and descriptions, depending on the
    !       model being used.
    !.....  Set variable names, units, and descriptions from
    !                   arrays defined in the MOBILE6 include file

    NEFS = MXM6EFS

    !.....  Allocate memory for and set the public variables for the
    !       EF names, units, and descriptions...
    ALLOCATE( EFSNAM( NEFS ),       &
              EFSUNT( NEFS ),       &
              EFSDSC( NEFS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EFSNAM...EFSDSC', PROGNAME )
    EFSNAM = ' '              ! array
    EFSUNT = ' '
    EFSDSC = ' '

    K = 0

    !.....  Loop over all emission processes, all pollutants

    DO I = 1, MXM6EPR
    DO J = 1, MXM6POLS

        !.....  Check if this is a valid pollutant/process combo
        IF( M6POL2EF( I,J ) == -1 ) CYCLE

        K = K + 1

        !.....  Pull process name and description from include file
        EFSNAM( K ) = M6PROCS( I ) // ETJOIN
        EFSDSC( K ) = 'EFs for ' // M6PRCDSC( I )

        !.....  If pollutant is HC, append specified volatile pollutant
        !       name; otherwise, use name and desc from include file
        IF( M6POLS( J ) == 'HC' ) THEN
            EFSNAM( K ) = TRIM( EFSNAM( K ) ) // VOLNAM
            EFSDSC( K ) = TRIM( EFSDSC( K ) ) // ' ' // VOLNAM
        ELSE
            EFSNAM( K ) = TRIM( EFSNAM( K ) ) // M6POLS( J )
            EFSDSC( K ) = TRIM( EFSDSC( K ) ) // ' ' // M6POLDSC( J )
        END IF

        !.....  Store units from include file
        EFSUNT( K ) = M6UNIT

    END DO          ! pollutant loop
    END DO          ! emission process loop

    RETURN

END SUBROUTINE EFSETUP
