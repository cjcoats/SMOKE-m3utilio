
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

    LOGICAL       :: EFLAG    = .FALSE.      ! true: processing error found
    LOGICAL, SAVE :: FIRSTIME = .TRUE.       ! true: first time routine called
    LOGICAL       :: IFLAG    = .FALSE.      ! true: input file available

    CHARACTER(300)                  MESG             ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'EFSETUP'     ! program name

    !***********************************************************************
    !   begin body of subroutine EFSETUP

    !.....  Set flag for no input file available
    IFLAG = ( FNAME .EQ. 'NONE' )

    !.....  Process for MOBILE6 model
    IF ( MODELNAM .EQ. 'MOBILE6' ) THEN

    !.....  For new file, get environment variable for the volatile pol
    !               for mobile sources. For now, this routine only knows MOBILE6
        IF( IFLAG ) THEN
            IF( FIRSTIME ) THEN
                MESG = 'Volatile pollutant type'
                CALL ENVSTR( 'MB_HC_TYPE', MESG, 'TOG', VOLNAM, IOS)
            END IF

    !.....  Create message about the volatile pollutant that is being used
            IF ( IOS .LT. 0 ) THEN
                MESG = 'WARNING: Default mobile volatile ' //       &
                       'pollutant "'// TRIM( VOLNAM ) // '" used.'

            ELSE
                MESG = 'NOTE: Using mobile volatile pollutant "' // TRIM( VOLNAM ) // '".'
            END IF

    !.....  For existing file, confirm MOBILE6 names
        ELSE
    ! note: add here

        END IF

    !.....  Abort if emission factor model name is not known
    ELSE
        MESG = 'Model name "' // TRIM( MODELNAM ) //        &
               '" not recognized by program ' // PROGNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.....  Confirm that the pollutant of interest is valid for the model
    !           selected.
    IF( IFLAG ) THEN
        SELECT CASE( MODELNAM )

          CASE( 'MOBILE6' )

            L = LEN_TRIM( VOLNAM )
            PINDX = INDEX1( VOLNAM( 1:4 ), NM6VPOL, M6VPOLS )

            IF( PINDX .LE. 0 ) THEN

                EFLAG = .TRUE.
                MESG  = 'ERROR: Volatile pollutant type "' // TRIM( VOLNAM ) //  &
                        '" is invalid for the ' // TRIM( MODELNAM ) // ' model'
            END IF

        END SELECT
    END IF

    !.....  Write out previously prepared message about the pollutant of
    !           interest.  This could have been set 3 sections above.
    IF( IFLAG ) CALL M3MSG2( MESG )

    !.....  Abort if there has been a fatal problem up to this point
    IF( EFLAG ) THEN
        MESG = 'Problem configuring for ' // MODELNAM( 1:ML ) // ' emission factor model.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Store variable names, units, and descriptions, depending on the
    !           model being used.
    SELECT CASE( MODELNAM )

      CASE( 'MOBILE6' )

        IF( IFLAG ) THEN

    !.....  Set variable names, units, and descriptions from
    !                   arrays defined in the MOBILE6 include file

            NEFS = MXM6EFS

    !.....  Allocate memory for and set the public variables for the
    !                   EF names, units, and descriptions...
            ALLOCATE( EFSNAM( NEFS ),       &
                      EFSUNT( NEFS ),       &
                      EFSDSC( NEFS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EFSNAM...EFSDSC', PROGNAME )
            EFSNAM = ' '      ! array
            EFSUNT = ' '
            EFSDSC = ' '

            K = 0

    !.....  Loop over all emission processes
            DO I = 1, MXM6EPR

    !.....  Loop over all pollutants
                DO J = 1, MXM6POLS

    !.....  Check if this is a valid pollutant/process combo
                    IF( M6POL2EF( I,J ) == -1 ) CYCLE

                    K = K + 1

    !.....  Pull process name and description from include file
                    EFSNAM( K ) = M6PROCS( I ) // ETJOIN
                    EFSDSC( K ) = 'EFs for ' // M6PRCDSC( I )

    !.....  If pollutant is HC, append specified volatile pollutant
    !                           name; otherwise, use name and desc from include file
                    IF( M6POLS( J ) == 'HC' ) THEN
                        EFSNAM( K ) = TRIM( EFSNAM( K ) ) // VOLNAM
                        EFSDSC( K ) = TRIM( EFSDSC( K ) ) // ' ' // VOLNAM
                    ELSE
                        EFSNAM( K ) = TRIM( EFSNAM( K ) ) // M6POLS( J )
                        EFSDSC( K ) = TRIM( EFSDSC( K ) ) // ' ' // M6POLDSC( J )
                    END IF

    !.....  Store units from include file
                    EFSUNT( K ) = M6UNIT

                END DO      ! pollutant loop
            END DO      ! emission process loop

    !.....  If file available
        ELSE
    ! note: add here

        END IF

    END SELECT

    FIRSTIME = .FALSE.

    RETURN

END SUBROUTINE EFSETUP
