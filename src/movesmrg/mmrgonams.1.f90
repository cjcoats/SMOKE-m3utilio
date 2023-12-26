
SUBROUTINE MRGONAMS

!***********************************************************************
!  subroutine MRGONAMS body starts at line
!
!  DESCRIPTION:
!      The purpose of this subroutine is to set the default names of the
!      output files. It sets them regardless of whether they will be used
!      or not.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 2/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: MREPNAME, MONAME, SGINLNNAME, SUBOUTNAME,&
    &                    NUNITS, GRDUNIT, NGRPS, IUGRPNUM, SUBSECFLAG

    IMPLICIT NONE

    INTEGER         I, J, N, K    ! indices and counters

    INTEGER         IOS           ! tmp I/O status

    LOGICAL      :: KFLAG = .FALSE.          ! true: simplied output file names
    LOGICAL      :: MOLEFLAG = .FALSE.       ! true: outputting moles

    CHARACTER(300)  MESG    ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGONAMS' ! program name

!***********************************************************************
!   begin body of subroutine MRGONAMS

!.........  Retrieve variable to indicate whether to use annual or average day data
    MESG = 'Use customized MOVESMRG output file names'
    KFLAG = ENVYN( 'MOVESMRG_CUSTOM_OUTPUT', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "MOVESMRG_CUSTOM_OUTPUT"', 2 )
    END IF

!.........  Set default output file name(s) depending on inputs.  Set both
!           report file(s) and gridded file(s)...

!.........  Initialize - everything will be gridded
    N = 1
    IF( SUBSECFLAG ) N = NGRPS
    ALLOCATE( SUBOUTNAME( N ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SUBOUTNAME', PROGNAME )
    SUBOUTNAME = 'SUBOUT'

    DO K = 1, N

        IF( KFLAG ) THEN

            MREPNAME = 'REPMG'
            MONAME = 'MOUT'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IUGRPNUM( K )

        ELSE

            MREPNAME = 'REPMG'
            MONAME = 'MG'
            SGINLNNAME = 'SGINLN'
            IF( SUBSECFLAG ) WRITE( SUBOUTNAME( K ), '(A,I2.2)' ) 'SUBOUT', IUGRPNUM( K )

            CALL TRIM_AND_CONCAT( MREPNAME,   'T' )
            CALL TRIM_AND_CONCAT( MONAME,     'T' )
            CALL TRIM_AND_CONCAT( SGINLNNAME, 'T' )
            CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'T' )

            CALL TRIM_AND_CONCAT( MREPNAME,   'S' )
            CALL TRIM_AND_CONCAT( MONAME,     'S' )
            CALL TRIM_AND_CONCAT( SGINLNNAME, 'S' )
            CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), 'S' )

!.............  Now append mass or mole for I/O API files, depending on which
!               inputs were used

!.............  Set flag if any of the output species are mole-based
            DO I = 1, NUNITS

                J = INDEX( GRDUNIT( I ), 'mole' )
                IF( J .GT. 0 ) THEN
                    MOLEFLAG = .TRUE.
                    EXIT
                END IF

            END DO

!.............  Get output file names depending on if there are moles in units
            IF( MOLEFLAG ) THEN

                CALL TRIM_AND_CONCAT( MREPNAME,   '_L' )
                CALL TRIM_AND_CONCAT( MONAME,     '_L' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, '_L' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_L' )

            ELSE

                CALL TRIM_AND_CONCAT( MREPNAME,   '_S' )
                CALL TRIM_AND_CONCAT( MONAME,     '_S' )
                CALL TRIM_AND_CONCAT( SGINLNNAME, '_S' )
                CALL TRIM_AND_CONCAT( SUBOUTNAME( K ), '_S' )

            END IF

        END IF

    END DO

    RETURN

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram trims and concatonates the first
!               string with the second string
    SUBROUTINE TRIM_AND_CONCAT( PART1, PART2 )

!.............  Subprogram arguments
        CHARACTER(*)  PART1
        CHARACTER(*)  PART2

!----------------------------------------------------------------------

        PART1 = TRIM( PART1 ) // PART2

    END SUBROUTINE TRIM_AND_CONCAT

END SUBROUTINE MRGONAMS
