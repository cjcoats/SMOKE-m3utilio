
SUBROUTINE OPENPHYS( CALLER, FNAME, STATUS, PHYSNAME, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine sets the physical file name to a logical file
    !      name and opens the i/o api FILESET file.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ??
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !****************************************************************************/
    !
    ! Project Title: EDSS Tools Library
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
    USE MODFILESET, ONLY:  ALLFILES

    IMPLICIT NONE

    !......   INCLUDES
    INCLUDE 'IOCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !......  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT(IN   ) :: CALLER        ! calling routine
    CHARACTER(*), INTENT(IN   ) :: FNAME         ! logical file name
    LOGICAL     , INTENT(IN   ) :: STATUS        ! OPENSET status to use
    CHARACTER(*), INTENT(IN   ) :: PHYSNAME      ! physical file name
    LOGICAL     , INTENT(INOUT) :: EFLAG         ! true: error occurred

    !......  Other local variables
    CHARACTER(256)          MESG        !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENPHYS'     ! program name

    !***********************************************************************
    !   begin body of subroutine OPENPHYS

    !......  Set environment variable
    IF( .NOT. SETENVVAR( FNAME, PHYSNAME ) ) THEN
        MESG = 'ERROR: Could not set logical name for file '    &
               // CRLF() // BLANK10 // TRIM( PHYSNAME )
        CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

        !......  Open file for current pollutant
    ELSE IF( .NOT. OPENSET( FNAME, STATUS, CALLER ) ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: Could not open file: ' //            &
                   CRLF() // BLANK10 // TRIM( PHYSNAME )
        CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

        !......  Get description for file
    ELSE IF( .NOT. DESCSET( FNAME, ALLFILES ) ) THEN
        EFLAG = .TRUE.
        MESG = 'Could not get description of file: ' //     &
               CRLF()//BLANK10// TRIM( PHYSNAME )
        CALL M3EXIT( CALLER, 0, 0, MESG, 2 )

    END IF

    RETURN

END SUBROUTINE OPENPHYS
