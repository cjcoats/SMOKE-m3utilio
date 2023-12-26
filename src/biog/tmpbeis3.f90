
PROGRAM TMPBEIS3

    !***********************************************************************
    !  program body starts at line  187
    !
    !  DESCRIPTION:
    !       Computes hourly time stepped gridded biogenic emissions using
    !       normalized gridded emissions from Normbeis3 and postprocessed MM5
    !       meteorology.
    !
    !  PRECONDITIONS REQUIRED:
    !       Postprocessed MM5 meteorology that contains temperature,
    !       solar radiation, and pressure data.
    !       Normalized gridded emissions B3GRD from Normbeis3
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !       HRBIO, PREBMET
    !
    !  REVISION  HISTORY:
    !      3/01: Prototype by Jeff Vukovich
    !            Tested only on 36km Lambert domain
    !            Summer/winter switch file option not tested
    !      8/04: Incorporated BEIS v3.12 by C. Seppanen
    !      4/06: changed to BEIS3.13 by G. Pouliot
    !      3/08: changed to BEIS3.14 by G. Pouliot
    !      ?/14: changed to BEIS3.60 by G. Pouliot
    !      7/15  Changed from BEIS3.60 to BEIS3.61 by Baek
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
    !***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......  PARAMETERs
    CHARACTER(16), PARAMETER :: PNAME = 'TMPBEIS3'       !  program name
    CHARACTER(50), PARAMETER :: CVSW  = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.......  Logical names and unit numbers
    INTEGER         LDEV        !  unit number for log device

    !.......  Other local variables
    INTEGER         IOS         !  I/O status

    CHARACTER(16)   BEISVER     !  version of BEIS3 to use
    CHARACTER(300)  MESG        !  message buffer for M3EXIT()

    !***********************************************************************
    !   begin body of program TMPBEIS3

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PNAME )

    !.......  Get the BEIS3 model version to use
    MESG = 'Version of BEIS3 to use'
    CALL ENVSTR( 'BEIS_VERSION', MESG, '3.7', BEISVER, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BEIS_VERSION"', 2 )
    END IF

    SELECT CASE( BEISVER )
      CASE( '3.7' )
        CALL TMPBEIS360( CVSW )
      CASE( '3.61' )
        CALL TMPBEIS360( CVSW )
      CASE( '3.14' )
        CALL TMPBEIS314( CVSW )
      CASE DEFAULT
        MESG = 'ERROR: Unrecognized BEIS_VERSION setting; valid ' //&
               'settings are 3.14, 3.61 and 3.7'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END SELECT

    !.......   End of program
    CALL M3EXIT( PNAME, 0, 0, 'Completion of ' // PNAME, 0 )

END PROGRAM TMPBEIS3

