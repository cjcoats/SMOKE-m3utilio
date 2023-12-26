
PROGRAM NORMBEIS3

!***********************************************************************
!
!  DESCRIPTION:  Produces normalized biogenic emissions for use with
!                SMOKE-BEIS versions 3.14 and 3.6
!
!  SUBROUTINES AND FUNCTIONS CALLED: Calls Normbeis314 or Normbeis360
!
!  REVISION  HISTORY:
!       3/00 Prototype, Jeff Vukovich
!       8/04 Integrated v3.12, C. Seppanen
!       4/06 Changed Beis3.12 to BEIS3.13 G. Pouliot
!       3/08 Changed Beis3.13 to BEIS3.14 G. Pouliot
!       ?/14 Changed Beis3.14 to BEIS3.60 G. Pouliot
!       7/15 Changed BEIS3.60 to BEIS3.61 by Baek
!       8/20 Added BEIS3.70 with BELD5 and BEISFACS
!       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
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
!***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!.........  Local parameters
    CHARACTER(50), PARAMETER :: CVSW = '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!.........  Logical names and unit numbers
    INTEGER         LDEV    !  unit number for log device

!.........  Other local variables
    INTEGER         IOS     !  I/O status

    CHARACTER(16)   BEISVER !  version of BEIS3 to use
    CHARACTER(300)  MESG    !  message buffer

    CHARACTER(16), PARAMETER :: PNAME = 'NORMBEIS3'   !  program name

!***********************************************************************
!   begin body of program NORMBEIS3

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, CVSW, PNAME )

!.........  Get the BEIS3 model version to use
    MESG = 'Version of BEIS3 to use'
    CALL ENVSTR( 'BEIS_VERSION', MESG, '3.7', BEISVER, IOS )
    IF ( IOS .NE. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "BEIS_VERSION"', 2 )
    END IF

    SELECT CASE( BEISVER )
      CASE( '3.7' )
        CALL NORMBEIS370( CVSW )
      CASE( '3.61' )
        CALL NORMBEIS360( CVSW )
      CASE( '3.14' )
        CALL NORMBEIS312( CVSW )
      CASE DEFAULT
        MESG = 'ERROR: Unrecognized BEIS_VERSION setting; valid ' //&
        &       'settings are 3.14, 3.61, and 3.7'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END SELECT

!.........  End of program
    CALL M3EXIT( PNAME, 0, 0, ' ', 0 )

END PROGRAM NORMBEIS3

