
SUBROUTINE CHKISIZ( FILNAM, FILDESC, COMPDESC, NSRC, STATUS )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      Checks the number of sources and sets error status (0 okay, 1 bad)
!
!  PRECONDITIONS REQUIRED:
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

    IMPLICIT NONE

!...........   INCLUDE FILES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN   ) :: FILNAM         ! Name of file being checked
    CHARACTER(*), INTENT (IN   ) :: FILDESC        ! Description of file being checked
    CHARACTER(*), INTENT (IN   ) :: COMPDESC       ! Description of comparison value
    INTEGER     , INTENT (IN   ) :: NSRC           ! Number of sources comparing against
    INTEGER,      INTENT (  OUT) :: STATUS         ! Exit status

!...........   Other local variables

    CHARACTER(256)  MESG

!***********************************************************************
!   begin body of subroutine CHKISIZ

    STATUS = 0

    IF( .NOT. DESC3( FILNAM ) ) THEN
        STATUS = 1
        MESG = 'Could not read description for "' //&
        &      TRIM( FILNAM ) // '"'
        CALL M3MSG2( MESG )
    ENDIF

!.............  Compare the number of sources to the NSRC value
    IF( NCOLS3D .NE. NSRC ) THEN
        STATUS = 1
        WRITE( MESG,94010 )&
        &       'Number of inventory sources mismatch.' //&
        &       CRLF() // BLANK16 //&
        &       COMPDESC // ':', NSRC,&
        &       CRLF() // BLANK16 //&
        &       FILDESC // ':', NCOLS3D

        CALL M3MSG2( MESG )
    ENDIF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CHKISIZ

