
SUBROUTINE CHKSRCNO( CATDESC, FILNAM, NTEST, NSRC, EFLAG )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine compares the number of rows in an I/O API header
!      (NROWS3D) with the input number of sources for the purpose of checking
!      the number of sources between two files
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
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

!.........  INCLUDE FILES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    CHARACTER(*)    CATDESC    ! source category description
    CHARACTER(*)    FILNAM     ! name of file being checked
    INTEGER         NTEST      ! number of sources to check
    INTEGER         NSRC       ! number of sources to compare against
    LOGICAL         EFLAG      ! error flag

!.........  Other local variables

    CHARACTER(300)  MESG

    CHARACTER(16) :: PROGNAME = 'CHKSRCNO' ! program name

!***********************************************************************
!   begin body of subroutine CHKSRCNO

!.........  Check the number of sources
    IF( NTEST .NE. NSRC ) THEN

        EFLAG = .TRUE.
        WRITE( MESG,94010 )&
        &       'ERROR: Dimension mismatch. Source number in ' //&
        &       FILNAM // ' file: ', NTEST, ', ' // CRLF() //&
        &       BLANK10 // 'but in ' // CATDESC //&
        &       ' source inventory file: ', NSRC
        CALL M3MSG2( MESG )

    END IF

    RETURN
!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CHKSRCNO

