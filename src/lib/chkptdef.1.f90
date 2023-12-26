
SUBROUTINE CHKPTDEF( CHKNCHAR, CHKJSCC )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine compares the source definition that is input through
!      the subroutine arguments to the defition from the inventory file
!      that is stored in MODINFO.  This is used to check to ensure the
!      cross-reference files and the inventory use the same definition of
!      point sources.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Create 3/99 by M. Houyoux
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

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NCHARS, CATDESC, JSCC

    IMPLICIT NONE

!.........  INCLUDE FILES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: CHKNCHAR     ! number of chars in source def
    INTEGER, INTENT (IN) :: CHKJSCC      ! position of SCC in source def

!.........  Other local variables
    LOGICAL      :: EFLAG = .FALSE.      ! true: error found

    CHARACTER(300)  MESG

    CHARACTER(16) :: PROGNAME = 'CHKPTDEF' ! program name

!***********************************************************************
!   begin body of subroutine CHKPTDEF

!.........  Compare the input total number of sources with the inventory
    IF( CHKNCHAR .NE. NCHARS ) THEN

        EFLAG = .TRUE.
        WRITE( MESG,94010 )&
        &       'ERROR: Source definition mismatch. Number of ' //&
        &       'source definition' // CRLF() // BLANK10 //&
        &       'characteristics in file: ', CHKNCHAR,&
        &       ', but in ' // CATDESC //&
        &       ' source inventory file: ', NCHARS
        CALL M3MSG2( MESG )

    END IF

!.........  Compare the input position of SCCs with the inventory
    IF( CHKJSCC .NE. JSCC ) THEN

        EFLAG = .TRUE.
        WRITE( MESG,94010 )&
        &       'ERROR: Source definition mismatch. Position of ' //&
        &       'SCC in file source definition: ', CHKJSCC, ', ' //&
        &       CRLF() // BLANK10 // 'but in ' // CATDESC //&
        &       ' source inventory file: ', JSCC
        CALL M3MSG2( MESG )

    END IF

!.........  If there is an error, exit
    IF( EFLAG ) THEN

        MESG = 'Inconsistent input files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CHKPTDEF

