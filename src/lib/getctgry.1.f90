
SUBROUTINE GETCTGRY

!***********************************************************************
!  subroutine body starts at line 85
!
!  DESCRIPTION:
!     This subroutine retrieves the SMK_SOURCE environment variable to
!     determine what source category is being processed.  It returns
!     the first letter, category, and category description based on that
!     environment variable setting.
!
!  PRECONDITIONS REQUIRED:
!
!
!  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
!
!  REVISION  HISTORY:
!       Created 3/99 by M Houyoux
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

!...........   MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CRL, CATEGORY, CATDESC, CATLEN

    IMPLICIT NONE

!...........   INCLUDES:
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   LOCAL PARAMETER:
    INTEGER  , PARAMETER :: NLIST   = 5
    CHARACTER, PARAMETER :: LETLIST( NLIST ) =&
    &                          ( / 'A', 'B', 'M', 'P', 'E' / )

    CHARACTER(6), PARAMETER :: CATTYPE( NLIST ) =&
    &                          ( / 'AREA  ', 'BIOGEN', 'MOBILE',&
    &                              'POINT ', 'EVERY ' / )

    CHARACTER(6), PARAMETER :: CATLDSC( NLIST ) =&
    &                          ( / 'Area  ', 'Biogen', 'Mobile',&
    &                              'Point ', 'Every ' / )

!...........   LOCAL VARIABLES their descriptions:

    INTEGER      :: IOS = 0    ! i/o status
    INTEGER         J          ! index

    CHARACTER(16) :: STRBLK = ' '
    CHARACTER(16)    STRVAL
    CHARACTER(200)   MESG       ! Message buffer

    CHARACTER(16) :: PNAME = 'GETCTGRY'    ! Program name

!***********************************************************************
!   begin body of subroutine GETCTGRY

!.........  Retrieve environment variable that indicates the source of interest
    MESG = 'Control for which source category for controls'
    CALL ENVSTR( 'SMK_SOURCE', MESG, STRBLK, STRVAL, IOS )
    IF ( IOS .NE. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "SMK_SOURCE"', 2 )
    END IF

    CRL = ADJUSTL( STRVAL )
    J = INDEX1( CRL, NLIST, LETLIST )

    IF( J .LE. 0 ) THEN

        MESG = 'ERROR: Do not recognize SMK_SOURCE environment '//&
        &       'variable setting "' // CRL // '"'
        CALL M3MSG2( MESG )

        MESG = 'Problem with environment variable settings'
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )

    ELSE

        CATEGORY = CATTYPE( J )
        CATDESC  = CATLDSC( J )

        CATLEN   = LEN_TRIM( CATEGORY )

    ENDIF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I6, :, 2X ) )

END SUBROUTINE GETCTGRY
