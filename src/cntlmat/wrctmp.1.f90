
SUBROUTINE WRCTMP( IDEV, POLID, IDX, VIDX, LOUTANY )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine writes control packet data table indices
!      to a temporary file.
!
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
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
!***************************************************************************

!.........  MODULES for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NIPPA

    IMPLICIT NONE

!...........   SUBROUTINE ARGUMENTS:

    INTEGER     , INTENT (IN) :: IDEV           ! logical file name
    INTEGER     , INTENT (IN) :: POLID          ! pollutant number
    INTEGER     , INTENT (IN) :: IDX ( NSRC )   ! index to data tables
    INTEGER     , INTENT (IN) :: VIDX( NIPPA )  ! pollutant/act flags
    LOGICAL     , INTENT(OUT) :: LOUTANY        ! true: at least one pollutant output

!...........   Other local variables

    INTEGER   S   ! indices

    LOGICAL, SAVE :: FIRSTTIME = .TRUE.

!***********************************************************************
!   Begin body of subroutine WRCTMP

    IF ( FIRSTTIME ) THEN
        LOUTANY   =  .FALSE.
        FIRSTTIME = .FALSE.
    ENDIF

!.............. Write indices to control factor packets to a temporary file
!               for only those pollutants that have controls
    IF ( VIDX( POLID ) .EQ. 1 ) THEN
        DO S = 1, NSRC

            WRITE( IDEV, '(I8)' ) IDX( S )

        END DO   ! end source loop
        LOUTANY = .TRUE.

    END IF

    RETURN

END SUBROUTINE WRCTMP
