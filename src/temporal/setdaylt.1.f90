
SUBROUTINE SETDAYLT

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine allocates and populates the by-source array that
!      indicates which sources get daylight time and which do not.
!
!  PRECONDITIONS REQUIRED:
!      Number of sources set
!      Region codes read from inventory file.
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 10/2000 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!************************************************************************
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

!...........   MODULES for public variables
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: FLTRDAYL, CIFIP

!.........  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTY, CNTYCOD, USEDAYLT

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

!...........   Other local variables

    INTEGER         J, S        ! counters and indices

    INTEGER         FLTR        ! tmp filter value
    INTEGER         IOS         ! i/o status

    CHARACTER(FIPLEN3) CFIP     ! tmp region code
    CHARACTER(FIPLEN3) PFIP     ! region code from previous iteration
    CHARACTER(300)  MESG        ! message buffer

    CHARACTER(16) :: PROGNAME = 'SETDAYLT' ! program name

!***********************************************************************
!   begin body of subroutine SETDAYLT

!.........  Allocate memory for daylight time filter
    ALLOCATE( FLTRDAYL( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FLTRDAYL', PROGNAME )

!.........  Initialize all sources as using daylight time
    FLTRDAYL = 1    ! array

!.........  Use region codes daylight time assignment information (from MODSTCY)
!           to set sources that do not use daylight time
!.........  No error if region code is not found in list, because if it is
!           not in the list, then the daylight status can remain in use
    PFIP = ' '
    DO S = 1, NSRC

        CFIP = CIFIP( S )

        IF( CFIP .NE. PFIP ) THEN

            J = FINDC( CFIP, NCOUNTY, CNTYCOD )

            FLTR = 1
            IF( J .GT. 0 ) THEN

                IF( .NOT. USEDAYLT( J ) ) FLTR = 0

            END IF

            PFIP = CFIP

        END IF

        FLTRDAYL( S ) = FLTR

    END DO

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94000 FORMAT( I2.2 )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE SETDAYLT

