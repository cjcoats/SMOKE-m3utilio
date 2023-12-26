
SUBROUTINE WRDAYMSG( JDATE, MESG )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      Writes a text message to stdout and log file about which day is
!      being processed
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!      Created 9/99 by M. Houyoux
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

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: JDATE    ! Julian date
    CHARACTER(*), INTENT(OUT) :: MESG     ! message buffer

!...........   Local variables
    INTEGER         DAY          !  day of week number

!***********************************************************************
!   begin body of subroutine WRDAYMSG

    DAY = WKDAY( JDATE )

    MESG= 'Processing '// TRIM( DAYS( DAY ) )// ' '// MMDDYY( JDATE )
    CALL M3MSG2( MESG )

    RETURN

END SUBROUTINE WRDAYMSG

