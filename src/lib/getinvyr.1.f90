
INTEGER FUNCTION GETINVYR( LINE )

!***********************************************************************
!  function body starts at line
!
!  DESCRIPTION:
!      This function returns the inventory year from a string if the
!      string contains the INVYEAR packet.
!
!  PRECONDITIONS REQUIRED:

!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created by M. Houyoux 12/98
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
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

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*) LINE    !  description of source category

    INTEGER         INY
    INTEGER         L1, L2

    CHARACTER(300)  MESG        !  message buffer

    CHARACTER(16) :: PROGNAME = 'GETINVYR' ! program name

!***********************************************************************
!   begin body of function GETINVYR

    L1 = INDEX( LINE, 'INVYEAR' )
    L2 = LEN_TRIM( LINE )

!.........  Process INVYEAR packet
    IF( L1 .GT. 0 ) THEN

        INY = STR2INT( LINE( L1+7:L2 ) )

        IF( INY .LE. 0 ) THEN

            MESG = 'Incorrectly set year using ' //&
            &       'INVYEAR packet.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSEIF( INY .LT. 1970 ) THEN

            MESG = 'INVYEAR packet has set 4-digit year ' //&
            &       'below 1970 minimum in PTINV file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

        GETINVYR = INY

    ELSE

        GETINVYR = IMISS3

    ENDIF


    RETURN

END FUNCTION GETINVYR
