
SUBROUTINE CHECKMEM( MSTATUS, ONVAR, CALLER )

!***********************************************************************
!  subroutine body starts at line  105
!
!  DESCRIPTION:
!       Reports an error and exits if memory status flag is non-zero.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created  ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***********************************************************************
!
! Project Title: EDSS Tools Library
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

!...........   ARGUMENTS and their descriptions:

    INTEGER     , INTENT( IN ) ::  MSTATUS !  ALLOCATE function exit status
    CHARACTER(*), INTENT( IN ) ::  ONVAR   !  Variable name(s) of previous ALLOCATE statement
    CHARACTER(*), INTENT( IN ) ::  CALLER  !  Name of calling program

!...........   Local variables

    CHARACTER(256)  MESG
    CHARACTER(32)   PNAME

!***********************************************************************
!   begin body of function CHECKMEM

!.........  Abort if memory status is non-zero

    IF( MSTATUS .GT. 0 ) THEN
        PNAME = TRIM( CALLER ) // ':' // 'CHECKMEM'
        WRITE( MESG, '( 3 A, I10 )' )&
        &       'Failure allocating memory for "', TRIM( ONVAR ),&
        &       '":  STATUS=', MSTATUS
        CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

END SUBROUTINE CHECKMEM

