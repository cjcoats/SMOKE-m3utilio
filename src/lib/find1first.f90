
INTEGER FUNCTION FIND1FIRST( KEY, N, LIST )

!***********************************************************************
!  function body starts at line 55
!
!  DESCRIPTION:
!       Returns first instance of INTEGER / REAL / CHARACTER KEY in LIST
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED: FIND1
!
!  REVISION  HISTORY:
!       10/01: Created by C. Seppanen
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
!***********************************************************************

    USE M3UTILIO
    IMPLICIT NONE

!.........  Function arguments
    INTEGER, INTENT (IN) :: KEY        ! key to search for
    INTEGER, INTENT (IN) :: N          ! number of entries in LIST
    INTEGER, INTENT (IN) :: LIST( N )  ! table to be searched

!.........  Local function variables
    INTEGER INDEX

!***********************************************************************
!   begin body of function FIND1FIRST

!.........  Use FIND1 to get location of key
    INDEX = FIND1( KEY, N, LIST )

!.........  If the key is found, search backward until the first entry is reached
    IF( INDEX > 0 ) THEN
        DO
            IF( INDEX < 1 ) EXIT
            IF( LIST( INDEX ) /= KEY ) EXIT
            INDEX = INDEX - 1
        END DO

        INDEX = INDEX + 1
    END IF

    FIND1FIRST = INDEX

    RETURN

END FUNCTION FIND1FIRST

!..........................................................

INTEGER FUNCTION FINDR1FIRST( KEY, N, LIST )

    USE M3UTILIO
    IMPLICIT NONE

!.........  Function arguments
    REAL,    INTENT (IN) :: KEY        ! key to search for
    INTEGER, INTENT (IN) :: N          ! number of entries in LIST
    REAL,    INTENT (IN) :: LIST( N )  ! table to be searched

!.........  Local function variables
    INTEGER INDEX

!***********************************************************************
!   begin body of function FINDR1FIRST

!.........  Use FINDR1 to get location of key
    INDEX = FINDR1( KEY, N, LIST )

!.........  If the key is found, search backward until the first entry is reached
    IF( INDEX > 0 ) THEN
        DO
            IF( INDEX < 1 ) EXIT
            IF( LIST( INDEX ) /= KEY ) EXIT
            INDEX = INDEX - 1
        END DO

        INDEX = INDEX + 1
    END IF

    FINDR1FIRST = INDEX

    RETURN

END FUNCTION FINDR1FIRST

!..........................................................

INTEGER FUNCTION FINDCFIRST( KEY, N, LIST )

    USE M3UTILIO
    IMPLICIT NONE

!.........  Function arguments
    CHARACTER(*), INTENT (IN) :: KEY        ! key to search for
    INTEGER,      INTENT (IN) :: N          ! number of entries in LIST
    CHARACTER(*), INTENT (IN) :: LIST( N )  ! table to be searched

!.........  Local function variables
    INTEGER INDEX

!***********************************************************************
!   begin body of function FINDCFIRST

!.........  Use FINDC to get location of key
    INDEX = FINDC( KEY, N, LIST )

!.........  If the key is found, search backward until the first entry is reached
    IF( INDEX > 0 ) THEN
        DO
            IF( INDEX < 1 ) EXIT
            IF( LIST( INDEX ) /= KEY ) EXIT
            INDEX = INDEX - 1
        END DO

        INDEX = INDEX + 1
    END IF

    FINDCFIRST = INDEX

    RETURN

END FUNCTION FINDCFIRST

