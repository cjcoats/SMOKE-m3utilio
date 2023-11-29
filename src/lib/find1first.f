
        INTEGER FUNCTION FIND1FIRST( KEY, N, LIST )

C***********************************************************************
C  function body starts at line 55
C
C  DESCRIPTION:
C       Returns first instance of INTEGER / REAL / CHARACTER KEY in LIST
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED: FIND1
C
C  REVISION  HISTORY:
C       10/01: Created by C. Seppanen
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***********************************************************************

        USE M3UTILIO
        IMPLICIT NONE

C.........  Function arguments
        INTEGER, INTENT (IN) :: KEY        ! key to search for
        INTEGER, INTENT (IN) :: N          ! number of entries in LIST
        INTEGER, INTENT (IN) :: LIST( N )  ! table to be searched

C.........  Local function variables
        INTEGER INDEX

C***********************************************************************
C   begin body of function FIND1FIRST

C.........  Use FIND1 to get location of key
        INDEX = FIND1( KEY, N, LIST )

C.........  If the key is found, search backward until the first entry is reached
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

C..........................................................

        INTEGER FUNCTION FINDR1FIRST( KEY, N, LIST )

        USE M3UTILIO
        IMPLICIT NONE

C.........  Function arguments
        REAL,    INTENT (IN) :: KEY        ! key to search for
        INTEGER, INTENT (IN) :: N          ! number of entries in LIST
        REAL,    INTENT (IN) :: LIST( N )  ! table to be searched

C.........  Local function variables
        INTEGER INDEX

C***********************************************************************
C   begin body of function FINDR1FIRST

C.........  Use FINDR1 to get location of key
        INDEX = FINDR1( KEY, N, LIST )

C.........  If the key is found, search backward until the first entry is reached
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

C..........................................................

        INTEGER FUNCTION FINDCFIRST( KEY, N, LIST )

        USE M3UTILIO
        IMPLICIT NONE

C.........  Function arguments
        CHARACTER(*), INTENT (IN) :: KEY        ! key to search for
        INTEGER,      INTENT (IN) :: N          ! number of entries in LIST
        CHARACTER(*), INTENT (IN) :: LIST( N )  ! table to be searched

C.........  Local function variables
        INTEGER INDEX

C***********************************************************************
C   begin body of function FINDCFIRST

C.........  Use FINDC to get location of key
        INDEX = FINDC( KEY, N, LIST )

C.........  If the key is found, search backward until the first entry is reached
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

