
       SUBROUTINE ALOCPTBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the speed
C      cross-reference tables that contain the speed profile numbers, and
C      it initializes these to missing.  The subroutine argument is
C      an array that contains the dimensions for each of the different
C      groups of the cross-reference.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 1/03 by C. Seppanen
C
C****************************************************************************/
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
C***************************************************************************

        USE M3UTILIO

C...........   This module is for cross reference tables
        USE MODXREF, ONLY: ISPD01, ISPD02, ISPD03, ISPD04, ISPD05,
     &                     ISPD06, ISPD07, ISPD08, ISPD09

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       IOS              ! i/o status

        CHARACTER(16) :: PROGNAME = 'ALOCPTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCPTBL

C.........  First deallocate if these have previously been allocated
        IF( ALLOCATED( ISPD02 ) ) THEN

            DEALLOCATE( ISPD02, ISPD03, ISPD04, ISPD05 )
            DEALLOCATE( ISPD06, ISPD07, ISPD08, ISPD09 )

        END IF
                      ! SCC=left, FIP=0
        ALLOCATE( ISPD02( ICSIZE( 2 ) ),
     &            ISPD03( ICSIZE( 3 ) ),
     &            ISPD04( ICSIZE( 4 ) ),
     &            ISPD05( ICSIZE( 5 ) ),
     &            ISPD06( ICSIZE( 6 ) ),
     &            ISPD07( ICSIZE( 7 ) ),
     &            ISPD08( ICSIZE( 8 ) ),
     &            ISPD09( ICSIZE( 9 ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISPD09', PROGNAME )

        ISPD01 = IMISS3     ! scalar
        ISPD02 = IMISS3     ! array
        ISPD03 = IMISS3     ! array
        ISPD04 = IMISS3     ! array
        ISPD05 = IMISS3     ! array
        ISPD06 = IMISS3     ! array
        ISPD07 = IMISS3     ! array
        ISPD08 = IMISS3     ! array
        ISPD09 = IMISS3     ! array

        RETURN

        END SUBROUTINE ALOCPTBL
