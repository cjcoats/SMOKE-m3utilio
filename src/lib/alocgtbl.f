
        SUBROUTINE ALOCGTBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the gridding
C      cross-reference tables that contain the gridding surrogate code numbers,
C      and it initializes these to missing.  The subroutine argument is
C      an array that contains the dimensions for each of the different groups
C      of the cross-reference.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C     Created 3/99 by M. Houyoux
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
        USE MODXREF, ONLY: ISRG01, ISRG02, ISRG03, ISRG04, ISRG05,
     &                     ISRG06, ISRG07, ISRG08, ISRG09

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       IOS              ! i/o status

        CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCGTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCGTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( ISRG02 ) ) THEN

            DEALLOCATE( ISRG02, ISRG03, ISRG04, ISRG05,
     &                  ISRG06, ISRG07, ISRG08, ISRG09 )

        END IF

        ALLOCATE( ISRG02( ICSIZE( 2 ) ),
     &            ISRG03( ICSIZE( 3 ) ),
     &            ISRG04( ICSIZE( 4 ) ),
     &            ISRG05( ICSIZE( 5 ) ),
     &            ISRG06( ICSIZE( 6 ) ),
     &            ISRG07( ICSIZE( 7 ) ),
     &            ISRG08( ICSIZE( 8 ) ),
     &            ISRG09( ICSIZE( 9 ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ISRG02...ISRG09', PROGNAME )

        ISRG01 = IMISS3     ! is scalar
        ISRG02 = IMISS3     ! array
        ISRG03 = IMISS3     ! array
        ISRG04 = IMISS3     ! array
        ISRG05 = IMISS3     ! array
        ISRG06 = IMISS3     ! array
        ISRG07 = IMISS3     ! array
        ISRG08 = IMISS3     ! array
        ISRG09 = IMISS3     ! array

        RETURN

        END SUBROUTINE ALOCGTBL
