
        SUBROUTINE ALOCMTBL( ICSIZE )

C***********************************************************************
C  subroutine body starts at line
C
C  DESCRIPTION:
C      This subroutine allocates memory for the portion of the VMT mix
C      table that contain the various mobile-source characteristics. The
C      subroutine argument is an array that contains the dimensions for each
C      of the different groups of the cross-reference.
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Created 2/2000 by M. Houyoux
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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
        USE MODXREF, ONLY: IMVS01, IMVS02, IMVS03, IMVS04, IMVS05,
     &                     IMVS06, IMVS07, IMVS08, IMVS09, IMVS10,
     &                     IMVS11, IMVS12

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        INTEGER, INTENT(IN) :: ICSIZE( * )  ! size of x-ref groups

C...........   Other local variables
        INTEGER       IOS              ! i/o status

        CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCMTBL' ! program name

C***********************************************************************
C   begin body of subroutine ALOCMTBL

C.........  First deallocate if these have previously been allocated
        IF ( ALLOCATED( IMVS02 ) ) THEN

            DEALLOCATE( IMVS02, IMVS03, IMVS04, IMVS05 )
            DEALLOCATE( IMVS06, IMVS07, IMVS08, IMVS09 )

        END IF

        ALLOCATE( IMVS02( ICSIZE( 2 ) ),
     &            IMVS03( ICSIZE( 3 ) ),
     &            IMVS04( ICSIZE( 4 ) ),
     &            IMVS05( ICSIZE( 5 ) ),
     &            IMVS06( ICSIZE( 6 ) ),
     &            IMVS07( ICSIZE( 7 ) ),
     &            IMVS08( ICSIZE( 8 ) ),
     &            IMVS09( ICSIZE( 9 ) ),
     &            IMVS10( ICSIZE(10 ) ),
     &            IMVS11( ICSIZE(11 ) ),
     &            IMVS12( ICSIZE(12 ) ), STAT=IOS )
        CALL CHECKMEM( IOS, 'IMVS02...IMVS12', PROGNAME )

        IMVS02 = IMISS3    ! array
        IMVS03 = IMISS3    ! array
        IMVS04 = IMISS3    ! array
        IMVS05 = IMISS3    ! array
        IMVS06 = IMISS3    ! array
        IMVS07 = IMISS3    ! array
        IMVS08 = IMISS3    ! array
        IMVS09 = IMISS3    ! array
        IMVS09 = IMISS3    ! array
        IMVS09 = IMISS3    ! array
        IMVS09 = IMISS3    ! array

       RETURN

        END SUBROUTINE ALOCMTBL
