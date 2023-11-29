
        CHARACTER(16) FUNCTION VERCHAR( BUFFER )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      This function returns the version number from a string filled in
C      by SCCS operations
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version 11/2023 by CJC
C
C**************************************************************************
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

        IMPLICIT NONE

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT( IN ) :: BUFFER         ! string FIPS code

C...........   Other local variables
        INTEGER         J

C***********************************************************************
C   begin body of function VERCHAR

        VERCHAR = ADJUSTL( BUFFER )

        IF( VERCHAR( 1:1 ) .EQ. '%' ) THEN
            VERCHAR = 'Dev'
        ELSE
            J  = INDEX  ( VERCHAR, ' ' )
            VERCHAR = ADJUSTL( VERCHAR( J+1: ) )
        ENDIF

        RETURN

        END FUNCTION VERCHAR

