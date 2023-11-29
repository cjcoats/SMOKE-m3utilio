
        LOGICAL FUNCTION USEEXPGEO()

C***********************************************************************
C  function USEEXPGEO body starts at line 50
C
C  DESCRIPTION:
C      The first time this function is called it will get the value
C      of the environment variable USE_EXP_GEO_CODES to indicate if
C      expanded geographic codes should be used. Subsequent calls
C      will return the saved setting.
C
C  REVISION  HISTORY:
C       Created 12/13 by C. Seppanen
C       Version 11/2023 by CJC:  USE M3UTILIO and related changes
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
C****************************************************************************

        USE M3UTILIO

C...........   Other local variables
        INTEGER        IOS    ! i/o status

        LOGICAL, SAVE :: FIRSTTIME  = .TRUE.
        LOGICAL, SAVE :: EXPGEOFLAG = .FALSE.

        CHARACTER(16), PARAMETER :: PNAME = 'USEEXPGEO' ! program name

C****************************************************************************
C   begin body of function USEEXPGEO

        IF ( FIRSTTIME ) THEN
            EXPGEOFLAG = ENVYN( 'USE_EXP_GEO_CODES',
     &                          'Use expanded geographic codes', .FALSE., IOS )
            IF ( IOS .GT. 0 ) THEN
                CALL M3EXIT( PNAME,0,0, 'Bad env vble "USE_EXP_GEO_CODES"', 2 )
            END IF

            FIRSTTIME = .FALSE.
        END IF

        USEEXPGEO = EXPGEOFLAG

        RETURN

        END FUNCTION USEEXPGEO
