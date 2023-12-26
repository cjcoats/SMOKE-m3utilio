
LOGICAL FUNCTION USEEXPGEO()

    !***********************************************************************
    !  function USEEXPGEO body starts at line 50
    !
    !  DESCRIPTION:
    !      The first time this function is called it will get the value
    !      of the environment variable USE_EXP_GEO_CODES to indicate if
    !      expanded geographic codes should be used. Subsequent calls
    !      will return the saved setting.
    !
    !  REVISION  HISTORY:
    !       Created 12/13 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" format, and related changes
    !
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !                System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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
    !****************************************************************************

    USE M3UTILIO

    !...........   Other local variables
    INTEGER        IOS    ! i/o status

    LOGICAL, SAVE :: FIRSTTIME  = .TRUE.
    LOGICAL, SAVE :: EXPGEOFLAG = .FALSE.

    CHARACTER(16), PARAMETER :: PNAME = 'USEEXPGEO' ! program name

    !****************************************************************************
    !   begin body of function USEEXPGEO

    IF ( FIRSTTIME ) THEN
        EXPGEOFLAG = ENVYN( 'USE_EXP_GEO_CODES',            &
                            'Use expanded geographic codes', .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "USE_EXP_GEO_CODES"', 2 )
        END IF

        FIRSTTIME = .FALSE.
    END IF

    USEEXPGEO = EXPGEOFLAG

    RETURN

END FUNCTION USEEXPGEO
