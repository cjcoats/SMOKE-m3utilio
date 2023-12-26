
SUBROUTINE GETBASYR( NSRC, BASEYEAR )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     This subroutine examines the years for all of the sources, and determines
    !     the base year from the most frequent of these base years.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:  M3EXIT
    !
    !  REVISION  HISTORY:
    !       Created 4/99 by M Houyoux
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
    !****************************************************************************

    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: INVYR

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   ARGUMENTS and their descriptions:
    INTEGER, INTENT (IN) :: NSRC            ! number of sources
    INTEGER, INTENT(OUT) :: BASEYEAR        ! base inventory year

    !.......   Argument dimensioned source arrrays
    INTEGER     INDX( NSRC )       ! sorting index for processing years
    INTEGER    YRGRP( NSRC )       ! group number of year

    !.......   LOCAL VARIABLES their descriptions:
    INTEGER         J, S           ! indices
    INTEGER         ICNT           ! counter
    INTEGER         IOS            ! i/o status
    INTEGER         GRPMAX         ! count of most prevalent year
    INTEGER         NYEARS         ! number of different years in inventory
    INTEGER         PGRP           ! year group from previous loop iteration
    INTEGER         PYEAR          ! year from previous loop iteration

    CHARACTER(300)  MESG           ! Message buffer

    !***********************************************************************
    !   begin body of subroutine GETBASYR

    !.....  Create sorting index for years
    DO S = 1, NSRC
        INDX( S ) = S
    END DO

    !.....  Sort years (from MODSOURC)
    CALL SORTI1( NSRC, INDX, INVYR )

    !.....  Count the number of different years in inventory
    PYEAR = 0
    ICNT  = 0
    DO S = 1, NSRC

        J = INDX( S )

        IF( INVYR( J ) .NE. PYEAR ) THEN

            ICNT = ICNT + 1

            PYEAR      = INVYR( J )

        END IF

        YRGRP( S ) = ICNT

    END DO

    !.....  Post-process list of year groups to find most prevalent year
    PGRP   = 1
    ICNT   = 0
    GRPMAX = 0
    NYEARS = 0
    IF ( NSRC .EQ. 1) THEN

        BASEYEAR = PYEAR
        NYEARS = 1
        PGRP   = YRGRP( 1 )
        PYEAR = INVYR( INDX( 1 ) )

    ELSE
        DO S = 1, NSRC

            IF( YRGRP( S ) .NE. PGRP .OR. S .EQ. NSRC ) THEN

                IF( ICNT .GT. GRPMAX ) THEN
                    GRPMAX = ICNT
                    BASEYEAR = PYEAR
                END IF

                NYEARS = NYEARS + 1
                ICNT   = 0
                PGRP   = YRGRP( S )

            END IF

            PYEAR = INVYR( INDX( S ) )
            ICNT  = ICNT + 1

        END DO
    ENDIF

    IF( NYEARS .EQ. 1 ) THEN

        WRITE( MESG,94010 ) 'NOTE: Inventory base year set to ', BASEYEAR

    ELSE

        WRITE( MESG,94010 ) 'NOTE: The number of inventory years was ', NYEARS, &
               'and the inventory base' // CRLF() // BLANK10 //                 &
               'year was set to the most frequent year, ', BASEYEAR

    END IF

    CALL M3MSG2( MESG )
    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.....94xxx

94010 FORMAT( 10( A, :, I6, :, 2X ) )

END SUBROUTINE GETBASYR
