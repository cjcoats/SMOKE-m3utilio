
SUBROUTINE GETINAME( CATEGORY, ENAME, ANAME )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     This subroutine sets the inventory file names based on the source category
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 3/99 by M Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    IMPLICIT NONE

    !.......   ARGUMENTS and their descriptions:
    CHARACTER(*), INTENT (IN) :: CATEGORY      ! category name
    CHARACTER(*), INTENT(OUT) :: ENAME         ! i/o api inventory file name
    CHARACTER(*), INTENT(OUT) :: ANAME         ! ASCII inventory file name

    !.......   LOCAL VARIABLES their descriptions:

    CHARACTER(300)  MESG                       ! Message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETINAME'      ! Program name

    !***********************************************************************
    !   begin body of subroutine GETINAME

    SELECT CASE( CATEGORY )

      CASE( 'AREA' )

        ENAME = 'AREA'
        ANAME = 'ASRC'

      CASE( 'BIOGEN' )

        ENAME = 'BGRD'
        ANAME = '    '

      CASE( 'MOBILE' )

        ENAME = 'MOBL'
        ANAME = 'MSRC'

      CASE( 'POINT' )

        ENAME = 'PNTS'
        ANAME = 'PSRC'

      CASE DEFAULT
        MESG = 'Category ' // CATEGORY // ' not known in program.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END SELECT

    RETURN

END SUBROUTINE GETINAME
