
INTEGER FUNCTION CVTRDTYPE( SMKROAD, RLASAFLAG, ULASAFLAG )

    !***********************************************************************
    !  function body starts at line 68
    !
    !  DESCRIPTION:
    !       Converts inventory road types into MOBILE6 type
    !
    !  PRECONDITIONS REQUIRED: none
    !
    !  SUBROUTINES AND FUNCTIONS CALLED: none
    !
    !  REVISION  HISTORY:
    !       10/01: Created by C. Seppanen
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
    !***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !....  Includes
    INCLUDE 'M6CNST3.h90'      ! MOBILE6 constants

    !....  Function arguments
    INTEGER, INTENT (IN) :: SMKROAD       ! road type in SMOKE code ( 1 - 19 )
    LOGICAL, INTENT (IN) :: RLASAFLAG     ! true: treat rural local roads as arterial
    LOGICAL, INTENT (IN) :: ULASAFLAG     ! true: treat urbal local roads as arterial

    !....  Road type parameters
    INTEGER, PARAMETER :: RURALINTERSTATE = 1       ! rural interstate
    INTEGER, PARAMETER :: RURALPRINCART   = 2       ! rural principle arterial
    INTEGER, PARAMETER :: RURALMINORART   = 6       ! rural minor arterial
    INTEGER, PARAMETER :: RURALMAJORCOLL  = 7       ! rural major collector
    INTEGER, PARAMETER :: RURALMINORCOLL  = 8       ! rural minor collector
    INTEGER, PARAMETER :: RURALLOCAL      = 9       ! rural local
    INTEGER, PARAMETER :: URBANINTERSTATE = 11      ! urban interstate
    INTEGER, PARAMETER :: URBANFREEWAY    = 12      ! urban freeway
    INTEGER, PARAMETER :: URBANPRINCART   = 14      ! urban principle arterial
    INTEGER, PARAMETER :: URBANMINORART   = 16      ! urban minor arterial
    INTEGER, PARAMETER :: URBANCOLL       = 17      ! urban collector
    INTEGER, PARAMETER :: URBANLOCAL      = 19      ! urban local

    !....  LOCAL VARIABLES and their descriptions:
    LOGICAL :: EFLAG      = .FALSE.       ! true: error found

    CHARACTER(300)          MESG          !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'CVTRDTYPE'     ! program name

    !***********************************************************************
    !   begin body of function CVTRDTYPE

    SELECT CASE( SMKROAD )
      CASE( RURALINTERSTATE )
        CVTRDTYPE = M6FREEWAY
      CASE( RURALPRINCART:RURALMINORCOLL )        ! rural arterials through collectors
        CVTRDTYPE = M6ARTERIAL
      CASE( RURALLOCAL )
        IF( RLASAFLAG ) THEN
            CVTRDTYPE = M6ARTERIAL
        ELSE
            CVTRDTYPE = M6LOCAL
        END IF
      CASE( URBANINTERSTATE:URBANFREEWAY )        ! urban interstate and freeway
        CVTRDTYPE = M6FREEWAY
      CASE( URBANPRINCART:URBANCOLL )             ! urban arterials through collectors
        CVTRDTYPE = M6ARTERIAL
      CASE( URBANLOCAL )
        IF( ULASAFLAG ) THEN
            CVTRDTYPE = M6ARTERIAL
        ELSE
            CVTRDTYPE = M6LOCAL
        END IF
      CASE DEFAULT
        EFLAG = .TRUE.

        WRITE( MESG, '(A,I8,A)' ) 'ERROR: Road type ', SMKROAD,&
                             ' is not recognized.'
        CALL M3MESG( MESG )
    END SELECT

    IF( EFLAG ) THEN
        MESG = 'Problem converting road type'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

END FUNCTION CVTRDTYPE
