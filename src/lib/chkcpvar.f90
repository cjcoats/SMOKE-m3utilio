
INTEGER FUNCTION CHKCPVAR( CPVNAM, NIPOL, EINAM )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      Uses the rules of naming control and projection variables to compare
    !      a control or projection variable name to the list of inventory pollutant
    !      names.  Returns an integer value based on the comparison:
    !         Returns -3 if variable should be skipped
    !         Returns -2 if error
    !         Returns -1 if CPVNAM does not apply to any pollutants in EINAM
    !         Returns  0 if CPVNAM applies to all pollutants
    !         Returns  number of inventory pollutant (1 to NIPOL) if CPVNAM applies
    !            to a single inventory pollutant
    !      Rules of control/projection variable nameing:
    !         1) A variable name of "all" applies to ALL pollutants in the
    !            inventory. How it is applied depends on the type of control matrix.
    !         2) A variable name that matches an inventory pollutant name applies
    !            only to that inventory pollutant
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !**************************************************************************
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
    !***************************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......   INCLUDE FILES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: CPVNAM             ! Name of control/proj variable
    INTEGER     , INTENT (IN) :: NIPOL              ! Number of inventory pollutant names
    CHARACTER(*), INTENT (IN) :: EINAM( NIPOL )     ! List of inventory pollutant names

    !.......   Other local variables

    INTEGER         I, L1, L2

    LOGICAL         EFLAG

    CHARACTER(NAMLEN3)   CBUF
    CHARACTER(256)       MESG

    !***********************************************************************
    !   begin body of function CHKCPVAR

    IF( CPVNAM .EQ. 'all' ) THEN
        CHKCPVAR = 0
        RETURN
    ELSEIF( CPVNAM .EQ. 'scc' ) THEN
        CHKCPVAR = -3
        RETURN
    ENDIF

    CHKCPVAR = -1     ! Initialize

    EFLAG = .FALSE.
    DO I = 1, NIPOL
        IF( CHKCPVAR .EQ. -1 .AND.      &
            CPVNAM   .EQ. EINAM( I ) ) THEN
            CHKCPVAR = I
            CBUF = EINAM( I )

        ELSEIF( CHKCPVAR .NE. -1 .AND.  &
                CPVNAM   .EQ. EINAM( I ) ) THEN
            EFLAG = .TRUE.

            MESG = 'ERROR: Control/projection variable "' // TRIM( CPVNAM ) //  &
                   'applies to inventory pollutants "' // TRIM( CBUF ) //       &
                  '" and "' // TRIM( EINAM( I ) ) // '"'
            CALL M3MSG2( MESG )

        ENDIF

    ENDDO

    !.....  Handle error condition
    IF( EFLAG ) THEN
        CHKCPVAR = -2

        MESG = 'ERROR: Control/projection variables may only' //    &
               'be applied to one or all inventory pollutants'
        CALL M3MSG2( MESG )

    ENDIF

    RETURN

END FUNCTION CHKCPVAR

