
        INTEGER FUNCTION CHKCPVAR( CPVNAM, NIPOL, EINAM )

C***********************************************************************
C  function body starts at line
C
C  DESCRIPTION:
C      Uses the rules of naming control and projection variables to compare
C      a control or projection variable name to the list of inventory pollutant
C      names.  Returns an integer value based on the comparison:
C         Returns -3 if variable should be skipped
C         Returns -2 if error
C         Returns -1 if CPVNAM does not apply to any pollutants in EINAM
C         Returns  0 if CPVNAM applies to all pollutants
C         Returns  number of inventory pollutant (1 to NIPOL) if CPVNAM applies
C            to a single inventory pollutant
C      Rules of control/projection variable nameing:
C         1) A variable name of "all" applies to ALL pollutants in the
C            inventory. How it is applied depends on the type of control matrix.
C         2) A variable name that matches an inventory pollutant name applies
C            only to that inventory pollutant
C
C  PRECONDITIONS REQUIRED:
C
C  SUBROUTINES AND FUNCTIONS CALLED:
C
C  REVISION  HISTORY:
C       Version ??/???? by ???
C       Version 11/2023 by CJC:  bug-fix, INTENT, USE M3UTILIO and related changes
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

        USE M3UTILIO

        IMPLICIT NONE

C...........   INCLUDE FILES:
        INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

C...........   SUBROUTINE ARGUMENTS
        CHARACTER(*), INTENT (IN) :: CPVNAM         ! Name of control/proj variable
        INTEGER     , INTENT (IN) :: NIPOL          ! Number of inventory pollutant names
        CHARACTER(*), INTENT (IN) :: EINAM( NIPOL ) ! List of inventory pollutant names

C...........   Other local variables

        INTEGER         I, L1, L2

        LOGICAL         EFLAG

        CHARACTER(NAMLEN3)   CBUF
        CHARACTER(256)       MESG

        CHARACTER(16) :: PROGNAME = 'CHKCPVAR' ! program name

C***********************************************************************
C   begin body of function CHKCPVAR

        IF( CPVNAM .EQ. 'all' ) THEN
            CHKCPVAR = 0
            RETURN
        ELSEIF( CPVNAM .EQ. 'scc' ) THEN
            CHKCPVAR = -3
            RETURN
        ENDIF

        CHKCPVAR = -1 ! Initialize

        EFLAG = .FALSE.
        DO I = 1, NIPOL
            IF( CHKCPVAR .EQ. -1 .AND.
     &          CPVNAM   .EQ. EINAM( I ) ) THEN
                CHKCPVAR = I
                CBUF = EINAM( I )

            ELSEIF( CHKCPVAR .NE. -1 .AND.
     &              CPVNAM   .EQ. EINAM( I ) ) THEN
                EFLAG = .TRUE.

                MESG = 'ERROR: Control/projection variable "' //
     &                 TRIM( CPVNAM ) //
     &                 'applies to inventory pollutants "' //
     &                 TRIM( CBUF ) // '" and "' //
     &                 TRIM( EINAM( I ) ) // '"'
                CALL M3MSG2( MESG )

            ENDIF

        ENDDO

C.........  Handle error condition
        IF( EFLAG ) THEN
            CHKCPVAR = -2

            MESG = 'ERROR: Control/projection variables may only' //
     &             'be applied to one or all inventory pollutants!'
            CALL M3MSG2( MESG )

        ENDIF

        RETURN

        END FUNCTION CHKCPVAR

