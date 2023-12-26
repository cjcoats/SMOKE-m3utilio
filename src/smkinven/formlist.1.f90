
SUBROUTINE FORMLIST

!**************************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine is designed to read SMKINVEN_FORMULA and determine how many formulas and list of variables
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 5/2012 by B. Baek
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!**************************************************************************
    USE M3UTILIO

!...........   Modules for public variables
!...........   This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVTBL, ITNAMA, ITCASA

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NMAP,NIACT, NIPPA, NPPOL, EANAM, ACTVTY,&
    &                   NCOMP, VAR_FORMULA, CHKPLUS, CHKMINUS,&
    &                   FORMULAS, VIN_A, VIN_B, VNAME
    IMPLICIT NONE

!.........   Subroutine arguments

!.........   Local variables
    INTEGER         I, J, F, L, N, V1, V2, VA, VB    ! indices and counters

    INTEGER         IOS                   ! i/o status
    INTEGER         LEQU                  ! position of '=' in formula
    INTEGER         LDIV                  ! position of '-' or '+' in formula
    INTEGER         LMNS                  ! position of '-' in formula
    INTEGER         LPLS                  ! position of '+' in formula

    LOGICAL,SAVE :: FIRSTIME = .TRUE.     ! true: first time
    LOGICAL      :: TFLAG    = .FALSE.    ! true: current formula has error
    LOGICAL      :: EFLAG    = .FALSE.    ! true: error found

    CHARACTER(300)  MESG                  ! Message buffer

    CHARACTER(16) :: PROGNAME = 'FORMLIST'    !  program name

!***********************************************************************
!   Begin body of subroutine FORMULAS

!.........  Figure out how many variables there are based on the
!           number of commas found in the string.
    NCOMP = 1
    L = LEN_TRIM( VAR_FORMULA )
    DO I = 1, L
        IF( VAR_FORMULA( I:I ) == ',' ) NCOMP = NCOMP + 1
    ENDDO

    NMAP = NMAP + NCOMP  ! NCOMP more variable(s) to map

!.........  Allocate array to store formulas
    IF( FIRSTIME ) THEN
        ALLOCATE( CHKPLUS( NCOMP ),&
        &         CHKMINUS( NCOMP ),&
        &         FORMULAS( NCOMP ),&
        &            VIN_A( NCOMP ),&
        &            VIN_B( NCOMP ),&
        &            VNAME( NCOMP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CHKPLUS...VNAME', PROGNAME )
        FIRSTIME = .FALSE.
    END IF
!.........  Split out formulas in string to array
    CALL PARSLINE( VAR_FORMULA, NCOMP, FORMULAS )

!.........  Loop through formulas
    DO F = 1, NCOMP

!.............  Make sure formula makes sense
        LEQU = INDEX( FORMULAS( F ), '=' )
        LPLS = INDEX( FORMULAS( F ), '+' )
        LMNS = INDEX( FORMULAS( F ), '-' )

        CHKPLUS( F )  = ( LPLS .GT. 0 )
        CHKMINUS( F ) = ( LMNS .GT. 0 )

        LDIV = LPLS
        IF( CHKMINUS( F ) ) LDIV = LMNS

        IF( LEQU .LE. 0 .OR.&
        &  ( .NOT. CHKPLUS(F) .AND. .NOT. CHKMINUS(F) ) ) THEN
            MESG = 'Could not interpret formula for extra ' //&
            &       'pollutant from ' // 'SMKINVEN_FORMULAS=' //&
            &       TRIM( FORMULAS( F ) )
            CALL M3MSG2( MESG )
            TFLAG = .TRUE.
        END IF

!.............  Extract formula variable names
        L      = LEN_TRIM( FORMULAS( F ) )
        VNAME( F )= ADJUSTL( FORMULAS( F )(      1:LEQU-1 ) )
        VIN_A( F )= ADJUSTL( FORMULAS( F )( LEQU+1:LDIV-1 ) )
        VIN_B( F )= ADJUSTL( FORMULAS( F )( LDIV+1:L      ) )

!.............  Find formula inputs in existing variable list
        J = INDEX1( VIN_A( F ), NINVTBL, ITCASA )
        IF( J < 1 ) THEN
            L = LEN_TRIM( VIN_A( F ) )
            MESG = 'Variable "'// VIN_A( F )( 1:L ) //&
            &   '" from formula was not found in inventory ' //&
            &   'pollutant code (CAS nubmer)'
            CALL M3MSG2( MESG )

        ELSE
            VIN_A( F ) = ITNAMA( J )

        END IF

        J = INDEX1( VIN_B( F ), NINVTBL, ITCASA )
        IF( J < 1 ) THEN
            L = LEN_TRIM( VIN_B( F ) )
            MESG = 'Variable "'// VIN_B( F )( 1:L ) //&
            &   '" from formula was not found in inventory ' //&
            &   'pollutant code (CAS nubmer)'
            CALL M3MSG2( MESG )

        ELSE
            VIN_B( F ) = ITNAMA( J )

        END IF

        VA = INDEX1( VIN_A( F ), NIPPA, EANAM )
        VB = INDEX1( VIN_B( F ), NIPPA, EANAM )

        IF( VA .LE. 0 ) THEN
            TFLAG = .TRUE.
            L = LEN_TRIM( VIN_A( F ) )
            MESG = 'Variable "'// VIN_A( F )( 1:L ) //&
            &       '" from formula was not found in inventory.'
            CALL M3MSG2( MESG )
        END IF

        IF( VB .LE. 0 ) THEN
            TFLAG = .TRUE.
            L = LEN_TRIM( VIN_B( F ) )
            MESG = 'Variable "'// VIN_B( F )( 1:L ) //&
            &       '" from formula was not found in inventory.'
            CALL M3MSG2( MESG )
        END IF

        V1 = INDEX1( VIN_A( F ), NIACT, ACTVTY )
        V2 = INDEX1( VIN_B( F ), NIACT, ACTVTY )

        IF( V1 .GT. 0 ) THEN
            TFLAG = .TRUE.
            L = LEN_TRIM( VIN_A( F ) )
            MESG = 'ERROR: Variable "'//VIN_A(F)(1:L)//'" is an'//&
            &       'activity, which is not allowed in a formula.'
            CALL M3MSG2( MESG )
        END IF

        IF( V2 .GT. 0 ) THEN
            TFLAG = .TRUE.
            L = LEN_TRIM( VIN_B( F ) )
            MESG = 'ERROR: Variable "'//VIN_B(F)(1:L)//'" is an'//&
            &       'activity, which is not allowed in a formula.'
            CALL M3MSG2( MESG )
        END IF

        IF( TFLAG ) THEN
            WRITE( MESG,94010 ) 'ERROR: Problem processing '//&
            &       'formula', F, ': "'//TRIM(FORMULAS(F)) // '"'
            CALL M3MSG2( MESG )
        END IF

        IF ( TFLAG ) EFLAG = .TRUE.

    END DO

    IF ( EFLAG ) THEN
        MESG = 'ERROR: Problem processing formulas. ' //&
        &       'See previous error messages.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!******************  FORMAT  STATEMENTS   ******************************

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE FORMLIST
