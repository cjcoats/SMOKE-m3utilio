
SUBROUTINE MRGUNITS

    !***********************************************************************
    !  subroutine MRGUNITS body starts at line 93
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/99 by M. Houyoux
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

    !.......  MODULES for public variables
    !.......  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: SMATCHK, GRDFAC, TOTFAC,    &
                        GRDUNIT, TOTUNIT, NMSPC, NIPPA, NUNITS

    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: SPCUNIT_L, SPCUNIT_S, GRDENV, TOTENV

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(16), EXTERNAL :: MULTUNIT
    REAL         , EXTERNAL :: UNITFAC

    !.......   Other local variables

    INTEGER         IOS          ! tmp I/O status
    INTEGER         L, V            ! counter

    REAL            EMFAC, FAC1, FAC2

    CHARACTER(300)  BUFFER       ! text buffer
    CHARACTER(300)  MESG         ! message buffer

    CHARACTER(NAMLEN3) EMUNIT          ! converted emissions unit

    CHARACTER(NAMLEN3) GRDUNIT_I       ! initialized gridded outputs units
    CHARACTER(NAMLEN3) GDEN_I          ! initialized gridded denominator
    CHARACTER(NAMLEN3) GNUM_I          ! initialized gridded numerator
    CHARACTER(NAMLEN3) GDEN            ! work gridded denominator
    CHARACTER(NAMLEN3) GNUM            ! work gridded numerator
    CHARACTER(NAMLEN3) GRDBUF          ! work gridded output units

    CHARACTER(NAMLEN3) TOTUNIT_I       ! initialized output totals units
    CHARACTER(NAMLEN3) TDEN_I          ! initialized totals denominator
    CHARACTER(NAMLEN3) TNUM_I          ! initialized totals numerator
    CHARACTER(NAMLEN3) TDEN            ! work totals denominator
    CHARACTER(NAMLEN3) TNUM            ! work totals numerator
    CHARACTER(NAMLEN3) TOTBUF          ! work output totals  units

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGUNITS'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGUNITS

    !.......  Allocate memory for units conversion factors and units.

    !.......  If using speciation, allocate arrays to number of species;
    !           otherwise, use total number of activities and pollutants.
    !           This approach assumes that each species has the same units
    !           regardless of which pollutant it comes from, i.e.
    !           ALD2 from VOC must have the same units as ALD2 from ACETALD
    NUNITS = NMSPC

    ALLOCATE( GRDFAC( NUNITS ),         &
             GRDUNIT( NUNITS ),         &
              TOTFAC( NUNITS+NIPPA ),   &
             TOTUNIT( NUNITS+NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GRDFAC...TOTUNIT', PROGNAME )

    !.......  Initialize all
    GRDFAC  = 1.       ! array
    TOTFAC  = 1.       ! array
    GRDUNIT = ' '      ! array
    TOTUNIT = ' '      ! array

    !.......  Loop through pollutants, and create output units accordingly
    !.......  For speciation, use the first unit for the speciation units from
    !           a given pollutants speciation factors
    DO V = 1, NUNITS

        !.......  Initialize the output units
        CALL UNITMATCH( SPCUNIT_L( V ) )
        CALL UNITMATCH( SPCUNIT_S( V ) )

        !.......  Convert emissions units (g/hr) to tons/hr
        EMUNIT = 'tons/hr'
        GRDUNIT_I = MULTUNIT( SPCUNIT_L( V ), EMUNIT )
        TOTUNIT_I = MULTUNIT( MULTUNIT( SPCUNIT_S( V ), EMUNIT ), 'hr/day' )

        !.......  Set the trial units
        GRDBUF = GRDUNIT_I
        TOTBUF = TOTUNIT_I
        IF( GRDENV .NE. ' ' ) GRDBUF = GRDENV
        IF( TOTENV .NE. ' ' ) TOTBUF = TOTENV

        !.......  Set the numerators and denominators
        L = INDEX( GRDUNIT_I, '/' )
        GNUM_I = ADJUSTL( GRDUNIT_I(   1:L-1     ) )
        GDEN_I = ADJUSTL( GRDUNIT_I( L+1:NAMLEN3 ) )

        L = INDEX( TOTUNIT_I, '/' )
        TNUM_I = ADJUSTL( TOTUNIT_I(   1:L-1     ) )
        TDEN_I = ADJUSTL( TOTUNIT_I( L+1:NAMLEN3 ) )

        L = INDEX( GRDBUF, '/' )
        GNUM = ADJUSTL( GRDBUF(   1:L-1     ) )
        GDEN = ADJUSTL( GRDBUF( L+1:NAMLEN3 ) )

        L = INDEX( TOTBUF, '/' )
        TNUM = ADJUSTL( TOTBUF(   1:L-1     ) )
        TDEN = ADJUSTL( TOTBUF( L+1:NAMLEN3 ) )

        !.......  Get factor for the numerators for the gridded outputs...
        FAC1 = UNITFAC( SPCUNIT_L( V ), GRDBUF, .TRUE. )          ! speciation

        !.......  Get factor for the denominators for the gridded outputs
        FAC2 = UNITFAC( EMUNIT, GRDBUF, .FALSE. )

        !.......  In case E.V. setting was bogus rebuild output units based
        !               on which factors were valid
        !.......  Also set negative factors (from unknown conversions) to 1.
        CALL CORRECT_UNITS( GNUM_I, GDEN_I, GNUM, GDEN, GRDBUF )

        !.......  When SMAT is used to calcuate model species
        IF( SMATCHK ) THEN
            EMFAC = UNITFAC( 'g/hr', EMUNIT, .TRUE. )
            FAC1 = FAC1 * EMFAC
        END IF

        !.......  Set factors for gridded outputs
        GRDFAC( V ) = FAC1 / FAC2

        !.......  Get conversion factor for the numerators for totals
        FAC1 = UNITFAC( SPCUNIT_S( V ), TOTBUF, .TRUE. )          ! speciation

        !.......  Get factors for the denominators for the totals.  Note that
        !               the hourly data are output as daily totals.
        FAC2 = UNITFAC( TOTUNIT_I, TOTBUF, .FALSE. )

        !.......  In case E.V. setting was bogus rebuild output units based
        !               on which factors were valid
        !.......  Also set negative factors (from unknown conversions) to 1.
        CALL CORRECT_UNITS( TNUM_I, TDEN_I, TNUM, TDEN, TOTBUF )

        !.......  Set factors for totaled outputs
        TOTFAC( V ) = FAC1 / FAC2

        !.......  Set the output units per pollutant/activity
        GRDUNIT( V ) = GRDBUF
        TOTUNIT( V ) = TOTBUF

    END DO

    !.......  Set up units for non-speciated emissions report (tons/day for now)
    EMFAC = UNITFAC( 'g/hr', 'tons/day', .TRUE. )
    DO V = 1, NIPPA
        TOTFAC ( NUNITS+V ) = EMFAC
        TOTUNIT( NUNITS+V ) = 'tons/day'
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.......94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram corrects badly formatted units
    SUBROUTINE CORRECT_UNITS( NUM_I, DEN_I, NUM, DEN, OUTUNIT )

        !.......  Subprogram arguments
        CHARACTER(NAMLEN3), INTENT(IN   ) :: NUM_I
        CHARACTER(NAMLEN3), INTENT(IN   ) :: DEN_I
        CHARACTER(NAMLEN3), INTENT(INOUT) :: NUM
        CHARACTER(NAMLEN3), INTENT(INOUT) :: DEN
        CHARACTER(NAMLEN3), INTENT(  OUT) :: OUTUNIT

        !----------------------------------------------------------------------

        IF( FAC1 .LT. 0. ) THEN
            NUM = NUM_I
            FAC1 = 1.
        END IF
        IF( FAC2 .LT. 0. ) THEN
            DEN  = DEN_I
            FAC2 = 1.
        END IF

        OUTUNIT =  TRIM( NUM ) // '/' // TRIM( DEN )

    END SUBROUTINE CORRECT_UNITS

END SUBROUTINE MRGUNITS
