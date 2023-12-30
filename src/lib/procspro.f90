
SUBROUTINE PROCSPRO( NMSPC, SPCNAM )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      Abridge full pollutants table to a data structure that can be used
    !      in ASGNSPRO
    !
    !  PRECONDITIONS REQUIRED:
    !      Expects cross-reference tables to be set to IMISS3 if not defined
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
    !****************************************************************************/
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

    !......   MODULES for public variables
    !......   This module contains the speciation profile tables
    USE MODSPRO, ONLY: NSPFUL, NSPROF, NSPECIES, SPROFN,    &
                       IDXSPRO, IDXSSPEC, INPRF, SPECID

    IMPLICIT NONE

    !......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !......  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NMSPC               ! No. of species for 1 pol
    CHARACTER(*), INTENT (IN) :: SPCNAM( NMSPC )     ! Species names for 1 pol

    !......  Other local variables
    INTEGER          I, J, K, L1          !  counters and indices
    INTEGER          IOS                  !  i/o status

    CHARACTER(300)     MESG           ! message buffer

    CHARACTER(NAMLEN3) SBUF           ! species name buffer
    CHARACTER(SPNLEN3) PCODE          ! current speciation profile code
    CHARACTER(SPNLEN3) PREVCODE       ! previous speciation profile code

    CHARACTER(16) :: PROGNAME = 'PROCSPRO'     ! program name

    !***********************************************************************
    !   begin body of subroutine PROCSPRO

    !......  Loop through sorted, unprocessed speciation profiles table to count
    !           the number of unique profiles

    PREVCODE = EMCMISS3
    J = 0
    DO I = 1, NSPFUL

        PCODE = INPRF( I )

        IF( PCODE .NE. PREVCODE ) THEN
            J = J + 1
            PREVCODE = PCODE
        ENDIF

    END DO            !  end loop on speciation profile table

    NSPROF = J

    !......  Allocate memory for unique profiles arrays
    ALLOCATE(  SPROFN( NSPROF ),        &
              IDXSPRO( NSPROF ),        &
             NSPECIES( NSPROF ),        &
             IDXSSPEC( NSPROF,NMSPC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IDXSSPEC', PROGNAME )

    !......  Initialize processed speciation profiles counter array
    NSPECIES = 0       ! array

    !......  Loop through sorted, unprocessed speciation profiles table for
    !           current pollutant to process it.  This algorithm works because we
    !           have stored (in INPRF) only those profiles for the current
    !           pollutant.

    PREVCODE = EMCMISS3
    J = 0
    DO I = 1, NSPFUL

        PCODE = INPRF ( I )
        SBUF  = SPECID( I )

        IF( PCODE .NE. PREVCODE ) THEN

            J = J + 1

            SPROFN  ( J ) = PCODE
            IDXSPRO ( J ) = I

            PREVCODE = PCODE

        END IF

        K = FINDC( SBUF, NMSPC, SPCNAM )

        IF( K .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: model species "' //             &
                   TRIM( SBUF ) // '" not found in sorted ' //      &
                   'names list developed for this pollutant'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ELSE

            NSPECIES( J ) = NSPECIES( J ) + 1
            IDXSSPEC( J,NSPECIES( J ) ) = K

        END IF

    END DO            !  end loop on speciation profile table

    RETURN

END SUBROUTINE PROCSPRO
