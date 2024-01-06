
SUBROUTINE RDINVPOL( FILNAM, NCNT, VCNT, VNAMES, VTYPES, SRCID, POLDAT, EFLAG )

    !***********************************************************************
    !  subroutine body starts at line 82
    !
    !  DESCRIPTION:
    !      Reads inventory pollutant-specific data for variables listed in VNAMES
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
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

    !.......  MODULES for public variables

    !.......  This module is required by the FileSetAPI
    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDE FILES:
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: FILNAM               ! Logical file name
    INTEGER     , INTENT (IN) :: NCNT                 ! Number of records
    INTEGER     , INTENT (IN) :: VCNT                 ! No. vars other than SRCID
    CHARACTER(*), INTENT (IN) :: VNAMES( VCNT )       ! Variable names
    INTEGER     , INTENT (IN) :: VTYPES( VCNT )       ! Variable types
    INTEGER     , INTENT(OUT) :: SRCID ( NCNT )       ! Data
    REAL        , INTENT(OUT) :: POLDAT( NCNT,VCNT )     ! Data
    LOGICAL     , INTENT(OUT) :: EFLAG                ! true: error found

    !.......   Other local variables

    INTEGER         IREAD ( NCNT )      ! integer read array
    INTEGER         K, LV, V      ! counters and indices
    INTEGER         IOS           ! i/o status

    CHARACTER(NAMLEN3)   VARBUF
    CHARACTER(256)       MESG

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDINVPOL'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDINVPOL

    EFLAG = .FALSE.

    IF( .NOT. READSET( FILNAM, 'SRCID', ALLAYS3, 1, 0, 0, SRCID ) ) THEN
        EFLAG = .TRUE.
        MESG = 'Error reading "SRCID" from file: ' // CRLF()// BLANK10// TRIM( FILNAM )
        CALL M3MSG2( MESG )
    END IF

    !.......  Read variables by type, and store as REAL in POLDAT
    DO V = 1, VCNT

        VARBUF = VNAMES( V )

        IF( VTYPES( V ) .EQ. M3INT ) THEN

            IF( .NOT. READSET( FILNAM, VARBUF, ALLAYS3, ALLFILES, 0, 0, IREAD )) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Could not read "' // TRIM( VARBUF ) // '" from file.'
                CALL M3MSG2( MESG )

            ELSE

                POLDAT( :,V ) = REAL( IREAD )       ! Array

            END IF

        ELSE IF( .NOT. READSET( FILNAM, VARBUF, ALLAYS3, ALLFILES, 0, 0, POLDAT(1,V) ) ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Could not read "' // TRIM( VARBUF ) // '" from file.'
            CALL M3MSG2( MESG )

        END IF

    END DO

    RETURN

END SUBROUTINE RDINVPOL

