
SUBROUTINE INITSTCY

    !***********************************************************************
    !  subroutine INITSTCY body starts at line
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to initialize the necessary fields
    !      for performing state and county totals.  The first call sets up the
    !      indices from each source to each county.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 8/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
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

    !.......  MODULES for public variables
    !.......  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: NMSRC, MCFIP,     &    ! no. of sources by category
                        LREPSTA,          &    ! output state total emissions flag
                        LREPCNY,          &    ! output county total emissions flag
                        LREPSCC,          &    ! output SCC total emissions flag
                        MEBSTA, &    ! state total speciated emissions
                        MEBCNY, &    ! county total speciated emissions
                        MEBSUM, &    ! source total speciated emissions total
                        MEBSRC, &    ! source total speciated emissions by hour
                        MEBSCC, &    ! SCC total speciated emissions
                        MEBSTC      ! state-SCC total speciated emissions

    !.......  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: MICNY, NCOUNTY, CNTYCOD

    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: CIFIP, CSCC

    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: MISCC, DISFLTYP, DISFL, GASFLTYP, GASFL, ETHFLTYP, ETHFL,&
                         CNGFLTYP, CNGFL, LPGFLTYP, LPGFL
    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVSCC, INVSCC, NINVIFIP, INVCFIP

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   Other local variables

    INTEGER          IOS          ! i/o status
    INTEGER          J, K, S      ! counter
    INTEGER          IDISFL        ! diesel fuel type
    INTEGER          IGASFL        ! gasoline fuel type
    INTEGER          IETHFL        ! ethanol fuel type
    INTEGER          ICNGFL        ! CNG fuel type
    INTEGER          ILPGFL        ! LPG fuel type

    LOGICAL, SAVE :: FIRSTIME = .TRUE.     ! true: first time routine called

    CHARACTER(FIPLEN3) CFIP          ! tmp cy/st/co code
    CHARACTER(FIPLEN3) PFIP          ! previous cy/st/co code

    CHARACTER(300)   MESG         ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'INITSTCY'     ! program name

    !***********************************************************************
    !   begin body of subroutine INITSTCY

    IF( FIRSTIME ) THEN

    !.......  Allocate memory for indices from Co/st/cy codes to counties
    !.......  Allocate memory for index from master list of SCCs to source SCC
        ALLOCATE( MICNY( NMSRC ),&
                  MISCC( NMSRC ),&
                  MCFIP( NMSRC ),&
                  GASFL( NMSRC ),&
                  DISFL( NMSRC ),&
                  CNGFL( NMSRC ),&
                  LPGFL( NMSRC ),&
                  ETHFL( NMSRC ), STAT=IOS )           ! arry for ethnol fuel type or not
        CALL CHECKMEM( IOS, 'MICNY...ETHFL', PROGNAME )
        DISFL = .FALSE.
        GASFL = .FALSE.
        CNGFL = .FALSE.
        LPGFL = .FALSE.
        ETHFL = .FALSE.

        IGASFL = ENVINT( 'GASOLINE_FUEL_CODE', 'Gasoline fuel type code [ex: 1]', 1, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "GASOLINE_FUEL_CODE"', 2 )
        END IF
        WRITE( GASFLTYP,'( I2.2)' ) IGASFL

        IDISFL = ENVINT( 'DIESEL_FUEL_CODE', 'Diesel fuel type code [ex: 2]', 2, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "DIESEL_FUEL_CODE"', 2 )
        END IF
        WRITE( DISFLTYP,'( I2.2)' ) IDISFL

        ICNGFL = ENVINT( 'CNG_FUEL_CODE', 'CNG fuel type code [ex: 3]', 3, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "CNG_FUEL_CODE"', 2 )
        END IF
        WRITE( CNGFLTYP,'( I2.2)' ) ICNGFL

        ILPGFL = ENVINT( 'LPG_FUEL_CODE', 'LPG fuel type code [ex: 4]', 4, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "LPG_FUEL_CODE"', 2 )
        END IF
        WRITE( LPGFLTYP,'( I2.2)' ) ILPGFL

        IETHFL = ENVINT( 'ETHANOL_FUEL_CODE', 'Ethanol fuel type code [ex: 5]', 5, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ETHANOL_FUEL_CODE"', 2 )
        END IF
        WRITE( ETHFLTYP,'( I2.2)' ) IETHFL

    !.......  Create indices to counties from Co/st/cy codes and for SCCs
        PFIP = ' '
        DO S = 1, NMSRC

            CFIP = CIFIP( S )

            IF( CFIP .NE. PFIP ) THEN

                J = MAX( FINDC( CFIP, NCOUNTY, CNTYCOD ), 0 )
                K = MAX( FINDC( CFIP, NINVIFIP,INVCFIP ), 0 )
                PFIP = CFIP

            END IF

            MICNY( S ) = J
            MCFIP( S ) = K
            MISCC( S ) = MAX( FINDC( CSCC( S ), NINVSCC, INVSCC ), 0 )

    !......   Check whether diseel fuel type or not
            IF( CSCC(S)( 11:12 ) == '22' ) THEN
                IF( CSCC(S)( 13:14 ) == GASFLTYP ) GASFL( S ) = .TRUE.
                IF( CSCC(S)( 13:14 ) == DISFLTYP ) DISFL( S ) = .TRUE.
                IF( CSCC(S)( 13:14 ) == CNGFLTYP ) CNGFL( S ) = .TRUE.
                IF( CSCC(S)( 13:14 ) == LPGFLTYP ) LPGFL( S ) = .TRUE.
                IF( CSCC(S)( 13:14 ) == ETHFLTYP ) ETHFL( S ) = .TRUE.
            END IF

        END DO

        FIRSTIME = .FALSE.

    END IF

    !.......  Initialize totals to zero...
    !.......  SCC totals...
    IF( LREPSCC ) THEN
        MEBSCC = 0.
    END IF

    !.......  State totals...
    IF( LREPSTA ) THEN
        MEBSTA = 0.
        IF( LREPSCC ) THEN
            MEBSTC = 0.
        END IF
    END IF

    !.......  County totals...
    IF( LREPCNY ) THEN
        MEBCNY = 0.
    END IF

    !.......  Source totals...
    MEBSUM = 0.

    RETURN

END SUBROUTINE INITSTCY
