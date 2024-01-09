 
SUBROUTINE SETREFCNTY

    !***********************************************************************
    !  subroutine SETREFCNTY body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine sets up reference county information for Movesmrg.
    !      It reads the county cross-reference file, fuel month file, and
    !      reference county emission factors list.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 4/10 by C. Seppanen
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

    !.......  MODULES for public variables
    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: XDEV, MDEV, FDEV, NSRCCELLS, NREFSRCS, REFSRCS

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC

    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: CIFIP

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, INVCFIP

    !.......  This module is used for reference county information
    USE MODMBSET, ONLY: NINVC, NREFC, MCREFSORT, MCREFIDX

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  LOCAL VARIABLES and their descriptions:

    !.......  Array to store counties inside grid
    CHARACTER(FIPLEN3)  GRDFIP( NSRC )

    !.......  Arrays for building list of reference county sources
    INTEGER           SRCREFIDX( NSRC )        ! ref. county index for each source
    CHARACTER(FIPLEN3)  INVCNTY( NINVC )        ! inventory counties sorted by reference county
    CHARACTER(FIPLEN3)  REFCNTY( NREFC )        ! list of sorted reference counties

    !.......  Other local variables
    INTEGER   I, J, K, S      ! indexes and counters
    INTEGER   MXNREFSRCS      ! max. no. sources per reference county
    INTEGER   NGRDFIP         ! no. counties inside grid
    INTEGER   REFIDX          ! current reference county index
    INTEGER   IOS             ! error status

    LOGICAL :: SKIPFIP = .FALSE.       ! true: county is not in cross-reference

    CHARACTER(FIPLEN3) PFIP        ! prev fips
    CHARACTER(FIPLEN3) REFFIP      ! prev fips
    CHARACTER(300)     MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'SETREFCNTY'     ! program name

    !***********************************************************************
    !   begin body of subroutine SETREFCNTY

    !.......  Build list of counties within the grid
    GRDFIP  = ' '       ! array
    PFIP    = ' '
    NGRDFIP = 0
    DO S = 1, NSRC

        IF( NSRCCELLS( S ) == 0 ) CYCLE

        IF( CIFIP( S ) .NE. PFIP ) THEN
            NGRDFIP = NGRDFIP + 1
            GRDFIP( NGRDFIP ) = CIFIP( S )
            PFIP = CIFIP( S )
        END IF

    END DO

    !.......  Read county cross-reference file
    CALL RDMXREF( XDEV, NGRDFIP, GRDFIP )

    !.......  Build list of inventory counties sorted by reference county
    ALLOCATE( NREFSRCS( NREFC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'NREFSRCS', PROGNAME )

    DO I = 1, NINVC
        INVCNTY( I ) = MCREFSORT( I,1 )
    END DO


    DO I = 1, NREFC
        REFCNTY( I ) = MCREFIDX( I,1 )
    END DO

    !.......  Build list of sources for each reference county
    !           Start by counting number of sources for each reference county
    NREFSRCS = 0          ! array
    SRCREFIDX = 0       ! array

    PFIP = ' '
    REFFIP = ' '
    REFIDX = 0
    DO S = 1, NSRC

        IF( NSRCCELLS( S ) == 0 ) CYCLE

        IF( CIFIP( S ) .NE. PFIP ) THEN
            SKIPFIP = .FALSE.
            PFIP = CIFIP( S )

            J = INDEX1( CIFIP( S ), NINVC, INVCNTY )
            IF( J <= 0 ) THEN
                SKIPFIP = .TRUE.
                WRITE( MESG, 94010 ) 'WARNING: No emissions will '//&
                  'be calculated for inventory county ' //CIFIP( S )//&
                  ' because it is not listed in the county '//&
                  'cross-reference file'
                CALL M3MESG( MESG )
                CYCLE
            END IF

            IF( MCREFSORT( J,2 ) .NE. REFFIP ) THEN
                REFFIP = MCREFSORT( J,2 )

                REFIDX = INDEX1( REFFIP, NREFC, REFCNTY )
                IF( REFIDX <= 0 ) THEN
                    WRITE( MESG, 94010 ) 'INTERNAL ERROR: ' //&
                      'Problem with reference county mapping for '&
                      // 'county ' // CIFIP( S )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF
        END IF

        IF( .NOT. SKIPFIP ) THEN

            NREFSRCS( REFIDX ) = NREFSRCS( REFIDX ) + 1

            SRCREFIDX( S ) = REFIDX

        END IF

    END DO

    !.......  Get maximum number of sources per reference county
    MXNREFSRCS = 0
    DO I = 1, NREFC

        IF( NREFSRCS( I ) > MXNREFSRCS ) MXNREFSRCS = NREFSRCS( I )

    END DO

    ALLOCATE( REFSRCS( NREFC, MXNREFSRCS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'REFSRCS', PROGNAME )
    REFSRCS = 0       ! array
    NREFSRCS = 0      ! array

    REFIDX = 0
    DO S = 1, NSRC

        IF( NSRCCELLS( S ) == 0 ) CYCLE

        REFIDX = SRCREFIDX( S )
        IF( REFIDX == 0 ) CYCLE

        K = NREFSRCS( REFIDX ) + 1
        REFSRCS( REFIDX, K ) = S
        NREFSRCS( REFIDX ) = K

    END DO

    !.......  Read fuel month reference file
    CALL RDFMREF( MDEV )

    !.......  Read emission factors file list
    CALL RDMRCLIST( FDEV )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.......94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE SETREFCNTY
