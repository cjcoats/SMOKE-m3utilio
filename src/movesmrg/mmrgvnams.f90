
SUBROUTINE MRGVNAMS

    !***********************************************************************
    !  subroutine MRGVNAMS body starts at line
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to merge the pollutant-to-species,
    !      pollutant, and species names from all of the open source categories,
    !      and populate the global unique lists of these names.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
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
    !*************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: PDEV, MNIPPA, NIPPA, MEANAM, EINAM, EMNAM,      &
                        EANAM, EMIDX, NSMATV, TSVDESC, NMSPC

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVSTAT, INVDCOD, INVDVTS

    !.......  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: NHAP, HAPNAM, MNSMATV_L, MSVDESC_L,            &
                         MSVUNIT_L, MSVUNIT_S, SPCUNIT_L, SPCUNIT_S,    &
                         GRDENV, TOTENV

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE

    !.......   Temporary arrays for building sorted pol-to-species names list
    INTEGER             INDXA  ( MNSMATV_L )     ! sorting index
    CHARACTER(PLSLEN3)  TVSORTA( MNSMATV_L )     ! basis for sorting
    CHARACTER(PLSLEN3)  TVDESCA( MNSMATV_L )     ! pol-to-spec concat

    !.......   Other local variables
    INTEGER         I, J, J1, J2, L, K, L1, L2, M, N, V      ! counters and indices
    INTEGER         LJ                  ! string length of emis types joiner
    INTEGER         IOS                 ! i/o error status
    INTEGER         NCNT                ! counter

    LOGICAL         EFLAG               ! error flag
    LOGICAL         FOUND               ! true: found HAP

    CHARACTER(NAMLEN3)   CPOL         ! tmp pol/act buffer
    CHARACTER(NAMLEN3)   BNUM         ! tmp biogenic units numerator
    CHARACTER(300)       MESG         ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGVNAMS'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGVNAMS
    EFLAG = .FALSE.
    FOUND = .FALSE.
    !.......  Read, sort, and store pollutant codes/names file
    CALL RDCODNAM( PDEV )

    !.......  Check if list of pollutants from MEPROC file contains
    !         any HAPs
    DO I = 1, MNIPPA
        J = INDEX1( MEANAM( I ), MXIDAT, INVDNAM )

        IF( J .GT. 0 ) THEN
            IF( INVDVTS( J ) == 'V' .OR.&
                INVDVTS( J ) == 'T' ) THEN
                FOUND = .TRUE.
                EXIT
            END IF
        END IF

    END DO

    !.......  If any HAP was found, update TOG to NONHAPTOG
    IF( FOUND ) THEN
        DO I = 1, MNIPPA
            IF( MEANAM( I ) == 'TOG' ) THEN
                MEANAM( I ) = 'NONHAP'//TRIM( MEANAM(I) )
            END IF
        END DO
    END IF

    !.......  Loop through emission process/pollutant combinations
    !         and update status of entry in master list.
    !         Also count number of pollutants to subtract when calculating
    !         NONHAPTOG
    NHAP = 0
    DO I = 1, MNIPPA

        CPOL = MEANAM( I )
        J = INDEX1( CPOL, MXIDAT, INVDNAM )
        IF( J .GT. 0 ) THEN
            INVSTAT( J ) = INVSTAT( J ) * 2

            IF( INVDVTS( J ) == 'V' .OR.    &
                INVDVTS( J ) == 'T' ) THEN
                NHAP = NHAP + 1
            END IF
        ELSE
            EFLAG = .TRUE.
            MESG = 'ERROR: Variable "' // TRIM( CPOL ) //   &
                   '" from MEPROC file is not in ' //       &
                   CRLF() // BLANK10 // 'inventory table.'
            CALL M3MSG2( MESG )
        END IF

    END DO

    !.......  If any process/pollutant names are not found in
    !         master list, exit
    IF( EFLAG ) THEN
        MESG = 'Make sure that master pollutants and/or ' //    &
               'activities files are' //                        &
               CRLF() // BLANK10 // 'consistent with those '//  &
               'used to process the inventory.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Otherwise, convert INVSTAT back from 2/-2 to 1/-1, and 
    !         zero for not used
    ELSE

        INVSTAT = INVSTAT / 2      ! integer math for array

    END IF

    !.......  Allocate memory for list of HAPS
    ALLOCATE( HAPNAM( NHAP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'HAPNAM', PROGNAME )

    !.......  Loop through master list count the actual number of
    !         process/pollutant combinations
    NIPPA = 0
    NHAP = 0
    DO I = 1, MXIDAT
        IF( INVSTAT( I ) .GT. 0 ) THEN
            NIPPA = NIPPA + 1

            IF( INVDVTS( I ) == 'V' .OR.    &
                INVDVTS( I ) == 'T' ) THEN
                NHAP = NHAP + 1
                HAPNAM( NHAP ) = INVDNAM( I )
            END IF
        END IF
    END DO

    !.......  Allocate memory for array of sorted process/pollutant names and
    !         for pollutants only
    ALLOCATE( EANAM( NIPPA ),    &
              EINAM( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EANAM,EINAM', PROGNAME )

    !.......  Initialize all
    EANAM   = ' '       ! array
    EINAM   = ' '       ! array

    !.......  Create array of process/pollutant names matching order of INVTABLE
    !.......  Also store pollutants-only array
    J1 = 0
    DO I = 1, MXIDAT

        IF( INVSTAT( I ) .GT. 0 ) THEN
            J1 = J1 + 1
            EANAM( J1 ) = INVDNAM( I )
            EINAM( J1 ) = EANAM( J1 )
        END IF

    END DO

    !.......  Loop through variable descriptions, find position of
    !         of pollutant for output, and store concatenated position with
    !         species name.  Make sure the same name is not stored twice.
    NCNT = 0
    CALL BUILD_VDESC_UNSORT( NCNT, MNSMATV_L, MSVDESC_L )

    NSMATV = NCNT

    !.......  Make sure some species match some of the pollutants
    IF( NSMATV .EQ. 0 ) THEN
        MESG = 'ERROR: No speciation factors match the inventory'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Sort cumulative list
    CALL SORTIC( NSMATV, INDXA, TVSORTA )

    !.......  Allocate memory for sorted list
    ALLOCATE( TSVDESC( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TSVDESC', PROGNAME )

    !.......  Store sorted list
    DO I = 1, NSMATV
        J = INDXA( I )
        TSVDESC( I ) = TVDESCA( J )
    END DO

    !.......  Create array of sorted unique species, sorted in order of their
    !         associated pollutants, and then in alphabetical order by species...
    !
    !.......  Allocate memory with the number of variables in the speciation
    !         matrices, which will always be >= needed space
    !.......  Also allocate memory for the index between the master species names
    !         and the master pollutant names
    ALLOCATE( EMNAM( NSMATV ),    &
              EMIDX( NSMATV ),    &
          SPCUNIT_L( NSMATV ),    &
          SPCUNIT_S( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMNAM...SPCUNIT_S', PROGNAME )
    EMNAM     = ' '      ! array
    SPCUNIT_L = ' '      ! array
    SPCUNIT_S = ' '      ! array

    !.......  Call subprogram to store species names in appropriate order
    CALL BUILD_SPECIES_ARRAY(  NSMATV, TSVDESC,  NMSPC,  EMNAM )

    !.......  Create index between master species names and master inventory names
    !.......  If there are multiple pollutants per species, the last pollutant
    !         will be put in the index.
    !.......  Also create the speciation units per pollutant
    DO I = 1, NSMATV

        !.......  Get positions of pollutants and model species
        L1 = INDEX( TSVDESC( I ), SPJOIN )
        L2 = LEN_TRIM( TSVDESC( I ) )
        J  = INDEX1( TSVDESC( I )(    1:L1-1 ), NIPPA, EANAM )
        K  = INDEX1( TSVDESC( I )( L1+1:L2   ), NMSPC, EMNAM )

        !.......  Store position of pollutant for each species
        EMIDX( K ) = J

        !.......  Find speciation name in one of the speciation matrices and
        !         set units accordingly.  Set it based on the first one found.
        IF( SPCUNIT_L( K ) .EQ. ' ' ) THEN

            M = INDEX1( TSVDESC( I ), MNSMATV_L, MSVDESC_L )
            IF( M .GT. 0 ) THEN
                SPCUNIT_L( K ) = MSVUNIT_L( M )
                SPCUNIT_S( K ) = MSVUNIT_S( M )
            END IF

        END IF

    END DO

    !.......  Switch speciation matrix between moles and mass based on
    !         the setting of  MRG_GRDOUT_UNIT and MRG_TOTOUT_UNIT
    IF( INDEX( GRDENV, 'mole' ) < 1 ) THEN
        SPCUNIT_L = SPCUNIT_S
    END IF

    !.......  Deallocate temporary arrays
    DEALLOCATE( INVDCOD, INVDNAM, INVSTAT )

    !.......  Abort if error(s) found
    IF( EFLAG ) THEN
        MESG = 'Problem processing variable names.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats.......94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This internal subprogram builds the unsorted list of unique
    !         pollutant-to-species speciation variable descriptions.

    SUBROUTINE BUILD_VDESC_UNSORT( NCNT, NVARS, VDESCS )

        !.......  Subprogram arguments

        INTEGER     , INTENT(IN OUT) :: NCNT                    ! running count
        INTEGER     , INTENT    (IN) :: NVARS                   ! no of var descs
        CHARACTER(*), INTENT    (IN) :: VDESCS( NVARS )         ! spec var descs

        !.......  Local subprogram arrays
        INTEGER,            ALLOCATABLE, SAVE :: EMBIN ( : )
        CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: EMLIST( : )

        !.......  Local subprogram variables
        INTEGER, SAVE :: EMCNT = 0                    ! count of species
        INTEGER          I, K1, K2, K3, L, L2         ! counters and indices
        INTEGER          LJ, LS                       ! string lengths of separators

        CHARACTER(5)           CBIN           ! tmp pollutant bin
        CHARACTER(5)           CPOL           ! tmp pollutant index
        CHARACTER(5)           IBUF           ! tmp variable position to INVDNAM
        CHARACTER(NAMLEN3)     POLNAM         ! tmp pollutant name
        CHARACTER(NAMLEN3)     SPCNAM         ! tmp species name
        CHARACTER(IODLEN3)     TDESC          ! tmp combined pollutant  species

        !----------------------------------------------------------------------

        !.......  Allocate local memory
        IF ( .NOT. ALLOCATED( EMBIN ) ) THEN
            ALLOCATE(  EMBIN( NVARS ),      &
                      EMLIST( NVARS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMLIST', PROGNAME )
            EMBIN  = 0
            EMLIST = ' '
        END IF

        LS = LEN_TRIM( SPJOIN )
        DO I = 1, NVARS

            TDESC  = VDESCS( I )
            L      = INDEX( TDESC, SPJOIN )          ! find species separator
            IF( L .LE. 0 ) THEN
                MESG = 'ERROR: Variable descriptions do not ' //    &
                       'contain proper separator in spec matrix.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            POLNAM = TDESC( 1:L-1 )

            L2     = LEN_TRIM( TDESC )
            SPCNAM = TDESC( L+LS:L2  )                 ! extract species name

            K1 = INDEX1( POLNAM, NIPPA, EINAM )              ! find pol in sorted list
            K2 = INDEX1( TDESC, NCNT, TVDESCA )              ! find combo

            !.......  Look for species name in temporary list. If found,
            !         assign pollutant bin number; if not, initialize pollutant bin
            !         number.
            K3 = INDEX1( SPCNAM, EMCNT, EMLIST )

            IF( K3 .GT. 0 ) THEN
                K3 = EMBIN( K3 )

            ELSE IF ( EMCNT + 1 <= NVARS ) THEN
                EMCNT = EMCNT + 1
                EMLIST( EMCNT ) = SPCNAM
                EMBIN ( EMCNT ) = K1
                K3 = K1

            ELSE
                MESG = 'INTERNAL ERROR: Array overflow for EMBIN.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

            !.......  When pollutant is found, store sorting variable
            !.......  Sorting variable is designed to sort by pollutant (e.g.CO),
            !         then species, and then emission type (e.g. EXH_CO).
            !.......  If multiple pollutants contribute to the same species, the
            !         sorting must ensure the species are together. An example
            !         of this is ROG and VOC creating the same species. The
            !         CBIN variable ensures that this is the case.
            IF( K1 .GT. 0 .AND. K2 .LE. 0 ) THEN
                NCNT = NCNT + 1

                WRITE( CBIN, '(I5.5)' ) K3
                WRITE( CPOL, '(I5.5)' ) K1
                WRITE( IBUF, '(I5.5)' ) I
                INDXA  ( NCNT ) = NCNT
                TVSORTA( NCNT ) = CBIN // SPCNAM // CPOL // IBUF
                TVDESCA( NCNT ) = TDESC

            ELSE IF( K1 .LE. 0 ) THEN
                L1 = LEN_TRIM( TDESC( 1:L-1 ) )
                L2 = LEN_TRIM( SPCNAM )
                MESG = 'NOTE: Speciation factor for "' //&
                        TDESC( 1:L1 ) // '-to-' // TRIM( SPCNAM )//&
                       '" was not associated with ' //&
                       CRLF()// BLANK10// 'an inventory pollutant.'
                CALL M3MSG2( MESG )

            END IF

        END DO

    END SUBROUTINE BUILD_VDESC_UNSORT

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  Build list of species
    !.......  To do this, must only condense the pol-to-species list, in case
    !         multiple pollutants are creating the same species.  Condense by
    !         removing later-appearing duplicates

    SUBROUTINE BUILD_SPECIES_ARRAY( LSMATV,LOCDESC,LMSPC,LOCNAM)

        !.......  Subprogram arguments

        INTEGER     , INTENT (IN) :: LSMATV                     ! input count
        CHARACTER(*), INTENT (IN) :: LOCDESC( LSMATV )          ! spec var descs
        INTEGER     , INTENT(OUT) :: LMSPC                      ! output count
        CHARACTER(*), INTENT(OUT) :: LOCNAM ( LSMATV )          ! species list

        !.......  Local subprogram variables
        INTEGER      I, J, K, N, NCNT, L1, L2              ! counters and indices
        INTEGER      IPTR                                  ! tmp position in array

        !----------------------------------------------------------------------

        !.......  Populate array by looping through master list of pol-to-species in
        !         output order, and if that entry is there for local subroutine
        !         arguments, and if the species hasn't been added yet, then add it.
        NCNT = 0
        DO I = 1, NSMATV

            !.......  Check if master pol-to-species is in local list. If not cycle
            J = INDEX1( TSVDESC( I ), LSMATV, LOCDESC )
            IF( J .LE. 0 ) CYCLE
            N = NSMATV - I

            !.......  Search remaining list of pollutants-to-species for current
            !         iteration's value.
            IPTR = MIN( I+1,NSMATV )
            J = INDEX1( TSVDESC( I ), N, TSVDESC( IPTR ) )

            !.......  See if species is not already in the list of species
            L1 = INDEX( TSVDESC( I ), SPJOIN )
            L2 = LEN_TRIM( TSVDESC( I ) )
            K  = INDEX1( TSVDESC( I )( L1+1:L2 ), NCNT, LOCNAM )

            !.......  If pollutant-to- species is not found, make sure that
            !         species only is not
            IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                NCNT = NCNT + 1
                LOCNAM( NCNT ) = TSVDESC( I )( L1+1:L2 )
            ENDIF
        END DO
        LMSPC = NCNT

    END SUBROUTINE BUILD_SPECIES_ARRAY

END SUBROUTINE MRGVNAMS
