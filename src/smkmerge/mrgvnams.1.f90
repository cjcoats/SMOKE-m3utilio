
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
!       Created 2/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!*************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,&
    &                    ARFLAG, MRFLAG, PRFLAG, SFLAG, PDEV,&
    &                    ANIPOL, BNIPOL, PNIPOL, NIPOL,&
    &                    MNIPPA, NIPPA,&
    &                    AEINAM, BEINAM, MEANAM, PEINAM, EINAM,&
    &                    AEMNAM, BEMNAM, MEMNAM, PEMNAM, EMNAM,&
    &                    EANAM, EMIDX, SPCUNIT,&
    &                    AONAMES, BONAMES, MONAMES, PONAMES, TONAMES,&
    &                    AOUNITS, BOUNITS, MOUNITS, POUNITS, TOUNITS,&
    &                    ANSMATV, BNSMATV, MNSMATV, PNSMATV, NSMATV,&
    &                    ARNMSPC,          MRNMSPC, PRNMSPC,&
    &                    ANRMATV,          MNRMATV, PNRMATV,&
    &                    ASVDESC, BSVDESC, MSVDESC, PSVDESC, TSVDESC,&
    &                    ANMSPC,  BNMSPC,  MNMSPC,  PNMSPC,  NMSPC,&
    &                    ARVDESC,          MRVDESC, PRVDESC,&
    &                    ASVUNIT,          MSVUNIT, PSVUNIT

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: MXIDAT, INVDNAM, INVSTAT, INVDCOD

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Temporary arrays for building sorted pol-to-species names list
    INTEGER                            MXPTOSPC     ! max for dimension
    INTEGER           , ALLOCATABLE :: INDXA  ( : ) ! sorting index
    CHARACTER(PLSLEN3), ALLOCATABLE :: TVSORTA( : ) ! basis for sorting
    CHARACTER(PLSLEN3), ALLOCATABLE :: TVDESCA( : ) ! pol-to-spec concat

!...........   Temporary arrays for sorting biogenic species names
    INTEGER           , ALLOCATABLE :: TMPBIDXA( : ) ! position in EMNAM
    INTEGER           , ALLOCATABLE :: TMPBIDXB( : ) ! sorting index
    CHARACTER(NAMLEN3), ALLOCATABLE :: TMPBNAM ( : ) ! extra spc names

!...........   Other local variables
    INTEGER         I, J, J1, J2, L, K, L1, L2, M, N  ! counters and indices
    INTEGER         LJ              ! string length of emis types joiner
    INTEGER         IOS             ! i/o error status
    INTEGER         JA, JM, JP      ! position in rctvty var desc of pol-spc
    INTEGER         MXSMATV         ! maximum size for unsorted pol-to-spc list
    INTEGER         NCNT            ! counter
    INTEGER      :: NIACT  = 0      ! actual no activites

    LOGICAL      :: EFLAG = .FALSE. ! error flag

    CHARACTER(NAMLEN3)   CPOL     ! tmp pol/act buffer
    CHARACTER(NAMLEN3)   BNUM     ! tmp biogenic units numerator
    CHARACTER(300)       MESG     ! message buffer

    CHARACTER(16) :: PROGNAME = 'MRGVNAMS' ! program name

!***********************************************************************
!   begin body of subroutine MRGVNAMS

!.........  Read, sort, and store pollutant codes/names file
    CALL RDCODNAM( PDEV )

!.........  Loop through pollutants and/or activities in each source category
!           and update status of pollutants and activities in master list.
    IF( AFLAG )&
    &  CALL SET_MASTER_POL_STAT( 'area'    , ANIPOL, AEINAM, EFLAG )
    IF( BFLAG )&
    &  CALL SET_MASTER_POL_STAT( 'biogenic', BNIPOL, BEINAM, EFLAG )
    IF( MFLAG )&
    &  CALL SET_MASTER_POL_STAT( 'mobile'  , MNIPPA, MEANAM, EFLAG )
    IF( PFLAG )&
    &  CALL SET_MASTER_POL_STAT( 'point'   , PNIPOL, PEINAM, EFLAG )

!.........  If any pollutants/activities from matrices not found in master
!           list, exit
    IF( EFLAG ) THEN
        MESG = 'Make sure that master pollutants and/or ' //&
        &       'activities files are' // CRLF() // BLANK10 //&
        &       'consistent with those used to process the ' //&
        &       'inventory.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!.........  Otherwise, convert INVSTAT back from 2/-2 to 1/-1, and zero for
!           not used
    ELSE

        INVSTAT = INVSTAT / 2  ! integer math for array

    END IF

!.........  Loop through master list count the actual number of total pollutants
!           and activities
    NIPOL = 0
    NIACT = 0
    DO I = 1, MXIDAT
        IF( INVSTAT( I ) .GT. 0 ) NIPOL = NIPOL + 1
        IF( INVSTAT( I ) .LT. 0 ) NIACT = NIACT + 1
    END DO
    NIPPA = NIPOL + NIACT

!.........  Allocate memory for array of sorted pollutants and activities and
!           for pollutants only, and for the input variable names and units
!.........  NOTE - the input variable names could be different from EANAM if
!           user is using emissions from the inventory file, and has selected
!           average day emissions instead of annual emissions.
    ALLOCATE( EANAM( NIPPA ),&
    &          EINAM( NIPOL ),&
    &        TONAMES( NIPPA ),&
    &        TOUNITS( NIPPA ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TOUNITS', PROGNAME )

!.........  Initialize all
    EANAM   = ' '   ! array
    EINAM   = ' '   ! array
    TONAMES = ' '   ! array
    TOUNITS = ' '   ! array

!.........  Create array of sorted unique pollutants and activities
!.........  Also determine number of pollutants and store pollutants-only array
    J1 = 0
    J2 = 0
    LJ = LEN_TRIM( ETJOIN )
    DO I = 1, MXIDAT

        IF( INVSTAT( I ) .NE. 0 ) THEN
            J1 = J1 + 1
            EANAM( J1 ) = INVDNAM( I )

            K  = INDEX( EANAM( J1 ), ETJOIN )
            L2 = LEN_TRIM( EANAM( J1 ) )
            IF( K .GT. 0 ) THEN
                CPOL = EANAM( J1 )( K+LJ:L2 )
            ELSE
                CPOL = EANAM( J1 )
            END IF
        END IF

        IF( INVSTAT( I ) .GT. 0 ) THEN
            J2 = J2 + 1
            EINAM( J2 ) = CPOL
        END IF

    END DO

!.........  Create array of input variable names and units that combines
!           this information from anthropogenic source categories.
    IF( AFLAG )&
    &  CALL BUILD_NAMES_UNITS( ANIPOL, AEINAM, AONAMES, AOUNITS,&
    &                          TONAMES, TOUNITS )
    IF( BFLAG )&
    &  CALL BUILD_NAMES_UNITS( BNIPOL, BEINAM, BONAMES, BOUNITS,&
    &                          TONAMES, TOUNITS )
    IF( MFLAG )&
    &  CALL BUILD_NAMES_UNITS( MNIPPA, MEANAM, MONAMES, MOUNITS,&
    &                          TONAMES, TOUNITS )
    IF( PFLAG )&
    &  CALL BUILD_NAMES_UNITS( PNIPOL, PEINAM, PONAMES, POUNITS,&
    &                          TONAMES, TOUNITS )

!.........  Create array of sorted unique pol-to-species, sorted in order of
!           pollutants, and then in alphabetical order by species...

!.........  Allocate maximum possible memory needed by summing number of
!           variables from all source categories
!.........  Allocate memory for sorted list

    MXPTOSPC = ANSMATV + ARNMSPC + BNSMATV + MNSMATV + MRNMSPC +&
    &           PRNMSPC + PNSMATV

    ALLOCATE( INDXA( MXPTOSPC ),&
    &        TVSORTA( MXPTOSPC ),&
    &        TVDESCA( MXPTOSPC ),&
    &        TSVDESC( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TSVDESC', PROGNAME )

!.........  Loop through source-specific variable descriptions, find position of
!           of pollutant for output, and store concatenated position with
!           species name.  Make sure the same name is not stored twice.
    NCNT = 0
    JA   = ANRMATV - ARNMSPC + 1
    JM   = MNRMATV - MRNMSPC + 1
    JP   = PNRMATV - PRNMSPC + 1
    MXSMATV = ANSMATV + BNSMATV + MNSMATV + PNSMATV
    IF( ARFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, ARNMSPC, MXSMATV, ARVDESC(JA) ) ! area rctvty
    IF( MRFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, MRNMSPC, MXSMATV, MRVDESC(JM) ) ! mobile rctvty
    IF( PRFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, PRNMSPC, MXSMATV, PRVDESC(JP) ) ! point rctvty
    IF( AFLAG .AND. SFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, ANSMATV, MXSMATV, ASVDESC ) ! area spectn
    IF( BFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, BNSMATV, MXSMATV, BSVDESC ) ! biogenics
    IF( MFLAG .AND. SFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, MNSMATV, MXSMATV, MSVDESC ) ! mobile spectn
    IF( PFLAG .AND. SFLAG )&
    &  CALL BUILD_VDESC_UNSORT( NCNT, PNSMATV, MXSMATV, PSVDESC ) ! point spectn

    NSMATV = NCNT

!.........  Make sure some species match some of the pollutants
    IF( SFLAG .AND. NSMATV .EQ. 0 ) THEN
        MESG = 'ERROR: No speciation factors match the inventory!'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Sort cumulative list
    CALL SORTIC( NSMATV, INDXA, TVSORTA )

!.........  Store sorted list
    DO I = 1, NSMATV
        J = INDXA( I )
        TSVDESC( I ) = TVDESCA( J )
    END DO

!.........  Create array of sorted unique species, sorted in order of their
!           associated pollutants, and then in alphabetical order by species...
!
!.........  Allocate memory with the number of variables in the speciation
!           matrices, which will always be >= needed space
!.........  Also allocate memory for the index between the master species names
!           and the master pollutant names
    ALLOCATE( AEMNAM( ANSMATV ),&
    &          BEMNAM( BNSMATV ),&
    &          MEMNAM( MNSMATV ),&
    &          PEMNAM( PNSMATV ),&
    &           EMNAM( NSMATV ),&
    &           EMIDX( NSMATV ),&
    &         SPCUNIT( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPCUNIT', PROGNAME )
    AEMNAM  = ' '  ! array
    BEMNAM  = ' '  ! array
    MEMNAM  = ' '  ! array
    PEMNAM  = ' '  ! array
    EMNAM   = ' '  ! array
    SPCUNIT = ' '  ! array

!.........  Call subprogram to store species names in appropriate order
!           for all source categories and the total.
    IF( AFLAG .AND. SFLAG )&
    &  CALL BUILD_SPECIES_ARRAY( ANSMATV, ASVDESC, ANMSPC, AEMNAM )
    IF( BFLAG )&
    &  CALL BUILD_SPECIES_ARRAY( BNSMATV, BSVDESC, BNMSPC, BEMNAM )
    IF( MFLAG .AND. SFLAG )&
    &  CALL BUILD_SPECIES_ARRAY( MNSMATV, MSVDESC, MNMSPC, MEMNAM )
    IF( PFLAG .AND. SFLAG )&
    &  CALL BUILD_SPECIES_ARRAY( PNSMATV, PSVDESC, PNMSPC, PEMNAM )
    IF( SFLAG )&
    &  CALL BUILD_SPECIES_ARRAY(  NSMATV, TSVDESC, NMSPC, EMNAM )

!.........  Resort biogenic species names to be consistent with global species
!           list
    IF( BFLAG ) THEN
        ALLOCATE( TMPBNAM( BNMSPC ),&
        &         TMPBIDXA( BNMSPC ),&
        &         TMPBIDXB( BNMSPC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'TMPBIDXB', PROGNAME )

        DO I = 1, BNMSPC
            J = INDEX1( BEMNAM( I ), NMSPC, EMNAM )
            TMPBIDXA( I ) = J
            TMPBIDXB( I ) = I
        END DO

        CALL SORTI1( BNMSPC, TMPBIDXB, TMPBIDXA )

        TMPBNAM = BEMNAM  ! array

        DO I = 1, BNMSPC
            J = TMPBIDXB( I )
            BEMNAM( I ) = TMPBNAM( J )
        END DO

    END IF

!.........  Create index between master species names and master inventory names
!.........  If there are multiple pollutants per species, the last pollutant
!           will be put in the index.
!.........  Also create the speciation units per pollutant
    IF( SFLAG ) THEN

!.............  Create temporary name for numerator of biogenic units
!.............  Assumes that the units will be the same for the biogenic input
        IF ( BFLAG ) THEN
            L2 = INDEX( BOUNITS( 1 ), '/' )
            BNUM = BOUNITS( 1 )( 1:L2-1 )
        END IF

        DO I = 1, NSMATV

!.................  Get positions of pollutants and model species
            L1 = INDEX( TSVDESC( I ), SPJOIN )
            L2 = LEN_TRIM( TSVDESC( I ) )
            J  = INDEX1( TSVDESC( I )(    1:L1-1 ), NIPPA, EANAM )
            K  = INDEX1( TSVDESC( I )( L1+1:L2   ), NMSPC, EMNAM )

!.................  Store position of pollutant for each species
            EMIDX( K ) = J

!.................  Find speciation name in one of the speciation matrices and
!                   set units accordingly.  Set it based on the first one found.
            IF( SPCUNIT( K ) .EQ. ' ' ) THEN

                IF( AFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), ANSMATV, ASVDESC )
                    IF( M .GT. 0 ) SPCUNIT( K ) = ASVUNIT( M )
                END IF

                IF( BFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), BNSMATV, BSVDESC )
                    L = LEN_TRIM( BNUM )
                    IF( M .GT. 0 ) SPCUNIT( K ) = BNUM( 1:L ) // '/'&
                    &                           // BNUM( 1:L )
                END IF

                IF( MFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), MNSMATV, MSVDESC )
                    IF( M .GT. 0 ) SPCUNIT( K ) = MSVUNIT( M )
                END IF

                IF( PFLAG ) THEN
                    M = INDEX1( TSVDESC( I ), PNSMATV, PSVDESC )
                    IF( M .GT. 0 ) SPCUNIT( K ) = PSVUNIT( M )
                END IF

            END IF

        END DO

    END IF      ! End if speciation

!.........  Deallocate temporary arrays
    DEALLOCATE( INDXA, TVSORTA, TVDESCA, INVDCOD, INVDNAM, INVSTAT )
    IF( BFLAG ) DEALLOCATE( TMPBNAM, TMPBIDXA, TMPBIDXB )

!.........  Abort if error(s) found
    IF( EFLAG ) THEN
        MESG = 'Problem processing variable names.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram updates the status for the master
!               list of inventory pollutants
    SUBROUTINE SET_MASTER_POL_STAT( CATDESC, N, PNAMES, EFLAG )

!.............  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: CATDESC     ! src category desc
        INTEGER     , INTENT (IN) :: N           ! number of src-spec pols
        CHARACTER(*), INTENT (IN) :: PNAMES( N ) ! names of src-spec pols
        LOGICAL     , INTENT(OUT) :: EFLAG       ! error flag

!.............  Local subprogram variables
        INTEGER      I, J, L      ! counters and indices

        CHARACTER(NAMLEN3) VBUF

!----------------------------------------------------------------------

        DO I = 1, N

            VBUF = PNAMES( I )
            J = INDEX1( VBUF, MXIDAT, INVDNAM )

!................. Adjust data status if actual pollutant or activity is being
!                  input to Smkmerge
            IF( J .GT. 0 ) THEN

                INVSTAT( J ) = INVSTAT( J ) * 2

            ELSE

                EFLAG = .TRUE.
                L = LEN_TRIM( VBUF )
                MESG = 'ERROR: Variable "' // VBUF( 1:L ) //&
                &       '" from ' // CATDESC //&
                &       ' inventory file is not in ' // CRLF() //&
                &       BLANK10 // 'inventory table.'
                CALL M3MSG2( MESG )

            END IF

        END DO

    END SUBROUTINE SET_MASTER_POL_STAT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.............  This internal subprogram builds the sorted list of input
!               emissions/activity variable names and their units
    SUBROUTINE BUILD_NAMES_UNITS( NINVARS, REFNAMS, INNAMS,&
    &                              INUNITS, OUTNAMS, OUTUNITS )

!.............  Subprogram arguments

        INTEGER     , INTENT (IN) :: NINVARS            ! no. input vars
        CHARACTER(*), INTENT (IN) :: REFNAMS( NINVARS ) ! reference names
        CHARACTER(*), INTENT (IN) :: INNAMS ( NINVARS ) ! input var names
        CHARACTER(*), INTENT (IN) :: INUNITS( NINVARS ) ! input units
        CHARACTER(*), INTENT(OUT) :: OUTNAMS( NIPPA   ) ! out variable names
        CHARACTER(*), INTENT(OUT) :: OUTUNITS( NIPPA  ) ! out units

!.............  Local subprogram variables
        INTEGER        K, L, V     ! counters and indices
        CHARACTER(300) MESG        ! message buffer

!----------------------------------------------------------------------

!.............  Loop through input variable names
        DO V = 1, NINVARS

!.................  Look for reference variable name in master list

            K = INDEX1( REFNAMS( V ), NIPPA, EANAM )

!.................  When output names and units have already been set, check
!                   for consistency.
            IF( OUTNAMS( K ) .NE. ' ' ) THEN

                IF( OUTUNITS( K ) .NE. INUNITS( V ) ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Units for input variable "' //&
                    &       TRIM( INNAMS( V ) ) // '" are not ' //&
                    &       'consistent among the input files.'
                    CALL M3MSG2( MESG )

                END IF

!.................  Otherwise, set the output names
            ELSE
                OUTNAMS ( K ) = INNAMS ( V )
                OUTUNITS( K ) = INUNITS( V )
            END IF

        END DO

    END SUBROUTINE BUILD_NAMES_UNITS

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.............  This internal subprogram builds the unsorted list of unique
!               pollutant-to-species speciation variable descriptions.
    SUBROUTINE BUILD_VDESC_UNSORT( NCNT, NVARS, MXTMP, VDESCS )

        INCLUDE 'PARMS3.EXT'    !  I/O API parameters

!.............  Subprogram arguments

        INTEGER     , INTENT(IN OUT) :: NCNT            ! running count
        INTEGER     , INTENT    (IN) :: NVARS           ! no of var descs
        INTEGER     , INTENT    (IN) :: MXTMP           ! max no of var descs
        CHARACTER(*), INTENT    (IN) :: VDESCS( NVARS ) ! spec var descs

!.............  Local subprogram arrays
        INTEGER,            ALLOCATABLE, SAVE :: EMBIN ( : )
        CHARACTER(NAMLEN3), ALLOCATABLE, SAVE :: EMLIST( : )

!.............  Local subprogram variables
        INTEGER, SAVE :: EMCNT = 0            ! count of species
        INTEGER          I, K1, K2, K3, L, L2 ! counters and indices
        INTEGER          LJ, LS               ! string lengths of separators

        CHARACTER(5)           CBIN   ! tmp pollutant bin
        CHARACTER(5)           CPOL   ! tmp pollutant index
        CHARACTER(5)           IBUF   ! tmp variable position to INVDNAM
        CHARACTER(NAMLEN3)     POLNAM ! tmp pollutant name
        CHARACTER(NAMLEN3)     SPCNAM ! tmp species name
        CHARACTER(IODLEN3)     TDESC  ! tmp combined pollutant & species

!----------------------------------------------------------------------

!.............  Allocate local memory
        IF ( .NOT. ALLOCATED( EMBIN ) ) THEN
            ALLOCATE( EMBIN( MXTMP ),&
            &  EMLIST( MXTMP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'EMLIST', PROGNAME )
            EMBIN  = 0
            EMLIST = ' '
        END IF

        LJ = LEN_TRIM( ETJOIN )
        LS = LEN_TRIM( SPJOIN )
        DO I = 1, NVARS

            TDESC  = VDESCS( I )
            L      = INDEX( TDESC, SPJOIN )  ! find species separator
            IF( L .LE. 0 ) THEN
                MESG = 'ERROR: Variable descriptions do not ' //&
                &       'contain proper separator in spec matrix.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            L2     = INDEX( TDESC, ETJOIN )  ! find emission type separator

!.................  Extract pollutant name (without emissions type, if applies)
            IF( L2 .GT. 0 ) THEN
                POLNAM = TDESC( L2+LJ:L-1 )
            ELSE
                POLNAM = TDESC( 1:L-1 )
            END IF

            L2     = LEN_TRIM( TDESC )
            SPCNAM = TDESC( L+LS:L2  )     ! extract species name

            K1 = INDEX1( POLNAM, NIPOL, EINAM )  ! find pol in sorted list
            K2 = INDEX1( TDESC, NCNT, TVDESCA )  ! find combo

!.................  Look for species name in temporary list. If found, assign
!                   pollutant bin number, if not, initialize pollutant bin
!                   number.
            K3 = INDEX1( SPCNAM, EMCNT, EMLIST )
            IF( K3 .GT. 0 ) THEN
                K3 = EMBIN( K3 )

            ELSE IF ( EMCNT + 1 <= MXTMP ) THEN
                EMCNT = EMCNT + 1
                EMLIST( EMCNT ) = SPCNAM
                EMBIN ( EMCNT ) = K1
                K3 = K1

            ELSE
                MESG = 'INTERNAL ERROR: Array overflow for EMBIN.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END IF

!.................  When pollutant is found, store sorting variable
!.................  Sorting variable is designed to sort by pollutant (e.g.CO),
!                   then species, and then emission type (e.g. EXH_CO).
!.................  If multiple pollutants contribute to the same species, the
!                   sorting must ensure the species are together. An example
!                   of this is ROG and VOC creating the same species. The
!                   CBIN variable ensures that this is the case.
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
                &        TDESC( 1:L1 ) // '-to-' // SPCNAM( 1:L2 )//&
                &       '" was not associated with ' //&
                &       CRLF()// BLANK10// 'an inventory pollutant.'
                CALL M3MSG2( MESG )

            END IF

        END DO

    END SUBROUTINE BUILD_VDESC_UNSORT

!----------------------------------------------------------------------
!----------------------------------------------------------------------

!.........  Build list of species
!.........  To do this, must only condense the pol-to-species list, in case
!           multiple pollutants are creating the same species.  Condense by
!           removing later-appearing duplicates
    SUBROUTINE BUILD_SPECIES_ARRAY( LSMATV,LOCDESC,LMSPC,LOCNAM)

!.............  Subprogram arguments

        INTEGER     , INTENT (IN) :: LSMATV             ! input count
        CHARACTER(*), INTENT (IN) :: LOCDESC( LSMATV )  ! spec var descs
        INTEGER     , INTENT(OUT) :: LMSPC              ! output count
        CHARACTER(*), INTENT(OUT) :: LOCNAM ( LSMATV )  ! species list

!.............  Local subprogram variables
        INTEGER      I, J, K, N, NCNT, L1, L2      ! counters and indices
        INTEGER      IPTR                          ! tmp position in array

!----------------------------------------------------------------------

!.............  Populate array by looping through master list of pol-to-species in
!               output order, and if that entry is there for local subroutine
!               arguments, and if the species hasn't been added yet, then add it.
        NCNT = 0
        DO I = 1, NSMATV

!.................  Check if master pol-to-species is in local list. If not cycle
            J = INDEX1( TSVDESC( I ), LSMATV, LOCDESC )
            IF( J .LE. 0 ) CYCLE

            N = NSMATV - I

!.................  Search remaining list of pollutants-to-species for current
!                   iteration's value.
            IPTR = MIN( I+1,NSMATV )
            J = INDEX1( TSVDESC( I ), N, TSVDESC( IPTR ) )

!.................  See if species is not already in the list of species
            L1 = INDEX( TSVDESC( I ), SPJOIN )
            L2 = LEN_TRIM( TSVDESC( I ) )
            K  = INDEX1( TSVDESC( I )( L1+1:L2 ), NCNT, LOCNAM )

!.................  If pollutant-to- species is not found, make sure that
!                   species only is not
            IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                NCNT = NCNT + 1
                LOCNAM( NCNT ) = TSVDESC( I )( L1+1:L2 )
            ENDIF
        END DO
        LMSPC = NCNT

    END SUBROUTINE BUILD_SPECIES_ARRAY

END SUBROUTINE MRGVNAMS
