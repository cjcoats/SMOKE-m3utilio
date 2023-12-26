
SUBROUTINE GENUSLST

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine processes generates the source-category-specific lists
!      of inventory characteristics in the MODLISTS module.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!****************************************************************************/
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id$
!
! COPYRIGHT (C) 2005, Environmental Modeling for Policy Development
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

!.........  MODULES for public variables
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC, CIFIP, CSCC, CISIC, CINTGR, CMACT,&
    &                    CORIS, CBLRID, CPDESC, CNAICS, CVTYPE

!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP, NINVSCC, NINVSCL, NINVSIC,&
    &                    NINVSIC2, NINVMACT, NINVORIS,&
    &                    INVCFIP, INVSCC, INVSCL, INVSIC,&
    &                    INVSIC2, INVMACT, INVORIS,&
    &                    INVORFP, IORSMTCH, INVODSC, ORISBLR,&
    &                    OBSRCBG, OBSRCNT, NORISBLR, NOBLRSRC,&
    &                    OBSRCNM, ORISFLAG, NINVNAICS, INVNAICS,&
    &                    NINVVTYP, INVVTYP

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, LSCCEND

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL:: SETSCCTYPE, CHKEXPSCC, CHKEXPSIC

!...........   Sorting index
    INTEGER     INDX( NSRC )

!...........   Local allocateable arrays for ORIS lists
    INTEGER, ALLOCATABLE :: FOIDXA  ( : )  ! sorting index for oris
    INTEGER, ALLOCATABLE :: OBIDXA  ( : )  ! sorting index for oris//blr

!...........   Other local variables
    INTEGER          I, J, J1, L1, L2, N, NS, S
    INTEGER       :: J2 = 0
    INTEGER          IOS                 ! allocate i/o status

    LOGICAL       :: EFLAG    = .FALSE.  ! true: error has occurred
    LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first call to subroutine
    LOGICAL, SAVE :: FIRSTORS = .TRUE.   ! true: first run of ORIS arrays
    LOGICAL          SCCFLAG             ! true: SCC type is different from previous

    CHARACTER(300)     MESG          ! message buffer
    CHARACTER(FIPLEN3) CFIP          ! current cntry/st/co code
    CHARACTER(FIPLEN3) PCFIP        ! previous iteration cntry/st/co
    CHARACTER(VTPLEN3) PVTYP         ! previous vehicle type
    CHARACTER(VTPLEN3) TVTYP         ! tmp vehicle type
    CHARACTER(SCCLEN3) PSCC          ! previous iteration SCC
    CHARACTER(SCCLEN3) PSCCL         ! previous iteration left SCC
    CHARACTER(SCCLEN3) SCCL          ! tmp left SCC
    CHARACTER(SCCLEN3) TSCC          ! tmp SCC
    CHARACTER(SICLEN3) CSIC          ! tmp char SIC
    CHARACTER(SICLEN3) CSIC2         ! tmp 2-digit SIC
    CHARACTER(SICLEN3) PCSIC         ! previous char SIC
    CHARACTER(SICLEN3) PCSIC2        ! previous 2-char SIC
    CHARACTER(MACLEN3) TMACT         ! tmp char MACT code
    CHARACTER(MACLEN3) PMACT         ! previous char MACT code
    CHARACTER(NAILEN3) TNAICS        ! tmp char NAICS code
    CHARACTER(NAILEN3) PNAICS        ! previous char NAICS code
    CHARACTER(BLRLEN3) BLID          ! tmp boiler ID
    CHARACTER(BLRLEN3) PBLID         ! previous boiler ID
    CHARACTER(ORSLEN3) CORS          ! tmp DOE plant ID
    CHARACTER(ORSLEN3) PCORS         ! previous DOE plant ID
    CHARACTER(OBRLEN3) PCORSBLR      ! previous DOE plant ID // boiler
    CHARACTER(DSCLEN3) PDSC          ! tmp plant description

    CHARACTER(16) :: PROGNAME = 'GENUSLST' ! program name

!***********************************************************************
!   begin body of subroutine GENUSLST

!.........  Include all processing in a firstime loop, because this information
!           only needs to be created once per program run, and all of the
!           outputs are stored in MODLISTS

    IF( FIRSTIME ) THEN

        MESG = 'Generating unique lists from inventory data...'
        CALL M3MSG2( MESG )

!.............  Check if CIFIP is allocated.
!.............  If it is, generate unique list of country/state/county codes
        IF( ASSOCIATED( CIFIP ) ) THEN

!.................  Count number of unique codes
            PCFIP = ' '
            J1 = 0
            DO S = 1, NSRC

                CFIP = CIFIP( S )
                IF( CFIP .NE. PCFIP ) J1 = J1 + 1
                PCFIP = CFIP

            END DO
            NINVIFIP = J1

!.................  Allocate memory for country/state/county lists
            ALLOCATE( INVCFIP( NINVIFIP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVCFIP', PROGNAME )
            INVCFIP = ''

!.................  Create unique country/state/county codes list
            PCFIP = ' '
            J1 = 0
            DO S = 1, NSRC

                CFIP = CIFIP( S )
                IF( CFIP .NE. PCFIP ) THEN
                    J1 = J1 + 1
                    INVCFIP( J1 ) = CFIP
                    PCFIP = CFIP
                END IF

            END DO

        END IF   ! End of CIFIP allocated or not

!.............  Check if CVTYPE is allocated
!.............  If it is, generate unique list of vehicle types
        IF( ALLOCATED( CVTYPE ) ) THEN

!.................  Initialize vehicle type sorting index
            DO S = 1, NSRC
                INDX( S ) = S
            END DO

!.................  Sort vehicle types
            CALL SORTIC( NSRC, INDX, CVTYPE )

!.................  Count number of unique vehicle types
            PVTYP = '-9'
            J1 = 0
            DO S = 1, NSRC
                J = INDX( S )
                TVTYP = CVTYPE( J )

                IF( TVTYP /= PVTYP ) J1 = J1 + 1

                PVTYP = TVTYP
            END DO
            NINVVTYP = J1

!.................  Allocate memory for vehicle type list
            ALLOCATE( INVVTYP( NINVVTYP ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVVTYP', PROGNAME )

!.................  Create list of unique vehicle types
            PVTYP = '-9'
            J1 = 0
            DO S = 1, NSRC
                J = INDX( S )
                TVTYP = CVTYPE( J )

                IF( TVTYP /= PVTYP ) THEN
                    J1 = J1 + 1
                    INVVTYP( J1 ) = TVTYP
                    PVTYP = TVTYP
                END IF
            END DO
        END IF  ! end vehicle type processing

!.............  Check if CSCC is allocated.
!.............  If it is, generate unique list of SCCs and left SCCs
        IF( ASSOCIATED( CSCC ) ) THEN

!.................  Initialize SCC sorting index
            DO S = 1, NSRC
                INDX( S ) = S
            END DO

!.................  Sort all SCCs in the point sources inventory in increasing
!                   order
            CALL SORTIC( NSRC, INDX, CSCC )

!.................  Count number of unique SCCs
            PSCC = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TSCC = CSCC( J )

                IF( TSCC .NE. PSCC ) J1 = J1 + 1

                PSCC = TSCC

            END DO
            NINVSCC = J1

!.................  Allocate memory for SCC lists
            ALLOCATE( INVSCC( NINVSCC ),&
            &          INVSCL( NINVSCC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVSCC,INVSCL', PROGNAME )

!.................  Create unique SCCs lists
            PSCC = '-9'
            PSCCL = '-9'
            J1 = 0        ! to skip TSCC = 0
            J2 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TSCC  = CSCC( J )

!.....................  Set type of SCC
                SCCFLAG = SETSCCTYPE( TSCC )
                SCCL  = TSCC( 1:LSCCEND )

                IF( TSCC .NE. PSCC ) THEN
                    J1 = J1 + 1
                    INVSCC( J1 ) = TSCC
                    PSCC = TSCC
                END IF

!.....................  Don't include expanded SCCs in left SCC list
                IF( CHKEXPSCC( TSCC ) ) CYCLE

                IF( SCCL .NE. PSCCL ) THEN
                    J2 = J2 + 1
                    INVSCL( J2 ) = SCCL
                    PSCCL = SCCL
                END IF

            END DO

        END IF   ! End SCC processing
        NINVSCL =  J2

!.............  Check if CISIC is allocated.
!.............  If it is, generate unique list of SICs
        IF( ASSOCIATED( CISIC ) ) THEN

!.................  Initialize SIC sorting index
            DO S = 1, NSRC
                INDX( S ) = S
            END DO

!.................  Sort all SICs in the inventory in increasing order
            CALL SORTIC( NSRC, INDX, CISIC )

!.................  Count number of unique SICs
            PCSIC = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                CSIC = CISIC( J )

                IF( CSIC .NE. PCSIC ) J1 = J1 + 1

                PCSIC = CSIC

            END DO
            NINVSIC = J1

!.................  Allocate memory for SIC lists
            ALLOCATE( INVSIC( NINVSIC ),&
            &         INVSIC2( NINVSIC ), STAT=IOS )
            CALL CHECKMEM( IOS,' INVSIC,INVSIC2', PROGNAME )

!.................  Create unique SIC lists
            PCSIC = '-9'
            PCSIC2 = '-9'
            J1 = 0
            J2 = 0
            DO S = 1, NSRC

                J = INDX( S )

                CSIC = CISIC( J )
                CSIC2 = CSIC( 1:SICLEN3-2 )

                IF( CSIC .NE. PCSIC ) THEN
                    J1 = J1 + 1
                    INVSIC( J1 ) = CSIC
                    PCSIC = CSIC
                END IF

!.....................  Don't include expanded SICs in 2-digit SIC list
                IF( CHKEXPSIC( CSIC ) ) CYCLE

                IF( CSIC2 .NE. PCSIC2 ) THEN
                    J2 = J2 + 1
                    INVSIC2( J2 ) = CSIC2
                    PCSIC2 = CSIC2
                END IF

            END DO

        END IF   ! End SIC processing
        NINVSIC2 = J2

!.............  Check if CMACT is allocated.
!.............  If it is, generate unique list of MACT codes
        IF( ASSOCIATED( CMACT ) ) THEN

!.................  Initialize MACT sorting index
            DO S = 1, NSRC
                INDX( S ) = S
            END DO

!.................  Sort all MACTs in the inventory in increasing order
            CALL SORTIC( NSRC, INDX, CMACT )

!.................  Count number of unique MACTs
            PMACT = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TMACT = CMACT( J )

                IF( TMACT /= PMACT ) J1 = J1 + 1

                PMACT = TMACT

            END DO
            NINVMACT = J1

!.................  Allocate memory for MACT lists
            ALLOCATE( INVMACT( NINVMACT ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVMACT', PROGNAME )

!.................  Create unique MACTs list
            PMACT = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TMACT = CMACT( J )

                IF( TMACT /= PMACT ) THEN
                    J1 = J1 + 1
                    INVMACT( J1 ) = TMACT
                    PMACT = TMACT
                END IF

            END DO

        END IF    ! End MACT processing

!.............  Check if CNAICS is allocated.
!.............  If it is, generate unique list of NAICS codes
        IF( ASSOCIATED( CNAICS ) ) THEN

!.................  Initialize NAICS sorting index
            DO S = 1, NSRC
                INDX( S ) = S
            END DO

!.................  Sort all NAICS in the inventory in increasing order
            CALL SORTIC( NSRC, INDX, CNAICS )

!.................  Count number of unique NAICS
            PNAICS = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TNAICS = CNAICS( J )

                IF( TNAICS /= PNAICS ) J1 = J1 + 1

                PNAICS = TNAICS

            END DO
            NINVNAICS = J1

!.................  Allocate memory for NAICS lists
            ALLOCATE( INVNAICS( NINVNAICS ), STAT=IOS )
            CALL CHECKMEM( IOS, 'INVNAICS', PROGNAME )

!.................  Create unique NAICS list
            PNAICS = '-9'
            J1 = 0
            DO S = 1, NSRC

                J = INDX( S )

                TNAICS = CNAICS( J )

                IF( TNAICS /= PNAICS ) THEN
                    J1 = J1 + 1
                    INVNAICS( J1 ) = TNAICS
                    PNAICS = TNAICS
                END IF

            END DO

        END IF    ! End NAICS processing

        FIRSTIME = .FALSE.

    END IF  ! End of firstime

!.........  Create list of FIPs, ORIS IDs and boiler IDs from the inventory
!           and how many sources for each...
    IF ( ORISFLAG             .AND.&
    &     FIRSTORS             .AND.&
    &     ALLOCATED ( CORIS  ) .AND.&
    &     ALLOCATED ( CBLRID )       ) THEN

        MESG = 'Generating ORIS lists...'
        CALL M3MSG2( MESG )

!.............  First, count the number of unique records
        NINVORIS = 0
        NORISBLR = 0
        NOBLRSRC = 0
        PCORS    = ' '
        PBLID    = ' '
        DO S = 1, NSRC

            CORS = CORIS ( S )
            BLID = CBLRID( S )

!.................  Skip missing ORIS IDs
            IF ( CORS == ' ' ) CYCLE

!.................  Count unique ORIS IDs; ORIS IDs are not guaranteed to be in
!                   blocks, but this will count the maximum number
            IF ( CORS /= PCORS ) THEN
                NINVORIS = NINVORIS + 1
            END IF

!.................  Skip missing boiler IDs
            IF ( BLID == ' ' ) CYCLE

            NOBLRSRC = NOBLRSRC + 1

!.................  Count unique ORIS // boiler combos
            IF ( CORS /= PCORS .OR.&
            &     BLID /= PBLID      ) THEN
                NORISBLR = NORISBLR + 1
            END IF

            PCORS = CORS
            PBLID = BLID

        END DO

!.............  Allocate memory for unsorted ORIS lists
        ALLOC TE( FOIDXA( NINVORIS ),&
        &         INVORFP( NINVORIS ),&
        &         INVODSC( NINVORIS ),&
        &        IORSMTCH( NINVORIS ),&
        &        INVORISA( NINVORIS ),&
        &        INVORFPA( NINVORIS ),&
        &        INVODSCA( NINVORIS ),&
        &          OBIDXA( NORISBLR ),&
        &         ORISBLR( NORISBLR ),&
        &         OBSRCNT( NORISBLR ),&
        &         OBSRCBG( NORISBLR ),&
        &         OBSRCNM( NOBLRSRC ),  STAT=IOS )
        CALL CHECKMEM( IOS, 'FOIDXA...OBSRCNM', PROGNAME )

!.............  Store arrays
        INVODSCA = ' '
        NINVORIS = 0
        NORISBLR = 0
        PCORS    = ' '
        PBLID    = ' '
        DO S = 1, NSRC

            CFIP = CIFIP ( S )
            CORS = CORIS ( S )
            BLID = CBLRID( S )
            PDSC = CPDESC( S )

!.................  Skip missing ORIS IDs
            IF ( CORS == ' ' ) CYCLE

!.................  Unsorted ORIS arrays
            IF ( CORS /= PCORS ) THEN
                J = INDEX1( CORS, NINVORIS, INVORISA )

                IF( J <= 0 ) THEN
                    NINVORIS = NINVORIS + 1
                    FOIDXA ( NINVORIS ) = NINVORIS
                    INVORIS( NINVORIS ) = CORS
                    INVORFP( NINVORIS ) = CFIP
                    INVODSC( NINVORIS ) = PDSC
                ELSE
                    IF( INVORFP( J ) /= CFIP ) THEN
                        MESG = 'WARNING: Different FIPS codes ' //&
                        &      'found for ORIS ID ' // CORS&
                        &      // '.  Will use ' // CFIP //&
                        &      ' for reporting.'
                        CALL M3MESG( MESG )
                    END IF

                    IF( INVODSCA( J ) /= PDSC ) THEN
                        MESG = 'WARNING: Different plant ' //&
                        &  'descriptions found for ORIS ID ' // CORS&
                        &  // '.  Will use ' // INVODSCA( J ) //&
                        &  ' for reporting.'
                        CALL M3MESG( MESG )
                    END IF
                END IF
            END IF

!.................  Skip missing boiler IDs
            IF ( BLID == ' ' ) CYCLE

!.................  Unsorted oris/boiler array
            J = INDEX1( CORS // BLID, NORISBLR, ORISBLRA )

            IF( J <= 0 ) THEN
                NORISBLR = NORISBLR + 1
                OBIDXA ( NORISBLR ) = NORISBLR
                ORISBLR( NORISBLR ) = CORS // BLID
                OBSRCNT( NORISBLR ) = 1
            ELSE
                OBSRCNT( J ) = OBSRCNT( J ) + 1
            END IF

            PCORS = CORS
            PBLID = BLID

        END DO

        IF( EFLAG ) THEN
            MESG = 'Problem counting ORIS IDs'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Sort arrays
        CALL SORTI(   NINVORIS, FOIDXA, INVORIS )
        CALL PERMUTI( NINVORIS, FOIDXA, INVORIS, INVORFP, INVODSC )

        CALL SORTI(   NORISBLR, OBIDXA, ORISBLR )
        CALL PERMUTI( NORISBLR, OBIDXA, ORISBLR, OBSRCNT )

!.............  Build ORIS/boiler source start array
        OBSRCBG( 1 ) = 1
        DO I = 2, NORISBLR
            OBSRCBG( I ) = OBSRCBG( I-1 ) + OBSRCNT( I-1 )
        END DO

!.............  Create list of sources for each ORIS/boiler combo
        OBSRCNM = 0

        DO S = 1, NSRC

            CORS = CORIS ( S )
            BLID = CBLRID( S )

!.................  Skip missing ORIS or boiler IDs
            IF( CORS == ' ' .OR. BLID == ' ' ) CYCLE

!.................  Find combination in sorted list
            J = FINDC( CORS // BLID, NORISBLR, ORISBLR )

            N  = OBSRCBG( J )
            NS = N + OBSRCNT( J )
            DO
                IF( OBSRCNM( N ) == 0 ) THEN
                    OBSRCNM( N ) = S
                    EXIT
                ELSE
                    N = N + 1
                END IF
            END DO

        END DO
        FIRSTORS = .FALSE.

    END IF   ! End boiler processing

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94300 FORMAT( A, I2.2, A, I2.2, A )

END SUBROUTINE GENUSLST
