
SUBROUTINE SETNONHAP( NRAWBP )

!**************************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine integrates criteria and toxics pollutant
!      emissions by creating NONHAP values.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!      Created 11/02 by C. Seppanen
!      modified 1/2/06 by B. Baek
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
!***************************************************************************
    USE M3UTILIO

!...........   MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: POLVAL, IPOSCOD, NPCNT, CSOURC, CINTGR

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: INVSTAT, MXIDAT, INVDNAM, INVDVTS,&
    &                    ITMSPC, ITEXPL, ITNAMA, NINVTBL

!...........   This module contains the cross-reference tables
    USE MODXREF, ONLY: LNONHAP, PROC_HAPS

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NPPOL, NCHARS, NEM, NDY

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

!.........  Pollutant names
    CHARACTER(11),      PARAMETER :: NONHAPCOMP  = 'NONHAP_TYPE'
    CHARACTER(6),       PARAMETER :: NONHAPREFIX = 'NONHAP'
    CHARACTER(4),       PARAMETER :: NOIEND = '_NOI'

!.........   SUBROUTINE ARGUMENTS
    INTEGER , INTENT (INOUT)      :: NRAWBP  ! no. raw records by pollutant

!.........   Local allocatable arrays
    INTEGER, ALLOCATABLE :: VOCPOS( : )      ! location of VOC or TOG entry in srcs array
    INTEGER, ALLOCATABLE :: NVOCPOS( : )     ! location of VOC or TOG entry in srcs array for VOC_INV
    INTEGER, ALLOCATABLE :: VNMPOS( : )      ! position of VOC or TOG in poll name array

    INTEGER, ALLOCATABLE :: TMPIDX( : )      ! sorting index
    INTEGER, ALLOCATABLE :: TMPNPCNT( : )    ! number of pol per source
    INTEGER, ALLOCATABLE :: TMPPOSCOD( : )   ! tmp pollutant code array

    REAL,    ALLOCATABLE :: HAPEANN( : )     ! summed annual VOC/TOG HAPs emissions
    REAL,    ALLOCATABLE :: HAPEDAY( : )     ! summed average day VOC/TOG HAPs emissions

    REAL,    ALLOCATABLE :: TMPPOLVAL( :,: ) ! tmp emissions array

    CHARACTER(NAMLEN3),ALLOCATABLE :: NONHAPNAM( : ) ! NONHAPS pol name array
    CHARACTER(3)      ,ALLOCATABLE :: NONHAPMOD( : ) ! NONHAPS modes array

!.........   Other local variables
    INTEGER  H,I,II,J,K,L,LL,L2,M,N,S   ! counters and indices

    INTEGER  IDXSIZE      ! size of sorting index
    INTEGER  IOS          ! I/O status
    INTEGER  CPOL         ! current pollutant number
    INTEGER  PPOL         ! previous pollutant number
    INTEGER  CPOLRAW      ! pollutant number in raw arrays
    INTEGER  CURRPOS      ! current position in POLVAL and IPOSCOD arrays
    INTEGER  NEWSRCPOS    ! current position in POLVAL and IPOSCOD arrays with VOC_INV poll
    INTEGER  NHAP         ! number of pollutants that needs NONHAP calculation in src array
    INTEGER  NNONHAP      ! number of pollutants that needs NONHAP calculation
    INTEGER  NONVPOS      ! position of NONHAP[VOC|TOG] in poll name array
    INTEGER  NVPOS        ! tmp location of VOC or TOG entry in srcs array for VOC_INV
    INTEGER  VPOS         ! tmp location of VOC or TOG entry in srcs array
    INTEGER  TOXPOS       ! position of NOI toxic in poll name array
    INTEGER  MXWARN       ! maximum number of warnings
    INTEGER  STIDX        ! starting index into pollutant array for current source
    INTEGER  ENDIDX       ! ending index into pollutant array
    INTEGER::NWARN = 0    ! current number of warnings

    LOGICAL  FNDPOL       ! true: found toxic pollutant to be processed
    LOGICAL  FNDVOC       ! true: found VOC entry in inventory
    LOGICAL  MFLAG        ! true: treat all sources integrated
    LOGICAL  LASTFLAG     ! true: entry is last for current source
    LOGICAL  NEEDSORT     ! true: need to resort pollutants for current source
    LOGICAL  NFLAG        ! true: start non-HAP calculation
    LOGICAL::EFLAG = .FALSE. ! true: error occured
    LOGICAL::PROCVOC=.FALSE. ! true: process VOC pollutants
    LOGICAL::PROCTOG=.FALSE. ! true: process TOG pollutants

    CHARACTER(3)       VOC_TOG  ! VOC or TOG
    CHARACTER(14)      NONHAP   ! name of NONHAP[VOC|TOG]
    CHARACTER(NAMLEN3) POLNAM   ! temporary pollutant name
    CHARACTER(NAMLEN3) CURPOL   ! temporary current pollutant name
    CHARACTER(NAMLEN3) VOCPOL   ! temporary current pollutant name
    CHARACTER(300)     BUFFER   ! message buffer
    CHARACTER(256)     MESG     ! message buffer

    CHARACTER(10) :: PROGNAME = 'SETNONHAP' ! program name

!***********************************************************************
!       begin body of subroutine SETNONHAP

!.........  Treat all or partially treat sources as integrate.
    SELECT CASE( PROC_HAPS )
      CASE( 'ALL' )
        ALLOCATE( LNONHAP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LNONHAP', PROGNAME )
        LNONHAP = .TRUE.   ! array
        CINTGR  = 'Y'      ! array
        NFLAG = .TRUE.

      CASE( 'NONE' )
        ALLOCATE( LNONHAP( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'LNONHAP', PROGNAME )
        LNONHAP = .FALSE.   ! array
        CINTGR  = 'N'       ! array
        NFLAG = .TRUE.

      CASE( 'PARTIAL' )   ! use LNONHAP setting from ASGNNHAP.f
        NFLAG = .TRUE.
        WHERE( LNONHAP )
            CINTGR = 'Y'   ! array
        ELSEWHERE
            CINTGR = 'N'   ! array
        ENDWHERE

      CASE DEFAULT
        NFLAG = .FALSE.  ! Don't combine CAP VOC and HAPs.
        CINTGR = ' '

    END SELECT
    IF( .NOT. NFLAG ) RETURN   ! skip combining CAP VOC and HAPs

    CALL ENVSTR( NONHAPCOMP, MESG, 'VOC', VOC_TOG, IOS )
    MESG ='Processing NONHAP'// VOC_TOG // ' calculation.'
    CALL M3MSG2( MESG )

    IF( IOS .NE. 0 ) THEN
        MESG = 'ERROR: NONHAP_TYPE environment variable ' //&
        &       'is not defined for NONHAPVOC calculation'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( VOC_TOG /= 'VOC' .AND. VOC_TOG /= 'TOG' ) THEN
        MESG = 'ERROR: Define NONHAP_TYPE environment variable ' //&
        &       'to either VOC or TOG only.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Set target processing NONHAP pollutant name (NONHAPVOC or NONHAPTOG)
    NONHAP = NONHAPREFIX // VOC_TOG

!.........  Get maximum number of warnings
    MXWARN = ENVINT( WARNSET, ' ', 100, IOS )
    NWARN = 0

!.........  Define size of arrays for NONHAP calculation
    NNONHAP = 0
    DO I = 1, MXIDAT

!.............  Counts total numbers of names for mobile emission processing modes
        IF( INVSTAT( I ) < 0 ) CYCLE

        L  = INDEX( INVDNAM( I ) , ETJOIN )

        IF( L > 0 ) THEN
            LL = LEN( INVDNAM( I ) )

            IF( INVDNAM( I )( L+2:LL ) == VOC_TOG ) THEN
                NNONHAP = NNONHAP + 1
            END IF

!.............  Counts total numbers of pollutatns that need NONHAP calculation for VOC and TOG
        ELSE IF( INVDNAM( I ) == VOC_TOG ) THEN
            NNONHAP = NNONHAP + 1

        END IF

    END DO

!.........  Check that either VOC or TOG is in the inventory
    IF( NNONHAP == 0 ) THEN
        MESG = 'WARNING: No ' // NONHAP // ' calculattion was ' //&
        &       ' processed because there are no ' // VOC_TOG //&
        &       ' & ' // NONHAP // ' found in the INVTABLE file.'//&
        &       CRLF() // BLANK5 // ':: Skipping ' // NONHAP //&
        &      ' calculation and NHAPEXCLUDE processing steps.'
        CALL M3MSG2( MESG )
        CINTGR = ' '    ! reset it to blank for wrinvchr.f
        RETURN
    END IF

!.........  Determine maximum size for allocating memory
!.........  Max number of souces with pol/act including a new species (VOC_INV)
    NRAWBP = NRAWBP + ( NNONHAP * NSRC )

    ALLOCATE(    VNMPOS( NNONHAP ),&
    &          NONHAPNAM( NNONHAP ),&
    &          NONHAPMOD( NNONHAP ),&
    &             VOCPOS( NNONHAP ),&
    &            NVOCPOS( NNONHAP ),&
    &            HAPEANN( NNONHAP ),&
    &            HAPEDAY( NNONHAP ),&
    &           TMPNPCNT( NSRC ),&
    &          TMPPOSCOD( NRAWBP ),&
    &          TMPPOLVAL( NRAWBP,NPPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TMPPOLVAL', PROGNAME )

    VNMPOS = 0
    NONHAPNAM = ' '
    NONHAPMOD = ' '
    TMPNPCNT  = 0
    TMPPOSCOD = 0
    TMPPOLVAL = BADVAL3

!.........  Store the poistion of pollutants for HAPs and NONHAPs
    NHAP = 0
    DO I = 1, MXIDAT

!.............  Skip if pol is activity pol as of VMT
        IF( INVSTAT( I ) < 0 ) CYCLE

        L  = INDEX( INVDNAM( I ) , ETJOIN )

!.............  Store positions of names for mobile emission processing modes
        IF( L > 0 ) THEN
            LL = LEN( INVDNAM( I ) )

            IF( INVDNAM( I )( L+2:LL ) == VOC_TOG ) THEN
                NHAP = NHAP + 1
                VNMPOS( NHAP ) = I
                NONHAPMOD( NHAP ) = INVDNAM( I )( 1:L-1 )
                NONHAPNAM( NHAP ) = INVDNAM( I )( 1:L+1 ) //&
                &               NONHAPREFIX // INVDNAM( I )( L+2:LL )

!.....................  look for NONHAP pol names in a list of pols
                NONVPOS=INDEX1( NONHAPNAM( NHAP ), MXIDAT, INVDNAM )
                IF( NONVPOS < 1 ) THEN
                    MESG = 'ERROR: Pollutant ' // TRIM&
                    &       ( NONHAPNAM( NHAP ))//' is not found '&
                    &       // 'in the INVTABLE file.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

            END IF

!.............  Store position of pollutatns that need NONHAP calculation for VOC and TOG
        ELSE IF( INVDNAM( I ) == VOC_TOG ) THEN
            NHAP = NHAP + 1
            VNMPOS( NHAP ) = I     ! position of VOC or TOG pol in a list of pols
            NONHAPNAM( NHAP ) = NONHAPREFIX // INVDNAM( I )

!.................  Look for NONHAP pol names in a list of pols
            NONVPOS = INDEX1( NONHAPNAM( NHAP ), MXIDAT, INVDNAM )
            IF( NONVPOS < 1 ) THEN
                MESG = 'ERROR: Pollutant ' //&
                &       TRIM( NONHAPNAM( NHAP  ))//' is not found '&
                &       // 'in the INVTABLE file.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

    END DO     ! eloop of a list of inv pols

!.........  Set up logical flags based on found pollutants (VOC or TOG)
    PROCVOC = .TRUE.
    IF( VOC_TOG == 'TOG' ) THEN
        PROCTOG = .TRUE.
        PROCVOC = .FALSE.
    END IF

!.........  Check if any toxics pollutants are to be subtracted
    FNDPOL = .FALSE.
    DO I = 1, MXIDAT
        IF( INVDVTS( I ) /= 'N' ) THEN
            IF( .NOT. PROCTOG .AND. INVDVTS( I ) == 'T' ) CYCLE
            IF( INVSTAT( I ) /= 0 ) THEN
                FNDPOL = .TRUE.
                EXIT
            END IF
        END IF
    END DO

    IF( .NOT. FNDPOL ) THEN
        MESG = 'WARNING: No toxics were selected from ' //&
        &       'the inventory as a subtratacter for ' //&
        &       'the NHAPEXCLUDE processing.'
        CALL M3MSG2( MESG )
        CINTGR = ' '    ! reset it to blank for wrinvchr.f
        RETURN
    END IF

    CURRPOS = 0
    NEWSRCPOS = 0

!.........  Loop through sources
    DO I = 1, NSRC

!.............  Initialize values for this source
        VOCPOS = 0
        NVOCPOS = 0
        HAPEANN = 0.
        HAPEDAY = 0.

        FNDVOC = .FALSE.
        FNDPOL = .FALSE.
        LASTFLAG = .FALSE.

!.............  store org number of pollunt to tmp array
        TMPNPCNT( I ) = NPCNT( I )

!.............  Process source if it is not integrated
        IF( ALLOCATED( LNONHAP ) ) THEN
            IF( .NOT. LNONHAP( I ) ) THEN

!.................  Loop through pollutants for this source
                DO J = 1, NPCNT( I )


!.....................  Increment current position in arrays
                    CURRPOS = CURRPOS + 1

!.....................  Store pollutant for this position
                    CPOL = IPOSCOD( CURRPOS )

!.....................  Increment current new position in arrays
                    NEWSRCPOS = NEWSRCPOS + 1

!.....................  Store to tmp arrays of IPOSCOD and POLVAL
                    TMPPOSCOD( NEWSRCPOS )   = IPOSCOD( CURRPOS )
                    TMPPOLVAL( NEWSRCPOS,: ) = POLVAL ( CURRPOS,: )

!.....................  If current pollutant is VOC(_INV) entry, store position and emissions
                    DO H = 1, NNONHAP

                        IF( CPOL == VNMPOS( H ) ) THEN

!.............................  Check to see a new pol VOC_INV listed in INVTABLE file
                            VOCPOL = TRIM( INVDNAM( CPOL ) ) // '_INV'
                            NONVPOS = INDEX1( VOCPOL, MXIDAT, INVDNAM )

                            IF( NONVPOS < 1 ) THEN
                                MESG = 'ERROR: Pollutant '//TRIM(VOCPOL)&
                                & //' is not found in the INVTABLE file.'
                                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                            END IF

                        END IF

                    END DO

!.....................  If pollutant is not part of VOC or TOG, cycle
                    IF( INVDVTS( CPOL ) == 'N' ) CYCLE

!.....................  Find pollutant position in raw list
                    CPOLRAW = INDEX1( INVDNAM( CPOL ), NINVTBL, ITNAMA )

!.....................  If pollutant is not a model species, set it to zero
                    IF( .NOT. ITMSPC( CPOLRAW ) ) THEN
                        POLVAL( CURRPOS,NEM ) = 0.0
                        POLVAL( CURRPOS,NDY ) = 0.0

!..................... Otherwise, if pollutant is not an explicit species, rename to NOI
                    ELSE IF( .NOT. ITEXPL( CPOLRAW ) ) THEN

!.........................  Create NOI name
                        POLNAM = INVDNAM( CPOL )
                        IF( LEN_TRIM( POLNAM ) > 11 ) THEN
                            POLNAM = POLNAM(1:12)
                        END IF
                        POLNAM = TRIM( POLNAM ) // NOIEND

!.........................  Find NOI name in pollutant names array
                        TOXPOS = INDEX1( POLNAM, MXIDAT, INVDNAM )

!.........................  If found, set pollutant for current source to NOI pollutant
                        IF( TOXPOS > 0 ) THEN
                            IPOSCOD( CURRPOS ) = TOXPOS
                            INVSTAT( TOXPOS ) = 2
                        ELSE
                            MESG = 'ERROR: Non-integrated toxic ' //&
                            &       'pollutant ' // TRIM( POLNAM ) //&
                            &       ' was not found in the ' //&
                            &       'INVTABLE file.'
                            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                        END IF

                    END IF  ! check if pollutant is model or explicit species

!.....................  Store to tmp arrays of IPOSCOD and POLVAL due to the name change to _NOI
                    TMPPOSCOD( NEWSRCPOS )   = IPOSCOD( CURRPOS )
                    TMPPOLVAL( NEWSRCPOS,: ) = POLVAL ( CURRPOS,: )

                END DO  ! loop through pollutants

!.................  Skip rest of loop since we're done with this source
                CYCLE

            END IF
        END IF

!.............  Loop through pollutants for this source
        DO J = 1, NPCNT( I )

!.................  Increment current position in arrays
            CURRPOS   = CURRPOS + 1

!.................  Store pollutant for this position
            CPOL   = IPOSCOD( CURRPOS )
            POLNAM = INVDNAM( CPOL )

            IF( VOC_TOG == TRIM( POLNAM ) ) FNDVOC = .TRUE.

!.....................  Increment current new position in arrays
            NEWSRCPOS = NEWSRCPOS + 1

!.................  Store to tmp arrays of IPOSCOD and POLVAL
            TMPPOSCOD( NEWSRCPOS )   = IPOSCOD( CURRPOS )
            TMPPOLVAL( NEWSRCPOS,: ) = POLVAL ( CURRPOS,: )

!.................  Check VOC or TOG existence
            L = INDEX( POLNAM, ETJOIN )
            IF( L > 0 ) THEN
                LL = LEN( POLNAM )
                IF( VOC_TOG == TRIM( POLNAM( L+2:LL ) ) ) FNDVOC = .TRUE.
            ELSE
                IF( VOC_TOG == TRIM( POLNAM ) ) FNDVOC = .TRUE.
            END IF

!.................  If current pollutant is VOC(_INV) entry, save position and emissions
            DO H = 1, NNONHAP

                IF( CPOL == VNMPOS( H ) ) THEN

!.........................  Store original and new locations
                    VOCPOS ( H ) = CURRPOS     ! store org VOC location
                    NVOCPOS( H ) = NEWSRCPOS   ! store new VOC location

!.........................  Check to see a new pol VOC_INV listed in INVTABLE file
                    VOCPOL = TRIM( POLNAM ) // '_INV'
                    NONVPOS = INDEX1( VOCPOL, MXIDAT, INVDNAM )

                    IF( NONVPOS < 1 ) THEN
                        MESG = 'ERROR: Pollutant '//TRIM( VOCPOL )&
                        &     //' is not found in the INVTABLE file.'
                        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                    END IF

                END IF

            END DO

!.................  Check if this is last pollutant for this source
            IF( J == NPCNT( I ) ) THEN
                LASTFLAG = .TRUE.

!.................  Otherwise, if not part of VOC or TOG, skip rest of loop
            ELSE IF( INVDVTS( CPOL ) == 'N' ) THEN
                CYCLE

            END IF

!.................  Check emission processing modes
            IF( L > 0 ) THEN
                M = INDEX1( POLNAM( 1:L-1 ), NNONHAP, NONHAPMOD )
            ELSE
                M = INDEX1( ' ', NNONHAP, NONHAPMOD )
            END IF


!.................  Sum toxic emissions for this source
!                   INVDVTS = 'V' => part of VOC and TOG
!                   INVDVTS = 'T' => part of TOG only
            IF( PROCVOC .AND. INVDVTS( CPOL ) == 'V' ) THEN
                HAPEANN( M ) = HAPEANN( M ) + POLVAL( CURRPOS,NEM )
                HAPEDAY( M ) = HAPEDAY( M ) + POLVAL( CURRPOS,NDY )
                FNDPOL = .TRUE.
            END IF

            IF( PROCTOG ) THEN
                HAPEANN( M ) = HAPEANN( M ) + POLVAL( CURRPOS,NEM )
                HAPEDAY( M ) = HAPEDAY( M ) + POLVAL( CURRPOS,NDY )
                FNDPOL = .TRUE.

            END IF

!.................  If this is not the last entry for source, cycle
!                   Otherwise, check conditions and subtract toxic
!                   emissions from criteria values
            IF( .NOT. LASTFLAG ) CYCLE

!.................  Format information for this source
            CALL FMTCSRC( CSOURC( I ), NCHARS, BUFFER, L2 )

!.................  Check the status of VOC and HAPs prior to NONHAPVOC calculation
            IF( .NOT. FNDVOC .AND. .NOT. FNDPOL ) THEN
                MESG = 'WARNING: Both ' // VOC_TOG //&
                &    ' and toxics are NOT found for the source:' //&
                &    CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )
            END IF

            IF( FNDVOC .AND. .NOT. FNDPOL ) THEN
                MESG = 'ERROR: Found ' // VOC_TOG //&
                &    ' but no toxics found for the source:' //&
                &    CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF

            IF( .NOT. FNDVOC .AND. FNDPOL ) THEN
                MESG = 'ERROR: Found toxics but no ' // VOC_TOG //&
                &    ' found for the source:' //&
                &    CRLF() // BLANK5 // BUFFER( 1:L2 )
                CALL M3MESG( MESG )
                EFLAG = .TRUE.
            END IF

!.................  Loop through pollutants that need NONHAP calculation
            DO K = 1, NNONHAP

!.....................  Restore position of VOC or TOG pollutants
                VPOS   = VOCPOS( K )
                NVPOS  = NVOCPOS( K )
                CURPOL = INVDNAM( VNMPOS( K ) )   ! VOC, EXH__VOC,,,

!....................  Check if this source had no criteria and no toxics
!                      This happens for activity only sources
                IF( VPOS == 0 .AND. FNDPOL ) CYCLE
                IF( VPOS == 0 .AND. .NOT. FNDPOL ) CYCLE
                IF( VPOS /= 0 .AND. .NOT. FNDPOL ) CYCLE

!.....................  Subtract toxic emissions from criteria emissions
                POLVAL( VPOS,NEM ) =&
                &        POLVAL( VPOS,NEM ) - HAPEANN( K )
                POLVAL( VPOS,NDY ) =&
                &        POLVAL( VPOS,NDY ) - HAPEDAY( K )

!.....................  Check that annual NONHAP value is not negative
                IF( POLVAL( VPOS,NEM ) < 0. ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Total annual toxic ' //&
                        &   'emissions greater than annual ' //&
                        &   TRIM( CURPOL )// ' emissions for source:'&
                        &   // CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

                    POLVAL( VPOS,NEM ) = .0
                END IF

!....................  Check that average day NONHAP value is not negative
                IF( POLVAL( VPOS,NDY ) < 0. ) THEN
                    IF( NWARN <= MXWARN ) THEN
                        MESG = 'WARNING: Total average day ' //&
                        &       'toxic emissions greater than ' //&
                        &       'average day ' // TRIM( CURPOL ) //&
                        &       ' emissions for source:' // CRLF() //&
                        &       BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        NWARN = NWARN + 1
                    END IF

                    POLVAL( VPOS,NDY ) = .0
                END IF

!.....................  Rename VOC to NONHAPVOC
!.....................  increment a number of poll for NONHAP[VOC|TOG]
                NONVPOS = INDEX1( NONHAPNAM( K ), MXIDAT, INVDNAM )
                IPOSCOD  ( VPOS )   = NONVPOS
                TMPPOSCOD( NVPOS )  = NONVPOS
                TMPPOLVAL( NVPOS,: )= POLVAL( VPOS,: )
                INVSTAT  ( NONVPOS )= 2

            END DO  !loop through pollutants that need NONHAP calculation

        END DO  ! loop through pollutants

    END DO  ! loop through sources

!........   Error message
    IF( EFLAG ) THEN
        MESG = 'ERROR: There is a problem during combining VOC and HAPs'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    NRAWBP = NEWSRCPOS

!.........  Update NPCNT, IPOSCOD and POLVAL including VOC_INV to retain org VOC
    DEALLOCATE( IPOSCOD, POLVAL )

    ALLOCATE( POLVAL( NRAWBP,NPPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'POLVAL', PROGNAME )
    ALLOCATE( IPOSCOD( NRAWBP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IPOSCOD', PROGNAME )

    POLVAL = 0.0
    IPOSCOD = 0

!.........  Reinitialize arrays with updated arrays including VOC_INV
    CURRPOS = 0
    DO S = 1, NSRC

        NPCNT( S ) = TMPNPCNT( S )

        DO J = 1, NPCNT( S )
!................  Increment current position in arrays
            CURRPOS   = CURRPOS + 1

!................  Store pollutant for this position
            IPOSCOD( CURRPOS ) = TMPPOSCOD( CURRPOS )
            POLVAL( CURRPOS,:) = TMPPOLVAL( CURRPOS,: )
            POLNAM = INVDNAM( IPOSCOD( CURRPOS) )
        ENDDO

    ENDDO

    DEALLOCATE( TMPNPCNT, TMPPOLVAL, TMPPOSCOD )

!.........  Sort POLVAL and IPOSCOD to put new pollutants (NONHAP and NOI) in correct order

!.........  Determine maximum size for sorting index, then allocate memory
    IDXSIZE = MAXVAL( NPCNT )

    ALLOCATE( TMPIDX( IDXSIZE ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TMPIDX', PROGNAME )
    ALLOCATE( TMPPOLVAL( IDXSIZE,NPPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TMPPOLVAL', PROGNAME )
    ALLOCATE( TMPPOSCOD( IDXSIZE ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TMPPOSCOD', PROGNAME )

    STIDX = 1

!.........  Loop through sources
    DO I = 1, NSRC

!.............  Set ending index for this source
        ENDIDX = STIDX + NPCNT( I ) - 1

!.............  Initialize sorting index and check if any sorting needs to be done
        PPOL = 0
        NEEDSORT = .FALSE.

        DO J = 1, NPCNT( I )
            CPOL = IPOSCOD( STIDX + J - 1 )

!.................  If current pollutant is lower than previous, need to resort
            IF( CPOL < PPOL ) THEN
                NEEDSORT = .TRUE.
            END IF

            PPOL = CPOL
            TMPIDX( J ) = J
        END DO

!.............  Make sure this source needs to be resorted
        IF( NEEDSORT ) THEN

!.................  Store current values in temporary arrays
            TMPPOLVAL( 1:NPCNT( I ),: ) = POLVAL ( STIDX:ENDIDX,: )
            TMPPOSCOD( 1:NPCNT( I ) )   = IPOSCOD( STIDX:ENDIDX )

!.................  Sort section of IPOSCOD corresponding to this source
            CALL SORTI1( NPCNT( I ), TMPIDX( 1:NPCNT( I ) ),&
            &             IPOSCOD( STIDX:ENDIDX ) )

!.................  Loop through pollutants for current source
            DO J = 1, NPCNT( I )

                K = TMPIDX( J )

                POLVAL ( STIDX + J - 1,: ) = TMPPOLVAL( K,: )
                IPOSCOD( STIDX + J - 1 )   = TMPPOSCOD( K )

            END DO

        END IF  ! check if pollutants need sorting

!.............  Increment starting index
        STIDX = ENDIDX + 1

    END DO

!.........  Deallocate local memory
    DEALLOCATE( TMPIDX, TMPPOLVAL, TMPPOSCOD )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

93000 FORMAT( A )
94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE SETNONHAP

