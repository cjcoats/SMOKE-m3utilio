
SUBROUTINE XREFTBL( OPTYPE, NXREF )

    !***********************************************************************
    !  subroutine body starts at line  71
    !
    !  DESCRIPTION:
    !      This subroutine will do the logic for grouping the cross references for
    !      any source category.  In the logic phase, it will merely store the group
    !      number and the count in the group for each x-ref entry. Then, depending on
    !      the type of operation (e.g., temporal, speciation, reactivity packet),
    !      it will allocate memory for an populate the appropriate grouped
    !      tables.  The unsorted x-ref data and grouped tables will be passed by a
    !      module for cross-reference data.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !     Created 1/99 by M. Houyoux
    !     Version 11/2023 by CJC:  USE M3UTILIO; conversion to ".f90"
    !
    !*************************************************************************
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

    !.........  MODULES for public variables
    !.........  This module is for cross reference tables
    USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, ISPTA, CMACTA, CISICA,   &
                       TXCNT, NXTYPES, XDUPCHK

    !.........  This module contains the information about the source category
    USE MODINFO, ONLY: MXCHRS, NIPPA, SCCLEV1, SCCLEV2,                 &
                       SCCLEV3, LSCCEND, SC_ENDP, SC_BEGP,              &
                       JSCC, RSCCBEG, PLTIDX, EANAM, NIACT, NCHARS

    !.........  This module contains the speciation profiles
    USE MODSPRO, ONLY: CMBDEX

    IMPLICIT NONE

    !...........   INCLUDES

    INCLUDE 'EMCNST3.h90'   !  emissions constant parameters

    !...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: OPTYPE ! operation type (tmprl,spec,ctg...)
    INTEGER     , INTENT (IN) :: NXREF  ! no. ungrouped x-ref entries

    !.........  EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: SETSCCTYPE, CHKEXPSCC, CHKEXPSIC, USEEXPGEO

    !...........   Local parameters
    INTEGER, PARAMETER :: SNFLEN3 = SRCLEN3 - FIPLEN3

    !...........  Local allocatable arrays...

    !...........   Arrays for intial pass through x-ref to determine degree of
    !              each record and store the type and count of that type
    INTEGER     XTYPE( NXREF )      ! group number of x-ref entry
    INTEGER     XTCNT( NXREF )      ! position in group of x-ref entry
    LOGICAL     DEFAULT( NIPPA+1 )  ! true: if default entry in x-ref

    !...........  Arrays for counting number of x-ref records in each matching
    !             degree and comparing current record with previous record as part
    !             of the counting. Note that although all arrays below are
    !             available for each degree of matching, only the array elements
    !             that are appropriate for a given degree are actually populated.
    INTEGER            N     ( 0:NXTYPES ) ! cnt for degree of matching

    CHARACTER(FIPLEN3) PCFIP ( NXTYPES )  ! previous co/st/cy code
    CHARACTER(SCCLEN3) PTSCC ( NXTYPES )  ! previous SCC
    CHARACTER(SICLEN3) PCSIC ( NXTYPES )  ! previous SIC
    CHARACTER(SRCLEN3) PCSRC ( NXTYPES )  ! previous CSRC
    CHARACTER(SS5LEN3) PCSSC ( NXTYPES )  ! previous CSRC(part) // SCC
    CHARACTER(MACLEN3) PCMCT ( NXTYPES )  ! previous MACT

    !...........   Array of source characeristics
    CHARACTER(300)          CHARS( MXCHRS )

    !...........   Other local variables
    INTEGER       I, J, J1, J2, K, L  ! counter and indices

    INTEGER       IDUM             ! dummy integer
    INTEGER       IOS              ! i/o status
    INTEGER       ISP              ! temporary pollutant position in EANAM
    INTEGER       LOPT             ! length of OPTYPE
    INTEGER       NCHKCHR          ! position of last non-SCC src char
    INTEGER       NT               ! code for specificity of x-ref entry
    INTEGER    :: PISP = IMISS3    ! previous iteration ISP

    LOGICAL    :: CFLAG = .FALSE.  ! true: operation type is control cntls
    LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred
    LOGICAL    :: FFLAG = .FALSE.  ! true: operation type is emis. factors
    LOGICAL    :: GFLAG = .FALSE.  ! true: operation type is ctg cntls
    LOGICAL    :: IFLAG = .FALSE.  ! true: operation type is gridding
    LOGICAL    :: JFLAG = .FALSE.  ! true: operation type is projection
    LOGICAL    :: LFLAG = .FALSE.  ! true: operation type is allowable cntls
    LOGICAL    :: MFLAG = .FALSE.  ! true: operation type is VMT mix
    LOGICAL    :: NFLAG = .TRUE.   ! true: has pol or activity-specific
    LOGICAL    :: OFLAG = .FALSE.  ! true: use extra SCC matches
    LOGICAL    :: PFLAG = .FALSE.  ! true: operation type is speeds
    LOGICAL    :: POADFLT          ! true: okay to have pol/act-spec dfaults
    LOGICAL    :: RFLAG = .FALSE.  ! true: operation type is reactivty cntls
    LOGICAL    :: SAMEFLAG = .FALSE. ! true: same x-ref chars as previous iter
    LOGICAL       SCCFLAG          ! true: SCC type is different from previous
    LOGICAL    :: SFLAG = .FALSE.  ! true: operation type is speciation
    LOGICAL    :: TFLAG = .FALSE.  ! true: operation type is temporal
    LOGICAL    :: TAGFLAG = .FALSE. ! true: operation type is tagging
    LOGICAL    :: XFLAG = .FALSE.  ! true: operation type is nonhapVOC exclusion
    LOGICAL    :: YFLAG = .FALSE.  ! true: operation type is area-to-point
    LOGICAL, SAVE :: FIRSTIME = .TRUE.  ! true: first time subroutine called
    LOGICAL, SAVE :: FULLSCC  = .FALSE. ! true: use only full SCC entries

    CHARACTER          CDUM          ! dummy character string
    CHARACTER(256)     BUFFER        ! source definition buffer
    CHARACTER(256)     MESG          ! message buffer

    CHARACTER(SCCLEN3) SCCL          ! left digits of TSCC
    CHARACTER(SCCLEN3) SCCR          ! RHS digits of TSCC (old method)
    CHARACTER(SCCLEN3) SCCR_A        ! RHS for level 1 of SCC
    CHARACTER(SCCLEN3) SCCR_B        ! RHS for level 2 of SCC
    CHARACTER(SCCLEN3) SCCR_C        ! RHS for level 3 of SCC
    CHARACTER(SCCLEN3) SCRZERO       ! buf for 0 right 5 digits of TSCC
    CHARACTER(SNFLEN3) CNFIP         ! characterstics without FIPS code
    CHARACTER(SRCLEN3) CSRC          ! temporary source char string
    CHARACTER(FPLLEN3) CFPL         ! temporary FIPS // plant ID
    CHARACTER(FIPLEN3) CFIP          ! temporary (character) FIPS code
    CHARACTER(MACLEN3) CMCT          ! temporary MACT code
    CHARACTER(FIPLEN3) FIPZERO       ! buffer for zero FIPS code
    CHARACTER(MACLEN3) MCTZERO       ! buffer for zero MACT code
    CHARACTER(SCCLEN3) PSCC          ! previous SCC
    CHARACTER(SCCLEN3) TSCC          ! temporary SCC
    CHARACTER(SCCLEN3) SCCZERO       ! buffer for zero SCC
    CHARACTER(SCCLEN3) SCCZ_A        ! buffer for level 1 zero string
    CHARACTER(SCCLEN3) SCCZ_B        ! buffer for level 2 zero string
    CHARACTER(SCCLEN3) SCCZ_C        ! buffer for level 3 zero string
    CHARACTER(SS5LEN3) CSRCSCC       ! buffer for source // SCC
    CHARACTER(SICLEN3) CSIC          ! buffer for SIC
    CHARACTER(SICLEN3) CSICL         ! buffer for left 2-digit SIC
    CHARACTER(SICLEN3) CSICR         ! buffer for right SIC
    CHARACTER(SICLEN3) SICZERO       ! buffer for SIC zero
    CHARACTER(SICLEN3) SICRZERO      ! buffer for right SIC zero

    CHARACTER(16), PARAMETER :: PNAME = 'XREFTBL' ! program name

    !***********************************************************************
    !   begin body of subroutine XREFTBL

    !.........  For first time routine is called ...
    IF( FIRSTIME ) THEN

        MESG = 'Use only full SCC matches'
        FULLSCC = ENVYN ( 'FULLSCC_ONLY', MESG, .FALSE., I )
        IF ( I .GT. 0 ) THEN
            CALL M3EXIT( PNAME,0,0, 'Bad env vble "FULLSCC_ONLY"', 2 )
        END IF
        FIRSTIME = .FALSE.

    ENDIF

    LOPT = LEN_TRIM( OPTYPE )

    !.........  Check for valid operation type
    SELECT CASE( OPTYPE )

      CASE( 'ALLOWABLE' )
        POADFLT = .TRUE.
        LFLAG   = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'AR2PT' )
        POADFLT = .FALSE.
        YFLAG   = .TRUE.
        NFLAG   = .FALSE.
      CASE( 'CONTROL' )
        POADFLT = .TRUE.
        CFLAG   = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'CTG' )
        POADFLT = .TRUE.
        GFLAG   = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'EMISFACS' )
        POADFLT = .TRUE.
        FFLAG   = .TRUE.
      CASE( 'GRIDDING' )
        POADFLT = .FALSE.
        IFLAG   = .TRUE.
        NFLAG   = .FALSE.
      CASE( 'MACT' )
        POADFLT = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'NONHAP' )
        POADFLT = .TRUE.
        XFLAG   = .TRUE.
        NFLAG   = .FALSE.
        OFLAG   = .TRUE.
        PLTIDX  = 2     ! for point (added for source-specific NHAPEXCLUDE)
      CASE( 'PROJECTION' )
        POADFLT = .TRUE.
        JFLAG   = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'REACTIVITY' )
        POADFLT = .FALSE.
        RFLAG   = .TRUE.
        OFLAG   = .TRUE.
      CASE( 'SPEED' )
        POADFLT = .FALSE.
        PFLAG   = .TRUE.
        NFLAG   = .FALSE.
      CASE( 'SPECIATION' )
        POADFLT = .TRUE.
        SFLAG   = .TRUE.
      CASE( 'TAGGING' )
        POADFLT = .FALSE.
        TAGFLAG = .TRUE.
      CASE( 'TEMPORAL' )
        POADFLT = .TRUE.
        TFLAG   = .TRUE.
      CASE( 'VMTMIX' )
        POADFLT = .FALSE.
        MFLAG   = .TRUE.
        NFLAG   = .FALSE.

      CASE DEFAULT

        MESG = 'INTERNAL ERROR: Operation type "' // OPTYPE( 1:LOPT ) //        &
               '" not known in subroutine ' // PNAME
        CALL M3MSG2( MESG )
        CALL M3EXIT( PNAME, 0, 0, ' ', 2 )

    END SELECT

    !.........  Set up zero strings for FIPS code of zero and SCC code of zero
    FIPZERO  = REPEAT( '0', FIPLEN3 )
    SCCZERO  = REPEAT( '0', SCCLEN3 )
    SCRZERO  = REPEAT( '0', SCCLEN3 - LSCCEND )
    SCCZ_A   = REPEAT( '0', SCCLEN3 - SCCLEV1 )
    SCCZ_B   = REPEAT( '0', SCCLEN3 - SCCLEV2 )
    SCCZ_C   = REPEAT( '0', SCCLEN3 - SCCLEV3 )
    SICZERO  = REPEAT( '0', SICLEN3 )
    SICRZERO = REPEAT( '0', SICLEN3 - 2 )
    MCTZERO  = REPEAT( '0', MACLEN3 )

    !.........  Initialize default array
    DEFAULT = .FALSE.   ! array

    !.........  Initialize arrays for counting number of x-ref records in each
    !           degree of matching
    N      = 0    ! arrays
    PCFIP  = ' '
    PCSRC  = ' '
    PCSSC  = ' '
    PTSCC  = ' '

    !.........  Initialize source characteristics
    CHARS = ' '  ! array

    !.........  Initialize index check
    NCHKCHR = NCHARS
    IF( JSCC .GT. 0 ) NCHKCHR = NCHARS - 1

    !.........  Loop through and count entries of each type. Store type.
    !.........  For CSRC, don't include pollutant for grouping.
    ISP    = 0
    PSCC   = ' '
    DO I = 1, NXREF

        J = INDXTA( I )
        CSRC    = CSRCTA( J )( 1:SC_ENDP( NCHARS ) )
        TSCC    = CSCCTA( J )

        IF( NFLAG ) ISP = ISPTA ( J )  ! no pollutants for gridding

        IF( ALLOCATED( CMACTA ) ) THEN
            CMCT = CMACTA( J )
        ELSE
            CMCT = MCTZERO
        END IF

        IF( ALLOCATED( CISICA ) ) THEN
            CSIC = CISICA( J )
        ELSE
            CSIC = SICZERO
        END IF

        DO J = 1, NCHARS
            CHARS( J ) = CSRC( SC_BEGP( J ):SC_ENDP( J ) )
        END DO

        !.............  Rearrange CHARS if SCC is a part of the source definition
        !               because we have now stored SCC separately.  It will be
        !               much easier to group the source characteristics this way.
        !               Still go back to using original definition when storing
        !               in tables (i.e., use CSRC)
        IF( JSCC .GT. 0 ) THEN
            DO J = JSCC, NCHARS - 1
                CHARS( J ) = CHARS( J + 1 )
            END DO
            CHARS( NCHARS ) = ' '
        END IF

        !.............  Set up partial strings for checking country/state/county
        CFIP    = CHARS( 1 )
        CFPL   = CSRC( SC_BEGP( 1 ):SC_ENDP( PLTIDX ) )
        IF( PLTIDX /= 0 .AND. PLTIDX <= NCHARS ) THEN
            CNFIP = CSRC( SC_BEGP( PLTIDX ):SC_ENDP( NCHARS ) )
        ELSE
            CNFIP = ' '
        END IF

        IF( CSIC /= SICZERO ) THEN      !!  If SIC given, setup SIC fields
            TSCC  = SCCZERO
            CSICL = CSIC( 1:SICLEN3-2 )
            CSICR = CSIC( SICLEN3-3:SICLEN3 )

        ELSE            !!If SIC *not* included, setup SCC fields
        
            CSICR = SICRZERO

            !.................  Set type of SCC
            SCCFLAG = SETSCCTYPE( TSCC )

            !.................  If SCC type has changed, reset zero strings
            IF( SCCFLAG ) THEN
                SCRZERO  = REPEAT( '0', SCCLEN3 - LSCCEND )
                SCCZ_A   = REPEAT( '0', SCCLEN3 - SCCLEV1 )
                SCCZ_B   = REPEAT( '0', SCCLEN3 - SCCLEV2 )
                SCCZ_C   = REPEAT( '0', SCCLEN3 - SCCLEV3 )
            END IF

            !.................  Standard strings for SCC left and right matching
            SCCL    = TSCC(       1:LSCCEND )
            SCCR    = TSCC( RSCCBEG:SCCLEN3 )

            !.................  More partial strings for special SCC levels matching
            SCCR_A = TSCC( SCCLEV1 + 1:SCCLEN3 )
            SCCR_B = TSCC( SCCLEV2 + 1:SCCLEN3 )
            SCCR_C = TSCC( SCCLEV3 + 1:SCCLEN3 )

        END IF

        !.............  Reset flag for identifying when a record with the same
        !               info but different pollutant is encountered
        SAMEFLAG = .FALSE.

        !.............  Select cases
        !.............  Note that since these are sorted in order of increasing FIPS
        !               code, SCC, pollutant index, etc., that the entries with zero for
        !               these characteristics will appear earlier in the sorted list
        IF( CFIP .EQ. FIPZERO ) THEN                           ! FIPS code is default

            IF( CMCT .NE. MCTZERO ) THEN                       ! have valid MACT code

                IF( TSCC .EQ. SCCZERO ) THEN                      ! SCC is default

                    NT = 32
                    IF( CMCT .NE. PCMCT( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCMCT( NT ) = CMCT

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                                             ! Complete SCC

                    NT = 33
                    IF( CMCT .NE. PCMCT( NT ) .OR.          &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCMCT( NT ) = CMCT
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

        !.....................  Set SCC to zero to avoid lower level SCC checks
                TSCC = REPEAT( '0', SCCLEN3 )

            ELSE IF( CSICR .NE. SICRZERO .OR.                 &    ! Full SIC defined
                      ( CSIC .NE. SICZERO .AND. CHKEXPSIC( CSIC ) ) ) THEN

                NT = 27
                IF( CSIC .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCSIC( NT ) = CSIC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( CSIC .NE. SICZERO ) THEN                ! Left SIC defined

                NT = 26
                IF( CSICL .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCSIC( NT ) = CSICL

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( TSCC .EQ. SCCZERO ) THEN                  ! SCC is default

                IF( POADFLT ) THEN                          ! Pollutant-specific permitted

                    IF( ISP .EQ. 0 ) THEN
                        NT = 1
                        N( NT ) = 1

                    ELSEIF( ISP .NE. 0 .AND. .NOT.DEFAULT( ISP ) ) THEN

                        DEFAULT( ISP ) = .TRUE.
                        NT = 1
                        N( NT ) = 1

                    ELSE                                    ! Report and skip
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                         ! Pollutant-specific not permitted

                    IF( ISP .EQ. 0 .AND.                &    ! Pollutant not specified
                        .NOT.DEFAULT( 1 ) ) THEN

                        DEFAULT( 1 ) = .TRUE.
                        NT = 1
                        N( NT ) = 1

                    ELSE IF( ISP .NE. 0 ) THEN                ! Report and skip
                        MESG = 'Cannot use pollutant-specific ultimate-default'
                        CALL REPORT_INVALID_XREF( MESG )
                        NT = 0

                    ELSE IF( DEFAULT( 1 ) ) THEN              ! Report and skip
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

            ELSE IF( .NOT. FULLSCC .AND. SCCR .EQ. SCRZERO .AND.&
                     .NOT. TAGFLAG .AND. .NOT. CHKEXPSCC( TSCC ) ) THEN            ! Left SCC

                NT = 2
                IF( TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            ELSE                                             ! Complete SCC

                NT = 3
                IF( TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            END IF

        !.................  Section for special SCC levels for controls. This is not an
        !                   efficient way to implement this, but it's needed so long
        !                   as the old Right-left method is still needed.
            IF ( .NOT. FULLSCC .AND. OFLAG .AND. NT .NE. 1 .AND.&
                  TSCC .NE. SCCZERO .AND. .NOT. CHKEXPSCC( TSCC ) )THEN

                IF( SCCR_A .EQ. SCCZ_A ) THEN             !  Level-1 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 17
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN        !  Level-2 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 18
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN        !  Level-3 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 19
                    IF( TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PTSCC( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

            END IF

        ELSEIF( .NOT. USEEXPGEO() .AND.                         &
                 CFIP( STALEN3+1:FIPLEN3 ) == '000' ) THEN      ! County code is default

            IF( CMCT .NE. MCTZERO ) THEN                        ! have valid MACT code

                IF( TSCC .EQ. SCCZERO ) THEN                    ! SCC is default

                    NT = 34
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        CMCT .NE. PCMCT( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP( NT ) = CFIP
                        PCMCT( NT ) = CMCT

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                                             ! Complete SCC

                    NT = 35
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) .OR.              &
                        CMCT .NE. PCMCT( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP( NT ) = CFIP
                        PTSCC( NT ) = TSCC
                        PCMCT( NT ) = CMCT

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

        !.....................  Set SCC to zero to avoid lower level SCC checks
                TSCC = REPEAT( '0', SCCLEN3 )

            ELSE IF( CSICR .NE. SICRZERO .OR.                   &    ! Full SIC defined
                     ( CSIC .NE. SICZERO .AND. CHKEXPSIC( CSIC ) ) ) THEN

                NT = 29
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    CSIC .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP( NT ) = CFIP
                    PCSIC( NT ) = CSIC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( CSIC .NE. SICZERO ) THEN                ! Left SIC defined

                NT = 28
                IF( CFIP  .NE. PCFIP( NT ) .OR.             &
                    CSICL .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP( NT ) = CFIP
                    PCSIC( NT ) = CSICL

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( TSCC .EQ. SCCZERO ) THEN             ! SCC code is default

                NT = 4
                IF( CFIP .NE. PCFIP( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP( NT ) = CFIP

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSEIF( .NOT. FULLSCC .AND. SCCR .EQ. SCRZERO .AND.         &
                    .NOT. TAGFLAG .AND. .NOT. CHKEXPSCC( TSCC ) ) THEN             ! left SCC

                NT = 5
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0

                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            ELSE                                         ! Complete SCC

                NT = 6
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0

                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            END IF

            !.................  Section for special SCC levels for controls. This is not an
            !                   efficient way to implement this, but it's needed so long
            !                   as the old Right-left method is still needed.

            IF ( .NOT. FULLSCC .AND. OFLAG .AND. NT .NE. 4 .AND.                &
                 TSCC .NE. SCCZERO .AND. .NOT. CHKEXPSCC( TSCC ) ) THEN

                IF( SCCR_A .EQ. SCCZ_A ) THEN             !  State/Level-1 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 20
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN        !  State/Level-2 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 21
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN        !  State/Level-3 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 22
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

            END IF

        ELSEIF( CNFIP .EQ. ' ' .OR.             &
                CNFIP .EQ. TSCC   ) THEN        ! Country/St/Co code is complete

            IF( CMCT .NE. MCTZERO ) THEN                       ! have valid MACT code

                IF( TSCC .EQ. SCCZERO ) THEN                      ! SCC is default

                    NT = 36
                    IF( CFIP .NE. PCFIP( NT ) .OR.          &
                        CMCT .NE. PCMCT( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP( NT ) = CFIP
                        PCMCT( NT ) = CMCT

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE                                             ! Complete SCC

                    NT = 37
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) .OR.              &
                        CMCT .NE. PCMCT( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP( NT ) = CFIP
                        PTSCC( NT ) = TSCC
                        PCMCT( NT ) = CMCT

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

                !.....................  Set SCC to zero to avoid lower level SCC checks
                TSCC = REPEAT( '0', SCCLEN3 )

            ELSEIF( CSICR .NE. SICRZERO .OR.                 &    ! Full SIC defined
                    ( CSIC .NE. SICZERO .AND. CHKEXPSIC( CSIC ) ) ) THEN

                NT = 31
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    CSIC .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PCSIC( NT ) = CSIC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( CSIC .NE. SICZERO ) THEN                ! Left SIC defined

                NT = 30
                IF( CFIP  .NE. PCFIP( NT ) .OR.             &
                    CSICL .NE. PCSIC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PCSIC( NT ) = CSICL

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSE IF( TSCC .EQ. SCCZERO ) THEN             ! SCC code is default

                NT = 7
                IF( CFIP .NE. PCFIP( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP( NT ) = CFIP

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0
                END IF

            ELSEIF( .NOT. FULLSCC .AND. .NOT. CHKEXPSCC( TSCC) .AND.    &
                     SCCR .EQ. SCRZERO ) THEN            ! Left SCC

                NT = 8
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0

                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            ELSE                                             ! Complete SCC

                NT = 9
                IF( CFIP .NE. PCFIP( NT ) .OR.              &
                    TSCC .NE. PTSCC( NT ) ) THEN
                    N( NT ) = N( NT ) + 1
                    PCFIP ( NT ) = CFIP
                    PTSCC( NT ) = TSCC

                ELSEIF( ISP .EQ. PISP ) THEN
                    CALL REPORT_DUP_XREF
                    NT = 0

                ELSE
                    SAMEFLAG = .TRUE.
                END IF

            END IF                                                  ! End SCC

            !.................  Section for special SCC levels for controls. This is not an
            !                   efficient way to implement this, but it's needed so long
            !                   as the old Right-left method is still needed.
            IF ( .NOT. FULLSCC .AND. OFLAG .AND. NT .NE. 7 .AND.        &
                 TSCC .NE. SCCZERO .AND. .NOT. CHKEXPSCC( TSCC ) ) THEN

                IF( SCCR_A .EQ. SCCZ_A ) THEN             !  FIPS/Level-1 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 23
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_B .EQ. SCCZ_B ) THEN        !  FIPS/Level-2 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 24
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                ELSE IF( SCCR_C .EQ. SCCZ_C ) THEN        !  FIPS/Level-3 SCC

                    IF( .NOT. SAMEFLAG )    N( NT ) = N( NT ) - 1
                    NT = 25
                    IF( CFIP .NE. PCFIP( NT ) .OR.              &
                        TSCC .NE. PTSCC( NT ) ) THEN
                        N( NT ) = N( NT ) + 1
                        PCFIP ( NT ) = CFIP
                        PTSCC ( NT ) = TSCC

                    ELSEIF( ISP .EQ. PISP ) THEN
                        CALL REPORT_DUP_XREF
                        NT = 0
                    END IF

                END IF

            END IF

        ELSE                                            ! Plant is specified

            IF( TSCC .EQ. SCCZERO ) THEN               ! SCC code is default

                !.....................  Loop through plant-specific characteristics. Only the
                !                       plant is permitted to not have an SCC not specified
                !.....................  Process NT 16 through 12, and 10
                NT = 9 + MXCHRS
                DO J = MXCHRS, 2, -1

                    IF( CHARS( J ) .NE. ' ' ) THEN

                        IF( CSRC .NE. PCSRC( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PCSRC( NT ) = CSRC
                            EXIT                          ! End loop with NT

                        ELSEIF( ISP .NE. PISP ) THEN
                            EXIT                          ! End loop with NT

                        ELSE
                            CALL REPORT_DUP_XREF
                            NT = 0
                            EXIT                          ! End loop with NT
                        END IF

                    END IF

                    NT = NT - 1

                    !......................  Adjust NT for indexing quirk.  Since NT=11 must have
                    !                        an SCC match, we're really processing for NT=10
                    !                        since in this section, SCC=0
                    IF( NT .EQ. 11 ) NT = 10

                END DO               ! End loop on plant characteristics

                !.....................  Check for Plant-MACT combination
                IF( CMCT .NE. MCTZERO ) THEN                           ! Plant/MACT specified

                    !.........................  Give warning and skip if other fields besides plant are included
                    IF( CFPL .NE. CSRC ) THEN
                        MESG = 'Plant-MACT entry with other source characteristics is skipped'
                        CALL REPORT_INVALID_XREF( MESG )
                        NT = 0

                    ELSE

                        NT = 38
                        IF( CFIP .NE. PCFIP( NT ) .OR.              &
                            CFPL .NE. PCSRC( NT ) .OR.              &
                            CMCT .NE. PCMCT( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PCFIP( NT ) = CFIP
                            PCSRC( NT ) = CFPL
                            PCMCT( NT ) = CMCT

                        ELSEIF( ISP .EQ. PISP ) THEN
                            CALL REPORT_DUP_XREF
                            NT = 0
                        END IF
                    END IF

                END IF

            ELSEIF( .NOT. FULLSCC .AND. .NOT. CHKEXPSCC( TSCC ) .AND.   &
                     SCCR .EQ. SCRZERO ) THEN             ! Left SCC

                MESG = 'Partial SCC "' // TSCC // '" is given instead of full SCC'
                CALL REPORT_INVALID_XREF( MESG )
                NT = 0

            ELSE                                              ! Complete SCC
                !.....................  Loop through plant-specific characteristics,
                !                       and store the most specific entries first.
                !.....................  Only the most specific and plant-only can have
                !                       full TSCC assignment.
                !.....................  Process NT 16 through 11
                NT = 9 + MXCHRS
                DO J = MXCHRS, 2, -1
                    IF( ( J .EQ. NCHKCHR .OR. J .EQ. 2  ) .AND.     &
                        CHARS( J ) .NE. ' ' ) THEN

                        CSRCSCC = CSRC( 1:PTENDL3( J ) ) // TSCC
                        IF( CSRCSCC .NE. PCSSC( NT ) ) THEN
                            N( NT ) = N( NT ) + 1
                            PCSSC( NT ) = CSRCSCC
                            EXIT                          ! End loop with NT

                        ELSEIF( ISP .NE. PISP ) THEN
                            EXIT                          ! End loop with NT

                        ELSE
                            CALL REPORT_DUP_XREF
                            NT = 0
                            EXIT                          ! End loop with NT
                        END IF

                    ELSE IF( CHARS( J ) .NE. ' ' .AND.              &
                            .NOT. XFLAG ) THEN
                        MESG = 'Non-zero SCC is not allowed to be specified'
                        CALL REPORT_INVALID_XREF( MESG )
                        NT = 0
                        EXIT                          ! End loop with NT

                    ELSE IF( CHARS( J ) .NE. ' ' .AND. XFLAG ) THEN

                        CYCLE

                    END IF

                    NT = NT - 1

                ENDDO               ! End loop on plant characteristics

            END IF     ! End SCC

        END IF     ! End degree of Country/State/County code, or plant specified

        PISP = ISP

        XTYPE( I ) = NT
        XTCNT( I ) = N( NT )             ! Dimensioned from 0 so NT can = 0

    END DO                           ! End Loop on sorted x-ref entries

    !.........  If processing tagging data, call routine to process it and
    !           then skip to the end of the program
    IF( TAGFLAG ) THEN

        CALL TAGTABLE( N( 1 ), NXREF, XTYPE, XTCNT )
        GO TO 300

    END IF

    !.........  Allocate the memory for the source-characteristics portion of the
    !           grouped cross-reference tables
    CALL ALOCCHRT( N( 1 ) )

    !.........  Populate the grouped tables of cross-reference source-characteristics
    CALL FILLCHRT( NXREF, XTYPE, XTCNT )

    !.........  Depending on the operation type, first allocate memory for the tables
    !           and initialize profile codes where needed.
    !.........  Then, populate the tables from the sorted and post-processed
    !           cross-reference tables

    !.........  Temporal x-ref tables
    IF( TFLAG ) THEN

        CALL ALOCTTBL( NIPPA, N( 1 ) )
        CALL FILLTTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  Speciation x-ref tables
    ELSE IF( SFLAG ) THEN

        CALL ALOCSTBL( NIPPA, N( 1 ) )
        CALL FILLSTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  Gridding x-ref tables
    ELSE IF( IFLAG ) THEN
        CALL ALOCGTBL( N( 1 ) )
        CALL FILLGTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  Emission factor x-ref tables
    ELSE IF( FFLAG ) THEN
        CALL ALOCETBL( NIACT, N( 1 ) )
        CALL FILLETBL( NIACT, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  Vehicle mix
    ELSE IF( MFLAG ) THEN
        CALL ALOCMTBL( N( 1 ) )
        CALL FILLMTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  Speeds
    ELSE IF( PFLAG ) THEN
        CALL ALOCPTBL( N( 1 ) )
        CALL FILLPTBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  non-HAP inclusion/exclusions file
    ELSE IF( XFLAG ) THEN
    ! Do nothing, because all that is needed is the CHRT* arrays

    !.........  Area-to-point factors assignment
    ELSE IF( YFLAG ) THEN
        CALL ALOCATBL( N( 1 ) )
        CALL FILLATBL( NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    !.........  All control x-ref tables
    ELSE

        CALL ALOCCTBL( NIPPA, N( 1 ) )
        CALL FILLCTBL( NIPPA, NXREF, N( 1 ), XTYPE, XTCNT( 1 ) )

    END IF

300 CONTINUE

    !.........  Store count of records in each group in final variable
    DO I = 1, NXTYPES
        TXCNT( I ) = N( I )
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.............  This internal subprogram writes a warning message for
    !               duplicate entries in the cross-reference file.
    !               Disabled if MODXREF/XDUPCHK is .FALSE.

    SUBROUTINE REPORT_DUP_XREF

        !.............  Local variables
        INTEGER        :: L1
        INTEGER        :: L2
        CHARACTER(300) :: BUFFER
        CHARACTER(300) :: MESG

        !......................................................................

        IF ( .NOT.XDUPCHK ) RETURN

        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )

        MESG = 'WARNING: Duplicate entry in ' // OPTYPE( 1:LOPT ) //    &
               ' x-ref file:' // CRLF() // BLANK10 // BUFFER

        IF( TSCC .NE. SCCZERO ) THEN
            MESG = TRIM( MESG ) // ' TSCC: ' // TSCC
        END IF

        IF( ISP .GT. 0 ) THEN
            MESG = TRIM( MESG ) // ' POL: ' // EANAM( ISP )
        END IF

        IF( CMCT .NE. MCTZERO ) THEN
            MESG = TRIM( MESG ) // ' MACT: ' // CMCT
        END IF

        IF( CSIC .NE. SICZERO ) THEN
            MESG = TRIM( MESG ) // ' SIC: ' // CSIC
        END IF

        CALL M3MSG2( MESG )

    END SUBROUTINE REPORT_DUP_XREF

    !----------------------------------------------------------------------
    !.............  This internal subprogram writes a warning message for
    !               invalid cross-reference entries
    SUBROUTINE REPORT_INVALID_XREF( INMESG )

        !.............  Subprogram arguments
        CHARACTER(*), INTENT( IN ) :: INMESG           ! Input message to include on line 2

        !.............  Local parameters
        CHARACTER(48), PARAMETER :: PART1 = 'WARNING: Skipping invalid cross-reference entry.'

        !.............  Local variables
        INTEGER        :: L1
        INTEGER        :: L2
        CHARACTER(300) :: BUFFER
        CHARACTER(300) :: MESG

        !......................................................................

        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L1 )

        MESG = PART1  //                                        &
               CRLF() // BLANK10 // TRIM( INMESG ) // ':' //    &
               CRLF() // BLANK10 // BUFFER

        IF( CMCT .NE. MCTZERO ) THEN
            MESG = TRIM( MESG ) // ' MACT:' // CMCT
        END IF

        IF( ISP .GT. 0 ) THEN
            MESG = TRIM( MESG ) // ' POA:' // EANAM( ISP )
        END IF

        CALL M3MESG( MESG )

    END SUBROUTINE REPORT_INVALID_XREF

    !----------------------------------------------------------------------

END SUBROUTINE XREFTBL
