
SUBROUTINE TAGTABLE( ICSIZE, NXREF, XTYPE, XTCNT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine allocates memory and populates cross-reference tables
    !      for the tagging operation.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 1/2009 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and 
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !         System
    ! File: @(#)$Id$
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***************************************************************************

    USE M3UTILIO

    !.......  This module is for cross reference tables
    USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, CMACTA, CISICA, ISPTA,&
                       CTAGNA, NXTYPES

    !.......  This module contains the speciation profiles
    USE MODSPRO, ONLY: SPCLIST

    !.......   This module is for cross reference tables for tagging
    USE MODTAG, ONLY: TAGXCNT, TAGCHRT03, TAGCHRT04, TAGCHRT06,&
            TAGCHRT07, TAGCHRT09, TAGCHRT10, TAGCHRT11, TAGCHRT26,&
            TAGCHRT27, TAGCHRT28, TAGCHRT29, TAGCHRT30, TAGCHRT31,&
            TAGCHRT32, TAGCHRT33, TAGCHRT34, TAGCHRT35, TAGCHRT36,&
            TAGCHRT37, TAGT03, TAGT04, TAGT06, TAGT07,TAGT09,TAGT10,&
            TAGT11, TAGT26, TAGT27, TAGT28, TAGT29, TAGT30, TAGT31,&
            TAGT32, TAGT33, TAGT34, TAGT35, TAGT36, TAGT37,&
            NTAGSALL, TAGSPECIES

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: ICSIZE( * )         ! size of x-ref groups
    INTEGER, INTENT (IN) :: NXREF               ! no. ungrpd x-ref entries
    INTEGER, INTENT (IN) :: XTYPE ( NXREF )     ! group no. of x-ref entry
    INTEGER, INTENT (IN) :: XTCNT ( NXREF )     ! pos. in x-ref group

    !.......   Local field position array
    INTEGER     ENDLEN( MAX( MXARCHR3, MXMBCHR3, MXPTCHR3 ) )

    !.......   Other local variables
    INTEGER       I, J, L, K, T, V         ! counter and indices
    INTEGER       ISP                ! temporary species position in SPCLIST
    INTEGER       IOS                  ! i/o status

    LOGICAL       EFLAG                ! true: error has occurred

    CHARACTER(512)     MESG        ! message buffer

    CHARACTER(STALEN3) CSTA        ! temporary (character) state code
    CHARACTER(FIPLEN3) CFIP        ! temporary (character) FIPS code
    CHARACTER(TAGLEN3) CTAG        ! temporary tag label
    CHARACTER(SRCLEN3) CSRC        ! temporary source char string
    CHARACTER(SCCLEN3) TSCC        ! temporary SCC
    CHARACTER(SCCLEN3) SCCZERO     ! buffer for zero SCC
    CHARACTER(SICLEN3) SICZERO     ! buffer for zero SCC
    CHARACTER(SICLEN3) CSIC        ! buffer for SIC
    CHARACTER(SICLEN3) CSICL       ! buffer for left 2-digit SIC
    CHARACTER(MACLEN3) CMCT        ! buffer for MACT code

    CHARACTER(16), PARAMETER :: PROGNAME = 'TAGTABLE'     ! program name

    !***********************************************************************
    !   begin body of subroutine TAGTABLE
    EFLAG = .FALSE.
    !.......  First deallocate if these have previously been allocated
    IF ( ALLOCATED( TAGCHRT03 ) ) THEN

        DEALLOCATE( TAGCHRT03, TAGCHRT04, TAGCHRT06 )
        DEALLOCATE( TAGCHRT07, TAGCHRT09, TAGCHRT10, TAGCHRT11 )
        DEALLOCATE( TAGCHRT26, TAGCHRT27, TAGCHRT28 )
        DEALLOCATE( TAGCHRT29, TAGCHRT30, TAGCHRT31, TAGCHRT32 )
        DEALLOCATE( TAGCHRT33, TAGCHRT34, TAGCHRT35, TAGCHRT36 )
        DEALLOCATE( TAGCHRT37 )

    END IF

    !.......  Set up zero strings for SCC and SIC codes of zero
    SCCZERO = REPEAT( '0', SCCLEN3 )
    SICZERO = REPEAT( '0', SICLEN3 )

    !.......  Set the local field position array based on the source category
    ENDLEN = 1      ! array
    SELECT CASE ( CATEGORY )
      CASE( 'AREA' )
        ENDLEN( 1:MXARCHR3 ) = ARENDL3( 1:MXARCHR3 )

      CASE( 'MOBILE' )
        ENDLEN( 1:MXMBCHR3 ) = MBENDL3( 1:MXMBCHR3 )

      CASE( 'POINT' )
        ENDLEN( 1:MXPTCHR3 ) = PTENDL3( 1:MXPTCHR3 )

    END SELECT

    !.......  Allocate tables for all of the valid tagging combinations
    J = MAX( 1, ICSIZE( 3 ) )                         ! SCC=all, FIP=0
    ALLOCATE( TAGCHRT03( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT03', PROGNAME )
    ALLOCATE( TAGT03( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT03', PROGNAME )
    TAGT03 = EMCMISS3

    J = MAX( 1, ICSIZE( 4 ) )                         ! SCC=0, FIP=state
    ALLOCATE( TAGCHRT04( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT04', PROGNAME )
    ALLOCATE( TAGT04( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT04', PROGNAME )
    TAGT04 = EMCMISS3

    J = MAX( 1, ICSIZE( 6 ) )                         ! SCC=all, FIP=state
    ALLOCATE( TAGCHRT06( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT06', PROGNAME )
    ALLOCATE( TAGT06( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT06', PROGNAME )
    TAGT06 = EMCMISS3

    J = MAX( 1, ICSIZE( 7 ) )                         ! SCC=0, FIP=all
    ALLOCATE( TAGCHRT07( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT07', PROGNAME )
    ALLOCATE( TAGT07( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT07', PROGNAME )
    TAGT07 = EMCMISS3

    J = MAX( 1, ICSIZE( 9 ) )                         ! SCC=all, FIP=all
    ALLOCATE( TAGCHRT09( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT09', PROGNAME )
    ALLOCATE( TAGT09( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT09', PROGNAME )
    TAGT09 = EMCMISS3

    J = MAX( 1, ICSIZE( 10 ) )                        ! PLANT=non-blank, SCC=0
    ALLOCATE( TAGCHRT10( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT10', PROGNAME )
    ALLOCATE( TAGT10( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT10', PROGNAME )
    TAGT10 = EMCMISS3

    J = MAX( 1, ICSIZE( 11 ) )                        ! PLANT=non-blank, SCC=all
    ALLOCATE( TAGCHRT11( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT11', PROGNAME )
    ALLOCATE( TAGT11( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT11', PROGNAME )
    TAGT11 = EMCMISS3

    J = MAX( 1, ICSIZE( 26 ) )                         ! SIC=2-digit, FIP=0
    ALLOCATE( TAGCHRT26( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT26', PROGNAME )
    ALLOCATE( TAGT26( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT26', PROGNAME )
    TAGT26 = EMCMISS3

    J = MAX( 1, ICSIZE( 27 ) )                         ! SIC=all, FIP=0
    ALLOCATE( TAGCHRT27( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT27', PROGNAME )
    ALLOCATE( TAGT27( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT27', PROGNAME )
    TAGT27 = EMCMISS3

    J = MAX( 1, ICSIZE( 28 ) )                         ! SIC=2-digit, FIP=state
    ALLOCATE( TAGCHRT28( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT28', PROGNAME )
    ALLOCATE( TAGT28( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT28', PROGNAME )
    TAGT28 = EMCMISS3

    J = MAX( 1, ICSIZE( 29 ) )                         ! SIC=all, FIP=state
    ALLOCATE( TAGCHRT29( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT29', PROGNAME )
    ALLOCATE( TAGT29( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT29', PROGNAME )
    TAGT29 = EMCMISS3

    J = MAX( 1, ICSIZE( 30 ) )                         ! SIC=2-digit, FIP=all
    ALLOCATE( TAGCHRT30( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT30', PROGNAME )
    ALLOCATE( TAGT30( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT30', PROGNAME )
    TAGT30 = EMCMISS3

    J = MAX( 1, ICSIZE( 31 ) )                         ! SIC=all, FIP=all
    ALLOCATE( TAGCHRT31( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT31', PROGNAME )
    ALLOCATE( TAGT31( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT31', PROGNAME )
    TAGT31 = EMCMISS3

    J = MAX( 1, ICSIZE( 32 ) )                         ! MACT=all, FIP=0, SCC=0
    ALLOCATE( TAGCHRT32( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT32', PROGNAME )
    ALLOCATE( TAGT32( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT32', PROGNAME )
    TAGT32 = EMCMISS3

    J = MAX( 1, ICSIZE( 33 ) )                         ! MACT=all, FIP=0, SCC=all
    ALLOCATE( TAGCHRT33( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT33', PROGNAME )
    ALLOCATE( TAGT33( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT33', PROGNAME )
    TAGT33 = EMCMISS3

    J = MAX( 1, ICSIZE( 34 ) )                         ! MACT=all, FIP=state, SCC=0
    ALLOCATE( TAGCHRT34( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT34', PROGNAME )
    ALLOCATE( TAGT34( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT34', PROGNAME )
    TAGT34 = EMCMISS3

    J = MAX( 1, ICSIZE( 35 ) )                         ! MACT=all, FIP=state, SCC=all
    ALLOCATE( TAGCHRT35( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT35', PROGNAME )
    ALLOCATE( TAGT35( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT35', PROGNAME )
    TAGT35 = EMCMISS3

    J = MAX( 1, ICSIZE( 36 ) )                         ! MACT=all, FIP=all, SCC=0
    ALLOCATE( TAGCHRT36( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT36', PROGNAME )
    ALLOCATE( TAGT36( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT36', PROGNAME )
    TAGT36 = EMCMISS3

    J = MAX( 1, ICSIZE( 37 ) )                         ! MACT=all, FIP=all, SCC=all
    ALLOCATE( TAGCHRT37( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGCHRT37', PROGNAME )
    ALLOCATE( TAGT37( -1:J,NTAGSALL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGT37', PROGNAME )
    TAGT37 = EMCMISS3

    !.......  Loop through sorted tagging cross-reference table populate tables
    DO I = 1, NXREF
        J      = INDXTA( I )
        CSRC   = CSRCTA( J )
        TSCC   = CSCCTA( J )
        CTAG   = CTAGNA( J )
        ISP    = ISPTA ( J )

        V = INDEX1( SPCLIST( ISP ), NTAGSALL, TAGSPECIES )

        IF( ISP == 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Species code was 0 ' //&
                   'for tagging cross-reference entry ' //&
                   CRLF()// BLANK10// 'in program ' // PROGNAME
            CALL M3MESG( MESG )
        END IF

        IF( ALLOCATED( CMACTA ) ) THEN
            CMCT = CMACTA( J )
        ELSE
            CMCT = ' '
        END IF

        IF( ALLOCATED( CISICA ) ) THEN
            CSIC = CISICA( J )
        ELSE
            CSIC = SICZERO
        END IF

        !.......  Set up partial strings for country/state/county
        CFIP   = CSRC( 1:FIPLEN3 )
        CSTA   = CSRC( 1:STALEN3 )

        !.......  If SIC given, setup SIC fields
        IF( CSIC /= SICZERO ) THEN
            TSCC  = SCCZERO
            CSICL = CSIC( 1:SICLEN3-2 )
        END IF

        T      = XTYPE ( I )      ! extract what group this entry is in
        K      = XTCNT ( I )      ! extract position in that group

        SELECT CASE ( T )

          CASE( 3 )
            TAGCHRT03( K ) = TSCC
            TAGT03( K,V ) = CTAG
          CASE( 4 )
            TAGCHRT04( K ) = CSTA
            TAGT04( K,V ) = CTAG
          CASE( 6 )
            TAGCHRT06( K ) = CSTA // TSCC
            TAGT06( K,V ) = CTAG
          CASE( 7 )
            TAGCHRT07( K ) = CFIP
            TAGT07( K,V ) = CTAG
          CASE( 9 )
            TAGCHRT09( K ) = CFIP // TSCC
            TAGT09( K,V ) = CTAG
          CASE( 10 )
            TAGCHRT10( K ) = CSRC( 1:ENDLEN( 2 ) )
            TAGT10( K,V ) = CTAG
          CASE( 11 )
            TAGCHRT11( K ) = CSRC( 1:ENDLEN( 2 ) ) // TSCC
            TAGT11( K,V ) = CTAG
          CASE( 26 )
            TAGCHRT26( K ) = CSICL
            TAGT26( K,V ) = CTAG
          CASE( 27 )
            TAGCHRT27( K ) = CSIC
            TAGT27( K,V ) = CTAG
          CASE( 28 )
            TAGCHRT28( K ) = CSTA //CSICL
            TAGT28( K,V ) = CTAG
          CASE( 29 )
            TAGCHRT29( K ) = CSTA // CSIC
            TAGT29( K,V ) = CTAG
          CASE( 30 )
            TAGCHRT30( K ) = CFIP // CSICL
            TAGT30( K,V ) = CTAG
          CASE( 31 )
            TAGCHRT31( K ) = CFIP // CSIC
            TAGT31( K,V ) = CTAG

    !.......  MACT based cases
          CASE( 32 )
            TAGCHRT32( K ) = CMCT
            TAGT32( K,V ) = CTAG
          CASE( 33 )
            TAGCHRT33( K ) = TSCC // CMCT
            TAGT33( K,V ) = CTAG
          CASE( 34 )
            TAGCHRT34( K ) = CSTA // CMCT
            TAGT34( K,V ) = CTAG
          CASE( 35 )
            TAGCHRT35( K ) = CSTA // TSCC // CMCT
            TAGT35( K,V ) = CTAG
          CASE( 36 )
            TAGCHRT36( K ) = CFIP // CMCT
            TAGT36( K,V ) = CTAG
          CASE( 37 )
            TAGCHRT37( K ) = CFIP // TSCC // CMCT
            TAGT37( K,V ) = CTAG

          CASE DEFAULT

            EFLAG = .TRUE.
            WRITE( MESG,94010 )                                     &
                   'INTERNAL ERROR: Cross-reference category', T,   &
                   'not known in subroutine ' // PROGNAME
            CALL M3MESG( MESG )

        END SELECT

    END DO

    !.......  If error flag, then abort
    IF( EFLAG ) THEN
        MESG = 'Problem processing tagging records.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Save table sizes
    ALLOCATE( TAGXCNT( NXTYPES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'TAGXCNT', PROGNAME )
    DO I = 1, NXTYPES
        TAGXCNT( I ) = ICSIZE( I )
    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE TAGTABLE
