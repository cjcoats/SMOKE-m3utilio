
SUBROUTINE ASGNTAG( SNAM, NSRCIN, NTAG, TAGNAMLOC, TAGIDX )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !    For each source and current species, find the most specific tag
    !    that applies to each source. Do this using the grouped tables of
    !    tagging cross references from RDTAG.  The hierarchical order is
    !    defined in this subroutine, and can be determined from the in-source
    !    comments below. Once a profile code has been identified for all sources,
    !    return to the calling program.
    !
    !  PRECONDITIONS REQUIRED:
    !    Expects cross-reference tables to be set to EMCMISS3 if not defined
    !
    !  SUBROUTINE AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !    Created 2/2009 by M. Houyoux, US EPA 
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90", and
    !       related changes
    !***************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !    System
    ! File: @(#)$Id$
    !
    ! COPYRIGHT (C) 2004, Environmental Protection Agency
    !
    ! smoke@unc.edu
    !
    ! Pathname: $Source$
    ! Last updated: $Date$
    !
    !***************************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......   This module contains the source arrays
    USE MODSOURC, ONLY: CSOURC, CSCC, CISIC, CMACT

    !.......   This module contains the cross-reference tables for tagging
    USE MODTAG, ONLY: TAGXCNT, TAGCHRT03, TAGCHRT04, TAGCHRT06,     &
            TAGCHRT07, TAGCHRT09, TAGCHRT10, TAGCHRT11, TAGCHRT26,  &
            TAGCHRT27, TAGCHRT28, TAGCHRT29, TAGCHRT30, TAGCHRT31,  &
            TAGCHRT32, TAGCHRT33, TAGCHRT34, TAGCHRT35, TAGCHRT36,  &
            TAGCHRT37, TAGT03, TAGT04,                              &
            TAGT06, TAGT07, TAGT09, TAGT10,                         &
            TAGT11, TAGT26, TAGT27, TAGT28, TAGT29, TAGT30, TAGT31, &
            TAGT32, TAGT33, TAGT34, TAGT35, TAGT36, TAGT37,         &
            NTAGSALL, TAGSPECIES

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT    (IN) :: SNAM        ! species name of interest
    INTEGER     , INTENT    (IN) :: NSRCIN      ! number of sources
    INTEGER     , INTENT    (IN) :: NTAG        ! number of tags for this species
    CHARACTER(*), INTENT    (IN) :: TAGNAMLOC( NTAG )     ! tag names for current species/pollutant
    INTEGER     , INTENT   (OUT) :: TAGIDX( NSRCIN )    ! tag indices

    !.......  Other local variables
    INTEGER         S, V        !  counters and indices

    INTEGER          F0, F1, F2, F3, F4, F5, F6      ! tmp find indices
    INTEGER       :: F0B = 0          ! extra find index for mobile
    INTEGER       :: F2B = 0          ! extra find index for mobile
    INTEGER       :: F4B = 0          ! extra find index for mobile
    INTEGER          IOS              ! i/o status
    INTEGER       :: NWARN=0          ! current number of warnings of each type to write
    INTEGER, SAVE :: MXWARN           ! maximum number of warnings of each type to write

    REAL             CNVFAC           ! tmp pol-to-pol conversion factor

    LOGICAL       :: EFLAG    = .FALSE.     ! true: error detected
    LOGICAL       :: LRDCOMBO               ! true: for first COMBO for a pollutant/emis type
    LOGICAL, SAVE :: FIRSTIME = .TRUE.      ! true: first time subrtn called
    LOGICAL, SAVE :: MACTFLAG = .FALSE.     ! true: MACT codes available in inventory
    LOGICAL, SAVE :: LNOZERO  = .FALSE.     ! true: skip speication for zero emission
    LOGICAL, SAVE :: SICFLAG  = .FALSE.     ! true: SIC available in inventory

    CHARACTER(300)       BUFFER      ! source fields buffer
    CHARACTER(300)       MESG        ! message buffer
    CHARACTER(FIPLEN3)   CFIP        ! tmp (character) FIPS code
    CHARACTER(STALEN3)   CSTA        ! tmp Country/state code
    CHARACTER(SRCLEN3)   CSRC        ! tmp source chars string
    CHARACTER(SCCLEN3)   TSCC        ! tmp 10-digit SCC
    CHARACTER(SCCLEN3)   TSCCINIT    ! tmp initial 10-digit SCC
    CHARACTER(SS0LEN3):: CHK11=' '   ! tmp FIPS // Plant // SCC
    CHARACTER(FPLLEN3):: CHK10=' '   ! tmp FIPS code // plant id
    CHARACTER(FPSLEN3):: CHK09=' '   ! tmp FIPS code // SCC
    CHARACTER(STSLEN3):: CHK06=' '   ! tmp Country/state code // SCC
    CHARACTER(STILEN3):: CHK28=' '   ! tmp Country/state code // left SIC
    CHARACTER(STILEN3):: CHK29=' '   ! tmp Country/state code // SIC
    CHARACTER(FPILEN3):: CHK30=' '   ! tmp FIPS code // left SIC
    CHARACTER(FPILEN3):: CHK31=' '   ! tmp FIPS code // SIC
    CHARACTER(MSCLEN3):: CHK33=' '   ! tmp SCC // MACT
    CHARACTER(MSTLEN3):: CHK34=' '   ! tmp Country/state code // MACT
    CHARACTER(MSSLEN3):: CHK35=' '   ! tmp Country/state code // SCC // MACT
    CHARACTER(MFPLEN3):: CHK36=' '   ! tmp FIPS code // MACT
    CHARACTER(MFSLEN3):: CHK37=' '   ! tmp FIPS code // SCC // MACT
    CHARACTER(MACLEN3)   CMCT        ! tmp MACT code
    CHARACTER(SICLEN3)   CSIC        ! tmp SIC code
    CHARACTER(SICLEN3)   CSICL       ! tmp left SIC code
    CHARACTER(VIDLEN3)   CVID        ! tmp vehicle type

    CHARACTER(16), PARAMETER :: PROGNAME = 'ASGNTAG'     ! program name

    !***********************************************************************
    !   begin body of subroutine ASGNTAG

    !.......  For first time routine is called in all cases,
    IF( FIRSTIME ) THEN

    !.......  Figure out if SIC and/or MACT codes are available
        IF ( ASSOCIATED ( CISIC ) ) SICFLAG  = .TRUE.
        IF ( ASSOCIATED ( CMACT ) ) MACTFLAG = .TRUE.

        FIRSTIME = .FALSE.

    ENDIF

    !.......  Initialize source tagging index to 0.
    TAGIDX = 0        ! array

    !.......  Find index in complete list of species that are getting tags
    V  = INDEX1( SNAM, NTAGSALL, TAGSPECIES )

    DO S = 1, NSRCIN

        CSRC  = CSOURC( S )
        CFIP  = CSRC( 1:FIPLEN3 )
        CSTA  = CFIP( 1:STALEN3 )
        TSCC  = CSCC( S )

        !.......  Set type of SCC
        CHK09 = CFIP // TSCC
        CHK06 = CSTA // TSCC

        IF( SICFLAG ) THEN
            CSIC = CISIC( S )
            CSICL = CSIC( SICEXPLEN3+1:SICEXPLEN3+2 )
            CHK28 = CSTA // CSICL
            CHK29 = CSTA // CSIC
            CHK30 = CFIP // CSICL
            CHK31 = CFIP // CSIC
        END IF

        IF( MACTFLAG ) THEN
            CMCT  = CMACT( S )
            CHK33 = TSCC // CMCT
            CHK34 = CSTA // CMCT
            CHK35 = CSTA // TSCC // CMCT
            CHK36 = CFIP // CMCT
            CHK37 = CFIP // TSCC // CMCT
        END IF

        TSCCINIT = TSCC

        !.......  Create selection
        IF( CATEGORY == 'POINT' ) THEN

            CHK11   = CSRC( 1:PTENDL3( 2 ) ) // TSCC
            CHK10   = CSRC( 1:PTENDL3( 2 ) )

        END IF

        !.................
        !.......  Now find and apply speciation profiles data
        !.................

        !.......  In the tables used in the following heirarchy, all cross-
        !    reference entries are by definition, pollutant- specific.
        !    The cross-reference tables (e.g,, TAGCHRT03 come from MODTAG)

        !.......  Try for pollutant-specific PLANT non-blank// SCC match; then
        !    pollutant-specific PLANT non-blank       match
        F1 = FINDC( CHK11, TAGXCNT( 11 ), TAGCHRT11 )
        F0 = FINDC( CHK10, TAGXCNT( 10 ), TAGCHRT10 )

        IF( F1 .GT. 0 .AND. TAGT11(F1,V) /= EMCMISS3 ) THEN
            TAGIDX(S) = INDEX1( TAGT11( F1,V ), NTAG, TAGNAMLOC )
            CYCLE                           !  to end of sources-loop

        ELSEIF( F0 .GT. 0 .AND. TAGT10(F0,V) /= EMCMISS3 ) THEN
            TAGIDX(S) = INDEX1( TAGT10( F0,V ), NTAG, TAGNAMLOC )
            CYCLE                           !  to end of sources-loop

        END IF

        !.......  If MACT available in inventory...
        IF ( MACTFLAG ) THEN

            !.......  Try for pollutant-specific FIPS code, SCC match,  MACT code; then
            !           pollutant-specific FIPS code  MACT code; then
            !           pollutant-specific Cy/st code, SCC match,  MACT code; then
            !           pollutant-specific Cy/st code  MACT code; then
            !           pollutant-specific SCC match  MACT code; then
            !           pollutant-specific MACT code
            F5 = FINDC( CHK37, TAGXCNT( 37 ), TAGCHRT37 )
            F4 = FINDC( CHK36, TAGXCNT( 36 ), TAGCHRT36 )
            F3 = FINDC( CHK35, TAGXCNT( 35 ), TAGCHRT35 )
            F2 = FINDC( CHK34, TAGXCNT( 34 ), TAGCHRT34 )
            F1 = FINDC( CHK33, TAGXCNT( 33 ), TAGCHRT33 )
            F0 = FINDC( CMCT , TAGXCNT( 32 ), TAGCHRT32 )

            IF( F5 .GT. 0 .AND. TAGT37( F5,V ) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT37( F5,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F4 .GT. 0 .AND. TAGT36(F4,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT36( F4,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F3 .GT. 0 .AND. TAGT35(F3,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT35( F3,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F2 .GT. 0 .AND. TAGT34(F2,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT34( F2,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F1 .GT. 0 .AND. TAGT33(F1,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT33( F1,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F0 .GT. 0 .AND. TAGT32(F0,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT32( F0,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            END IF
        END IF

        !.......  If SIC available in inventory...
        IF ( SICFLAG ) THEN

            !.......  Try for pollutant-specific FIPS code  SIC match; then
            !           pollutant-specific FIPS code  left SIC match; then
            !           pollutant-specific Cy/st code  SIC match; then
            !           pollutant-specific Cy/st code  left SIC match; then
            !           pollutant-specific SIC match; then
            !           pollutant-specific left SIC match

            F5 = FINDC( CHK31, TAGXCNT( 31 ), TAGCHRT31 )
            F4 = FINDC( CHK30, TAGXCNT( 30 ), TAGCHRT30 )
            F3 = FINDC( CHK29, TAGXCNT( 29 ), TAGCHRT29 )
            F2 = FINDC( CHK28, TAGXCNT( 28 ), TAGCHRT28 )
            F1 = FINDC( CSIC , TAGXCNT( 27 ), TAGCHRT27 )
            F0 = FINDC( CSICL, TAGXCNT( 26 ), TAGCHRT26 )

            IF( F5 .GT. 0 .AND. TAGT31( F5,V ) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT31( F5,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F4 .GT. 0 .AND. TAGT30(F4,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT30( F4,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F3 .GT. 0 .AND. TAGT29(F3,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT29( F3,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F2 .GT. 0 .AND. TAGT28(F2,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT28( F2,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F1 .GT. 0 .AND. TAGT27(F1,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT27( F1,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            ELSEIF(F0 .GT. 0 .AND. TAGT26(F0,V) /= EMCMISS3 ) THEN
                TAGIDX(S)= INDEX1( TAGT26( F0,V ), NTAG, TAGNAMLOC )
                CYCLE                           !  to end of sources-loop

            END IF
        END IF

        !.......  Try for pollutant-specific FIPS code  SCC match; then
        !           pollutant-specific Cy/st code  SCC match; then
        !           pollutant-specific SCC match

        F5 = FINDC( CHK09, TAGXCNT( 9 ), TAGCHRT09 )
        F3 = FINDC( CHK06, TAGXCNT( 6 ), TAGCHRT06 )
        F1 = FINDC( TSCC , TAGXCNT( 3 ), TAGCHRT03 )

        IF( F5 .GT. 0 .AND. TAGT09(F5,V) /= EMCMISS3 ) THEN
            TAGIDX(S)= INDEX1( TAGT09( F5,V ), NTAG, TAGNAMLOC )
            CYCLE                           !  to end of sources-loop

        ELSEIF( F3 .GT. 0 .AND. TAGT06(F3,V) /= EMCMISS3 ) THEN
            TAGIDX(S)= INDEX1( TAGT06( F3,V ), NTAG, TAGNAMLOC )
            CYCLE                           !  to end of sources-loop

        ELSEIF( F1 .GT. 0 .AND. TAGT03(F1,V) /= EMCMISS3 ) THEN
            TAGIDX(S)= INDEX1( TAGT03( F1,V ), NTAG, TAGNAMLOC )
            CYCLE                               !  to end of sources-loop

        END IF

        !.......  Try for any FIPS code match
        F0 = FINDC( CFIP, TAGXCNT( 7 ), TAGCHRT07 )

        IF( F0 .GT. 0 .AND. TAGT07(F0,V) /= EMCMISS3 ) THEN
            TAGIDX(S)= INDEX1( TAGT07( F0,V ), NTAG, TAGNAMLOC )
            CYCLE                               !  to end of sources-loop
        END IF

        !.......  Try for any country/state code match (not, pol-specific)
        F0 = FINDC( CSTA, TAGXCNT( 4 ), TAGCHRT04 )

        IF( F0 .GT. 0 .AND. TAGT04(F0,V) /= EMCMISS3 ) THEN
            TAGIDX(S)= INDEX1( TAGT04( F0,V ), NTAG, TAGNAMLOC )
            CYCLE                           !  to end of sources-loop
        END IF

    END DO            !  end loop on source, S

    IF( EFLAG ) THEN
        MESG = 'Problem assigning tags to sources'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94300 FORMAT( A, I2.2, A, I2.2, A )

END SUBROUTINE ASGNTAG
