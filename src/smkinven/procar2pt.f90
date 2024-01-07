
SUBROUTINE PROCAR2PT( NRAWBP )

    !**************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine process the area-to-point sources, adding X and Y
    !      locations and adjusting the annual and average day emissions.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !      Split from adjustinv.f 1/03 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
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

    !.......   MODULES for public variables
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: CIFIP, NPCNT, IPOSCOD, TPFLAG, INVYR,   &
                        POLVAL, CSOURC, CSCC, CINTGR,           &
                        XLOCA, YLOCA, CELLID, CEXTORL,          &
                        CISIC, CSRCTYP, CMACT, CNAICS, CSHAPE

    !.......   This module contains the cross-reference tables
    USE MODXREF, ONLY: AR2PTIDX, AR2PTTBL, AR2PTCNT

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NEM, NDY, NEF, NRP, NPPOL

    !.......  This module contains the arrays for the area-to-point x-form
    USE MODAR2PT, ONLY: AR2PTABL, NCONDSRC, REPAR2PT

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constat parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER , INTENT (INOUT) :: NRAWBP      ! no. raw records by pollutant

    !.......   Local pointers
    INTEGER, POINTER :: OLDNPCNT  ( : )       !  number of pollutants per source
    INTEGER, POINTER :: OLDIPOSCOD( : )       !  positn of pol in INVPCOD
    INTEGER, POINTER :: OLDTPFLAG ( : )       ! temporal profile types
    INTEGER, POINTER :: OLDINVYR  ( : )       ! inventory year

    REAL   , POINTER :: OLDPOLVAL( :,: )       ! emission values

    CHARACTER(FIPLEN3), POINTER :: OLDCIFIP ( : )     ! FIPS code
    CHARACTER(ALLLEN3), POINTER :: OLDCSOURC( : )     ! concat src
    CHARACTER(SCCLEN3), POINTER :: OLDCSCC  ( : )     ! scc code
    CHARACTER(SICLEN3), POINTER :: OLDCISIC ( : )     ! SIC code

    CHARACTER(MACLEN3), POINTER :: OLDCMACT  ( : )      ! MACT code
    CHARACTER(NAILEN3), POINTER :: OLDCNAICS ( : )      ! NAICS code
    CHARACTER(STPLEN3), POINTER :: OLDCSRCTYP( : )      ! source type code
    CHARACTER(EXTLEN3), POINTER :: OLDCEXTORL( : )      ! extended orl
    CHARACTER(INTLEN3), POINTER :: OLDCINTGR ( : )      ! integrate status
    CHARACTER(SHPLEN3), POINTER :: OLDCSHAPE ( : )      ! area-ff10 SHAPE-ID

    !.......   Local allocatable arrays
    INTEGER, ALLOCATABLE :: REPIDX( : )          ! index for sorting
    INTEGER, ALLOCATABLE :: REPLSCC( : )         ! left half of scc
    INTEGER, ALLOCATABLE :: REPRSCC( : )         ! right half of scc
    INTEGER, ALLOCATABLE :: REPSTA( : )          ! state fips code
    INTEGER, ALLOCATABLE :: REPPOL( : )          ! pollutant number
    REAL   , ALLOCATABLE :: REPORIGEMIS( : )     ! original emissions
    REAL   , ALLOCATABLE :: REPSUMEMIS ( : )     ! split and summed emissions

    !.......   Other local variables
    INTEGER         I,J,K,S         ! counters
    INTEGER         IOS             ! I/O error status
    INTEGER         LSTA            ! last state code
    INTEGER         LPOL            ! last pollutant code
    INTEGER         NA2PSRCS        ! no. of area-to-point sources to add
    INTEGER         NA2PRECS        ! no. of sources with pollutants to add
    INTEGER         NEWSRCPOS       ! position in new source arrays
    INTEGER         NEWRECPOS       ! position in new source w/ pollutant arrays
    INTEGER         OLDRECPOS       ! position in old source w/ pollutant arrays
    INTEGER         NREPSRCS        ! total no. of area-to-point sources
    INTEGER         TBLE            ! current area-to-point table number
    INTEGER         TSTA            ! temporary state code
    INTEGER         TPOL            ! temporary pollutant code
    INTEGER         REPPOS          ! position in reporting arrays
    INTEGER         TMPPOS          ! temporary position in reporting arrays
    INTEGER         ROW             ! current area-to-point row
    INTEGER         OLDNSRC         ! old number of srcs

    REAL            FACTOR          ! factor for area-to-point conversion

    CHARACTER(SRCLEN3)  CSRC            !  temporary source information
    CHARACTER(SCCLEN3)  LSCC            !  last scc code
    CHARACTER(SCCLEN3)  TSCC            !  temporary scc code
    CHARACTER(256)      MESG            !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'PROCAR2PT'     ! program name

    !***********************************************************************
    !   begin body of subroutine PROCAR2PT

    !.......  Determine total number of sources to be added; if source only
    !         has one location, can use existing arrays and don't need to
    !         create additional memory for it
    NA2PSRCS = 0
    NA2PRECS = 0
    NREPSRCS = 0

    DO S = 1, NSRC
        IF( AR2PTTBL( S ) /= 0 ) THEN
            NA2PSRCS = NA2PSRCS +  AR2PTCNT( S ) - 1
            NA2PRECS = NA2PRECS + (AR2PTCNT( S ) - 1)*NPCNT( S )
            NREPSRCS = NREPSRCS + NPCNT( S )
        END IF
    END DO

    !.......  Check if any sources need to be processed
    IF( NREPSRCS == 0 ) RETURN

    !.......  Allocate memory for reporting
    ALLOCATE( REPIDX ( NREPSRCS ),    &
              REPLSCC( NREPSRCS ),    &
              REPRSCC( NREPSRCS ),    &
               REPSTA( NREPSRCS ),    &
               REPPOL( NREPSRCS ),    &
          REPORIGEMIS( NREPSRCS ),    &
           REPSUMEMIS( NREPSRCS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'REPIDX...REPSUMEMIS', PROGNAME )

    !.......  Store old number of sources
    OLDNSRC = NSRC

    IF( NA2PSRCS > 0 ) THEN

        !.......  Update total number of sources and records
        NSRC   = NSRC   + NA2PSRCS
        NRAWBP = NRAWBP + NA2PRECS

        !.......  Associate temporary pointers with sorted arrays
        OLDCIFIP   => CIFIP
        OLDCISIC   => CISIC
        OLDNPCNT   => NPCNT
        OLDIPOSCOD => IPOSCOD
        OLDTPFLAG  => TPFLAG
        OLDINVYR   => INVYR

        OLDPOLVAL  => POLVAL

        OLDCSOURC  => CSOURC
        OLDCSCC    => CSCC

        OLDCSRCTYP => CSRCTYP

        IF( ASSOCIATED( CMACT ) ) THEN
            OLDCMACT   => CMACT
            OLDCNAICS  => CNAICS
        END IF

        IF( ASSOCIATED( CINTGR ) ) THEN
            OLDCINTGR   => CINTGR
        END IF

        IF( ASSOCIATED( CEXTORL ) ) THEN
            OLDCEXTORL   => CEXTORL
        END IF

        IF( ASSOCIATED( CSHAPE ) ) THEN
            OLDCSHAPE   => CSHAPE
        END IF

        !.......  Nullify original sorted arrays
        NULLIFY( CIFIP, CISIC, NPCNT, IPOSCOD, TPFLAG, INVYR,    &
                 POLVAL, CSOURC, CSCC, CMACT, CSRCTYP, CNAICS,    &
                CEXTORL, CINTGR, CSHAPE )

        !.....  Deallocate original X and Y location arrays
        !       Don't need to store old values since they aren't set
        DEALLOCATE( XLOCA, YLOCA, CELLID  )

        !.......  Allocate memory for larger sorted arrays
        ALLOCATE( CIFIP( NSRC ),    &
                  CISIC( NSRC ),    &
                  NPCNT( NSRC ),    &
                 TPFLAG( NSRC ),    &
                  INVYR( NSRC ),    &
                IPOSCOD( NRAWBP ),    &
                 POLVAL( NRAWBP,NPPOL ),    &
                 CSOURC( NSRC ),    &
                   CSCC( NSRC ),    &
                CSRCTYP( NSRC ),    &
                  XLOCA( NSRC ),    &
                  YLOCA( NSRC ),    &
                 CELLID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CIFIP...CELLID', PROGNAME )

        CSRCTYP = ' '       ! array

        IF( ASSOCIATED( OLDCMACT ) ) THEN
            ALLOCATE( CMACT( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CMACT', PROGNAME )
            ALLOCATE( CNAICS( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CNAICS', PROGNAME )

            CMACT   = ' '       ! array
            CNAICS  = ' '       ! array
        END IF

        IF( ASSOCIATED( OLDCINTGR ) ) THEN
            ALLOCATE( CINTGR( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CINTGR', PROGNAME )

            CINTGR  = ' '       ! array
        END IF

        IF( ASSOCIATED( OLDCEXTORL ) ) THEN
            ALLOCATE( CEXTORL( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CEXTORL', PROGNAME )

            CEXTORL = ' '       ! array
        END IF

        IF( ASSOCIATED( OLDCSHAPE ) ) THEN
            ALLOCATE( CSHAPE( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'CSHAPE', PROGNAME )

            CSHAPE = ' '       ! array
        END IF

        CISIC = ' '           ! array
        NPCNT = 0             ! array

        POLVAL = BADVAL3      ! array

        XLOCA = BADVAL3       ! array
        YLOCA = BADVAL3       ! array
        CELLID = 0            ! array

    END IF

    NEWSRCPOS = 0
    NEWRECPOS = 0
    OLDRECPOS = 1
    REPPOS = 0

    !.......  Loop through original sources
    DO S = 1, OLDNSRC

        !.......  Check if current source is to be processed
        IF( AR2PTTBL( S ) /= 0 ) THEN

            !.......  Set area-to-point table and row for current source
            TBLE = AR2PTTBL( S )
            ROW  = AR2PTIDX( S )

            !.......  If no sources are actually added, just modify existing values
            IF( NA2PSRCS == 0 ) THEN
                XLOCA( S ) = AR2PTABL( ROW,TBLE )%LON
                YLOCA( S ) = AR2PTABL( ROW,TBLE )%LAT

                FACTOR = AR2PTABL( ROW,TBLE )%ALLOC

                !.......  Loop through pollutants for this source
                DO K = OLDRECPOS, OLDRECPOS + NPCNT( S ) - 1

                    !.......  Store information for reporting
                    REPPOS = REPPOS + 1
                    REPIDX( REPPOS )  = REPPOS
                    REPLSCC( REPPOS ) = STR2INT( CSCC( S )( SCCEXPLEN3+1:SCCEXPLEN3+5 ) )
                    REPRSCC( REPPOS ) = STR2INT( CSCC( S )( SCCEXPLEN3+6:SCCLEN3 ) )
                    REPSTA( REPPOS )  = STR2INT( CSOURC( S )( FIPEXPLEN3+1:FIPEXPLEN3+3 ) )
                    REPPOL( REPPOS )  = IPOSCOD( K )
                    REPORIGEMIS( REPPOS ) = POLVAL( K,NEM )

                    !.....  If not adding sources, factor should be 1.0, so skip
                    !       adjusting emissions
                    IF( FACTOR /= 1. ) THEN
                        POLVAL( K,NEM ) = POLVAL( K,NEM ) * FACTOR
                        POLVAL( K,NDY ) = POLVAL( K,NDY ) * FACTOR
                    END IF

                    REPSUMEMIS( REPPOS ) = POLVAL( K,NEM )

                END DO

                OLDRECPOS = OLDRECPOS + NPCNT( S )

            ELSE

                !.......  Loop through all locations for this source
                DO J = 0, AR2PTCNT( S ) - 1

                    !.......  Increment source position and copy source info
                    NEWSRCPOS = NEWSRCPOS + 1

                    CIFIP( NEWSRCPOS ) = OLDCIFIP( S )
                    CISIC( NEWSRCPOS ) = OLDCISIC( S )
                    NPCNT( NEWSRCPOS ) = OLDNPCNT( S )
                    TPFLAG( NEWSRCPOS ) = OLDTPFLAG( S )
                    INVYR( NEWSRCPOS ) = OLDINVYR( S )
                    CSOURC( NEWSRCPOS ) = OLDCSOURC( S )
                    CSCC( NEWSRCPOS ) = OLDCSCC( S )
                    CSRCTYP( NEWSRCPOS ) = OLDCSRCTYP( S )

                    IF( ASSOCIATED( OLDCMACT ) ) THEN
                        CMACT  ( NEWSRCPOS ) = OLDCMACT  ( S )
                        CNAICS ( NEWSRCPOS ) = OLDCNAICS ( S )
                    END IF

                    IF( ASSOCIATED( OLDCINTGR ) ) THEN
                        CINTGR( NEWSRCPOS ) = OLDCINTGR( S )
                    END IF

                    IF( ASSOCIATED( OLDCEXTORL ) ) THEN
                        CEXTORL( NEWSRCPOS ) = OLDCEXTORL( S )
                    END IF

                    IF( ASSOCIATED( OLDCSHAPE ) ) THEN
                        CSHAPE( NEWSRCPOS ) = OLDCSHAPE( S )
                    END IF

                    !.......  Store X and Y locations
                    XLOCA( NEWSRCPOS ) = AR2PTABL( ROW+J,TBLE )%LON
                    YLOCA( NEWSRCPOS ) = AR2PTABL( ROW+J,TBLE )%LAT

                    !.......  Set factor for this location
                    FACTOR = AR2PTABL( ROW+J,TBLE )%ALLOC

                    !.......  Loop through pollutants for this source
                    DO K = OLDRECPOS, OLDRECPOS + OLDNPCNT( S ) - 1

                        !.......  Increment record position and store pollutant code
                        NEWRECPOS = NEWRECPOS + 1
                        IPOSCOD( NEWRECPOS ) = OLDIPOSCOD( K )

                        !.......  Adjust annual and average day emissions based on ar2pt factor
                        POLVAL( NEWRECPOS, NEM ) =  OLDPOLVAL( K,NEM ) * FACTOR

                        POLVAL( NEWRECPOS,NDY )  =  OLDPOLVAL( K,NDY ) * FACTOR

                        !.......  Copy remaining values to new array
                        POLVAL( NEWRECPOS,NEF:NRP ) = OLDPOLVAL( K,NEF:NRP )

                        !.......  Store information for reporting
                        IF( J == 0 ) THEN
                            REPPOS = REPPOS + 1
                            REPIDX( REPPOS )  = REPPOS
                            REPLSCC( REPPOS ) = STR2INT( OLDCSCC( S )( SCCEXPLEN3+1:SCCEXPLEN3+5 ) )
                            REPRSCC( REPPOS ) = STR2INT( OLDCSCC( S )( SCCEXPLEN3+6:SCCLEN3 ) )
                            REPSTA( REPPOS )  = STR2INT( OLDCSOURC( S )( FIPEXPLEN3+1:FIPEXPLEN3+3 ) )
                            REPPOL( REPPOS )  = OLDIPOSCOD( K )
                            REPORIGEMIS( REPPOS ) = OLDPOLVAL( K,NEM )
                            REPSUMEMIS( REPPOS )  = POLVAL( NEWRECPOS,NEM )
                        ELSE
                            TMPPOS = REPPOS - OLDNPCNT( S ) + K - OLDRECPOS + 1
                            REPSUMEMIS( TMPPOS ) = REPSUMEMIS( TMPPOS ) + POLVAL( NEWRECPOS,NEM )
                        END IF

                    END DO      ! loop through pollutants

                END DO      ! loop through area-to-point locations

                OLDRECPOS = OLDRECPOS + OLDNPCNT( S )

            END IF              ! end check if any sources are added

        ELSE

            !.....  Not processing current source, but if we're adding sources,
            !       then need to copy information to new arrays
            IF( NA2PSRCS > 0 ) THEN

                NEWSRCPOS = NEWSRCPOS + 1

                CIFIP( NEWSRCPOS )  = OLDCIFIP( S )
                CISIC( NEWSRCPOS )  = OLDCISIC( S )
                NPCNT( NEWSRCPOS )  = OLDNPCNT( S )
                TPFLAG( NEWSRCPOS ) = OLDTPFLAG( S )
                INVYR( NEWSRCPOS )  = OLDINVYR( S )
                CSOURC( NEWSRCPOS ) = OLDCSOURC( S )
                CSCC( NEWSRCPOS )   = OLDCSCC( S )
                CSRCTYP( NEWSRCPOS ) = OLDCSRCTYP( S )

                IF( ASSOCIATED( OLDCMACT ) ) THEN
                    CMACT  ( NEWSRCPOS ) = OLDCMACT  ( S )
                    CNAICS ( NEWSRCPOS ) = OLDCNAICS ( S )
                END IF

                IF( ASSOCIATED( OLDCINTGR ) ) THEN
                    CINTGR( NEWSRCPOS ) = OLDCINTGR( S )
                END IF

                IF( ASSOCIATED( OLDCEXTORL ) ) THEN
                    CEXTORL( NEWSRCPOS ) = OLDCEXTORL( S )
                END IF

                IF( ASSOCIATED( OLDCSHAPE ) ) THEN
                    CSHAPE( NEWSRCPOS ) = OLDCSHAPE( S )
                END IF

                DO K = OLDRECPOS, OLDRECPOS + OLDNPCNT( S ) - 1

                    NEWRECPOS = NEWRECPOS + 1
                    IPOSCOD( NEWRECPOS ) = OLDIPOSCOD( K )
                    POLVAL( NEWRECPOS,: ) = OLDPOLVAL( K,: )

                END DO

                OLDRECPOS = OLDRECPOS + OLDNPCNT( S )

            END IF

        END IF      ! check if source is processed

    END DO      ! loop through sources

    !.......  Deallocate old source and emissions arrays
    IF( NA2PSRCS > 0 ) THEN
        NULLIFY( OLDCIFIP, OLDCISIC, OLDNPCNT, OLDIPOSCOD,    &
                 OLDTPFLAG, OLDINVYR, OLDPOLVAL, OLDCSOURC,   &
                 OLDCSCC, OLDCSRCTYP )

        IF( ASSOCIATED( OLDCMACT ) ) THEN
            NULLIFY( OLDCMACT, OLDCNAICS )
        END IF

        IF( ASSOCIATED( OLDCINTGR ) ) THEN
            NULLIFY( OLDCINTGR )
        END IF

        IF( ASSOCIATED( OLDCEXTORL ) ) THEN
            NULLIFY( OLDCEXTORL )
        END IF

        IF( ASSOCIATED( OLDCSHAPE ) ) THEN
            NULLIFY( OLDCSHAPE )
        END IF

    END IF

    !.......  Sort source information for reporting
    CALL SORTI4( NREPSRCS, REPIDX, REPSTA, REPLSCC, REPRSCC, REPPOL )

    !.......  Determine total number of sources accounting for multiple
    !           pollutants and counties
    LSTA = 0
    LSCC = EMCMISS3
    LPOL = 0
    NCONDSRC = 0

    DO I = 1,NREPSRCS
        J = REPIDX( I )

        TSTA = REPSTA( J )
        WRITE( TSCC, '(2I5.5)' ) REPLSCC( J ), REPRSCC( J )
        CALL PADZERO( TSCC )
        TPOL = REPPOL( J )

        IF( TSTA /= LSTA .OR.    &
            TSCC /= LSCC .OR.    &
            TPOL /= LPOL      ) THEN
            NCONDSRC = NCONDSRC + 1
        END IF

        LSTA = TSTA
        LSCC = TSCC
        LPOL = TPOL
    END DO

    !.......  Allocate final size for reporting array
    ALLOCATE( REPAR2PT( NCONDSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'REPAR2PT', PROGNAME )
    REPAR2PT%NFIPS    = 0
    REPAR2PT%ORIGEMIS = 0.
    REPAR2PT%SUMEMIS  = 0.

    !.......  Loop through and store final reporting information
    LSTA = 0
    LSCC = EMCMISS3
    LPOL = 0
    K    = 0

    DO I = 1,NREPSRCS
        J = REPIDX( I )

        TSTA = REPSTA( J )
        WRITE( TSCC, '(2I5.5)' ) REPLSCC( J ), REPRSCC( J )
        CALL PADZERO( TSCC )
        TPOL = REPPOL( J )

        IF( TSTA /= LSTA .OR.    &
            TSCC /= LSCC .OR.    &
            TPOL /= LPOL      ) THEN

            K = K + 1
            REPAR2PT( K )%STATE = TSTA
            REPAR2PT( K )%SCC   = TSCC
            REPAR2PT( K )%POLL  = TPOL
            REPAR2PT( K )%NFIPS = 1

            REPAR2PT( K )%ORIGEMIS = REPORIGEMIS( J )
            REPAR2PT( K )%SUMEMIS  = REPSUMEMIS ( J )
        ELSE
            REPAR2PT( K )%NFIPS = REPAR2PT( K )%NFIPS + 1

            REPAR2PT( K )%ORIGEMIS = REPAR2PT( K )%ORIGEMIS + REPORIGEMIS( J )
            REPAR2PT( K )%SUMEMIS  = REPAR2PT( K )%SUMEMIS  + REPSUMEMIS( J )
        END IF

        LSTA = TSTA
        LSCC = TSCC
        LPOL = TPOL

    END DO

    !.......  Deallocate temporary reporting arrays
    DEALLOCATE( REPIDX, REPSTA, REPLSCC, REPRSCC, REPPOL, REPORIGEMIS, REPSUMEMIS )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE PROCAR2PT
