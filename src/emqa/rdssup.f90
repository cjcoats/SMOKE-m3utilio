
SUBROUTINE RDSSUP ( PDEV )

    !***********************************************************************
    !  subroutine body starts at line  84
    !
    !  DESCRIPTION:
    !       Read new-format SSUP files and places content into MODSOURC
    !       sparse matrix profiles-and-fractins data structures
    !       <NSPFRC,SPPOLS(:),SPPNS(:,:),SPPROF(:),SPFRAC(:)>
    !
    !  PRECONDITIONS REQUIRED:
    !       setenv [AMP]SSUP   <path>
    !
    !  REVISION  HISTORY:
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
    !***********************************************************************

    USE M3UTILIO

    !.......   MODULES for public variables
    !.......  This module contains the speciation-profiles matrix, among other things.
    USE MODSOURC, ONLY : NSPFRC, SPPNLO, SPPNHI, SPPROF, SPFRAC,&
                         NGSPRO, GSPROID

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CRL, NSRC, NIPPA

    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: RPT_, NSPCPOL, SPCPOL, LSPCPOL


    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !......  External Function
    INTEGER,EXTERNAL :: GETFLINE

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: PDEV       ! unit no.: speciation supplemental

    !.......  Local variables
    INTEGER             IOS, IS1, IS2
    CHARACTER(NAMLEN3)  ::   CBUF

    CHARACTER(512) ::   LINE
    CHARACTER(512) ::   MESG

    CHARACTER(32) :: TPROF( 50 ), TFACS( 50 )

    INTEGER     IREC, PCNT
    INTEGER     I, J, K, L, L1, L2, N, N1, N2, NS, M, S, V

    CHARACTER(16), PARAMETER :: PNAME = 'RDSSUP'     ! subroutine name

    !***********************************************************************
    !   begin body of subroutine RDSSUP

    !.......  Unique list of speciation profiles
    ALLOCATE( GSPROID( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GSPROID', PNAME )
    GSPROID = ''

    !.......  Open SSUP speciation supplemental file
    MESG = 'Supplemental speciation file'
    N = GETFLINE( PDEV, MESG )

    IREC   = 0
    NSPFRC = 0
    NGSPRO = 0

    DO  I = 1, N

        READ( PDEV, '(A)', IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG,94010 ) 'I/O error', IOS, 'reading ' //&
                  'supplemental speciation file at line', IREC
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        !.......  See if this line is a pollutant name
        L1 = INDEX( LINE, '"' )                ! Find start quote

        !.......  If pollutant name, figure out which pollutant index and
        !               reset source counter to 0.
        IF ( L1 .GT. 0 ) THEN

            L2 = LEN_TRIM( LINE )             ! Find end quote
            CBUF = LINE( L1+1:L2-1 )

            !.......  Check if this pollutant is one selected for reporting
            V = INDEX1( CBUF, NSPCPOL, SPCPOL )
            IF ( V .GT. 0 ) LSPCPOL( V ) = .TRUE.

        !.......  If not pollutant name, then continue to read in the
        !         pollutant codes and store them by source
        ELSE IF ( V .GT. 0 ) THEN

            N1 = INDEX( LINE, 'NFRAC=' )
            N2 = N1 + LEN( 'NFRAC=' )
            READ( LINE(N2: ), *, IOSTAT=IS1 ) PCNT
            IF( N1 > 0 ) THEN
                S = S + 1
                NSPFRC = NSPFRC + PCNT

            ELSE
                READ( LINE, * ) ( TPROF( K ), TFACS( K ), K = 1, PCNT )

                !.......  Count the unique list of GSPRO
                DO K = 1, PCNT
                    IF( INDEX1( TPROF( K ), NSRC, GSPROID ) < 1 ) THEN
                        NGSPRO = NGSPRO + 1
                        GSPROID( NGSPRO ) = TPROF( K )
                    END IF
                END DO
            END IF

            IF( S .EQ. NSRC ) V = 0

        END IF

    END DO

        !    !......  Allocate output arrays
    DEALLOCATE( GSPROID )
    ALLOCATE( SPPNLO( NSRC ),&
              SPPNHI( NSRC ),&
              SPPROF( NSPFRC ),&
              SPFRAC( NSPFRC ),&
             GSPROID( NGSPRO ), STAT = IOS )
    CALL CHECKMEM( IOS, 'SPPNLO...SPFRAC', PNAME )
    SPFRAC = 1.0
    GSPROID = ''

        !    !......  Read/Build output arrays for fractions, fractional profiles
    REWIND( PDEV )

    V    = 0
    M    = 0
    S    = 0
    NS   = 0
    IREC = 0

    DO I = 1, N

        READ( PDEV, '(A)', IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            WRITE( MESG,94010 )&
                 'I/O error', IOS, 'reading supplemental  speciation file at line', IREC
            CALL M3EXIT( PNAME, 0, 0, MESG, 2 )
        END IF

        !.......  See if this line is a pollutant name
        L1 = INDEX( LINE, '"' )                ! Find start quote

        !.......  If pollutant name, figure out which pollutant index and
        !               reset source counter to 0.
        IF ( L1 .GT. 0 ) THEN

            L2 = LEN_TRIM( LINE )             ! Find end quote
            CBUF = LINE( L1+1:L2-1 )

            !.......  Check if this pollutant is one selected for reporting
            V = INDEX1( CBUF, NSPCPOL, SPCPOL )
            IF ( V .GT. 0 ) THEN
                LSPCPOL( V ) = .TRUE.
                S = 0
            END IF

        !.......  If not pollutant name, then continue to read in the
        !                   pollutant codes and store them by source
        ELSE IF ( V .GT. 0 ) THEN

            N1 = INDEX( LINE, 'NFRAC=' )
            N2 = N1 + LEN( 'NFRAC=' )
            READ( LINE(N2: ), *, IOSTAT=IS1 ) PCNT
            IF( N1 > 0 ) THEN
                S = S + 1
                SPPNLO( S ) = M + 1
                SPPNHI( S ) = M + PCNT
                M = M + PCNT

            ELSE
                READ( LINE, * ) ( SPPROF( K ), SPFRAC( K ), K = SPPNLO(S), SPPNHI(S) )

                !.......  Develop the unique list of GSPRO
                DO K = SPPNLO(S), SPPNHI(S)
                    IF( INDEX1( SPPROF( K ), NGSPRO, GSPROID ) < 1 ) THEN
                        NS = NS + 1
                        GSPROID( NS ) = SPPROF( K )
                    END IF
                END DO

            END IF

            IF( S == NSRC+1 ) V = 0

        END IF

    END DO

    REWIND( PDEV )

    RETURN

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE RDSSUP
