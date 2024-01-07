
SUBROUTINE PROCPKTS( PDEV, CDEV, GDEV, LDEV, MDEV, WDEV, CPYEAR,    &
                     PKTTYP, ENAME, LPSASGN, USEPOL, SFLAG, LPTMP, LCTMP )

    !***********************************************************************
    !  subroutine body starts at line 116
    !
    !  DESCRIPTION:
    !      This subroutine is responsible for processing control packet data
    !      that has already been grouped into the control cross-reference tables
    !      and control data tables.  For control packets that do not depend
    !      on pollutants, it calls routines to assign the control to the
    !      appropriate sources and generate the appropriate matrix.  For control
    !      packets that do depend on pollutants, it creates the pollutant
    !      group structure and processes the packet for all pollutant groups.  The
    !      index to the control data tables is stored for all pollutants in
    !      temporary files for later use, after the output files have been
    !      opened.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Started 3/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !************************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !       System
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

    !.....  MODULES for public variables
    !.....  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC

    !.....  This module is for cross reference tables
    USE MODXREF, ONLY: ASGNINDX

    !.....  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: PNAMMULT, PNAMPROJ, FACTOR,                     &
                        BACKOUT, DATVAL, GRPINDX, GRPFLAG, GRPSTIDX,    &
                        GRPCHAR, GRPINEM, GRPOUTEM, POLSFLAG,           &
                        NVPROJ, NVCMULT, PCTLFLAG

    !.....  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NSRC, NIPPA, NIPOL, NIACT, NPPOL,      &
                       EINAM, EANAM, ACTVTY

    IMPLICIT NONE

    !.....   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....   SUBROUTINE ARGUMENTS:

    INTEGER     , INTENT (IN) :: PDEV          ! file unit no. for tmp PROJ file
    INTEGER     , INTENT (IN) :: CDEV          ! file unit no. for tmp CTL file
    INTEGER     , INTENT (IN) :: GDEV          ! file unit no. for tmp CTG file
    INTEGER     , INTENT (IN) :: LDEV          ! file unit no. for tmp ALW file
    INTEGER     , INTENT (IN) :: MDEV          ! file unit no. for tmp MACT file
    INTEGER     , INTENT (IN) :: WDEV          ! file unit no. for warnings/error file
    INTEGER     , INTENT (IN) :: CPYEAR        ! year to project to
    CHARACTER(*), INTENT (IN) :: PKTTYP        ! packet type
    CHARACTER(*), INTENT (IN) :: ENAME         ! inventory file name
    LOGICAL     , INTENT (IN) :: LPSASGN       ! true: matrix needs to be by pollutant
    LOGICAL , INTENT (IN OUT) :: USEPOL( NIPPA )     ! true: pol in current pkt
    LOGICAL     , INTENT(OUT) :: SFLAG         ! true: at least one packet done
    LOGICAL     , INTENT(OUT) :: LPTMP         ! true: projection tmp file written
    LOGICAL     , INTENT(OUT) :: LCTMP         ! true: control tmp file written

    !.....  Reshaped inventory pollutants and associated variables
    !       INTEGER         NGRP                    ! number of pollutant groups
    !       INTEGER         NGSZ                    ! number of pollutants per group
    !       INTEGER               , ALLOCATABLE:: IPSTAT ( : )       ! pol status (0|1)

    !.....   Other local variables

    INTEGER         I, J, K, L, L1, L2, S, V          ! counters and indices

    INTEGER         IOS                       ! i/o error status
    INTEGER         NGRP                      ! number of reporting groups
    INTEGER         SCCBEG                    ! begining of SCC in CSOURC string
    INTEGER         SCCEND                    ! end of SCC in CSOURC string
    INTEGER         VIDXMULT( NIPPA )         ! pollutant flags for control
    INTEGER         VIDXPROJ( NIPPA )         ! pollutant flags for projection

    LOGICAL       :: DSFLAG   = .FALSE.       ! true: this call to ASGNCNTL has data-specific match
    LOGICAL          EFLAG                    ! error flag
    LOGICAL, SAVE :: FIRSTIME = .TRUE.        ! true: first time routine called
    LOGICAL, SAVE :: OFLAG(NPACKET) = .FALSE.       ! true: tmp file has not been opened
    LOGICAL       :: CCHEK                    ! temporary control tmp file status per packet

    CHARACTER(5)    CPOS            ! tmp sorted position of pol
    CHARACTER(256)  LINE            ! read buffer for a line
    CHARACTER(256)  MESG            ! message buffer

    CHARACTER(NAMLEN3), SAVE :: RPOL     ! pol name for reactivity controls
    CHARACTER(FPLLEN3) CPLT           ! tmp point src info through plant
    CHARACTER(FPLLEN3) PPLT           ! previous CPLT
    CHARACTER(STALEN3) CSTA           ! tmp char state
    CHARACTER(STALEN3) PSTA           ! previous char state
    CHARACTER(SCCLEN3) TSCC           ! tmp SCC
    CHARACTER(SCCLEN3) PSCC           ! previous SCC

    CHARACTER(16), PARAMETER :: PROGNAME = 'PROCPKTS'     ! program name

    !***********************************************************************
    !   Begin body of subroutine PROCPKTS

    IF( FIRSTIME ) THEN

        SFLAG = .FALSE.         ! Initialize status as no packets applied
        FIRSTIME = .FALSE.

    END IF
    
    EFLAG    = .FALSE.
    !.....  Reactivity packet...
    IF( PKTTYP .EQ. 'REACTIVITY' ) THEN

        !.....  Get environment variable setting for reactivity pollutant
        MESG = 'Pollutant for creating reactivity matrix'
        CALL ENVSTR( 'REACTIVITY_POL', MESG, 'VOC', RPOL, IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "REACTIVITY_POL"', 2 )
        END IF

        !.....  Make sure that the pollutant for the reactivity packet is
        !       in the inventory
        J = INDEX1( RPOL, NIPOL, EINAM )
        IF( J .LE. 0 ) THEN

            MESG = 'Environment variable "REACTIVITY_POL" is set '//    &
                   'to pollutant "' // TRIM( RPOL ) // '",'  //         &
                   CRLF() // BLANK10 // 'but this pollutant is not in the inventory'
            CALL M3MSG2( MESG )

            MESG = 'Reactivity matrix creation skipped.'
            CALL M3MSG2( MESG )
            RETURN

        END IF

        !.....  Generate reactivity matrices
        USEPOL = .TRUE.          ! array
        CALL GENREACT( CPYEAR, ENAME, RPOL, USEPOL )

        SFLAG = .TRUE.

    !.....  If projection or control packet, then allocate memory for
    !       needed arrays, pollutant-specific arrays, and output matrices.
    ELSE

        !......  Allocate arrays for managing output pollutants/activities
        !......  Create array for names of pollutants that receive projections

        !.....  Create array of flags indicating which controls are
        !       applied to each pollutant receiving at least one type
        !       of control or projection
        IF( .NOT. ALLOCATED( PNAMMULT ) ) THEN

        !......  Create array for names of pollutants that receive controls
            ALLOCATE( PNAMMULT( NIPPA ),        &
                      PNAMPROJ( NIPPA ),        &
                      PCTLFLAG( NIPPA, 4 ), STAT=IOS )
            CALL CHECKMEM( IOS, 'PCTLFALG', PROGNAME )
            PNAMMULT = ' '        ! array
            PNAMPROJ = ' '        ! array
            PCTLFLAG = .FALSE.        ! array

        END IF

        !.....  Allocate for work output factor for both projection and control
        IF( .NOT. ALLOCATED( FACTOR ) ) THEN
            ALLOCATE( ASGNINDX( NSRC ),     &
                        FACTOR( NSRC ), STAT=IOS )
            CALL CHECKMEM( IOS, 'FACTOR', PROGNAME )
        END IF
        ASGNINDX = 0
        FACTOR = 1.          !  array

        !.....  Allocate for work arrays for multiplicative controls
        IF( PKTTYP .NE. 'PROJECTION' ) THEN

            !.....  Allocate main work arrays
            IF( .NOT. ALLOCATED( BACKOUT ) ) THEN
                ALLOCATE( BACKOUT( NSRC ),&
                          DATVAL( NSRC,NPPOL ), STAT=IOS )
                CALL CHECKMEM( IOS, 'DATVAL', PROGNAME )
            END IF
            BACKOUT = 0.               ! array
            DATVAL  = 0.               ! array

            !.....  Allocate first set of reporting arrays
            IF( .NOT. ALLOCATED( GRPINDX ) ) THEN
                ALLOCATE( GRPINDX( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'GRPINDX', PROGNAME )
            END IF
            IF( .NOT. ALLOCATED( GRPSTIDX ) ) THEN
                ALLOCATE( GRPSTIDX( NSRC ),&
                           GRPCHAR( NSRC ), STAT=IOS )
                CALL CHECKMEM( IOS, 'GRPCHAR', PROGNAME )
                GRPINDX  = 0              ! array
            END IF

            !.....  If haven't already, get set up for group reporting...
            IF ( .NOT. ALLOCATED( GRPFLAG ) ) THEN

                !......  Count the number of groups in the inventory
                IF( CATEGORY .EQ. 'POINT' ) THEN

                    PPLT = ' '
                    NGRP = 0
                    DO S = 1, NSRC
                        CPLT = CSOURC( S )( 1:FPLLEN3 )
                        IF( CPLT .NE. PPLT ) THEN
                            NGRP = NGRP + 1
                            PPLT = CPLT
                        END IF
                        GRPINDX ( S ) = NGRP
                        GRPSTIDX( S ) = S         ! Needed for loops, but not used to sort
                    END DO

                ELSE

                    IF( CATEGORY .EQ. 'AREA' ) THEN
                        SCCBEG = ARBEGL3( 2 )
                        SCCEND = ARENDL3( 2 )
                    ELSE                ! MOBILE
                        SCCBEG = MBBEGL3( 5 )
                        SCCEND = MBENDL3( 5 )
                    END IF

                    !......  Build and sort source array for SCC-state grouping
                    DO S = 1, NSRC
                        CSTA = CSOURC( S )( 1     :STALEN3 )
                        TSCC = CSOURC( S )( SCCBEG:SCCEND  )

                        GRPSTIDX( S ) = S
                        GRPCHAR ( S ) = CSTA // TSCC
                    END DO

                    CALL SORTIC( NSRC, GRPSTIDX, GRPCHAR )

                    !......  Count the number of state/SCCs in the domain
                    PSTA = ' '
                    PSCC = ' '
                    SCCBEG = STALEN3 + 1
                    SCCEND = STALEN3 + SCCLEN3
                    DO S = 1, NSRC
                        J = GRPSTIDX( S )
                        CSTA = GRPCHAR( J )( 1     :STALEN3 )
                        TSCC = GRPCHAR( J )( SCCBEG:SCCEND  )
                        IF( CSTA .NE. PSTA .OR.     &
                            TSCC .NE. PSCC       ) THEN
                            NGRP = NGRP + 1
                            PSTA = CSTA
                            PSCC = TSCC
                        END IF
                        GRPINDX( J ) = NGRP
                    END DO

                END IF        ! If point or non-point sources

                IF( .NOT. ALLOCATED( GRPFLAG ) ) THEN
                    ALLOCATE( GRPFLAG( NGRP ),          &
                              GRPINEM( NGRP, NIPPA ),   &
                             GRPOUTEM( NGRP, NIPPA ), STAT=IOS )
                    CALL CHECKMEM( IOS, 'GRPOUTEM', PROGNAME )
                END IF
            END IF        ! If group information not previously allocated

            GRPINEM  = 0.     ! array
            GRPOUTEM = 0.     ! array
            GRPFLAG  = .FALSE.      ! array

        END IF           ! For multiplicative control packets

        !.....  Create array for indicating the status of pollutants at each iteration

    END IF           ! For reactivity or not

    !.....  If the packet does not have pol/act-specific entries, then
    !       create using the "all" keyword.
    !.....  For projection matrix only for now...
    IF ( PKTTYP .EQ. 'PROJECTION' .AND. .NOT. LPSASGN ) THEN

        CALL ASGNCNTL( NSRC, WDEV, PKTTYP, 'all', DSFLAG, ASGNINDX )

        IF ( .NOT. OFLAG(1) ) THEN
            CALL OPENCTMP( PKTTYP, PDEV )
            OFLAG(1) = .TRUE.
        END IF
        CALL WRCTMP( PDEV, 1, ASGNINDX, 1, LPTMP )

        SFLAG = .TRUE.

    !.....  For projection for specific pollutants and multiplicative controls...
    ELSE IF ( PKTTYP .NE. 'REACTIVITY' ) THEN

        !.....  Loop through the pollutant groups...

        !.....  Apply the current packets to appropriate sources and pollutants
        !       in the inventory.
        !.....  Write temporary ASCII files containing the indices to the
        !       control data tables.  This is because we only want to write out
        !       the I/O API control matrices for the pollutants that are
        !       actually affected by controls, but don't know which pollutants
        !       to open the output file(s) with until all of the pollutants have
        !       been processed.

        DO V = 1, NIPPA

            !.....  Skip pollutants that do not apply to this packet.
            IF( .NOT. USEPOL( V ) ) CYCLE

            SELECT CASE( PKTTYP )

              CASE( 'PROJECTION' )

                !......  Generate projection matrix
                CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), DSFLAG, ASGNINDX )
                IF( DSFLAG ) POLSFLAG = .TRUE.

                VIDXPROJ = 0         ! array
                CALL UPDATE_POLLIST( V, ASGNINDX, 0, VIDXPROJ, NVPROJ, PNAMPROJ )

                IF ( .NOT. OFLAG(1) ) THEN
                    CALL OPENCTMP( PKTTYP, PDEV )
                    OFLAG(1) = .TRUE.
                END IF
                CALL WRCTMP( PDEV, V, ASGNINDX, VIDXPROJ, LPTMP )

                SFLAG = .TRUE.

              CASE( 'CTG' )

                VIDXMULT = 0     ! array
                CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), DSFLAG, ASGNINDX )

                CALL UPDATE_POLLIST( V, ASGNINDX, 2, VIDXMULT, NVCMULT, PNAMMULT )
                IF ( .NOT. OFLAG(2) ) THEN
                    CALL OPENCTMP( PKTTYP, GDEV )
                    OFLAG(2) = .TRUE.
                END IF

                CALL WRCTMP( GDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                LCTMP = ( LCTMP .OR. CCHEK )

                SFLAG = .TRUE.

              CASE( 'CONTROL' )

                !.....  Skip activities because they
                !       do not have the base-year control effectiveness
                J = INDEX1( EANAM( V ), NIACT, ACTVTY )
                IF ( J .GT. 0 ) THEN
                    MESG = 'Skipping activity "' //                 &
                               TRIM( EANAM( V ) )// '" since '//    &
                               'CONTROL packet cannot apply.'
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

                VIDXMULT = 0     ! array
                CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), DSFLAG, ASGNINDX )

                CALL UPDATE_POLLIST( V, ASGNINDX, 1, VIDXMULT, NVCMULT, PNAMMULT )
                IF ( .NOT. OFLAG(3) ) THEN
                    CALL OPENCTMP( PKTTYP, CDEV )
                    OFLAG(3) = .TRUE.
                END IF
                CALL WRCTMP( CDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                LCTMP = ( LCTMP .OR. CCHEK )

                SFLAG = .TRUE.

              CASE( 'ALLOWABLE' )

                VIDXMULT = 0     ! array
                CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), DSFLAG, ASGNINDX )

                CALL UPDATE_POLLIST( V, ASGNINDX, 3, VIDXMULT, NVCMULT, PNAMMULT )
                IF ( .NOT. OFLAG(4) ) THEN
                    CALL OPENCTMP( PKTTYP, LDEV )
                    OFLAG(4) = .TRUE.
                END IF
                CALL WRCTMP( LDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                LCTMP = ( LCTMP .OR. CCHEK )

                SFLAG = .TRUE.

              CASE( 'MACT' )

                VIDXMULT = 0     ! array
                CALL ASGNCNTL( NSRC, WDEV, PKTTYP, EANAM( V ), DSFLAG, ASGNINDX )

                CALL UPDATE_POLLIST( V, ASGNINDX, 4, VIDXMULT, NVCMULT, PNAMMULT )

                IF ( .NOT. OFLAG(5) ) THEN
                    CALL OPENCTMP( PKTTYP, MDEV )
                    OFLAG(5) = .TRUE.
                END IF
                CALL WRCTMP( MDEV, V, ASGNINDX, VIDXMULT, CCHEK )
                LCTMP = ( LCTMP .OR. CCHEK )

                SFLAG = .TRUE.

            END SELECT

        END DO       ! End loop on pollutants and activities

    END IF           ! End select on pol-specific packet or not

    !.....   Rewind tmp files
    IF( PDEV .GT. 0 ) REWIND( PDEV )
    IF( CDEV .GT. 0 ) REWIND( CDEV )
    IF( GDEV .GT. 0 ) REWIND( GDEV )
    IF( LDEV .GT. 0 ) REWIND( LDEV )
    IF( MDEV .GT. 0 ) REWIND( MDEV )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.....   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.....   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.....  This internal subprogram flags pollutants that have any
    !       controls applied.
    SUBROUTINE UPDATE_POLLIST( POLID, IDX, CFLG, VIDX, NPCNT, PNAMES )

        !.....  Subprogram arguments
        INTEGER     , INTENT (IN) :: POLID                    ! pollutant index
        INTEGER     , INTENT (IN) :: IDX(NSRC)          ! index to data tables
        INTEGER     , INTENT (IN) :: CFLG                    ! control flag
        INTEGER , INTENT (IN OUT) :: VIDX(NIPPA)             ! pollutant flags
        INTEGER , INTENT (IN OUT) :: NPCNT                   ! pollutant/act count
        CHARACTER(*), INTENT (IN OUT) :: PNAMES(NIPPA)           ! pollutant/act names

        !.....  Local variables
        INTEGER   I, S           ! counters and indices
        LOGICAL   ::      SRCLOOP= .TRUE.           ! true: pollutant 'I' does not
                ! have controls applied

        !----------------------------------------------------------------------

        !.....  For current pollutant, loop through sources until a source with
        !       controls is encountered. Terminate loop when all sources have
        !       been examined.
        SRCLOOP = .TRUE.
        S = 0
        DO WHILE( SRCLOOP .AND. S .LT. NSRC )

            S = S + 1
            IF ( IDX(S) .GT. 0 ) SRCLOOP = .FALSE.          ! controls encountered,
                    ! exit source loop
        END DO          ! end source loop

        !.....  Check to see if current pollutant has controls applied
        IF ( SRCLOOP ) THEN               ! no controls
            VIDX( POLID ) = 0

        ELSE                                  ! controls

            !......  See if pollutant is already in list
            I = INDEX1( EANAM( POLID ), NPCNT, PNAMES )

            !...... If not, add it
            IF( I .LE. 0 ) THEN
                NPCNT  = NPCNT + 1
                PNAMES( NPCNT ) = EANAM( POLID )
                I = NPCNT
            END IF

            !......  Update pol/act flag for this packet
            VIDX( POLID )  = 1

            !......  Update control type flag if control type is set
            IF( CFLG .GT. 0 ) PCTLFLAG( I, CFLG ) = .TRUE.

        END IF

        RETURN

    END SUBROUTINE UPDATE_POLLIST

END SUBROUTINE PROCPKTS
