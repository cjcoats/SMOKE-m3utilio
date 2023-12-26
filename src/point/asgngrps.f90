
SUBROUTINE ASGNGRPS( NSP, NCRIT, MXCHK, CRITVALS, CRITYPES, NINVGRP )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     This routine assigns groups to the inventory sources based on comparison
    !     criteria provided as subroutine arguments.
    !
    !  PRECONDITIONS REQUIRED:
    !     Comparison criteria must be in a specific format in which the "ORs" are
    !     on separate rows and the "ANDs" are on the same row.  Modules MODINFO
    !     and MODSOURC must be populated, and MODELEV unallocated.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !     EVALCRIT is responsible for testing the formulas
    !
    !  REVISION  HISTORY:
    !       Written 7/2001 by M. Houyoux
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
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: CSOURC, CIFIP, STKHT, STKDM, STKTK, STKVE,    &
                        XLOCA, YLOCA

    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: GINDEX, GROUPID, GRPGID, GRPLAT, GRPLON,    &
                       GRPDM, GRPHT, GRPTK, GRPVE, GRPFL, GRPCNT,GRPFIP

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NCHARS, SC_BEGP, SC_ENDP

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   INCLUDES:
    INCLUDE 'EMCNST3.h90'       ! emissions constant parameters
    INCLUDE 'CONST3.EXT'        ! physical and mathematical constants

    !.......   ARGUMENTS and their descriptions:
    INTEGER     , INTENT (IN) :: NSP         ! no. of variables in formulas
    INTEGER     , INTENT (IN) :: NCRIT       ! no. of comparison criteria
    INTEGER     , INTENT (IN) :: MXCHK       ! max. no. checks per NSP and NCRIT
    REAL        , INTENT (IN) :: CRITVALS( NCRIT, MXCHK, NSP )     ! formula values
    CHARACTER(*), INTENT (IN) :: CRITYPES( NCRIT, MXCHK, NSP )     ! formula types
    INTEGER     , INTENT(OUT) :: NINVGRP     ! no. of groups

    !.......  .h90ERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: EVALCRIT

    !.......   LOCAL PARAMETERS and their descriptions:
    INTEGER,       PARAMETER :: MXLOCGRP = 1000      ! Max number groups per facility
    INTEGER,       PARAMETER :: MXSPGRP  = 1000      ! Max number sources per group
    CHARACTER(16), PARAMETER :: PROGNAME = 'ASGNGRPS'       !  program name

    !.......   Local allocatable arrays
    INTEGER             TGRPCNT( NSRC )     ! tmp no. sources in group
    REAL                TGRPHT ( NSRC )     ! tmp stack heights [m]
    REAL                TGRPDM ( NSRC )     ! tmp stack diameters [m]
    REAL                TGRPTK ( NSRC )     ! tmp stack exit temperatures [K]
    REAL                TGRPVE ( NSRC )     ! tmp stack exit velocities [m/s]
    REAL                TGRPFL ( NSRC )     ! tmp stack flow rates [m^3/s]

    REAL                VALS   ( NSP )     ! tmp source stack parameter array
    REAL                REFS   ( NSP )     ! tmp group stack parameter array
    CHARACTER(PLTLEN3)  CHRS( NSP )            ! dummy test character strings

    LOGICAL             GPSTAT( NSP, NCRIT, MXCHK )    ! not used except EVALCRIT call
    CHARACTER(PLTLEN3)  COMCHRS( NSP, NCRIT, MXCHK )    ! dummy formula strings

    !.......   Local fixed arrays
    INTEGER     GSRC  ( MXLOCGRP, MXSPGRP )     ! source IDs for local groups
    INTEGER     GCNT  ( MXLOCGRP )      ! local group count
    INTEGER     LOCGID( MXLOCGRP )      ! local group ID
    REAL        G_DM  ( MXLOCGRP )      ! local group diameter
    REAL        G_FL  ( MXLOCGRP )      ! local group exit flow rate
    REAL        G_HT  ( MXLOCGRP )      ! local group stack height
    REAL        G_TK  ( MXLOCGRP )      ! local group exit temperature
    REAL        G_VE  ( MXLOCGRP )      ! local group exit velocity

    !.......   OTHER LOCAL VARIABLES and their descriptions:

    INTEGER     G, I, J, L2, S, S2          ! counters and indices

    INTEGER     GLM                  ! local G max value
    INTEGER     IOS                  ! i/o status
    INTEGER     LOCG                 ! local G
    INTEGER  :: MXGCNT = 0               ! max of all GCNT entries
    INTEGER  :: NLOCGRP = 1          ! number of local (facility) groups
    INTEGER     PRVG                 ! G for previous interation
    INTEGER, SAVE :: NWARN = 0       ! warning count
    INTEGER, SAVE :: MXWARN          ! max no. warnings

    REAL        DM                   ! tmp stack diameter [m]
    REAL        FL                   ! tmp stack exit flow [m]
    REAL        HT                   ! tmp stack height [m]
    REAL        TK                   ! tmp stack exit temperature [K]
    REAL        VE                   ! tmp stack exit velocity [m/s]
    REAL        W1, W2               ! tmp multipler weights

    LOGICAL :: LANYGRP = .FALSE.      ! true: one or more groups in inventory
    LOGICAL :: LGROUP  = .FALSE.      ! true: group found for this source
    LOGICAL :: LGRPALL = .FALSE.      ! true: group found for previous facility
    LOGICAL :: STATUS  = .FALSE.      ! true: tmp evaluation status

    CHARACTER(FIPLEN3) CFIP          ! tmp FIPS code
    CHARACTER(FIPLEN3) PFIP          ! previous FIPS code
    CHARACTER(300)  BUFFER           ! src description buffers
    CHARACTER(400)  MESG             ! msg buffer

    CHARACTER(CHRLEN3) PLT       ! tmp facility ID
    CHARACTER(CHRLEN3) PPLT      ! tmp facility ID

    !***********************************************************************
    !   begin body of subroutine ASGNGRPS

    !.......  Write status message
    MESG = 'Assigning stack groups...'
    CALL M3MSG2( MESG )

    !.......   Get maximum number of warnings
    MXWARN = ENVINT( WARNSET , ' ', 100, I )

    TGRPCNT = 0
    TGRPHT = 0.
    TGRPDM = 0.
    TGRPTK = 0.
    TGRPVE = 0.
    TGRPFL = 0.
    CHRS    = ' '       ! array
    COMCHRS = ' '       ! array

    !.......  Loop through sources to establish groups
    G    = 0
    PFIP = ' '
    PPLT = ' '
    DO S = 1, NSRC

        CFIP= CIFIP ( S )
        PLT = CSOURC( S )( SC_BEGP( 2 ):SC_ENDP( 2 ) )
        HT  = STKHT ( S )
        DM  = STKDM ( S )
        TK  = STKTK ( S )
        VE  = STKVE ( S )
        FL  = 0.25 * PI * DM * DM * VE

        GINDEX( S ) = S           ! Set sorting index

        !.......  For a new facility, assume that a group will form, but store
        !               source ID to reset group if none materialized.
        IF( CFIP .NE. PFIP .OR. PLT .NE. PPLT ) THEN

            !.......  If there was one or more groups for the previous facility,
            !                   compute the group numbers and store info by source
            IF( LGRPALL ) THEN

                GLM = 0
                DO I = 1, NLOCGRP

                    !...........  If there are multiple stacks in the local group,
                    !                           then it's a global group.  Update global arrays.
                    IF ( GCNT( I ) .GT. 1 ) THEN

                        DO J = 1, GCNT( I )

                            S2 = GSRC( I,J )
                            GROUPID( S2 ) = G + LOCGID( I )
                            TGRPCNT( S2 ) = GCNT( I )
                            TGRPHT ( S2 ) = G_HT( I )
                            TGRPDM ( S2 ) = G_DM( I )
                            TGRPTK ( S2 ) = G_TK( I )
                            TGRPVE ( S2 ) = G_VE( I )
                            TGRPFL ( S2 ) = G_FL( I )

                            GLM = MAX( GLM, GROUPID( S2 ) )              ! Max for facility

                        END DO

                    END IF

                END DO

                !.......  Update global counter. When get to the end, we'll have
                !                       to reset the groups so that there aren't any holes.
                G = GLM

            END IF

            !.......  Initialize info for possible current group
            PFIP    = CFIP
            PPLT    = PLT
            LGRPALL = .FALSE.
            NLOCGRP = 1
            GSRC  ( 1,1 ) = S
            GCNT  ( 1 ) = 1
            LOCGID( 1 ) = 1
            G_HT  ( 1 ) = HT
            G_DM  ( 1 ) = DM
            G_TK  ( 1 ) = TK
            G_VE  ( 1 ) = VE
            G_FL  ( 1 ) = FL

        !.......  For the same facility, compare stack parameters to the
        !               group stack parameters using the tolerances.
        !.......  If there is a match in stack parameters, set flag current source
        !               with the group number,
        ELSE

            VALS( HT_IDX ) = HT
            VALS( DM_IDX ) = DM
            VALS( TK_IDX ) = TK
            VALS( VE_IDX ) = VE
            VALS( FL_IDX ) = FL

            !.......  Compare this source to all other groups/sources for this
            !                   plant
            LGROUP = .FALSE.
            DO I = 1, NLOCGRP

                IF( I .LE. MXLOCGRP ) THEN
                    REFS( HT_IDX ) = G_HT( I )
                    REFS( DM_IDX ) = G_DM( I )
                    REFS( TK_IDX ) = G_TK( I )
                    REFS( VE_IDX ) = G_VE( I )
                    REFS( FL_IDX ) = G_FL( I )
                END IF

                !.......  Check tolerances. Use REFS for RANK field since RANK
                !                       can't be used to group stacks (it makes no sense)
                STATUS = EVALCRIT( NSP, NCRIT, MXCHK, VALS, REFS,    &
                                   REFS, CHRS, CRITVALS, COMCHRS,    &
                                   CRITYPES, GPSTAT)

                !.......  If stack parameters meet the criteria, then recompute
                !                       the group stack parameters as the weighted average.
                !.......  Exit from loop if a match has been found.
                IF( STATUS ) THEN

                    GCNT( I ) = GCNT( I ) + 1
                    W1 = REAL( GCNT( I ) - 1 ) / REAL( GCNT( I ) )
                    W2 = 1. / REAL( GCNT( I ) )
                    G_HT( I ) = G_HT( I ) * W1 + HT * W2
                    G_DM( I ) = G_DM( I ) * W1 + DM * W2
                    G_TK( I ) = G_TK( I ) * W1 + TK * W2
                    G_VE( I ) = G_VE( I ) * W1 + VE * W2
                    G_FL( I ) = G_FL( I ) * W1 + FL * W2

                    IF( GCNT( I ) .LE. MXSPGRP ) THEN
                        GSRC( I,GCNT( I ) ) = S
                    END IF

                    MXGCNT = MAX( MXGCNT,GCNT( I ) )

                    LGROUP  = .TRUE.
                    LGRPALL = .TRUE.
                    LANYGRP = .TRUE.
                    EXIT

                END IF          ! If stack parms are meeting group criteria

            END DO          ! loop through temporary groups at the plant

            !.......  If a match was not found for this source, add its
            !         stack parameters to the temporary group list
            IF( .NOT. LGROUP ) THEN
                NLOCGRP = NLOCGRP + 1

                IF ( NLOCGRP .LE. MXLOCGRP ) THEN
                    G_HT  ( NLOCGRP ) = HT
                    G_DM  ( NLOCGRP ) = DM
                    G_TK  ( NLOCGRP ) = TK
                    G_VE  ( NLOCGRP ) = VE
                    G_FL  ( NLOCGRP ) = FL
                    GCNT  ( NLOCGRP ) = 1
                    LOCGID( NLOCGRP ) = NLOCGRP

                    IF( GCNT( NLOCGRP ) .LE. MXSPGRP ) THEN
                        GSRC( NLOCGRP,GCNT( NLOCGRP ) ) = S
                    END IF

                    MXGCNT = MAX( MXGCNT,GCNT( NLOCGRP ) )

                END IF

            END IF

        END IF          ! If same facility or not

    END DO              ! loop through sources

    !.......  Ensure that we didn't run out of space for the per-facility groups
    IF ( NLOCGRP .GT. MXLOCGRP ) THEN

        WRITE( MESG,94010 ) 'INTERNAL ERROR: The maximum number ' //    &
               'of groups per plant (', MXLOCGRP, ') was ' //    &
               'exceeded with a number of ', NLOCGRP
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Ensure that we had enough space for the sources per group
    IF( MXGCNT .GT. MXSPGRP ) THEN
        WRITE( MESG,94010 ) 'INTERNAL ERROR: The maximum number ' //    &
               'of sources per group (', MXSPGRP, ') was ' //    &
               'exceeded with a number of ', MXGCNT
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.......  Postprocess group IDs to remove gaps in number and to count the
    !           number of inventory groups.  If LANYGRP is false, then GROUPID
    !           will be zeros, so do not want to sort.  GINDEX will be maintained
    !           as the source number.
    IF( LANYGRP ) THEN
        CALL SORTI1( NSRC, GINDEX, GROUPID )
    ELSE
        MESG = 'NOTE: No stack groups assigned.'
        CALL M3MSG2( MESG )
    END IF

    PRVG = 0
    G    = 0
    DO S = 1, NSRC

        I    = GINDEX( S )
        LOCG = GROUPID( I )

        IF ( LOCG .NE. 0 ) THEN

            !.......  If a new group, increment group ID
            IF ( LOCG .NE. PRVG ) G = G + 1

            !.......  Store global group number
            GROUPID( I ) = G
            PRVG = LOCG

        END IF

    END DO

    !.......  Store the actual number of groups present
    NINVGRP = G

    !.......  Allocate memory for group information and populate group arrays
    ALLOCATE( GRPGID( NINVGRP ),    &
              GRPLAT( NINVGRP ),    &
              GRPLON( NINVGRP ),    &
               GRPDM( NINVGRP ),    &
               GRPHT( NINVGRP ),    &
               GRPTK( NINVGRP ),    &
               GRPVE( NINVGRP ),    &
               GRPFL( NINVGRP ),    &
              GRPCNT( NINVGRP ),    &
              GRPFIP( NINVGRP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GRPGID...GRPFIP', PROGNAME )

    GRPGID  = 0
    GRPLAT  = BADVAL3
    GRPLON  = BADVAL3
    GRPDM   = BADVAL3
    GRPHT   = BADVAL3
    GRPTK   = BADVAL3
    GRPVE   = BADVAL3
    GRPFL   = BADVAL3
    GRPCNT  = 0
    GRPFIP  = ' '

    !.......  Store the group information more succinctly
    PRVG = 0
    DO S = 1, NSRC
        I = GINDEX( S )
        G = GROUPID( I )

        !.......  Skip sources that aren't in a group
        IF ( G .EQ. 0 ) CYCLE

        !.......  Store. Location will be from first source encountered in group
        IF ( G .NE. PRVG ) THEN
            GRPGID( G ) = G
            GRPLAT( G ) = YLOCA  ( I )
            GRPLON( G ) = XLOCA  ( I )
            GRPDM ( G ) = TGRPDM ( I )
            GRPHT ( G ) = TGRPHT ( I )
            GRPTK ( G ) = TGRPTK ( I )
            GRPVE ( G ) = TGRPVE ( I )
            GRPFL ( G ) = TGRPFL ( I )
            GRPCNT( G ) = TGRPCNT( I )
            GRPFIP( G ) = CIFIP  ( I )
        !.......  Otherwise, give warnings if stacks in the same group do not have
        !               the same stack locations.
        ELSE

            CALL FMTCSRC( CSOURC( I ), NCHARS, BUFFER, L2 )

            NWARN = NWARN + 1

            IF ( FLTERR( GRPLAT( G ), YLOCA( I ) ) ) THEN
                WRITE( MESG,94078 ) 'WARNING: Source latitude',    &
                       YLOCA( I ), 'is inconsistent with group' //    &
                       CRLF() // BLANK10 // 'latitude', GRPLAT(G),    &
                       'taken from the first source in the group '    &
                       // 'for source:' // CRLF() // BLANK10 //    &
                       BUFFER( 1:L2 )
                IF( NWARN <= MXWARN ) CALL M3MESG( MESG )
            END IF

            IF ( FLTERR( GRPLON( G ), XLOCA( I ) ) ) THEN
                WRITE( MESG,94078 ) 'WARNING: Source longitude',    &
                       XLOCA( I ), 'is inconsistent with group' //    &
                       CRLF() // BLANK10 // 'longitude', GRPLON(G),    &
                       'taken from the first source in the group '    &
                       // 'for source:' // CRLF() // BLANK10 //    &
                       BUFFER( 1:L2 )
                IF( NWARN <= MXWARN )CALL M3MESG( MESG )
            END IF

        END IF

        PRVG = G

    END DO

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

94078 FORMAT( A, 1X, F8.2, 1X, A, 1X, F8.2, 1X, A )

END SUBROUTINE ASGNGRPS

