
SUBROUTINE MRGELEV( NSRC, NMAJOR, NPING, KEY1, KEY2, KEY4, CNV, LINIT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine multiplies a source-emissions vector with optionally
    !      a speciation array and multiplicative control array.
    !      The first time this routine is called, a PinG- and elevated-source-
    !      specific set of arrays are allocated for storing and processing the
    !      PinG and elevated emissions.  This routine computes the elevated
    !      emissions for UAM-style processing and the PinG grouped emissions
    !      for CMAQ PinG files.
    !
    !  PRECONDITIONS REQUIRED:
    !      Assumes any source that will be a PinG source is also a major (elevated)
    !      source.
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ????
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !****************************************************************************/
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
    !.....  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: ELEVFLAG, PEMSRC, PSMATX, PRINFO, INLINEFLAG,   &
                        SRCGRPFLAG

    !.....  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: PCUMATX

    !.....  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: ELEVSIDX, PINGGIDX, NGROUP, GRPGID, PGRPEMIS,    &
                       ELEVEMIS, LMAJOR, LPING, GROUPID,                &
                       ELEVGRPID, EMELEVGRP

    IMPLICIT NONE

    !.....  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC        ! number of source
    INTEGER     , INTENT (IN) :: NMAJOR      ! no. elevated sources
    INTEGER     , INTENT (IN) :: NPING       ! no. plume-in-grid sources
    INTEGER     , INTENT (IN) :: KEY1        ! inven emissions index
    INTEGER     , INTENT (IN) :: KEY2        ! mult controls index
    INTEGER     , INTENT (IN) :: KEY4        ! speciation index
    REAL        , INTENT (IN) :: CNV         ! units conversion factor
    LOGICAL     , INTENT (IN) :: LINIT       ! true: initialize ELEVEMIS

    !.....  Other local variables
    INTEGER         I, J, K, L, S       ! counters and indicies
    INTEGER         IDX             ! index to list of counties in grid
    INTEGER         IOS             ! i/o status

    REAL(8)         SUM1            ! sum for GOUT1
    REAL(8)         SUM2            ! sum for GOUT2
    REAL(8)         MULT            ! tmp value with multiplictv controls
    REAL(8)         REAC            ! tmp value with reactivity controls
    REAL(8)         VAL             ! tmp value
    REAL(8)         VMP             ! tmp market penetration value

    LOGICAL, SAVE:: FIRSTIME = .TRUE.      ! true: first time routine called

    CHARACTER(300)  MESG            ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'MRGELEV'     ! program name

    !***********************************************************************
    !   begin body of subroutine MRGELEV

    !.....  For the first time this routine is called, process the plume-in-
    !       grid indicator array to allocate and generate the necessary
    !       indices
    IF( FIRSTIME ) THEN

        !.....  Allocate indices used for processing in this routine
        ALLOCATE( ELEVSIDX( NSRC ),     &
                  PINGGIDX( NSRC ),     &
                  ELEVEMIS( NMAJOR ),   &
                  PGRPEMIS( NGROUP ), STAT=IOS )        ! grouped PinG
        CALL CHECKMEM( IOS, 'ELEVSIDX...PGRPEMIS', PROGNAME )
        ELEVSIDX = 0           ! array
        PINGGIDX = 0           ! array

        !.....  Create indices used for processing in this routine
        I = 0
        J = 0
        DO S = 1, NSRC

            !.....  Set elevated sources index
            IF( ELEVFLAG .AND. LMAJOR( S ) ) THEN
                I = I + 1
                IF( I .GT. NMAJOR ) CYCLE

                !.....  Store index to major count
                ELEVSIDX( S ) = I

            END IF

            !.....  Set elevated sources index
            IF( INLINEFLAG .AND. LMAJOR( S ) .AND. (.NOT. LPING(S)) ) THEN

                K = FIND1( GROUPID( S ), NGROUP, GRPGID )
                IF( K .GT. 0 ) ELEVSIDX( S ) = K
            END IF

            !.....  Set PinG indicies
            IF( LPING( S ) ) THEN

                !.....  Find group ID in list
                K = FIND1( GROUPID( S ), NGROUP, GRPGID )

                !.....  Store index to groups
                IF( K .GT. 0 ) PINGGIDX( S ) = K

            END IF

        END DO

        !.....  Abort if dimensions incorrect
        IF( ELEVFLAG .AND. I .NE. NMAJOR ) THEN
            WRITE( MESG,94010 ) 'INTERNAL ERROR: Number of elevated sources I=', I,    &
                    'unequal to dimension NMAJOR=', NMAJOR
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF



        FIRSTIME = .FALSE.

    END IF      ! end of firstime section

    !.....  Initialize emissions values if calling routine indicates
    !       that this is a new species being processed.
    IF ( LINIT ) THEN
        PGRPEMIS = 0      ! array
        ELEVEMIS = 0      ! array
        IF( SRCGRPFLAG ) EMELEVGRP = 0.
    END IF

    !.....  Check if this is a valid inventory pollutant for this call
    IF( KEY1 .GT. 0 ) THEN

        !..... If multiplicative controls  speciation
        IF( KEY2 .GT. 0 .AND. KEY4 .GT. 0 ) THEN

            DO S = 1, NSRC

                !.....  Skip if source is not elevated and not PinG
                IF( ELEVSIDX( S ) .EQ. 0 .AND.    &
                    PINGGIDX( S ) .EQ. 0       ) CYCLE

                VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 )                 ! Spec value
                MULT = VAL * PCUMATX( S,KEY2 )                          ! w/ control

                VMP  = PRINFO( S,2 )
                VAL = ( MULT * (1.-VMP) + PRINFO( S,1 ) * VMP )                 ! w/ reactivity control

                !....  NOTE - apply units conversion after application
                !       of reactivity matrix.
                IDX = PINGGIDX( S )
                IF( IDX .GT. 0 )  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                IDX = ELEVSIDX( S )
                IF( IDX .GT. 0 )  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                IF( SRCGRPFLAG ) THEN
                    IDX = ELEVGRPID( S )
                    EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                END IF

            END DO

        !..... If multiplicative controls only
        ELSE IF( KEY2 .GT. 0 ) THEN

            DO S = 1, NSRC

                IF( ELEVSIDX( S ) .EQ. 0 .AND.    &
                    PINGGIDX( S ) .EQ. 0       ) CYCLE

                VAL = PEMSRC( S,KEY1 ) * PCUMATX( S,KEY2 )

                IDX = PINGGIDX( S )
                IF( IDX .GT. 0 )  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                IDX = ELEVSIDX( S )
                IF( IDX .GT. 0 )  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                IF( SRCGRPFLAG ) THEN
                    IDX = ELEVGRPID( S )
                    EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                END IF

            END DO

        !.....  If speciation only
        ELSE IF( KEY4 .GT. 0 ) THEN

            DO S = 1, NSRC

                IF( ELEVSIDX( S ) .EQ. 0 .AND.    &
                    PINGGIDX( S ) .EQ. 0       ) CYCLE

                VAL  = PEMSRC ( S,KEY1 ) * PSMATX( S,KEY4 )
                VMP  = PRINFO( S,2 )
                VAL = ( VAL * (1.-VMP) + PRINFO( S,1 ) * VMP )

                IDX = PINGGIDX( S )
                IF( IDX .GT. 0 )  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                IDX = ELEVSIDX( S )
                IF( IDX .GT. 0 )  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                IF( SRCGRPFLAG ) THEN
                    IDX = ELEVGRPID( S )
                    EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                END IF

            END DO

        !.....  If inventory pollutant only
        ELSE

            DO S = 1, NSRC

                IF( ELEVSIDX( S ) .EQ. 0 .AND.    &
                    PINGGIDX( S ) .EQ. 0       ) CYCLE

                VAL = PEMSRC( S,KEY1 )

                IDX = PINGGIDX( S )
                IF( IDX .GT. 0 )  PGRPEMIS( IDX ) = PGRPEMIS( IDX ) + VAL * CNV

                IDX = ELEVSIDX( S )
                IF( IDX .GT. 0 )  ELEVEMIS( IDX ) = ELEVEMIS( IDX ) + VAL * CNV

                IF( SRCGRPFLAG ) THEN
                    IDX = ELEVGRPID( S )
                    EMELEVGRP( IDX ) = EMELEVGRP( IDX ) + VAL * CNV
                END IF

            END DO

        END IF      ! End which of controls and speciation

    END IF      ! End if no inventory emissions

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

94000 FORMAT( A )

94010 FORMAT( 10 ( A, :, I8, :, 2X ) )

END SUBROUTINE MRGELEV
