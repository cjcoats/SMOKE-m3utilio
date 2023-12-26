
SUBROUTINE SELECTSRC( RCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!    The SELECTSRC routine is responsible for selecting the sources from the
!    inventory records based on the settings for a given report.  If no
!    groups are used in the current report, then all sources are selected.
!    Selected sources will have INDEXA( S ) = 1, while unselected sources will
!    have INDEXA( S ) = 0.
!
!  PRECONDITIONS REQUIRED:
!    The inventory data are already required to have been read in prior to
!    this routine being called
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 4/2001 by M Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: CIFIP, SPPNLO, SPPNHI, SPPROF, SPFRAC, INDEXA

!.........  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: RPT_, LREGION, NREGNGRP, REGNNAM,&
    &                    PINGOUT3, ELEVOUT3, NOELOUT3, PSFLAG,&
    &                    ALLRPT, SLFLAG, SSFLAG, NREGREC, EXCLDRGN

!.........  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: NOUTREC, NSRCDROP, OUTSRC, OUTBIN, OUTSPRO, OUTSFAC

!.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: LMAJOR, LPING

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC

    USE MODSPRO, ONLY: MXSPEC

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: RCNT    ! current report number

!...........  Local variables
    INTEGER         J, L, N, S      ! counters and indices
    INTEGER         IOS             ! i/o status
    INTEGER         REGNIDX         ! index to list of region groups for this rpt

    LOGICAL      :: EFLAG = .FALSE.  ! True: error has been detected

    CHARACTER(FIPLEN3) CFIP                ! tmp country/state/county code
    CHARACTER(256)     MESG                ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'SELECTSRC' ! program name

!***********************************************************************
!   begin body of subroutine SELECTSRC

!.........  Set report-specific local settings
    RPT_    = ALLRPT( RCNT )
    LREGION = .FALSE.

!.........  Determine current report has any groups
    REGNIDX = 0
    IF( RPT_%REGNNAM .NE. ' ' ) THEN
        REGNIDX = INDEX1( RPT_%REGNNAM, NREGNGRP, REGNNAM( : ) )
        LREGION = .TRUE.

        IF( REGNIDX .LE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Group "'// TRIM( RPT_%REGNNAM )//&
            &   '"' // CRLF() // BLANK10 // 'not found in list '//&
            &   'of groups.'
            CALL M3MSG2( MESG )
        END IF

    END IF

!.........  NOTE - if adding group types, do it here and in errors below

!.........  Report internal errors if group name is not found
    IF ( EFLAG ) THEN
        MESG = 'Problems using groups'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Deallocate arrays if they've been allocated before
    ALLOCATE( INDEXA( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'INDEXA', PROGNAME )

!......... Initialize list of selected sources
    DO S = 1, NSRC
        INDEXA( S ) = 1
    END DO

!.........  If need to select sources, loop through sources to identify
!.........  the ones not to output
    IF( LREGION .OR. RPT_%ELEVSTAT .GT. 0 ) THEN

        DO S = 1, NSRC

!.................  If using a region group, search for FIPS code in list
            IF( LREGION ) THEN
                CFIP = CIFIP( S )
                J = FINDC( CFIP, NREGREC( REGNIDX ),&
                &           EXCLDRGN ( 1,REGNIDX )   )

                IF ( J .GT. 0 ) INDEXA( S ) = 0

            END IF

!.................  If selecting elevated srcs...
            SELECT CASE( RPT_%ELEVSTAT )
              CASE( PINGOUT3 )
                IF ( .NOT. LPING( S ) ) INDEXA( S ) = 0

              CASE( ELEVOUT3 )
                IF ( .NOT. LMAJOR( S ) ) INDEXA( S ) = 0

              CASE( NOELOUT3 )
                IF ( LMAJOR( S ) ) INDEXA( S ) = 0

            END SELECT

!.................  NOTE - this should always be at end of select statements
!.................  Keep a count of the number of selected sources
            IF( INDEXA( S ) .EQ. 0 ) NSRCDROP = NSRCDROP + 1

        END DO   ! End loop through sources

    END IF

!.........  Now create compressed list of sources, but leave INDEXA as is
!           for use by REPMRGGRD

    IF( PSFLAG ) THEN
        NOUTREC = 0
        DO S = 1, NSRC

            IF ( INDEXA( S ) .EQ. 1 ) THEN
                NOUTREC = NOUTREC + SPPNHI( S ) - SPPNLO( S ) + 1
            END IF

        END DO

    ELSE

        NOUTREC = NSRC - NSRCDROP

    END IF

    IF ( ALLOCATED( OUTSRC  ) )  DEALLOCATE( OUTSRC )
    IF ( ALLOCATED( OUTBIN  ) )  DEALLOCATE( OUTBIN )
    IF ( ALLOCATED( OUTSPRO ) )  DEALLOCATE( OUTSPRO )
    IF ( ALLOCATED( OUTSFAC ) )  DEALLOCATE( OUTSFAC )

    ALLOCATE( OUTSRC( NOUTREC ),&
    &          OUTBIN( NOUTREC ),&
    &         OUTSPRO( NOUTREC ),&
    &         OUTSFAC( NOUTREC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'OUTSRC,OUTSFAC', PROGNAME )

!.........  Now create compressed list of sources, but leave INDEXA as is
!           for use by REPMRGGRD

    IF( PSFLAG ) THEN

        J = 0
        DO S = 1, NSRC

            IF( INDEXA( S ) .EQ. 0 ) CYCLE

            DO N = SPPNLO( S ), SPPNHI( S )
                J = J + 1
                OUTSRC(  J ) = S
                OUTSPRO( J ) = SPPROF( N )
                OUTSFAC( J ) = SPFRAC( N )
            END DO

        END DO

    ELSE

        J = 0
        DO S = 1, NSRC

            IF( INDEXA( S ) .EQ. 0 ) CYCLE

            J = J + 1
            OUTSRC(  J ) = S
            OUTSPRO( J ) = 'NA'
            OUTSFAC( J ) = 1.0

        END DO

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE SELECTSRC
