
SUBROUTINE APPLREAC( NSRC, NREAC, KEY1, KEY2, PRJFLAG, LMKTPON,&
&                     EMIS, EMIST, IDX, REPEM, PRJFC, MKTPEN,&
&                     RMAT, RCTINFO )

!***********************************************************************
!  subroutine body starts at line 59
!
!  DESCRIPTION:
!      This subroutine applies the reactivity matrix for the species of
!      interest (indicated by the KEY, which has already been checked to
!      ensure species is in the reactivity matrix).  The replacement
!      emissions are applied to the sources in the array of inventory
!      emissions.  The base emissions from the reactivity matrix are used
!      to generate a temporal factor from the ratio of base emissions to
!      the inventory emissions, which could be temporalized.  If the inventory
!      emissions are already projected, the projection factor from the
!      reactivity matrix is applied.
!
!      It also populates a market penetration array for all sources.  The
!      default value is 0 for sources that don't get reactivity controls.  If
!      the flag is set to ignore market penetration, then the penetration for
!      reactivity sources is 1.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!     ??/???? by ?????
!
!****************************************************************************/
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

    IMPLICIT NONE

!.........  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC             ! number inventory sources
    INTEGER     , INTENT (IN) :: NREAC            ! no. in react matrix
    INTEGER     , INTENT (IN) :: KEY1             ! position in EMIS
    INTEGER     , INTENT (IN) :: KEY2             ! position in RMAT
    LOGICAL     , INTENT (IN) :: PRJFLAG          ! true: emissions are proj
    LOGICAL     , INTENT (IN) :: LMKTPON          ! true: use mkt pen
    REAL        , INTENT (IN) :: EMIS   ( NSRC,* )! inventory emissions
    REAL        , INTENT (IN) :: EMIST  ( NSRC,* )! inv or hourly emissions
    INTEGER     , INTENT (IN) :: IDX    ( NREAC ) ! src idx of react. matrix
    REAL        , INTENT (IN) :: REPEM  ( NREAC ) ! replacement emissions
    REAL        , INTENT (IN) :: PRJFC  ( NREAC ) ! stored orig. speciation
    REAL        , INTENT (IN) :: MKTPEN ( NREAC ) ! market penetration
    REAL        , INTENT (IN) :: RMAT   ( NREAC,* ) ! reactivity facs
    REAL        , INTENT(OUT) :: RCTINFO( NSRC,2 )! per-source react. info

!.........  Other local variables
    INTEGER         J, S        !  counters and indices

    REAL            EM          !  tmp emission value
    REAL            PFAC        !  tmp projection factor

    CHARACTER(16), PARAMETER :: PROGNAME = 'APPLREAC' ! program name

!***********************************************************************
!   begin body of subroutine APPLREAC

!.........  Initialize reactivity arrays for all sources
    RCTINFO( :,1 ) = EMIS( :,KEY1 )
    RCTINFO( :,2 ) = 0.

!.........  Exit from subroutine if reactivity matrix does not apply for this
!           species
    IF( KEY2 .LE. 0 ) RETURN

!.........  Note that inventory and hourly (if present) have already been
!           checked to insure the same base and projection years.
!.........  Compute the reactivity-controlled emissions by source and set the
!           market penetration.
    PFAC = 1.
    DO J = 1, NREAC

!.............  Retrieve the source index from the reactivity matrix
        S = IDX( J )

!.............  If the inventory emissions are already projected, apply the
!               reactivity projection factor
        IF( PRJFLAG ) PFAC = PRJFC( J )

!.............  Compute output emissions as the product of the replacement
!               emissions, times the projection factor, times the temporal
!               adjustment (or 1. if EMIS=EMIST). Screen for divide by zero.
        EM = 0.
        IF( EMIS( S,KEY1 ) .GT. 0 ) THEN

            EM = REPEM( J ) * PFAC *&
            &     EMIST( S,KEY1 ) / EMIS( S,KEY1 )

        END IF

        RCTINFO( S,1 ) = EM * RMAT( J,KEY2 )

!.............  If using market penetration, then set to value in matrix
        IF( LMKTPON ) THEN

            RCTINFO( S,2 ) = MKTPEN( J )

!.............  If ignoring market penetration, then set to unity
        ELSE

            RCTINFO( S,2 ) = 1.

        END IF

    END DO

    RETURN

END SUBROUTINE APPLREAC
