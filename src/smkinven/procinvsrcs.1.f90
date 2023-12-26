
SUBROUTINE PROCINVSRCS( NRAWSRCS )

!**************************************************************************
!  subroutine body starts at line 114
!
!  DESCRIPTION:
!      This subroutine sorts and stores the unique source information.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!      Created 2/03 by C. Seppanen based on procinven.f
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

!...........   MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: SRCIDA, CSOURCA, CSOURC, CIFIP,&
    &                    CSCC, XLOCA, YLOCA, CELLID, IRCLAS,&
    &                    IVTYPE, CLINK, CVTYPE

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NSRC

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constat parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER , INTENT (IN) :: NRAWSRCS ! no. raw srcs

!...........   Other local variables
    INTEGER         I, J, S            ! counter and indices
    INTEGER         IOS                ! I/O status

    CHARACTER(ALLLEN3) TSRC        !  tmp source information

    CHARACTER(16) :: PROGNAME = 'PROCINVSRCS' ! program name

!***********************************************************************
!   begin body of subroutine PROCINVSRCS

!.........  Allocate memory for sorted inventory arrays
    ALLOCATE( CIFIP( NSRC ),&
    &           CSCC( NSRC ),&
    &         CSOURC( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSOURC', PROGNAME )

    SELECT CASE( CATEGORY )
      CASE( 'AREA' )
        ALLOCATE( XLOCA( NSRC ),&
        &          YLOCA( NSRC ),&
        &         CELLID( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CELLID', PROGNAME )
      CASE( 'MOBILE' )
        ALLOCATE( IRCLAS( NSRC ),&
        &          IVTYPE( NSRC ),&
        &           CLINK( NSRC ),&
        &          CVTYPE( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'CVTYPE', PROGNAME )
      CASE( 'POINT' )
    END SELECT

!.........  Loop through sources to store sorted arrays
!           for output to I/O API file.
!.........  Keep case statement outside the loops to speed processing
    SELECT CASE ( CATEGORY )
      CASE( 'AREA' )

        DO I = 1, NRAWSRCS

            S = SRCIDA( I )
            TSRC = CSOURCA( I )

            CIFIP( S ) = TSRC( 1:FIPLEN3 )
            CSCC( S )  = TSRC( SCCPOS3:SCCPOS3+SCCLEN3-1 )
            CSOURC( S ) = TSRC

        END DO

        XLOCA = BADVAL3   ! array
        YLOCA = BADVAL3   ! array
        CELLID = 0        ! array

      CASE( 'MOBILE' )

        DO I = 1, NRAWSRCS

            S = SRCIDA( I )
            TSRC = CSOURCA( I )

            CIFIP( S )  = TSRC( 1:FIPLEN3 )
            IRCLAS( S ) =&
            &    STR2INT( TSRC( RWTPOS3:RWTPOS3+RWTLEN3-1 ) )
            IVTYPE( S ) =&
            &    STR2INT( TSRC( VIDPOS3:VIDPOS3+VIDLEN3-1 ) )
            CSCC( S )   = TSRC( MSCPOS3:MSCPOS3+SCCLEN3-1 )
            CLINK( S )  = TSRC( LNKPOS3:LNKPOS3+LNKLEN3-1 )
            CSOURC( S ) = TSRC
            CVTYPE( S ) = 'MOVES'

        END DO

      CASE( 'POINT' )

        DO I = 1, NRAWSRCS

            S = SRCIDA( I )
            TSRC = CSOURCA( I )

            CIFIP( S ) = TSRC( 1:FIPLEN3 )
            CSCC( S ) = TSRC( CH4POS3:CH4POS3+SCCLEN3-1 )

            CSOURC( S ) = TSRC
        END DO

    END SELECT

!.........  Deallocate per-source unsorted arrays
    DEALLOCATE( CSOURCA )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE PROCINVSRCS
