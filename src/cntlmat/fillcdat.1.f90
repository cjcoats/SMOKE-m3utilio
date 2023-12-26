
SUBROUTINE FILLCDAT( PKTTYP, JT, PKTINFO )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine populates the control data table depending on the
!      packet type being processed when the routine is called (PKTTYP)
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Started 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***************************************************************************
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

!.........  MODULES for public variables
!.........  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: FACCTG, CUTCTG, FACMACT, FACRACT,&
    &                    ICTLEQUP, CCTLSIC, FACCEFF, FACREFF,&
    &                    FACRLPN, CALWSIC, FACALW, EMCAPALW,&
    &                    EMREPALW, EMREPREA, PRJFCREA, MKTPNREA,&
    &                    CSCCREA, CSPFREA, CPRJSIC, PRJFC, CTLRPLC,&
    &                    MACEXEFF, MACNWEFF, MACNWFRC, CMACSRCTYP,&
    &                    CTGCOMT, CTLCOMT, ALWCOMT,&
    &                    REACOMT, PRJCOMT, MACCOMT

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'CPKTDAT.EXT'   !  control packet contents

!...........   SUBROUTINE ARGUMENTS:

    CHARACTER(*), INTENT (IN) :: PKTTYP    ! packet type
    INTEGER     , INTENT (IN) :: JT        ! index to control data tables
    TYPE(CPACKET),INTENT (IN) :: PKTINFO   ! packet information

    CHARACTER(16), PARAMETER :: PROGNAME = 'FILLCDAT' ! program name

!***********************************************************************
!   Begin body of subroutine FILLCDAT

    SELECT CASE( PKTTYP )

      CASE( 'CTG' )
        FACCTG ( JT ) = PKTINFO%FAC1
        CUTCTG ( JT ) = PKTINFO%FAC2
        FACMACT( JT ) = PKTINFO%FAC3
        FACRACT( JT ) = PKTINFO%FAC4
        CTGCOMT( JT ) = PKTINFO%COMMENT

      CASE( 'CONTROL' )
        ICTLEQUP( JT ) = INT    ( PKTINFO%FAC1 )
        CCTLSIC ( JT ) = PKTINFO%CSIC
        FACCEFF ( JT ) = PKTINFO%FAC2
        FACREFF ( JT ) = PKTINFO%FAC3
        FACRLPN ( JT ) = PKTINFO%FAC4
        IF( PKTINFO%REPFLAG == 'R' .OR.&
        &    PKTINFO%REPFLAG == 'Y'      ) THEN
            CTLRPLC( JT ) = .TRUE.
        END IF
        CTLCOMT( JT ) = PKTINFO%COMMENT

      CASE( 'ALLOWABLE' )
        CALWSIC ( JT ) = PKTINFO%CSIC
        FACALW  ( JT ) = PKTINFO%FAC1
        EMCAPALW( JT ) = PKTINFO%FAC2
        EMREPALW( JT ) = PKTINFO%FAC3
        ALWCOMT( JT ) = PKTINFO%COMMENT

      CASE( 'REACTIVITY' )
        EMREPREA( JT ) = PKTINFO%FAC1
        PRJFCREA( JT ) = PKTINFO%FAC2
        MKTPNREA( JT ) = PKTINFO%FAC3
        CSCCREA ( JT ) = PKTINFO%NSCC
        CSPFREA ( JT ) = PKTINFO%TMPPRF
        REACOMT( JT ) = PKTINFO%COMMENT

      CASE( 'PROJECTION' )
        CPRJSIC( JT ) = PKTINFO%CSIC
        PRJFC  ( JT ) = PKTINFO%FAC1
        PRJCOMT( JT ) = PKTINFO%COMMENT

      CASE( 'MACT' )
        MACEXEFF( JT ) = PKTINFO%FAC1
        MACNWEFF( JT ) = PKTINFO%FAC2
        MACNWFRC( JT ) = PKTINFO%FAC3

        CMACSRCTYP( JT ) = PKTINFO%CSTYP
        CALL PADZERO( CMACSRCTYP( JT ) )
        MACCOMT( JT ) = PKTINFO%COMMENT

!.............  Make sure src type is only 00, 01, or 02
        IF( CMACSRCTYP( JT ) /= '01' .AND.&
        &    CMACSRCTYP( JT ) /= '02'       ) THEN
            CMACSRCTYP( JT ) = '00'
        END IF

    END SELECT

    RETURN

END SUBROUTINE FILLCDAT
