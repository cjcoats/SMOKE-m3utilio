
SUBROUTINE FILLSTBL( NIPOL, NXREF, ICSIZE, XTYPE, XTCNT )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine populates the speciation profile codes part of the
!      grouped speciation cross-reference tables.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 3/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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
    USE M3UTILIO

!...........   This module is for cross reference tables
    USE MODXREF, ONLY: CSPT01, CSPT02, CSPT03, CSPT04, CSPT05,&
    &        CSPT06, CSPT07, CSPT08, CSPT09, CSPT10,&
    &        CSPT11, CSPT12, CSPT13, CSPT14, CSPT15, CSPT16,&
    &        CSPT26, CSPT27, CSPT28, CSPT29, CSPT30, CSPT31,&
    &        CSPT32, CSPT33, CSPT34, CSPT35, CSPT36, CSPT37,&
    &        INDXTA, ISPTA, CSPRNA

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NIPOL           ! no. pollutants
    INTEGER     , INTENT (IN) :: NXREF           ! no. ungrpd x-ref entries
    INTEGER     , INTENT (IN) :: ICSIZE( * )     ! size of x-ref groups
    INTEGER     , INTENT (IN) :: XTYPE ( NXREF ) ! group no. of x-ref entry
    INTEGER     , INTENT (IN) :: XTCNT ( NXREF ) ! pos. in x-ref group

!...........   Other local variables
    INTEGER       I, J, K, T       ! counter and indices
    INTEGER       ISP              ! temporary pollutant position in EINAM

    LOGICAL    :: EFLAG = .FALSE.  ! true: error has occurred

    CHARACTER(300)     MESG    ! message buffer
    CHARACTER(SPNLEN3) SPCODE  ! tmp for speciation profile code

    CHARACTER(16) :: PROGNAME = 'FILLSTBL' ! program name

!***********************************************************************
!   begin body of subroutine FILLSTBL

!.........  Store the temporal profile codes for each x-ref entry, depending
!           on the group (XTYPE) and the position in that group (XTCNT)

    DO I = 1, NXREF

        J      = INDXTA( I )
        ISP    = ISPTA ( J )
        SPCODE = CSPRNA( J )

        T      = XTYPE ( I )
        K      = XTCNT ( I )

        IF( ISP .EQ. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'INTERNAL ERROR: Pollutant code was 0 ' //&
            &       'for speciation cross-reference entry ' //&
            &       CRLF()// BLANK10// 'in program ' // PROGNAME
            CALL M3MESG( MESG )
        ENDIF

!.................  Populate tables depending on type. Note that all profiles
!                   are, by definition, pollutant-specific.
        SELECT CASE ( T )

          CASE( 0 )  ! Skip this x-ref because it is invalid or duplicate

          CASE( 1 )
            CSPT01( ISP )   = SPCODE

          CASE( 2 )
            CSPT02( K,ISP ) = SPCODE

          CASE( 3 )
            CSPT03( K,ISP ) = SPCODE

          CASE( 4 )
            CSPT04( K,ISP ) = SPCODE

          CASE( 5 )
            CSPT05( K,ISP ) = SPCODE

          CASE( 6 )
            CSPT06( K,ISP ) = SPCODE

          CASE( 7 )
            CSPT07( K,ISP ) = SPCODE

          CASE( 8 )
            CSPT08( K,ISP ) = SPCODE

          CASE( 9 )
            CSPT09( K,ISP ) = SPCODE

          CASE( 10 )
            CSPT10( K,ISP ) = SPCODE

          CASE( 11 )
            CSPT11( K,ISP ) = SPCODE

          CASE( 12 )
            CSPT12( K,ISP ) = SPCODE

          CASE( 13 )
            CSPT13( K,ISP ) = SPCODE

          CASE( 14 )
            CSPT14( K,ISP ) = SPCODE

          CASE( 15 )
            CSPT15( K,ISP ) = SPCODE

          CASE( 16 )
            CSPT16( K,ISP ) = SPCODE

!.............  SIC specific entries
          CASE( 26 )
            CSPT26( K,ISP ) = SPCODE

          CASE( 27 )
            CSPT27( K,ISP ) = SPCODE

          CASE( 28 )
            CSPT28( K,ISP ) = SPCODE

          CASE( 29 )
            CSPT29( K,ISP ) = SPCODE

          CASE( 30 )
            CSPT30( K,ISP ) = SPCODE

          CASE( 31 )
            CSPT31( K,ISP ) = SPCODE

!.............  MACT specific entries
          CASE( 32 )
            CSPT32( K,ISP ) = SPCODE

          CASE( 33 )
            CSPT33( K,ISP ) = SPCODE

          CASE( 34 )
            CSPT34( K,ISP ) = SPCODE

          CASE( 35 )
            CSPT35( K,ISP ) = SPCODE

          CASE( 36 )
            CSPT36( K,ISP ) = SPCODE

          CASE( 37 )
            CSPT37( K,ISP ) = SPCODE

          CASE DEFAULT

        END SELECT

    ENDDO                            ! End Loop on sorted x-ref entries

    IF( EFLAG ) THEN
        MESG = 'Problem processing cross-reference records.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    RETURN

END SUBROUTINE FILLSTBL
