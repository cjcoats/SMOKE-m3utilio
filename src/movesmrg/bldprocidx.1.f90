
SUBROUTINE BLDPROCIDX

!***********************************************************************
!  subroutine BLDPROCIDX body starts at line
!
!  DESCRIPTION:
!      This subroutine builds a mapping from the list of MOVES emission
!      processes to the pollutant-species combinations.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 4/10 by C. Seppanen
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
!****************************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains data structures and flags specific to Movesmrg
    USE MODMVSMRG, ONLY: RPDFLAG, RPVFLAG, EMPOLIDX

!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: NSMATV, TSVDESC

    IMPLICIT NONE

!.........  INCLUDES:
    INCLUDE 'EMCNST3.EXT'    !  emissions constant parameters
    INCLUDE 'MVSCNST3.EXT'   !  MOVES constants

!.........  Other local variables
    INTEGER   J, K, L1, V  ! counters and indexes
    INTEGER   IOS       ! error status

    CHARACTER(NAMLEN3) :: CPROC ! tmp process buffer
    CHARACTER(PLSLEN3) :: SVBUF ! tmp speciation name buffer
    CHARACTER(300) :: MESG      ! message buffer

    CHARACTER(16) :: PROGNAME = 'BLDPROCIDX' ! program name

!***********************************************************************
!   begin body of subroutine BLDPROCIDX

    ALLOCATE( EMPOLIDX( NSMATV ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMPOLIDX', PROGNAME )
    EMPOLIDX = 0   ! array

!.........  Loop through pollutant-species combos
    DO V = 1, NSMATV

        SVBUF = TSVDESC( V )
        L1 = INDEX( SVBUF, ETJOIN )
        CPROC = SVBUF( 1:L1-1 )

!.............  Find emission process in master MOVES list
        IF( RPDFLAG ) THEN
            J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
        ELSE IF( RPVFLAG ) THEN
            J = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
        ELSE
            J = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
        END IF

        IF( J .LE. 0 ) THEN

!.................  Check if process is handled by other modes
            IF( RPDFLAG ) THEN
                J = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
                K = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
            ELSE IF( RPVFLAG ) THEN
                J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
                K = INDEX1( CPROC, MXMVSPPROCS, MVSPPROCS )
            ELSE
                J = INDEX1( CPROC, MXMVSDPROCS, MVSDPROCS )
                K = INDEX1( CPROC, MXMVSVPROCS, MVSVPROCS )
            END IF

            IF( J .LE. 0 .AND. K .LE. 0 ) THEN
                MESG = 'ERROR: Requested emission process ' //&
                &  TRIM( CPROC ) // ' is not output by MOVES.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

    END DO

END SUBROUTINE BLDPROCIDX
