
SUBROUTINE RDMAPPOL( NSRC, NVARS, NPVAR, VARNAMS, OUTVALS )

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!     Reads in pollutant or activity data and associated
!     variables from map-formatted or old-format inventory files.
!     Stores the data by source instead of by sparse record.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!************************************************************************
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
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NMAP, MAPNAM, MAPFIL

    USE MODFILESET

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: NSRC                 ! no. inven sources
    INTEGER     , INTENT (IN) :: NVARS                ! number of pol/act
    INTEGER     , INTENT (IN) :: NPVAR                ! number of vars per pol/act
    CHARACTER(*), INTENT (IN) :: VARNAMS( NVARS )     ! vars to map
    REAL        , INTENT(OUT) :: OUTVALS( NSRC, NVARS*NPVAR )! output data

!...........   Local allocatable arrays
    INTEGER, ALLOCATABLE :: SRCID( :   )  ! source IDs in sparse pol files
    REAL   , ALLOCATABLE :: ETMP ( :,: )  ! emissions in sparse pol files

!...........   Other local variables
    INTEGER          I, J, L, M, N, S, V   ! counters and indices

    INTEGER          ADJINDX      ! adjustment index for special average day read
    INTEGER          IOS          ! i/o status
    INTEGER          NSPARSE      ! number of rows in sparse pol input file
    INTEGER          NSRCLOC      ! number of rows in sparse pol input file

    LOGICAL       :: EFLAG = .FALSE.  ! true: error found
    LOGICAL       :: OFLAG = .FALSE.  ! true: old inventory format
    LOGICAL       :: SFLAG = .FALSE.  ! true: error found on read of inven data

    CHARACTER(16)   :: RNAME = 'IOAPI_DAT' ! logical name for reading pols
    CHARACTER(256)     MESG   ! message buffer
    CHARACTER(NAMLEN3) VBUF   ! tmp variable name buffer

    CHARACTER(16) :: PROGNAME = 'RDMAPPOL'   !  program name

!***********************************************************************
!   begin body of subroutine RDMAPPOL

!........  Initialize output arrays to zero
    OUTVALS = 0.   ! array

    DO V = 1, NVARS

        ADJINDX = 2
        VBUF = VARNAMS( V )

!............  Check if average day prefix is in the name
        IF( VARNAMS( V )( 1:3 ) .EQ. AVEDAYRT ) THEN
            L = LEN_TRIM( VARNAMS( V ) )
            VBUF = VARNAMS( V )( 4:L )
            ADJINDX = 3

!...............  Give error if reading more than a single variable
            IF( NPVAR .GT. 1 ) THEN
                MESG = 'INTERNAL ERROR: Cannot specify average ' //&
                &       'day value and have NPVAR > 1'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

        END IF

!............  Find variable in map
        M = INDEX1( VBUF, NMAP, MAPNAM )
        IF( M .LE. 0 ) THEN
            MESG = 'INTERNAL ERROR: Could not find variable "'//&
            &       TRIM( VBUF ) // '" in map file'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!............  Open physical file name for this pollutant or activity
!............  Also, get description
        CALL OPENPHYS( PROGNAME, RNAME, FSREAD3, MAPFIL( M ),&
        &               EFLAG )

        NSPARSE = NROWS3D

!............  Allocate memory for local read arrays
        IF( ALLOCATED( ETMP ) ) DEALLOCATE( ETMP, SRCID )
        ALLOCATE(  ETMP( NSPARSE, NPVAR ),&
        &          SRCID( NSPARSE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETMP,SRCID', PROGNAME )

!.............  Use subroutine to read in data. This subroutine
!               already handles integer and real
        J = ADJINDX
        CALL RDINVPOL( RNAME, NSPARSE, NPVAR, VNAMESET( J ),&
        &               VTYPESET( J ), SRCID, ETMP, SFLAG     )

        IF( SFLAG ) EFLAG = .TRUE.
        IF( EFLAG ) CYCLE

!...........  Transfer sparse-storage pollutant to output arrays
        DO I = 1, NPVAR

            J = I + ( V-1 ) * NPVAR
            DO N = 1, NSPARSE
                S = SRCID( N )
                OUTVALS( S,J ) = ETMP( N,I )
            END DO

        END DO     ! End loop over variable per pol or act (if any)

!...............  Close output file for this variable
        IF( .NOT. CLOSESET( RNAME ) ) THEN
            MESG = 'Could not close file:'//CRLF()//BLANK10//&
            &           TRIM( MAPFIL( M ) )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END DO         ! End loop over pols or acts

!.........  Abort if an error has been found
    IF( EFLAG ) THEN

        MESG = 'Problem reading pollutant data'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    IF( ALLOCATED( SRCID ) ) DEALLOCATE( SRCID, ETMP )

    RETURN

END SUBROUTINE RDMAPPOL


