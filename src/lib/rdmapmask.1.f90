
SUBROUTINE RDMAPMASK( ENAME, NMAPVAR, MAPVARS, MAPFILES, NSRC,&
&                     NVAR, NLOOP, INVARS, RDVARS, MASK, OUTVALS)

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!     Reads in pollutants from map-formatted inventory files.
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

    IMPLICIT NONE

!...........   INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: ENAME                ! i/o api logical file name
    INTEGER     , INTENT (IN) :: NMAPVAR              ! no. map-file variables
    CHARACTER(*), INTENT (IN) :: MAPVARS ( NMAPVAR )  ! map-file variables
    CHARACTER(*), INTENT (IN) :: MAPFILES( NMAPVAR )  ! map-file files
    INTEGER     , INTENT (IN) :: NSRC                 ! no. inven sources
    INTEGER     , INTENT (IN) :: NVAR                 ! dimension 2 of OUTVALS
    INTEGER     , INTENT (IN) :: NLOOP                ! no. vars to read
    CHARACTER(*), INTENT (IN) :: INVARS ( NLOOP )     ! var names
    CHARACTER(*), INTENT (IN) :: RDVARS ( NLOOP )     ! vars to read
    INTEGER     , INTENT (IN) :: MASK   ( NLOOP )     ! variable index to dim2 of OUTVALS
    REAL        , INTENT(OUT) :: OUTVALS( NSRC, NVAR )! output emissions

!...........   Local allocatable arrays
    INTEGER, ALLOCATABLE :: SRCID( : )  ! source IDs in sparse pol files
    REAL   , ALLOCATABLE :: ETMP ( : )  ! emissions in sparse pol files

!...........   Other local variables
    INTEGER          I, K, N, S   ! counters and indices

    INTEGER          IOS          ! i/o status
    INTEGER          JIDX         ! mask index
    INTEGER          NSPARSE      ! number of rows in sparse pol input file

    LOGICAL       :: EFLAG = .FALSE.  ! true: error found
    LOGICAL       :: SFLAG = .FALSE.  ! true: error found on read of inven data

    CHARACTER(256)     MESG   ! message buffer
    CHARACTER(NAMLEN3) PNAME  ! logical file name for data files
    CHARACTER(NAMLEN3) VBUF   ! tmp variable name

    CHARACTER(16) :: PROGNAME = 'RDMAPMASK'   !  program name

!***********************************************************************
!   begin body of subroutine RDMAPMASK

    PNAME = 'IOAPI_DAT'

!........  Initialize output emissions
    OUTVALS = 0.   ! array

    DO I = 1, NLOOP

!.............  Do not read variable if mask is zero
        JIDX = MASK( I )
        IF( JIDX .LE. 0 ) CYCLE

        VBUF = RDVARS( I )

!.............  Bounds check
        IF( JIDX .GT. NVAR ) THEN

            MESG = 'INTERNAL ERROR: Prevented overflow for read '//&
            &       'of variable "' // TRIM( VBUF ) // '"'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

!.............  Find variable to read in map list of available vars
        K = INDEX1( VBUF, NMAPVAR, MAPVARS )

!.............  If variable is invalid, give error and go to get iteration
        IF( K .LE. 0 ) THEN
            EFLAG = .TRUE.
            MESG = 'ERROR: Variable "' // TRIM( VBUF ) //&
            &       '" requested by program but not found in ' //&
            &       'map-formatted inventory.'
            CALL M3MSG2( MESG )
            CYCLE
        END IF

!.............  Set temporary environment variable to use for reading
!               pollutant file. Also, get the header description.
        CALL OPENPHYS( PROGNAME, PNAME, FSREAD3, MAPFILES( K ),&
        &               EFLAG )

        NSPARSE = NROWS3D

!................  Allocate memory for local read arrays
        IF( ALLOCATED( ETMP ) ) DEALLOCATE( ETMP, SRCID )
        ALLOCATE( ETMP( NSPARSE ),&
        &         SRCID( NSPARSE ), STAT=IOS )
        CALL CHECKMEM( IOS, 'ETMP,SRCID', PROGNAME )

        CALL RDINVPOL( PNAME, NSPARSE, 1, VBUF, M3REAL, SRCID,&
        &               ETMP, SFLAG )

        IF( SFLAG ) EFLAG = .TRUE.
        IF( EFLAG ) CYCLE

!...........  Transfer sparse-storage pollutant to output arrays
        DO N = 1, NSPARSE
            S = SRCID( N )
            OUTVALS( S,JIDX ) = ETMP( N )
        END DO

!...........  If map-formatted inventory, close pollutant file
        IF( .NOT. CLOSESET( PNAME ) ) THEN
            MESG = 'Could not close file:'// CRLF()// BLANK10//&
            &       TRIM( MAPFILES( K ) )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END DO      ! end loop through input variables

!.........  Abort if an error has been found
    IF( EFLAG ) THEN

        MESG = 'Problem reading pollutant data'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END SUBROUTINE RDMAPMASK


