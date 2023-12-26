
SUBROUTINE WPINGSTK( FNAME, SDATE, STIME, LFLAG )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !       Opens the output files for the Elevpoint program
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 8/99 by M Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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
    !***********************************************************************
    USE M3UTILIO

    !.......   MODULES for public variables
    !.......  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: NGROUP, GRPIDX, GRPGIDA, GRPCNT, GRPCOL,    &
                       GRPROW, GRPDM, GRPFL, GRPHT, GRPLAT, GRPLON,    &
                       GRPTK, GRPVE, GRPXL, GRPYL, GRPFIP,GRPLMAJOR,    &
                       GRPLPING, GRPACRES, FFLAG

    IMPLICIT NONE

    !......    Subroutine arguments and their descriptions
    CHARACTER(*), INTENT (IN) :: FNAME       ! i/o api inventory file
    INTEGER     , INTENT (IN) :: SDATE       ! Julian start date
    INTEGER     , INTENT (IN) :: STIME       ! start time
    LOGICAL     , INTENT (IN) :: LFLAG       ! true: write lat/lon

    !.......   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'WPINGSTK'       !  subroutine name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.......   Local arrays
    !.......   These are for sorting groups and outputting in sorted order
    INTEGER ::    LOCGID( NGROUP )
    INTEGER ::    LOCCNT( NGROUP )
    INTEGER ::    LOCCOL( NGROUP )
    INTEGER ::    LOCROW( NGROUP )
    INTEGER ::    LOCFIP( NGROUP )
    INTEGER :: LOCLMAJOR( NGROUP )
    INTEGER ::  LOCLPING( NGROUP )

    REAL    ::     LOCDM( NGROUP )
    REAL    ::     LOCFL( NGROUP )
    REAL    ::     LOCHT( NGROUP )
    REAL    ::    LOCLAT( NGROUP )
    REAL    ::    LOCLON( NGROUP )
    REAL    ::     LOCTK( NGROUP )
    REAL    ::     LOCVE( NGROUP )
    REAL    ::     LOCXL( NGROUP )
    REAL    ::     LOCYL( NGROUP )
    REAL    ::  LOCACRES( NGROUP )

    !.......   Other local variables
    INTEGER         I, J        ! indices and counters

    CHARACTER(300)  MESG

    !***********************************************************************
    !   begin body of subroutine WPINGSTK

    !.......  Store sorted information
    DO I = 1, NGROUP
        J = GRPIDX( I )

        LOCGID( I ) = GRPGIDA( J )
        LOCCNT( I ) = GRPCNT( J )
        LOCCOL( I ) = GRPCOL( J )
        LOCROW( I ) = GRPROW( J )
        LOCDM ( I ) = GRPDM ( J )
        LOCFL ( I ) = GRPFL ( J )
        LOCHT ( I ) = GRPHT ( J )
        LOCLAT( I ) = GRPLAT( J )
        LOCLON( I ) = GRPLON( J )
        LOCTK ( I ) = GRPTK ( J )
        LOCVE ( I ) = GRPVE ( J )
        LOCXL ( I ) = GRPXL ( J )
        LOCYL ( I ) = GRPYL ( J )
        LOCFIP( I ) = STR2INT( GRPFIP( J ) )
        LOCLMAJOR( I ) = GRPLMAJOR( J )
        LOCLPING( I ) = GRPLPING( J )
        IF (FFLAG) LOCACRES( I ) = GRPACRES( J)
    END DO

    MESG = 'Error writing to output file "' //TRIM( FNAME )// '"'

    IF ( .NOT. WRITE3( FNAME, 'ISTACK', SDATE, STIME, LOCGID )) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( LFLAG ) THEN         ! skip writing Lat and Lon variables.
        IF( .NOT. WRITE3( FNAME,'LATITUDE',SDATE,STIME,LOCLAT )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IF( .NOT. WRITE3( FNAME,'LONGITUDE',SDATE,STIME,LOCLON )) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKDM', SDATE, STIME, LOCDM ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKHT', SDATE, STIME, LOCHT ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKTK', SDATE, STIME, LOCTK ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKVE', SDATE, STIME, LOCVE ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKFLW', SDATE, STIME, LOCFL ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'STKCNT', SDATE, STIME, LOCCNT )) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'ROW', SDATE, STIME, LOCROW ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'COL', SDATE, STIME, LOCCOL ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'XLOCA', SDATE, STIME, LOCXL ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME, 'YLOCA', SDATE, STIME, LOCYL ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF


    IF ( .NOT. WRITE3( FNAME, 'IFIP', SDATE, STIME, LOCFIP ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. WRITE3( FNAME,'LMAJOR',SDATE,STIME,LOCLMAJOR ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF


    IF ( .NOT. WRITE3( FNAME,'LPING',SDATE, STIME, LOCLPING ) ) THEN
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF (FFLAG) THEN
        IF ( .NOT. WRITE3( FNAME,'ACRESBURNED',SDATE, STIME, LOCACRES ) ) THEN
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

    RETURN

END SUBROUTINE WPINGSTK

