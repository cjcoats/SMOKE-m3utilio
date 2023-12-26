
PROGRAM METCOMBINE

        !***********************************************************************
        !  program body starts at line
        !
        !  DESCRIPTION:
        !       This program combines 2D and the 1st layer of 3D meteorology
        !       files to create custom files for mobile processing with MOBILE6.
        !
        !  PRECONDITIONS REQUIRED:
        !       Postprocessed MM5 meteorology that contains temperature data
        !
        !  SUBROUTINES AND FUNCTIONS CALLED:
        !       CHKGRID
        !
        !  REVISION  HISTORY:
        !       Created 2/2004 by C. Seppanen
        !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
        !       related changes
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

    IMPLICIT NONE

        !.......  LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'METCOMBINE'         ! program name
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'          ! CVS release tag

        !.......  LOCAL VARIABLES

        !.......  Arrays
    INTEGER        VARFOUND( MXVARS3 )          ! stores which file each variable is in
    CHARACTER(16)   VARLIST( MXVARS3 )         ! list of variables to read
    CHARACTER(16)   VARUNIT( MXVARS3 )          ! stores the units for each variable
    CHARACTER(80)   VARDESC( MXVARS3 )          ! stores the description for each variable
    INTEGER         VARTYPE( MXVARS3 )          ! stores the type of each variable

    REAL, ALLOCATABLE :: VARDATA( : )          ! generic array to store variable data

        !.......  File units and logical names
    INTEGER         LDEV                            ! unit number for log file

    CHARACTER(16)   M1NAME                          ! logical name for 1st met file
    CHARACTER(16)   M2NAME                          ! logical name for 2nd met file
    CHARACTER(16)   INAME                           ! tmp name for met file
    CHARACTER(16)   ONAME                           ! logical name for output file

        !.......  Other local variables
    INTEGER         I, J, T                         ! counters
    INTEGER         IOS                             ! I/O status
    INTEGER         JDATE                           ! loop date
    INTEGER         JTIME                           ! loop time
    INTEGER         NSTEPS                          ! number of time steps
    INTEGER         NVARS                           ! actual number of variables
    INTEGER         SDATE                           ! start date
    INTEGER         STIME                           ! start time
    INTEGER         TSTEP                           ! time step

    LOGICAL      :: GRID_ERR = .FALSE.              ! true: error found in grid settings

    CHARACTER(256)  MESG                            ! message buffer


    !***********************************************************************
    !   begin body of program METCOMBINE

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Get environment variable settings...
    !.......  Get name of temperature variable to extract
    MESG = 'Meteorology variable names'
    IF( .NOT. STRLIST( 'VARLIST', MESG, MXVARS3,&
                        NVARS, VARLIST ) ) THEN
        MESG = 'Could not read list of meteorology variable names'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    VARFOUND = 0
    VARUNIT = ' '
    VARDESC = ' '
    VARTYPE = 0

    !.......  Open input meteorology files
    MESG = 'Enter logical name for the first meteorology file'
    M1NAME = 'METFILE1'
    M1NAME = PROMPTMFILE( MESG, FSREAD3, M1NAME, PROGNAME )

    MESG = 'Enter logical name for the second meteorology file'
    M2NAME = 'METFILE2'
    M2NAME = PROMPTMFILE( MESG, FSREAD3, M2NAME, PROGNAME )

    !.......  Get headers from files and compare grids and dates

    !.......  Get description of first file
    IF( .NOT. DESC3( M1NAME ) ) THEN
        MESG = 'Could not get description of file ' // TRIM( M1NAME )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Initialize grid settings
    CALL CHKGRID( M1NAME, 'GRID', 0, GRID_ERR )

        !.......  Save date and time settings
    SDATE = SDATE3D
    STIME = STIME3D
    TSTEP = TSTEP3D
    NSTEPS = MXREC3D

    !.......  Check that requested variables are in file
    DO I = 1, NVARS
        J = INDEX1( VARLIST( I ), NVARS3D, VNAME3D )
        IF( J > 0 ) THEN
            VARFOUND( I ) = 1
            VARUNIT( I ) = UNITS3D( J )
            VARDESC( I ) = VDESC3D( J )
            VARTYPE( I ) = VTYPE3D( J )
        END IF
    END DO

        !.......  Get description of second file
    IF( .NOT. DESC3( M2NAME ) ) THEN
        MESG = 'Could not get description of file ' // TRIM( M2NAME )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

        !.......  Check grid settings
    CALL CHKGRID( M2NAME, 'GRID', 0, GRID_ERR )

    IF( GRID_ERR ) THEN
        MESG = 'Grids are inconsistent between meteorology files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

        !.......  Check date and time settings
    IF( SDATE3D /= SDATE ) THEN
        MESG = 'Start dates are inconsistent between meteorology files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( STIME3D /= STIME ) THEN
        MESG = 'Start times are inconsistent between meteorology files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( TSTEP3D /= TSTEP ) THEN
        MESG = 'Time steps are inconsistent between meteorology files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( NSTEPS /= MXREC3D ) THEN
        MESG = 'Number of time steps is inconsistent between meteorology files.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

        !.......  Check for requested variables
    DO I = 1, NVARS
        J = INDEX1( VARLIST( I ), NVARS3D, VNAME3D )
        IF( J > 0 ) THEN
            IF( VARFOUND( I ) == 1 ) THEN
                MESG = 'Variable ' // TRIM( VARLIST( I ) ) //&
                       ' is in both meteorology files.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            ELSE
                VARFOUND( I ) = 2
                VARUNIT( I ) = UNITS3D( J )
                VARDESC( I ) = VDESC3D( J )
                VARTYPE( I ) = VTYPE3D( J )
            END IF
        END IF
    END DO

        !.......  Check that all variables are in the files
    DO I = 1, NVARS
        IF( VARFOUND( I ) == 0 ) THEN
            MESG = 'Variable ' // TRIM( VARLIST(I) ) // ' not in either meteorology file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END DO

    !.......  Open new output file

    !.......  Initialize file header with first meteorology file
    IF( .NOT. DESC3( M1NAME ) ) THEN
        MESG = 'Could not get description of file ' // TRIM( M1NAME )
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    NLAYS3D = 1
    NVARS3D = NVARS

    VNAME3D = ' '
    VTYPE3D = 0
    UNITS3D = ' '
    VDESC3D = ' '

    DO I = 1, NVARS
        VNAME3D( I ) = VARLIST( I )
        UNITS3D( I ) = VARUNIT( I )
        VDESC3D( I ) = VARDESC( I )
        VTYPE3D( I ) = VARTYPE( I )
    END DO

    MESG = 'Enter logical name for output meterology file'
    ONAME = 'OUTFILE'
    ONAME = PROMPTMFILE( MESG, FSUNKN3, ONAME, PROGNAME )

    !.......  Allocate memory for storing data
    ALLOCATE( VARDATA( NCOLS3D*NROWS3D ), STAT=IOS )
    CALL CHECKMEM( IOS, 'VARDATA', PROGNAME )

    !.......  Loop through time steps of met data
    JDATE = SDATE
    JTIME = STIME

    DO T = 1, NSTEPS

        DO I = 1, NVARS

            IF( VARFOUND( I ) == 1 ) THEN
                INAME = M1NAME
            ELSE
                INAME = M2NAME
            END IF

            IF( .NOT. READ3( INAME, VARLIST( I ), 1, JDATE, JTIME, VARDATA ) ) THEN
                MESG = 'Could not read variable "' // TRIM( VARLIST( I ) ) //    &
                       '" from file ' // TRIM( INAME )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

            IF( .NOT. WRITE3( ONAME, VARLIST( I ), JDATE, JTIME, VARDATA ) ) THEN
                MESG = 'Could not write variable "' // TRIM( VARLIST( I ) ) //   &
                       '" to file ' // TRIM( ONAME )
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            END IF

        END DO

        CALL NEXTIME( JDATE, JTIME, TSTEP )

    END DO

        !.......  End program successfully
    CALL M3EXIT( PROGNAME, 0, 0, 'Successful completion', 0 )

        !******************  FORMAT  STATEMENTS   ******************************

        !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

END PROGRAM METCOMBINE

