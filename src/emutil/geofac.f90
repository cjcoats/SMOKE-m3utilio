
PROGRAM GEOFAC

    !***********************************************************************
    !
    !  DESCRIPTION: Takes gridded SMOKE emissions file and applies
    !               a user-supplied factor for each individual species
    !               and applies it for a certain geographical region
    !               obtained from input mask file..
    !               The resulting emissions are output to a I/O API file.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Prototype 10/2000 :  by JMV
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
    !*************************************************************************
    USE M3UTILIO

    USE MODFILESET

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'         ! Emissions constants
    INCLUDE 'SETDECL.h90'       !  FileSetAPI function declarations

    !.......  EXTERNAL FUNCTIONS and their descriptions:

    INTEGER,EXTERNAL :: GETFLINE

    !.......   PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'GEOFAC'       !  program name
    CHARACTER(50), PARAMETER :: CVSW =&
         '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.......  LOCAL VARIABLES

    REAL  FEMIS                                  ! temporary value for emissions
    REAL, ALLOCATABLE :: EMIS( :, : , : )            ! emissions
    REAL, ALLOCATABLE :: SPCFAC( : )             ! species factors

    INTEGER  TSTEP                               ! time step
    INTEGER  I, J, K , L, M, N                      ! counters
    INTEGER, ALLOCATABLE ::  IMASK ( :, : )      ! mask values
    INTEGER  LDEV                                ! log file unit number
    INTEGER  RDEV                                ! species factors unit number
    INTEGER  NSPECS                              ! number of species
    INTEGER  HR                                  ! hour loop counter
    INTEGER  NLINES                              ! number of species factors
    INTEGER  NSTEPS                              ! number of time steps
    INTEGER  IOS                                 ! iostat
    INTEGER  SDATE                               ! start date
    INTEGER  STIME                               ! start time
    INTEGER  NLAYS                               ! number of layers in emis file

    CHARACTER(16), ALLOCATABLE :: SPCNAM ( : )     ! species names with facs
    CHARACTER(16)  ENAME                           ! logical name for gridded emis input file
    CHARACTER(16)  ONAME                           ! logical name for output file
    CHARACTER(16)  MNAME                           ! logical name for mask file
    CHARACTER(300) MESG                            ! message buffer for M3EXIT()


    !***********************************************************************
    !   begin body of program

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.

    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Prompt for name of I/O API input file

    ENAME = PROMPTSET( 'Enter logical name for SMOKE gridded input (I/O API) file',&
                       FSREAD3, 'INFILE', PROGNAME )

    IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN

        MESG = 'Could not get description of file "' // TRIM( ENAME ) // '"'
        CALL M3EXIT( PROGNAME , 0, 0, MESG, 2 )

    ENDIF

    !.......  Assign local variables

    TSTEP  = TSTEP3D
    SDATE  = SDATE3D
    STIME  = STIME3D
    NSPECS = NVARS3D
    NSTEPS = MXREC3D
    NLAYS  = NLAYS3D

    !.......   Get the mask file name

    MNAME = PROMPTMFILE( 'Enter logical name for mask input (I/O API) file',&
                         FSREAD3, 'GEOMASK', PROGNAME )

    !.......   Get the species fac file name

    RDEV = PROMPTFFILE( 'Enter logical name for factors file',  &
                        .TRUE., .TRUE., 'SPECFACS', PROGNAME )

    !.......  Read mask file header

    ALLOCATE( IMASK ( NCOLS3D, NROWS3D ), STAT=IOS )
    CALL CHECKMEM( IOS, 'IMASK', PROGNAME )


    IF ( .NOT. READ3( MNAME, 'IFAC', 1, 0, 0, IMASK ) ) THEN
        CALL M3EXIT( PROGNAME , 0, 0, 'Error reading factors from file ' // MNAME, 2 )
    ENDIF

    !.......  Read species factors

    NLINES = GETFLINE( RDEV, 'Species factors file' )

    ALLOCATE( SPCNAM( NLINES ),     &
              SPCFAC( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPCFAC', PROGNAME )

    !.......  Create file descriptions for output file

    FDESC3D  = '  '
    FDESC3D( 1 ) = '/' // PROGNAME // '/'

    !......   Read in factors for species

    DO N = 1, NLINES
        READ( RDEV, 93000 ) SPCNAM ( N ) , SPCFAC ( N )
        M = N + 1
        WRITE ( FDESC3D( M ) , 94000 ) SPCNAM( N ), SPCFAC( N )
    ENDDO

    !.......  Open output file

    ONAME = PROMPTSET( 'Enter logical name for OUTPUT gridded file ',&
                       FSUNKN3, 'OUTFILE',  PROGNAME      )

    !.......  Allocate memory for emissions

    ALLOCATE( EMIS( NCOLS3D, NROWS3D, NLAYS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMIS', PROGNAME )

    !.......  Loop over time steps

    DO HR = 1, NSTEPS

        !.......  Write to screen because WRITE3 only writes to LDEV

        DO  L = 1, NSPECS

            EMIS = 0.0           !  array

            !.......  Read input file for time and species of interest

            IF ( .NOT. READSET( ENAME, VNAMESET( L ), ALLAYS3,&
                                ALLFILES, SDATE, STIME, EMIS   ) ) THEN
                CALL M3EXIT( PROGNAME , 0, 0,&
                               'Error reading ' // VNAMESET( L ) //&
                               ' from file ' // ENAME , 2 )
            ENDIF

            DO M = 1, NLINES

                !.......   Find out if factor for this species exists

                IF ( VNAMESET ( L ) .EQ. SPCNAM ( M )  ) THEN

                    !.......  Loop through cells and apply factor to cells in mask

                    DO K = 1, NLAYS
                    DO J = 1, NROWS3D
                    DO I = 1, NCOLS3D

                        !.......   Find out if grid cell is in geographical region of interest

                        IF ( IMASK ( I,J ) .EQ. 1 ) THEN
                            EMIS ( I,J,K ) = EMIS ( I,J,K ) * SPCFAC( M )
                        END IF

                    END DO
                    END DO
                    END DO

                END IF

                !.......  Finished for this species

            ENDDO

            !......  Write out new emissions

            IF ( .NOT. WRITESET( ONAME, VNAMESET( L ), -1, SDATE, STIME, EMIS ) ) THEN

                CALL M3EXIT( PROGNAME , SDATE, STIME,&
                           'Could not write "' // VNAMESET( L ) // '" to ' // ONAME, 2 )

            END IF             !  if writeset() failed

        ENDDO

    !......  Next time step

        CALL NEXTIME( SDATE, STIME, TSTEP)

    ENDDO

    !.......   End of program:

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Error and warning message formats..... 91xxx

    !.......   Informational (LOG) message formats... 92xxx

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A, F6.3 )

    !.......   Internal buffering formats...... 94xxx

94000 FORMAT ( '/SPECIES FAC/ ', A, F6.3 )

END PROGRAM  GEOFAC

