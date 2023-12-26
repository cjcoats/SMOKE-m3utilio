
SUBROUTINE OPENEOUT( NGROUP, SDATE, STIME, ENAME, VFLAG, LFLAG, PDEV, MNAME )

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
    !       Create 8/99 by M Houyoux
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
    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CRL, CATDESC
    USE MODELEV, ONLY: FFLAG
    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !......    Subroutine arguments and their descriptions
    INTEGER     , INTENT (IN) :: NGROUP      ! number of ping groups
    INTEGER     , INTENT (IN) :: SDATE       ! start date of episode
    INTEGER     , INTENT (IN) :: STIME       ! start time of episode
    CHARACTER(*), INTENT (IN) :: ENAME       ! i/o api inventory file
    LOGICAL     , INTENT (IN) :: VFLAG       ! true: using variable grid
    LOGICAL     , INTENT (IN) :: LFLAG       ! true: write lat/lon
    INTEGER     , INTENT (OUT):: PDEV        ! ASCII file for major/ping src IDs
    CHARACTER(*), INTENT (OUT):: MNAME       ! logical name of ping srcs groups

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL      , EXTERNAL :: DSCM3GRD
    CHARACTER(50), EXTERNAL :: GETCFDSC
    CHARACTER(16), EXTERNAL :: VERCHAR

    !.......   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENEOUT'       !  subroutine name
    CHARACTER(50), PARAMETER ::  CVSW    = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag

    !.......   Other local variables
    INTEGER         J
        ! indices and counters
    LOGICAL         EFLAG                   !  true: error found

    CHARACTER(80)   GDESC                   !  grid description
    CHARACTER(300)  MESG

    CHARACTER(NAMLEN3) NAMBUF               ! file name buffer
    CHARACTER(NAMLEN3) COORD3D
    CHARACTER(NAMLEN3) COORUN3D
    CHARACTER(MXDLEN3) IFDESC2, IFDESC3     ! fields 2  3 from inven FDESC

    !***********************************************************************
    !   begin body of subroutine OPENEOUT

    !.......  Open ASCII output file
    EFLAG = .FALSE.
    PDEV  = PROMPTFFILE( 'Enter name for ELEVATED POINT SOURCE output file',    &
                         .FALSE., .TRUE., CRL // 'ELV', PROGNAME )

    !.......  If needed, set up and open plume-in-grid i/o api output file
    IF( NGROUP .GT. 0 ) THEN

        !.......  Get header information from inventory file

        IF ( .NOT. DESCSET( ENAME,-1 ) ) THEN
            MESG = 'Could not get description of file "' // TRIM( ENAME ) // '".'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        IFDESC2 = GETCFDSC( FDESC3D, '/FROM/', .TRUE. )
        IFDESC3 = GETCFDSC( FDESC3D, '/VERSION/', .TRUE. )

        !.......  Set up for opening I/O API output file header
        CALL HDRMISS3          ! Initialize for emissions

        !.......  Get grid description for setting most grid parameters for
        !               the output file
        IF( .NOT. DSCM3GRD(    &
                  GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,   &
                  P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,  &
                  XORIG3D, YORIG3D, XCELL3D, YCELL3D,           &
                  NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

            MESG = 'Could not get Models-3 grid description.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        COORUN3D = 'METERS '
        IF ( GDTYP3D .EQ. 1 ) COORUN3D = 'DEGREES '

        !....... Finalize i/o api header fields
        !....... NOTE - this is a time-independent file, but the Plume Dynamics
        !              Model that reads this file needs to have the start date and time
        !              of the Met file in here.
        FDESC3D( 1 ) = CATDESC // ' source stack groups file'
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS3D
        WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS3D
        FDESC3D( 11 ) = '/INVEN FROM/ ' // IFDESC2
        FDESC3D( 12 ) = '/INVEN VERSION/ ' // IFDESC3

        IF( VFLAG ) THEN         ! when processing variable grid
            FDESC3D( 13 ) = '/VARIABLE GRID/ ' // GDNAM3D
        END IF

        NVARS3D = 14       ! Start with 14, but add more depending on options

        NROWS3D = NGROUP
        NCOLS3D = 1
        NLAYS3D = 1
        SDATE3D = SDATE
        STIME3D = STIME
        TSTEP3D = 10000

        !.......  Set the file variables
        J = 1
        VNAME3D( J ) = 'ISTACK'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Stack group number'
        J = J + 1

        IF( LFLAG ) THEN        ! skip writing Lat and Lon variables
            VNAME3D( J ) = 'LATITUDE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Latitude'
            J = J + 1

            VNAME3D( J ) = 'LONGITUDE'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'degrees'
            VDESC3D( J ) = 'Longitude'
            J = J + 1

            NVARS3D = NVARS3D + 2
        END IF

        VNAME3D( J ) = 'STKDM'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Inside stack diameter'
        J = J + 1

        VNAME3D( J ) = 'STKHT'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm'
        VDESC3D( J ) = 'Stack height above ground surface'
        J = J + 1

        VNAME3D( J ) = 'STKTK'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'degrees K'
        VDESC3D( J ) = 'Stack exit temperature'
        J = J + 1

        VNAME3D( J ) = 'STKVE'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm/s'
        VDESC3D( J ) = 'Stack exit velocity'
        J = J + 1

        VNAME3D( J ) = 'STKFLW'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = 'm**3/s'
        VDESC3D( J ) = 'Stack exit flow rate'
        J = J + 1

        VNAME3D( J ) = 'STKCNT'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Number of stacks in group'
        J = J + 1

        VNAME3D( J ) = 'ROW'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Grid row number'
        J = J + 1

        VNAME3D( J ) = 'COL'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Grid column number'
        J = J + 1

        VNAME3D( J ) = 'XLOCA'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = COORUN3D
        VDESC3D( J ) = 'Projection x coordinate'
        J = J + 1

        VNAME3D( J ) = 'YLOCA'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = COORUN3D
        VDESC3D( J ) = 'Projection y coordinate'

        J = J + 1

        VNAME3D( J ) = 'IFIP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'FIPS CODE'

        J = J + 1

        VNAME3D( J ) = 'LMAJOR'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = '1= MAJOR SOURCE in domain, 0=otherwise'

        J = J + 1

        VNAME3D( J ) = 'LPING'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = '1=PING SOURCE in domain, 0=otherwise'

        !....... If fire, then add acres burned to dataset
        IF ( FFLAG ) THEN
            J = J + 1
            NVARS3D = NVARS3D + 1

            VNAME3D( J ) = 'ACRESBURNED'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'acres/day'
            VDESC3D( J ) = 'number of acres burned for a fire in one day'
        ENDIF

        MESG = 'Enter logical name for ELEVATED STACK GROUPS file'

        !.......  Open file. Use NAMBUF for HP.
        MNAME  = 'STACK_GROUPS'
        NAMBUF = PROMPTMFILE( MESG, FSUNKN3, MNAME, PROGNAME )
        MNAME  = NAMBUF

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************


    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10 ( A, :, I8, :, 2X  ) )

94015 FORMAT( A, 1X, I6, 1X, A, 1X, I5.5, 1X, A, 1X, I8.8, 1X,    &
            A, 1X, I6, 1X, A, 1X, I6,   1X, A )

94020 FORMAT( 5X, 'H[m]:', 1X, F6.2, 1X, 'D[m]:'  , 1X, F4.2, 1X,    &
                'T[K]:', 1X, F7.1, 1X, 'V[m/s]:', 1X, F10.1 )

END SUBROUTINE OPENEOUT

