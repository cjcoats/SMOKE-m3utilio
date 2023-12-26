
SUBROUTINE OPENMRGOUT( NGRP )

!***********************************************************************
!  subroutine OPENMRGOUT body starts at line 86
!
!  DESCRIPTION:
!      The purpose of this subroutine is to open all of the necessary
!      output files for the merge routine (both I/O API and ASCII files)
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Created 2/99 by M. Houyoux
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
!.........  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: SDATE, STIME, TSTEP, BYEAR, PYEAR,&
    &        LAVEDAY, LGRDOUT,&
    &        AFLAG, ANIPOL, ANMSPC, AEINAM, AEMNAM, AONAME,&
    &        BFLAG, BNIPOL, BNMSPC, BEINAM, BEMNAM, BONAME,&
    &        MFLAG, MNMSPC, NMSPC, MNIPPA, MEANAM, EMNAM, MEMNAM,&
    &        MONAME, EMLAYS,&
    &        PFLAG, PNIPOL, PNMSPC, PEINAM, PEMNAM, PONAME,&
    &        XFLAG, NIPPA, EANAM, TONAME, PINGFLAG,&
    &        INLINEFLAG, PINGNAME, INLINENAME,&
    &        ELEVFLAG, EVDEV, PELVNAME, LREPSTA, LREPCNY,&
    &        AREPNAME, BREPNAME, MREPNAME, PREPNAME, TREPNAME,&
    &        ARDEV, BRDEV, MRDEV, PRDEV, TRDEV,&
    &        VGRPCNT, SIINDEX, SPINDEX, GRDUNIT, VARFLAG,&
    &        SRCGRPFLAG, NSGOUTPUT, SRCGRPNAME, SGINLNNAME,&
    &        SUBSECFLAG, NGRPS, SUBOUTNAME

!.........  This module contains arrays for plume-in-grid and major sources
    USE MODELEV, ONLY: NGROUP, SGFIREFLAG

!.........  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: GRDNM, NCOLS, NROWS, P_ALP, P_BET, P_GAM,&
    &                   XCENT, YCENT, XORIG, YORIG, XCELL, YCELL,&
    &                   GDTYP, VGTYP, VGTOP, VGLVS

    USE MODFILESET

    IMPLICIT NONE

!.........  INCLUDES:

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........  SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: NGRP     ! Actual number of groups

!.........  EXTERNAL FUNCTIONS and their descriptions:

    CHARACTER(16), EXTERNAL :: MULTUNIT
    CHARACTER(16), EXTERNAL :: VERCHAR

!...........  Local parameters
    CHARACTER(16), PARAMETER :: PROGNAME = 'OPENMRGOUT' ! program name
    CHARACTER(50), PARAMETER :: CVSW =&
    &    '$Name SMOKEv5.0_Jun2023$' ! CVS release tag

!.........  Base and future year per

!.........  Other local variables

    INTEGER         I, J, K, L, L1, L2, L3, L4, LD, LJ, N, V

    INTEGER         IOS, IOS1, IOS2        ! i/o ENVSTR statuses

    LOGICAL      :: FDFLAG= .FALSE.   ! true: output file meta desc
    LOGICAL      :: EFLAG = .FALSE.   ! true: error is found

    CHARACTER(16)   MRGFDESC     ! name for met data description
    CHARACTER(50)   RTYPNAM      ! name for type of state/county report file
    CHARACTER(50)   UNIT         ! tmp units buffer

    CHARACTER(128)  OUTSCEN      ! output scenario name
    CHARACTER(128)  METADESC     ! output meta description
    CHARACTER(128)  FILEDESC     ! output file description
    CHARACTER(256)  MESG         ! message buffer
    CHARACTER(256)  OUTDIR       ! output path
    CHARACTER(IODLEN3) DESCBUF ! variable description buffer


!***********************************************************************
!   begin body of subroutine OPENMRGOUT

!.........  Evaluate environment variables to possibly use to build output
!           physical file names if logical file names are not defined.
    CALL ENVSTR( 'OUTPUT', 'Output directory', ' ', OUTDIR, IOS1 )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "OUTPUT"', 2 )
    END IF
    IF ( IOS1 .EQ. 0 ) L3 = LEN_TRIM( OUTDIR )

    CALL ENVSTR( 'ESCEN', 'Emission scenario', ' ', OUTSCEN, IOS2 )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "ESCEN"', 2 )
    END IF
    IF ( IOS2 .EQ. 0 ) L4 = LEN_TRIM( OUTSCEN )

!.........  Set default output file names
    CALL MRGONAMS

!.........  Set the EV name for met data description
    MESG = 'Setting for the environment variable name for file ' //&
    &         'meta description for output file'
    CALL ENVSTR( 'MRG_FILEDESC', MESG, ' ', MRGFDESC, IOS  )
    IF( IOS >= 0 ) THEN
        FDFLAG = .TRUE.
        MESG = 'Use this file meta description for output file'
        CALL ENVSTR( MRGFDESC, MESG, ' ', METADESC, IOS )
        IF( IOS < 0 ) FDFLAG = .FALSE.
    END IF

!.........  Set up header for I/O API output files
    FTYPE3D = GRDDED3
    P_ALP3D = DBLE( P_ALP )
    P_BET3D = DBLE( P_BET )
    P_GAM3D = DBLE( P_GAM )
    XCENT3D = DBLE( XCENT )
    YCENT3D = DBLE( YCENT )
    XORIG3D = DBLE( XORIG )
    YORIG3D = DBLE( YORIG )
    XCELL3D = DBLE( XCELL )
    YCELL3D = DBLE( YCELL )
    SDATE3D = SDATE
    STIME3D = STIME
    TSTEP3D = TSTEP
    NCOLS3D = NCOLS
    NROWS3D = NROWS
    NTHIK3D = 1
    GDTYP3D = GDTYP
    VGTYP3D = VGTYP
    VGTOP3D = VGTOP
    GDNAM3D = GRDNM

    FDESC3D = ' '   ! array
    FDESC3D( 2 ) = '/FROM/ '    // PROGNAME
    FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
    WRITE( FDESC3D( 5 ),94010 ) '/BASE YEAR/ ', BYEAR
    IF( PYEAR .NE. BYEAR )&
    &    WRITE( FDESC3D( 6 ),94010 ) '/PROJECTED YEAR/ ', PYEAR
    IF( VARFLAG ) THEN
        FDESC3D( 7 ) = '/VARIABLE GRID/' // GRDNM
    END IF

!.........  Set average day description buffer
    DESCBUF = ' '
    IF( LAVEDAY ) DESCBUF = 'Average day value, '
    LD = MAX( LEN_TRIM( DESCBUF ), 1 )

!.........  Set up and open I/O API output file
    IF( LGRDOUT ) THEN

!.............  Prompt for and gridded open file(s)
        IF( AFLAG ) THEN
            CALL SETUP_VARIABLES( ANIPOL, ANMSPC, AEINAM, AEMNAM )
            NLAYS3D = 1
            FDESC3D( 1 ) = 'Area source emissions data'
            IF( FDFLAG ) FDESC3D( 1 ) = METADESC

!.................  Open by logical name or physical name
            FILEDESC = 'AREA-SOURCE GRIDDED OUTPUT file'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC,'NETCDF',AONAME,I )

!.................  Open by sub-sector output files
            IF( SUBSECFLAG ) THEN
                FILEDESC = 'Sub-sector groups gridded output file'
                DO I = 1, NGRPS
                    CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'NETCDF', SUBOUTNAME(I), K )
                END DO
            END IF

        END IF

        IF( BFLAG ) THEN
            CALL SETUP_VARIABLES( BNIPOL, BNMSPC, BEINAM, BEMNAM )
            NLAYS3D = 1
            FDESC3D( 1 ) = 'Biogenic source emissions data'
            IF( FDFLAG ) FDESC3D( 1 ) = METADESC

!.................  Open by logical name or physical name
            FILEDESC = 'BIOGENIC GRIDDED OUTPUT file'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC,'NETCDF',BONAME,I )

        END IF

        IF( MFLAG ) THEN

!cs...............  Check if the number of mobile species exceeds the
!cs                 number of total species (can happen if there are
!cs                 more species in the speciation matrix than desired
!cs                 in the output)
            IF( MNMSPC > NMSPC ) THEN
                CALL SETUP_VARIABLES( MNIPPA, NMSPC, MEANAM, EMNAM )
            ELSE
                CALL SETUP_VARIABLES( MNIPPA, MNMSPC,MEANAM, MEMNAM)
            ENDIF
            NLAYS3D = 1
            FDESC3D( 1 ) = 'Mobile source emissions data'
            IF( FDFLAG ) FDESC3D( 1 ) = METADESC

!.................  Open by logical name or physical name
            FILEDESC = 'MOBILE-SOURCE GRIDDED OUTPUT file'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC,'NETCDF',MONAME,I )


        END IF

        IF( PFLAG ) THEN
            CALL SETUP_VARIABLES( PNIPOL, PNMSPC, PEINAM, PEMNAM )
            NLAYS3D = EMLAYS
            IF( ALLOCATED( VGLVS ) ) THEN
                J = LBOUND( VGLVS3D,1 )
                DO V = 0, EMLAYS
                    VGLVS3D( J ) = VGLVS( V )
                    J = J + 1
                END DO
            ELSE
                VGLVS3D = 0  ! array
            ENDIF
            FDESC3D( 1 ) = 'Point source emissions data'
            IF( FDFLAG ) FDESC3D( 1 ) = METADESC

!.................  Open by logical name or physical name
            FILEDESC = 'POINT-SOURCE GRIDDED OUTPUT file'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC,'NETCDF',PONAME,I )

        END IF

        IF( XFLAG ) THEN
            CALL SETUP_VARIABLES( NIPPA, NMSPC, EANAM, EMNAM )
            NLAYS3D = EMLAYS
            IF( ALLOCATED( VGLVS ) ) THEN
                J = LBOUND( VGLVS3D,1 )
                DO V = 0, EMLAYS
                    VGLVS3D( J ) = VGLVS( V )
                    J = J + 1
                END DO
            ELSE
                VGLVS3D = 0  ! array
            ENDIF
            FDESC3D( 1 ) = 'Multiple category emissions data'
            IF( FDFLAG ) FDESC3D( 1 ) = METADESC

!.................  Open by logical name or physical name
            FILEDESC = 'MULTI-SOURCE GRIDDED OUTPUT file'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC,'NETCDF',TONAME,I )

        END IF

    END IF  ! End of gridded output

!.........  Open plume-in-grid output
    IF( PINGFLAG ) THEN

!.............  Get variable names for I/O API file
        CALL SETUP_VARIABLES( PNIPOL, PNMSPC, PEINAM, PEMNAM )

!.............  Override gridded file settings
        NCOLS3D = 1
        NROWS3D = NGROUP
        NLAYS3D = 1
        GDTYP3D = GDTYP
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

        FDESC3D = ' '   ! array

        PINGNAME = PROMPTSET(&
        &               'Enter name for PING EMISSIONS OUTPUT file',&
        &               FSUNKN3, PINGNAME, PROGNAME )
    END IF


!.........  Open plume-in-grid output
    IF( INLINEFLAG ) THEN

!.............  Get variable names for I/O API file
        CALL SETUP_VARIABLES( PNIPOL, PNMSPC, PEINAM, PEMNAM )

!............. Manual correction for HFLUX units for fires

        DO V = 1, NVARSET
            IF (VNAMESET(V) == 'HFLUX') THEN
                VUNITSET(V) = 'BTU/hr'
            END IF
        END DO

!.............  Override gridded file settings
        NCOLS3D = 1
        NROWS3D = NGROUP
        NLAYS3D = 1
        GDTYP3D = GDTYP
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

        FDESC3D = ' '   ! array

        INLINENAME = PROMPTSET(&
        &               'Enter name for INLINE EMISSIONS OUTPUT file',&
        &               FSUNKN3, INLINENAME, PROGNAME )
    END IF

!.........  Open plume-in-grid output
    IF( ELEVFLAG .AND. NGRP .GT. 1 ) THEN

        WRITE( MESG, 94010 ) 'The number of processing groups is ',&
        &       NGRP, ', but it cannot be greater than' // CRLF() //&
        &       BLANK10 // '1 for ASCII elevated output.  Not ' //&
        &       'enough memory.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ELSE IF( ELEVFLAG ) THEN

        EVDEV = PROMPTFFILE(&
        &        'Enter name for ASCII ELEVATED SOURCES file',&
        &        .FALSE., .TRUE., PELVNAME, PROGNAME )

    END IF

!.........  Open source apportionment output files
    IF( SRCGRPFLAG ) THEN

!.............  Override gridded file settings
        NCOLS3D = 1
        NROWS3D = NSGOUTPUT
        NLAYS3D = 1
        GDTYP3D = GDTYP
        VGTYP3D = IMISS3
        VGTOP3D = BADVAL3

        FDESC3D = ' '   ! array
        IF( AFLAG ) THEN
            FDESC3D( 1 ) = 'Area source groups file'
        ELSE IF( BFLAG ) THEN
            FDESC3D( 1 ) = 'Biogenic source groups file'
        ELSE IF( MFLAG ) THEN
            FDESC3D( 1 ) = 'Mobile source groups file'
        ELSE IF( PFLAG ) THEN
            FDESC3D( 1 ) = 'Point source groups file'
        END IF
        IF( FDFLAG ) FDESC3D( 1 ) = METADESC
        FDESC3D( 2 ) = '/FROM/ ' // PROGNAME
        FDESC3D( 3 ) = '/VERSION/ ' // VERCHAR( CVSW )
        WRITE( FDESC3D(5), 94010 ) '/NCOLS3D/ ', NCOLS
        WRITE( FDESC3D(6), 94010 ) '/NROWS3D/ ', NROWS

!.............  Build list of variables for stack groups file
        J = 1
        VNAME3D( J ) = 'IGROUP'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Source group number'

        J = J + 1
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
        VNAME3D( J ) = 'GRPCNT'
        VTYPE3D( J ) = M3INT
        UNITS3D( J ) = 'none'
        VDESC3D( J ) = 'Number of sources in group'

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
        UNITS3D( J ) = ''
        VDESC3D( J ) = 'Projection x coordinate'

        J = J + 1
        VNAME3D( J ) = 'YLOCA'
        VTYPE3D( J ) = M3REAL
        UNITS3D( J ) = ''
        VDESC3D( J ) = 'Projection y coordinate'

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

        IF( SGFIREFLAG ) THEN
            J = J + 1
            VNAME3D( J ) = 'ACRESBURNED'
            VTYPE3D( J ) = M3REAL
            UNITS3D( J ) = 'acres/day'
            VDESC3D( J ) = 'number of acres burned for a fire in one day'
        END IF

        NVARS3D = J

        SRCGRPNAME = PROMPTMFILE(&
        &               'Enter name for OUTPUT SOURCE GROUPS file',&
        &               FSUNKN3, 'SRCGROUPS_OUT', PROGNAME )

!.............  Set up variables for emissions output file
        CALL SETUP_VARIABLES( NIPPA, NMSPC, EANAM, EMNAM )

        SGINLNNAME = PROMPTSET(&
        &               'Enter name for INLINE EMISSIONS OUTPUT file',&
        &               FSUNKN3, SGINLNNAME, PROGNAME )
    END IF

!.........  Open report file(s)

    IF( LREPSTA .OR. LREPCNY ) THEN

        IF( LREPSTA .AND. LREPCNY ) THEN
            RTYPNAM = 'STATE AND COUNTY'
        ELSE IF ( LREPSTA ) THEN
            RTYPNAM = 'STATE'
        ELSE IF ( LREPCNY ) THEN
            RTYPNAM = 'COUNTY'
        END IF

        L = LEN_TRIM( RTYPNAM )

        IF( AFLAG ) THEN
            FILEDESC = 'AREA-SOURCE ' // RTYPNAM( 1:L ) // ' REPORT'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'ASCII',&
            &                          AREPNAME, ARDEV )
        END IF

        IF( BFLAG ) THEN
            FILEDESC = 'BIOGENIC ' // RTYPNAM( 1:L ) // ' REPORT'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'ASCII',&
            &                          BREPNAME, BRDEV )
        END IF

        IF( MFLAG ) THEN
            FILEDESC = 'MOBILE-SOURCE '// RTYPNAM( 1:L )// ' REPORT'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'ASCII',&
            &                          MREPNAME, MRDEV )
        END IF

        IF( PFLAG ) THEN
            FILEDESC = 'POINT-SOURCE '// RTYPNAM( 1:L )// ' REPORT'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'ASCII',&
            &                          PREPNAME, PRDEV )
        END IF

        IF( XFLAG ) THEN
            FILEDESC = 'TOTAL '// RTYPNAM( 1:L )// ' REPORT'
            CALL OPEN_LNAME_OR_PNAME( FILEDESC, 'ASCII',&
            &                          TREPNAME, TRDEV )
        END IF

    END IF  ! End of state and/or county output

!.........  Abort if problem with output files.
    IF( EFLAG ) THEN

        MESG = 'Problem opening output file(s)'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats.............94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94050 FORMAT( A, 1X, I2.2, A, 1X, A, 1X, I6.6, 1X,&
    &        A, 1X, I3.3, 1X, A, 1X, I3.3, 1X, A   )

!*****************  INTERNAL SUBPROGRAMS  ******************************

CONTAINS

!.............  This internal subprogram uses WRITE3 and exits gracefully
!               if a write error occurred
    SUBROUTINE SETUP_VARIABLES( NIPPA_L, NMSPC_L,&
    &                            EANAM_L, EMNAM_L  )

!.............  MODULES for public variables
!.............  This module contains the major data structure and control flags
        USE MODMERGE, ONLY: SFLAG, NIPOL, EINAM

!.............  Internal subprogram arguments
        INTEGER     , INTENT (IN) :: NIPPA_L
        INTEGER     , INTENT (IN) :: NMSPC_L
        CHARACTER(*), INTENT (IN) :: EANAM_L( NIPPA_L )
        CHARACTER(*), INTENT (IN) :: EMNAM_L( NMSPC_L )

!.............  Local subprogram varibles
        INTEGER     I, J, K, LJ, M, V

        INTEGER     IOS    ! allocate status

        CHARACTER(NAMLEN3) CBUF

!------------------------------------------------------------------------

!.............  Set constants number and values for variables
!.............  Do this regardless of whether we have outputs or not
!.............  For speciation...
        IF( SFLAG ) THEN
            NVARSET = NMSPC_L

            IF( ALLOCATED( VTYPESET ) )&
            &    DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )
            IF( ALLOCATED( VARS_PER_FILE ) )&
            &    DEALLOCATE( VARS_PER_FILE )
            ALLOCATE( VTYPESET( NVARSET ),&
            &          VNAMESET( NVARSET ),&
            &          VUNITSET( NVARSET ),&
            &          VDESCSET( NVARSET ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )

            K = 0
            LJ = -1
!.................  Loop through global species index
            DO N = 1, NGRP
                DO V = 1, VGRPCNT( N )

!.........................  Access global indices
                    I = SIINDEX( V,N )
                    IF( I .EQ. 0 ) CYCLE     ! Skip unused pollutants
                    J = SPINDEX( V,N )
                    IF( J .EQ. LJ ) CYCLE    ! Do not repeat species

!.........................  Make sure current species is in local array
                    CBUF = EMNAM( J )
                    M = INDEX1( CBUF, NMSPC_L, EMNAM_L )
                    IF( M .LE. 0 ) CYCLE

                    DESCBUF= DESCBUF(1:LD)//' Model species '// CBUF

                    K = K + 1
                    VNAMESET( K ) = CBUF
                    VUNITSET( K ) = GRDUNIT( J )
                    VDESCSET( K ) = ADJUSTL( DESCBUF )
                    VTYPESET( K ) = M3REAL

                    LJ = J

                END DO
            END DO

!.............  For no speciation...
        ELSE

            NVARSET = NIPPA_L

            IF( ALLOCATED( VTYPESET ) )&
            &    DEALLOCATE( VTYPESET, VNAMESET, VUNITSET, VDESCSET )
            IF( ALLOCATED( VARS_PER_FILE ) )&
            &    DEALLOCATE( VARS_PER_FILE )
            ALLOCATE( VTYPESET( NVARSET ),&
            &          VNAMESET( NVARSET ),&
            &          VUNITSET( NVARSET ),&
            &          VDESCSET( NVARSET ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VDESCSET', PROGNAME )

            K = 0
!.................  Loop through global list
            DO V = 1, NIPPA

!.....................  Determine if variable is an emission type, pollutant,
!                       or activity, and set variable description accordingly
                I = INDEX1( EANAM( V ), NIPOL, EINAM )
                J = INDEX ( EANAM( V ), ETJOIN )
                IF( J .GT. 0 ) THEN
                    DESCBUF = DESCBUF(1:LD) //&
                    &         ' Emission type ' // EANAM( V )
                ELSE IF( I .GT. 0 ) THEN
                    DESCBUF = DESCBUF(1:LD) //&
                    &          ' Pollutant ' // EANAM( V )
                ELSE
                    DESCBUF = DESCBUF(1:LD) //&
                    &          ' Activity ' // EANAM( V )
                END IF

!.....................  Make sure current data variable is in local array
                CBUF = EANAM( V )
                M = INDEX1( CBUF, NIPPA_L, EANAM_L )
                IF( M .LE. 0 ) CYCLE

                K = K + 1

!.....................  Define variable information
                VNAMESET( K ) = EANAM  ( V )
                VUNITSET( K ) = GRDUNIT( V )
                VDESCSET( K ) = ADJUSTL( DESCBUF )
                VTYPESET( K ) = M3REAL

            END DO

        END IF

        RETURN

    END SUBROUTINE SETUP_VARIABLES

!------------------------------------------------------------------------
!------------------------------------------------------------------------

!.............  Routine to open by physical name and give error if problem
    SUBROUTINE OPEN_LNAME_OR_PNAME( FILEDESC, FTYPE,&
    &                                LNAME, FDEV      )

!.............  Subroutine arguments
        CHARACTER(*), INTENT  ( IN ) :: FILEDESC  ! file description
        CHARACTER(*), INTENT  ( IN ) :: FTYPE     ! NETCDF or ASCII
        CHARACTER(*), INTENT(IN OUT) :: LNAME     ! logical file name
        INTEGER     , INTENT   (OUT) :: FDEV      ! file unit (if any)

!.............  Local variables
        INTEGER         L, L2
        INTEGER         IOS

        CHARACTER(128)  DUMSTR       ! dummy string
        CHARACTER(256)  MESG         ! output string buffer
        CHARACTER(512)  PHYSNAME     ! output path & file buffer

        LOGICAL ::      ENVFLAG = .FALSE. ! true: could not set env variable

!------------------------------------------------------------------------

        FDEV = 0
        ENVFLAG = .FALSE.

!.............  Check if logical output file name is defined
        CALL ENVSTR( LNAME, FILEDESC, ' ', DUMSTR, IOS )

!.............  If not defined, give note, build physical output file name
!              and set environment variable
        IF( IOS1 .EQ. 0 .AND.&
        &    IOS2 .EQ. 0 .AND.&
        &    IOS  .NE. 0       ) THEN

            PHYSNAME = OUTDIR( 1:L3 ) // '/' // TRIM( LNAME ) //&
            &           '.' // OUTSCEN( 1:L4 ) // '.ncf'

            MESG = 'NOTE: Opening "' // TRIM( LNAME ) //&
            &       '" file using Smkmerge-defined physical ' //&
            &       'file name.'
            CALL M3MSG2( MESG )

            IF( .NOT. SETENVVAR( LNAME, PHYSNAME ) ) THEN
                ENVFLAG = .TRUE.
                MESG = 'ERROR: Could not set logical file name ' //&
                &       CRLF() // BLANK10 // 'for file ' //&
                &       TRIM( PHYSNAME )
                CALL M3MSG2( MESG )
            END IF
        END IF

        IF( .NOT. ENVFLAG ) THEN

!.............  Logical name is defined, open file.
            SELECT CASE( FTYPE )
              CASE( 'NETCDF' )
                LNAME = PROMPTSET( 'Enter name for '// FILEDESC,&
                &                    FSUNKN3, LNAME, PROGNAME    )
              CASE( 'ASCII' )
                FDEV  = PROMPTFFILE(&
                &        'Enter name for  '// FILEDESC,&
                &       .FALSE., .TRUE., LNAME, PROGNAME )

            END SELECT

        ELSE

            EFLAG = .TRUE.

        END IF

        RETURN

    END SUBROUTINE OPEN_LNAME_OR_PNAME

END SUBROUTINE OPENMRGOUT
