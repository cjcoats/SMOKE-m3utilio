
PROGRAM INLINETO2D

!***********************************************************************
!  program body starts at line
!
!  DESCRIPTION:
!    Program INLINETO2D reads STACK_GROUPS and INLN I/O API files and
!    converts them into a 2-D emission file for QA
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Original by G. Pouliot 2/25/2010
!       Revised 11/16/2011 to use GRIDDESC information for output file
!       and ignore grid info from input files
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!***********************************************************************
!
! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
!                System
! File: @(#)$Id: inlineto2d.f,v 1.19 2007/07/11 19:30:43 bbaek Exp $
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
!***********************************************************************
    USE M3UTILIO

!.........  MODULES for public variables
!.........  This module contains the global variables for the 3-d grid
    USE MODSOURC, ONLY: XLOCA, YLOCA

    USE MODGRID, ONLY: GRDNM, COORD, GDTYP, P_ALP, P_BET, P_GAM,&
    &                   XCENT, YCENT, XORIG, YORIG, XCELL, YCELL,&
    &                   NCOLS, NROWS, NGRID

    USE MODGRDLIB

    IMPLICIT NONE

!...........   INCLUDES:
    INCLUDE 'EMCNST3.EXT'

!...........   EXTERNAL FUNCTIONS
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: DSCM3GRD

!.........  LOCAL PARAMETERS and their descriptions:

    CHARACTER(16), PARAMETER :: PROGNAME = 'INLINETO2D' ! program name
    CHARACTER(50), PARAMETER ::&
    &CVSW = '$Name: SMOKEv5.0_Jun2023$' ! CVS release tag

!...........   LOCAL VARIABLES and their descriptions:
!...........   Emissions arrays
    REAL,      ALLOCATABLE :: EIN ( : ,: )
    REAL,      ALLOCATABLE :: EOUT( : ,:, :)
    REAL,      ALLOCATABLE :: LAT(:), LON(:)
    INTEGER,   ALLOCATABLE :: LMAJOR(:)
    INTEGER,   ALLOCATABLE :: ROW(:), COL(:), ISTACK(:)
    INTEGER,   ALLOCATABLE :: NEW_COL(:), NEW_ROW(:)

!...........   Input file descriptors
    INTEGER       DURATA, DURATB ! no. time steps
    INTEGER       NCOLSA, NCOLSB! no. columns file a, file b, output file
    INTEGER       NROWSA, NROWSB  ! no. rows
    INTEGER       NVARSA, NVARSB ! no. variables
    INTEGER       SDATEA, SDATEB   ! start date
    INTEGER       STIMEA, STIMEB   ! start time
    INTEGER       NLAYSA, NLAYSB ! number of layers in the file
    INTEGER       NCOLS_N, NROWS_N

    REAL, ALLOCATABLE :: EMIS_SUM(:,:)
    CHARACTER(16)   VNAMEA( MXVARS3 ),  VNAMEB( MXVARS3 ) ! variable names
    CHARACTER(16)   VUNITA( MXVARS3 ),  VUNITB( MXVARS3 ) ! variable units
    CHARACTER(80)   VDESCA( MXVARS3 ),  VDESCB( MXVARS3 ) ! var descrip
    INTEGER         VTYPEA( MXVARS3 ),  VTYPEB( MXVARS3 ) ! var type

!...........   Other local variables
    INTEGER       C, F, J, K, L, L1, L2, NL, V, T, S ! pointers and counters
    INTEGER       ISTART, IEND, ICT,V_R, V_I
    INTEGER       J_LOOP
    INTEGER       DUMMY                      ! dummy value for use
    INTEGER       EDATE                      ! ending julian date
    INTEGER       ETIME                      ! ending time HHMMSS
    INTEGER       NSRC                       ! no of sources
    INTEGER       IOS                        ! i/o status
    INTEGER       IREC                       ! line number count
    INTEGER       JDATE                      ! iterative julian date
    INTEGER       JTIME                      ! iterative time HHMMSS
    INTEGER       LB                         ! leading blanks counter
    INTEGER       LE                         ! location of end of string
    INTEGER       MXNFIL_1                   ! max no. of stack groups files
    INTEGER       MXNFIL_2                   ! max no. of emission files
    INTEGER       MXNFIL                     ! max no. of both
    INTEGER       NFILE                      ! no. of 2-d input files
    INTEGER       NSTEPS                     ! no. of output time steps
    INTEGER       NVOUT(2)                   ! no. of output variables
    INTEGER       NVOUT_I(2)                 ! no. of integer type output variables
    INTEGER       NVOUT_R(2)                 ! no. of real type output variables
    INTEGER       RDATE                      ! reference date
    INTEGER       SAVLAYS                    ! number of layers
    INTEGER       SDATE                      ! starting julian date
    INTEGER       SECS                       ! tmp seconds
    INTEGER       SECSMAX                    ! seconds maximum
    INTEGER       SECSMIN                    ! seconds minimum
    INTEGER       STIME                      ! starting time HHMMSS
    INTEGER       STEPS                      ! tmp number of steps
    INTEGER       TIMET                      ! tmp time from seconds
    INTEGER       TSTEP                      ! time step
    INTEGER       VLB                        ! VGLVS3D lower bound
    INTEGER       TMAX                       ! max time counter
    INTEGER       LDEV, S_CT

    CHARACTER(16)  ONAME                     ! ONAME
    CHARACTER(16)  FDESC                     ! tmp file description
    CHARACTER(16)  NAM                       ! tmp file name
    CHARACTER(16)  VNM                       ! tmp variable name
    CHARACTER(256) LINE                      ! input buffer
    CHARACTER(256) MESG                      ! message field
    CHARACTER(15)  RPTCOL                    ! single column in report line
    CHARACTER(300) RPTLINE                   ! line of report file

    LOGICAL    :: EFLAG   = .FALSE.   ! error flag
    LOGICAL    :: FIRST3D = .TRUE.    ! true: first 3-d file not yet input
    LOGICAL    :: LFLAG   = .FALSE.   ! true iff 3-d file input
    LOGICAL    :: TFLAG   = .FALSE.   ! true: grid didn't match
    LOGICAL    :: LATLON_IN_FILE

    CHARACTER(80) :: GDESC        ! grid description
    CHARACTER(NAMLEN3) COORD3D    ! coordinate system name
    CHARACTER(NAMLEN3) COORUN3D   ! coordinate system projection units
    CHARACTER(80) :: GDNAM

!***********************************************************************
!   begin body of program INLINETO2D

    LDEV = INIT3()

!.........  Write out copyright, version, web address, header info, and prompt
!           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    LATLON_IN_FILE = .TRUE.

    IF ( .NOT. OPEN3( 'INLN', FSREAD3, PROGNAME )) THEN
        MESG = 'Could not open file INLN '
        CALL M3MSG2( MESG )
        MESG = 'Ending program "INLINETO2D".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. DESC3( 'INLN' ) ) THEN
        MESG = 'Could not get description of file INLN '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        NROWSA = NROWS3D
        NCOLSA = NCOLS3D
        NLAYSA = NLAYS3D
        NVARSA = NVARS3D
        SDATEA = SDATE3D
        STIMEA = STIME3D
        DURATA = MXREC3D
        TSTEP = TSTEP3D

        DO V = 1, NVARSA
            VNAMEA( V ) = VNAME3D( V )
            VUNITA( V ) = UNITS3D( V )
            VDESCA( V ) = VDESC3D( V )
            VTYPEA( V ) = VTYPE3D( V )
        END DO
    END IF

    NSRC = NROWSA
    ALLOCATE(      EIN( NSRC, NVARSA),&
    &          EMIS_SUM(NVARSA,3), STAT=IOS )
    CALL CHECKMEM( IOS, 'EIN', PROGNAME )

    IF ( .NOT. OPEN3( 'STACK_GROUPS', FSREAD3, PROGNAME )) THEN
        MESG = 'Could not open file INLN '
        CALL M3MSG2( MESG )
        MESG = 'Ending program "INLINETO2D".'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF ( .NOT. DESC3( 'STACK_GROUPS' ) ) THEN
        MESG = 'Could not get description of file STACK_GROUPS '
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        NROWSB = NROWS3D
        NCOLSB = NCOLS3D
        NLAYSB = NLAYS3D
        NVARSB = NVARS3D
        SDATEB = SDATE3D
        STIMEB = STIME3D
        DURATB = MXREC3D

        DO V = 1, NVARSB
            VNAMEB( V ) = VNAME3D( V )
            VUNITB( V ) = UNITS3D( V )
            VDESCB( V ) = VDESC3D( V )
            VTYPEB( V ) = VTYPE3D( V )
        END DO

    END IF

    ALLOCATE( COL( NSRC ),&
    &          ROW( NSRC ),&
    &        XLOCA( NSRC ),&
    &        YLOCA( NSRC ),&
    &       ISTACK( NSRC ),&
    &          LAT( NSRC ),&
    &          LON( NSRC ),&
    &       LMAJOR( NSRC ),&
    &      NEW_COL( NSRC ),&
    &      NEW_ROW( NSRC ), STAT=IOS )
    CALL CHECKMEM( IOS, 'COL...NEW_ROW', PROGNAME )

    IF( .NOT. READ3( 'STACK_GROUPS', 'COL', 1, SDATEB, STIMEB, COL)) THEN
        MESG = 'Could not read  COL from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

    IF( .NOT. READ3( 'STACK_GROUPS', 'ROW', 1, SDATEB, STIMEB, ROW)) THEN
        MESG = 'Could not read  ROW from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

    IF( .NOT. READ3( 'STACK_GROUPS', 'XLOCA', 1,SDATEB, STIMEB, XLOCA)) THEN
        MESG = 'Could not read  XLOCA from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

    IF( .NOT. READ3( 'STACK_GROUPS', 'YLOCA', 1, SDATEB, STIMEB, YLOCA)) THEN
        MESG = 'Could not read  YLOCA from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

    IF (LATLON_IN_FILE) THEN
        IF( .NOT. READ3( 'STACK_GROUPS', 'LATITUDE', 1, SDATEB, STIMEB, LAT)) THEN
            MESG = 'Could not read  LAT from STACK_GROUPS'
            CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF

        IF( .NOT. READ3( 'STACK_GROUPS', 'LONGITUDE', 1, SDATEB, STIMEB, LON)) THEN
            MESG = 'Could not read  LON from STACK_GROUPS'
            CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
        ENDIF
    ENDIF ! latlon in file check

    IF( .NOT. READ3( 'STACK_GROUPS', 'ISTACK', 1, SDATEB, STIMEB, ISTACK)) THEN
        MESG = 'Could not read  ISTACK from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

    IF( .NOT. READ3( 'STACK_GROUPS', 'LMAJOR', 1, SDATEB, STIMEB, LMAJOR)) THEN
        MESG = 'Could not read  LMAJOR from STACK_GROUPS'
        CALL M3EXIT( PROGNAME, 0, 0,MESG, 2 )
    ENDIF

!***** GET GRID INFORMATION FOR OUTPUT FILE AND OPEN OUTPUT FILE
    IF( .NOT. DSCM3GRD( GDNAM3D, GDESC, COORD3D, GDTYP3D, COORUN3D,&
    &              P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,&
    &              XORIG3D, YORIG3D, XCELL3D, YCELL3D,&
    &              NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ELSE
        GDTYP = GDTYP3D
        P_ALP = P_ALP3D
        P_BET = P_BET3D
        P_GAM = P_GAM3D
        XCENT = XCENT3D
        YCENT = YCENT3D
        NCOLS = NCOLS3D
        NROWS = NROWS3D
        NGRID = NCOLS * NROWS

    END IF

    IF (LATLON_IN_FILE) THEN
        XLOCA(1:NSRC) = LON(1:NSRC)
        YLOCA(1:NSRC) = LAT(1:NSRC)

        GDNAM = GDNAM3D

        CALL CONVRTXY( NSRC, GDTYP, GDNAM, P_ALP, P_BET, P_GAM,&
        &               XCENT, YCENT, XLOCA, YLOCA )

        DO S = 1, NSRC
            IF( .NOT. INGRID( XLOCA(S), YLOCA(S), NCOLS_N, NROWS_N,&
            &                  COL(S), ROW(S) ) ) THEN
                CYCLE  ! To end of loop
            ENDIF
        ENDDO

    ENDIF    ! latlon in file check

    ALLOCATE (EOUT(NCOLS, NROWS, NVARSA))
    EOUT(1:NCOLS,1:NROWS,1:NVARSA) = 0.0

!.........  Set up layer structure for output file
    SDATE3D = SDATEA
    STIME3D = STIMEA
    NVARS3D = NVARSA
    GDTYP3D = GDTYP
    P_ALP3D = P_ALP
    P_BET3D = P_BET
    P_GAM3D = P_GAM
    XCENT3D = XCENT
    YCENT3D = YCENT
    NCOLS3D = NCOLS
    NROWS3D = NROWS

!.........  Set up layer structure for output file
    NLAYS3D = 1
    VGLVS3D = 0.     ! initialize array
    DO V = 1, NVARSA
        VNAME3D( V ) = VNAMEA( V )
        UNITS3D( V ) = VUNITA( V )
        VDESC3D( V ) = VDESCA( V )
        VTYPE3D( V ) = VTYPEA( V )
    END DO

    ONAME = PROMPTMFILE(&
    &        'Enter logical name for OUTPUT file',&
    &        FSUNKN3, 'OUTFILE', PROGNAME )

    JDATE = SDATEA
    JTIME = STIMEA
    TMAX  = DURATA

    EMIS_SUM(1:NVARSA,1:3) = 0.0

    DO T = 1, TMAX

        EOUT(1:NCOLS,1:NROWS,1:NVARSA) = 0.0

        DO V = 1, NVARSA
            IF( .NOT. READ3( 'INLN', VNAMEA(V), 1, JDATE,&
            &                  JTIME, EIN(1,V))) THEN
                MESG = 'Could not read "' // VNAMEA(V) //&
                &       '" from file INLN "'
                CALL M3EXIT( PROGNAME, RDATE, JTIME, MESG, 2 )
            ENDIF

            S_CT = 0
            DO S = 1, NSRC
                IF (COL(S) .lt. 1) CYCLE
                IF (ROW(S) .lt. 1) CYCLE
                IF (COL(S) .gt. ncols) CYCLE
                IF (ROW(S) .gt. nrows) CYCLE
                S_CT = S_CT + 1
                EOUT(COL(S),ROW(S),V) = EOUT(COL(S),ROW(S),V) + EIN(S,V)
                EMIS_SUM(V,3) = EMIS_SUM(V,3) + EIN(S,V)
            END DO

            EMIS_SUM(V,1) = EMIS_SUM(V,1) + SUM(EIN(1:NSRC,V))
            EMIS_SUM(V,2) = EMIS_SUM(V,2) + SUM(EOUT(1:NCOLS,1:NROWS,V))

            IF( .NOT. WRITE3( ONAME, VNAMEA(V), JDATE,&
            &                  JTIME, EOUT(1,1,V))) THEN
                MESG = 'Could not write "' // VNAMEA(V) //&
                &       '" to file OUTFILE "'
                CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )
            ENDIF

        END DO   ! loop through variables

        CALL NEXTIME( JDATE, JTIME, TSTEP )

    END DO       ! loop through timesteps

    CALL M3EXIT( PROGNAME, 0, 0, ' ', 0)

!******************  FORMAT  STATEMENTS   ******************************

!...........   Informational (LOG) message formats... 92xxx

92000 FORMAT( 5X, A )

!...........   Formatted file I/O formats............ 93xxx


93000 FORMAT(  A )

93010 FORMAT( A15 )

93020 FORMAT( I15 )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I7, :, 1X ) )

94020 FORMAT( A, :, I3, :, 1X, 10 ( A, :, F8.5, :, 1X ) )

END PROGRAM INLINETO2D
