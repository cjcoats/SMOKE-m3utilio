
LOGICAL FUNCTION DSCM3GRD( GNAME, GDESC, CNAME, CTYPE, PUNIT,       &
                           P_ALP, P_BET, P_GAM, XCENT, YCENT,       &
                           XORIG, YORIG, XCELL, YCELL, NCOLS,       &
                           NROWS, NTHIK )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns true it successfully reads the properties of the
    !      models-3 grid, and false otherwise.  It also returns those properties
    !      using the subroutine arguments.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !       I/O API routines DSCGRID()
    !
    !  REVISION  HISTORY:
    !       Copied and modified from DSCGRID, version 12/16/97
    !       Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",
    !           Bug-fixes in  , GRDINFO_CHECKER, GRDINFO_DEFINED and related changes
    !**************************************************************************
    !
    ! Project Title: EDSS Tools Library
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

    INCLUDE 'IOCNST3.h90'
    INCLUDE 'IOSTRG3.h90'

    !.......   ARGUMENTS and their descriptions.  All are output variables:

    CHARACTER(*), INTENT(  OUT) :: GNAME         !  grid name
    CHARACTER(*), INTENT(  OUT) :: GDESC         !  grid description
    CHARACTER(*), INTENT(  OUT) :: CNAME         !  coord sys name
    INTEGER     , INTENT(  OUT) :: CTYPE         !  coord sys type (I/O API code number)
    CHARACTER(*), INTENT(  OUT) :: PUNIT         !  projection units (e.g., meters or degrees)
    REAL(8)     , INTENT(  OUT) :: P_ALP         !  first, second, third map
    REAL(8)     , INTENT(  OUT) :: P_BET         !  projection descriptive
    REAL(8)     , INTENT(  OUT) :: P_GAM         !  parameters
    REAL(8)     , INTENT(  OUT) :: XCENT         !  lon for coord-system X=0
    REAL(8)     , INTENT(  OUT) :: YCENT         !  lat for coord-system Y=0
    REAL(8)     , INTENT(  OUT) :: XORIG         !  X-coordinate origin of grid (map units)
    REAL(8)     , INTENT(  OUT) :: YORIG         !  Y-coordinate origin of grid
    REAL(8)     , INTENT(  OUT) :: XCELL         !  X-coordinate cell dimension
    REAL(8)     , INTENT(  OUT) :: YCELL         !  Y-coordinate cell dimension
    INTEGER     , INTENT(  OUT) :: NCOLS         !  number of grid columns
    INTEGER     , INTENT(  OUT) :: NROWS         !  number of grid rows
    INTEGER     , INTENT(  OUT) :: NTHIK         !  BOUNDARY:  perimeter thickness (cells)

    !.......   EXTERNAL FUNCTIONS:

    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL
    INTEGER, EXTERNAL :: GETNLIST

    !.......   Local parameters
    INTEGER, PARAMETER :: MXGRDTYP = 13

    !.......   Grid types and names arrays
    INTEGER      :: GRDTYPES( MXGRDTYP ) = (/  LATGRD3    &
                                             , LATGRD3    &
                                             , LATGRD3    &
                                             , LAMGRD3    &
                                             , LAMGRD3    &
                                             , MERGRD3    &
                                             , MERGRD3    &
                                             , STEGRD3    &
                                             , STEGRD3    &
                                             , POLGRD3    &
                                             , POLGRD3    &
                                             , UTMGRD3    &
                                             , UTMGRD3 /)

    CHARACTER(15):: GRDNAMES( MXGRDTYP ) = (/  'LAT-LON        '    &
                                             , 'GEOGRAPHIC     '    &
                                             , 'LATGRD3        '    &
                                             , 'LAMBERT        '    &
                                             , 'LAMGRD3        '    &
                                             , 'MERCATOR       '    &
                                             , 'MERGRD3        '    &
                                             , 'STEREOGRAPHIC  '    &
                                             , 'STEGRD3        '    &
                                             , 'POLAR          '    &
                                             , 'POLGRD3        '    &
                                             , 'UTM            '    &
                                             , 'UTMGRD3        ' /)

    !.......   Local arrays (note- tried to make these allocatable, but
    !              this caused unexplainable crashing on SGI).
    CHARACTER(32) :: SEGMENT( 32 )
    CHARACTER(32) :: UPCSGMT( 32 )

    !.......   File units and logical/physical names:
    INTEGER, SAVE :: IDEV        !  unit number of grid information file
    CHARACTER(16), SAVE :: LNAME     !  logical name for grid information file

    !.......   Scratch local variables and their descriptions:

    INTEGER         I, J, L, L2, N      !  indices and string lengths
    INTEGER         IOS          !  I/O status return
    INTEGER         IOS2         !  I/O status return for ENVSTR setting
    INTEGER      :: IDUM = 0     !  integer dummy variable
    INTEGER         IREC         !  record counter
    INTEGER         LEGAL        !  valid length of string
    INTEGER         NDEV         !  next file number to be used

    REAL(8)      :: RDUM = 0     !  double precision dummy variable

    LOGICAL         EFLAG                 ! true: error detected
    LOGICAL      :: FIRSTIME = .TRUE.     ! true: first time routine called
    LOGICAL,SAVE :: GRIDPATH = .TRUE.     ! true: for G_GRIDPATH fmt or 1st read

    CHARACTER(16)   CDUM         !  character dummy buffer
    CHARACTER(300)  BUFFER       !  multi-purpose buffer
    CHARACTER(300)  GNBUF        !  grid name buffer
    CHARACTER(300)  GDBUF        !  grid description buffer
    CHARACTER(300)  LINE         !  line of file
    CHARACTER(300)  MESG         !  message buffer
    CHARACTER(300)  UPCLINE      !  upper case line of file

    CHARACTER(16), PARAMETER :: PROGNAME = 'DSCM3GRD'     ! program name

    !....  single grid file format intergrity flags
    LOGICAL      :: GNAME_FL = .FALSE.     !  true: found GNAME
    LOGICAL      :: GDESC_FL = .FALSE.     !  true: found GDESC
    LOGICAL      :: CTYPE_FL = .FALSE.     !  true: found CTYPE
    LOGICAL      :: PUNIT_FL = .FALSE.     !  true: found PUNIT
    LOGICAL      :: P_ALP_FL = .FALSE.     !  true: found P_ALP
    LOGICAL      :: P_BET_FL = .FALSE.     !  true: found P_BET
    LOGICAL      :: P_GAM_FL = .FALSE.     !  true: found P_GAM
    LOGICAL      :: XCENT_FL = .FALSE.     !  true: found XCENT
    LOGICAL      :: YCENT_FL = .FALSE.     !  true: found YCENT
    LOGICAL      :: XORIG_FL = .FALSE.     !  true: found XORIG
    LOGICAL      :: YORIG_FL = .FALSE.     !  true: found YORIG
    LOGICAL      :: XCELL_FL = .FALSE.     !  true: found XCELL
    LOGICAL      :: YCELL_FL = .FALSE.     !  true: found YCELL
    LOGICAL      :: NCOLS_FL = .FALSE.     !  true: found NCOLS
    LOGICAL      :: NROWS_FL = .FALSE.     !  true: found NROWS
    LOGICAL      :: NTHIK_FL = .FALSE.     !  true: found NTHIK
    LOGICAL         ALLFOUND     !  true: G_GRIDPATH format file detected

    !***********************************************************************
    !   begin body of function DSCM3GRD

    !.....  Get value of Models-3 environment variable for grid and check
    !           exit status

    EFLAG = .FALSE.
    IF ( FIRSTIME ) THEN
        MESG = 'Checking for defined environment variable'
        CALL ENVSTR( 'GRIDDESC', MESG, ' ', BUFFER, IOS )
        IF ( IOS .EQ. 0 ) THEN
            LNAME = 'GRIDDESC'
        ELSE
            LNAME = 'G_GRIDPATH'
        END IF
        FIRSTIME = .FALSE.
    END IF

    !.....  Try to read G_GRIDPATH file the first time or skip
    !           if we know this is a GRIDDESC file

    IF ( GRIDPATH ) THEN

        !.....  Open grid file
        IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

        IF( IDEV .LE. 0 ) THEN

            MESG = 'Could not open file "' // TRIM( LNAME ) // '"!'
            DSCM3GRD = .FALSE.
            RETURN

        END IF

        !.....  Initialize grid setting to missing
        GNAME   = ' '
        GDESC   = ' '
        CNAME   = ' '
        CTYPE   = IMISS3
        PUNIT   = ' '
        P_ALP   = DBLE( BADVAL3 )
        P_BET   = DBLE( BADVAL3 )
        P_GAM   = DBLE( BADVAL3 )
        XCENT   = DBLE( BADVAL3 )
        YCENT   = DBLE( BADVAL3 )
        XORIG   = DBLE( BADVAL3 )
        YORIG   = DBLE( BADVAL3 )
        XCELL   = DBLE( BADVAL3 )
        YCELL   = DBLE( BADVAL3 )
        NCOLS   = IMISS3
        NROWS   = IMISS3
        NTHIK   = IMISS3

        !.....  Make sure read it at the start of the file
        REWIND( IDEV )

        !.....  Read grid information file
        IREC = 0
        DO

            !.....  Read whole line from the file
            READ( IDEV, 93000, END=111, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF( IOS .GT. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,    &
                      'reading grid information file at line', IREC

                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.....  Skip blank lines
            IF( LINE .EQ. ' ' ) CYCLE

            !......  Adjust line to left and create upper case line
            LINE = ADJUSTL( LINE )
            UPCLINE = LINE
            CALL UPCASE( UPCLINE )
            L2 = LEN_TRIM( LINE )

            !.....  Skip line if it starts with a number
            IF( UPCLINE( 1:1 ) .GE. '0' .AND.           &
                UPCLINE( 1:1 ) .LE. '9'       ) CYCLE

            !.....  Parse the lower and upper case lines
            N = GETNLIST( L2, UPCLINE )
            IF( N .GT. 32 ) THEN
                WRITE( MESG,94010 ) 'WARNING: Line', IREC,              &
                       'in GRIDDESC file skipped because more ' //      &
                       'than 32 fields were found.'
                CALL M3MSG2( MESG )
                CYCLE
            END IF
            ! note: add allocatable array
            SEGMENT = ' '
            UPCSGMT = ' '
            CALL PARSLINE( LINE   , N, SEGMENT )
            CALL PARSLINE( UPCLINE, N, UPCSGMT )

            !.....  Search for keywords.  If the keyword exists, extract the
            !                   value(s) and cycle to the next line...

            !.....  Grid name
            SELECT CASE ( UPCSGMT( 1 ) )
              CASE( 'GDNAME_GD' )
                GNBUF = UPCSGMT( 3 )
                GNAME_FL = .TRUE.

            !.....  Grid description
              CASE( 'GDDESC_GD' )
                I = INDEX( LINE, '0' )
                GDBUF = ADJUSTL( LINE( I+1:L2 ) )
                GDESC_FL = .TRUE.

            !.....  Coordinate system type
              CASE( 'GDTYP_GD' )
                CNAME = UPCSGMT( 3 )

            !.....  Look for coordinate system type in known types
                J = INDEX1( CNAME, MXGRDTYP, GRDNAMES )
                IF( J .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    MESG = 'ERROR: Coordinate system type "' //         &
                            CNAME // '" is not known in grid ' //       &
                            'information reader.'
                    CALL M3MSG2( MESG )
                ELSE

                    CTYPE = GRDTYPES( J )
                    CTYPE_FL = .TRUE.

                END IF

            !.....  Grid units
              CASE( 'GDUNT_GD' )
                PUNIT = SEGMENT( 3 )
                PUNIT_FL = .TRUE.

            !.....  Alpha value
              CASE( 'P_ALP_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_ALP, IDUM )
                P_ALP_FL = .TRUE.

            !.....  Beta value
              CASE( 'P_BET_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_BET, IDUM )
                P_BET_FL = .TRUE.

            !.....  Gamma value
              CASE( 'P_GAM_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, P_GAM, IDUM )
                P_GAM_FL = .TRUE.

            !.....  X-center value
              CASE( 'XCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XCENT, IDUM )
                XCENT_FL = .TRUE.

            !.....  Y-center value
              CASE( 'YCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YCENT, IDUM )
                YCENT_FL = .TRUE.

            !.....  X-center value
              CASE( 'P_XCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XCENT, IDUM )
                XCENT_FL = .TRUE.

            !.....  Y-center value
              CASE( 'P_YCENT_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YCENT, IDUM )
                YCENT_FL = .TRUE.

            !.....  X-origin value
              CASE( 'XORIG_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XORIG, IDUM )
                XORIG_FL = .TRUE.

            !.....  Y-origin value
              CASE( 'YORIG_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YORIG, IDUM )
                YORIG_FL = .TRUE.

            !.....  Delta X value
              CASE( 'XCELL_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, XCELL, IDUM )
                XCELL_FL = .TRUE.

            !.....  Delta Y value
              CASE( 'YCELL_GD' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3DBLE, YCELL, IDUM )
                YCELL_FL = .TRUE.

            !.....  Number of columns
              CASE( 'NCOLS' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NCOLS )
                NCOLS_FL = .TRUE.

            !.....  Number of rows
              CASE( 'NROWS' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NROWS )
                NROWS_FL = .TRUE.

            !.....  Number of boundary cells
              CASE( 'NTHIK' )
                CALL GRDINFO_CHECKER( UPCSGMT, M3INT, RDUM, NTHIK )
                NTHIK_FL = .TRUE.

            END SELECT

        END DO

    !.....  Exit from read loop
111     CONTINUE

    !.....  Close grid file
        CLOSE( IDEV )

    END IF      ! end G_GRIDPATH file read condition

    !.....  Check that all G_GRIDPATH file variables were found
    ALLFOUND = ( GNAME_FL .AND. GDESC_FL .AND. CTYPE_FL .AND.       &
                 PUNIT_FL .AND. P_ALP_FL .AND. P_BET_FL .AND.       &
                 P_GAM_FL .AND. XCENT_FL .AND. YCENT_FL .AND.       &
                 XORIG_FL .AND. YORIG_FL .AND. XCELL_FL .AND.       &
                 YCELL_FL .AND. NCOLS_FL .AND. NROWS_FL .AND.       &
                 NTHIK_FL )

    !.....  Read G_GRIDPATH file if all variables were found
    IF ( ALLFOUND .AND. GRIDPATH ) THEN

    !.....  Check length of character variables that matter
        L = LEN_TRIM( GNBUF )
        IF ( L .GT. NAMLEN3 ) THEN

            WRITE( MESG,94010 ) 'Grid name "' // GNBUF( 1:L ) //        &
                   '" has maximum allowable length of', NAMLEN3,        &
                   '.' // CRLF() // BLANK10 // 'Truncating to "' //     &
                   GNBUF( 1:NAMLEN3 ) // '".'
            CALL M3WARN( PROGNAME, 0, 0, MESG )

        END IF
        GNAME = GNBUF( 1:NAMLEN3 )

        !.....  Set grid description, and truncate
        GDESC = GDBUF

        !.....  Make sure everything important has been defined
        CALL GRDINFO_DEFINED( 'GDNAME_GD', GNAME, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'GDDESC_GD', GDESC, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'GDTYP_GD', CNAME, RDUM, IDUM, M3CHAR )
        CALL GRDINFO_DEFINED( 'P_ALP_GD', CDUM, P_ALP, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'P_BET_GD', CDUM, P_BET, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'P_GAM_GD', CDUM, P_GAM, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XCENT_GD', CDUM, XCENT, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YCENT_GD', CDUM, YCENT, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XORIG_GD', CDUM, XORIG, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YORIG_GD', CDUM, YORIG, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'XCELL_GD', CDUM, XCELL, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'YCELL_GD', CDUM, YCELL, IDUM, M3DBLE )
        CALL GRDINFO_DEFINED( 'NCOLS',    CDUM, RDUM, NCOLS, M3INT )
        CALL GRDINFO_DEFINED( 'NROWS',    CDUM, RDUM, NROWS, M3INT )
        CALL GRDINFO_DEFINED( 'NTHIK',    CDUM, RDUM, NTHIK, M3INT )

    ELSE

        !.....  Intialize grid description to NA for DSCGRID
        GDESC = 'No description available'          !  not used in DSCGRID

        !.....  set GNAME environment variable for DSCGRID
        MESG = 'grid system name for DSCGRID.F'
        CALL ENVSTR( 'IOAPI_GRIDNAME_1', MESG, ' ', GNAME, IOS2 )

        IF( IOS2 .NE. 0 ) THEN
            WRITE( MESG,94010 ) 'I/O error setting env vble "GNAME": IOSTAT=', IOS2
            CALL M3MESG( MESG )
            DSCM3GRD = .FALSE.
            RETURN
        END IF

        !.....  check GRIDDESC file format found

        IF ( DSCGRID( GNAME, CNAME, CTYPE, P_ALP, P_BET,        &
                      P_GAM, XCENT, YCENT, XORIG, YORIG,        &
                      XCELL, YCELL, NCOLS, NROWS, NTHIK ) ) THEN

            ! note: GRIDDESC file is still open from DSCGRID.F function
            GRIDPATH = .FALSE.
            EFLAG    = .FALSE.

        ELSE

            !.....  Both G_GRIDPATH and GRIDDESC file formats failed
            EFLAG = .TRUE.

        END IF      ! end GRIDDESC file format check

    END IF

    DSCM3GRD = ( .NOT. EFLAG )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats.... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

    !********************** INTERNAL SUBPROGRAMS ****************************

CONTAINS

    SUBROUTINE GRDINFO_CHECKER( STRINGS, OUTTYPE, ROUT, IOUT )

        !.....  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: STRINGS( 3 )
        INTEGER     , INTENT (IN) :: OUTTYPE
        REAL(8)     , INTENT(OUT) :: ROUT             ! real output value
        INTEGER     , INTENT(OUT) :: IOUT             ! integer output value

        !.....  Local variables
        INTEGER   I, L

        CHARACTER(300) MESG

        !..........

        IF( OUTTYPE .EQ. M3INT ) THEN

            IF( .NOT. CHKINT( STRINGS( 3 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // TRIM( STRINGS( 1 ) ) //     &
                       ' value is invalid in grid information file.'
                CALL M3MSG2( MESG )
            ELSE
                IOUT = STR2INT( STRINGS( 3 ) )
            END IF

        ELSE IF( OUTTYPE .EQ. M3DBLE ) THEN

            IF( .NOT. CHKREAL( STRINGS( 3 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // TRIM( STRINGS( 3 ) ) //     &
                       ' value is invalid in grid information file.'
                CALL M3MSG2( MESG )
            ELSE
                ROUT = DBLE( STR2REAL( STRINGS( 3 ) ) )
            END IF

        END IF

        RETURN

    END SUBROUTINE GRDINFO_CHECKER

    !..........
    !..........

    SUBROUTINE GRDINFO_DEFINED( KEYWORD, CVAL, RVAL, IVAL, INTYPE )

        !.....  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: KEYWORD
        CHARACTER(*), INTENT (IN) :: CVAL
        REAL(8)     , INTENT (IN) :: RVAL
        INTEGER     , INTENT (IN) :: IVAL
        INTEGER     , INTENT (IN) :: INTYPE

        !.....  Local variables
        CHARACTER(300) MESG

        !..........

        IF( ( INTYPE .EQ. M3CHAR .AND. CVAL .EQ. ' '    ) .OR.      &
            ( INTYPE .EQ. M3DBLE .AND. RVAL .LT. AMISS3 ) .OR.      &
            ( INTYPE .EQ. M3INT  .AND. IVAL .EQ. IMISS3 ) ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Keyword "' // KEYWORD //         &
                   '" not found in grid information file.'
            CALL M3MSG2( MESG )

        END IF

    END SUBROUTINE GRDINFO_DEFINED

END FUNCTION DSCM3GRD

