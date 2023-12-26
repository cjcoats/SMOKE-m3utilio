
LOGICAL FUNCTION DSCM3LAY( NLAYS, VGTYP, VGTOP, VGLVS )

!***********************************************************************
!  function body starts at line
!
!  DESCRIPTION:
!      This function returns true it successfully reads the properties of the
!      models-3 layer structure, and false otherwise.  It also returns those
!      properties using the subroutine arguments.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Copied and modified from DSCM3LAY, version 1.10
!
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!       Bug-fixes in GRDINFO_CHECKER, GRDINFO_DEFINED
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

    INCLUDE 'IOCNST3.EXT'

!...........   ARGUMENTS and their descriptions.  All are output variables:

    INTEGER, INTENT( OUT ) :: NLAYS      ! number of layers
    INTEGER, INTENT( OUT ) :: VGTYP      ! type of vertical coordinates
    REAL   , INTENT( OUT ) :: VGTOP      ! model-top, for sigma coord types
    REAL   , INTENT( OUT ) :: VGLVS( MXLAYS3 + 1 ) ! vertical coord values

!...........   EXTERNAL FUNCTIONS:
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL
    INTEGER, EXTERNAL :: GETNLIST

!...........   Local parameters
    INTEGER, PARAMETER :: MXLAYTYP = 6

!...........   Layer structure names arrays
    CHARACTER(7):: LAYNAMES( MXLAYTYP ) = ( / 'VGSGPH3'&
    &                                        , 'VGSGPN3'&
    &                                        , 'VGSIGZ3'&
    &                                        , 'VGPRES3'&
    &                                        , 'VGZVAL3'&
    &                                        , 'VGHVAL3' / )

!...........   Local arrays (note- tried to make these allocatable, but
!              this caused unexplainable crashing on SGI).
    CHARACTER(32) :: SEGMENT( 32 )
    CHARACTER(32) :: UPCSGMT( 32 )

!...........   File units and logical/physical names:
    INTEGER, SAVE :: IDEV    !  unit number of grid information file
    CHARACTER(16)    LNAME   !  logical name for grid information file

!...........   Scratch local variables and their descriptions:

    INTEGER         I, J, L, L2, LB, N  !  indices and string lengths
    INTEGER         IOS      !  I/O status return
    INTEGER         IREC     !  record counter
    INTEGER         NLAYREAD !  number of layers for reading sigma levels

    REAL(8)         DMISS3   !  double precision missing value

    LOGICAL      :: EFLAG = .FALSE.   ! true: error detected
    LOGICAL      :: FIRSTIME = .TRUE. ! true: first time routine called

    CHARACTER(300)  BUFFER   !  multi-purpose buffer
    CHARACTER(300)  LINE     !  line of file
    CHARACTER(300)  MESG     !  message buffer
    CHARACTER(300)  UPCLINE  !  upper case line of file

    CHARACTER(16) :: PROGNAME = 'DSCM3LAY' ! program name

!***********************************************************************
!   begin body of function DSCM3LAY

!.........  Get value of Models-3 environment variable for grid and check
!           exit status
    FIRSTIME = .FALSE.
    LNAME = 'G_GRIDPATH'

!.........  Open grid file
    IDEV = GETEFILE( LNAME, .TRUE., .TRUE., PROGNAME )

    IF( IDEV .LE. 0 ) THEN

        MESG = 'Could not open file "' // TRIM( LNAME ) // '"!'
        DSCM3LAY = .FALSE.
        RETURN

    END IF

!.........  Initialize grid setting to missing
    NLAYS   = IMISS3
    VGTYP   = IMISS3
    DMISS3  = DBLE( BADVAL3 )

!.........  Make sure read it at the start of the file
    REWIND( IDEV )

!.........  Read grid information file
    IREC = 0
    DO

!.............  Read whole line from the file
        READ( IDEV, 93000, END=111, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,&
            &       'reading grid information file at line', IREC

            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Skip blank lines
        IF( LINE .EQ. ' ' ) CYCLE

!.............  Adjust line to left and create upper case line
        LINE = ADJUSTL( LINE )
        UPCLINE = LINE
        CALL UPCASE( UPCLINE )
        L2 = LEN_TRIM( LINE )

!.............  Skip line if it starts with a number
        IF( UPCLINE( 1:1 ) .GE. '0' .AND.&
        &    UPCLINE( 1:1 ) .LE. '9'       ) CYCLE

!.............  Parse the lower and upper case lines
        N = GETNLIST( L2, UPCLINE )
        IF( N .GT. 32 ) THEN
            WRITE( MESG,94010 ) 'WARNING: Line', IREC,&
            &       'in G_GRIDPATH file skipped because ' //&
            &       CRLF() // BLANK10 // 'more than 32' //&
            &       'fields were found.'
            CALL M3MSG2( MESG )
            CYCLE
        END IF
! note: add allocatable array
        SEGMENT = ' '
        UPCSGMT = ' '
        CALL PARSLINE( LINE   , N, SEGMENT )
        CALL PARSLINE( UPCLINE, N, UPCSGMT )

!.............  Search for keywords.  If the keyword exists, extract the
!               value(s) and cycle to the next line...

!.............  Grid name
        SELECT CASE ( UPCSGMT( 1 ) )
          CASE( 'NLAYS' )

!.................  Check to ensure the layer number is an integer
            IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                &       'ERROR: NLAYS value is not an integer '//&
                &       'in grid description ' // CRLF() // BLANK10//&
                &       'file at line', IREC
                CALL M3MSG2( MESG )

!.................  For integer input...
            ELSE
                NLAYS = STR2INT( SEGMENT( 3 ) )

            END IF

!.............  Grid description
          CASE( 'VG_TYP_GD' )

!.................  Check to ensure the layer number is an integer, and if
!                   not, see if layer type name is being given
            IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN

                VGTYP = INDEX1( SEGMENT( 3 ), MXLAYTYP, LAYNAMES )

                IF ( VGTYP .LE. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )&
                    &   'ERROR: VG_TYP_GD value is not an integer '//&
                    &   'or does not match the '// CRLF()// BLANK10//&
                    &   'I/O API keywords in grid description ' //&
                    &   'file at line', IREC
                    CALL M3MSG2( MESG )
                END IF

!.................  For integer input...
            ELSE
                VGTYP = STR2INT( SEGMENT( 3 ) )

            END IF

!.............  Coordinate system type
          CASE( 'VGTOP_GD' )

!.................  Check to ensure the layer number is an floating point value
            IF ( .NOT. CHKREAL( SEGMENT( 3 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                &       'ERROR: VGTOP_GD value cannot be read '//&
                &       'from grid description'// CRLF()// BLANK10//&
                &       'file at line', IREC
                CALL M3MSG2( MESG )

!.................  For real input...
            ELSE
                VGTOP = STR2REAL( SEGMENT( 3 ) )
            END IF

!.............  Layer levels
          CASE( 'VGLVS_GD' )

!.................  Check to ensure the layer number is an integer
            IF ( .NOT. CHKINT( SEGMENT( 3 ) ) ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                &       'ERROR: VGLVS_GD layer number is not '//&
                &       'and integer in grid '// CRLF()// BLANK10//&
                &       'description file ' //&
                &       'at line', IREC
                CALL M3MSG2( MESG )

!.................  For integer input...
            ELSE

!.....................  Store the number of levels in the file
                NLAYREAD = STR2INT( SEGMENT( 3 ) )

!.....................  Check to ensure that the number of layers available
!                       is consistent with the number of layers needed in file
                IF ( NLAYREAD .GT. MXLAYS3 + 1 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )&
                    &   'ERROR: Number of layers (', NLAYREAD,&
                    &   ') exceeds maximum (', MXLAYS3+1,&
                    &   ') at line', IREC
                    CALL M3MSG2( MESG )
                    CYCLE
                END IF

!.....................  Check to ensure that the number of layers available
!                       is consistent with the number specified elsewhere in the
!                       file
                IF ( NLAYREAD-1 .NE. NLAYS ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )&
                    &   'ERROR: Number of layers (', NLAYREAD-1,&
                    &   ') differs from NLAYS in file (', NLAYS,&
                    &   ') at line', IREC
                    CALL M3MSG2( MESG )
                    CYCLE

                END IF

!.....................  Read the layer structure
                IF( LB .EQ. 0 )&
                &    READ( IDEV, *, ERR=999 )&
                &        ( VGLVS( I ), I=1, NLAYREAD )

            END IF

        END SELECT

    END DO

!.........  Exit from read loop
111 CONTINUE

    DSCM3LAY = ( .NOT. EFLAG )

!.........  Close grid file
    CLOSE( IDEV )

    RETURN

999 DSCM3LAY = .FALSE.
    MESG = 'ERROR: Could not read vertical layer structure ' //&
    &       'from grid description file.'
    CALL M3MSG2( MESG )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

!********************** INTERNAL SUBPROGRAMS ****************************

CONTAINS

    SUBROUTINE GRDINFO_CHECKER( STRINGS, OUTTYPE,&
    &                            ROUT, IOUT )

!.............  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: STRINGS( 3 )
        INTEGER     , INTENT (IN) :: OUTTYPE
        REAL(8)     , INTENT(OUT) :: ROUT     ! real output value
        INTEGER     , INTENT(OUT) :: IOUT     ! integer output value

!.............  Local variables
        INTEGER   I, L

        CHARACTER(300) MESG

!..........................................................................

        IF( OUTTYPE .EQ. M3INT ) THEN

            IF( .NOT. CHKINT( STRINGS( 3 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // TRIM( STRINGS( 1 ) ) //&
                &       ' value is invalid in grid information' //&
                &       ' file.'
                CALL M3MSG2( MESG )
            ELSE
                IOUT = STR2INT( STRINGS( 3 ) )
            END IF

        ELSE IF( OUTTYPE .EQ. M3DBLE ) THEN

            IF( .NOT. CHKREAL( STRINGS( 3 ) ) ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: ' // TRIM( STRINGS( 3 ) ) //&
                &       ' value is invalid in grid information' //&
                &       ' file.'
                CALL M3MSG2( MESG )
            ELSE
                ROUT = DBLE( STR2REAL( STRINGS( 3 ) ) )
            END IF

        END IF

        RETURN

    END SUBROUTINE GRDINFO_CHECKER

!..........................................................................
!..........................................................................

    SUBROUTINE GRDINFO_DEFINED( KEYWORD, CVAL, RVAL,&
    &                            IVAL, INTYPE )

!.............  Subroutine arguments
        CHARACTER(*), INTENT (IN) :: KEYWORD
        CHARACTER(*), INTENT (IN) :: CVAL
        REAL(8)     , INTENT (IN) :: RVAL
        INTEGER     , INTENT (IN) :: IVAL
        INTEGER     , INTENT (IN) :: INTYPE

!.............  Local variables
        CHARACTER(300) MESG

!..........................................................................

        IF( ( INTYPE .EQ. M3CHAR .AND. CVAL .EQ. ' '    ) .OR.&
        &    ( INTYPE .EQ. M3DBLE .AND. RVAL .LT. AMISS3 ) .OR.&
        &    ( INTYPE .EQ. M3INT  .AND. IVAL .EQ. IMISS3 ) ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: Keyword "' // KEYWORD //&
            &       '" not found in grid information file.'
            CALL M3MSG2( MESG )

        END IF

    END SUBROUTINE GRDINFO_DEFINED

END FUNCTION DSCM3LAY

