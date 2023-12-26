
SUBROUTINE RDSRGHDR( VFLAG, FDEV, SRGFMT )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine reads the header of the spatial surrogates file.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
    !       Version 11/2023 by CJC:  USE M3UTILIO".f90" free source format, and
    !       related changes
    !**************************************************************************
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
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !...........   Subroutine arguments
    LOGICAL      , INTENT  (IN) :: VFLAG          ! true: using variable grid
    INTEGER      , INTENT  (IN) :: FDEV           ! File unit number
    CHARACTER(*) , INTENT  (OUT):: SRGFMT         ! Format of surrogates file

    !...........   EXTERNAL FUNCTIONS and their descriptions:

    LOGICAL, EXTERNAL :: DSCM3GRD
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: CHKREAL

    !...........   Local parameters

    INTEGER, PARAMETER :: MXSEG = 16              ! # of potential line segments
    INTEGER, PARAMETER :: MXGRDTYP = 13

    !...........   Grid types and names arrays
    INTEGER, PARAMETER:: GRDTYPES( MXGRDTYP ) = (/  LATGRD3     &
                                                  , LATGRD3     &
                                                  , LATGRD3     &
                                                  , LAMGRD3     &
                                                  , LAMGRD3     &
                                                  , MERGRD3     &
                                                  , MERGRD3     &
                                                  , STEGRD3     &
                                                  , STEGRD3     &
                                                  , POLGRD3     &
                                                  , POLGRD3     &
                                                  , UTMGRD3     &
                                                  , UTMGRD3 /)

    CHARACTER(16), PARAMETER :: GRDNAMES( MXGRDTYP ) = (/  'LAT-LON        '    &
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

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDSRGHDR'         ! Program name

    !...........   Other arrays

    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines

    !...........   Local variables

    INTEGER       I, J, L        ! counters and indices
    INTEGER       IOS            ! i/o status
    INTEGER       IREC           ! record counter
    INTEGER       NTHIK          ! boundary thickness

    LOGICAL    :: EFLAG = .FALSE.               ! Error flag

    CHARACTER(7)       PROJUNIT       ! Projection units
    CHARACTER(80)      MSGEND         ! tmp end of message
    CHARACTER(320)     LINE           ! Read buffer for a line
    CHARACTER(300)     MESG           ! Message buffer
    CHARACTER(NAMLEN3) COORD3D        ! coordinate system name
    CHARACTER(NAMLEN3) COORUN3D       ! coordinate system projection units
    CHARACTER(NAMLEN3) PROJTYPE       ! coordinate system projection name

    !***********************************************************************
    !   Begin body of subroutine RDSRGHDR

    IREC  = 0
    EFLAG = .FALSE.

    REWIND( FDEV )

    !.........  Loop through lines of file until the header line is encountered
    DO

        READ ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010)                                         &
              'I/O error', IOS, 'reading gridding information '//       &
              'file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        CALL UPCASE( LINE )

        IF ( LINE .EQ. ' ' ) CYCLE          ! skip all blank lines

        !.............  Determine if current line is the header
        IF( VFLAG ) THEN
            I = INDEX( LINE, '#VARIABLE_GRID' )
        ELSE
            I = INDEX( LINE, '#GRID' )
        END IF

        IF ( I .GT. 0 ) THEN
            SRGFMT = 'MODELS3'      ! set format to 'MODELS3'
            IF( IREC .NE. 1 ) THEN
                MESG = 'First line of the file did not contain the header line.'
                CALL M3WARN( PROGNAME, 0, 0, MESG )
            END IF

            EXIT                    ! exit read loop with header found

        END IF

    END DO

    SELECT CASE( SRGFMT )

      CASE ( 'MODELS3' )

        !.............  Parse the line of data into segments based on the rules
        !               for "list-formatted" in fortran, but not requiring
        !               quotes around the text strings

        CALL PARSLINE( LINE, MXSEG, SEGMENT )

        !.............  Make sure appropriate segments of header line appear to
        !               be either INTEGER or REAL as expected

        IF ( .NOT. CHKREAL( SEGMENT( 3 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 4 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 5 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 6 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKINT ( SEGMENT( 7 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKINT ( SEGMENT( 8 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKINT ( SEGMENT( 9 ) ) )  EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 12 ) ) ) EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 13 ) ) ) EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 14 ) ) ) EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 15 ) ) ) EFLAG = .TRUE.
        IF ( .NOT. CHKREAL( SEGMENT( 16 ) ) ) EFLAG = .TRUE.

        IF ( EFLAG ) THEN

            WRITE( MESG, 94010 )        &
              'Unexpected data type encountered in header of the file'
            CALL M3MESG( MESG )

        ELSE

            GDNAM3D =            SEGMENT ( 2 )
            XORIG3D  = STR2DBLE( SEGMENT ( 3 ) )
            YORIG3D  = STR2DBLE( SEGMENT ( 4 ) )
            XCELL3D  = STR2DBLE( SEGMENT ( 5 ) )
            YCELL3D  = STR2DBLE( SEGMENT ( 6 ) )
            NCOLS3D  = STR2INT ( SEGMENT ( 7 ) )
            NROWS3D  = STR2INT ( SEGMENT ( 8 ) )
            NTHIK3D  = STR2INT ( SEGMENT ( 9 ) )
            PROJTYPE =           SEGMENT ( 10 )
            PROJUNIT =           SEGMENT ( 11 )
            P_ALP3D  = STR2DBLE( SEGMENT ( 12 ) )
            P_BET3D  = STR2DBLE( SEGMENT ( 13 ) )
            P_GAM3D  = STR2DBLE( SEGMENT ( 14 ) )
            XCENT3D  = STR2DBLE( SEGMENT ( 15 ) )
            YCENT3D  = STR2DBLE( SEGMENT ( 16 ) )

        END IF

    END SELECT

    !.........  Set project code based on projection type
    J =  INDEX1( PROJTYPE, MXGRDTYP, GRDNAMES )
    IF ( J .LE. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'Projection type "' // TRIM( PROJTYPE ) // '" is not recognized.'
        CALL M3MSG2( MESG )

    !.........  Otherwise, set the grid type code number
    !.........  Initialize grid information based on the surrogates file
    ELSE

        GDTYP3D = GRDTYPES( J )
        CALL CHKGRID( GDNAM3D, 'SURROGATES', 0, EFLAG )

    END IF

    !.........  Abort if an error is found
    IF ( EFLAG ) THEN

        MESG = 'Problem with processing grid information'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    RETURN

    !.........  If get here, then the header was never found    !
999 MESG = 'Surrogates file is missing a header line'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDSRGHDR





