
SUBROUTINE CONVRTLL( NSRC, CTYPE, GDNAM, P_ALP, P_BET,  &
                     P_GAM, XCENT, YCENT, XVALS, YVALS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine converts coordinates from the grid
    !      defined by the subroutine arguments to lat-lon.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created based on CONVRTXY 07/13 by C. Seppanen
    !
    !       Version 10/2016 by C. Coats:  USE M3UTILIO
    !       Version 11/2023 by CJC:  conversion to ".f90" source format
    !****************************************************************************/
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

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....  SUBROUTINE ARGUMENTS
    INTEGER,            INTENT (IN) :: NSRC              !  actual source count
    INTEGER,            INTENT (IN) :: CTYPE             !  coord sys type
    CHARACTER(NAMLEN3), INTENT (IN) :: GDNAM             !  grid name
    REAL(8),            INTENT (IN) :: P_ALP             !  first, second, third map
    REAL(8),            INTENT (IN) :: P_BET             !  projection descriptive
    REAL(8),            INTENT (IN) :: P_GAM             !  parameters
    REAL(8),            INTENT (IN) :: XCENT             !  lon for coord-system X=0
    REAL(8),            INTENT (IN) :: YCENT             !  lat for coord-system Y=0
    REAL(8),            INTENT(IN OUT) :: XVALS( NSRC )     !  x location (input grid coord)
    REAL(8),            INTENT(IN OUT) :: YVALS( NSRC )     !  y location (input grid coord)

    !.......   Other local variables
    INTEGER     UZONE                   ! UTM zone
    INTEGER     S
    REAL        ALP, BET, GAM, XC, YC
    REAL        XLOC, XX, YLOC, YY      ! tmp x and y coordinates

    LOGICAL     EFLAG                   ! true: error detected

    CHARACTER(NAMLEN3) TMPGDNAM         ! temporary grid name
    CHARACTER(256)     MESG             ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'CONVRTLL'     ! program name

    !***********************************************************************
    !   begin body of subroutine CONVRTLL

    !.....  Copy input grid name to temporary variable since I/O API
    !           routines may change it to the coordinate system name

    TMPGDNAM = GDNAM
    EFLAG    =.FALSE.
    ALP      = P_ALP
    BET      = P_BET
    GAM      = P_GAM
    XC       = XCENT
    YC       = YCENT

    IF ( CTYPE .EQ. LATGRD3 ) THEN
        RETURN

    ELSE IF ( CTYPE .EQ. UTMGRD3 ) THEN

        UZONE = NINT( P_ALP )

        DO S = 1, NSRC

            XLOC = XVALS( S )
            YLOC = YVALS( S )

            IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

            CALL UTM2LL( XLOC, YLOC, UZONE, XX, YY )
            XVALS( S ) = XX
            YVALS( S ) = YY

        ENDDO

    ELSE IF ( CTYPE .EQ. LAMGRD3 ) THEN

        IF( .NOT.LAMBERT( TMPGDNAM, ALP, BET, GAM, XC, YC ) ) THEN
            MESG = 'ERROR: Could not initialize Lambert grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DO S = 1, NSRC

            XLOC = XVALS( S )
            YLOC = YVALS( S )

            IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

            IF ( .NOT. LAM2LL( XLOC, YLOC, XX, YY ) ) THEN
                EFLAG = .TRUE.

            ELSE

                XVALS( S ) = XX
                YVALS( S ) = YY

            END IF

        ENDDO

        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem converting coordinates from Lambert to Lat-lon'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

    ELSE IF ( CTYPE .EQ. POLGRD3 ) THEN

        IF( .NOT. POLSTE( TMPGDNAM, ALP, BET, GAM, XC, YC ) ) THEN
            MESG = 'ERROR: Could not initialize Polar Stereographic grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DO S = 1, NSRC

            XLOC = XVALS( S )
            YLOC = YVALS( S )

            IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

            IF ( .NOT. POL2LL( XLOC, YLOC, XX, YY ) ) THEN
                EFLAG = .TRUE.

            ELSE

                XVALS( S ) = XX
                YVALS( S ) = YY

            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem converting coordinates from Polar Stereographic to Lat-lon'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2)

        END IF

    ELSE                !  error
        WRITE( MESG,94010 ) 'ERROR: Do not know how to convert from coordinate type number', CTYPE, 'to Lat-lon'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF              !  if coord type UTM, Lambert, or Polar Stereographic

    RETURN


    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats.... 93xxx

93030 FORMAT( I8.8 )

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CONVRTLL
