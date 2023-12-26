
SUBROUTINE CONVRTXY( NSRC, CTYPE, GDNAM, P_ALP, P_BET,&
                     P_GAM, XCENT, YCENT, XVALS, YVALS )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine converts coordinates from lat-lon to the grid
    !      defined by the subroutine arguments.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Added POLGRD3 as a supported coord sys type - E. Giroux CNRC 03/2004
    !
    !       Version 10/2016 by C. Coats:  USE M3UTILIO
    !
    !       Version 11/2023 by CJC:  Conversion to ".f90",  and related changes
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

    !.....  SUBROUTINE ARGUMENTS
    INTEGER,      INTENT    (IN) :: NSRC              !  actual source count
    INTEGER,      INTENT    (IN) :: CTYPE             !  coord sys type
    CHARACTER(*), INTENT    (IN) :: GDNAM             !  grid name
    REAL(8),      INTENT    (IN) :: P_ALP             !  first, second, third map
    REAL(8),      INTENT    (IN) :: P_BET             !  projection descriptive
    REAL(8),      INTENT    (IN) :: P_GAM             !  parameters
    REAL(8),      INTENT    (IN) :: XCENT             !  lon for coord-system X=0
    REAL(8),      INTENT    (IN) :: YCENT             !  lat for coord-system Y=0
    REAL(8),      INTENT(IN OUT) :: XVALS( NSRC )     !  x location (input lon)
    REAL(8),      INTENT(IN OUT) :: YVALS( NSRC )     !  y location (input lat)

    !.......   Other local variables
    INTEGER     UZONE                   ! UTM zone
    INTEGER     S

    REAL        ALP, BET, GAM, XC, YC
    REAL        XLOC, XX, YLOC, YY      ! tmp x and y coordinates

    LOGICAL     EFLAG                   ! true: error detected

    CHARACTER(256) MESG                 !  message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'CONVRTXY'     ! program name

    !***********************************************************************
    !   begin body of subroutine CONVRTXY

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

            CALL LL2UTM( XLOC, YLOC, UZONE, XX, YY )
            XVALS( S ) = XX
            YVALS( S ) = YY

        ENDDO

    ELSE IF ( CTYPE .EQ. LAMGRD3 ) THEN

        IF( .NOT. LAMBERT( GDNAM, ALP, BET, GAM, XC, YC ) ) THEN
            MESG = 'ERROR: Could not initialize Lambert grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DO S = 1, NSRC

            XLOC = XVALS( S )
            YLOC = YVALS( S )

            IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

            IF ( .NOT. LL2LAM( XLOC, YLOC, XX, YY ) ) THEN
                EFLAG = .TRUE.

            ELSE

                XVALS( S ) = XX
                YVALS( S ) = YY

            END IF

        ENDDO

        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem converting coordinates ' //  &
                   'from Lat-lon to Lambert'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        ENDIF

    ELSE IF ( CTYPE .EQ. POLGRD3 ) THEN

        IF( .NOT. POLSTE( GDNAM, ALP, BET, GAM, XC, YC ) ) THEN
            MESG = 'ERROR: Could not initialize Polar ' //      &
                   'Stereographic grid.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        DO S = 1, NSRC

            XLOC = XVALS( S )
            YLOC = YVALS( S )

            IF( XLOC .LT. AMISS3 .OR. YLOC .LT. AMISS3 ) CYCLE

            IF ( .NOT. LL2POL( XLOC, YLOC, XX, YY ) ) THEN
                EFLAG = .TRUE.

            ELSE

                XVALS( S ) = XX
                YVALS( S ) = YY

            END IF

        END DO

        IF( EFLAG ) THEN
            MESG = 'ERROR: Problem converting coordinates ' //      &
                   'from Lat-lon to Polar Stereographic'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2)

        END IF

    ELSE                !  error
        WRITE( MESG,94010 ) 'ERROR: Do not know how to convert ' // &
               'from Lat-lon to coordinate type number', CTYPE

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF              !  if coord type UTM, Lambert, or Polar Stereographic

    RETURN

    !.......   Internal buffering formats.... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE CONVRTXY
