
SUBROUTINE RDB4FAC( NSEF, FDEV, NLINES, VGID, LINDX,&
&                    LFAC, WNTF, LWGT, FACS  )

!***********************************************************************
!  subroutine body starts at line XX
!
!  DESCRIPTION:
!       Reads in the BEIS3 emissions factors from the BFAC file.
!
!  PRECONDITIONS REQUIRED:
!
!  REVISION  HISTORY:
!       03/01 protoype by J. Vukovich
!       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
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
!***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!...........   ARGUMENTS and their descriptions: actually-occurring ASC table


    INTEGER, INTENT (IN)  :: NSEF    !  no. biogenic emission factors
    INTEGER, INTENT (IN)  :: FDEV    !  unit number for elev srcs file
    INTEGER, INTENT (IN)  :: NLINES  !  no. veg types

    CHARACTER(16), INTENT (OUT)  :: VGID( NLINES )       ! veg ids
    INTEGER, INTENT (OUT)        :: LINDX( NLINES )      ! leaf area index
    REAL, INTENT (OUT)           :: LFAC( NLINES )       ! leaf biomass
    REAL, INTENT (OUT)           :: WNTF( NLINES )       ! winter factor
    REAL, INTENT (OUT)           :: LWGT( NLINES )       ! specific leaf wgt
    REAL, INTENT (OUT)           :: FACS( NLINES, NSEF ) ! emis facs

    LOGICAL      :: EFLAG = .FALSE.  !  error flag
    INTEGER      :: MXSEG            ! # of potential line segments

    INTEGER       I, J               !  counters
    INTEGER       ISTAT              !  iostat error

    CHARACTER(50)   SEGMENT( NSEF + 7 )   ! Segments of parsed lines
    CHARACTER(300)  MESG             !  message buffer
    CHARACTER(300)  LINE             !  buffer for variables

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDB4FAC' ! program name

!***********************************************************************
!   begin body of subroutine RDB3FAC

!.........  Set number of potential line segments
    MXSEG = NSEF + 7

!.......... Read in emissions factors for each veg id

    DO I = 1, NLINES



        READ( FDEV, 93030, IOSTAT=ISTAT )  LINE
        IF ( I .LE. 2) CYCLE

        IF ( ISTAT .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'Error', ISTAT,&
            &    'reading EMISSION FACTOR file at line', I
            CALL M3MESG( MESG )
        END IF

!.............  Separate the line of data into each part
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

        VGID ( I ) = SEGMENT( 2 )
        LINDX( I ) = STR2INT ( SEGMENT ( 3 ) )
        LFAC ( I ) = STR2REAL( SEGMENT ( 4 ) )
        WNTF ( I ) = STR2REAL( SEGMENT ( 5 ) )
        LWGT ( I ) = STR2REAL( SEGMENT ( 6 ) )

!          write (*,*) I, VGID(I), LINDX(I), LFAC(I), WNTF(I), LWGT(I)

        DO J = 1, NSEF

            FACS( I , J ) = STR2REAL( SEGMENT ( J + 6 )  )

        ENDDO

    ENDDO

    IF( EFLAG ) THEN
        MESG = 'Problem reading biogenic emissions factors file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93010 FORMAT( 8X, A16, A )
93020 FORMAT( A16, A )
93030 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I10, :, 2X ) )

END SUBROUTINE RDB4FAC
