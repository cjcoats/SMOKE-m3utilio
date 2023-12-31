
        SUBROUTINE RDB3FAC( B309FLAG, NSEF, FDEV, NLINES, VGID, LINDX,
     &                      LFAC, WNTF, LWGT, FACS  )

C***********************************************************************
C  subroutine body starts at line XX
C
C  DESCRIPTION:
C       Reads in the BEIS3 emissions factors from the BFAC file.
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C       03/2001 protoype by J. Vukovich
C       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
C
C***********************************************************************
C
C Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
C                System
C File: @(#)$Id$
C
C COPYRIGHT (C) 2004, Environmental Modeling for Policy Development
C All Rights Reserved
C
C Carolina Environmental Program
C University of North Carolina at Chapel Hill
C 137 E. Franklin St., CB# 6116
C Chapel Hill, NC 27599-6116
C
C smoke@unc.edu
C
C Pathname: $Source$
C Last updated: $Date$
C
C***********************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   ARGUMENTS and their descriptions: actually-occurring ASC table
        LOGICAL, INTENT (IN)  :: B309FLAG  ! true: using v3.09 of BEIS

        INTEGER, INTENT (IN)  :: NSEF    !  no. biogenic emission factors
        INTEGER, INTENT (IN)  :: FDEV    !  unit number for elev srcs file
        INTEGER, INTENT (IN)  :: NLINES  !  no. veg types

        CHARACTER(16), INTENT (OUT)  ::  VGID( NLINES )       ! veg ids
        INTEGER, INTENT (OUT)        :: LINDX( NLINES )      ! leaf area index
        REAL,    INTENT (OUT)        ::  LFAC( NLINES )       ! leaf biomass
        REAL,    INTENT (OUT)        ::  WNTF( NLINES )       ! winter factor
        REAL,    INTENT (OUT)        ::  LWGT( NLINES )       ! specific leaf wgt
        REAL,    INTENT (OUT)        ::  FACS( NLINES, NSEF ) ! emis facs

        LOGICAL      :: EFLAG            !  error flag
        INTEGER      :: MXSEG            ! # of potential line segments

        INTEGER       I, J               !  counters
        INTEGER       ISTAT              !  iostat error

        CHARACTER(20)   SEGMENT( NSEF + 4 )   ! Segments of parsed lines
        CHARACTER(300)  MESG             !  message buffer
        CHARACTER(300)  LINE             !  buffer for variables

        CHARACTER(16), PARAMETER :: PROGNAME = 'RDB3FAC' ! program name

C***********************************************************************
C   begin body of subroutine RDB3FAC

C.........  Set number of potential line segments
        MXSEG = NSEF + 4
        EFLAG = .FALSE.

C.......... Read in emissions factors for each veg id

        DO I = 1, NLINES

C...........  Factors files have slightly different formats
          IF( B309FLAG ) THEN
             READ( FDEV, 93010, IOSTAT=ISTAT ) VGID( I ) , LINE
          ELSE
             READ( FDEV, 93020, IOSTAT=ISTAT ) VGID( I ) , LINE
          END IF

          IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', ISTAT,
     &              'reading EMISSION FACTOR file at line', I
               CALL M3MESG( MESG )
          END IF

C.............  Separate the line of data into each part
          CALL PARSLINE( LINE, MXSEG, SEGMENT )

          LINDX( I ) = STR2INT ( SEGMENT ( 1 ) )
          LFAC ( I ) = STR2REAL( SEGMENT ( 2 ) )
          WNTF ( I ) = STR2REAL( SEGMENT ( 3 ) )
          LWGT ( I ) = STR2REAL( SEGMENT ( 4 ) )

          DO J = 1, NSEF

            FACS( I , J ) = STR2REAL( SEGMENT ( J + 4 )  )

          ENDDO

        ENDDO

        IF( EFLAG ) THEN
            MESG = 'Problem reading biogenic emissions factors file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx

93010   FORMAT( 8X, A16, A )
93020   FORMAT( A16, A )

C...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10 ( A, :, I10, :, 2X ) )

        END SUBROUTINE RDB3FAC
