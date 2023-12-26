
SUBROUTINE FMTCSRC( INSTRING, NCIN, OUTBUFF, LENOUT )

!***********************************************************************
!  subroutine body starts at line 78
!
!  DESCRIPTION:
!      This subroutine formats the INSTRING assuming that is a string of
!      source characteristics that has proper segment lengths.  It adds
!      labels and newlines as necessary and returns the formatted OUTBUFF
!      and its length.
!
!  PRECONDITIONS REQUIRED:
!      INSTRING with correct segement lengths
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Functions: I/O API functions
!
!  REVISION  HISTORY:
!       Created 11/98 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!...........   Modules for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NCHARS, SC_BEGP, SC_ENDP

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: INSTRING ! Input string
    INTEGER     , INTENT (IN) :: NCIN     ! No. of chars to format (1 to 8)
    CHARACTER(*), INTENT(OUT) :: OUTBUFF  ! Formatted output string
    INTEGER     , INTENT(OUT) :: LENOUT   ! Length of output string

!...........   Local paramaters
    INTEGER, PARAMETER :: MXPTLBL = 8
    INTEGER, PARAMETER :: MXARLBL = 3
    INTEGER, PARAMETER :: MXMBLBL = 6

!.........  Output labels (note: these could be dynamic if in MODINFO)
    CHARACTER(6) :: PTLABEL( MXPTLBL ) =  &!  message buffer
    &              ( / 'Region', 'Plant ', 'Char1 ', 'Char2 ',&
    &                  'Char3 ', 'Char4 ', 'Char5 ', 'Data  ' / )

    CHARACTER(6) :: ARLABEL( MXARLBL ) =  &!  message buffer
    &              ( / 'Region', 'SCC   ', 'Data  ' / )

    CHARACTER(9) :: MBLABEL( MXMBLBL ) =  &!  message buffer
    &              ( / 'Region   ', 'Road Type', 'Link     ',&
    &                  'Vtype    ', 'SCC      ', 'Data     ' / )

    CHARACTER(9), SAVE :: LABEL( MXPTLBL )

!...........   Other local variables
    INTEGER         I, L, K, L1, L2, L3, L4  ! counters and indices

    INTEGER, SAVE :: TMPNUM      ! temporary number of chars
    INTEGER          NLOOP       ! number of loop iterations

    LOGICAL, SAVE :: FIRSTIME = .TRUE.

    CHARACTER(300)  BUFFER  !  string buffer

    CHARACTER(16) :: PROGNAME = 'FMTCSRC'  ! program name

!***********************************************************************
!   begin body of subroutine FMTCSRC

    IF( FIRSTIME ) THEN

        FIRSTIME = .FALSE.

        SELECT CASE( CATEGORY )
          CASE( 'AREA' )
            TMPNUM = MXARLBL
            LABEL( 1:MXARLBL ) = ARLABEL  ! array

          CASE( 'MOBILE' )
            TMPNUM = MXMBLBL
            LABEL( 1:MXMBLBL ) = MBLABEL  ! array

          CASE( 'POINT' )
            TMPNUM = MXPTLBL
            LABEL( 1:MXPTLBL ) = PTLABEL  ! array

          CASE DEFAULT
            TMPNUM = NCIN

        END SELECT

    END IF

!.........  Make sure not to exceed NCHARS legal value
    NLOOP = MIN( NCIN, TMPNUM )

!.........  Initialize output buffer
    OUTBUFF = ' '
    IF( NCHARS .GE. 1 ) THEN
        L1 = SC_BEGP( 1 )
        L2 = SC_ENDP( 1 )
        L4 = LEN_TRIM( LABEL( 1 ) )

        WRITE( OUTBUFF, 94900 ) LABEL( 1 )(1:L4), INSTRING( L1:L2 )

    ENDIF

!.........  Loop through the remaining source chars, writing the output string
!           for each populated source characterstic
    K = LEN_TRIM( OUTBUFF )
    DO I = 2, NLOOP

        L = LEN_TRIM( OUTBUFF )

        L1 = SC_BEGP( I )    ! retrieve stored field lengths
        L2 = SC_ENDP( I )

        BUFFER = ADJUSTL( INSTRING( L1:L2 ) )   ! Left-justifying
        L3 = LEN_TRIM( BUFFER )
        K = K + 8 + L3

!.............  Continue to contribute to buffer if not blank
        IF( BUFFER .NE. ' ' ) THEN

            L4 = LEN_TRIM( LABEL( I ) )

!.................  Include a return and indent if line gets too long
            IF( K .GT. EMOUTLN3 ) THEN
                K = 18 + L3
                WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) //&
                &     CRLF() // BLANK10 // LABEL( I )( 1:L4 ) ,&
                &       BUFFER( 1:L3 )
            ELSE
                WRITE( OUTBUFF, 94900 ) OUTBUFF( 1:L ) // ' ' //&
                &       LABEL( I )( 1:L4 ) , BUFFER( 1:L3 )

            ENDIF
        ENDIF

    ENDDO

    LENOUT = LEN_TRIM( OUTBUFF )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94900 FORMAT( A, ': ', A )

END SUBROUTINE FMTCSRC





