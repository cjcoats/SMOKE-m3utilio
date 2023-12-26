
SUBROUTINE RDTZONE( FDEV, NDIM, NZS, NZF, TZONE0,&
&                    TZONST, TFIPST, TZONEF, TFIPEF )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine reads the time zones file, sorts it, and returns
!      the sorted data
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    INTEGER      FDEV              !  time zones file unit number
    INTEGER      NDIM              !  dimension for time zone arrays
    INTEGER      NZS               !  no of state-specific
    INTEGER      NZF               !  no of county-specific
    INTEGER      TZONE0            !  fallback zone
    INTEGER      TZONST( NDIM )    !  state-specific time zones
    INTEGER      TFIPST( NDIM )    !  state FIPS codes (2 digit)
    INTEGER      TZONEF( NDIM )    !  FIPS-specific time zones
    INTEGER      TFIPEF( NDIM )    !  state/county FIPS codes (5 digits)

!...........   Unsorted time zone records
    INTEGER      INDXSA( NDIM )
    INTEGER      TZONSA( NDIM )
    INTEGER      TFIPSA( NDIM )
    INTEGER      INDXFA( NDIM )
    INTEGER      TZONFA( NDIM )
    INTEGER      TFIPFA( NDIM )

!...........   Other local variables
    INTEGER         FIP, TZONE          ! tmp fips code and time zone
    INTEGER         I, J                ! counters and indices
    INTEGER         IOS                 ! i/o status
    INTEGER         IREC                ! record counter

    LOGICAL      :: EFLAG = .FALSE.     ! error flag

    CHARACTER(300)  LINE    !  Input line from POINT file
    CHARACTER(300)  MESG    !  message buffer

    CHARACTER(16) :: PROGNAME = 'RDTZONE' ! program name

!***********************************************************************
!   begin body of subroutine RDTZONE

    TZONE0 = 5      !  default:  EST
    NZS    = 0
    NZF    = 0
    IREC   = 0
    DO

        READ( FDEV, *, END=12, IOSTAT=IOS ) FIP, TZONE
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )&
            &       'I/O error', IOS,  'reading time zones file ' //&
            &       'at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        IF ( FIP .EQ. 0 ) THEN              !  fallback -- all sources

            TZONE0 = TZONE

        ELSE IF ( MOD( FIP,100 ) .EQ. 0 ) THEN     !  state-specific zone

            NZS = NZS + 1
            IF ( NZS .LE. NDIM ) THEN
                INDXSA( NZS ) = NZS
                TFIPSA( NZS ) = FIP / 1000
                TZONSA( NZS ) = TZONE
            END IF

        ELSE                                        !  county-specific zone

            NZF = NZF + 1
            IF ( NZF .LE. NDIM ) THEN
                INDXFA( NZF ) = NZF
                TFIPFA( NZF ) = FIP
                TZONFA( NZF ) = TZONE
            END IF

        END IF      !  if fip zero, or nn000, or not.

    ENDDO

12  CONTINUE        !  exit from loop reading ZDEV

    IF ( NZS .GT. NDIM .OR. NZF .GT. NDIM ) THEN
        WRITE( MESG,94010 )&
        &  'Number of state-only records  :', NZS, CRLF()// BLANK5//&
        &  'Number of state&county records:', NZF, CRLF()// BLANK5//&
        &  'Memory allocated:', NDIM
        CALL M3MSG2( MESG )

        MESG = 'INTERNAL ERROR: Insufficient memory allocated ' //&
        &       'for time zones tables'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
    END IF

    CALL SORTI1( NZS, INDXSA, TFIPSA )
    DO I = 1, NZS
        J = INDXSA( I )
        TZONST( I ) = TZONSA( J )
        TFIPST( I ) = TFIPSA( J )
    ENDDO

    CALL SORTI1( NZF, INDXFA, TFIPFA )
    DO I = 1, NZF
        J = INDXFA( I )
        TZONEF( I ) = TZONFA( J )
        TFIPEF( I ) = TFIPFA( J )
    ENDDO

!.........  Rewind file

    REWIND( FDEV )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************
!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94120 FORMAT( I5.5 )

END  SUBROUTINE RDTZONE
