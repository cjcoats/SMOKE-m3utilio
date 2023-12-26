
SUBROUTINE BLDENAMS( CATEGORY_L, NIPPA_L, NPPOA, EANAM_L,&
&                     OUTNAMES, OUTUNITS, OUTTYPES, OUTDESCS )

!***********************************************************************
!  subroutine body starts at line 154
!
!  DESCRIPTION:
!      This subroutine builds names for pol/act-specific inventory variables
!      such as rule effectiveness and control efficiency, while ensuring that
!      no names are duplicated.  Since inventory pol/act names are allowed
!      to be 16 characters, there needs to be a truncation and check for
!      duplicates, and then an insertion of the prefixes.  If the UNITS
!      header option has been used in an inventory input file, it assigns
!      variable units based on these settings.
!
!  PRECONDITIONS REQUIRED:
!      Source CATEGORY_L specified correctly
!      Pollutant/activity list and number provided
!      Memory allocated for output arrays
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!     Subroutines: Models-3 subroutines
!     Functions: Models-3 functions
!
!  REVISION  HISTORY:
!       Created 12/98 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!*************************************************************************
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

!.........  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: INVDNAM, INVDDSC, MXIDAT

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: CATEGORY_L               ! source CATEGORY
    INTEGER     , INTENT (IN) :: NIPPA_L                  ! no. of pol/act
    INTEGER     , INTENT (IN) :: NPPOA                    ! no. vars per
    ! pol or activity
    CHARACTER(*), INTENT (IN) :: EANAM_L ( NIPPA_L )      ! pol/act names
    CHARACTER(*), INTENT(OUT) :: OUTNAMES( NIPPA_L,NPPOA )! var names
    CHARACTER(*), INTENT(OUT) :: OUTUNITS( NIPPA_L,NPPOA )! var units
    INTEGER     , INTENT(OUT) :: OUTTYPES( NIPPA_L,NPPOA )! var typ:int/real
    CHARACTER(*), INTENT(OUT) :: OUTDESCS( NIPPA_L,NPPOA )! var descriptions

!.........  Data arrays for variable names...

!.........  Area source variable name parameters
    CHARACTER(CPRTLEN3), PARAMETER :: ARPREFIX( 2:NARPPOL3 ) =&
    &     (/ AVEDAYRT,   EMISFCRT,   CTLEFFRT, RULEFFRT, RULPENRT /)

    CHARACTER(NAMLEN3), PARAMETER :: ARUNITS( NARPPOL3 ) =&
    &     (/ 'tons/yr  ', 'tons/day ', 'SCC units',&
    &        '%        ', '%        ', '%        ' /)

    INTEGER, PARAMETER :: ARTYPES( NARPPOL3 ) =&
    &     (/ M3REAL, M3REAL, M3REAL, M3REAL, M3REAL, M3REAL /)

    CHARACTER(MXDLEN3), PARAMETER :: ARDESCS( NARPPOL3 )=&
    &     (/&
    &     'Annual Emissions                                      ',&
    &     'Average Day Emissions                                 ',&
    &     'Emission Factors                                      ',&
    &     'Control efficiency (in [0,1], or "MISSING": < -9.0E36)',&
    &     'Rule effectiveness (in [0,1], or "MISSING": < -9.0E36)',&
    &     'Rule penetration (in [0,1],   or "MISSING": < -9.0E36)'&
    &     /)

!.........  Mobile source variable name parameters
    CHARACTER(CPRTLEN3), PARAMETER :: MBPREFIX( 2:NMBPPOL3 ) = (/ AVEDAYRT /)

    CHARACTER(NAMLEN3), PARAMETER :: MBUNITS( NMBPPOL3 ) = (/&
    &     'tons/yr ', 'tons/day' /)

    INTEGER MBTYPES( NMBPPOL3 )
    DATA MBTYPES / M3REAL, M3REAL /

    CHARACTER(MXDLEN3), PARAMETER :: MBDESCS( NMBPPOL3 ) =&
    &     (/ 'Annual Data     ',&
    &        'Average Day Data' /)

!.........  Point source variable name parameters
    CHARACTER(CPRTLEN3) PTPREFIX( 2:NPTPPOL3 )
    DATA PTPREFIX / AVEDAYRT, CTLEFFRT, RULEFFRT,&
    &                EMISFCRT, CECOD1RT, CECOD2RT /

    CHARACTER(NAMLEN3), PARAMETER :: PTUNITS( NPTPPOL3 )=&
    &            (/ 'tons/yr  ', 'tons/day ',&
    &               '%        ', '%        ',&
    &               'SCC units', 'n/a      ',&
    &               'n/a      '                    /)

    INTEGER, PARAMETER :: PTTYPES( NPTPPOL3 ) =&
    &    (/ M3REAL, M3REAL, M3REAL, M3REAL, M3REAL, M3INT, M3INT /)

    CHARACTER(MXDLEN3), PARAMETER :: PTDESCS( NPTPPOL3 ) =&
    &     (/&
    &     'Annual Emissions                                      ',&
    &     'Average Day Emissions                                 ',&
    &     'Control efficiency (in [0,1], or "MISSING": < -9.0E36)',&
    &     'Rule effectiveness (in [0,1], or "MISSING": < -9.0E36)',&
    &     'Emission Factors                                      ',&
    &     'Primary Control Equipment Code                        ',&
    &     'Secondary Control Equipment Code                      '&
    &     /)

!...........   Unsorted pollutant/activity records
    INTEGER       INDEXA ( NIPPA_L )
    CHARACTER(NAMLEN3) ABRNAMA( NIPPA_L )    !  pollutant/activity names

!...........   Other local variables
    INTEGER         COD       !  tmp for pollutant/activity code
    INTEGER         I, J, K   !  counters and indices
    INTEGER         IDIF      !  abridged name length
    INTEGER         L, L2, LU !  length indices
    INTEGER         LCNT      !  same-abrigded-name counter
    INTEGER         POLLNUM   !  position of current pollutant in INVDNAM

    CHARACTER(NAMLEN3)  LNAM  !  previous pollutant/activity name
    CHARACTER(NAMLEN3)  NAM   !  current pollutant/activity name
    CHARACTER(10)   BUFFER    !  tmp string buffer
    CHARACTER(300)  MESG      !  message buffer

    CHARACTER(16) :: PROGNAME = 'BLDENAMS' ! program name

!***********************************************************************
!   begin body of subroutine BLDENAMS

!.........  Truncate variable names based on original length and length of
!.........  fields to be inserted

    IDIF = NAMLEN3 - CPRTLEN3
    DO I = 1, NIPPA_L
        INDEXA ( I ) = I
        L = MIN( IDIF, LEN_TRIM( EANAM_L( I ) ) )
        ABRNAMA( I ) = EANAM_L( I )( 1:L )
    END DO

!.........  Sort to set up for duplicates search
    CALL SORTIC( NIPPA_L, INDEXA, ABRNAMA )

!.........  Search for duplicates and rename them
    LNAM = EMCMISS3
    LCNT = 0
    DO I = 1, NIPPA_L

        J = INDEXA( I )
        NAM = ABRNAMA( J )

        IF( NAM .EQ. LNAM ) THEN
            LCNT = LCNT + 1
            WRITE( BUFFER, '(I10)' ) LCNT  ! turn in to string
            BUFFER = ADJUSTL( BUFFER )     ! left justify

            L  = LEN_TRIM( BUFFER )         ! get length of trailer
            L2 = MIN( IDIF -L -1, LEN_TRIM( NAM ) )  ! remaining name
            ABRNAMA( J ) = NAM( 1:L2 ) // '_' // BUFFER( 1:L ) ! store name
        ELSE
            LCNT = 0
            LNAM = NAM
        ENDIF

    ENDDO

!.........  Build output names
    DO I = 1, NIPPA_L

!.............  If pol/act name is in inventory table, append its description
!               to the annual data description
        NAM = EANAM_L( I )
        IF( ALLOCATED( INVDNAM ) ) THEN
            POLLNUM = INDEX1( NAM, MXIDAT, INVDNAM )
        ELSE
            POLLNUM = 0
        END IF

        OUTNAMES( I,1 ) = EANAM_L( I )  ! Annual emis name is same as pol/act

        L  = LEN_TRIM( ABRNAMA( I ) )

        SELECT CASE( CATEGORY_L )
          CASE( 'AREA' )
            OUTUNITS( I,1 ) = ARUNITS( 1 )
            OUTTYPES( I,1 ) = ARTYPES( 1 )

            IF( POLLNUM > 0 ) THEN
                OUTDESCS( I,1 ) = TRIM( ARDESCS( 1 ) ) // ' for ' //&
                &                  INVDDSC( POLLNUM )
            ELSE
                OUTDESCS( I,1 ) = ARDESCS( 1 )
            END IF

            DO J = 2, NPPOA
                OUTNAMES( I,J ) = ARPREFIX( J )//ABRNAMA( I )( 1:L )
                OUTUNITS( I,J ) = ARUNITS( J )
                OUTTYPES( I,J ) = ARTYPES( J )
                OUTDESCS( I,J ) = ARDESCS( J )
            ENDDO

          CASE( 'MOBILE' )

            OUTUNITS( I,1 ) = MBUNITS ( 1 )
            OUTTYPES( I,1 ) = MBTYPES ( 1 )

            IF( POLLNUM > 0 ) THEN
                OUTDESCS( I,1 ) = TRIM( MBDESCS( 1 ) ) // ' for ' //&
                &                  INVDDSC( POLLNUM )
            ELSE
                OUTDESCS( I,1 ) = MBDESCS ( 1 )
            END IF

            DO J = 2, NPPOA
                OUTNAMES( I,2 ) = MBPREFIX( 2 )//ABRNAMA( I )( 1:L )
                OUTUNITS( I,2 ) = MBUNITS ( 2 )
                OUTTYPES( I,2 ) = MBTYPES ( 2 )
                OUTDESCS( I,2 ) = MBDESCS ( 2 )
            ENDDO

          CASE( 'POINT' )

            OUTUNITS( I,1 ) = PTUNITS( 1 )
            OUTTYPES( I,1 ) = PTTYPES( 1 )

            IF( POLLNUM > 0 ) THEN
                OUTDESCS( I,1 ) = TRIM( PTDESCS( 1 ) ) // ' for ' //&
                &                  INVDDSC( POLLNUM )
            ELSE
                OUTDESCS( I,1 ) = PTDESCS( 1 )
            END IF

            DO J = 2, NPPOA
                OUTNAMES( I,J ) = PTPREFIX( J )//ABRNAMA( I )( 1:L )
                OUTUNITS( I,J ) = PTUNITS( J )
                OUTTYPES( I,J ) = PTTYPES( J )
                OUTDESCS( I,J ) = PTDESCS( J )
            ENDDO

          CASE DEFAULT
            MESG = 'INTERNAL ERROR: Do not know how to build ' //&
            &       'names for CATEGORY ' // TRIM( CATEGORY_L )
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END SELECT

    END DO

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END  SUBROUTINE BLDENAMS
