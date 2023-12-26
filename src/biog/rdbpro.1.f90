
SUBROUTINE RDBPRO( FDEV, SPPRO, NIPOL, EINAM, MSPCS, EMSPC,&
&                   MLFAC, MSFAC )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!
!  Find speciation profile to use for speciation of biogenic emissions.
!  Calculate mole and mass factors for each species based on this profile.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       11/99 by Jeff Vukovich
!       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
!
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
!***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!...........   Include files

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   EXTERNAL FUNCTIONS and their descriptions:
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

!...........   Subroutine arguments (note- outputs MXSPFUL, MXSPEC, and SPCNAMES
!              passed via module MODSPRO)

    INTEGER     , INTENT  (IN) :: FDEV            ! file unit number
    INTEGER     , INTENT  (IN) :: NIPOL           ! number of pollutants
    INTEGER     , INTENT  (IN) :: MSPCS           ! max number of species
    CHARACTER(*), INTENT  (IN) :: EINAM( NIPOL )  ! pollutant names
    CHARACTER(*), INTENT  (IN) :: SPPRO           ! biogenic profile to find
    CHARACTER(*), INTENT  (IN) :: EMSPC( MSPCS )

    REAL, INTENT (OUT) :: MLFAC( MSPCS, NIPOL )
    REAL, INTENT (OUT) :: MSFAC( MSPCS, NIPOL )

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDBPRO' ! program name
    INTEGER,       PARAMETER :: MXSEG = 6        ! # of potential line segments

!...........   Other arrays

    CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines

!...........   Local variables

    INTEGER        I, J, K    ! counters and indices
    INTEGER        INPRFTP    ! tmp. profile number
    INTEGER        IOS        ! i/o status
    INTEGER        IREC       ! record counter
    INTEGER        NLINES     ! number of lines in data file
    INTEGER        PPOS       ! tmp position (from INDEX1) of pol in POLNAMA
    INTEGER        SPOS       ! tmp position (from INDEX1) of pol in SPECNMA

    REAL           SPLTFAC, SDIV, SMFAC ! tmp speciation profile factors
    LOGICAL        FOUND_FLAG             ! flag for requested profile
    LOGICAL        EFLAG             ! error flag

    CHARACTER(SPNLEN3)   TMPPRF     ! tmp profile number
    CHARACTER(16)  POLNAM     ! pollutant name
    CHARACTER(16)  SPECNM     ! tmp species name
    CHARACTER(300) LINE       ! buffer for profile data
    CHARACTER(256) MESG       ! message buffer

!***********************************************************************
!   Begin body of subroutine RDBPRO

!...........  Make sure routine arguments are valid
    IF( FDEV .LE. 0 .OR. NIPOL .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: Invalid subroutine arguments'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    EFLAG = .FALSE.
    FOUND_FLAG = .FALSE.
    MLFAC = 0.0   ! array
    MSFAC = 0.0   ! array

!...........   Determine length of input file

    NLINES = GETFLINE( FDEV, 'Speciation profile file' )

!...........   Initialize species count per pollutant and flag for indicating
!              true molar conversions (NOTE - for some pollutants like PM10,
!              there is no mole-based factor and outputu should be in units
!              of gm/mole into the mole-base speciation matrix)

!...........   Read through input file to determine the total number
!              of pollutants in the input file, to determine the
!              number of profiles per pollutant, to store the unique
!              species names, and to store the units for mass-based and
!              mole-based conversions
    IREC   = 0
    DO I = 1, NLINES

        READ( FDEV,93100,END=999,IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF( IOS .GT. 0 ) THEN
            WRITE( MESG, 94010 )&
            &    'I/O error', IOS, 'reading speciation profile '//&
            &    'file at line', IREC
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Skip blank and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Separate the line of data into each part
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

!.............  Left-justify character strings and convert factors to reals
        TMPPRF = ADJUSTL ( SEGMENT( 1 ) )

        IF( TMPPRF .NE. SPPRO ) CYCLE

        FOUND_FLAG = .TRUE.

        POLNAM  = ADJUSTL ( SEGMENT( 2 ) )
        SPECNM  = ADJUSTL ( SEGMENT( 3 ) )
        SPLTFAC = STR2REAL( SEGMENT( 4 ) )
        SDIV    = STR2REAL( SEGMENT( 5 ) )
        SMFAC   = STR2REAL( SEGMENT( 6 ) )

!.............  Search for pollutant in list of valid names, and go to the end
!               of the loop if not found (skip entry)

        J = INDEX1( POLNAM, NIPOL, EINAM )

        IF( J .EQ. 0 ) CYCLE

        SPOS = INDEX1( SPECNM, MSPCS, EMSPC )

        IF ( SPOS .GT. 0 ) THEN

            IF ( SDIV .EQ. 0 ) THEN
                EFLAG = .TRUE.
                CYCLE
            ENDIF

            MLFAC( SPOS, J ) = SPLTFAC / SDIV
            MSFAC( SPOS, J ) = SMFAC * GM2TON

        END IF


    END DO

    IF( .NOT. FOUND_FLAG ) THEN
        MESG =  SPPRO // ' profile not found in GSPRO'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( EFLAG ) THEN
        MESG = 'At least one of the divisors was zero in GSPRO.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF


    RETURN

!......... Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of speciation' // CRLF() // BLANK5 //&
    &       'profile file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )


!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93100 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDBPRO
