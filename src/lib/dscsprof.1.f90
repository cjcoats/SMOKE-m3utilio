
SUBROUTINE DSCSPROF( FDEV, NIPOL, EINAM )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine determines the maximum number of species per pollutant,
!      determines the maximum number of profile table entries per pollutant,
!      and creates a table that stores the species name per pollutant.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ????
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
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
!***************************************************************************
    USE M3UTILIO

!...........   Modules for public variables
!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CRL

!...........   This module contains the speciation profile tables
    USE MODSPRO, ONLY: MXSPFUL, MXSPEC, SPCNAMES, MOLUNITS

    IMPLICIT NONE

!...........   Include files

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   Subroutine arguments (note- outputs MXSPFUL, MXSPEC, and SPCNAMES
!              passed via module MODSPRO)

    INTEGER     , INTENT  (IN) :: FDEV            ! file unit number
    INTEGER     , INTENT  (IN) :: NIPOL           ! number of pollutants
    CHARACTER(*), INTENT  (IN) :: EINAM( NIPOL )  ! pollutant names

!...........   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

!.........  Local parameters
    INTEGER, PARAMETER :: MXSEG = 6        ! # of potential line segments
    INTEGER, PARAMETER :: TMPNSPEC = 5000  ! # tmp number of species per pollutant

!...........   Arrays for getting pollutant-specific information from file
    INTEGER       NENTRA ( NIPOL )   ! number of table entries per pollutant
    INTEGER       NSPECA ( NIPOL )   ! number of species per pollutant
    CHARACTER(16) POLNAMA( NIPOL )   ! unsorted pollutant names

!...........   Arrays for getting species-specific information from file
    INTEGER,            ALLOCATABLE :: INDX1A ( : )    ! sorting index for SPECNMA
    CHARACTER(NAMLEN3), ALLOCATABLE :: SPECNMA ( : )   ! unsort spcs names
    LOGICAL,            ALLOCATABLE :: LMOLAR ( : )    ! true: moles conversion is not mass
    CHARACTER(NAMLEN3)              :: TMPNAMES( TMPNSPEC,NIPOL ) ! unsort names per pollutant
    INTEGER        IPOS( 10 )       ! position in input pollutant list

!...........   Other arrays
    CHARACTER(20) SEGMENT( MXSEG )             ! Segments of parsed lines

!...........   Local variables

    INTEGER        I, J, K, L, M, N ! counters and indices
    INTEGER        ICOUNT     ! tmp counter while populating SPCNAMES
    INTEGER        INPRFTP    ! tmp. profile number
    INTEGER        IOS        ! i/o status
    INTEGER        IPOL       ! pollutant counter
    INTEGER        IREC       ! record counter
    INTEGER        ISP        ! species names counter
    INTEGER        NIPOS      ! number of pollutant matches
    INTEGER        NLINES     ! number of lines in data file
    INTEGER        PPOS       ! tmp position (from INDEX1) of pol in POLNAMA
    INTEGER        SPOS       ! tmp position (from INDEX1) of pol in SPECNMA

    REAL           FAC1, FAC2, FAC3 ! tmp speciation profile factors

    LOGICAL     :: EFLAG    = .FALSE.   ! true: error found
    LOGICAL     :: INHEADER = .FALSE.   ! true: in file header
    LOGICAL     :: BFLAG    = .FALSE.   ! true: running for biogenics

    CHARACTER(256) LINE       ! read buffer for a line
    CHARACTER(256) MESG       ! message buffer

    CHARACTER(SPNLEN3)  TMPPRF     ! tmp profile number
    CHARACTER(NAMLEN3)  POLNAM     ! tmp pollutant name
    CHARACTER(NAMLEN3)  SPECNM     ! tmp species name
    CHARACTER(SPNLEN3)  SPPRO      ! biogenics speciation profile to use

    CHARACTER(16) :: PROGNAME = 'DSCSPROF' ! program name

!***********************************************************************
!   Begin body of subroutine DSCSPROF

!...........  Make sure routine arguments are valid
    IF( FDEV .LE. 0 .OR. NIPOL .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: Invalid subroutine arguments'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!...........  Get source category, transfered by MODINFO module
    CALL GETCTGRY

!...........  For biogenics, evaluate speciation profile code to use
    IF( CRL == 'B' ) THEN
        BFLAG = .TRUE.
        MESG = 'Speciation profile to use for biogenics'
        CALL ENVSTR( 'BIOG_SPRO', MESG, ' ', SPPRO, IOS )

        IF( IOS .NE. 0 ) THEN
            MESG = 'ERROR: Variable BIOG_SPRO needs to be set ' //&
            &       'to run'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF
    END IF

!...........   Determine length of input file and allocate memory for
!              a temporary species names array, an array that
!              associates a pollutant with each species name, an
!              index array, and an array to determine output units

    NLINES = GETFLINE( FDEV, 'Speciation profile file' )

    ALLOCATE( SPECNMA( NLINES ),&
    &           INDX1A( NLINES ),&
    &           LMOLAR( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'LMOLAR', PROGNAME )

!...........   Initialize species count per pollutant and flag for indicating
!              true molar conversions (NOTE - for some pollutants like PM10,
!              there is no mole-based factor and outputu should be in units
!              of gm/mole into the mole-base speciation matrix)
    NENTRA   = 0        ! array
    NSPECA   = 0.       ! array
    POLNAMA  = ' '      ! array
    TMPNAMES = ' '      ! array
    LMOLAR   = .FALSE.  ! array

!...........   Read through input file to determine the total number
!              of pollutants in the input file, to determine the
!              number of profiles per pollutant, to store the unique
!              species names, and to store the units for mass-based and
!              mole-based conversions
    ICOUNT = 1
    IPOL   = 0
    IREC   = 0
    ISP    = 0
    DO I = 1, NLINES

        READ( FDEV,93000,END=999, IOSTAT=IOS ) LINE

        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG, 94010 )&
            &    'I/O error', IOS, 'reading speciation profile '//&
            &    'file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

!.............  Skip blank and comment lines

        IF( BLKORCMT( LINE ) ) CYCLE

!.............  Separate the line of data into each part
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

!.............  Left-justify character strings and convert factors to reals
        TMPPRF = ADJUSTL ( SEGMENT( 1 ) )

!.............  For biogenics, skip any entry that is not needed
        IF( BFLAG .AND. TMPPRF .NE. SPPRO ) CYCLE

        POLNAM = ADJUSTL ( SEGMENT( 2 ) )
        SPECNM = ADJUSTL ( SEGMENT( 3 ) )
        FAC1   = STR2REAL( SEGMENT( 4 ) )
        FAC2   = STR2REAL( SEGMENT( 5 ) )
        FAC3   = STR2REAL( SEGMENT( 6 ) )

!.............  Make sure divsor factor is not zero
        IF( FAC2 .EQ. 0. ) THEN
            WRITE( MESG,94010 ) 'WARNING: Zero divisor found at '//&
            &       'line ', IREC, '. Setting to 1.'
            CALL M3MESG( MESG )
            FAC2 = 1.
        END IF

!.............  Check width of character fields of fixed width
        L = LEN_TRIM( TMPPRF )
        IF( L .GT. SPNLEN3 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Speciation profile code ' //&
            &       'exceeds max width of', SPNLEN3, 'at line', IREC
            CALL M3MESG( MESG )
        END IF

        L = LEN_TRIM( POLNAM )
        IF( L .GT. NAMLEN3 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Pollutant name ' //&
            &       'exceeds max width of', NAMLEN3, 'at line', IREC
            CALL M3MESG( MESG )
        END IF

        L = LEN_TRIM( SPECNM )
        IF( L .GT. NAMLEN3 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Species name ' //&
            &       'exceeds max width of', NAMLEN3, 'at line', IREC
            CALL M3MESG( MESG )
        END IF

!.............  Search for pollutant in list of valid names, and go to the end
!               of the loop if none found (skip entry).  Record number
!               and position of all matches.
        M    = 0
        IPOS = 0   ! local array
        DO N = 1, NIPOL
            IF( POLNAM .EQ. EINAM( N ) ) THEN
                M = M + 1
                IF( M .LE. 10 ) THEN
                    IPOS( M ) = N
                ELSE
                    MESG = 'INTERNAL ERROR: IPOS array overflow '//&
                    &       'prevented in ' // PROGNAME
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF
            END IF
        END DO
        NIPOS = M

        IF ( MAXVAL( IPOS ) .EQ. 0 ) CYCLE

!.............  Search for pollutant unique list of all pollutants
        PPOS = INDEX1( POLNAM, IPOL, POLNAMA )

        IF ( PPOS .EQ. 0 ) THEN       ! if current POLNAM is not in
        !    POLNAMA, then
            IPOL = IPOL + 1
            POLNAMA( IPOL ) = POLNAM    ! add POLNAM to POLNAMA
            NENTRA ( IPOL ) = 1         ! init for first entry per pol

            PPOS = IPOL   ! Set for storing species count, below

        ELSE     ! if current POLNAM is already in POLNAMA, then

!.................  If a new profile number, then add to count of table entries
!                   for this pollutant
            NENTRA( PPOS ) = NENTRA( PPOS ) + 1

        END IF

        SPOS = INDEX1( SPECNM, ISP, SPECNMA )

        IF ( SPOS .LE. 0 ) THEN    ! if current SPECNM is not in
        ! SPECNMA, then
            ISP = ISP + 1
            INDX1A ( ISP )  = ISP
            SPECNMA( ISP )  = SPECNM

!.................  If mole-based = mass based, then use molar transform
            IF( FAC1/FAC2 .NE. FAC3 ) LMOLAR( ISP ) = .TRUE.

        END IF

!.............  Check if species is already stored for current pollutant, and
!               if not, increment species-per-pollutant counter and
!               add species to list.
        DO M = 1, NIPOS

            K = NSPECA( IPOS( M ) )
            J = INDEX1( SPECNM, K, TMPNAMES( 1,IPOS( M ) ) )

            IF( J .LE. 0 ) THEN

                K = K + 1

                IF( K .LE. TMPNSPEC ) THEN
                    TMPNAMES( K, IPOS( M ) ) = SPECNM

                ELSE

                    WRITE(MESG,94010)&
                    &  'INTERNAL ERROR: The', TMPNSPEC, 'species '//&
                    &  'per pollutant limit was exceeded by ' //&
                    &  'pollutant '//TRIM(POLNAM)//' in ' // PROGNAME
                    CALL M3MSG2( MESG )
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

                END IF

                NSPECA( IPOS( M ) ) = K

            END IF
        END DO

    END DO              ! End loop over speciation profile input lines

    IF( IPOL .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No pollutants found in speciation '//&
        &       'profiles match the inventory!'
        CALL M3MSG2( MESG )
    END IF

    IF( ISP .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No species found in speciation profile!'
        CALL M3MSG2( MESG )
    END IF

    IF( EFLAG ) THEN
        MESG = 'Problem(s) with speciation profiles file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!...........   Determine the max species per full table, and the max species
!              per pollutant
    MXSPFUL = MAXVAL( NENTRA )
    MXSPEC  = MAXVAL( NSPECA )

!...........   Allocate memory for species names array and units to use for
!              mole-based transformations. Also initialize.
    ALLOCATE( SPCNAMES( MXSPEC,NIPOL ),&
    &          MOLUNITS( MXSPEC,NIPOL ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MOLUNITS', PROGNAME )

    SPCNAMES = ' '   ! array
    MOLUNITS = ' '   ! array

!...........   Sort master species names
    CALL SORTIC( ISP , INDX1A, SPECNMA )  ! sort on SPECNMA

!...........   Initialize species names table and counter per pollutant
    SPCNAMES = ' ' ! array

!...........   Cycle through count of all valid pollutants (NIPOL) and all
!              species associated with these pollutants (ISP).  Check if species
!              is valid for the current pollutant, and if so, store in the
!              output species name list.
    DO I = 1, NIPOL

        ICOUNT = 0
        DO J = 1, ISP

!.................  Process species in sorted order
            K = INDX1A( J )

!.................  Find species in list of valid species per pollutant
            N = INDEX1( SPECNMA( K ), NSPECA( I ),&
            &            TMPNAMES( 1,I )            )

            IF ( N .GT. 0 ) THEN

                ICOUNT = ICOUNT + 1
                SPCNAMES( ICOUNT, I ) = SPECNMA( K )

!......................  When the species does not have molar factors, store
!                        the molar units as mass units
                IF( LMOLAR( K ) ) THEN
                    MOLUNITS( ICOUNT, I ) = SMOLUNIT
                ELSE
                    MOLUNITS( ICOUNT, I ) = SMASUNIT
                END IF

            END IF

        END DO

    END DO

!........  Rewind file

    REWIND( FDEV )

    RETURN

!......... Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
    &       'Check format of speciation' // CRLF() // BLANK5 //&
    &       'profile file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )


!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE DSCSPROF
