
SUBROUTINE RDGREF( FDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     Reads the gridding cross-reference file for area or mobile sources.  It
    !     allocates memory (locally) for reading the unsorted x-refs. It sorts the
    !     x-refs for processing. It allocates memory for the appropriate x-ref
    !     tables and populates the tables (passed via modules).
    !
    !  PRECONDITIONS REQUIRED:
    !     File unit FDEV already is opened... MORE
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 4/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
    !       related changes
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

    !.......  MODULES for public variables
    !.......  This module is for cross reference tables
    USE MODXREF, ONLY: INDXTA, CSRCTA, CSCCTA, ISRGCDA

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NCHARS, SC_ENDP

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER, INTENT (IN) :: FDEV       ! cross-reference file unit no.

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: CHKINT
    LOGICAL, EXTERNAL :: BLKORCMT
    INTEGER, EXTERNAL :: GETFLINE

    !.......   Local parameters
    INTEGER, PARAMETER :: AREATYP  = 1
    INTEGER, PARAMETER :: BIOGTYP  = 2
    INTEGER, PARAMETER :: MOBILTYP = 3

    CHARACTER(6), PARAMETER :: LOCCATS( 3 ) =&
                           (/ 'AREA  ', 'BIOG  ', 'MOBILE' /)

    !.......   Array of input fields
    CHARACTER(SCCLEN3)  FIELDARR( 3 )

    !.......   Other local variables
    INTEGER         I, J, J1, J2, L, N        !  counters and indices

    INTEGER         FIP         !  temporary FIPS code
    INTEGER         IDUM        !  dummy variable
    INTEGER         IOS         !  i/o status
    INTEGER         IREC        !  record counter
    INTEGER         ISRG        !  tmp surrogate ID
    INTEGER         JS          !  position of SCC in source chars in x-ref file
    INTEGER         LINTYPE     !  temporary source category code
    INTEGER         LPCK        !  length of point definition packet
    INTEGER      :: NCP = 0     !  input point source header parm
    INTEGER         NLINES      !  number of lines
    INTEGER         NXREF       !  number of valid x-ref entries
    INTEGER         THISTYP     !  index in LOCCATS for CATEGORY
    INTEGER         VTYPE       !  temporary vehicle type number

    LOGICAL      :: EFLAG = .FALSE.       !  true: error found
    LOGICAL      :: LDUM  = .FALSE.       !  dummy
    LOGICAL      :: SKIPREC = .FALSE.     !  true: record skipped in x-ref file

    CHARACTER(2)       SCC2         !  1st  2nd character of SCC
    CHARACTER(300)     LINE         !  line buffer
    CHARACTER(300)     MESG         !  message buffer
    CHARACTER(ALLLEN3) CSRCALL      !  buffer for source char, incl pol
    CHARACTER(FIPLEN3) CFIP         !  buffer for CFIPS code
    CHARACTER(RWTLEN3) CRWT         !  buffer for roadway type
    CHARACTER(SICLEN3) CDUM         !  dummy buffer for SIC code
    CHARACTER(MACLEN3) CDUM2        !  dummy buffer for MACT code
    CHARACTER(SCCLEN3) CHKZERO      !  buffer to check for zero SCC
    CHARACTER(SCCLEN3) TSCC         !  temporary SCC or roadway type
    CHARACTER(VIDLEN3) CVID         !  buffer for vehicle type ID

    CHARACTER(16) :: PROGNAME = 'RDGREF'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDGREF

    !.......  Ensure that the CATEGORY is valid
    !.......  Use THISTYP to  et LINTYPE when aren't sure if the record is for
    !           current source category or not.
    THISTYP = INDEX1( CATEGORY, 3, LOCCATS )

    IF( THISTYP .LE. 0 ) THEN
        MESG = 'INTERNAL ERROR: category "' // TRIM( CATEGORY ) //&
               '" is not valid in routine ' // PROGNAME
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ENDIF

    !.......  Get the number of lines in the file
    NLINES = GETFLINE( FDEV, 'Gridding cross reference file' )

    !.......  Allocate memory for unsorted data used in all source categories
    ALLOCATE( ISRGCDA( NLINES ),        &
               CSCCTA( NLINES ),        &
               CSRCTA( NLINES ),        &
               INDXTA( NLINES ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ISRGCDA...INDXTA', PROGNAME )

    !.......  Set up constants for loop.

    !.......  Second pass through file: read lines and store unsorted data for
    !           the source category of interest
    IREC   = 0
    N      = 0
    CDUM  = ' '
    CDUM2 = ' '
    DO I = 1, NLINES

        READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF ( IOS .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,           &
                'reading GRIDDING X-REF file at line', IREC
            CALL M3MESG( MESG )
            CYCLE
        END IF

        !.......  Skip blank or comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

        !.......  Depending on source category, transfer line to temporary
        !               fields.  In cases where blanks are allowed, do not use
        !               STR2INT to prevent warning messages.
        SELECT CASE( CATEGORY )

          CASE( 'AREA' )
            CALL PARSLINE( LINE, 3, FIELDARR )

            CFIP = FIELDARR( 1 )
            TSCC = FIELDARR( 2 )

            !.......  Post-process x-ref information to scan for '-9', pad
            !                   with zeros, compare SCC version master list.
            CALL FLTRXREF( CFIP, CDUM, TSCC, ' ', CDUM2,        &
                           IDUM, IDUM, IDUM, LDUM, SKIPREC )

          CASE( 'MOBILE' )
            CALL PARSLINE( LINE, 3, FIELDARR )

            L = LEN_TRIM( FIELDARR( 2 ) )

            CFIP = FIELDARR( 1 )
            TSCC = FIELDARR( 2 )

            !.......  Post-process x-ref information to scan for '-9', pad
            !                   with zeros.  Do not include SCC in call below because
            !                   right SCC will not work.
            CALL FLTRXREF( CFIP, CDUM, TSCC, ' ', CDUM2,        &
                           IDUM, IDUM, IDUM, LDUM, SKIPREC )

        END SELECT

        !.......  Make sure that the spatial surrogates code is an integer
        IF( .NOT. CHKINT( FIELDARR( 3 ) ) ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Spatial surrogates ' // &
                   'code is not an integer at line', IREC
            CALL M3MESG( MESG )
        END IF

        !.......  If this record is in error, or should be skipped because
        !               it doesn't match any sources, go to next iteration
        IF( EFLAG .OR. SKIPREC ) CYCLE

        !.......  Convert surrogate code to an integer
        ISRG = STR2INT( FIELDARR( 3 ) )

        N = N + 1
        IF( N .GT. NLINES ) CYCLE          ! Ensure no overflow

        !.......  Store case-indpendent fields
        INDXTA ( N ) = N
        ISRGCDA( N ) = ISRG
        CSCCTA ( N ) = TSCC

        !.......  Store sorting criteria as right-justified in fields
        CSRCALL = ' '
        IF( CATEGORY .EQ. 'AREA' ) THEN
            CALL BLDCSRC( CFIP, TSCC, CHRBLNK3, CHRBLNK3,       &
                          CHRBLNK3, CHRBLNK3, CHRBLNK3,         &
                          POLBLNK3, CSRCALL   )
        ELSE
            CALL BLDCSRC( CFIP, RWTBLNK3, CHRBLNK3, CHRBLNK3,   &
                          CHRBLNK3, CHRBLNK3, CHRBLNK3,         &
                          POLBLNK3, CSRCALL   )
        END IF

        CSRCTA( N ) = CSRCALL( 1:SC_ENDP( NCHARS ) ) // TSCC

    END DO          ! End of loop on I for reading in speciation x-ref file

    !.......  Reset number of cross-reference entries in case some were dropped
    NXREF = N

    !.......  Write errors for problems with input
    IF( NXREF .EQ. 0 ) THEN
        EFLAG = .TRUE.
        MESG = 'ERROR: No valid gridding cross-reference entries    !'
        CALL M3MSG2( MESG )

    ELSEIF( NXREF .GT. NLINES ) THEN
        EFLAG = .TRUE.
        WRITE( MESG,94010 ) 'INTERNAL ERROR: dimension for ' //     &
               'storing gridding cross-reference was', NLINES,      &
               CRLF() // BLANK10 // 'but actually needed', NXREF
        CALL M3MSG2( MESG )

    ENDIF

    !.......  Check for errors reading XREF file, and abort
    IF( EFLAG ) THEN
        MESG = 'Problem reading gridding cross-reference file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    CALL M3MSG2( 'Processing gridding cross-reference file...' )

    CALL SORTIC( NXREF, INDXTA, CSRCTA )

    !.......  Group cross-reference data into tables for different groups
    CALL XREFTBL( 'GRIDDING', NXREF )

    !.......  Deallocate other temporary unsorted arrays
    DEALLOCATE( ISRGCDA, CSRCTA, CSCCTA, INDXTA )

    !.......  Rewind file
    REWIND( FDEV )

    RETURN

    !.......  Error message for reaching the end of file too soon
999 MESG = 'End of file reached unexpectedly. ' //&
           'Check format of gridding' // CRLF() // BLANK5 //&
           'cross-reference file.'
    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE RDGREF
