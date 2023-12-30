
SUBROUTINE GENPROJ( PDEV, PYEAR, ENAME )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine processes the projection data and writes out
    !      the projection matrix.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !************************************************************************
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
    !*************************************************************************
    USE M3UTILIO

    !.......  MODULES for public variables
    !.......  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC, CISIC, CMACT

    !.......  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: POLSFLAG, NVPROJ, FACTOR, RPTDEV,&
                        PNAMPROJ, PRJFC

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: MXCHRS, BYEAR, CATDESC, NCHARS, NSRC,&
                       SC_BEGP, SC_ENDP, CRL, CATEGORY

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN OUT) :: PDEV       ! temporary file
    INTEGER     , INTENT (IN) :: PYEAR      ! projection year
    CHARACTER(*), INTENT (IN) :: ENAME      ! emission inventory file name

    !.......  Local static arrays
    LOGICAL          LF   ( MXCHRS )          !  true: column should be output
    CHARACTER(20)    CHARS( MXCHRS )          !  source fields for output

    !.......   Local variables that depend on module variables
    CHARACTER(10) CHRHDRS( NCHARS )      ! Source characteristics headers

    !.......   Logical names and unit numbers
    !        INTEGER          ODEV           ! unit number of output tmp file
    CHARACTER(16)    PNAME          ! logical name for projection matrix

    !.......   Other local variables

    INTEGER          J, K, L, N, S, V        ! counters and indices
    INTEGER          IDUM              ! dummy integer
    INTEGER          IOS               ! i/o error status
    INTEGER          NC                ! local number src chars
    INTEGER          RDEV              ! Report unit number

    LOGICAL       :: EFLAG    = .FALSE.      ! true: error has occurred
    LOGICAL, SAVE :: FIRSTOUT = .TRUE.       ! true: output header line
    LOGICAL, SAVE :: APPLFLAG = .FALSE.      ! true: something has been applied
    LOGICAL       :: MACTFLAG = .FALSE.      ! true: MACT field is allocated

    CHARACTER(SICLEN3) :: CSIC        ! SIC code
    CHARACTER(200)     :: PATHNM      ! path name for tmp file
    CHARACTER(220)        FILENM      ! file name
    CHARACTER(512)        MESG        ! message buffer
    CHARACTER(512)        HDRBUF      ! report header buffer
    CHARACTER(NAMLEN3) :: PNAM        ! tmp pol/act name

    CHARACTER(16), PARAMETER :: PROGNAME = 'GENPROJ'     ! program name

    !***********************************************************************
    !   begin body of subroutine GENPROJ

    !.......  Get path for temporary files
    MESG = 'Path where temporary control files will be written'
    CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_TMPDIR"', 2 )
    END IF

    !.......  Open reports file
    RPTDEV( 4 ) = PROMPTFFILE( 'Enter logical name for PROJECTION REPORT',  &
                               .FALSE., .TRUE., CRL // 'PROJREP', PROGNAME )
    RDEV = RPTDEV( 4 )

    !.......  Open *output* temporary file
    ! NOTE: This is commented out because output tmp file is not used in wcntlrep.f
    !        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_proj_rep'
    !        ODEV = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )

    !.......  Set up and open output projection matrices
    CALL OPENPMAT( ENAME, BYEAR, PYEAR, PNAME )

    !......  Write header for report.
    WRITE( RDEV, 93000 ) '# Processed as '// CATDESC// ' sources'
    WRITE( RDEV, 93000 ) '# Projection factors applied with /PROJECTION/ packet'

    WRITE( RDEV, 93390 ) '#      from base year    ', BYEAR
    WRITE( RDEV, 93390 ) '#      to projected year ', PYEAR

    IF( POLSFLAG ) THEN
        WRITE( RDEV, 93000 ) '#      using pollutant-specific assignments'
    ELSE
        WRITE( RDEV, 93000 ) '#      to all pollutants uniformly'
    END IF

    !.......  Define source-category specific header
    !.......  NOTE that (1) will not be used and none will be for area sources
    CHRHDRS( 1 ) = 'Region'
    SELECT CASE( CATEGORY )
      CASE( 'AREA' )
        CHRHDRS( 2 ) = 'SCC'

      CASE( 'MOBILE' )
        CHRHDRS( 2 ) = 'Road '
        CHRHDRS( 3 ) = 'Link'
        CHRHDRS( 4 ) = 'Veh Type'
        CHRHDRS( 5 ) = 'SCC'

      CASE( 'POINT' )
        CHRHDRS( 2 ) = 'Plant ID'
        IF ( NCHARS .GE. 3 ) CHRHDRS( 3 ) = 'Char 1'
        IF ( NCHARS .GE. 4 ) CHRHDRS( 4 ) = 'Char 2'
        IF ( NCHARS .GE. 5 ) CHRHDRS( 5 ) = 'Char 3'
        IF ( NCHARS .GE. 6 ) CHRHDRS( 6 ) = 'Char 4'
        IF ( NCHARS .GE. 7 ) CHRHDRS( 7 ) = 'Char 5'
        CHRHDRS( NCHARS ) = 'SCC'

    END SELECT

    !.......  Check if MACT field is available from inventory
    MACTFLAG = ( ASSOCIATED( CMACT ) )

    !.......  Loop through all sources and store projection information for
    !           those that have it.  Otherwise, set projection factor=1.

    !.......  Initialize valid columns
    LF = .FALSE.      ! array
    DO J = 1, NCHARS
        LF( J ) = .TRUE.
    END DO

    !.......  If no pollutant-specific assignments, write out single pfac var
    IF ( .NOT. POLSFLAG ) NVPROJ = 1

    !.......  Loop through pollutants that are getting projections
    DO V = 1, NVPROJ

        IF( POLSFLAG ) THEN
            PNAM = PNAMPROJ( V )
        ELSE
            PNAM = 'pfac'
        END IF

        !.......  Loop through sources, retrieve projection packet index,
        !             set projection factor, and write out new tmp file info.
        DO S = 1, NSRC

            READ( PDEV, * ) K

            !......  If source has projection info...
            IF ( K .GT. 0 ) THEN

                !.......  Store projection factor for current pollutant/act
                FACTOR( S ) = PRJFC( K )

                !                    WRITE( ODEV,93300 ) 1, PNAM, FACTOR( S )
                APPLFLAG = .TRUE.

                !......  Write report
                !......  Format source characteristic information
                CALL PARSCSRC( CSOURC(S), MXCHRS, SC_BEGP, SC_ENDP, LF, NC, CHARS )
                NC = MIN( NC, NCHARS )

                !.......  Write out projection information for all sources
                !                   that are getting projected
                IF( ASSOCIATED( CISIC ) ) THEN
                    CSIC = CISIC( S )
                ELSE
                    CSIC = ''
                END IF

                IF( POLSFLAG ) THEN

                    !.......  Write header and use same header names as Svmkreport
                    !                           to facilitate EMF importing
                    WRITE( HDRBUF, 94015 ) '#       Variable', (CHRHDRS(N), N=1, NC), 'SIC'
                    WRITE( MESG, 94015 ) PNAM,&
                     ( CHARS(J)(1:SC_ENDP(J)-SC_BEGP(J)+1),J=1,NC ), CSIC
                ELSE
                    WRITE( HDRBUF, 94016 ) ( CHRHDRS(N), N=1, NC ), 'SIC'
                    WRITE( MESG, 94016 )&
                     ( CHARS(J)(1:SC_ENDP(J)-SC_BEGP(J)+1),J=1,NC ), CSIC
                END IF

                !.......  Add more to header and output if MACT is present
                IF( MACTFLAG ) THEN
                    L = LEN_TRIM( HDRBUF )
                    HDRBUF = HDRBUF( 1:L ) // '; MACT'

                    L = LEN_TRIM( MESG )
                    MESG = MESG( 1:L ) // '; ' // CMACT( S )
                ENDIF

                !.......  Complete header
                HDRBUF = TRIM( HDRBUF ) // '; Data value'

                IF( FIRSTOUT ) THEN
                    WRITE( RDEV, '(A)' ) TRIM( HDRBUF )
                    FIRSTOUT = .FALSE.
                ENDIF

                WRITE( RDEV, 94020 ) TRIM( MESG ), FACTOR( S )

            !.......  If source does not have projection info., set to 1.
            ELSE

                FACTOR( S ) = 1.0

                !......  Add to tmp file line
                !                    WRITE( ODEV, 93300 ) 0, 'D', 1.0

            END IF

        END DO              ! End loop through sources

        !.......  Write projection factors
        IF( .NOT. WRITESET( PNAME,PNAM,ALLFILES,0,0,FACTOR )) THEN

            MESG = 'Problem writing "'// TRIM( PNAM )// '" to "' // TRIM( PNAME )// '"'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END DO              ! End loop through pol/act

    IF( .NOT. APPLFLAG ) THEN

        MESG = 'WARNING: No PROJECTION packet entries match inventory.'
        CALL M3MSG2( MESG )

        MESG = 'WARNING: Projection matrix will not be created'
        CALL M3MSG2( MESG )

        !.......  Write not into report file
        WRITE( RDEV, 93000 ) '# No projection packet entries matched the inventory.'

        RETURN

    END IF

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx
94015 FORMAT( A16, ';', 1X, 10( A, :,';', 1X ) )

94016 FORMAT( 10( A, :, ';', 1X ) )

94020 FORMAT( A, ';', 1X, E13.5 )

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !93300   FORMAT( I2, 1X, '"', A, '"', 3( 1X, E12.5 ) )

93390 FORMAT( A, I4.4 )

    !******************  INTERNAL SUBPROGRAMS  *****************************

END SUBROUTINE GENPROJ
