
PROGRAM SMKREPORT

    !***********************************************************************
    !  subroutine body starts at line 129
    !
    !  DESCRIPTION:
    !    The SMKREPORT routine create emissions and activity reports for one
    !    major source category at a time (area, mobile, or point). It permits
    !    the user to control the columns and rows in the report through a series
    !    of instructions.  These reports allow users to quality assure emissions
    !    inventories by comparison and analysis of reports from different stages
    !    of SMOKE processing. The outputs can be read into the Java Analysis and
    !    Report Tool (JART).
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Version ??/???? by ???
    !       Revised  7/2003 by A. Holland
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   MODULES for public variables
    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: QAFMTL3, NMATX, SSFLAG, SLFLAG, NREPORT, RPT_, AFLAG, ALLRPT

    !.......  This module contains report arrays for each output bin
    USE MODREPBN, ONLY: NSPCIN

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, INVPIDX, NCHARS, JSCC

    IMPLICIT NONE

    !.......   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters
    INCLUDE 'SETDECL.h90'       !  FileSetAPI variables and functions

    !.......   LOCAL PARAMETERS
    CHARACTER(50), PARAMETER :: CVSW     = '$Name SMOKEv5.0_Jun2023$'     ! CVS release tag
    CHARACTER(16), PARAMETER :: PROGNAME = 'SMKREPORT'     ! program name

    !.......   Gridding Matrix
    INTEGER, ALLOCATABLE :: GMAT( : )     ! Contiguous gridding matrix

    !.......   Speciation matrices
    INTEGER, ALLOCATABLE :: SLMAT( :,: )     ! mole-based
    INTEGER, ALLOCATABLE :: SSMAT( :,: )     ! mass-based

    !.......   File units and logical/physical names
    INTEGER :: ADEV = 0       !  ASCII elevated file
    INTEGER :: CDEV = 0       !  reports configuration file
    INTEGER :: EDEV = 0       !  elevated source ID file
    INTEGER :: GDEV = 0       !  gridding supplemental file
    INTEGER :: LDEV = 0       !  log-device
    INTEGER :: NDEV = 0       !  SCC descriptions
    INTEGER :: NIDEV = 0      !  SIC descriptions
    INTEGER :: NMDEV = 0      !  MACT descriptions
    INTEGER :: NNDEV = 0      !  NAICS descriptions
    INTEGER :: NODEV = 0      !  ORIS descriptions
    INTEGER :: MODEV = 0      !  src mapping file
    INTEGER :: PDEV = 0       !  speciation supplemental file
    INTEGER :: NPDEV = 0      !  GSPRO descriptions
    INTEGER :: RDEV(3) = (/ 0,0,0 /)     !  ASCII reports from Cntlmat program
    INTEGER :: SDEV = 0       !  ASCII inven input file
    INTEGER :: TDEV = 0       !  temporal supplemental files
    INTEGER :: YDEV = 0       !  country/state/county names file

    INTEGER, ALLOCATABLE :: ODEV( : )       !  output file unit numbers

    CHARACTER(16)  :: ANAME  = ' '       !  logical name for ASCII inven input
    CHARACTER(16)  :: ENAME  = ' '       !  logical name for I/O API inven input
    CHARACTER(16)  :: CUNAME = ' '       !  multiplicative control matrix input
    CHARACTER(16)  :: GNAME  = ' '       !  gridding matrix input
    CHARACTER(16)  :: LNAME  = ' '       !  layer fractions input file
    CHARACTER(16)  :: PRNAME = ' '       !  projection matrix input
    CHARACTER(16)  :: SLNAME = ' '       !  speciation matrix input
    CHARACTER(16)  :: SSNAME = ' '       !  speciation matrix input
    CHARACTER(16)  :: TNAME  = ' '       !  hourly emissions input file

    CHARACTER(300) :: FNAME = ' '        !  output physical/logical file name
    CHARACTER(300) :: PNAME = ' '        !  previous output file name

    !.......   Other local variables
    INTEGER      I, J,  K, L, N           ! indices and counters

    INTEGER      HWID                     ! header width
    INTEGER      IOS                      ! i/o status
    INTEGER      EDIDX                    ! ending index of loop
    INTEGER   :: GDIM    = 0              ! dimension of contiguous gridding mat
    INTEGER   :: NSLIN   = 1              ! no. mole input speciation variables
    INTEGER   :: NSSIN   = 1              ! no. mass input speciation variables

    REAL         RNFILES                  ! real number of files per report
    REAL         RNSECT                   ! real number of sections per report

    LOGICAL       :: EFLAG    = .FALSE.         ! true: error found
    LOGICAL       :: ZEROFLAG = .FALSE.         ! true: report zero values

    CHARACTER(300)     MESG                 !  message buffer
    CHARACTER(QAFMTL3) OUTFMT               !  data output format string

    !***********************************************************************
    !   begin body of program SMKREPORT

    LDEV = INIT3()

    !.......  Write out copyright, version, web address, header info, and prompt
    !           to continue running the program.
    CALL INITEM( LDEV, CVSW, PROGNAME )

    !.......  Prompt for and open REPCONFIG file
    CDEV = PROMPTFFILE( 'Enter logical name for the REPORT CONFIGURATION file',&
                        .TRUE., .TRUE., 'REPCONFIG', PROGNAME )

    !.......  Scan report configuration file to determine input file types,
    !           get global file flags and settings, and determine maximum
    !           values for use in memory allocation.
    CALL SCANREPC( CDEV )

    !.......  Get environment variable settings
    ZEROFLAG = ENVYN( 'REPORT_ZERO_VALUES', 'Leave entries ' // &
                      'with values equal to zero in reports',   &
                      .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PNAME,0,0, 'Bad env vble "REPORT_ZERO_VALUES"', 2 )
    END IF

    !.......  Prompt for and open all other input files
    CALL OPENREPIN( ENAME, ANAME, CUNAME, GNAME, LNAME,         &
                    PRNAME, SLNAME, SSNAME, TNAME, RDEV,        &
                    SDEV, GDEV, PDEV, TDEV, EDEV, YDEV, NDEV,   &
                    NIDEV, NPDEV, ADEV, NMDEV, NNDEV, NODEV )

    !.......  Read and store all report instructions
    CALL RDRPRTS( CDEV )

    !.......  Build pollutant, activity, emis-type, and species indices to output
    !           data columns based on the selected output data from the reports
    !.......  Index arrays are stored in the bin module
    CALL BLDREPIDX( SLNAME, SSNAME )

    !.......  Allocate memory for gridding matrix (even if not used so that it
    !           can be passed through subroutines)
    GDIM = NGRID + 2*NMATX
    ALLOCATE( GMAT( GDIM ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GMAT', PROGNAME )

    !.......  Allocate memory for speciation matrices (even if no speciation
    !           so that arrays can be passed through subroutines).
    N = 1
    IF( SLFLAG .OR. SSFLAG ) N = NSRC

    IF( SLFLAG ) NSLIN = NSPCIN
    ALLOCATE( SLMAT( N, NSLIN ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SLMAT', PROGNAME )

    IF( SSFLAG ) NSSIN = NSPCIN
    ALLOCATE( SSMAT( N, NSSIN ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SSMAT', PROGNAME )

    !.......  Read one-time input file data
    CALL RDREPIN( NSLIN, NSSIN, RDEV, SDEV, GDEV, PDEV, TDEV,       &
                  EDEV, YDEV, NDEV, NIDEV, NPDEV, NMDEV, NNDEV,     &
                  NODEV, ADEV, ENAME, CUNAME, GNAME, LNAME,         &
                  PRNAME, SLNAME, SSNAME, GMAT( 1 ),                &
                  GMAT( NGRID+1 ), GMAT( NGRID+NMATX+1 ),           &
                  SSMAT, SLMAT )

    !.......  Preprocess the country/state/county data
    ! note: Could add routine to reduce list of co/st/cy data to just records
    !    n: selected.

    !.......  Preprocess the inventory data
    ! note: Could add routine to reduce source to just records
    !    n: selected across all groups.

    !.......  Read and store all group definitions
    CALL RDGRPS( CDEV )

    !.......  Loop through reports
    DO N = 1, NREPORT

        RPT_ = ALLRPT( N )

        !......  Determine number of output files/sections per report

        IF( RPT_%RPTNVAR .GT. RPT_%NUMDATA ) THEN
            RPT_%RPTNVAR = RPT_%NUMDATA
        END IF

        IF( RPT_%RPTMODE .EQ. 1 ) THEN

            RPT_%NUMSECT = 1

            RNFILES = REAL( RPT_%NUMDATA ) / REAL( RPT_%RPTNVAR )

            IF( RNFILES .LT. 1.0 ) THEN

                RPT_%NUMFILES = 1

            ELSE

                RPT_%NUMFILES = INT( RNFILES )

                IF( RNFILES .GT. RPT_%NUMFILES ) THEN
                    RPT_%NUMFILES = RPT_%NUMFILES + 1
                END IF

            END IF

        ELSE IF( RPT_%RPTMODE .EQ. 2 ) THEN

            RPT_%NUMFILES = 1

            RNSECT = REAL( RPT_%NUMDATA ) / REAL( RPT_%RPTNVAR )

            IF( RNSECT .LT. 1.0 ) THEN

                RPT_%NUMSECT = 1

            ELSE

                RPT_%NUMSECT = INT( RNSECT )

                IF( RNSECT .GT. RPT_%NUMSECT ) THEN
                    RPT_%NUMSECT = RPT_%NUMSECT + 1
                END IF

            END IF

        ELSE

            RPT_%NUMFILES = 1
            RPT_%NUMSECT = 1

        END IF

        ALLRPT( N )%NUMFILES = RPT_%NUMFILES
        ALLRPT( N )%NUMSECT  = RPT_%NUMSECT
        ALLRPT( N )%RPTNVAR  = RPT_%RPTNVAR


        WRITE( MESG,94010 ) '***** CHECKING INPUTS FOR REPORT', N, ' *****'
        CALL M3MSG2( MESG )

        !.......  QA reports configuration file settings
        CALL QAREPIN( N, IOS )

        !.......  Skip report if errors are found
        IF( IOS .GT. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) '***** SKIPPING REPORT', N, ' *****'
            CALL M3MSG2( MESG )
            CYCLE
        END IF

        !.......  Write message to log and standard output for report that is
        !         being processed

        !.......  Get file name
        FNAME = RPT_%OFILENAM

        !.......  If current file is different than previous
        IF( FNAME .NE. PNAME ) THEN

            IF( ALLOCATED( ODEV ) ) DEALLOCATE( ODEV )

            !......  Allocate output file number array
            ALLOCATE( ODEV( RPT_%NUMFILES ), STAT=IOS )
            CALL CHECKMEM( IOS, 'ODEV', PROGNAME )
            ODEV = 0

            !.......  When not first report...
            IF( N .GT. 1 ) THEN

            !.......  Add Metadata to previous file if current file number is
            !         different from previous file number
            !         CALL WRMETADAT( FDEV )
            !         note: Need to write this

            !.......  Close output file(s)
                DO I = 1, RPT_%NUMFILES
                    CLOSE( ODEV( I ) )
                END DO

            END IF

            !.......  Open new output file if current file number is different
            !         previous file number.
            CALL OPENREPOUT( FNAME, ODEV, MODEV )

        END IF


        MESG = BLANK10 // 'Selecting records...'
        CALL M3MSG2( MESG )

        !.......  Select inventory records
        CALL SELECTSRC( N )

        !.......  Apply gridding information
        IF( RPT_%USEGMAT ) THEN
            CALL REPMRGGRD( N, GMAT( 1 ), GMAT( NGRID+1 ), GMAT( NGRID+NMATX+1 ), EFLAG   )
        END IF

        !.......  Skip remainder of report if error found so far
        IF( EFLAG ) CYCLE

        MESG = BLANK10 // 'Aggregating output records...'
        CALL M3MSG2( MESG )

        !.......  Assign bin numbers to selected records
        CALL ASGNBINS( N )

        MESG = BLANK10 // 'Reading emissions data and writing report...'
        CALL M3MSG2( MESG )

        !.......  Update inventory input names and units, depending on status of
        !               average day emissions.
        INVPIDX = 0
        IF ( RPT_%AVEDAY ) INVPIDX = 1
        IF( .NOT. AFLAG ) CALL GETSINFO( ENAME )
        IF( JSCC .GT. 0 ) NCHARS = NCHARS - 1          ! duplicate of rdrepin.f

        !.......  Determine input units and create conversion factors
        CALL REPUNITS( N )

        !.......  If report is multifile or database
        IF( RPT_%RPTMODE .EQ. 1 .OR.    &
            RPT_%RPTMODE .EQ. 3 .OR.    &
            RPT_%RPTMODE .EQ. 0 ) THEN
            EDIDX = RPT_%NUMFILES

        !.......  If report is multisection
        ELSE
            EDIDX = RPT_%NUMSECT

        END IF

        DO I = 1, EDIDX

            J = I
            IF( RPT_%RPTMODE .EQ. 2 ) J = 1

            !.......  Write report header
            CALL WRREPHDR( ODEV( J ), N, I, HWID, OUTFMT )

            !.......  Loop through time steps (if any) and sum emissions into bins
            !               for the appropriate time resolution...

            !.......  For mole-based speciation...
            IF( RPT_%USESLMAT ) THEN
                CALL GENRPRT( ODEV( J ), N, ADEV, MODEV,  ENAME,&
                       TNAME, LNAME, OUTFMT, SLMAT, ZEROFLAG,   &
                       EFLAG )

            !.......  For mass-based and no speciation
            ELSE
                CALL GENRPRT( ODEV( J ), N, ADEV, MODEV, ENAME, &
                       TNAME, LNAME, OUTFMT, SSMAT, ZEROFLAG,   &
                       EFLAG )
            END IF

            !.......  Save file number to use in next iteration
            PNAME  = FNAME

        END DO        ! end loop over files/sections

    END DO       ! end loop over reports

    !.......  Completion with errors
    IF( EFLAG ) THEN
        MESG = 'Problem creating reports'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Normal completion
    ELSE
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 0 )
    END IF

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

END PROGRAM SMKREPORT


