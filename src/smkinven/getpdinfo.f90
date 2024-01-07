
SUBROUTINE GETPDINFO( FDEV, TZONE, INSTEP, OUTSTEP, TYPNAM,    &
                      FNAME, SDATE, STIME, NSTEPS, NPDVAR,     &
                      NPDVSP, MXPDSRC, EAIDX, SPIDX, CFLAG )

    !***************************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !      This subroutine gets the vital information from the day-specific or
    !      hour-specific input files so that memory can be allocated to read in
    !      the data.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !      Subroutines: I/O API subroutine
    !
    !  REVISION  HISTORY:
    !       Created 12/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***************************************************************************
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
    !.......   This module is the inventory arrays
    USE MODSOURC, ONLY: CINTGR, INTGRFLAG

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NIPPA, NSPDAT, EANAM, NCOMP, VAR_FORMULA,  VNAME

    !.......  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: NUNIQCAS, UCASNKEP, UNIQCAS, UCASIDX, ITNAMA, SCASIDX

    IMPLICIT NONE

    !.......  INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......  SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN):: FDEV              ! file unit no.
    INTEGER     , INTENT (IN):: TZONE             ! output time zone
    INTEGER     , INTENT (IN):: INSTEP            ! expected data time step HHMMSS
    INTEGER     , INTENT (IN):: OUTSTEP           ! output time step HHMMSS
    CHARACTER(*), INTENT (IN):: TYPNAM            ! name of processing type
    CHARACTER(*), INTENT (IN):: FNAME             ! logical file name
    INTEGER     , INTENT(OUT):: SDATE             ! Julian start date in TZONE
    INTEGER     , INTENT(OUT):: STIME             ! start time of data in TZONE
    INTEGER     , INTENT(OUT):: NSTEPS            ! no. time steps
    INTEGER     , INTENT(OUT):: NPDVAR            ! no. pol/act variables
    INTEGER     , INTENT(OUT):: NPDVSP            ! no. pol/act/special data
    INTEGER     , INTENT(OUT):: MXPDSRC           ! max. no. srcs over all times
    INTEGER     , INTENT(OUT):: EAIDX( NIPPA )    ! index to EANAM
    INTEGER     , INTENT(OUT):: SPIDX( MXSPDAT )    ! index to SPDATNAM
    LOGICAL     , INTENT(OUT):: CFLAG             ! true: CEM data processing

    !.......   EXTERNAL FUNCTIONS and their descriptions:

    INTEGER, EXTERNAL :: GETFLINE
    INTEGER, EXTERNAL :: GETFORMT

    !.......  Local parameters

    CHARACTER(16), PARAMETER :: PROGNAME = 'GETPDINFO'     !  program name
    CHARACTER(16), PARAMETER :: FORMEVNM = 'SMKINVEN_FORMULA'

    !.......  Local arrays...
    INTEGER         EASTAT( NIPPA )       ! true: act/pol present in data
    INTEGER         SPSTAT( MXSPDAT )         ! true: special data variable used

    INTEGER, ALLOCATABLE :: MXPDPT( : )         ! true: special data variable used

    !.......  Other local variables
    INTEGER         I, J, N, V, NV, IV, L, LL              ! counters and indices
    INTEGER         CIDX, NCIDX            ! tmp data index
    INTEGER         FILFMT               ! format code of files in list
    INTEGER         INVFMT               ! inventory format code
    INTEGER         IOS                  ! i/o status
    INTEGER         NLINE                ! number of lines
    INTEGER         NPPCAS               !  no. of pollutants per CAS number

    LOGICAL       :: DFLAG    = .FALSE.      ! true: day-specific processing

    CHARACTER(NAMLEN3) INVNAM       ! temporary pollutant name
    CHARACTER(NAMLEN3) POLNAM       ! temporary pollutant name
    CHARACTER(300)     MESG            !  message buffer

    !***********************************************************************
    !   begin body of program GETPDINFO

    EASTAT = 0      ! array
    SPSTAT = 0      ! array

    MESG = 'Determining number of time steps for ' // TYPNAM // '-specific files...'
    CALL M3MSG2( MESG )

    !.......  Check no of new calculated pollutants
    CALL ENVSTR( FORMEVNM, MESG, ' ', VAR_FORMULA, IOS )
    IF( LEN_TRIM( VAR_FORMULA ) > 0 ) CALL FORMLIST

    !.......  Perform case-specific settings
    SELECT CASE( TYPNAM )
      CASE( 'day' )
        DFLAG = .TRUE.

      CASE( 'hour' )
        DFLAG = .FALSE.

      CASE DEFAULT
        MESG = 'INTERNAL ERROR: Do not know type ' // TYPNAM
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END SELECT

    !.......  Determine whether combining VOC and HAPs together or not
    DO V = 1, NIPPA
        IF( INDEX( EANAM(V),'_NOI'  ) > 0 ) INTGRFLAG = .TRUE.
        IF( INDEX( EANAM(V),'NONHAP') > 0 ) INTGRFLAG = .TRUE.
    END DO

    !.......  Ensure that input file is a list-formatted file
    INVFMT = GETFORMT( FDEV, -1 )

    IF( INVFMT .NE. LSTFMT ) THEN
        MESG = TYPNAM// '-specific input file is not provided by '//    &
               'a list of files OR ' // CRLF() // BLANK10 //            &
               'files in list provided could not be found.'
        CALL M3MSG2( MESG )

        MESG = 'Problem reading inventory file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Get the dates (in the output time zone) from the files,
    !         flag the pollutants of interest, and flag the special variables
    !         contained in the file.
    CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG,    &
                   FNAME, SDATE, STIME, NSTEPS, FILFMT,             &
                   EASTAT, SPSTAT )

    !.......  Initialize for the maximum number of records per time step
    ALLOCATE( MXPDPT( NSTEPS ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MXPDPT', PROGNAME )
    MXPDPT = 0      ! array

    MESG = 'Determining number of sources for ' // TYPNAM // '-specific files...'
    CALL M3MSG2( MESG )

    !.......  Get the maximum number of records per time step - i.e., populate
    !         MXPDSRC
    CALL RDLOOPPD( FDEV, TZONE, INSTEP, OUTSTEP, MXPDSRC, DFLAG,    &
                   FNAME, SDATE, STIME, NSTEPS, FILFMT,             &
                   EASTAT, SPSTAT )

    !.......  Check whether processing CEM dataset or not
    NV = 0
    DO V = 1, NIPPA
        NV = NV + EASTAT( V )
    END DO
    IF( NV == NIPPA ) CFLAG = .TRUE.

    N = 0
    DO V = 1, NIPPA

        !.......  Processing CEM data
        IF( CFLAG ) THEN
            N = N + 1
            EAIDX( N ) = V
            CYCLE
        END IF

        !.......  Add multiple inventory pollutant(s) with same CAS name
        !         Find code corresponding to current pollutant before you add
        IF( EASTAT( V ) > 0 ) THEN
            N = N + 1
            EAIDX( N ) = V
            CIDX   = EASTAT( V )
            NPPCAS = UCASNKEP( CIDX )
            DO J = 2, NPPCAS
                NCIDX   = UCASIDX( CIDX ) + J - 1
                POLNAM = ITNAMA( SCASIDX( NCIDX ) )
                NV = INDEX1( POLNAM, NIPPA, EANAM )
                IF( INDEXINT1( NV, NIPPA, EAIDX ) < 1 .AND. NV > 0 ) THEN
                    N = N + 1
                    EAIDX( N ) = NV
                END IF
            END DO
        END IF

    END DO

    !......  Add new computed pollutants
    DO I = 1, NCOMP
        POLNAM = VNAME( I )
        NV = INDEX1( POLNAM, NIPPA, EANAM )
        IV = FIND1( NV, N, EAIDX )
        IF( NV > 0 .AND. IV < 1 ) THEN        ! only add if it doesn't exit in PDAY
            N = N + 1
            EAIDX( N ) = NV
        END IF
    END DO

    NPDVAR = N

    !.......  Create index to special data variable names for current data files
    !.......  The idex serves a different purpose from EAIDX and is constructed
    !           differently intentionally.
    N = 0
    DO V = 1, MXSPDAT

        IF( SPSTAT( V ) > 0 ) THEN
            N = N + 1
            SPIDX( V ) = N
        END IF

    END DO
    NSPDAT = N

    NPDVSP = NPDVAR + NSPDAT

    !.......  Compute the maximum number of sources per time step
    !.......  NOTE - MXPDPT is in the MODDAYHR module
    MXPDSRC = MAXVAL( MXPDPT )

    !.......  If no sources matched then error
    IF ( MXPDSRC .EQ. 0 ) THEN

        MESG = 'No ' // TYPNAM //'-specific sources matched ' //    &
               'the inventory for the time period processed.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    DEALLOCATE( MXPDPT )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE GETPDINFO


