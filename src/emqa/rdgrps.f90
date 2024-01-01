
SUBROUTINE RDGRPS( FDEV )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !     The RDGRPS routine reads the groups from the REPCONFIG file, interprets
    !     the entries, and stores the fully detailed group information in a table
    !     for each type of group.
    !      - Subgrid will store cell numbers that are included
    !      - Region will store Regions that are excluded.
    !C
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 7/2000 by M Houyoux
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

    !.......  MODULES for public variables
    !.......  This module contains Smkreport-specific settings
    USE MODREPRT, ONLY: NALLPCKT, NLINE_RC, REGNNAM, SUBGNAM,           &
                        NREGNGRP, NREGRAW, REG_IDX, NSUBGRID,           &
                        NSBGRAW, SBG_IDX, MXGRPREC, EXCLDRGN,           &
                        VALIDCEL, LENLAB3, PKT_IDX, LIN_DEFGRP,         &
                        GRP_LABEL, INREPORT, LIN_SUBREGN,               &
                        LIN_GROUP, GRPNRECS, GRP_INCLSTAT, RPT_IDX,     &
                        ALLRPT, RPT_, NREGREC, PKTCOUNT, LIN_SUBGRID

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, NCOLS, NROWS

    !.......  This module contains the arrays for state and county summaries
    USE MODSTCY, ONLY: NCOUNTRY, NSTATE, NCOUNTY, CTRYCOD, STATCOD, CNTYCOD

    IMPLICIT NONE

    !.......   INCLUDES
    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    INTEGER     , INTENT (IN) :: FDEV           ! output file unit number

    !....... EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL,EXTERNAL :: BLKORCMT

    !.......   Local allocatable arrays...

    !.......   Region group input allocatable arrays
    INTEGER, ALLOCATABLE :: REGNREC ( : )        ! no. records per group
    CHARACTER(FIPLEN3), ALLOCATABLE :: REGRAW  ( :,: )      ! raw codes from input file
    INTEGER, ALLOCATABLE :: REGTYPE ( :,: )      ! 0=bad,1=country,2=st,3=county
    LOGICAL, ALLOCATABLE :: REGSTAT ( :,: )      ! raw status (true=include)

    !.......   Subgrid input allocatable arrays
    INTEGER,        ALLOCATABLE :: SBGNREC ( : )       ! no. records per subgrid

    CHARACTER(100), ALLOCATABLE :: SBGRAW  ( :,: )     ! raw info from input file
    LOGICAL,        ALLOCATABLE :: SBGSTAT ( :,: )     ! raw status (true=incld)

    !.......   Per grid-cell local allocatable arrays
    LOGICAL     LCELSTAT( NGRID )
    LOGICAL     LCEL    ( NGRID )

    !.......   Per county allocatable arrays
    LOGICAL     LRGN( NCOUNTY )         ! true: county included

    !.......   Per input line local allocatable arrays
    INTEGER     LINECODE( NLINE_RC )     ! 1= in-line region; 2= in-line subgrid

    !.......   Other local arrays
    INTEGER   NCNT( NALLPCKT )       ! Count of groups defined by SELECT

    !.......   Other local variables
    INTEGER         C, I, J, K, N              ! counters and indices

    INTEGER         IC                      ! tmp partial region code
    INTEGER         IOS                     ! i/o status
    INTEGER         IREC                    ! line number
    INTEGER         LEVEL                   ! match level for region groups
    INTEGER      :: NINCL = 0               ! tmp number of included cells

    LOGICAL      :: EFLAG = .FALSE.         ! true: error found

    CHARACTER(200)  BUFFER       ! tmp label buffer
    CHARACTER(300)  MESG         ! tmp message buffer
    CHARACTER(FIPLEN3) CREGN     ! tmp partial region code

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDGRPS'     ! program name

    !***********************************************************************
    !   begin body of subroutine RDGRPS

    !.......  Allocate memory for reading and storing raw information from
    !         the input file.
    ALLOCATE( REGNNAM( NREGRAW ),       &
              SUBGNAM( NSBGRAW ), STAT=IOS )
    CALL CHECKMEM( IOS, 'REGNNAM...LRGN', PROGNAME )

    REGNNAM = ' '         ! array
    SUBGNAM = ' '         ! array
    LINECODE = 0           ! array
    LCELSTAT = .FALSE.     ! array
    LCEL     = .FALSE.     ! array
    LRGN     = .FALSE.     ! array

    !.......  Read in defined group labels
    CALL READ_GROUPS( FDEV, NGRID, 'DEFINED LABELS', LCELSTAT )

    !.......  Initialize number of groups of each type defined with SELECT
    NCNT = 0       ! array

    !.......  Read in inline group labels and compare to defined groups.  If not
    !         defined, try to match as country, state, or county and store
    !         additional names.
    CALL READ_GROUPS( FDEV, NGRID, 'SELECT LABELS', LCELSTAT )

    !.......  If error was found so far, abort
    IF( EFLAG ) THEN
        MESG = 'Problem with groups in input file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Set the number for each type of group
    NREGNGRP = NREGRAW + NCNT( REG_IDX )
    NSUBGRID = NSBGRAW + NCNT( SBG_IDX )

    !.......  If no groups are defined, leave the subroutine
    IF( NREGNGRP .EQ. 0 .AND.       &
        NSUBGRID .EQ. 0       ) RETURN

    !.......  The maximum number of records per group has been set in SCANREPC
    !         and needs to be at least 1
    MXGRPREC = MAX( MXGRPREC, 1 )

    !.......  Allocate memory for raw group information. The total number of
    !         groups is the sum of the defined groups and the unmatched and valid
    !         in-line groups.
    !.......  Reallocate memory for group labels so that labels can be reset with
    !           the inline labels as well.
    IF( NREGNGRP .GT. 0 ) THEN
        DEALLOCATE( REGNNAM )
        ALLOCATE( REGNNAM( NREGNGRP ),              &
                  REGNREC( NREGNGRP ),              &
                   REGRAW( MXGRPREC,NREGNGRP ),     &
                  REGSTAT( MXGRPREC,NREGNGRP ),     &
                  REGTYPE( MXGRPREC,NREGNGRP ),     &
                  NREGREC( NREGNGRP ),              &
                 EXCLDRGN( NCOUNTY,NREGNGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'REGNNAM...EXCLDRGN', PROGNAME )

        REGNNAM  = ' '         ! array
        REGNREC  = 0           ! array
        REGRAW   = ' '         ! array
        REGSTAT  = .TRUE.      ! array (default is include)
        REGTYPE  = 0           ! array (default is record invalid)
        EXCLDRGN = ' '         ! array

    END IF

    !.......  Same notes as above, but do for subgrids.
    IF( NSUBGRID .GT. 0 ) THEN
        DEALLOCATE( SUBGNAM )
        ALLOCATE( SUBGNAM( NSUBGRID ),              &
                  SBGNREC( NSUBGRID ),              &
                   SBGRAW( MXGRPREC,NSUBGRID ),     &
                  SBGSTAT( MXGRPREC,NSUBGRID ),     &
                 VALIDCEL( NGRID,NSUBGRID ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SUBGNAM...VALIDCEL', PROGNAME )

        SUBGNAM  = ' '         ! array
        SBGNREC  = 0           ! array
        SBGRAW   = ' '         ! 2-d array
        SBGSTAT  = .TRUE.      ! 2-d array (initialize to "include")
        VALIDCEL = 0           ! 2-d array

    END IF

    !.......  Read raw group information.  "Raw" means the entries as they
    !         appear in the input file, but not converted to the data structures
    !         needed for further processing.
    CALL READ_GROUPS( FDEV, NGRID, 'STORE', LCELSTAT )

    !.......  Convert region group entries to data structures needed for
    !           further processing
    IF( NREGNGRP .GT. 0 ) THEN

        !.......  Loop through different regions
        DO N = 1, PKTCOUNT( REG_IDX )

            !.......  Initialize region list indicator depending on first subgrid
            !                   entry
            IF( REGSTAT( 1,N ) ) THEN               ! include
                LRGN = .FALSE.                      ! array
            ELSE                                    ! exclude
                LRGN = .TRUE.                       ! array
            END IF

            !.......  Loop through records in each packet or in-line region group
            DO I = 1, REGNREC( N )

                !.......  Depending on type of the record (country, state, or
                !         county), set comparison length
                SELECT CASE( REGTYPE( I,N ) )
                  CASE( 1 )                                     ! country
                    IC = FIPEXPLEN3 + 1
                  CASE( 2 )                                     ! state
                    IC = STALEN3
                  CASE( 3 )                                     ! county
                    IC = FIPLEN3
                  CASE DEFAULT
                    CYCLE

                END SELECT

                CREGN = REGRAW( I,N )( 1:IC )

                !....... Loop through counties and determine which ones
                DO J = 1, NCOUNTY

                    IF( CNTYCOD( J )( 1:IC ) == CREGN ) THEN
                        LRGN( J ) = REGSTAT( I,N )
                    END IF

                END DO              ! End loop on counties

            END DO              ! End loop on entries in region group

            !.......  Create list of counties to exclude from group
            K = 0
            DO J = 1, NCOUNTY

                IF( .NOT. LRGN( J ) ) THEN
                    K = K + 1
                    EXCLDRGN( K,N ) = CNTYCOD( J )
                END IF
                NREGREC( N ) = K

            END DO          ! End loop on counties

        END DO              ! End loop on region groups

    END IF

    !.......  Convert subgrid entries to data structures needed for further
    !         processing.
    IF( NSUBGRID .GT. 0 ) THEN

        !.......  Loop through different subgrids
        DO N = 1, PKTCOUNT( SBG_IDX )

            !.......  Initialize cell list indicator depending on first subgrid entry
            IF( SBGSTAT( 1,N ) ) THEN               ! include
                LCEL = .FALSE.
            ELSE                                    ! exclude
                LCEL = .TRUE.
            END IF

            !.......  Loop through records in each packet or in-line subgrid
            DO I = 1, SBGNREC( N )

                !.......  Determine which cells current entry applies
                BUFFER = SBGRAW( I,N )
                NINCL = 0
                LCELSTAT = .FALSE.
                CALL PARSE_SUBGRID( BUFFER, NGRID, LCELSTAT, NINCL )

                !.......  If current entry is an include, then include cells
                IF( SBGSTAT( I,N ) ) THEN                  ! include

                    DO C = 1, NGRID
                        IF( LCELSTAT( C ) ) LCEL( C ) = .TRUE.
                    END DO

                !.......  If current entry is an exclude, then exclude cells
                ELSE                                       ! exclude

                    DO C = 1, NGRID
                        IF( LCELSTAT( C ) ) LCEL( C ) = .FALSE.
                    END DO

                END IF

            END DO

            !.......  Based on global cell status, create list of valid cells
            J = 0
            DO C = 1, NGRID
                IF( LCEL( C ) ) THEN
                    J = J + 1
                    VALIDCEL( J,N ) = C
                END IF
            END DO

        END DO

    END IF

    !.......  Deallocate raw group information
    IF( ALLOCATED( REGRAW ) ) DEALLOCATE( REGNREC, REGRAW, REGSTAT, REGTYPE )
    IF( ALLOCATED( SBGRAW ) ) DEALLOCATE( SBGRAW, SBGSTAT )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I10, :, 1X ) )

    !******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

    !.......  This subprogram provides reads the group definitions
    !               in various ways, as indicated by the main program
    !               argument.
    SUBROUTINE READ_GROUPS( FDEV, NGRID, READTYPE, LGRDSTAT )

        !.......  External functions
        LOGICAL      CHKINT
        CHARACTER(2) CRLF
        INTEGER      GETNLIST
        INTEGER      INDEX1
        INTEGER      STR2INT

       EXTERNAL     CHKINT, CRLF, GETNLIST, INDEX1, STR2INT

        !.......  Subprogram arguments
        INTEGER     , INTENT (IN) :: FDEV               ! input file unit
        INTEGER     , INTENT (IN) :: NGRID              ! no. grid cells
        CHARACTER(*), INTENT (IN) :: READTYPE           ! Reading type
        LOGICAL     , INTENT (IN) :: LGRDSTAT( NGRID )         ! true: report cell

        !.......   Local parameters
        INTEGER, PARAMETER :: MXSEG = 100

        !.......  Subprogram local arrays
        CHARACTER(256) SEGMENT( MXSEG )

        !.......  Local variables
        INTEGER       I, J, L, N               ! counters and indices

        INTEGER       IOS              ! i/o status
        INTEGER    :: NS = 1           ! no. segments in line
        INTEGER       RCNT             ! record count

        CHARACTER(FIPLEN3) CFIP         ! tmp region code
        CHARACTER(300) BUFFER           ! tmp line buffer as uppercase
        CHARACTER(300) LINE             ! tmp line buffer
        CHARACTER(300) MESG             ! mesg buffer

        CHARACTER(LENLAB3) :: PREGNNAM           ! previous region name
        CHARACTER(LENLAB3) :: PSUBGNAM           ! previous subgrid name

        !----------------------------------------------------------------------

        !.......  Rewind input file
        REWIND( FDEV )

        !.......  Loop though file to store local array of labeled group names
        SEGMENT  = ' '             ! array
        PREGNNAM = ' '
        PSUBGNAM = ' '
        IREC = 0
        DO I = 1, NLINE_RC

            READ( FDEV, 93000, END=999, IOSTAT=IOS ) LINE
            IREC = IREC + 1

            IF ( IOS .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 ) 'I/O error', IOS,               &
                  'reading report configuration file at line', IREC
                CALL M3MESG( MESG )
                CYCLE
            END IF

            !.......  Skip blank lines and comment lines
            IF( BLKORCMT( LINE ) ) CYCLE

            !.......  Screen for appended comments and remove them
            CALL RMCOMMNT( '##', LINE )

            !.......  Left-justify and convert line to upper case
            BUFFER = ADJUSTL( LINE )
            CALL UPCASE( BUFFER )

            !.......  Initialize segment from previous iteration
            SEGMENT( 1:NS ) = ' '

            !.......  Parse line into segments
            L = LEN_TRIM( BUFFER )
            NS = GETNLIST( L, BUFFER )
            IF( NS .GT. MXSEG ) NS = MXSEG
            CALL PARSLINE( BUFFER, NS, SEGMENT )

            !.......  Interpret line of code.  Set global variables in MODREPRT.
            CALL PRCLINRC( IREC, NS, BUFFER, SEGMENT )

            !.......  If a new group...
            SELECT CASE( READTYPE )

                !.......  Read the defined groups and store the labels
              CASE( 'DEFINED LABELS' )

                !.......  Get current count of current packet
                RCNT = PKTCOUNT( PKT_IDX )

                SELECT CASE( PKT_IDX )

                !.......  Store region group label
                  CASE( REG_IDX )
                    IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL

                !.......  Store subgrid label in local array
                  CASE( SBG_IDX )
                    IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL

                END SELECT

                !.......  Count SELECT statements that do not use defined regions.
                !         Determine those that coorespond to valid
                !         entries and those that should be ignored.
              CASE( 'SELECT LABELS' )

                !.......  Get current count of current packet
                RCNT = PKTCOUNT( PKT_IDX )

                SELECT CASE( PKT_IDX )

                !.......  Store region group label
                  CASE( REG_IDX )
                    IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL

                !.......  Store subgrid label in local array
                  CASE( SBG_IDX )
                    IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL

                END SELECT

                !.......  Skip if report section not started yet.
                IF( .NOT. INREPORT ) CYCLE

                !.......  A region is being selected
                IF( LIN_SUBREGN ) THEN

                    !.......  Search for region name in defined regions
                    J = INDEX1( RPT_%REGNNAM, NREGRAW, REGNNAM )

                    !.......  If region name is not found, then try to convert
                    !                           to a region code and compare with the inventory.
                    IF( J .LE. 0 ) THEN

                        !.......  Check if label is an integer
                        IF( .TRUE. .OR. CHKINT( RPT_%REGNNAM ) ) THEN

                            !.......  Convert label to region code
                            CFIP = RPT_%REGNNAM

                            !.......  Check if label is a valid region code
                            !.......  REGNTMP is a dummy argument at this stage
                            CALL CHECK_REGIONS( CFIP, LEVEL, IOS )

                            !.......  If code is valid, store line number of
                            !                                   record
                            IF( IOS .EQ. 0 ) THEN
                                N = NCNT( REG_IDX ) + 1
                                NCNT( REG_IDX ) = N
                                LINECODE( IREC ) = 1

                            !.......  Otherwise, give warning because the label
                            !                                   does not match and it's not an inline code
                            ELSE
                                WRITE( MESG,94010 ) 'WARNING: Region label "' //    &
                                  TRIM( RPT_%REGNNAM )// '" at line', IREC,         &
                                  'does not match any groups' // CRLF()//BLANK10 // &
                                  'or region codes in the inventory. '              &
                                  //'SELECT REGION will be ignored.'
                                CALL M3MSG2( MESG )

                            END IF

                        !.......  Label is not an integer, label is invalid
                        ELSE
                            L = LEN_TRIM( RPT_%REGNNAM )
                            WRITE( MESG,94010 ) 'WARNING: Region label "' //    &
                               TRIM( RPT_%REGNNAM ) // '" at line', IREC,       &
                               IREC, 'does not match any groups.'               &
                               // CRLF() // BLANK10 //                          &
                               'SELECT REGION will be ignored.'
                            CALL M3MSG2( MESG )

                        END IF              ! If label is an integer or not

                    END IF                      ! If region label not found in groups list

                !.......  Store current region label for use in next iteration
                    PREGNNAM = RPT_%REGNNAM

                END IF                         ! If region selected

                !.......  A subgrid is being selected
                IF( LIN_SUBGRID ) THEN

                    !.......  Search for subgrid names in defined subgrids
                    J = 0
                    IF( NSBGRAW .GT. 0 )   J = INDEX1( RPT_%SUBGNAM, NSBGRAW, SUBGNAM )

                    !.......  If subgrid name is not found, then try to interpret
                    !         entry as a subgrid definition.
                    IF( J .LE. 0 ) THEN

                        !....... Check if subgrid is defined in-line
                        NINCL    = 0
                        LCELSTAT = .FALSE.                           ! array
                        CALL PARSE_SUBGRID( RPT_%SUBGNAM, NGRID, LCELSTAT, NINCL )

                        !.......  If subgrid is valid, increase count and
                        !                               flag line as a in-line subgrid
                        IF( NINCL .GT. 0 ) THEN
                            N = NCNT( SBG_IDX ) + 1
                            NCNT( SBG_IDX ) = N
                            LINECODE( IREC ) = 2

                        !.......  If subgrid is invalid
                        ELSE
                            L = LEN_TRIM( RPT_%SUBGNAM )
                            WRITE( MESG,94010 )                         &
                               'WARNING: Subgrid definition "' //       &
                               RPT_%SUBGNAM( 1:L ) // '" at line',      &
                               IREC, 'is not defined and cannot' //     &
                               CRLF() // BLANK10 //                     &
                               'be interpreted. SELECT SUBGRID ' //     &
                               'will be ignored.'
                            CALL M3MSG2( MESG )

                        END IF

                    END IF                            ! If subgrid name not found in list

                    !.......  Store current subgrid label for use in next
                    !                           iteration
                    PSUBGNAM = RPT_%SUBGNAM

                END IF                            ! If subgrid selected on current line

                !.......  Store the raw information for the groups
              CASE( 'STORE' )

                !.......  If line has the group label
                IF( LIN_DEFGRP ) THEN

                    !.......  Get number of packet type
                    RCNT = PKTCOUNT( PKT_IDX )

                    SELECT CASE( PKT_IDX )

                    !.......  Store region group label
                      CASE( REG_IDX )
                        IF( LIN_DEFGRP ) REGNNAM( RCNT ) = GRP_LABEL

                    !.......  Store subgrid label in local array
                      CASE( SBG_IDX )
                        IF( LIN_DEFGRP ) SUBGNAM( RCNT ) = GRP_LABEL

                    END SELECT

                !.......  If line is inside a group packet
                ELSE IF( LIN_GROUP ) THEN

                    !.......  Get number of packet type
                    RCNT = PKTCOUNT( PKT_IDX )

                    !.......  Get records count for current packet
                    N = GRPNRECS

                    !.......  Choose group type
                    SELECT CASE( PKT_IDX )

                        !.......  Store region
                      CASE( REG_IDX )

                        !.......  Make sure region code is an integer
                        IF( .TRUE. .OR. CHKINT( SEGMENT( 1 ) ) ) THEN
                            CFIP = SEGMENT(1)

                            CALL CHECK_REGIONS( CFIP, LEVEL, IOS )

                            IF( IOS .EQ. 0 ) THEN

                                REGNREC( RCNT )  = N
                                REGRAW ( N,RCNT )= CFIP
                                REGSTAT( N,RCNT )= GRP_INCLSTAT
                                REGTYPE( N,RCNT )= LEVEL

                            END IF

                        !.......  Give an error if code not an integer
                        ELSE
                            EFLAG = .TRUE.
                            WRITE( MESG,94010 )                         &
                                   'ERROR: Region code not an ' //      &
                                   'integer in group definition ' //    &
                                   'at line', IREC
                            CALL M3MESG( MESG )

                        END IF

                        !.......  Store subgrid label in local array
                      CASE( SBG_IDX )

                        !.......  Make sure subgrid entries are specified properly
                        NINCL = 0
                        LCELSTAT = .FALSE.                   ! array
                        CALL PARSE_SUBGRID( BUFFER, NGRID, LCELSTAT, NINCL )

                        !.......  If the string could be parse as a subgrid, then
                        !         store unparsed string
                        IF( NINCL .GT. 0 ) THEN

                            SBGNREC( RCNT )   = GRPNRECS
                            SBGRAW ( N,RCNT ) = BUFFER
                            SBGSTAT( N,RCNT ) = GRP_INCLSTAT

                        END IF

                    END SELECT

                !.......  If not in group, is this line a valid Select-specified region?
                !         If so, make sure not a duplicate in-line and store it.
                ELSE IF( LINECODE( IREC ) .EQ. 1 ) THEN

                    !.......  Convert in-line code to region code and rename it.
                    CFIP = RPT_%REGNNAM
                    CALL CHECK_REGIONS( CFIP, LEVEL, IOS )
                    CALL RENAME_REGION( RPT_%REGNNAM, CFIP )

                    !.......  Get report number
                    N = PKTCOUNT( RPT_IDX )
                    ALLRPT( N )%REGNNAM = RPT_%REGNNAM

                    !.......  See if this name is already stored and if not,
                    !                           store it.
                    RCNT = PKTCOUNT( REG_IDX )
                    J = INDEX1( RPT_%REGNNAM, RCNT, REGNNAM )

                    !.......  Store without checking status because LINECODE = 1
                    !                           only if code has already been through CHECK_REGIONS
                    IF( J .LE. 0 ) THEN
                        RCNT = RCNT + 1
                        REGNNAM( RCNT ) = RPT_%REGNNAM
                        PKTCOUNT( REG_IDX ) = RCNT

                        REGNREC( RCNT )   = 1
                        REGRAW ( 1,RCNT ) = CFIP
                        REGSTAT( 1,RCNT ) = .TRUE.
                        REGTYPE( 1,RCNT ) = LEVEL
                    END IF

                    !.......  If not in group, is this line a valid Select-specified
                    !         subgrid?  If so, make sure not a duplicate in-line and
                    !         store it.
                ELSE IF( LINECODE( IREC ) .EQ. 2 ) THEN

                    !.......  Rename subgrid name
                    CALL RENAME_SUBGRID( RPT_%SUBGNAM )

                    !.......  Get report number
                    N = PKTCOUNT( RPT_IDX )
                    ALLRPT( N )%SUBGNAM = RPT_%SUBGNAM

                    !.......  See if this name is already stored and if not,
                    !        store it.
                    RCNT = PKTCOUNT( SBG_IDX )
                    J = INDEX1( RPT_%SUBGNAM, RCNT, SUBGNAM )

                    IF ( J .LE. 0 ) THEN
                        RCNT = RCNT + 1
                        SUBGNAM( RCNT ) = RPT_%SUBGNAM
                        PKTCOUNT( SBG_IDX ) = RCNT

                        SBGNREC( RCNT )   = 1
                        SBGRAW ( 1,RCNT ) = BUFFER
                        SBGSTAT( 1,RCNT ) = .TRUE.
                    END IF

                END IF              ! If not group entry or SELECT-specified

                !.......  No default case for calling internal subprogram
              CASE DEFAULT
                L = LEN_TRIM( READTYPE )
                MESG = 'INTERNAL ERROR: Can not call READ_GROUPS '//&
                       'subprogram with type "' // READTYPE( 1:L )//&
                       '"'
                CALL M3MSG2( MESG )
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

            END SELECT

        END DO

        !.......  Successful completion of subprogram
        RETURN

        !.......  Problem(s) reading input file...
999     WRITE( MESG,94010 ) 'INTERNAL ERROR: Unexpected end of ' //&
               'file at line', IREC
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        !......  FORMAT  STATEMENTS   ........

        !.......   Formatted file I/O formats...... 93xxx
93000   FORMAT( A )

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE READ_GROUPS

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This subprogram compares a region code to the valid country,
    !         state, and county codes and sets corresponding entries in a
    !         logical array aligning with the CNTYCOD array
    SUBROUTINE CHECK_REGIONS( REGN, LEVEL, STATUS )

        !.......  Exernal subroutines
        INTEGER    FINDC
       EXTERNAL   FINDC

        !.......  Subprogram arguments
        CHARACTER(FIPLEN3), INTENT (IN) :: REGN                ! region code
        INTEGER, INTENT (OUT):: LEVEL               ! sub-region level code
        INTEGER, INTENT (OUT):: STATUS              ! exit status

        !.......  Local variables
        INTEGER  K, L, N             ! counters and indices

        INTEGER  FIP                 ! tmp country/state/county code
        INTEGER  RCHK                ! region code for comparison

        CHARACTER(300) MESG          ! mesg buffer
        !----------------------------------------------------------------------

        !.......  Initialize output variables
        STATUS = 1
        LEVEL  = 0
        CALL PADZERO( REGN )

        !.......  Find in country list
        IF( REGN( FIPEXPLEN3+2:FIPLEN3 ) == '00000' ) THEN
            K = FINDC( REGN, NCOUNTRY, CTRYCOD )
            STATUS = 0
            LEVEL = 1

        !.......  Find in state list
        ELSE IF( REGN( STALEN3+1:FIPLEN3 ) == '000' ) THEN
            K = FINDC( REGN, NSTATE, STATCOD )
            STATUS = 0
            LEVEL = 2

        !.......  Find in county list
        ELSE
            K = FINDC( REGN, NCOUNTY, CNTYCOD )
            STATUS = 0
            LEVEL = 3

        END IF

        RETURN

        !......  FORMAT  STATEMENTS   ........

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE CHECK_REGIONS

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This subprogram compares interprets a subgrid definition
    !               and sets a logical array with one record per grid cell
    !               to true for the records in the subgrid.
    !.......  NINCL is > 0 when subgrid definition is valid
    !.......  LCELSTAT must be initialized by caller
    SUBROUTINE PARSE_SUBGRID( SGRANGE, NGRID, LCELSTAT, NINCL )

        !.......  External subroutines
        LOGICAL   CHKINT
        INTEGER   STR2INT

       EXTERNAL  CHKINT, STR2INT

        !.......  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: SGRANGE                   ! ASCII cell range
        INTEGER     , INTENT (IN) :: NGRID                     ! no. cells
        LOGICAL     , INTENT (OUT):: LCELSTAT( NGRID )         ! true: cell selected
        INTEGER     , INTENT (OUT):: NINCL                     ! no. cells selected

        !.......  Local variables
        INTEGER  C, C1, I, J, L, L1, L2             ! counters and indices

        INTEGER  X1, Y1                ! starting coordinate
        INTEGER  X2, Y2                ! ending corrdinate

        CHARACTER(100) BUFFER          ! input buffer
        CHARACTER(100) XBUF            ! tmp x-cell buffer
        CHARACTER(100) YBUF            ! tmp y-cell buffer
        CHARACTER(300) MESG            ! mesg buffer

        !----------------------------------------------------------------------

        !.......  Initialize number of valid cells
        NINCL = 0

        !.......  Transfer range to upper case
        BUFFER = SGRANGE
        CALL UPCASE( BUFFER )

        !.......  Find the first set of cell numbers
        L1 = INDEX( BUFFER, '(' )
        L2 = INDEX( BUFFER, ')' )
        C1 = INDEX( BUFFER, ',' )

        !.......  Extract the cell positions and compare to the valid ranges
        !.......  Give errors if bad entries for first corrdinate
        IF( L1 .GT. 0 .AND.&
            L2 .GT. 0 .AND.&
            C1 .GT. 0       ) THEN

            XBUF = BUFFER( L1+1:C1-1 )
            YBUF = BUFFER( C1+1:L2-1 )

            !.......  Ensure x-coordinate buffer is an integer
            IF( CHKINT( XBUF ) ) THEN
                X1 = STR2INT( XBUF )

                ! note:  This could be updated to only give the warning the first time
                !.......  Ensure x-coordinate range is valid
                !.......  Check minimum value and reset of out of range
                IF( X1 .LT. 1 ) THEN
                    L = LEN_TRIM( XBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting x-coordinate "' //&
                           XBUF( 1:L ) // '" at input line', IREC,&
                           'to minimum value of 1.'
                    CALL M3MESG( MESG )
                    X1 = 1

                !.......  Check maximum value and reset of out of range
                ELSE IF( X1 .GT. NCOLS ) THEN
                    L = LEN_TRIM( XBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting x-coordinate "' //&
                           XBUF( 1:L ) // '" at input line', IREC,&
                           'to grid maximum of', NCOLS
                    CALL M3MESG( MESG )
                    X1 = NCOLS

                END IF

            !.......  Give error if value is not an integer
            ELSE
                EFLAG = .TRUE.
                L = LEN_TRIM( XBUF )
                WRITE( MESG,94010 )&
                       'ERROR: Bad x-coordinate "' // XBUF( 1:L ) //&
                       '" in subgrid definition at line', IREC
                CALL M3MESG( MESG )

            END IF

            !.......  Ensure y-coordinate buffer is an integer
            IF( CHKINT( YBUF ) ) THEN
                Y1 = STR2INT( YBUF )

                !.......  Ensure y-coordinate range is valid
                !.......  Check minimum value and reset of out of range
                IF( Y1 .LT. 1 ) THEN
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting y-coordinate "' //&
                           YBUF( 1:L ) // '" at input line', IREC,&
                           'to minimum value of 1.'
                    CALL M3MESG( MESG )
                    Y1 = 1

                !.......  Check maximum value and reset of out of range
                ELSE IF( Y1 .GT. NROWS ) THEN
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting y-coordinate "' //&
                           YBUF( 1:L ) // '" at input line', IREC,&
                           'to grid maximum of', NROWS
                    CALL M3MESG( MESG )
                    Y1 = NROWS

                END IF

            !.......  Give error if value is not an integer
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                       'ERROR: Bad y-coordinate "' //TRIM( YBUF ) //&
                       '" in subgrid definition at line', IREC
                CALL M3MESG( MESG )

            END IF

        !.......  Otherwise, bad entry should be ignored
        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Could not find starting ' //&
                   'coordinate of subgrid at line', IREC
            CALL M3MESG( MESG )

        END IF

        !.......  Find the "TO" divider
        L1 = INDEX( BUFFER, 'TO' )

        !.......  If "TO" is not found, bad entry should be ignored
        IF( L1 .LE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Could not find "TO" ' //&
                   'separator of subgrid definition at line', IREC
            CALL M3MESG( MESG )

        END IF

        !.......  Remove the first cells and find the second ones
        L2 = LEN_TRIM( BUFFER )
        BUFFER = BUFFER( L1+2: L2 )

        L1 = INDEX( BUFFER, '(' )
        L2 = INDEX( BUFFER, ')' )
        C1 = INDEX( BUFFER, ',' )
        IF( L1 .GT. 0 .AND.&
            L2 .GT. 0 .AND.&
            C1 .GT. 0       ) THEN

            XBUF = BUFFER( L1+1:C1-1 )
            YBUF = BUFFER( C1+1:L2-1 )

            !.......  Ensure x-coordinate buffer is an integer
            IF( CHKINT( XBUF ) ) THEN
                X2 = STR2INT( XBUF )

                !.......  Ensure x-coordinate range is valid
                !.......  Check minimum value and reset of out of range
                IF( X2 .LT. 1 ) THEN
                    WRITE( MESG,94010 )&
                           'WARNING: resetting x-coordinate "' //&
                           TRIM( XBUF ) // '" at input line', IREC,&
                           'to minimum value of 1.'
                    CALL M3MESG( MESG )
                    X2 = 1

                !.......  Check maximum value and reset of out of range
                ELSE IF( X2 .GT. NCOLS ) THEN
                    WRITE( MESG,94010 )&
                           'WARNING: resetting x-coordinate "' //&
                           TRIM( XBUF ) // '" at input line', IREC,&
                           'to grid maximum of', NCOLS
                    CALL M3MESG( MESG )
                    X2 = NCOLS

                END IF

            !.......  Give error if value is not an integer
            ELSE
                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                       'ERROR: Bad x-coordinate "' // TRIM( XBUF ) //&
                       '" in subgrid definition at line', IREC
                CALL M3MESG( MESG )

            END IF

            !.......  Ensure y-coordinate buffer is an integer
            IF( CHKINT( YBUF ) ) THEN
                Y2 = STR2INT( YBUF )

                !.......  Ensure y-coordinate range is valid
                !.......  Check minimum value and reset of out of range
                IF( Y2 .LT. 1 ) THEN
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting y-coordinate "' //&
                           YBUF( 1:L ) // '" at input line', IREC,&
                           'to minimum value of 1.'
                    CALL M3MESG( MESG )
                    Y2 = 1

                !.......  Check maximum value and reset of out of range
                ELSE IF( Y2 .GT. NROWS ) THEN
                    L = LEN_TRIM( YBUF )
                    WRITE( MESG,94010 )&
                           'WARNING: resetting y-coordinate "' //&
                           YBUF( 1:L ) // '" at input line', IREC,&
                           'to grid maximum of', NROWS
                    CALL M3MESG( MESG )
                    Y2 = NROWS

                END IF

            ELSE
                EFLAG = .TRUE.
                L = LEN_TRIM( YBUF )
                WRITE( MESG,94010 )&
                       'ERROR: Bad y-coordinate "' // YBUF( 1:L ) //&
                       '" in subgrid definition at line', IREC
                CALL M3MESG( MESG )

            END IF

        !.......  If coordinate badly formed, bad entry should be ignored
        ELSE
            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ERROR: Could not find end ' //&
                   'coordinate of subgrid at line', IREC
            CALL M3MESG( MESG )

        END IF

        !.......  If error found, return with NINCL = 0
        IF( EFLAG ) THEN
            RETURN
        END IF

        !.......  Now loop through the x-cell range and y-cell range and set
        !         the status of the grid cells accordingly
        DO J = Y1, Y2
        DO I = X1, X2

            C = ( J - 1 ) * NCOLS + I
            LCELSTAT( C ) = .TRUE.
            NINCL = NINCL + 1

        END DO
        END DO

        RETURN

        !......  FORMAT  STATEMENTS   ........

        !.......   Internal buffering formats...... 94xxx

94010   FORMAT( 10( A, :, I10, :, 1X ) )

    END SUBROUTINE PARSE_SUBGRID

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This subprogram renames internal region names
    SUBROUTINE RENAME_REGION( NAM, CFIP )

        CHARACTER(*), INTENT (IN OUT) :: NAM
        CHARACTER(FIPLEN3), INTENT (IN    ) :: CFIP

        !----------------------------------------------------------------------

        WRITE( NAM, '(A,I12.12)' ) 'In-line region ', CFIP

    END SUBROUTINE RENAME_REGION

    !----------------------------------------------------------------------
    !----------------------------------------------------------------------

    !.......  This subprogram renames internal subgrid names
    SUBROUTINE RENAME_SUBGRID( NAM )

        CHARACTER(*), INTENT (IN OUT) :: NAM
        CHARACTER(LENLAB3) TMPNAM

    !----------------------------------------------------------------------
        TMPNAM = 'In-line region ' // TRIM( NAM )
        NAM = TMPNAM

    END SUBROUTINE RENAME_SUBGRID

END SUBROUTINE RDGRPS

