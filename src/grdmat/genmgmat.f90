
SUBROUTINE GENMGMAT( ENAME, GNAME, UNAME, FDEV, MXSCEL, MXCSRC,     &
                     MXCCL, NMSRC, NMATX, NMATXU, UFLAG, VFLAG,     &
                     DEFSRGID, FSGFLAG, NX, IX, CX, NU, IU, CU,     &
                     NCOEF, CMAX, CMIN, NCOEFU )

    !***********************************************************************
    !  subroutine body starts at line
    !
    !  DESCRIPTION:
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by M. Houyoux 5/99
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
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

    !.......   MODULES for public variables
    !.......   This module is the source inventory arrays
    USE MODSOURC, ONLY: VMT, XLOCA, YLOCA, CIFIP, CSOURC, CLINK,        &
                        IRCLAS, XLOC1, YLOC1, XLOC2, YLOC2

    !.......  This module contains the global variables for the 3-d grid
    USE MODGRID, ONLY: NGRID, GRDNM, COORD, NCOLS, NROWS, XOFF, YOFF

    !.......   This module contains the gridding surrogates tables
    USE MODSURG, ONLY: NCELLS, FIPCELL, NSRGS, SRGLIST, NSRGFIPS,       &
                       SRGFIPS, NTSRGDSC, SRGFNAM, SRGFCOD,             &
                       SRGFMT, SRGNCOLS, SRGNROWS, NTLINES,             &
                       SRGCSUM, SRGFRAC

    !.......  This module contains the lists of unique source characteristics
    USE MODLISTS, ONLY: NINVIFIP

    !.......  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NCHARS

    !.......   This module contains the cross-reference tables
    USE MODXREF, ONLY: ASRGID

    USE MODGRDLIB

    IMPLICIT NONE

    !.......   INCLUDES

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.......   SUBROUTINE ARGUMENTS
    CHARACTER(*), INTENT (IN) :: ENAME             ! inventory file name
    CHARACTER(*), INTENT (IN) :: GNAME             ! gridding mtx logical name
    CHARACTER(*), INTENT (IN) :: UNAME             ! ungridding mtx logical name
    INTEGER     , INTENT (IN) :: FDEV              ! surg codes report file
    INTEGER     , INTENT (IN) :: MXSCEL            ! max sources per cell
    INTEGER     , INTENT (IN) :: MXCSRC            ! max cells per source
    INTEGER     , INTENT (IN) :: MXCCL             ! max cells per county or link
    INTEGER     , INTENT (IN) :: NMSRC             ! no. mobile sources
    INTEGER     , INTENT (IN) :: NMATX             ! no. source-cell intersects
    INTEGER     , INTENT (IN) :: NMATXU            ! county-cell intrscts for all sources
    LOGICAL     , INTENT (IN) :: UFLAG             ! true: create gridding matrix
    LOGICAL     , INTENT (IN) :: VFLAG             ! true: using variable grid
    INTEGER     , INTENT (IN) :: DEFSRGID          ! default surrogate ID
    LOGICAL     , INTENT (IN) :: FSGFLAG           ! true: using default fallback surrogate
    INTEGER     , INTENT(OUT) :: NX  ( NGRID )     ! no. srcs per cell
    INTEGER     , INTENT(OUT) :: IX  ( NMATX )     ! src IDs
    REAL        , INTENT(OUT) :: CX  ( NMATX )     ! gridding coefficients
    INTEGER     , INTENT(OUT) :: NU  ( NMSRC )     ! no. cells per source
    INTEGER     , INTENT(OUT) :: IU  ( NMATXU )     ! cell numbers
    REAL        , INTENT(OUT) :: CU  ( NMATXU )     ! ungridding coefficients
    INTEGER     , INTENT(OUT) :: NCOEF             ! no. of gridding coeffs
    INTEGER     , INTENT(OUT) :: CMAX              ! max no. of sources per cell
    INTEGER     , INTENT(OUT) :: CMIN              ! min no. of sources per cell
    INTEGER     , INTENT(OUT) :: NCOEFU            ! no. of ungridding coeffs

    !.......   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: DSCM3GRD
    INTEGER, EXTERNAL :: GETFLINE
    LOGICAL, EXTERNAL :: BLKORCMT

    !.......   Local allocatable arrays...
    !.......   LOCAL VARIABLES and their descriptions:
    !.......   Local parameters
    INTEGER,       PARAMETER :: MXSEG = 10              ! # of potential line segments
    CHARACTER(16), PARAMETER :: PROGNAME = 'GENMGMAT'     ! program name

    !.......   Other arrays
    CHARACTER(20) SEGMENT( MXSEG )                 ! Segments of parsed lines
    CHARACTER(100), ALLOCATABLE :: TPLINE( : )     ! tmp line buffer

    !.......   Scratch Gridding Matrix (subscripted by source-within-cell, cell)

    INTEGER     IDXSRT ( NMATX )    ! sorting index
    REAL        CSJ    ( NMATX )    ! factors (with adjustments)
    INTEGER     IS     ( NMATX )    ! source IDs
    REAL        CS     ( NMATX )    ! facs (no adjustments), for ungrid
    INTEGER     NCL    ( NMATX )    ! position for county and/or link
    INTEGER     IDXSRT2( NMATX )    ! 2nd sorting index
    INTEGER     IC     ( NMATX )    ! cell IDs

    !.......   Scratch Ungridding Matrix information
    INTEGER     NNCL      ( NSRC )     ! no of county or link index by source
    INTEGER, ALLOCATABLE :: CLIDX     ( : )     ! county or link index by source
    INTEGER, ALLOCATABLE :: CNT_CL    ( : )     ! cell count per county or link
    INTEGER, ALLOCATABLE :: VMT_CELL( :,: )     ! cell numbers for county or link
    REAL   , ALLOCATABLE :: VMT_FRAC( :,: )     ! VMT fraction for cell/county or cell/link
    REAL   , ALLOCATABLE :: VMT_CL_INV( : )     ! inverse of VMT by county or link
    CHARACTER(FIPLEN3), ALLOCATABLE :: vmt_label ( : )     ! FIPS codes

    !.......   Temporary array for flagging sources that are outside the domain
    LOGICAL     INDOMAIN( NMSRC )     ! true: src is in the domain

    !.....   Arrays for links intersecting with cells
    !.....   Note that the NGRID dimension could conceivably be too small if
    !        a link winds through the whole domain, but this is a case that
    !        is not worth going to extra trouble for.
    INTEGER     ACEL( NGRID )        ! number of cell intersctns per src
    REAL        AFAC( NGRID )        ! fraction of link in cell

    CHARACTER(FIPLEN3)  FIPNOSRG( NINVIFIP )      ! cy/st/co codes w/o surrogates
    CHARACTER(SRCLEN3)    LKOGRD( NMSRC )     ! link srcs outside the grid

    !.......   Temporary arrays for storing surrogate codes to use
    INTEGER     SURGID1( NSRC )     ! primary surrogate code
    INTEGER     SURGID2( NSRC )     ! secondary surrogate code

    !.......   Other local variables

    INTEGER         C, CNT, F, I, II, IDX, J, JJ, K, KK, N, NT, S     !  indices and counters
    INTEGER         LC, LN                  !  "previous" values of these

    INTEGER         GDEV          !  for surrogate coeff file
    INTEGER         CNTMAX        !  counter for storing correct max dimensions
    INTEGER         COL           !  tmp column
    INTEGER         TCOL          ! tmp column
    INTEGER         ID1, ID2      !  tmp primary and secondary surg codes
    INTEGER         IOS           !  i/o status
    INTEGER         ISIDX         !  surrogate ID code index
    INTEGER         NFIP          !  no of FIPS list from sources
    INTEGER         ISDEF         !  default surrogate ID code index
    INTEGER         L2            !  string length
    INTEGER         LNKEND        !  width of sourc info to end of link ID
    INTEGER         NCEL          !  tmp number of cells
    INTEGER         NCOULNK       !  number of counties and links
    INTEGER         NLKOGRD       !  no. of link sources outside the grid
    INTEGER         NNOSRG        !  no. of cy/st/co codes with no surrogates
    INTEGER         ROW           ! tmp row
    INTEGER         TROW          ! tmp row
    INTEGER         RWT           !  tmp roadway type
    INTEGER         NTL           ! max no. of line buffers
    INTEGER         TGTSRG        ! target surrogates code
    INTEGER         SSC           ! surrogates code
    INTEGER      :: NLINES = 0    ! number of lines in input file

    REAL            ADJ, ADJC       ! tmp adjustment factors
    REAL            ALEN            ! link length
    REAL            FRAC            ! tmp surrogate fraction
    REAL            XBEG, YBEG      ! tmp X and Y link start coordinates
    REAL            XEND, YEND      ! tmp X and Y link end   coordinates
    REAL            SUM
    REAL*8          XX, YY

    LOGICAL      :: EFLAG = .FALSE.        ! true: error detected
    LOGICAL      :: IFLAG = .FALSE.        ! true: internal error detected
    LOGICAL      :: LFLAG = .FALSE.        ! true: location data available
    LOGICAL      :: XYSET = .FALSE.        ! true: X/Y available for src
    LOGICAL      :: LASTIME = .FALSE.      ! true: X/Y available for src
    LOGICAL      :: CFLAG   = .FALSE.      ! true: called by sizgmat, false: called by gen[a|m]gmat
    LOGICAL      :: WFLAG = .FALSE.        ! true: per iteration warning flag

    CHARACTER(FIPLEN3) CFIP       !  tmp country/state/county code
    CHARACTER(FIPLEN3) LFIP       !  cy/st/co code from previous iteration
    CHARACTER(FIPLEN3) LLFIP      !  cy/st/co code from previous iteration
    CHARACTER(200)  LINE          !  Read buffer for a line
    CHARACTER(16)   COORUNIT      !  coordinate system projection units
    CHARACTER(80)   GDESC         !  grid description
    CHARACTER(256)  BUFFER        !  source fields buffer
    CHARACTER(256)  MESG          !  message buffer
    CHARACTER(256)  NAMBUF        !  surrogate file name buffer
    CHARACTER(256)  NAMBUFT       !  tmp surrogate file name buffer
    CHARACTER(256)  TSRGFNAM      !  tmp surrogate file name buffer

    CHARACTER(LNKLEN3)    CLNK       ! tmp link ID
    CHARACTER(LNKLEN3)    LLNK       ! previous link ID
    CHARACTER(SRCLEN3)    CSRC       ! tmp source chars string

    !***********************************************************************
    !   begin body of subroutine GENMGMAT

    !.......  Get grid name from the environment and read grid parameters
    IF( .NOT. DSCM3GRD( GRDNM, GDESC, COORD, GDTYP3D, COORUNIT,         &
                        P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D,    &
                        XORIG3D, YORIG3D, XCELL3D, YCELL3D,             &
                        NCOLS3D, NROWS3D, NTHIK3D ) ) THEN

        MESG = 'Could not get Models-3 grid description'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    !.......  Store grid parameters for later processing
    ELSE

        IF( NCOLS3D * NROWS3D .NE. NGRID ) THEN

            MESG = 'INTERNAL ERROR: Number of cells in "' //        &
                   TRIM( PROGNAME ) // '" are inconsistent '//      &
                   'with calling program.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

        END IF

    END IF

    !.......  If ungridding matrix is needed, allocate memory for VMT,
    !           and read the VMT data.
    IF ( UFLAG ) THEN
        ALLOCATE( VMT( NSRC ), STAT=IOS )
        CALL CHECKMEM( IOS, 'VMT', PROGNAME )
        CALL RDMAPPOL( NSRC, 1, 1, 'VMT', VMT )
    END IF

    !.......  Initialize number of sources per cell counter
    NX = 0       ! Array

    !.......  Allocate memory for temporary gridding matrix and other
    IC  = 0
    IS  = 0
    NCL = 0
    NNCL =0
    CS  = 0.
    CSJ = 0.
    FIPNOSRG = ' '
    SURGID1 = 0       ! array
    SURGID2 = 0       ! array

    LLNK    = ' '
    LLFIP   = ' '
    LFIP    = ' '
    N       = 0

    !.......  Store the number and values of unfound cy/st/co codes
    !.......  Keep track of sources that are outside the domain
    DO I = 1, NSRC

        CFIP = CIFIP( I )
        CLNK = CLINK ( I )

        !.......  Count and store the number of county and links
        IF ( CFIP .NE. LLFIP .OR. CLNK .NE. LLNK ) N = N + 1
        NNCL( I ) = N

        !.......  storing st/cy/ct index within a domain
        F = FINDC( CFIP, NSRGFIPS, SRGFIPS )

        IF ( F .LE. 0 ) THEN

            IF( CFIP .NE. LFIP ) THEN
                NNOSRG = NNOSRG + 1
                FIPNOSRG( NNOSRG ) = CFIP
                LFIP = CFIP
            END IF

            INDOMAIN( I ) = .FALSE.
            CYCLE       ! To next source

        END IF

        LLFIP = CFIP
        LLNK  = CLNK

    END DO

    !.......  Set flag to indicate that XLOCA/YLOCA are available
    LFLAG = ALLOCATED( XLOCA )

    MESG = 'Computing gridding matrix and statistics...'
    CALL M3MSG2( MESG )

    J       = 0
    LFIP    = ' '
    LLNK    = ' '
    LNKEND  = VIDPOS3 - 1
    NNOSRG  = 0
    NLKOGRD = 0
    ADJ     = 1.        !  temporary
    ADJC    = 1.        !  temporary

    !.......  default fallback surrogate will run at last after
    !           re-assigned zero fraction surrogate

    DO II = 1, NSRGS      ! loop through only the surrogate code assigned by sources

        TGTSRG = SRGLIST( II )
        ISDEF  = FIND1( DEFSRGID, NSRGS, SRGLIST )
        KK = II

        IF( FSGFLAG ) THEN

            !.......  default fallback surrogate will run at last after
            !                   re-assigned zero fraction surrogate
            IF( II >= ISDEF ) THEN
                IF( II == NSRGS ) THEN
                    TGTSRG = DEFSRGID
                    KK = ISDEF
                ELSE
                    TGTSRG = SRGLIST( II + 1 )
                    KK = II + 1
                END IF
            END IF

        END IF

        NTL    = NTLINES( KK )

        !.......  Warning message when there are no surrogate available due to out of domain
        IF( NTL .EQ. 0 ) THEN
            WRITE( MESG,94010 ) 'WARNING: The surrogate', TGTSRG,       &
                ' does not exist within a domain range: emissions'//    &
                ' of sources assigned to', TGTSRG,                      &
                ' will be set to zero'
            CALL M3MSG2( MESG )

            SRGFRAC = 0.0
            SRGCSUM = 0.0

        ELSE

            WRITE( MESG,94010 ) 'Looping surrogate', TGTSRG,    &
                ' through sources to generate gridding matrix'
            CALL M3MSG2( MESG )

            !.......  Allocate memory for indices to surrogates tables for each source
            ALLOCATE( TPLINE( NTL ), STAT=IOS )
            CALL CHECKMEM( IOS, 'TPLINE', PROGNAME )

            !.......  If surrogates are needed, read and store the gridding surrogates,
            !         allocate memory for the surrogate assignments, and assign
            !         surrogates to each source.
            NT = 0
            TPLINE = ' '
            TSRGFNAM = ' '

            DO I = 1, NTSRGDSC      ! Open all surrogate files using the same srg code

                !.......  Prompt for and open I/O API output file(s)...
                CALL GETENV( 'SRGPRO_PATH', NAMBUF )
                NAMBUFT = TRIM( NAMBUF ) // SRGFNAM( I )

                IF( TGTSRG .NE. SRGFCOD( I ) ) CYCLE

                IF( NAMBUFT .NE. TSRGFNAM  ) THEN
                    CALL OPEN_SRGFILE

                    !.......  Reading surrogate files
                    DO JJ = 1, NLINES

                        READ ( GDEV, 93000, END=111, IOSTAT=IOS )LINE

                        IF ( BLKORCMT( LINE ) ) CYCLE

                        CALL PARSLINE( LINE, MXSEG, SEGMENT )
                        SSC    = STR2INT ( SEGMENT( 1 ) )
                        TCOL   = STR2INT ( SEGMENT( 3 ) )
                        TROW   = STR2INT ( SEGMENT( 4 ) )

                        !......  Check the value of the column number
                        IF( TCOL .LT. 0 .OR.  TCOL .GT. SRGNCOLS .OR.   &
                          ( TROW .EQ. 0 .AND. TCOL .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                        END IF

                        !......  Check the value of the row number
                        IF( TROW .LT. 0 .OR.  TROW .GT. SRGNROWS .OR.   &
                          ( TCOL .EQ. 0 .AND. TROW .NE. 0 ) ) THEN
                            WFLAG = .TRUE.
                            CALL M3MESG( MESG )

                        !......  Special treatment for cell (0,0) (skip for now)
                        ELSE IF( TROW .EQ. 0 .AND. TCOL .EQ. 0 ) THEN
                            CYCLE

                        END IF

                        !......  Adjust column and row for subgrid
                        TCOL = TCOL - XOFF
                        TROW = TROW - YOFF

                        !......  Skip entry after subgrid adjustment
                        IF( TCOL .LE. 0 .OR. TCOL .GT. NCOLS .OR.   &
                            TROW .LE. 0 .OR. TROW .GT. NROWS ) CYCLE

                        !......  Skip entry if rows and columns are out of range
                        IF( WFLAG ) CYCLE

                        !.......  Skip entry if SSC is not in the assigned SRGLIST by source
                        IF( SSC .EQ. TGTSRG ) THEN
                            NT = NT + 1
                            TPLINE( NT ) = LINE
                        END IF

111                 END DO

                    TSRGFNAM = NAMBUFT        ! store a previous surrogate file name buffer
                    CLOSE( GDEV )

                END IF          ! skip if surrogate file has the same srg file

            END DO                   ! loop over all surrogate files in SRGDESC file

            !.......  Populating assigned surrogated
            CALL RDSRG4GRD( NT, TPLINE, CFLAG )             ! populating surrogates

            DEALLOCATE( TPLINE )

        END IF           ! end of populating surrogates

        !.......    Compute gridding matrix:
        !.......       First case:   explicit link (ILINK > 0
        !.......       Second case:  some LNKDEF entry applies
        !.......       Third case:   FIP+roadtype cross-reference match
        !.......       Fourth case:  state+roadtype cross-reference match
        !.......       fifth case:   roadtype cross-reference match
        !.......       sixth case:   fallback default

        DO S = 1, NSRC

            CFIP = CIFIP ( S )
            RWT  = IRCLAS( S )
            CLNK = CLINK ( S )
            CSRC = CSOURC( S )
            SSC  = ASRGID( S )

            IF( SSC .NE. TGTSRG ) CYCLE

            !.......  Determine if x/y location is available
            XYSET = .FALSE.
            IF( LFLAG ) XYSET = ( XLOCA( S ) .GT. AMISS3 )

            !.......  Special case for source has an x/y location
            IF( XYSET ) THEN

                !.......  If source is in the domain, get cell number and store
                XX = XLOCA( S )
                YY = YLOCA( S )
                IF( INGRID( XX, YY, NCOLS, NROWS, COL, ROW  ) ) THEN

                    C = ( ROW-1 ) * NCOLS + COL
                    NX( C ) = NX( C ) + 1
                    J = J + 1

                    !.......  Check that the maximum number of sources per cell is ok
                    IF ( J .LE. NMATX ) THEN
                        IS ( J ) = S
                        IC ( J ) = C
                        NCL( J ) = NNCL( S )
                        CS ( J ) = 1.0
                    END IF

                !......  Otherwise, mark source as being outside domain
                ELSE
                    INDOMAIN( S ) = .FALSE.

                END IF

                CYCLE               ! To head of loop over sources

            END IF       ! End if assigned point location or not

            !.......  Find FIP/RoadWayType adjustment for FIP and RWT
            !         ADJ = ADJMV( NADJ1, FIP, RWT, ADJFIP, ADJRWT, ADJFAC1 )

            !.......  Process for link source...
            IF ( CLNK .NE. ' ' ) THEN

                !.......  Convert coordinate system from UTM to required system

                XBEG = XLOC1 ( S )
                YBEG = YLOC1 ( S )
                XEND = XLOC2 ( S )
                YEND = YLOC2 ( S )

                !.......  Compute the fractions of the link in each grid cell
                CALL LNK2GRD( NGRID, XBEG, YBEG, XEND, YEND,&
                              NCEL, ACEL, AFAC, ALEN, IFLAG  )

                !.......  Make sure that there was enough storage
                IF ( IFLAG ) THEN
                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    WRITE( MESG,94010 )&
                       'INTERNAL ERROR: Overflow for link to ' //&
                       'grid conversion for:' // CRLF()// BLANK10//&
                       BUFFER( 1:L2 )

                    CALL M3MESG( MESG )
                    CYCLE

                !.......  If link is outside the grid, store its information and
                !         go to next loop iteration
                ELSE IF( NCEL .EQ. 0 ) THEN

                    NLKOGRD = NLKOGRD + 1
                    LKOGRD( NLKOGRD ) = CSRC( 1:LNKEND )

                    CYCLE

                !.......  Write warning if the link has zero-length link
                ELSE IF( ALEN .EQ. 0.0 ) THEN

                    CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                    WRITE( MESG,94010 )&
                      'WARNING: Zero-length link for: ' // CRLF() //&
                       BLANK10 // BUFFER( 1:L2 ) // '.'// CRLF() //&
                       BLANK10 // 'Emissions put in cell', ACEL( 1 )
                    CALL M3MESG( MESG )

                END IF

                !.......  Loop through cells intersecting the current link
                DO K = 1, NCEL

                    C = ACEL( K )

                    !......  Find cell-based adjustment for cell C
                    !......  Warn when both types of adjustments are applied

                    !.......  Make sure cell ID is in valid range.
                    IF( C .GT. NGRID .OR. C .LT. 0 ) THEN
                        EFLAG = .TRUE.
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                        WRITE( MESG,94010 ) 'Tried to assign link to illegal cell =', C, 'for:' //       &
                               CRLF() // BLANK5 // BUFFER( 1:L2 )
                        CALL M3MESG( MESG )
                        CYCLE                        ! to next cell in link

                    !.......  Otherwise, increase the count of records and
                    !.......  store the source count for this cell
                    ELSE
                        J = J + 1
                        NX( C ) = NX( C ) + 1

                    END IF

                    !.......  Check that the maximum number of src-cell intersections
                    !         is okay
                    IF ( J .LE. NMATX ) THEN

                        IDXSRT( J ) = J
                        IS    ( J ) = S
                        IC    ( J ) = C
                        NCL   ( J ) = NNCL( S )
                        CS    ( J ) = AFAC( K )
                        CSJ   ( J ) = ADJ * ADJC * AFAC( K )
                    END IF

                END DO      ! end loop over cells for this link

                CYCLE       ! to next source

            END IF          ! end of link-specific processing

    !....... Start of non-link processing (will only get here if the
    !                  source is a non-link source)...

    !....... Look for source in the grid-by-link table
    !                S = FIND2( FIP, RWT, NUMLNK, LNKFIP, LNKRWT )      !  use LNKDEF factor
    !                IF ( S .GT. 0 ) THEN

    !                    DO K = 1, LNKCNT( S )
    !                        C = LNKCEL( J,S )
    !                        J = NX( C ) + 1

    !.......  Find cell-based adjustment for cell C
    !                        ADJC = ADJMV( NADJ2, C, 0, ADJCELL, ADJCELL,
    !                                     ADJFAC2 )

    !.......  Warn when both types of adjustments are applied
    !                        IF( ADJ .NE. 1. .AND. ADJC .NE. 1. ) THEN

    !                            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
    !                            WRITE( MESG,94010 ) 'WARNING: Both FIP/'//
    !                                  'RoadWayType and cell-specific adjustments'
    !                                  // ' applied for ' // CRLF() // BLANK10//
    !                                  BUFFER( 1:L2 ) // ' Cell:', C
    !                            CALL M3MSG2( MESG )

    !                        END IF

    !.......  Check that the maximum number of sources per cell is ok
    !                        IF ( J .LE. MXSCEL ) THEN
    !                            IS ( J,C ) = S
    !                            CS ( J,C ) = LNKFAC( K,S )
    !                            CSJ( J,C ) = ADJ * ADJC * LNKFAC( K,S )
    !                        END IF

    !.......  Store the count of sources for current cell
    !                        NX( C )   = J

    !                    END DO       ! end of loop on cells for this source

    !                    CYCLE        ! to next source

    !                END IF           ! end of grid-by-link processing

            !.......  Process for non-link, non-grid-by-link sources...

            !.......  Retrieve the indices to the surrogates tables
            ISIDX = 1
            F    = FINDC( CFIP, NSRGFIPS, SRGFIPS )

            !.......  Store the number and values of unfound cy/st/co codes
            !.......  Keep track of sources that are outside the domain
            IF ( F .LE. 0 ) THEN

            !.......  Re-assigning org assigned srg to default fallback srg
                IF( FSGFLAG ) ASRGID( S ) = DEFSRGID
            !                    IF ( FIP .NE. LFIP .OR. CLNK .NE. LLNK ) N = N - 1

                CYCLE               ! To next source

            END IF

            !.......  Loop through all of the cells intersecting this co/st/cy code.
            DO K = 1, NCELLS( F )

                C = FIPCELL( K,F )               ! Retrieve cell number

                    !.......  Set the surrogate fraction
                CALL SETFRAC( S, ISIDX, TGTSRG, K, F, 2,&
                             INDOMAIN( S ), CSRC, DEFSRGID, FSGFLAG,&
                             ID1, ID2, FRAC, CFLAG )

                    !.......  Re-assigning org assigned srg to default fallback srg
                IF( ID2 .EQ. DEFSRGID .AND. FSGFLAG ) THEN
                    ASRGID( S ) = DEFSRGID
                    CYCLE
                END IF

                    !.......  Store surg IDs for reporting
                SURGID1( S ) = ID1
                SURGID2( S ) = ID2

                IF( FRAC .GT. 0 ) THEN

                    J = J + 1                                  ! increment no. src-cell intersection
                    NX( C ) = NX( C ) + 1                      ! increment no. srcs per cell

                    !.......  Find cell-based adjustment for cell C
                    !.......  Warn when both types of adjustments are applied

                    !.......  Check that the maximum number of src-cell intersections is okay
                    IF ( J .LE. NMATX ) THEN
                        IDXSRT( J ) = J
                        IS    ( J ) = S
                        IC    ( J ) = C
                        NCL   ( J ) = NNCL( S )
                        CS    ( J ) = FRAC
                        CSJ   ( J ) = ADJ * ADJC * FRAC
                    END IF

                END IF      ! if surrogate fraction > 0.

            END DO        !  end of loop on cells K for this FIP

            LFIP = CFIP
            LLNK = CLNK

        END DO            !  end loop on sources S, computing gridding matrix.

    END DO            !  end loop on surrogate code.

    !.......  For first time routine is called in all cases,
    IF( LASTIME ) THEN
        MESG = 'WARNING: Some surrogates renormalized when ' //     &
               'total of surrogates by county were ' // CRLF() //   &
                BLANK10 // 'greater than 1.'
        CALL M3MSG2( MESG )
    ENDIF

    !.......  Abort if error
    IF( EFLAG .OR. IFLAG ) THEN
        MESG = 'Problem creating gridding matrix.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    IF( NLKOGRD .NE. 0 ) THEN

        WRITE( MESG,94010 ) 'NOTE: ', NLKOGRD, 'link sources were outside grid.'
        CALL M3MSG2( MESG )

    END IF

    NCOEF = J
    NCOULNK = N

    !.......  Abort if overflow occurred
    IF ( NCOEF .GT. NMATX ) THEN

        WRITE( MESG,94010 )                                                     &
          'INTERNAL ERROR: Gridding and ungridding matrices not written.' //    &
          CRLF() // BLANK10 // 'Arrays would have overflowed.' //               &
          CRLF() // BLANK10 // 'Current maximum cell-source intersections (NMATX) =',&
          NMATX, '.' //                                                         &
          CRLF() // BLANK10 // 'Actual  maximum cell-source intersections         =',&
          NCOEF, '.'
        CALL M3MSG2( MESG )
        CALL M3EXIT( PROGNAME, 0, 0, 'Error-exit', 2 )

    END IF

    !.......  Sort the scratch gridding matrix arrays to organize by cell and by
    !         source.
    CALL SORTI2( NCOEF, IDXSRT, IC, IS )

    !.......  Compress scratch gridding matrix into output representation
    !.......  and compute statistics.  Already NX part of gridding matrix

    CMAX = 0
    CMIN = 99999999
    DO C = 1, NGRID

        J = NX( C )
        IF (      J .GT. CMAX ) THEN
            CMAX = J
        ELSE IF ( J .GT. 0 .AND. J .LT. CMIN ) THEN
            CMIN = J
        END IF

    END DO

    DO K = 1, NCOEF      ! Loop through cell-src intersections

        J = IDXSRT( K )
        IX( K ) = IS ( J )
        CX( K ) = CSJ( J )

    END DO

    !.......  Write gridding matrix
    MESG = 'Writing out GRIDDING MATRIX file...'
    CALL M3MSG2( MESG )

    IF( .NOT. WRITE3( GNAME, 'ALL', 0, 0, NX ) ) THEN
        MESG = 'Error writing GRIDDING MATRIX file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.......  Write output surrogates codes
    DO S = 1, NSRC
        WRITE( FDEV,93360 ) S, SURGID1( S ), SURGID2( S )
    END DO

    !..........
    !......   Generate ungridding matrix ......
    !..........

    !......  If we need to create the ungridding matrix...
    IF( UFLAG ) THEN

        DO J = 1, NMATX
            IDXSRT2( J ) = J
        END DO

        !.......  Sort the scratch gridding matrix arrays to organize by county/link,
        !               cell, and source
        CALL SORTI3( NCOEF, IDXSRT2, NCL, IC, IS )

        !.......  Allocate memory for county/link VMT within grid
        ALLOCATE( CLIDX( NSRC ),            &
                 CNT_CL( NCOULNK ),         &
               VMT_CELL( NCOULNK,MXCCL ),   &
               VMT_FRAC( NCOULNK,MXCCL ),   &
             VMT_CL_INV( NCOULNK ),         &
              vmt_label( NCOULNK ), STAT=IOS )
        CALL CHECKMEM( IOS, 'vmt_label', PROGNAME )
        CLIDX       = 0        ! array
        CNT_CL      = 0        ! array
        VMT_CELL    = 0        ! array
        VMT_FRAC    = 0.       ! array
        VMT_CL_INV  = 0.           ! array
        vmt_label = ' '

        !.......  Compute county/link VMT total within grid
        DO K = 1, NCOEF
            J = IDXSRT2( K )
            S = IS ( J )
            N = NCL( J )
            IF ( N .EQ. 0 ) CYCLE           ! skip records outside the grid

            !.......  County and/or link total VMT
            if ( s .gt. 0 ) then
                vmt_label( n ) = cifip( s )
            else
                vmt_label( n ) = ' '
            endif
            VMT_CL_INV( N ) = VMT_CL_INV( N ) + VMT(S) * CS(J)

        END DO               ! end loop to compute VMT totals

        !.......  Inverse of county and/or link total VMT
        DO N = 1, NCOULNK
            IF( VMT_CL_INV( N ) .GT. 0. ) VMT_CL_INV( N ) =&
                                               1. / VMT_CL_INV( N )
        END DO

        !.......  Loop through cell-source intersections and compute county total
        !               VMT by cell/source over County total VMT
        LN = -1
        LC = -1
        CNTMAX = 0
        DO K = 1, NCOEF

            J   = IDXSRT2( K )
            C   = IC( J )
            S   = IS( J )
            N   = NCL( J )
            IF ( N .EQ. 0 .OR. C .LE. 0 ) CYCLE      ! skip records outside the grid

            CLIDX( S ) = N

            IF ( N .NE. LN ) THEN
                CNT = 1
            ELSE IF ( C .NE. LC ) THEN
                CNT = CNT + 1
            END IF

            IF( CNT .LE. MXCCL ) THEN
                CNT_CL  ( N )     = CNT
                VMT_CELL( N,CNT ) = C
                VMT_FRAC( N,CNT ) =  VMT_FRAC( N,CNT ) +&
                                     VMT(S) * CS(J) * VMT_CL_INV(N)
            END IF
            IF( CNT > CNTMAX ) THEN
                CNTMAX = CNT
            END IF
            LN = N
            LC = C
        END DO

        IF( CNTMAX .GT. MXCCL ) THEN
            WRITE( MESG,94010 )&
             'INTERNAL ERROR: Ungridding matrix not written'&
             // CRLF() // BLANK10 // 'Arrays would have overflowed.'&
             // CRLF() // BLANK10 // 'Current maximum cells per county/link (MXCCL) =',&
             MXCCL, '.' // CRLF() // BLANK10 //&
             'Actual  maximum cells per county/link    =',&
             CNTMAX, '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )
        END IF

        !.......  Create ungridding matrix. The factors have been created so that any
        !               non-link sources within a county have the same factors used for
        !               computing temperatures.  This means that a source (e.g., urban
        !               interstates) that appear in only some cells in the county will
        !               use temperatures that are based on all sources in the county).
        !               The link sources within a county have different factors than the
        !               non-link sources, and these use only the cells that the link intersects.

        K  = 0
        DO S = 1, NSRC
            N = CLIDX( S )

            !......  Skip sources that are outside the grid
            IF( N .EQ. 0 ) THEN
                NU( S ) = 0
                CYCLE

            !......  Store the number of cells per source as the number of cell
            !        intersections with the county or link
            ELSE
                NU( S ) = CNT_CL( N )
            END IF

            sum = 0
            !......  Store the cell numbers and VMT fractions into the ungridding matrix
            DO I = 1, CNT_CL( N )
                K = K + 1
                IU( K ) = VMT_CELL( N,I )
                CU( K ) = VMT_FRAC( N,I )
            END DO

        END DO
        NCOEFU = K

        !.......  Check for overflow
        IF( NCOEFU .GT. NMATXU ) THEN

            WRITE( MESG,94010 ) 'INTERNAL ERROR: Ungridding matrix not written' //  &
             CRLF() // BLANK10 // 'Arrays would have overflowed.' //                &
             CRLF() // BLANK10 // 'Current cell-source intersections (NMATXU) =', NMATXU, '.' // &
             CRLF() // BLANK10 // 'Actual  cell-source intersections          =', NCOEFU, '.'
            CALL M3MSG2( MESG )
            CALL M3EXIT( PROGNAME, 0, 0, ' ', 2 )

        END IF

        !.......  Write out ungridding matrix

        MESG = 'Writing out UNGRIDDING MATRIX file...'
        CALL M3MSG2( MESG )

        IF ( .NOT. WRITE3( UNAME, 'ALL', 0, 0, NU ) ) THEN
            MESG = 'Problem writing UNGRIDDING MATRIX file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

    END IF

    !.......  Report FIPS that don't have surrogate data
    !.......  Report links that are outside the grid
    !        CALL RPSRCOUT( NNOSRG, NLKOGRD, FIPNOSRG, LKOGRD )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats...... 93xxx

93000 FORMAT( A )

93360 FORMAT( I8, 1X, I8, 1X, I8 )

    !.......   Internal buffering formats...... 94xxx

94010 FORMAT( 10( A, :, I9, :, 1X ) )


    !******************  INTERNAL SUBPROGRAMS  *****************************
CONTAINS

    !.......  This internal subprogram opens individual surrogate file

    SUBROUTINE OPEN_SRGFILE

        !.......  Get the number of lines in the surrogate description file desription file
        GDEV = GETEFILE( NAMBUFT,.TRUE., .TRUE., PROGNAME )

        IF( GDEV .LT. 0 ) THEN
            MESG = 'Could not open input surrogate file' // TRIM( NAMBUFT )
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

        REWIND( GDEV )

        RETURN

    END SUBROUTINE OPEN_SRGFILE

END SUBROUTINE GENMGMAT

