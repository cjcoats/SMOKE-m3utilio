
SUBROUTINE WRPDEMIS( DAYFLAG, JDATE, JTIME, TIDX, NPDSRC, NVAR,&
&                     NVASP, FNAME, PFLAG, CFLAG, EAIDX, SPIDX,&
&                     LASTSTEP, PDIDX, PDDATA, EFLAG, MXPDSRC )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!      Subroutines: I/O API subroutine
!
!  REVISION  HISTORY:
!      Created 12/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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

!.........  MODULES for public variables
!...........   This module is the inventory arrays
    USE MODSOURC, ONLY: CSOURC, CINTGR, INTGRFLAG

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: INVSTAT, MXIDAT, INVDNAM, INVDVTS,&
    &                    ITMSPC, ITEXPL, ITNAMA, NINVTBL

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: NSRC, NCHARS, NIPPA, EANAM, NCOMP, VAR_FORMULA,&
    &                   CHKPLUS, CHKMINUS, VIN_A, VIN_B, VNAME

!.........  This module contains data for day- and hour-specific data
    USE MODDAYHR, ONLY: PDTOTL, NPDPT, IDXSRC, SPDIDA, CODEA,&
    &                    EMISVA, DYTOTA, CIDXA

!.........  This module contains the lists of unique inventory information
    USE MODLISTS, ONLY: FIREFLAG, NUNIQCAS, UNIQCAS, UCASIDX, SCASIDX,&
    &                    UCASNKEP, ITNAMA, ITFACA, ITKEEPA

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!.........  SUBROUTINE ARGUMENTS
    LOGICAL     , INTENT  (IN) :: DAYFLAG              ! true: day-, false: hour-spec
    INTEGER     , INTENT  (IN) :: JDATE                ! Julian date
    INTEGER     , INTENT  (IN) :: JTIME                ! time HHMMSS
    INTEGER     , INTENT  (IN) :: TIDX                 ! time index
    INTEGER     , INTENT  (IN) :: NPDSRC               ! no. part-day srcs
    INTEGER     , INTENT  (IN) :: NVAR                 ! no. pol/act vars
    INTEGER     , INTENT  (IN) :: NVASP                ! no. pol/act/special
    CHARACTER(*), INTENT  (IN) :: FNAME                ! output file name
    LOGICAL     , INTENT  (IN) :: PFLAG                ! true: gen profiles
    LOGICAL     , INTENT  (IN) :: CFLAG                ! true: CEM processing
    INTEGER     , INTENT  (IN) :: EAIDX( NVAR )        ! pol/act index
    INTEGER     , INTENT  (IN) :: SPIDX( MXSPDAT )     ! special var index
    LOGICAL     , INTENT  (IN) :: LASTSTEP             ! true: last timestep
    INTEGER     , INTENT (OUT) :: PDIDX ( NPDSRC )     ! sparse src index
    REAL        , INTENT (OUT) :: PDDATA( NPDSRC,NVASP)! sparse data storage
    LOGICAL     , INTENT (OUT) :: EFLAG                ! true: error found
    INTEGER     , INTENT  (IN) :: MXPDSRC              ! maximum period-specific sources

!...........   Local allocatable arrays
    INTEGER,       ALLOCATABLE, SAVE :: EAIDX2( : )    ! reverse index for EAIDX
    LOGICAL,       ALLOCATABLE, SAVE :: NOMISS( :,: )
    INTEGER,       ALLOCATABLE, SAVE :: SPIDX2( : )
    INTEGER,       ALLOCATABLE, SAVE :: NHAPPOS( : )   ! positions of NONHAPVOCs
    CHARACTER(3),  ALLOCATABLE, SAVE :: NHAPMOD( : )

!...........   Local arrays
    INTEGER         SIDX( NPDSRC )      ! start index of a source
    INTEGER         EIDX( NPDSRC )      ! end index of a source

!...........   LOCAL PARAMETERS
    CHARACTER(16), PARAMETER   :: FORMEVNM = 'SMKINVEN_FORMULA'

!...........   Other local variables
    INTEGER          I, J, K, LK, LN, M, N, S, V, V2, NC, NV, IV, CV
    INTEGER          II, JJ, VV, F, NP, L, LL, L2, LS

    INTEGER          IOS                  ! i/o status
    INTEGER, SAVE :: NWARN = 0            ! warning count
    INTEGER, SAVE :: NWARN1= 0            ! warning count
    INTEGER, SAVE :: MXEA                 ! maximum pol/var # in EAIDX
    INTEGER, SAVE :: MXWARN               ! max no. warnings
    INTEGER          NOUT                 ! tmp no. sources per time step
    INTEGER          NPPCAS               ! no. of pollutants per CAS number
    INTEGER          IDXA                 ! position of first variable in source index
    INTEGER          IDXB                 ! position of second in iable in source index
    INTEGER          NVPOS                ! position of calculated new variable in source index
    INTEGER, SAVE :: NNHV = 0             ! no of NONHAPVOCs
    INTEGER          NVOC                 ! no of NONHAPVOCs
    INTEGER          NHAP                 ! no of HAPs
    INTEGER          NHVPOS               ! position of NONHAPVOC in PTDAY
    INTEGER          NVRAW                ! raw position of variable in source index
    INTEGER          WARNCNT              ! number of times warnings

    LOGICAL       :: FND_VINA = .FALSE.   ! true: found formula pollutant
    LOGICAL       :: FND_VINB = .FALSE.  ! true: found formula pollutant
    LOGICAL       :: DUPFLAG  = .FALSE.  ! true: record is an actual duplicate
    LOGICAL, SAVE :: FIRSTIME = .TRUE.   ! true: first time routine called
    LOGICAL, SAVE :: HOURFLAG = .FALSE.  ! true: hour-spec
    LOGICAL, SAVE :: DFLAG    = .FALSE.  ! true: error on duplicates
    LOGICAL, SAVE :: LFLAG    = .FALSE.  ! true: iteration on special var

    CHARACTER(3  )   TMPMOD           ! tmp mode (EXH,EVP,,,,)
    CHARACTER(6  )   TYPE             ! "Hourly" or "Daily"
    CHARACTER(256)   BUFFER           ! src description buffer
    CHARACTER(500)   MESG             ! message buffer
    CHARACTER(NAMLEN3) INVNAM, POLNAM ! tmp SMOKE name

    CHARACTER(16) :: PROGNAME = 'WRPDEMIS' !  program name

!***********************************************************************
!   begin body of program WRPDEMIS

!.........  If there is one or more computed output variable, get set up
    SIDX = 0
    EIDX = 0

!.........  Reset FIRSTIME flag for either day-spec or hour-spec
    IF( .NOT. FIRSTIME .AND. .NOT. DAYFLAG ) FIRSTIME = .TRUE.

!.........  For the first time the routine is called...
    IF( FIRSTIME .AND. .NOT. HOURFLAG ) THEN

!.............  Get settings from the environment.
        DFLAG = ENVYN( 'RAW_DUP_CHECK',&
        &               'Error if duplicate inventory records',&
        &               .FALSE., IOS )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "RAW_DUP_CHECK"', 2 )
        END IF

!.............  Get maximum number of warnings
        MXWARN = ENVINT( WARNSET , ' ', 100, I )
        IF ( IOS .GT. 0 ) THEN
            CALL M3EXIT( PROGNAME,0,0, 'Bad env vble WARNSET', 2 )
        END IF

!.............  Allocate memory for flag for writing missing-data messages
        IF( ALLOCATED( NOMISS ) ) DEALLOCATE( NOMISS )
        ALLOCATE( NOMISS( NSRC,NVASP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'NOMISS', PROGNAME )
        NOMISS = .TRUE.  ! Array

!.............  Create reverse index for pollutants and activities
        IF (ALLOCATED (EAIDX2 )) DEALLOCATE (EAIDX2)
        MXEA = MAXVAL( EAIDX )
        ALLOCATE( EAIDX2( MXEA ), STAT=IOS )
        CALL CHECKMEM( IOS, 'EAIDX2', PROGNAME )
        EAIDX2 = 0

        NNHV = 0
        DO V = 1, NVAR
            POLNAM = EANAM( EAIDX( V ) )
            L = INDEX( POLNAM, 'NONHAP' )
            IF( L > 0 ) NNHV = NNHV + 1
            EAIDX2( EAIDX( V ) ) = V
        END DO

!.............  Create array to store original no of VOC/TOG values
        IF( INTGRFLAG ) THEN
            IF( ALLOCATED( NHAPMOD ) ) DEALLOCATE( NHAPMOD )
            ALLOCATE( NHAPMOD( NNHV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NHAPMOD', PROGNAME )
            IF( ALLOCATED( NHAPPOS ) ) DEALLOCATE( NHAPPOS )
            ALLOCATE( NHAPPOS( NNHV ), STAT=IOS )
            CALL CHECKMEM( IOS, 'NHAPPOS', PROGNAME )
            NHAPMOD = ''
            NHAPPOS = 0
        END IF

!..............  Build mode array
        NNHV = 0
        DO V = 1, NVAR
            POLNAM = EANAM( EAIDX( V ) )
            L = INDEX( POLNAM, 'NONHAP' )
            IF( L < 1 ) CYCLE

            L  = INDEX( POLNAM, ETJOIN )
            LL = LEN_TRIM ( POLNAM )
            IF( L > 0  ) THEN
                TMPMOD = POLNAM( 1:L-1 )
                INVNAM = POLNAM( L+2:LL )
            ELSE
                TMPMOD = 'TMP'
                INVNAM = POLNAM
            END IF

            IF( INDEX1( TMPMOD, NNHV, NHAPMOD ) < 1 ) THEN
                NNHV = NNHV + 1
                NHAPMOD( NNHV ) = TMPMOD
                NHAPPOS( NNHV ) = V
            END IF
        END DO

!.............  Create reverse index for special variables
        IF( ALLOCATED (SPIDX2)) DEALLOCATE (SPIDX2)
        ALLOCATE( SPIDX2( MXSPDAT ), STAT=IOS )
        CALL CHECKMEM( IOS, 'SPIDX2', PROGNAME )
        SPIDX2 = 0

        DO V = 1, MXSPDAT
            K = SPIDX( V )
            IF( K .GT. 0 ) SPIDX2( K ) = V
        END DO

        FIRSTIME = .FALSE.
        IF( .NOT. DAYFLAG ) HOURFLAG = .TRUE.

    END IF

    PDIDX  = 0        ! array (index)
    PDDATA = BADVAL3  ! array (emissions/activities)
    PDTOTL = BADVAL3  ! array (total daily emissions/activities)

!.........  Sort sources & inventory data codes for current time step
!.........  Added CIDXA into sorting in SMOKE 2.3.2 to handle now using
!           inventory data codes instead of SMOKE pollutant codes only.
!           Want to prevent false duplicates caused by Inventory Table renaming
!           multiple Inventory Data Code (i.e., CAS numbers) to the same
!           SMOKE name.  These are not really duplicates and should be
!           ignored by the warning messages below.

    CALL SORTI2( NPDPT( TIDX ), IDXSRC( 1,TIDX ),&
    &             SPDIDA( 1,TIDX ), CIDXA( 1,TIDX ) )

!.........  Store sorted records for this hour
    LS = 0  ! previous source
    LN = -9 ! previous Inventory Data Code position in UNIQCAS
    K  = 0
    LK = 0
    DO I = 1, MIN( NPDPT( TIDX ), MXPDSRC )

        J = IDXSRC( I,TIDX )
        S = SPDIDA( J,TIDX )
        V = CODEA ( J,TIDX )
        N = CIDXA ( J,TIDX )

        IF( S < 1 ) CYCLE  ! Skip if source is missing

        POLNAM = EANAM( V )

!...........  Add multiple inventory pollutant(s) with same CAS name
!             Find code corresponding to current pollutant before you
        NPPCAS = UCASNKEP( N )

        IF( CFLAG ) THEN
            NPPCAS = 1
            CV     = V
            NCOMP = 0
            INTGRFLAG = .FALSE.
        END IF

        DO NP = 0, NPPCAS - 1

            NC = UCASIDX( N ) + NP
            IF( NP > 0 ) POLNAM = ITNAMA( SCASIDX( NC ) )
            V = INDEX1( POLNAM, NIPPA, EANAM )
            IF( NP > 0 ) LN = 0
            IF( CFLAG  ) V = CV                  ! restore original poll idx
            IF( V < 1  ) CYCLE                   ! skip if it is not listed in ann inv poll

!.............  Intialize as not a special data variable (e.g., not flow rate)
            LFLAG = .FALSE.

!.............  Check for index for special variables and set V to be consistent
!               with the output structure of the file
            IF ( V .GT. CODFLAG3 ) THEN   ! CODFLAG3 = 9000 for special data types
                V2 = V - CODFLAG3      ! remove flag on index
                V = NVAR + SPIDX( V2 ) ! reset to condensed order from master
                LFLAG = .TRUE.         ! flag as a special data variable

!.............  Otherwise, set index for period-specific pollutant or activity
            ELSE
                V = EAIDX2( V )

            END IF

!.............  Initialize duplicates flag
            DUPFLAG = .FALSE.
            NHVPOS = 0

!.............  If current source is not equal to previous source
            IF( S .NE. LS  .OR. I .EQ. NPDPT( TIDX ) ) THEN

!.................  Count no of source
                IF( S .NE. LS ) THEN
                    K = K + 1
                    PDIDX( K ) = S
                END IF

!.................  Get the location of start and end index
                IF ( K .EQ. 1 ) THEN
                    SIDX( K ) = I
                ELSE IF ( K .GT. 1 .AND. I .LT. NPDPT( TIDX )) THEN
                    SIDX( K ) = I
                    EIDX( K-1 ) = I - 1
                    LK = K-1
!.................  Only specify the end index for last source
                ELSE IF ( I .EQ.  NPDPT( TIDX )) THEN
                    LK = K
                    EIDX( K ) = I
                END IF

!.................  Calculate formula if needed
                IF( NP < 1 ) THEN       ! Skip PMC and Integration calcuation
                    IF( NCOMP > 0 .AND. K > 1 ) THEN

                        DO F = 1, NCOMP
                            IDXA = 0
                            IDXB = 0
                            NV = INDEX1( VNAME(F), NIPPA, EANAM )
                            NVPOS = EAIDX2( NV )     ! output pollutant index

                            FND_VINA = .FALSE.
                            FND_VINB = .FALSE.

                            DO II = SIDX(LK), EIDX(LK)
                                JJ = IDXSRC( II,TIDX )
                                VV = CODEA ( JJ,TIDX )
                                IF ( EANAM(VV) == VIN_A(F) ) THEN
                                    IDXA = JJ
                                    FND_VINA = .TRUE.
                                ELSE IF ( EANAM(VV) == VIN_B(F) ) THEN
                                    IDXB = JJ
                                    FND_VINB = .TRUE.
                                END IF
                            END DO    ! end of variable search loop

!.........................  Skip if formula variables are missing
                            IF( .NOT. FND_VINA .OR. .NOT. FND_VINB ) CYCLE

!.........................  Calculate formula variables
                            IF( IDXA > 0 .AND. IDXB > 0 .AND. NVPOS > 0 ) THEN
                                IF( CHKPLUS(F) )  THEN
                                    IF( PDDATA( LK,NVPOS ) > 0.0 ) THEN
                                        PDDATA( LK,NVPOS ) =  PDDATA( LK,NVPOS ) +&
                                        &    EMISVA( IDXA, TIDX ) + EMISVA( IDXB, TIDX )
                                        PDTOTL( LK,NVPOS ) =  PDTOTL( LK,NVPOS ) +&
                                        &    DYTOTA( IDXA, TIDX ) + DYTOTA( IDXB, TIDX )
                                    ELSE
                                        PDDATA( LK,NVPOS ) =&
                                        &    EMISVA( IDXA, TIDX ) + EMISVA( IDXB, TIDX )
                                        PDTOTL( LK,NVPOS ) =&
                                        &    DYTOTA( IDXA, TIDX ) + DYTOTA( IDXB, TIDX )
                                    END IF

                                END IF

                                IF( CHKMINUS(F) )  THEN
                                    IF( PDDATA( LK,NVPOS ) > 0.0 ) THEN
                                        PDDATA( LK,NVPOS ) = PDDATA( LK,NVPOS )+&
                                        &    EMISVA( IDXA, TIDX ) - EMISVA( IDXB, TIDX )
                                        PDTOTL( LK,NVPOS ) = PDTOTL( LK,NVPOS )+&
                                        &    DYTOTA( IDXA, TIDX ) - DYTOTA( IDXB, TIDX )
                                    ELSE
                                        PDDATA( LK,NVPOS ) =&
                                        &    EMISVA( IDXA, TIDX ) - EMISVA( IDXB, TIDX )
                                        PDTOTL( LK,NVPOS ) =&
                                        &    DYTOTA( IDXA, TIDX ) - DYTOTA( IDXB, TIDX )
                                    END IF
                                END IF

!.............................  Check for negative values for daily value
                                IF( PDDATA( LK,NVPOS ) < 0.0 )  THEN
                                    WARNCNT = WARNCNT + 1
                                    IF ( WARNCNT .LE. MXWARN ) THEN
                                        CALL FMTCSRC( CSOURC( LS ), NCHARS, BUFFER, L2 )
                                        TYPE = 'Hourly'
                                        IF( DAYFLAG ) TYPE = 'Daily'
                                        WRITE( MESG,94020 ) 'WARNING: '//&
                                        & 'Resetting negative value of "'//&
                                        & TRIM(TYPE)//' "'//TRIM(VNAME(F))//&
                                        & '" from',PDDATA( LK,NVPOS ),'to 0. for '//&
                                        &'source:'//CRLF()//BLANK10//BUFFER(1:L2)
                                        MXWARN = MXWARN + 1
                                    END IF
!                                CALL M3MESG( MESG )
                                    PDDATA( LK,NVPOS ) = 0.0
                                    PDTOTL( LK,NVPOS ) = 0.0
                                END IF

                            ELSE IF( WARNCNT .LE. MXWARN ) THEN
                                CALL FMTCSRC( CSOURC( LS ), NCHARS, BUFFER, L2 )
                                EFLAG = .TRUE.
                                MESG = 'ERROR: no sources for calculating ' //&
                                &   'variable "'//TRIM(VNAME(F))//'" J '//&
                                &    CRLF() // BLANK10 // BUFFER( 1:L2 )
                                CALL M3MSG2( MESG )
                                WARNCNT = WARNCNT + 1
                                CYCLE

                            END IF

                        END DO    ! end of formula loop

                    END IF   ! end of formula calculation

!.................  Combine VOC + HAPs for integrate/non-integrate option
!.................  Determine no of VOCs and store VOC/NONHAPVOC positions
!                   for later NONHAPVOC calculation by substracting
                    IF( INTGRFLAG .AND. K > 1 ) THEN

!.....................  Count no of HAPs
                        NHAP = 0
                        NVOC = 0
                        DO II = SIDX(LK), EIDX(LK)
                            JJ = IDXSRC( II,TIDX )
                            VV = CODEA ( JJ,TIDX )
                            IV = EAIDX2( VV )

                            POLNAM = EANAM( VV )

!.........................  Store values when NONHAP[VOC|TOG] is found
                            IF( INDEX( POLNAM, 'NONHAP' ) > 0 ) THEN
                                NVOC = NVOC + 1
                                PDDATA( LK,IV ) = EMISVA( JJ,TIDX )
                                PDTOTL( LK,IV ) = DYTOTA( JJ,TIDX )
                                IF( I == NPDPT( TIDX ) ) THEN
                                    PDDATA( LK,IV ) = 0.0
                                    PDTOTL( LK,IV ) = 0.0
                                END IF
                            END IF

!.........................  Skip if it is activity data, skip
!.........................  Skip if pollutant is not part of VOC or TOG, cycle
                            NV = INDEX1( POLNAM, MXIDAT, INVDNAM )
                            IF( INVSTAT( NV ) <   0  ) CYCLE
                            IF( INVDVTS( NV ) == 'N' ) CYCLE

                            NHAP = NHAP + 1

                        END DO

!.....................  Process non-integrated sources
                        IF( CINTGR( LS ) == 'N' ) THEN

!.........................  Search for HAP to converted to HAP_NOI for non-integrated source
                            DO II = SIDX(LK), EIDX(LK)
                                JJ = IDXSRC( II,TIDX )
                                VV = CODEA ( JJ,TIDX )
                                IV = EAIDX2( VV )

                                POLNAM = EANAM( VV )

!.................................  If pollutant is not a model species, set it to zero
                                IF( INDEX( POLNAM, '_NOI' ) > 0 ) THEN
                                    INVNAM = POLNAM
                                    L  = INDEX( POLNAM,'_NOI' )
                                    IF( L > 0 ) INVNAM = POLNAM( 1:L-1 )
                                    NVRAW = INDEX1( INVNAM, NINVTBL, ITNAMA )
                                    IF( .NOT. ITMSPC( NVRAW ) ) EMISVA( JJ, TIDX ) = 0.0
                                    PDDATA( LK,IV ) = EMISVA( JJ,TIDX )
                                    PDTOTL( LK,IV ) = DYTOTA( JJ,TIDX )
                                END IF

                            END DO

!.....................  Process integrated sources : NONHAPVOC = VOC - all HAPs
                        ELSE IF( CINTGR( LS ) == 'Y' ) THEN

                            IF( NVOC > 0 .AND. NHAP > 0 ) THEN
                                DO II = SIDX(LK), EIDX(LK)
                                    JJ = IDXSRC( II,TIDX )
                                    VV = CODEA ( JJ,TIDX )

                                    POLNAM = EANAM( VV )
                                    NV = INDEX1( POLNAM, MXIDAT, INVDNAM )

!.........................  Skip if it is activity data, skip
                                    IF( INVSTAT( NV ) < 0 ) CYCLE

!.........................  Skip if pollutant is not part of VOC or TOG, cycle
                                    IF( INVDVTS( NV ) == 'N' ) CYCLE

                                    L  = INDEX( POLNAM, ETJOIN )
                                    IF( L > 0  ) THEN
                                        TMPMOD = POLNAM( 1:L-1 )
                                    ELSE
                                        TMPMOD = 'TMP'
                                    END IF

!.............................  Retrieve NONHAPVOC positions
                                    M = INDEX1( TMPMOD, NNHV, NHAPMOD )
                                    NHVPOS = NHAPPOS( M )

                                    IF( NHVPOS < 1 ) CYCLE

!.............................  Compute NONHAPVOCs
!                               Subtract all mode-specific HAPs from Original VOC
                                    PDDATA( LK,NHVPOS ) = PDDATA( LK,NHVPOS )&
                                    &                      - EMISVA( JJ,TIDX )
                                    PDTOTL( LK,NHVPOS ) = PDTOTL( LK,NHVPOS )&
                                    &                      - DYTOTA( JJ,TIDX )

                                END DO    ! end of variable search loop

!.....................  Error if VOC or HAP is missing for integration
                            ELSE IF( NVOC > 0 .AND. NHAP < 1 ) THEN
                                CALL FMTCSRC( CSOURC( LS ), NCHARS, BUFFER, L2 )
                                MESG = 'ERROR: Found VOC|TOG but no toxics found '//&
                                &     ' for the source:'//CRLF()//BLANK10//BUFFER(1:L2)
                                CALL M3MESG( MESG )
                                EFLAG = .TRUE.

                            ELSE IF( NVOC < 1 .AND. NHAP > 0 ) THEN
                                CALL FMTCSRC( CSOURC( LS ), NCHARS, BUFFER, L2 )
                                MESG = 'ERROR: Found toxics but no VOC|TOG found '//&
                                &     ' for the source:'//CRLF()//BLANK10//BUFFER(1:L2)
                                CALL M3MESG( MESG )
                                EFLAG = .TRUE.

                            ELSE IF( NVOC < 1 .AND. NHAP < 1 ) THEN
                                CALL FMTCSRC( CSOURC( LS ), NCHARS, BUFFER, L2 )
                                MESG = 'ERROR: Both VOC|TOG and toxics are not found '//&
                                &     ' for the source:'//CRLF()//BLANK10//BUFFER(1:L2)
                                CALL M3MESG( MESG )
                                EFLAG = .TRUE.
                            END IF

                        END IF   ! integrated sources only

                    END IF
                END IF   ! NP loop

!.................  set previous source id
                IF( S .NE. LS ) LS = S

!............. If source is the same, look for duplicate data variables
            ELSE IF ( N .EQ. LN ) THEN
                DUPFLAG = .TRUE.      !  This iteration has a duplicate

            END IF

!.............  If emissions are not yet set for current source and variable
            IF( PDDATA( K,V ) .LT. AMISS3 ) THEN  ! PDDATA is blank

                PDDATA( K,V ) = EMISVA( J,TIDX )
                PDTOTL( K,V ) = DYTOTA( J,TIDX )

!.............  Otherwise, sum emissions, report, and set error if needed
            ELSE

                IF( .NOT. LFLAG .AND. EMISVA( J,TIDX ) .GE. 0. )&
                &    PDDATA( K,V )  = PDDATA( K,V ) + EMISVA( J,TIDX )

                IF( .NOT. LFLAG .AND. DYTOTA( J,TIDX ) .GE. 0. )&
                &    PDTOTL( K,V )  = PDTOTL( K,V ) + DYTOTA( J,TIDX )

!.................  If the source is an actual duplicate (the same source
!                   and Inventory Data Name), proceed accordingly.  Do
!                   not write warnings for Inventory Data Names that
!                   are combined into the same SMOKE name.
                IF ( DUPFLAG ) THEN
                    CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )
                    NWARN1 = NWARN1 + 1

                    IF( NWARN1 <= MXWARN ) THEN
                        IF( DFLAG ) THEN
                            EFLAG = .TRUE.
                            MESG = 'ERROR: Duplicate source in inventory:'//&
                            &       CRLF() // BLANK10 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                            CYCLE
                        ELSE IF ( LFLAG ) THEN
                            MESG = 'WARNING: Duplicate source in ' //&
                            &       'inventory. Will store only one value '//&
                            &       'for '// CRLF()// BLANK10//SPDATDSC(V2)//&
                            &       ':' // CRLF() // BLANK10 // BUFFER(1:L2)
                            CALL M3MESG( MESG )
                        ELSE
                            MESG = 'WARNING: Duplicate source in ' //&
                            &       'inventory will have summed emissions:'&
                            &       //CRLF() // BLANK10 // BUFFER( 1:L2 )
                            CALL M3MESG( MESG )
                        END IF
                    END IF

                END IF   !  If duplicate or not

            END IF

            LN = N  ! Set LN for next iteration

        END DO    ! loop over no of poll share same CAS

    END DO  ! End loop over data for this time step

!.........  Set tmp variable for loops
    NOUT = K

!.........  Check if there are missing values and output errors, if the
!           flag is set to treat these as errors
!.........  Also use loop to create diurnal profiles from emission values,
!           if needed.
    DO V = 1, NVASP

        DO I = 1, NOUT

            S = PDIDX( I )

!.................  Format source information
            CALL FMTCSRC( CSOURC( S ), NCHARS, BUFFER, L2 )

!.................  Check for missing values
            IF ( PDDATA( I,V ) .LT. AMISS3 ) THEN

!.....................  For fires data, reset "missing" values
!                       to zero, since Temporal will not need to
!                       know the difference between missing and zero.
!                       This is because there are no "annual" values
!                       to fill in for missing daily/hourly data for fires.
                IF( FIREFLAG ) PDDATA( I,V ) = 0.

                IF( NOMISS( S,V ) ) THEN
                    NWARN = NWARN + 1

                    IF( NWARN <= MXWARN ) THEN
                        IF ( V .LE. NVAR ) THEN
                            MESG = 'WARNING: Data missing for: ' //&
                            &   CRLF()//BLANK10//BUFFER( 1:L2 )//&
                            &   ' VAR: '// EANAM( EAIDX( V ) )
                        ELSE
                            K = V - NVAR
                            MESG = 'WARNING: Data missing for: ' //&
                            &   CRLF()//BLANK10//BUFFER( 1:L2 )//&
                            &   ' VAR: '//SPDATNAM( SPIDX2( K ) )
                        END IF

                        CALL M3MESG( MESG )
                    END IF
                    NOMISS( S,V ) = .FALSE.
                END IF

            END IF

!.................  Skip next section if special data variable
            IF ( V .GT. NVAR ) CYCLE

!.................  If profiles need to be created instead of hourly emissions
!.................  Make sure totals data are available!
            IF ( PFLAG .AND. PDTOTL( I,V ) .LT. AMISS3 ) THEN
                EFLAG = .TRUE.
                MESG = 'ERROR: Cannot create profiles ' //&
                &       'because daily total data missing for:' //&
                &       CRLF() // BLANK10 // BUFFER( 1:L2 ) //&
                &       ' VAR: ' // EANAM( EAIDX( V ) )
                CALL M3MESG( MESG )
                CYCLE

!.................  Prevent divide by zero
!.................  If totals are zero, then profile values should be zero
            ELSE IF( PFLAG .AND. PDTOTL( I,V ) .GT. 0 ) THEN

                PDDATA( I,V ) = PDDATA( I,V ) / PDTOTL( I,V )

            END IF

        END DO
    END DO

!.........  Write emissions for this time step
    IF ( .NOT. WRITE3( FNAME, ALLVAR3, JDATE, JTIME, PDIDX ) ) THEN
        MESG= 'Error writing output file "' // TRIM( FNAME ) // '"'
        CALL M3EXIT( PROGNAME, JDATE, JTIME, MESG, 2 )

    END IF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************
!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )
94020 FORMAT( 10( A, :, E10.2, :, 1X ) )

END SUBROUTINE WRPDEMIS
