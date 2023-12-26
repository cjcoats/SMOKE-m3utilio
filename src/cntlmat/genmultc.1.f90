
SUBROUTINE GENMULTC( CDEV, GDEV, LDEV, MDEV, NCPE, PYEAR,&
&                     ENAME, MNAME, CFLAG, GFLAG, LFLAG,&
&                     MFLAG )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!      This subroutine computes the multiplicative control factors
!      and writes out the multiplicative control matrix.
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Version ??/???? by ???
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
!***************************************************************************
!
! Project Title: Sparse Matrix Operator Kernel gsions (SMOKE) Modeling
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
!.........  This module contains the inventory arrays
    USE MODSOURC, ONLY: CSOURC

!.........  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: FACCEFF, FACREFF, FACRLPN,&
    &                    NVCMULT, PNAMMULT, FACTOR,&
    &                    DATVAL, BACKOUT, GRPINDX, GRPFLAG,&
    &                    GRPINEM, GRPOUTEM, RPTDEV, PCTLFLAG,&
    &                    CUTCTG, FACCTG, FACMACT,&
    &                    FACRACT, EMCAPALW, EMREPALW, GRPSTIDX,&
    &                    GRPCHAR, MACEXEFF, MACNWEFF,&
    &                    MACNWFRC, CTLRPLC, CTGCOMT, CTLCOMT, ALWCOMT,&
    &                    REACOMT, PRJCOMT, MACCOMT

!.........  This module contains the information about the source category
    USE MODINFO, ONLY: CATEGORY, NSRC, NEM, NDY, NCE, NRE, NRP,&
    &                   NPPOL, BYEAR, NCHARS, CATDESC, CRL

    IMPLICIT NONE

!...........   INCLUDES

    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters
    INCLUDE 'SETDECL.EXT'   !  FileSetAPI variables and functions

!...........   EXTERNAL FUNCTIONS and their descriptions:
    LOGICAL, EXTERNAL :: FLTERR

!...........   SUBROUTINE ARGUMENTS

    INTEGER,       INTENT (IN OUT) :: CDEV   ! file unit no. for tmp CTL file
    INTEGER,       INTENT (IN OUT) :: GDEV   ! file unit no. for tmp CTG file
    INTEGER,       INTENT (IN OUT) :: LDEV   ! file unit no. for tmp ALW file
    INTEGER,       INTENT (IN OUT) :: MDEV   ! file unit no. for tmp MACT file
    INTEGER,       INTENT (IN) :: NCPE   ! no. of control packet entries
    INTEGER,       INTENT (IN) :: PYEAR  ! projected year, or missing
    CHARACTER(16), INTENT (IN) :: ENAME  ! logical name for i/o api
    ! inventory input file
    CHARACTER(16), INTENT (IN OUT) :: MNAME  ! logical name for mult. cntl. mat.
    LOGICAL,       INTENT (IN) :: CFLAG  ! true = apply CTL controls
    LOGICAL,       INTENT (IN) :: GFLAG  ! true = apply CTG controls
    LOGICAL,       INTENT (IN) :: LFLAG  ! true = apply ALW controls
    LOGICAL,       INTENT (IN) :: MFLAG  ! true = apply MACT controls

!...........   Local allocatable arrays
!        INTEGER, ALLOCATABLE :: ALWINDX ( :,: ) ! indices to ALW controls table
!        INTEGER, ALLOCATABLE :: CTGINDX ( :,: ) ! indices to CTG controls table
!        INTEGER, ALLOCATABLE :: CTLINDX ( :,: ) ! indices to CTL controls table
!        INTEGER, ALLOCATABLE :: GRPINDX ( : )   ! index from sources to groups
!        INTEGER, ALLOCATABLE :: GRPSTIDX( : )   ! sorting index

!        REAL   , ALLOCATABLE :: GRPINEM ( :,: ) ! initial emissions
!        REAL   , ALLOCATABLE :: GRPOUTEM( :,: ) ! controlled emissions

!        LOGICAL, ALLOCATABLE :: GRPFLAG ( : )   ! true: group controlled

!        CHARACTER(STALEN3+SCCLEN3), ALLOCATABLE :: GRPCHAR( : ) ! group chars
!        REAL   , ALLOCATABLE :: BACKOUT ( : )   ! factor used to account for pol
    ! specific control info that is
    ! already in the inventory
!        REAL   , ALLOCATABLE :: DATVAL  ( :,: ) ! emissions and control settings
!        REAL   , ALLOCATABLE :: FACTOR  ( : )   ! multiplicative controls

!.........   Local arrays
    INTEGER                 OUTTYPES( NVCMULT,6 ) ! var type:int/real
    INTEGER                 ODEV( 4 )             ! tmp output files
    CHARACTER(NAMLEN3)      OUTNAMES( NVCMULT,6 ) ! var names
    CHARACTER(IOULEN3)      OUTUNITS( NVCMULT,6 ) ! var units
    CHARACTER(MXDLEN3)      OUTDESCS( NVCMULT,6 ) ! var descriptions

!...........   Other local variables
    INTEGER          E, I, J, K, L, L2, S  ! counters and indices

    INTEGER       :: ALWINDX = 0 ! indices to ALW controls table
    INTEGER       :: CTGINDX = 0 ! indices to CTG controls table
    INTEGER       :: CTLINDX = 0 ! indices to CTL controls table
    INTEGER       :: MACINDX = 0 ! indices to MACT controls table
    INTEGER          IDX      ! group index
    INTEGER          IOS      ! input/output status
    INTEGER          PIDX     ! previous IDX
    INTEGER          RDEV     ! Report unit number
    INTEGER          RDEV2    ! Report unit number 2
    INTEGER          SCCBEG   ! begining of SCC in CSOURC string
    INTEGER          SCCEND   ! end of SCC in CSOURC string

    REAL             ALWEMIS  ! allowable emissions
    REAL             CAP      ! emissions cap
    REAL             CTGFAC   ! control technology control factor
    REAL             CTGFAC2  ! MAXACF or RSNACF
    REAL             CTLEFF   ! tmp control efficiency
    REAL             CUTOFF   ! CTG cutoff for application of control
    REAL             DENOM    ! denominator of control back-out factor
    REAL             E1, E2   ! tmp emissions values
    REAL             EXSTCEFF ! control efficiency for existing sources
    REAL             FAC      ! tmp control factor
    REAL             INVEFF   ! inventory baseline efficiency
    REAL             MACT     ! max. achievable cntrl tech. cntrl factor
    REAL             NEWCEFF  ! control efficiency for new sources
    REAL             NEWFAC   ! tmp for new 1-ceff*reff*rpen
    REAL             NEWFRAC  ! fraction of new sources vs. existing
    REAL             OLDFAC   ! tmp for old (base) 1-ceff*reff*rpen
    REAL             RACT     ! reasonably achiev. cntrl tech. cntrl factor
    REAL             REPLACE  ! replacement emissions
    REAL             RULEFF   ! tmp rule effectiveness
    REAL             RULPEN   ! tmp rule penetration

    REAL, ALLOCATABLE :: BASECEFF( : )  ! control eff. from base inv.
    REAL, ALLOCATABLE :: BASEREFF( : )  ! rule eff. from base inv.
    REAL, ALLOCATABLE :: BASERPEN( : )  ! rule penetration from base inv.
    REAL, ALLOCATABLE :: CEFF( : )   ! control efficiency per source
    REAL, ALLOCATABLE :: REFF( : )   ! rule effectiveness per source
    REAL, ALLOCATABLE :: RPEN( : )   ! rule penetration per source

    LOGICAL          LAVEDAY  ! true: use average day emissions
    LOGICAL          LCOMPARE ! true: during a "replacement" of controls,
!                                         make sure new controls are more reduction
!                                         than the old controls.
    LOGICAL          LREPORT  ! true: output a record to report
    LOGICAL       :: LREPHDR  = .TRUE.  ! true: output report header
    LOGICAL, SAVE :: APPLFLAG = .FALSE. ! true: something has been applied
    LOGICAL, SAVE :: OPENFLAG = .FALSE. ! true: output file has been opened

    CHARACTER(1)       REPCTYPE   ! control type for report
    CHARACTER(100)     OUTFMT     ! header format buffer
    CHARACTER(200)     PATHNM     ! path name for tmp file
    CHARACTER(220)     FILENM     ! file name
    CHARACTER(256)     BUFFER     ! source fields buffer
    CHARACTER(256)     MESG       ! message buffer
    CHARACTER(SRCLEN3) CSRC       ! tmp source chars
    CHARACTER(NAMLEN3) PNAM       ! tmp pollutant name
    CHARACTER(NAMLEN3) CBUF       ! pollutant name temporary buffer
    CHARACTER(NAMLEN3) EBUF       ! pollutant name temporary buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'GENMULTC' ! program name

!***********************************************************************
!   begin body of subroutine GENMULTC

!.........  Get environment variables that control program behavior
    MESG = 'Use annual or average day emissions'
    LAVEDAY = ENVYN( 'SMK_AVEDAY_YN', MESG, .FALSE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_AVEDAY_YN"', 2 )
    END IF

!.........  Get path for temporary files
    MESG = 'Path where temporary control files will be written'
    CALL ENVSTR( 'SMK_TMPDIR', MESG, '.', PATHNM, IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "SMK_TMPDIR"', 2 )
    END IF

!.........  Get environment variables that control program behavior
    MESG = 'Use annual or average day emissions'
    LCOMPARE = ENVYN( 'COMPARE_REPLACE_CONTROL', MESG, .TRUE., IOS )
    IF ( IOS .GT. 0 ) THEN
        CALL M3EXIT( PROGNAME,0,0, 'Bad env vble "COMPARE_REPLACE_CONTROL"', 2 )
    END IF

!.........  Open reports files
    RPTDEV( 1 ) = PROMPTFFILE(&
    &              'Enter logical name for MULTIPLICATIVE ' //&
    &              'CONTROLS REPORT',&
    &              .FALSE., .TRUE., CRL // 'CREP', PROGNAME )
    RDEV = RPTDEV( 1 )

    RPTDEV( 5 ) = PROMPTFFILE(&
    &              'Enter logical name for DETAILED ' //&
    &              'CONTROL PACKET REPORT',&
    &              .FALSE., .TRUE., CRL // 'CREP2', PROGNAME )
    RDEV2 = RPTDEV( 5 )

!.........  Open *output* temporary files depending on whether an input
!           temporary file exists - indicating that the packet is being used
    IF( CDEV .GT. 0 ) THEN
        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_ctl_rep'
        ODEV( 1 ) = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
    END IF
    IF( GDEV .GT. 0 ) THEN
        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_ctg_rep'
        ODEV( 2 ) = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
    END IF
    IF( LDEV .GT. 0 ) THEN
        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_alw_rep'
        ODEV( 3 ) = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
    ENDIF
    IF( MDEV .GT. 0 ) THEN
        FILENM = TRIM( PATHNM ) // '/cntlmat_tmp_mact_rep'
        ODEV( 4 ) = GETEFILE( FILENM, .FALSE., .TRUE., PROGNAME )
    END IF

!.........  Allocate index to reporting groups
!        ALLOCATE( GRPINDX( NSRC ), STAT=IOS )
!        CALL CHECKMEM( IOS, 'GRPINDX', PROGNAME )
!        ALLOCATE( GRPSTIDX( NSRC ), STAT=IOS )
!        CALL CHECKMEM( IOS, 'GRPSTIDX', PROGNAME )
!        ALLOCATE( GRPCHAR( NSRC ), STAT=IOS )
!        CALL CHECKMEM( IOS, 'GRPCHAR', PROGNAME )
!        GRPINDX  = 0  ! array

!.........  Get set up for group reporting...
!        IF( CATEGORY .EQ. 'POINT' ) THEN

!.............  Count the number of groups in the inventory
!           PPLT = ' '
!            NGRP = 0
!            DO S = 1, NSRC
!                CPLT = CSOURC( S )( 1:FPLLEN3 )
!                IF( CPLT .NE. PPLT ) THEN
!                    NGRP = NGRP + 1
!                    PPLT = CPLT
!                END IF
!                GRPINDX ( S ) = NGRP
!                GRPSTIDX( S ) = S     ! Needed for loops, but not used to sort
!            END DO

!        ELSE

    IF( CATEGORY .EQ. 'AREA' ) THEN
        SCCBEG = ARBEGL3( 2 )
        SCCEND = ARENDL3( 2 )
    ELSE            ! MOBILE
        SCCBEG = MBBEGL3( 5 )
        SCCEND = MBENDL3( 5 )
    END IF

!.............  Build and sort source array for SCC-state grouping
!            DO S = 1, NSRC
!                CSTA = CSOURC( S )( 1     :STALEN3 )
!                TSCC = CSOURC( S )( SCCBEG:SCCEND  )

!                GRPSTIDX( S ) = S
!                GRPCHAR ( S ) = CSTA // TSCC
!            END DO

!            CALL SORTIC( NSRC, GRPSTIDX, GRPCHAR )

!.............  Count the number of state/SCCs in the domain
!            PSTA = ' '
!            PSCC = ' '
!            SCCBEG = STALEN3 + 1
!            SCCEND = STALEN3 + SCCLEN3
!            DO S = 1, NSRC
!                J = GRPSTIDX( S )
!                CSTA = GRPCHAR( J )( 1     :STALEN3 )
!                TSCC = GRPCHAR( J )( SCCBEG:SCCEND  )
!                IF( CSTA .NE. PSTA .OR. TSCC .NE. PSCC ) THEN
!                    NGRP = NGRP + 1
!                    PSTA = CSTA
!                    PSCC = TSCC
!                END IF
!                GRPINDX( J ) = NGRP
!            END DO

!        END IF

!...........  Allocate memory for the number of groups for storing emissions
!          ALLOCATE( GRPFLAG( NGRP ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'GRPFLAG', PROGNAME )
!          ALLOCATE( GRPINEM( NGRP, NVCMULT ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'GRPINEM', PROGNAME )
!          ALLOCATE( GRPOUTEM( NGRP, NVCMULT ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'GRPOUTEM', PROGNAME )

!...........  Initialize
!          GRPINEM  = 0. ! array
!          GRPOUTEM = 0. ! array
!          GRPFLAG  = .FALSE.  ! array

!...........  Allocate local memory
!          ALLOCATE( BACKOUT( NSRC ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'BACKOUT', PROGNAME )
!          ALLOCATE( DATVAL( NSRC,NPPOL ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'DATVAL', PROGNAME )
!          ALLOCATE( FACTOR( NSRC ), STAT=IOS )
!          CALL CHECKMEM( IOS, 'FACTOR', PROGNAME )

!...........  For each pollutant that receives controls, obtain variable
!             names for control efficiency, rule effectiveness, and, in the
!             case of AREA sources, rule penetration. These variable names
!             will be used in reading the inventory file.

    CALL BLDENAMS( CATEGORY, NVCMULT, 6, PNAMMULT, OUTNAMES,&
    &               OUTUNITS, OUTTYPES, OUTDESCS )

    DO I = 1, 6
        IF( OUTNAMES( 1,I )(1:NAMLEN3)  .EQ. PNAMMULT(1) ) NEM = I
        IF( OUTNAMES( 1,I )(1:CPRTLEN3) .EQ. AVEDAYRT ) NDY = I
        IF( OUTNAMES( 1,I )(1:CPRTLEN3) .EQ. CTLEFFRT ) NCE = I
        IF( OUTNAMES( 1,I )(1:CPRTLEN3) .EQ. RULEFFRT ) NRE = I
        IF( OUTNAMES( 1,I )(1:CPRTLEN3) .EQ. RULPENRT ) NRP = I
    END DO

!...........  Ensure the temporary files, if opened, are rewound
    IF( CDEV .GT. 0 ) REWIND( CDEV )
    IF( GDEV .GT. 0 ) REWIND( GDEV )
    IF( LDEV .GT. 0 ) REWIND( LDEV )
    IF( MDEV .GT. 0 ) REWIND( MDEV )

!...........  Fractionalize control-packet information

    IF( MFLAG ) THEN
        MACEXEFF = MACEXEFF*0.01  ! array
        MACNWEFF = MACNWEFF*0.01  ! array
        MACNWFRC = MACNWFRC*0.01  ! array
    END IF

    IF( CFLAG ) THEN
        FACCEFF = FACCEFF*0.01  ! array
        FACREFF = FACREFF*0.01  ! array
        FACRLPN = FACRLPN*0.01  ! array
    END IF

!...........  Set emissions index depending on average day or not
    E = NEM
    IF( LAVEDAY ) E = NDY

!...........  Loop through pollutants that receive controls
    DO I = 1, NVCMULT

!............  Set tmp pollutant name
        PNAM = PNAMMULT( I )
        EBUF = OUTNAMES(I,1)          ! set annual emis name

!............  Initialize control factor array
        FACTOR = 1.0  ! array

!...........  Read in emissions data from inventory file...
!...........  From map-formatted inventory or old format
        CALL RDMAPPOL( NSRC, 1, NPPOL, EBUF, DATVAL )

!...........  Adjust emissions values and compute group values
        FAC = YR2DAY( BYEAR )  ! year to day factor for later
        DO S = 1, NSRC

!...............  Check for missing values and reset to zero
            IF( DATVAL( S,E ) .LT. AMISS3 ) DATVAL( S,E ) = 0.0

!...............  Divide annual emissions to get average day
            IF( .NOT. LAVEDAY ) DATVAL( S,E ) = DATVAL( S,E ) * FAC

!.................  Compute group emissions before controls
            J = GRPINDX( S )
            GRPINEM ( J,I ) = GRPINEM ( J,I ) + DATVAL( S,E )

        END DO ! end source loop

!...........  If CONTROL packet is present: For the current pollutant, read
!             in control efficiency, rule effectiveness, and, in the case of
!             AREA sources, rule penetration.
        IF ( CFLAG ) THEN

!...........  Allocate arrays to store control information from base inventory
            IF( ALLOCATED( BASECEFF ) ) THEN
                DEALLOCATE( BASECEFF, BASEREFF, BASERPEN )
            END IF
            ALLOCATE( BASECEFF( NSRC ) )
            ALLOCATE( BASEREFF( NSRC ) )
            ALLOCATE( BASERPEN( NSRC ) )
            BASECEFF = 0.
            BASEREFF = 1.
            BASERPEN = 1.

!...........  Then calculate the factor which will be used to account
!             for control information already in the inventory
            DO S = 1, NSRC

                CTLEFF = 0.
                RULEFF = 1.
                RULPEN = 1.
                IF ( NCE .GT. 0 ) CTLEFF = DATVAL( S,NCE )
                IF ( NRE .GT. 0 ) RULEFF = DATVAL( S,NRE )
                IF ( NRP .GT. 0 ) RULPEN = DATVAL( S,NRP )

!..................  Store control information from base inventory
                BASECEFF( S ) = CTLEFF
                BASEREFF( S ) = RULEFF
                BASERPEN( S ) = RULPEN

!..................  Perform division by zero check.

                DENOM = ( 1.0 - CTLEFF*RULEFF*RULPEN )

                IF ( FLTERR( DENOM, 0.0 ) ) THEN
                    BACKOUT( S ) = 1.0/DENOM
                ELSE
                    BACKOUT( S ) = 0.0
                END IF

            END DO ! end source loop

        END IF

!.............................................................................
!...........  Apply /MACT/ packet controls if present for the current
!             pollutant
!.............................................................................
        IF ( MFLAG .AND. PCTLFLAG( I, 4 ) ) THEN

!...............  Loop through sources
            DO S = 1, NSRC

                E1  = DATVAL( S, E )
                FAC = 1.

!................  Read temporary input file indices
                READ( MDEV,* ) MACINDX

!................  If MACT packet applies to this source, compute factor
                IF ( MACINDX .GT. 0 ) THEN

!......................  Calculate inventory baseline efficiency
                    CTLEFF = 0.
                    RULEFF = 1.
                    RULPEN = 1.
                    IF ( NCE .GT. 0 ) CTLEFF = DATVAL( S,NCE )
                    IF ( NRE .GT. 0 ) RULEFF = DATVAL( S,NRE )
                    IF ( NRP .GT. 0 ) RULPEN = DATVAL( S,NRP )
                    INVEFF = CTLEFF * RULEFF * RULPEN

                    EXSTCEFF = MACEXEFF( MACINDX )
                    NEWCEFF  = MACNWEFF( MACINDX )
                    NEWFRAC  = MACNWFRC( MACINDX )

!.....................  Compute factor for existing sources
                    IF( EXSTCEFF > INVEFF ) THEN
                        FAC = (1-NEWFRAC) * ((1-EXSTCEFF)/(1-INVEFF))
                    ELSE
                        FAC = (1-NEWFRAC)
                    END IF

!.....................  Compute factor for new sources and sum
!                       with factor for existing sources
                    IF( NEWCEFF > INVEFF ) THEN
                        FAC = FAC + NEWFRAC * ((1-NEWCEFF)/(1-INVEFF))
                    ELSE
                        FAC = FAC + NEWFRAC
                    END IF

!...................  Overwrite temporary file line with new info
                    E2 = E1 * FAC
                    FACTOR( S ) = FACTOR( S ) * FAC

                    WRITE( ODEV(4),93300 ) 1, PNAM, E1, E2, FAC
                    APPLFLAG = .TRUE.

                    IF( FAC .GT. 1. ) THEN
                        CSRC = CSOURC( S )
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                        WRITE( MESG,94110 ) 'WARNING: MACT ' //&
                        &       'packet record number', MACINDX,&
                        &       'is increasing ' // CRLF() // BLANK10 //&
                        &       'emissions by factor', FAC, 'for:'&
                        &       // CRLF() // BLANK10 //&
                        &       BUFFER( 1:L2 ) // ' POL:' // PNAM
                        CALL M3MESG( MESG )
                    END IF
                ELSE
                    WRITE( ODEV(4),93300 ) 0, 'D', 0., 0., 0.

                END IF

            END DO ! end source loop

        END IF

!.............................................................................
!...........  Apply /CONTROL/ packet controls if present for the current
!             pollutant
!.............................................................................
        IF ( CFLAG .AND. PCTLFLAG( I, 1 ) ) THEN

!................  Allocate arrays to store control eff., rule eff. and
!                  rule penetration for each source
            IF( ALLOCATED( CEFF ) ) THEN
                DEALLOCATE( CEFF, REFF, RPEN )
            END IF
            ALLOCATE( CEFF( NSRC ) )
            ALLOCATE( REFF( NSRC ) )
            ALLOCATE( RPEN( NSRC ) )
            CEFF = 0.
            REFF = 1.
            RPEN = 1.

!...............  Loop through sources
            DO S = 1, NSRC

!...............  Initialize control info. to base values
                CEFF( S ) = BASECEFF( S )
                REFF( S ) = BASEREFF( S )
                RPEN( S ) = BASERPEN( S )

                E1  = DATVAL( S,E ) * FACTOR( S )
                FAC = 1.

!................  Read temporary input file indices
                READ( CDEV,* ) CTLINDX

!................  If control packet applies to this source, compute factor
                IF ( CTLINDX .GT. 0 ) THEN
                    CTLEFF = FACCEFF( CTLINDX )
                    RULEFF = FACREFF( CTLINDX )
                    RULPEN = FACRLPN( CTLINDX )

                    NEWFAC = ( 1.0 - CTLEFF*RULEFF*RULPEN )
                    OLDFAC = 1.0
                    IF( BACKOUT( S ) .NE. 0. ) OLDFAC= 1.0/BACKOUT( S )

                    LREPORT = .FALSE.

!.....................  Check if this is a "replace" entry
                    IF( CTLRPLC( CTLINDX ) ) THEN

!.........................  Check if new reduction is more than old reduction
                        IF( LCOMPARE ) THEN

!.............................  If so, then backout old reduction and apply new
                            IF ( NEWFAC .LT. OLDFAC ) THEN
                                FAC = BACKOUT( S )* NEWFAC

!.............................  Store control information for this source
                                CEFF( S ) = CTLEFF
                                REFF( S ) = RULEFF
                                RPEN( S ) = RULPEN
                                LREPORT   = .TRUE.

!.............................  Otherwise, do nothing (keep with old reduction)
                            ELSE
                                CSRC = CSOURC( S )
                                CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                                WRITE( MESG,94110) 'WARNING: CONTROL '//&
                                &'packet record number', CTLINDX,&
                                &'control ignored for source because '//&
                                &'base year controls are greater for:'//&
                                &CRLF()// BLANK10//&
                                &BUFFER( 1:L2 ) // ' POL:' // PNAM
                                CALL M3MESG( MESG )

                                FAC = 1.0

!.............................  Store control information for this source
                                CEFF( S ) = BASECEFF( S )
                                REFF( S ) = BASEREFF( S )
                                RPEN( S ) = BASERPEN( S )

                            ENDIF

!.........................  If not checking, then always backout old reduction
!                           and apply new reduction
                        ELSE
                            FAC = BACKOUT( S )* NEWFAC

                            CEFF( S ) = CTLEFF
                            REFF( S ) = RULEFF
                            RPEN( S ) = RULPEN
                            LREPORT   = .TRUE.

                        END IF

                        FACTOR( S ) = FAC
                        E1 = DATVAL( S,E )

                        REPCTYPE = "R"

!.....................  Otherwise, apply as additional controls
                    ELSE
                        FAC = ( 1.0 - CTLEFF*RULEFF*RULPEN )
                        FACTOR( S ) = FACTOR( S ) * FAC

!.........................  For this calculation, put entire
!                           additional reduction (REFF, RPEN, CEFF)
!                           into new CEFF and use base REFF and RPEN
                        CEFF( S ) = 1 - FACTOR( S )  ! simplify, while awaiting correct formula
                        REFF( S ) = 1.  ! Set to 1, while awaiting correct formula
                        RPEN( S ) = 1.  ! Set to 1, while awaiting correct formula
                        LREPORT   = .TRUE.
                        REPCTYPE = "A"

                    END IF

!...................  Write out report entry if control packet entry
!                     has been applied
                    IF( LREPORT ) THEN

                        SELECT CASE( CATEGORY )
                          CASE( 'POINT' )
                            IF ( LREPHDR ) THEN
                                WRITE(RDEV2,93421)
                                LREPHDR = .FALSE.
                            END IF

                            WRITE( RDEV2, 93420 )&
                            &   TRIM(CSOURC(S)(PTBEGL3(1):PTENDL3(1))),&
                            &   TRIM(CSOURC(S)(PTBEGL3(2):PTENDL3(2))),&
                            &   TRIM(CSOURC(S)(PTBEGL3(3):PTENDL3(3))),&
                            &   TRIM(CSOURC(S)(PTBEGL3(4):PTENDL3(4))),&
                            &   TRIM(CSOURC(S)(PTBEGL3(5):PTENDL3(5))),&
                            &   TRIM(CSOURC(S)(PTBEGL3(6):PTENDL3(6))),&
                            &   TRIM(PNAM), 100.*(1.-OLDFAC),&
                            &   100.*CEFF(S)*REFF(S)*RPEN(S), REPCTYPE,&
                            &   TRIM( CTLCOMT( CTLINDX ) )
                          CASE( 'AREA' )
                            IF ( LREPHDR ) THEN
                                WRITE(RDEV2,93423)
                                LREPHDR = .FALSE.
                            END IF
                            WRITE( RDEV2, 93422 )&
                            &   TRIM(CSOURC(S)(ARBEGL3(1):ARENDL3(1))),&
                            &   TRIM(CSOURC(S)(ARBEGL3(2):ARENDL3(2))),&
                            &   TRIM(PNAM), 100.*(1.-OLDFAC),&
                            &   100.*CEFF(S)*REFF(S)*RPEN(S), REPCTYPE,&
                            &   TRIM( CTLCOMT( CTLINDX ) )

                          CASE( 'MOBILE' )
                            IF ( LREPHDR ) THEN
                                WRITE(RDEV2,93425)
                                LREPHDR = .FALSE.
                            END IF
                            WRITE( RDEV2, 93424 )&
                            &   TRIM(CSOURC(S)(MBBEGL3(1):MBENDL3(1))),&
                            &   TRIM(CSOURC(S)(MBBEGL3(2):MBENDL3(2))),&
                            &   TRIM(CSOURC(S)(MBBEGL3(3):MBENDL3(3))),&
                            &   TRIM(CSOURC(S)(MBBEGL3(4):MBENDL3(4))),&
                            &   TRIM(CSOURC(S)(MBBEGL3(5):MBENDL3(5))),&
                            &   TRIM(PNAM), 100.*(1.-OLDFAC),&
                            &   100.*CEFF(S)*REFF(S)*RPEN(S), REPCTYPE,&
                            &   TRIM( CTLCOMT( CTLINDX ) )

                        END SELECT

                    END IF

!...................  Overwrite temporary file line with new info
                    E2 = E1 * FAC

                    WRITE( ODEV(1),93300 ) 1, PNAM, E1, E2, FAC
                    APPLFLAG = .TRUE.

                    IF( FAC .GT. 1. ) THEN
                        CSRC = CSOURC( S )
                        CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                        WRITE( MESG,94110 ) 'WARNING: CONTROL '//&
                        &       'packet record number', CTLINDX,&
                        &       'is increasing ' // CRLF()// BLANK10//&
                        &       'emissions by factor', FAC, 'for:'&
                        &       //CRLF()// BLANK10//&
                        &       BUFFER( 1:L2 ) // ' POL:' // PNAM
                        CALL M3MESG( MESG )
                    END IF

!..................  Overwrite temporary file line with new info
                ELSE
                    WRITE( ODEV(1),93300 ) 0, 'D', 0., 0., 0.

                END IF

            END DO ! end source loop

        END IF

!.............................................................................
!............  Apply /CTG/ packet
!.............................................................................
        IF ( GFLAG .AND. PCTLFLAG( I, 2 ) ) THEN

!...............  Compute CTG factor
            DO S = 1, NSRC

                E1  = DATVAL( S,E ) * FACTOR( S )
                FAC = 1.

!................  Read CTG packet index from input temporary file
                READ( GDEV,* ) CTGINDX

!................  If CTG packet applies to this source, compute factor
                IF ( CTGINDX .GT. 0 ) THEN
                    CUTOFF = CUTCTG ( CTGINDX )
                    CTGFAC = FACCTG ( CTGINDX )
                    MACT   = FACMACT( CTGINDX )
                    RACT   = FACRACT( CTGINDX )

!.....................  Check to see if emissions exceed cutoff and if
!                       necessary, apply controls
!.....................  The comparison emission value already has controls
!                       from the /CONTROL/ packet
                    IF ( E1 .GT. CUTOFF ) THEN

!........................  Initialize CTG factor with base factor from packet
                        FAC = CTGFAC

!........................  Compute output emissions with /CONTROL/ and base CTG
                        E2  = E1*CTGFAC

!........................  If emissions still exceed cutoff, apply second
!                          CTG factor
                        IF ( E2 .GT. CUTOFF ) THEN

!...........................  Use MACT factor if it is defined
                            IF ( MACT .GT. 0 ) THEN
                                CTGFAC2 = MACT

!...........................  Otherwise, use RACT factor
                            ELSE IF ( RACT .GT. 0 ) THEN
                                CTGFAC2 = RACT

!...........................  Otherwise, set to cutoff value
                            ELSE
                                CTGFAC2 = CUTOFF / E2

                            END IF

                            FAC = FAC * CTGFAC2
                            E2  = E2  * CTGFAC2

                        END IF

!........................  Compute aggregate factor for current source
                        FACTOR( S ) = FACTOR( S ) * FAC

!........................  Write output temporary file line with new info
                        WRITE( ODEV(2),93300 ) 1, PNAM, E1, E2, FAC
                        APPLFLAG = .TRUE.

!.........................  Write warning if emissions have increased
                        IF( FAC .GT. 1. ) THEN
                            CSRC = CSOURC( S )
                            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                            WRITE( MESG,94110 ) 'WARNING: CTG '//&
                            &   'packet record number', CTGINDX,&
                            &   'is increasing ' // CRLF()// BLANK10//&
                            &   'emissions by factor', FAC, 'for:' //&
                            &   CRLF()// BLANK10//&
                            &   BUFFER( 1:L2 ) // ' POL:' // PNAM
                            CALL M3MESG( MESG )
                        END IF

!.....................  If no controls, then write output temporary line only
                    ELSE
                        WRITE( ODEV(2),93300 ) 0, 'D', 0., 0., 0.

                    END IF

                ELSE ! If K = 0 (for sources not having applied "CTG"

                    WRITE( ODEV(2),93300 ) 0, 'D', 0., 0., 0.

                END IF

            END DO ! end source loop

        END IF

!.............................................................................
!............  Apply /ALLOWABLE/ packet
!.............................................................................
        IF ( LFLAG .AND. PCTLFLAG( I, 3 ) ) THEN

!...........  Process ALW packet
            DO S = 1, NSRC

                E1 = DATVAL( S,E ) * FACTOR( S )

!.................  Read allowable packet index from input tmp file
                READ( LDEV,* ) ALWINDX

                FAC = 1.
!..................  If ALLOWABLE packet applies to this source, compute factor
                IF ( ALWINDX .GT. 0 ) THEN
                    CAP     = EMCAPALW( ALWINDX )
                    REPLACE = EMREPALW( ALWINDX )

!.....................  Both Cap value and Replace value are defined, then
!                       compare emissions to Cap and set factor with Replace.
                    IF ( CAP .GE. 0. .AND. REPLACE .GE. 0. ) THEN

                        IF ( DATVAL( S,E ) .GT. CAP ) THEN
                            FAC = REPLACE/DATVAL( S,E )
                        END IF

!.....................  Only Cap value is defined, then compare emissions to Cap
!                       set factor with Cap
                    ELSE IF ( CAP .GE. 0. .AND. REPLACE .LT. 0. ) THEN

                        IF ( DATVAL( S,E ) .GT. CAP ) THEN
                            FAC = CAP/DATVAL( S,E )
                        END IF

!.....................  Only Replace value is defined, then set factor with
!                       Replace.
                    ELSE IF ( CAP .LT. 0. .AND. REPLACE .GE. 0. ) THEN

                        IF ( DATVAL( S,E ) .GT. REPLACE ) THEN
                            FAC = REPLACE/DATVAL( S,E )
                        END IF

                    END IF

!.....................  Write output temporary file line with new info
                    E2 = E1 * FAC
                    FACTOR( S ) = FACTOR( S ) * FAC

                    IF( FAC .NE. 1. ) THEN
                        WRITE( ODEV(3),93300 ) 1, PNAM, E1, E2, FAC
                        APPLFLAG = .TRUE.

!.........................  Write warning if emissions have increased
                        IF( FAC .GT. 1. ) THEN
                            CSRC = CSOURC( S )
                            CALL FMTCSRC( CSRC, NCHARS, BUFFER, L2 )
                            WRITE( MESG,94110 ) 'WARNING: ALLOWABLE '//&
                            &   'packet record number', ALWINDX,&
                            &   'is increasing ' // CRLF()// BLANK10//&
                            &   'emissions by factor', FAC, 'for:' //&
                            &   CRLF()// BLANK10//&
                            &   BUFFER( 1:L2 ) // ' POL:' // PNAM
                            CALL M3MESG( MESG )
                        END IF

                    ELSE
                        WRITE( ODEV(3),93300 ) 0, 'D', 0., 0., 0.
                    END IF

!..................  If no controls, then write output temporary line only
                ELSE
                    WRITE( ODEV(3),93300 ) 0, 'D', 0., 0., 0.

                END IF

            END DO ! end source loop

        END IF

!.............  Store output emissions for groups
!.............  This must be in a separate loop to account for all possible
!               combinations of packets
        DO S = 1, NSRC

            E1 = DATVAL( S,E )
            E2 = E1 * FACTOR( S )
            J  = GRPINDX( S )
            GRPOUTEM( J,I ) = GRPOUTEM( J,I ) + E2

!.................  Flag group if emissions are different
            IF( E1 .NE. E2 ) GRPFLAG( J ) = .TRUE.

        END DO

!.............  Open control matrix, if needed and if not opened before
        IF( APPLFLAG .AND. .NOT. OPENFLAG ) THEN
            CALL OPENCMAT( ENAME, 'MULTIPLICATIVE', MNAME )
            OPENFLAG = .TRUE.
        END IF

!.............  Write multiplicative controls for current pollutant
        IF( OPENFLAG ) THEN
            IF( .NOT. WRITESET( MNAME, PNAM, ALLFILES, 0, 0,&
            &                    FACTOR )                     ) THEN
                MESG = 'Failed to write multiplicative control ' //&
                &       'factors for pollutant ' // PNAM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            L = LEN_TRIM( PNAM )
            IF( .NOT. WRITESET( MNAME, 'CE_'//PNAM(1:L), ALLFILES, 0, 0,&
            &                    CEFF )                     ) THEN
                MESG = 'Failed to write control efficiency ' //&
                &       'for pollutant ' // PNAM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( .NOT. WRITESET( MNAME, 'RE_'//PNAM(1:L), ALLFILES, 0, 0,&
            &                    REFF )                     ) THEN
                MESG = 'Failed to write rule effectivness ' //&
                &       'for pollutant ' // PNAM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF

            IF( .NOT. WRITESET( MNAME, 'RP_'//PNAM(1:L), ALLFILES, 0, 0,&
            &                    RPEN )                     ) THEN
                MESG = 'Failed to write rule penetration ' //&
                &       'for pollutant ' // PNAM
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
        END IF

    END DO ! end pollutant loop

!.........  Write out controlled facilities report for point sources
    IF( APPLFLAG ) THEN
        WRITE( RDEV, 93000 ) 'Processed as '// CATDESC// ' sources'
        IF ( PYEAR .GT. 0 ) THEN
            WRITE( RDEV, 93390 ) 'Projected inventory year ', PYEAR
        ELSE
            WRITE( RDEV, 93390 ) 'Base inventory year ', BYEAR
        END IF
        IF( MFLAG ) WRITE( RDEV, 93000 )&
        &            'Controls applied with /MACT/ packet'
        IF( CFLAG ) WRITE( RDEV, 93000 )&
        &            'Controls applied with /CONTROL/ packet'
        IF( GFLAG ) WRITE( RDEV, 93000 )&
        &            'Controls applied with /CTG/ packet'
        IF( LFLAG ) WRITE( RDEV, 93000 )&
        &            'Controls applied with /ALLOWABLE/ packet'
        IF( LAVEDAY ) THEN
            WRITE( RDEV,93000 ) 'Average day data basis in report'
        ELSE
            WRITE( RDEV,93000 ) 'Annual total data basis in report'
        END IF

        IF ( CATEGORY .EQ. 'POINT' ) THEN
            WRITE( RDEV, 93000 )&
            &   'Emissions by controlled facility before and ' //&
            &   'after controls'
        ELSE
            WRITE( RDEV, 93000 )&
            &   'Emissions by controlled state/SCC before and ' //&
            &   'after controls'
        END IF

        IF ( CATEGORY .EQ. 'POINT' ) THEN
            WRITE( RDEV, 93400 )&
            &     ( PNAMMULT( I ), PNAMMULT( I ), I=1,NVCMULT )
        ELSE
            WRITE( RDEV, 93402 )&
            &     ( PNAMMULT( I ), PNAMMULT( I ), I=1,NVCMULT )
        END IF

        PNAM = '[tons/day]'
        WRITE( RDEV, 93405 ) ( PNAM, PNAM, I=1,NVCMULT )

        I = 26 + NVCMULT * 43
        WRITE( RDEV,93000 ) REPEAT( '-', I )

        PIDX = 0
        DO S = 1, NSRC
            K   = GRPSTIDX( S )
            IDX = GRPINDX ( K )

            IF( IDX .NE. PIDX .AND.&
            &    GRPFLAG( IDX )      ) THEN

                J = FIPLEN3 + 1
                SELECT CASE( CATEGORY )
                  CASE ( 'AREA' , 'MOBILE' )
                    CSRC = GRPCHAR( S )
                    WRITE( RDEV, 93410 )&
                    &    CSRC( 1:STALEN3 ),&
                    &    CSRC( STALEN3+1:STALEN3+SCCLEN3 ),&
                    &    ( GRPINEM ( IDX,I ), GRPOUTEM( IDX,I ),&
                    &      I = 1, NVCMULT )
                  CASE ( 'POINT' )
                    CSRC = CSOURC( S )
                    WRITE( RDEV, 93412 )&
                    &    CSRC( 1:FIPLEN3 ), CSRC( J:FPLLEN3 ),&
                    &    ( GRPINEM ( IDX,I ), GRPOUTEM( IDX,I ),&
                    &      I = 1, NVCMULT )
                END SELECT

                PIDX = IDX
            END IF
        END DO

    END IF

    IF( .NOT. APPLFLAG ) THEN

        MESG = 'WARNING: No MACT, CONTROL, EMS_CONTROL, CTG, or ' //&
        &       'ALLOWABLE packet entries match inventory.'
        CALL M3MSG2( MESG )

        MESG = 'WARNING: Multiplicative control will not be output!'
        CALL M3MSG2( MESG )

        WRITE( RDEV, 93000 ) 'No MACT, CONTROL, EMS_CONTROL, ' //&
        &       'CTG, or ALLOWABLE packet entries matched inventory.'

    END IF

!........  Change input temporary file unit numbers to be output files for
!          report processing.
    IF( CDEV .GT. 0 ) CDEV = ODEV( 1 )
    IF( GDEV .GT. 0 ) GDEV = ODEV( 2 )
    IF( LDEV .GT. 0 ) LDEV = ODEV( 3 )
    IF( MDEV .GT. 0 ) MDEV = ODEV( 4 )

!.........  Deallocate local memory
    DEALLOCATE( BACKOUT, DATVAL, FACTOR )

    DEALLOCATE( GRPINDX )

    DEALLOCATE( GRPFLAG, GRPINEM, GRPOUTEM )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

93300 FORMAT( I2, 1X, '"', A, '"', 3( 1X, E12.5 ) )

93390 FORMAT( A, I4.4 )

93400 FORMAT( ' Region;', 5X, 'Facility ID;', 1X,&
    &        100( '  In ', A16, ';', 1X, 'Out ', A16, :, ';' ) )

93402 FORMAT( '  State;', 5X, '        SCC;', 1X,&
    &        100( '  In ', A16, ';', 1X, 'Out ', A16, :, ';' ) )

93405 FORMAT( 7X, ';', 16X, ';', 1X,&
    &        100( 2X, A16, 3X, ';', 1X, A16, :, 4X, ';' ))

93410 FORMAT( 4X, A3, ';', 6X, A10, ';', 1X,&
    &        100( 10X, E11.4, :, ';' ))

93412 FORMAT( 1X, A6, ';', 1X, A15, ';', 1X,&
    &        100( 10X, E11.4, :, ';' ))

93420 FORMAT( 7(A, ','), 2(F8.2,','),A,',',A )
93421 FORMAT( 'FIP,PLANTID,CHAR1,CHAR2,CHAR2,CHAR4,SCC,POL,',&
    &        'BASCRR,NEWCRR,REPLACE,COMMENT' )

93422 FORMAT( 3(A, ','), 2(F8.2,','),A,',',A )
93423 FORMAT( 'FIP,SCC,POL,BASCRR,NEWCRR,REPLACE,COMMENT' )

93424 FORMAT( 6(A, ','), 2(F8.2,','),A,',',A )
93425 FORMAT( 'FIP,ROADTYPE,LINK,VEHTYPE,SCC,NEWCRR,REPLACE,COMMENT' )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

94110 FORMAT( A, 1X, I8, 1X, 5( A, :, E11.3, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram writes and error message
!               and then terminates program execution if an error
!               is encountered reading control information from the
!               inventory file
    SUBROUTINE WRITE_MESG_EXIT( OUTNAME, PROGNAME )

!.............  Subprogram arguments
        CHARACTER(*),  INTENT (IN) :: OUTNAME   ! name of inventory
        ! variable that generated
        ! the error
        CHARACTER(16), INTENT (IN) :: PROGNAME  ! name of calling subroutine

!.............  Local variables
        CHARACTER(300)   MESG                   ! message buffer

!----------------------------------------------------------------------

        MESG = 'Error reading ' // OUTNAME // ' from inventory file'

        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END SUBROUTINE WRITE_MESG_EXIT

!----------------------------------------------------------------------

END SUBROUTINE GENMULTC
