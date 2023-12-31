
SUBROUTINE ALOCPKTS( FDEV, WDEV, INYEAR, CPYEAR, PKTCNT, PKTBEG,&
&                     XRFCNT, LPTMP, LCTMP )

!***********************************************************************
!  subroutine body starts at line 107
!
!  DESCRIPTION:
!      This subroutine allocates memory for the valid control packets that are
!      in the input file
!
!  PRECONDITIONS REQUIRED:
!
!  SUBROUTINES AND FUNCTIONS CALLED:
!
!  REVISION  HISTORY:
!       Started 3/99 by M. Houyoux
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
!.........  This module contains the control packet data and control matrices
    USE MODCNTRL, ONLY: CUTCTG, FACCTG, FACMACT, FACRACT,&
    &                    ICTLEQUP, CCTLSIC, FACCEFF, FACREFF,&
    &                    FACRLPN, CALWSIC, FACALW, EMCAPALW,&
    &                    EMREPALW, CREASIC, EMREPREA, PRJFCREA,&
    &                    MKTPNREA, CSCCREA, CSPFREA, CPRJSIC,&
    &                    PRJFC, CTLRPLC, MACEXEFF, MACNWEFF, MACNWFRC,&
    &                    CMACSRCTYP, CTGCOMT, CTLCOMT, ALWCOMT,&
    &                    REACOMT, PRJCOMT, MACCOMT

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'   !  emissions constant parameters

!...........   SUBROUTINE ARGUMENTS:
!...........   Note: LPTMP and LCTMP needed only because PKTLOOP call needs them.
!                    They will not be set until after tmp files are written after
!                    calling PROCPKTS later in the main routine.

    INTEGER     , INTENT (IN) :: FDEV      ! in file unit number
    INTEGER     , INTENT (IN) :: WDEV      ! errors/warning file
    INTEGER     , INTENT (IN) :: INYEAR    ! year to project from
    INTEGER     , INTENT(OUT) :: CPYEAR    ! year to project to
    INTEGER     , INTENT(OUT) :: PKTCNT( NPACKET ) ! count of packet recs
    INTEGER     , INTENT(OUT) :: PKTBEG( NPACKET ) ! 1st line of pkt in file
    INTEGER     , INTENT(OUT) :: XRFCNT( NPACKET ) ! count of x-ref entries
    LOGICAL     , INTENT(OUT) :: LPTMP     ! true: projection tmp file written
    LOGICAL     , INTENT(OUT) :: LCTMP     ! true: control tmp file written

!...........   EXTERNAL FUNCTIONS:
    LOGICAL, EXTERNAL :: BLKORCMT

!...........   Logical names and unit numbers

    INTEGER         PDEV      ! file unit no. for tmp PROJ file
    INTEGER         CDEV      ! file unit no. for tmp CTL file
    INTEGER         GDEV      ! file unit no. for tmp CTG file
    INTEGER         LDEV      ! file unit no. for tmp ALW file
    INTEGER         MDEV      ! file unit no. for tmp MACT file

!...........   Other local variables

    INTEGER         I, J, K, L, L1, L2      ! counters and indices

    INTEGER         IOS       ! i/o error status
    INTEGER         IREC      ! line number
    INTEGER         SYEAR     ! tmp start year for packet

    LOGICAL      :: EFLAG  = .FALSE.   ! error flag
    LOGICAL      :: INSIDE = .FALSE.   ! true: read loop inside a packet
    LOGICAL      :: VALID  = .FALSE.   ! true: current packet is valid
    LOGICAL      :: YRSPEC = .FALSE.   ! true: packet is year-specific

    CHARACTER(7)              ACTION      ! buffer for PKTLOOP action
    CHARACTER(300)            LINE        ! read buffer for a line
    CHARACTER(300)            MESG        ! message buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'ALOCPKTS' ! program name

!***********************************************************************
!   Begin body of subroutine ALOCPKTS

!.........  Initialize packet count for all valid packets
    PKTCNT = 0   ! array

!.........  Write status message
    MESG = 'Scanning control/projection packets input file...'
    CALL M3MSG2( MESG )

!.........  Loop through control packets file and scan for packets. For each
!           packet type, count the number of entries.
!.........  For the projection packet, only pay attention to those that
!           have the same INYEAR.
!.........  Multiple packets of the same type are not permitted.  For
!           projections, this means that multiple packets for the same INYEAR
!           are not permitted.

    IREC = 0
    DO                ! head of FDEV read loop

        READ( FDEV, 93000, END = 101, IOSTAT=IOS ) LINE
        IREC = IREC + 1

        IF( IOS .GT. 0 ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'ERROR: I/O error', IOS,&
            &    'reading control packets file at line', IREC
            CALL M3MESG( MESG )
            CYCLE

        END IF

!.............  Skip blank lines and comment lines
        IF( BLKORCMT( LINE ) ) CYCLE

!.............  If inside packet...
        IF( INSIDE ) THEN

            J = INDEX( LINE, '!' )
            IF ( J .LE. 0 ) J = LEN_TRIM( LINE )
            I = INDEX( LINE( 1:J ), '/END/' )
            J = INDEX( LINE( 1:J ), '/'     )

!.................  Check for /END/ of packet
            IF( I .GT. 0 ) THEN
                INSIDE = .FALSE.

!.................  Encountered slash but not /END/
!                   It may be part of a point source characteristic, so check if
!                   it's at the beginning of the line
            ELSEIF( I .LE. 0 .AND. J .GT. 0 ) THEN
                IF( J .EQ. 1 .OR. LINE( 1:J ) .EQ. ' ' ) THEN
                    WRITE( MESG,94010 ) 'Problem at line', IREC,&
                    &       'of control packets file.' //&
                    &       CRLF() // BLANK10 //&
                    &       'Encountered a "/" before /END/.'
                    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
                END IF

!.................  For valid packets, count records
            ELSEIF( VALID ) THEN

                PKTCNT( K ) = PKTCNT( K ) + 1   ! Increment for packet K

            END IF

!.............  If outside packet, look for next packet header
        ELSE

            L  = LEN_TRIM( LINE )
            L1 = INDEX( LINE, '/' )
            L2 = INDEX( LINE( L1+1:L ), '/' ) + 1

!.................  Well-formed packet
            IF( L1 .GT. 0 .AND. L2 .GT. L1 + 1 ) THEN

!.....................  Check if this packet is in list, with special treatment
!                       for reactivity and projection packets.  Return position
!                       K in parameter list of valid packet names, the control
!                       packet year (if applicable), and a year-specific flag.
                CALL CHECK_PACKET( LINE, L1, L2, INYEAR, K,&
                &                   CPYEAR, YRSPEC )

!.....................  Print message because this year-specific packet does
!                       not apply to the inventory file
                IF( K .LE. 0 .AND. YRSPEC ) THEN
                    WRITE( MESG,94010 )&
                    &       'WARNING: Packet ' // LINE( L1:L2 ) //&
                    &       'does not apply to current' //&
                    &       CRLF() // BLANK10 // 'inventory year (',&
                    &       INYEAR, '). It will be skipped.'
                    CALL M3MSG2( MESG )
                    INSIDE = .TRUE.
                    VALID = .FALSE.

!.....................  Otherwise when K < 0, packet is not known
                ELSE IF( K .LE. 0 ) THEN
                    MESG = 'WARNING: Packet ' // LINE( L1:L2 ) //&
                    &       ' is not recognized. It will be skipped.'
                    CALL M3MSG2( MESG )
                    INSIDE = .TRUE.
                    VALID = .FALSE.

!.....................  When K is found, make sure it's the first time
                ELSE IF( PKTCNT( K ) .GT. 0 ) THEN
                    EFLAG = .TRUE.
                    WRITE( MESG,94010 )&
                    &       'ERROR: Second packet '// LINE( L1:L2 )//&
                    &       'found at line', IREC, 'of control file.'
                    CALL M3MSG2( MESG )
                    INSIDE = .TRUE.
                    VALID = .FALSE.

!.....................  Otherwise, store the line number this packet starts on
                ELSE
                    INSIDE = .TRUE.
                    VALID = .TRUE.
                    PKTBEG( K ) = IREC

                END IF

!.................  Badly-formed packet
            ELSE IF( L1 .GT. 0 .AND. L2 .LE. L1 ) THEN

                WRITE( MESG,94010 )&
                &       'WARNING: Badly formed packet at line',&
                &       IREC, 'in control packets file.'
                INSIDE = .TRUE.
                VALID = .FALSE.

            END IF  ! End well-formed packet or not

        END IF      ! End inside a packet or not

    END DO          ! End loop to scan and count records

101 CONTINUE        ! Exit from read loop

!.........  Abort of there was an error
    IF( EFLAG ) THEN
        MESG = 'Problem reading control packets file'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Allocate memory for control packet inputs...
!.........  CTG packet
    J = PKTCNT( 1 )
    ALLOCATE( CUTCTG( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CUTCTG', PROGNAME )
    ALLOCATE( FACCTG( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACCTG', PROGNAME )
    ALLOCATE( FACMACT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACMACT', PROGNAME )
    ALLOCATE( FACRACT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACRACT', PROGNAME )
    ALLOCATE( CTGCOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    CUTCTG  = 0.  ! array
    FACCTG  = 0.  ! array
    FACMACT = 0.  ! array
    FACRACT = 0.  ! array
    CTGCOMT = " " ! array

!.........  CONTROL packet
    J = PKTCNT( 2 )
    ALLOCATE( ICTLEQUP( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'ICTLEQUP', PROGNAME )
    ALLOCATE( CCTLSIC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CCTLSIC', PROGNAME )
    ALLOCATE( FACCEFF( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACCEFF', PROGNAME )
    ALLOCATE( FACREFF( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACREFF', PROGNAME )
    ALLOCATE( FACRLPN( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACRLPN', PROGNAME )
    ALLOCATE( CTLRPLC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTLRPLC', PROGNAME )
    ALLOCATE( CTLCOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    ICTLEQUP = 0   ! array
    CCTLSIC  = ' ' ! array
    FACCEFF  = 0.  ! array
    FACREFF  = 0.  ! array
    FACRLPN  = 0.  ! array
    CTLRPLC = .FALSE.  ! array
    CTLCOMT = " " ! array

!.........  ALLOWABLE packet
    J = PKTCNT( 3 )
    ALLOCATE( CALWSIC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CALWSIC', PROGNAME )
    ALLOCATE( FACALW( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'FACALW', PROGNAME )
    ALLOCATE( EMCAPALW( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMCAPALW', PROGNAME )
    ALLOCATE( EMREPALW( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMREPALW', PROGNAME )
    ALLOCATE( ALWCOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    CALWSIC = ' ' ! array
    FACALW   = 0. ! array
    EMCAPALW = 0. ! array
    EMREPALW = 0. ! array
    ALWCOMT = " " ! array

!.........  REACTIVITY packet
    J = PKTCNT( 5 )
    ALLOCATE( CREASIC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CREASIC', PROGNAME )
    ALLOCATE( EMREPREA( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'EMREPREA', PROGNAME )
    ALLOCATE( PRJFCREA( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'PRJFCREA', PROGNAME )
    ALLOCATE( MKTPNREA( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MKTPNREA', PROGNAME )
    ALLOCATE( CSCCREA( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSCCREA', PROGNAME )
    ALLOCATE( CSPFREA( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CSPFREA', PROGNAME )
    ALLOCATE( REACOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    CREASIC  = ' ' ! array
    EMREPREA = 0.  ! array
    PRJFCREA = 0.  ! array
    MKTPNREA = 0.  ! array
    CSCCREA  = ' ' ! array
    CSPFREA  = ' ' ! array
    REACOMT = " " ! array

!.........  PROJECTION packet
    J = PKTCNT( 6 )
    ALLOCATE( CPRJSIC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CPRJSIC', PROGNAME )
    ALLOCATE( PRJFC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'PRJFC', PROGNAME )
    ALLOCATE( PRJCOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    CPRJSIC = ' ' ! array
    PRJFC   = 0.  ! array
    PRJCOMT = " " ! array

!.........  MACT packet
    J = PKTCNT( 7 )
    ALLOCATE( CMACSRCTYP( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CMACSRCTYP', PROGNAME )
    ALLOCATE( MACEXEFF( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MACEXEFF', PROGNAME )
    ALLOCATE( MACNWEFF( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MACNWEFF', PROGNAME )
    ALLOCATE( MACNWFRC( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'MACNWFRC', PROGNAME )
    ALLOCATE( MACCOMT( J ), STAT=IOS )
    CALL CHECKMEM( IOS, 'CTGCOMT', PROGNAME )
    CMACSRCTYP = '00'  ! array
    MACEXEFF = 0.  ! array
    MACNWEFF = 0.  ! array
    MACNWFRC = 0.  ! array
    MACCOMT = " " ! array

!.........  Make sure that at least one packet is defined
    J = 0
    DO I = 1, NPACKET
        J = J + PKTCNT( I )
    END DO

    IF( J .EQ. 0 ) THEN
        MESG = 'No valid packets found in control packets file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

!.........  Now determine memory needed the cross-referencing portion of the
!           packets for each packet type.
!.........  For point sources, the memory needed is not necessarily equal to the
!           number of lines in the packet, because the SIC codes need to be
!           expanded to SCC codes.  So, for point sources process the packets to
!           find SIC NE 0 and SCC EQ 0, and expand the memory requirements
!           accordingly.

    ACTION = 'COUNT'
    CALL PKTLOOP( FDEV, PDEV, CDEV, GDEV, LDEV, MDEV, WDEV,&
    &              CPYEAR, ACTION, BLANK5, PKTCNT, PKTBEG,&
    &              XRFCNT, LPTMP, LCTMP )

!.........  Rewind file
    REWIND( FDEV )

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93000 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

!******************  INTERNAL SUBPROGRAMS  *****************************

CONTAINS

!.............  This internal subprogram checks the packet header(s) to find
!               valid packets in the input file.
    SUBROUTINE CHECK_PACKET( LINE, L1, L2, INYEAR, PKTIDX,&
    &                         OUTYEAR, YFLAG )

!.............  Subprogram arguments
        CHARACTER(*), INTENT (IN) :: LINE     ! local size for dimensioning
        INTEGER     , INTENT (IN) :: L1       ! start of packet in LINE
        INTEGER     , INTENT (IN) :: L2       ! end of packet in LINE
        INTEGER     , INTENT (IN) :: INYEAR   ! year for valid packets
        INTEGER     , INTENT(OUT) :: PKTIDX   ! packet index in master list
        INTEGER     , INTENT(OUT) :: OUTYEAR  ! year generated by packet(s)
        LOGICAL     , INTENT(OUT) :: YFLAG    ! flag for year-specific pkt

!.............  Local variables
        INTEGER    I, J, L, S1, S2    ! indices
        INTEGER    IOS             ! i/o status

        INTEGER, SAVE :: SAVYEAR = 0  ! saved output year

        LOGICAL, SAVE :: SFLAG     ! true: there is a saved output year

        CHARACTER(40)  :: BUFFER    ! packet buffer
        CHARACTER(300) :: MESG      ! message buffer

!----------------------------------------------------------------------

!.............  Initialize outputs
        OUTYEAR = 0
        YFLAG   = .FALSE.

!.............  If spaces are between slashes, reset comparison accordingly
        S1 = L1 + 1
        S2 = L2 - 1
        BUFFER = LINE( S1:S2 )
        J = INDEX( BUFFER, ' ' )
        IF( J .GT. 0 ) BUFFER = BUFFER( 1:J-1 )

!.............  Find packet in LINE in list of packets
        PKTIDX = 0
        DO I = 1, NPACKET

            L = LEN_TRIM( PKTLIST( I ) )
            IF( BUFFER .EQ. PKTLIST( I )( 1:L ) ) THEN
                PKTIDX = I
                EXIT
            END IF

        END DO

!.............  If packet is not recognized, return
        IF( PKTIDX .LE. 0 ) RETURN

!.............  Otherwise, continue...

!.............  Check for year information...
        L = L1 + LEN_TRIM( PKTLIST( PKTIDX ) ) + 1
        IF( L2 .GT. L ) THEN

            READ( LINE( L+1:L2 ), *, IOSTAT=IOS ) SYEAR, OUTYEAR

!.................  Misformated year information
            IF( IOS .GT. 0 ) THEN

                EFLAG = .TRUE.
                WRITE( MESG,94010 )&
                &    'ERROR: I/O error', IOS,&
                &    'reading packet years at line', IREC  ! IREC from main
                CALL M3MESG( MESG )
                PKTIDX = -1
                RETURN

!.................  Start year of packet does not apply (warning)
            ELSE IF( SYEAR .NE. INYEAR ) THEN

! NOTE: It is not this simple, because there might be multiple inventory years
!      in the input inventory.  Instead of comparing to the base year, compare
!      to a list of base years generated in GETSINFO.  The projection packet
!      data tables should then contain the start year.
!.................  Get type of projection entries: with year or without it (EPS)
                YFLAG = ENVYN( 'PROJECTION_YR_SPEC',&
                &       'Projection entries in year-specific format',&
                &       .TRUE., IOS )
                IF( YFLAG ) THEN
                    PKTIDX = -1
                    RETURN
                END IF

                SAVYEAR = OUTYEAR

!.................  Start and end year are the same... (warning)
            ELSE IF( SYEAR .EQ. OUTYEAR ) THEN
                WRITE( MESG,94010 )&
                &    'WARNING: Input year and output year are '//&
                &    'the same.'
                CALL M3MSG2( MESG )
                RETURN

!.................  Start year is fine, but output year disagrees with other pkt
            ELSE IF( SFLAG .AND. OUTYEAR .NE. SAVYEAR ) THEN ! (warning)
                WRITE( MESG,94010 )&
                &       'WARNING: Output year', OUTYEAR, 'from ' //&
                &       LINE( S1:S2 ) // ' packet, does not agree '//&
                &       'with previously encountered output year',&
                &       SAVYEAR, CRLF() // BLANK10 //&
                &       'Packet will be skipped.'
                CALL M3MESG( MESG )
                PKTIDX = -1
                RETURN

!.................  Packet is good
            ELSE

                SFLAG   = .TRUE.    ! year-specific has been encountered
                SAVYEAR = OUTYEAR   ! store saved year
                YFLAG   = .TRUE.

            END IF   ! End checks on year-specific records

        END IF       ! End year-specific processing

        OUTYEAR = SAVYEAR

        RETURN

!------------------- SUBPROGRAM FORMAT STATEMENTS ----------------------

!...........   Internal buffering formats............ 94xxx

94010   FORMAT( 10( A, :, I8, :, 1X ) )

    END SUBROUTINE CHECK_PACKET

END SUBROUTINE ALOCPKTS
