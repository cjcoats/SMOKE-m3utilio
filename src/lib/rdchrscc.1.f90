
SUBROUTINE RDCHRSCC( FDEV, NSCC, SCCLIST )

!***********************************************************************
!  subroutine body starts at line
!
!  DESCRIPTION:
!       Reads character SCC codes
!
!  PRECONDITIONS REQUIRED:
!
!  REVISION  HISTORY:
!       Written  1/99 by M. Houyoux
!       Version 11/2023 by CJC:  USE M3UTILIO and related changes
!
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
!****************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

!...........   INCLUDES
    INCLUDE 'EMCNST3.EXT'

!...........   ARGUMENTS and their descriptions: actually-occurring ASC table

    INTEGER     , INTENT (IN   ) :: FDEV             !  unit number for actual-SCC file
    INTEGER     , INTENT (IN   ) :: NSCC             !  number of lines in actual-SCC file
    CHARACTER(*), INTENT (  OUT) :: SCCLIST( NSCC )  !  list of SCCs

!...........   SCRATCH LOCAL VARIABLES and their descriptions:

    INTEGER         I             !  counters and indices
    INTEGER         IOS           !  I/O Status
    INTEGER         IREC          !  input line (record) number

    LOGICAL      :: EFLAG         !  error flag

    CHARACTER(7)    FMTSCC        !  read format for SCCs
    CHARACTER(300)  MESG          !  message buffer
    CHARACTER(SCCLEN3) SBUFFA     !  read SCC buffer

    CHARACTER(16) :: PROGNAME = 'RDCHRSCC' ! program name

!***********************************************************************
!   begin body of subroutine  RDCHRSCC

!.........  Write SCC read format
    WRITE( FMTSCC, '(A,I2.2,A)' ) '(A', SCCLEN3, ')'

    EFLAG = .FALSE.
    IREC  = 0
    DO I = 1, NSCC    !  head of the FDEV-read loop

        READ( FDEV, FMTSCC, END=999, IOSTAT=IOS ) SBUFFA
        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'I/O error', IOS,&
            &    'reading SCC CODES file at line', IREC
            CALL M3MESG( MESG )
            CYCLE   !  to end of loop

        END IF      !  if i/o error; else if out-of-order

        SCCLIST( I ) = SBUFFA

    ENDDO       !  end of the FDEV-read loop

    IF( EFLAG ) THEN

        MESG = 'Problem reading actual SCC file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    ENDIF

    REWIND( FDEV )

    RETURN

!.........  Special exit from loop for end of file

999 WRITE( MESG,94010 ) 'Unexpected end-of-file at line', IREC,&
    &                    'of actual SCC file.'

    CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93020 FORMAT( I7, I3 )


!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I10, :, 2X ) )

END SUBROUTINE RDCHRSCC
