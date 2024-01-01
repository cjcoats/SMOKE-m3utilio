
SUBROUTINE  RDASCC( ADEV, NDIM, NASC, ASCA7, ASCA3 )

    !***********************************************************************
    !  subroutine body starts at line  72
    !
    !  DESCRIPTION:
    !       Reads formatted actual-ASC file.
    !
    !  PRECONDITIONS REQUIRED:
    !       Actual-ASC file opened on unit ADEV.
    !       Actual-ASC file is sorted, and formatted (I7,I3)
    !
    !  REVISION  HISTORY:
    !       Prototype  12/96 by CJC for area-source submodel in SMOKE
    !       Version 11/2023 by CJC:  USE M3UTILIO and related changes
    !       Version 11/2023 by CJC:  USE M3UTILIO, ".f90" source format, and
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
    !****************************************************************************
    USE M3UTILIO

    IMPLICIT NONE

    !.......   ARGUMENTS and their descriptions: actually-occurring ASC table

    INTEGER, INTENT (IN   ) :: ADEV              !  unit number for actual-ASC file
    INTEGER, INTENT (IN   ) :: NDIM              !  max dimensioned number of ASCs
    INTEGER, INTENT (  OUT) :: NASC              !  actual number of ASCs returned
    INTEGER, INTENT (  OUT) :: ASCA7( NDIM )     !  leading-7 digits
    INTEGER, INTENT (  OUT) :: ASCA3( NDIM )     !  trailing-3 digits

    !.......   SCRATCH LOCAL VARIABLES and their descriptions:

    INTEGER         IREC                !  input line (record) number
    INTEGER         IOS                 !  I/O Status
    INTEGER         I                   !  loop counter
    INTEGER         ID7,  ID3
    INTEGER         LID7, LID3
    LOGICAL         EFLAG       !  input error flag
    CHARACTER(300)  MESG        !  message buffer for M3MESG() and M3EXIT()

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDASCC'     ! program name

    !***********************************************************************
    !   begin body of subroutine  RDASCC

    CALL M3MSG2( 'Reading ACTUAL ASCs file...' )

    IREC  =  0
    I     =  0
    NASC  =  0
    LID3  = -1
    LID7  = -1
    EFLAG = .FALSE.

11  CONTINUE                            !  head of the ADEV-read loop

        READ( ADEV, 93020, END=99, IOSTAT=IOS ) ID7, ID3

        IREC = IREC + 1

        IF ( IOS .GT. 0 ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'I/O error', IOS,       &
                'reading ACTUAL ASC file at line', IREC
            CALL M3MESG( MESG )
            GO TO  11                       !  to head of loop

        ELSE IF ( ( LID7 .GT. ID7 ) .OR.                &
                  ( LID7 .EQ. ID7 .AND. LID3 .GE. ID3 ) ) THEN

            EFLAG = .TRUE.
            WRITE( MESG,94010 ) 'ASC table out of order at line', IREC
            CALL M3MESG( MESG )
            GO TO  11

        END IF                  !  if i/o error; else if out-of-order

        I = I + 1
        IF ( I .LE. NDIM ) THEN

            ASCA7( I ) = ID7
            ASCA3( I ) = ID3

        END IF                  !  if I in bounds

        LID7 = ID7
        LID3 = ID3
    GO TO  11                       !  to head of loop

99  CONTINUE                            !  end of the ADEV-read loop

    WRITE( MESG,94010 ) 'Dimensioned ASC TABLE size', NDIM, ' Actual size', I
    CALL M3MSG2( MESG )

    IF ( I .GT. NDIM ) THEN
        CALL M3EXIT( 'SPCAMAT', 0, 0, 'ACTUAL ASC table overflow', 2 )
    ELSE IF ( EFLAG ) THEN
        CALL M3EXIT( 'SPCAMAT', 0, 0, 'Error reading ACTUAL ASC file.', 2 )
    END IF

    NASC = I

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.......   Formatted file I/O formats............ 93xxx

93020 FORMAT( I7, I3 )


    !.......   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I10, :, 2X ) )


END  SUBROUTINE  RDASCC

