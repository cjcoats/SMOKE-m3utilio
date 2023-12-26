
SUBROUTINE RDB4FAC_CSV( BEISV, NLINES, NSEF, FDEV, NVEG, VGLIST,&
&                    BIOTYPES,LINDX, LFAC, WNTF, LWGT, FACS  )

!***********************************************************************
!  subroutine body starts at line XX
!
!  DESCRIPTION:
!       Reads in the BEIS3 emissions factors from the BFAC file formatted as a csv file
!       for BELD4 only!
!
!  PRECONDITIONS REQUIRED:
!
!  REVISION  HISTORY:
!       07/13 protoype by G. Pouliot
!       07/14 cleaned up the code and add cycle commands for warning cases
!       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
!
!***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

!...........   ARGUMENTS and their descriptions: actually-occurring ASC table

    CHARACTER(11), INTENT(IN) :: BEISV   !  version of BEIS
    INTEGER,       INTENT(IN) :: NSEF    !  no. biogenic emission factors
    INTEGER,       INTENT(IN) :: FDEV    !  unit number input file
    INTEGER,       INTENT(IN) :: NVEG  !  no. veg types
    INTEGER,       INTENT(IN) :: NLINES  !  no. lines in input file
    CHARACTER(16), INTENT(IN) :: VGLIST(NVEG)
    CHARACTER(5) , INTENT(IN) :: BIOTYPES(NSEF)

    INTEGER, INTENT (OUT)     :: LINDX( NVEG )          ! leaf area index
    REAL,    INTENT (OUT)     ::  LFAC( NVEG )          ! leaf biomass
    REAL,    INTENT (OUT)     ::  WNTF( NVEG )          ! winter factor
    REAL,    INTENT (OUT)     ::  LWGT( NVEG )          ! specific leaf wgt
    REAL,    INTENT (OUT)     ::  FACS( NVEG, NSEF )    ! emis facs

    INTEGER, PARAMETER :: MXSEG = 4   ! # of potential line segments

    LOGICAL     CHKHDR            ! check the header BEISFAC4BELD5
    LOGICAL     EFLAG             !  error flag

    INTEGER     I, J               !  counters
    INTEGER     ISTAT              !  iostat error
    INTEGER     INDEX, SINDEX
    REAL        VALU
    CHARACTER(16)   VTYP,VGID,UNIT
    CHARACTER(50)   SEGMENT( 4 )   ! Segments of parsed lines
    CHARACTER(300)  MESG             !  message buffer
    CHARACTER(300)  LINE             !  buffer for variables

    CHARACTER(16), PARAMETER :: PROGNAME = 'RDB4FAC' ! program name

!***********************************************************************
!   begin body of subroutine RDB4FAC_CSV
!.......... Read in emissions factors for each veg id
    CHKHDR = .FALSE.
    EFLAG  = .FALSE.

    DO I = 1, NLINES

        READ( FDEV, 93030, IOSTAT=ISTAT )  LINE

        IF ( ISTAT .NE. 0 ) THEN
            EFLAG = .TRUE.
            WRITE( MESG,94010 )&
            &    'Error', ISTAT,&
            &    'reading EMISSION FACTOR file at line', I
            CALL M3MESG( MESG )
        END IF

!............   Check the header
        IF( BEISV == 'NORMBEIS360' ) THEN   ! BEIS v3.6 does not need a header from B360FAC file
            IF( LINE( 2:14 ) == 'BEISFAC4BELD5' ) THEN
                MESG = 'Error: B360FAC emission factors input file '//&
                &       'is not compatible for BEIS v3.61 but for v3.7.'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            CHKHDR = .TRUE.
        ELSE
            IF( LINE( 2:14 ) == 'BEISFAC4BELD5' ) THEN
                CHKHDR = .TRUE.
                CYCLE
            END IF
        END IF

        IF( I > 1 .AND. .NOT. CHKHDR ) THEN
            MESG = 'ERROR: #BEISFAC4BELD5 header is missing ' //&
            &       'from BEISFAC input file for BEIS v3.7'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        END IF

!.............  Separate the line of data into each part
        CALL PARSLINE( LINE, MXSEG, SEGMENT )

        VGID = SEGMENT( 1 )
        VTYP = SEGMENT( 2 )
        UNIT = SEGMENT( 3 )
        VALU = STR2REAL( SEGMENT ( 4 ) )

        INDEX = 0
        DO J = 1, NVEG
            IF (TRIM(VGID) .eq. TRIM(VGLIST(J)) ) THEN
                INDEX = J

            ENDIF
        ENDDO
        IF (INDEX .eq. 0) THEN
            MESG = 'WARNING: VEGETATION NAME: '//TRIM(VGID)//&
            &   ' found in emission factor file but not in land use file'
            CALL M3MESG( MESG )
            CYCLE
        ENDIF

        SELECT CASE (TRIM(VTYP))

          CASE ('LAI')
            LINDX(INDEX) = VALU
          CASE ('LFBIO')
            LFAC(INDEX) = VALU

          CASE ('WFAC')
            WNTF(INDEX) = VALU
          CASE ('SLW')
            LWGT(INDEX) = VALU
          CASE DEFAULT
            SINDEX = 0
            DO J = 1, NSEF

                IF (TRIM(VTYP) .eq. TRIM(BIOTYPES(J)) ) THEN
                    SINDEX = J
                ENDIF
            ENDDO


            IF (SINDEX .eq. 0) THEN
                EFLAG = .TRUE.
                MESG = 'Error: SPECIES NAME: '//TRIM(VTYP)//&
                &    'not found'
                CALL M3MESG( MESG )
                CYCLE
            ENDIF


            FACS(INDEX,SINDEX) = VALU

        END SELECT

    ENDDO

    IF( EFLAG ) THEN
        MESG = 'Problem reading biogenic emissions factors file.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    ENDIF

    RETURN

!******************  FORMAT  STATEMENTS   ******************************

!...........   Formatted file I/O formats............ 93xxx

93010 FORMAT( 8X, A16, A )
93020 FORMAT( A16, A )
93030 FORMAT( A )

!...........   Internal buffering formats............ 94xxx

94010 FORMAT( 10 ( A, :, I10, :, 2X ) )
94020 FORMAT( A )

END SUBROUTINE RDB4FAC_CSV
