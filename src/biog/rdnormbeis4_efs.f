
        SUBROUTINE RDNORMBEIS4_EFS( NLINES, NSEF, FDEV, NVEG, VGLIST,
     &                      BIOTYPES,LINDX, WNTF, FACS, AG_YN  )

C***********************************************************************
C  subroutine body starts at line XX
C
C  DESCRIPTION:
C    Reads in the BEIS4 emissions factors from the BEIS_FAC file formatted as
C    a csv file. These factors will be used in calculate normalized
C    emissions.
C
C  PRECONDITIONS REQUIRED:
C
C  REVISION  HISTORY:
C       03/22 protoype by J. Vukovich
C       10/2023 Adapted to USE M3UTILIO by Carlie J. Coats, Jr., UNCIE
C
C***********************************************************************

        USE M3UTILIO

        IMPLICIT NONE

C...........   ARGUMENTS and their descriptions: actually-occurring CSV table

        INTEGER, INTENT (IN)  :: NSEF    !  no. biogenic emission factors
        INTEGER, INTENT (IN)  :: FDEV    !  unit number input file
        INTEGER, INTENT (IN)  :: NVEG    !  no. veg types
        INTEGER, INTENT (IN)  :: NLINES  !  no. lines in input file
        CHARACTER(16), INTENT(IN) :: VGLIST(NVEG)
        CHARACTER(5) , INTENT(IN) :: BIOTYPES(NSEF)

C........Assign output vars
        INTEGER, INTENT (OUT) :: LINDX( NVEG )       ! leaf area index
        REAL,    INTENT (OUT) ::  WNTF( NVEG )       ! winter factor
        REAL,    INTENT (OUT) ::  FACS( NVEG, NSEF ) ! emis facs
        LOGICAL, INTENT (OUT) :: AG_YN( NVEG )       ! ag veg or not

        INTEGER,       PARAMETER :: MXSEG    = 4        ! # of potential line segments
        CHARACTER(16), PARAMETER :: PROGNAME = 'RDNORMBEIS4_EFS' ! program name

        LOGICAL         CHKHDR              ! check the header BEISFAC4BELD5
        LOGICAL         EFLAG               !  error flag

        INTEGER         I, J               !  counters
        INTEGER         ISTAT              !  iostat error
        INTEGER         INDEX, SINDEX

        CHARACTER(16)   VTYP,VGID,VUNIT
        REAL            VALU
        CHARACTER(50)   SEGMENT( 4 )   ! Segments of parsed lines
        CHARACTER(300)  MESG             !  message buffer
        CHARACTER(300)  LINE             !  buffer for variables


C***********************************************************************
C   begin body of subroutine RDNORMBEIS4_EFS

C.........  Set number of potential line segments
        CHKHDR = .FALSE.
        EFLAG  = .FALSE.

C.......... Read in emissions factors for each veg id
        DO I = 1, NLINES

          READ( FDEV, 93030, IOSTAT=ISTAT )  LINE

          IF ( ISTAT .NE. 0 ) THEN
                EFLAG = .TRUE.
                WRITE( MESG,94010 )
     &              'Error', ISTAT,
     &              'reading EMISSION FACTOR file at line', I
               CALL M3MESG( MESG )
          END IF

C............. Finding header which must equal "#BELD6" on first line

          IF ( I .EQ. 1 ) THEN
             IF ( LINE( 1:6 ) .NE. "#BELD6" ) THEN
               MESG = "Header line not equal to #BELD6; check " //
     &          "if BEIS_FAC file is for BELD6 use"
                EFLAG = .TRUE.
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
             ENDIF
             CYCLE
          ENDIF

C.............  Separate the line of data into each part
          CALL PARSLINE( LINE, MXSEG, SEGMENT )

          VGID  = SEGMENT( 1 )
          VTYP  = SEGMENT( 2 )
          VUNIT = SEGMENT( 3 )
          IF ( VTYP .NE. "AG_YN" ) THEN
              VALU = STR2REAL( SEGMENT ( 4 ))
          ENDIF
          INDEX = INDEX1( VGID, NVEG, VGLIST )
          IF (INDEX .LE. 0) THEN
            MESG = 'WARNING: VEGETATION NAME: '//TRIM(VGID)//
     &       ' found in emission factor file but not in land use file'
            CALL M3MESG( MESG )
            CYCLE
          ENDIF

C...... Find biotype and assign factor to biotype
          SELECT CASE ( VTYP )

             CASE ('LAI')
                LINDX(INDEX) = VALU
             CASE ('WFAC')
                WNTF(INDEX) = VALU
             CASE ('AG_YN')
                AG_YN( INDEX ) = ( SEGMENT( 4 ) .EQ. "Y" )
             CASE DEFAULT
                SINDEX = INDEX1( VTYP, NSEF, BIOTYPES )

                IF (SINDEX .LE. 0) THEN
                     EFLAG = .TRUE.
                     MESG = 'ERROR: SPECIES NAME: '//TRIM(VTYP)// 'not found'
                     CALL M3MESG( MESG )	
                     CYCLE
                 ENDIF

                 FACS(INDEX,SINDEX) = VALU

          END SELECT

        ENDDO   !!  End loop on number of lines in EF file

        IF( EFLAG ) THEN
            MESG = 'Problem reading biogenic emissions factors file.'
            CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
        ENDIF

        RETURN

C******************  FORMAT  STATEMENTS   ******************************

C...........   Formatted file I/O formats............ 93xxx
93030   FORMAT( A )

C...........   Internal buffering formats............ 94xxx
94010   FORMAT( 10 ( A, :, I10, :, 2X ) )
94020   FORMAT( A )

        END SUBROUTINE RDNORMBEIS4_EFS
