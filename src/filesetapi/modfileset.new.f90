
MODULE MODFILESET

    !***********************************************************************
    !  Module body starts at line 38
    !
    !  DESCRIPTION:
    !     This module contains the public and private variables and arrays
    !     needed to use the FileSetAPI.
    !
    !  PUBLIC SUBROUTINES AND FUNCTIONS:
    !       OPENSET:   Opens a file set; can be already opened, existing on disk, or new
    !       DESCSET:   Get description of file set
    !       CLOSESET:  Closes a file set
    !       PROMPTSET: Prompts the user for a logical file name, then tries
    !                  to open specified file set
    !
    !  REVISION HISTORY:
    !       Created 6/02 by C. Seppanen
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !       Version 1/2024 by CJC:  CONTAINed routines
    !***************************************************************************
    !
    ! Project Title: FileSetAPI
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
    !*************************************************************************
    USE M3UTILIO

    IMPLICIT NONE
    SAVE
    
    PUBLIC OPENSET, DESCSET, CLOSESET, PROMPTSET
    
    LOGICAL, EXTERNAL :: READSET     ! read from a file set
    LOGICAL, EXTERNAL :: WRITESET    ! write to a file set

    INTEGER, PUBLIC, PARAMETER :: ALLFILES = -1

    !.......  File set information
    INTEGER, PUBLIC              :: NVARSET                 ! total number of variables in the file set
    INTEGER, PUBLIC              :: NFILESET                ! total number of files in the file set
    INTEGER, PUBLIC, ALLOCATABLE :: VARS_PER_FILE( : )      ! number of variables per file

    !.......  Arrays for storing variable information (dim: NVARSET)
    INTEGER      , PUBLIC, ALLOCATABLE :: VTYPESET( : )      ! variable types
    CHARACTER(16), PUBLIC, ALLOCATABLE :: VNAMESET( : )      ! variable names
    CHARACTER(16), PUBLIC, ALLOCATABLE :: VUNITSET( : )      ! variable units
    CHARACTER(80), PUBLIC, ALLOCATABLE :: VDESCSET( : )      ! variable descriptions

    !.......  Internal wrapper data
    TYPE :: CHAR_PTR_ARRAY
        LOGICAL                :: RDONLY            ! read-only status
        CHARACTER(16), POINTER :: LNAMES( : )       ! logical file names
        CHARACTER(16), POINTER :: VARS( :,: )       ! variable names
    END TYPE

    INTEGER               , PUBLIC :: NOPENSETS = 0             ! total number of open file sets
    CHARACTER(16)         , PUBLIC :: RNAMES( MXFILE3 )         ! logical file names for open file sets
    TYPE( CHAR_PTR_ARRAY ), PUBLIC :: FILE_INFO( MXFILE3 )      ! file information for open file sets


    PRIVATE     ! everything else !


CONTAINS

    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-
    !!-=-=-=-=-  PUBLIC ROUTINES FIRST  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-
    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-

    LOGICAL FUNCTION OPENSET( ROOTNAME, FSTATUS, PGNAME )

        CHARACTER(*), INTENT(IN) :: ROOTNAME      ! logical file name for file set
        INTEGER,      INTENT(IN) :: FSTATUS       ! file mode/status (read/write, unknown)
        CHARACTER(*), INTENT(IN) :: PGNAME        ! name of calling program

        !......  External functions
        INTEGER,       EXTERNAL :: RMFILE
        CHARACTER(10), EXTERNAL :: GETCFDSC

        !......  Local variables
        INTEGER            I               ! counter
        INTEGER            IOS             ! I/O status
        INTEGER            FILEIDX         ! file index
        INTEGER            VARPOS          ! position in variable information arrays
        INTEGER            NFILEINT        ! number of files as an integer
        INTEGER            NVARINT         ! number of variables as an integer

        LOGICAL, SAVE   :: INITIAL = .TRUE.     ! true: first time through function
        LOGICAL            FILEEXIST       ! true: individual file exists on disk
        LOGICAL         :: SETEXIST = .FALSE.       ! true: file set exists on disk
        LOGICAL            TEMPEXIST       ! true: temporary file name exists on disk
        LOGICAL            EFLAG           ! true: error occured

        CHARACTER(2)       NFILESTR        ! number of files as a string
        CHARACTER(4)       NVARSTR         ! number of variables as a string
        CHARACTER(2)       INTBUF          ! integer string buffer
        CHARACTER(16)      ROOTNAME16      ! fixed length root file name
        CHARACTER(16)      PGNAME16        ! fixed length calling program name
        CHARACTER(16)      LNAME           ! temporary logical name
        CHARACTER(256)     ENVNAME         ! environment variable value for ROOTNAME
        CHARACTER(256)     TEMPNAME        ! temporary physical file name
        CHARACTER(256)     MESG            ! message buffer

        CHARACTER(24), PARAMETER :: FUNCNAME = 'MODFILESET/OPENSET'      ! function name

        !---------------------------------
        !  Begin body of function OPENSET
        !---------------------------------

        !......  Initialize arrays first time
        IF( INITIAL ) THEN
            RNAMES = CMISS3

            DO I = 1, MXFILE3
                NULLIFY( FILE_INFO( I )%VARS, FILE_INFO( I )%LNAMES )
            END DO

            INITIAL = .FALSE.
        END IF

        !......  Check length of file name
        IF( LEN_TRIM( ROOTNAME ) > 16 ) THEN
            MESG = 'Max file name length (16) exceeded for "' // ROOTNAME // '"'
            CALL M3MSG2( MESG )
            OPENSET = .FALSE.
            RETURN
        END IF

        !......  Check if file is already open
        ROOTNAME16 = ROOTNAME
        FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )

        !......  If file is already open, check status
        IF( FILEIDX /= 0 ) THEN

            !......  Trying to open as a new file
            IF( FSTATUS == FSNEW3 ) THEN
                MESG = 'File ' // ROOTNAME16 // ' already opened; ' //&
                       'cannot subsequently create in "NEW".'
                CALL M3WARN( FUNCNAME, 0, 0, MESG )
                OPENSET = .FALSE.
                RETURN

            !......  Trying to open as read/write or unknown when already readonly
            ELSE IF( FILE_INFO( FILEIDX )%RDONLY .AND.          &
                   ( FSTATUS == FSRDWR3 .OR. FSTATUS == FSUNKN3 ) ) THEN
                MESG = 'File ' // ROOTNAME16 // ' already opened READONLY; ' //&
                       'cannot subsequently open it for READ/WRITE.'
                CALL M3WARN( FUNCNAME, 0, 0, MESG )
                OPENSET = .FALSE.
                RETURN

            !......  Trying to open as unknown (already checked that status is not readonly)
            ELSE IF( FSTATUS == FSUNKN3 ) THEN

                !......  Check consistency of file set description
                IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                    MESG = 'Bad file description for file set ' // ROOTNAME16
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Check description against file info
                IF( .NOT. CHKFILESET( FILEIDX ) ) THEN
                    MESG = 'File description does not match for ' // 'file set ' // ROOTNAME16
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Loop through individual files
                VARPOS = 1
                DO I = 1, NFILESET

                    !........  Set necessary file description values
                    NVARS3D = VARS_PER_FILE( I )
                    VTYPE3D( 1:NVARS3D ) =  VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
                    VNAME3D( 1:NVARS3D ) = VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )

                    !........  Try to open file
                    IF( .NOT. OPEN3( FILE_INFO( FILEIDX )%LNAMES( I ), FSTATUS, FUNCNAME ) ) THEN
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    VARPOS = VARPOS + NVARS3D
                END DO

                !......  Successfully opened files
                OPENSET = .TRUE.
                RETURN

            !......  Trying to open as create
            ELSE IF( FSTATUS == FSCREA3 ) THEN

                !......  Check consistency of file set description
                IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                    MESG = 'Bad file description for file set ' // ROOTNAME16
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Close already open file
                IF( CLOSESET( ROOTNAME16 ) ) THEN
                    MESG = 'File ' // ROOTNAME16 // ' already opened. ' //&
                           'Closing, deleting, and re-opening it'
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                ELSE
                    MESG = 'File ' // ROOTNAME16 // ' already opened. ' //&
                           'Could not close to reopen with status FSCREA3'
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                    OPENSET = .FALSE.
                    RETURN
                END IF

            !......  Checked all problematic combinations
            ELSE

                OPENSET = .TRUE.
                RETURN

            END IF

        !......  Otherwise, file is not open
        ELSE
            IF( FSTATUS == FSNEW3 .OR. FSTATUS == FSUNKN3 .OR. FSTATUS == FSCREA3 ) THEN

                !......  Check consistency of file set description
                IF( .NOT. CHKSETDESC( ROOTNAME16 ) ) THEN
                    MESG = 'Bad file description for file set ' // ROOTNAME16
                    CALL M3WARN( FUNCNAME, 0, 0, MESG )
                    OPENSET = .FALSE.
                    RETURN
                END IF
            END IF
        END IF

        !......  Find a place for the new file
        FILEIDX = INDEX1( CMISS3, MXFILE3, RNAMES )
        IF( FILEIDX == 0 ) THEN
            MESG = 'Could not open ' // ROOTNAME16 //       &
                   'Maximum number of files already have been opened.'
            CALL M3WARN( FUNCNAME, 0, 0, MESG )
            OPENSET = .FALSE.
            RETURN
        END IF

        !......  Get env. variable setting for ROOTNAME
        CALL NAMEVAL( ROOTNAME16, ENVNAME )

        !......  Find out if file already exists on disk
        !..      First try original env. variable setting
        INQUIRE( FILE = ENVNAME, EXIST = FILEEXIST )

        !......  If original file doesn't exist, try first file of set
        IF( .NOT. FILEEXIST ) THEN
            CALL APPENDNAME( ENVNAME, 1, TEMPNAME, EFLAG )
            IF( EFLAG ) THEN
                MESG = 'Could not generate individual file names ' //   &
                       'for file set "' // TRIM( ROOTNAME ) // '"'
                CALL M3WARN( FUNCNAME, 0, 0, MESG )
                CALL CLEANUP( FILEIDX )
                OPENSET = .FALSE.
                RETURN
            END IF
            INQUIRE( FILE = TEMPNAME, EXIST = SETEXIST )
        END IF

        !......  Store logical file name and readonly status
        RNAMES( FILEIDX ) = ROOTNAME16
        FILE_INFO( FILEIDX )%RDONLY = ( FSTATUS .EQ. FSREAD3 )

        !......  Open file based on status and if file or set exists
        IF( FSTATUS == FSREAD3 .OR. FSTATUS == FSRDWR3 ) THEN

            !......  If single file exists, we don't have a file set
            IF( FILEEXIST ) THEN

                !......  Try to open single file
                IF( .NOT. OPEN3( ROOTNAME16, FSTATUS, FUNCNAME ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Get file description
                IF( .NOT. DESC3( ROOTNAME16 ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Allocate memory for logical file names and variables
                ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( 1 ),     &
                          FILE_INFO( FILEIDX )%VARS( NVARS3D,2 ), STAT=IOS )
                CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

                !......  Store logical file names and variables
                FILE_INFO( FILEIDX )%LNAMES( 1 ) = ROOTNAME16
                FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,1 ) = VNAME3D( 1:NVARS3D )
                FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,2 ) = ROOTNAME16

            !......  Otherwise, if the file set exists, open the individual files
            ELSE IF( SETEXIST ) THEN

                !......  Create new logical file name and set the env. variable
                LNAME = TRIM( ROOTNAME16 ) // '1'
                IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Try to open the first file in the set
                IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Get file description
                IF( .NOT. DESC3( LNAME ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Get total number of files from file header
                NFILESTR = GETCFDSC( FDESC3D, '/NUMBER OF FILES/', .TRUE. )
                NFILEINT = STR2INT( NFILESTR )
                IF( NFILEINT == IMISS3 ) THEN
                    MESG = 'Invalid number of files in header of file set "' // TRIM( ROOTNAME ) // '"'
                    CALL M3MSG2( MESG )
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Get total number of variables from file header
                NVARSTR  = GETCFDSC( FDESC3D, '/NUMBER OF VARIABLES/', .TRUE. )
                NVARINT  = STR2INT( NVARSTR )
                IF( NVARINT == IMISS3 ) THEN
                    MESG = 'Invalid number of variables in header of file set "' // TRIM( ROOTNAME ) // '"'
                    CALL M3MSG2( MESG )
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF

                !......  Allocate memory for logical file names and variables
                ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( NFILEINT ),     &
                          FILE_INFO( FILEIDX )%VARS( NVARINT,2 ),   STAT=IOS )
                CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

                !......  Store logical file names and variables
                FILE_INFO( FILEIDX )%LNAMES( 1 ) = LNAME
                FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,1 ) = VNAME3D( 1:NVARS3D )
                FILE_INFO( FILEIDX )%VARS( 1:NVARS3D,2 ) = LNAME
                VARPOS = NVARS3D + 1

                !......  Loop through remaining files in set
                DO I = 2, NFILEINT

                    !........  Create new logical and physical file names
                    WRITE( INTBUF, '(I2)' ) I
                    INTBUF = ADJUSTL( INTBUF )
                    LNAME = TRIM( ROOTNAME16 ) // TRIM( INTBUF )
                    CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                    IF( EFLAG ) THEN
                        CALL CLEANUP( FILEIDX )
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    !........  Set new env. variable
                    IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                        CALL CLEANUP( FILEIDX )
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    !........  Try to open file
                    IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                        CALL CLEANUP( FILEIDX )
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    !........  Get file description
                    IF( .NOT. DESC3( LNAME ) ) THEN
                        CALL CLEANUP( FILEIDX )
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    !........  Store logical file names and variables
                    FILE_INFO(FILEIDX)%LNAMES( I ) = LNAME
                    FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,1 ) = VNAME3D( 1:NVARS3D )
                    FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,2 ) = LNAME
                    VARPOS = VARPOS + NVARS3D

                END DO          ! loop over remaining files in file set

            ELSE
            !......  File doesn't exist on disk
                MESG = ROOTNAME16 // ':' // TRIM( ENVNAME )
                CALL M3MSG2( MESG )
                CALL M3WARN( FUNCNAME, 0, 0, 'File not available.' )

                CALL CLEANUP( FILEIDX )
                OPENSET = .FALSE.
                RETURN
            END IF

        ELSE IF( FSTATUS == FSNEW3 ) THEN

            !......  For status new, neither file should exist
            IF( FILEEXIST .OR. SETEXIST ) THEN
                MESG = ROOTNAME16 // ':' // TRIM( ENVNAME )
                CALL M3MSG2( MESG )
                MESG = 'File already exists on disk, cannot open as "NEW".'
                CALL M3WARN( FUNCNAME, 0, 0, MESG )
                CALL CLEANUP( FILEIDX )
                OPENSET = .FALSE.
                RETURN

            !......  Otherwise, try to create new files
            ELSE
                IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, FSTATUS ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF
            END IF

        ELSE IF( FSTATUS == FSUNKN3 ) THEN

            !......  If file or set exists, try to open it
            IF( FILEEXIST .OR. SETEXIST ) THEN

                !......  Allocate memory for logical file names and variables
                ALLOCATE( FILE_INFO( FILEIDX )%LNAMES( NFILESET ),     &
                          FILE_INFO( FILEIDX )%VARS( NVARSET,2 ),   STAT=IOS )
                CALL CHECKMEM( IOS, 'FILE_INFO%VARS', FUNCNAME )

                !......  Loop through individual files in the set
                VARPOS = 1
                DO I = 1, NFILESET

                    !........  Set necessary file description values
                    NVARS3D = VARS_PER_FILE( I )
                    VTYPE3D( 1:NVARS3D ) = VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
                    VNAME3D( 1:NVARS3D ) = VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )

                    !..........  Create new logical and physical file names if needed
                    IF( NFILESET > 1 ) THEN
                        WRITE( INTBUF, '(I2)' ) I
                        INTBUF = ADJUSTL( INTBUF )
                        LNAME = TRIM( ROOTNAME16 ) // TRIM( INTBUF )
                        CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                        IF( EFLAG ) THEN
                            CALL CLEANUP( FILEIDX )
                            OPENSET = .FALSE.
                            RETURN
                        END IF

                        !..........  Set new env. variable
                        IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                            CALL CLEANUP( FILEIDX )
                            OPENSET = .FALSE.
                            RETURN
                        END IF
                    ELSE
                        LNAME = ROOTNAME16
                    END IF

                    !........  Try to open file
                    IF( .NOT. OPEN3( LNAME, FSTATUS, FUNCNAME ) ) THEN
                        CALL CLEANUP( FILEIDX )
                        OPENSET = .FALSE.
                        RETURN
                    END IF

                    !........  Store logical file names and variables
                    FILE_INFO(FILEIDX)%LNAMES( I ) = LNAME
                    FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,1 ) = VNAME3D( 1:NVARS3D )
                    FILE_INFO(FILEIDX)%VARS( VARPOS:VARPOS+NVARS3D-1,2 ) = LNAME
                    VARPOS = VARPOS + NVARS3D

                END DO          ! loop over individual files

            !......  Otherwise, try to create new files
            ELSE
                IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, FSTATUS ) ) THEN
                    CALL CLEANUP( FILEIDX )
                    OPENSET = .FALSE.
                    RETURN
                END IF
            END IF

        ELSE IF( FSTATUS == FSCREA3 ) THEN

            !......  If file or set exists, delete it
            IF( FILEEXIST .OR. SETEXIST ) THEN

                !......  Loop through files in set
                DO I = 1, NFILESET

                    !........  If more than one file, generate subsequent file names
                    IF( NFILESET > 1 ) THEN
                        CALL APPENDNAME( ENVNAME, I, TEMPNAME, EFLAG )
                        IF( EFLAG ) THEN
                            CALL CLEANUP( FILEIDX )
                            OPENSET = .FALSE.
                            RETURN
                        END IF
                    ELSE
                        TEMPNAME = ENVNAME
                    END IF

                    !........  Check that file exists
                    INQUIRE( FILE = TEMPNAME, EXIST = TEMPEXIST )

                    !........  Try to delete files
                    IF( TEMPEXIST ) THEN
                        IOS = RMFILE( TEMPNAME )
                        IF( IOS /= 0 ) THEN
                            WRITE( MESG, 93010 ) 'Error number', IOS, 'removing file ' // TEMPNAME
                            CALL M3WARN( FUNCNAME, 0, 0, MESG )
                            CALL CLEANUP( FILEIDX )
                            OPENSET = .FALSE.
                            RETURN
                        END IF
                    END IF

                END DO
            END IF

                !......  Try to create new files
            IF( .NOT. CREATESET( FILEIDX, ROOTNAME16, ENVNAME, FSTATUS ) ) THEN
                CALL CLEANUP( FILEIDX )
                OPENSET = .FALSE.
                RETURN
            END IF

        !......  Illegal FSTATUS value
        ELSE

            MESG = 'File opening error: illegal FSTATUS argument.'
            CALL M3WARN( FUNCNAME, 0, 0, MESG )
            MESG = 'Legal values: 1-READONLY, 2-READ/WRITE, 3-NEW, 4 -UNKNOWN'
            CALL M3MSG2( MESG )
            WRITE( MESG, 93010 ) 'Value supplied by caller:', FSTATUS
            CALL M3MSG2( MESG )

            CALL CLEANUP( FILEIDX )
            OPENSET = .FALSE.
            RETURN

        END IF

        !......  Update number of open file sets
        NOPENSETS = NOPENSETS + 1

        OPENSET = .TRUE.

        RETURN

93010   FORMAT ( 5 ( A, :, I9, :, 2X ) )

    END FUNCTION OPENSET


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    LOGICAL FUNCTION CLOSESET( ROOTNAME )

        !......  Function arguments
        CHARACTER(*), INTENT(IN) :: ROOTNAME      ! logical file name for file set

        !......  Local variables
        INTEGER         I               ! counter
        INTEGER         FILEIDX         ! file index

        CHARACTER(16)   ROOTNAME16      ! fixed length root file name
        CHARACTER(256)  MESG            ! message buffer

        !---------------------------------
        !  Begin body of function CLOSESET
        !---------------------------------

        !......  Check length of file name
        IF( LEN_TRIM( ROOTNAME ) > 16 ) THEN
            MESG = 'Max file name length (16) exceeded for "' // TRIM( ROOTNAME ) // '"'
            CALL M3MSG2( MESG )
            CLOSESET = .FALSE.
            RETURN
        END IF

        !......  Get file index
        ROOTNAME16 = ROOTNAME
        FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )

        !......  If file is not open, exit with error
        IF( FILEIDX == 0 ) THEN
            MESG = 'File set "' // TRIM( ROOTNAME ) // '" is not currently open'
            CALL M3MSG2( MESG )
            CLOSESET = .FALSE.
            RETURN
        END IF

        !......  Loop through individual files
        DO I = 1, SIZE( FILE_INFO( FILEIDX )%LNAMES )
            IF( .NOT. CLOSE3( FILE_INFO( FILEIDX )%LNAMES( I ) ) ) THEN
                CLOSESET = .FALSE.
                RETURN
            END IF
        END DO

        !......  Write message only if file set contains more than one file
        IF( SIZE( FILE_INFO( FILEIDX )%LNAMES ) > 1 ) THEN
            MESG = 'Closing file set "' // TRIM( ROOTNAME ) // '"'
            CALL M3MSG2( MESG )
        END IF

        CALL CLEANUP( FILEIDX )
        NOPENSETS = NOPENSETS - 1

        CLOSESET = .TRUE.
        RETURN

    END FUNCTION CLOSESET


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    LOGICAL FUNCTION DESCSET( ROOTNAME, FILENUM )

        !......  Function arguments
        CHARACTER(*), INTENT(IN) :: ROOTNAME      ! logical file name for file set
        INTEGER,      INTENT(IN) :: FILENUM       ! position of file in set or ALLFILES

        !......  External functions
        CHARACTER(10), EXTERNAL :: GETCFDSC

        CHARACTER(24), PARAMETER :: FUNCNAME = 'MODFILESET/DESCSET'      ! function name

        !......  Local variables
        INTEGER            I               ! counter
        INTEGER            IOS             ! I/O status
        INTEGER            FILEIDX         ! file index
        INTEGER            VARPOS          ! position in variable information arrays
        INTEGER            NFILEINT        ! number of files as an integer
        INTEGER            NVARINT         ! number of variables as an integer

        CHARACTER(2)   NFILESTR            ! number of files as a string
        CHARACTER(4)   NVARSTR             ! number of variables as a string
        CHARACTER(16)  ROOTNAME16          ! fixed length root file name
        CHARACTER(16)  LNAME               ! temporary logical name
        CHARACTER(256) MESG                ! message buffer

        !---------------------------------
        !  Begin body of function DESCSET
        !---------------------------------

        !......  Check length of file name
        IF( LEN_TRIM( ROOTNAME ) > 16 ) THEN
            MESG = 'Max file name length (16) exceeded for "' // ROOTNAME // '"'
            CALL M3MSG2( MESG )
            DESCSET = .FALSE.
            RETURN
        END IF

        !......  Get file index
        ROOTNAME16 = ROOTNAME
        FILEIDX = INDEX1( ROOTNAME16, MXFILE3, RNAMES )

        !......  If file is not open, exit with error
        IF( FILEIDX == 0 ) THEN
            MESG = 'File set "' // TRIM( ROOTNAME ) // '" is not currently open'
            CALL M3MSG2( MESG )
            DESCSET = .FALSE.
            RETURN
        END IF

        !......  Deallocate variable information arrays if necessary
        IF( ALLOCATED( VARS_PER_FILE ) ) THEN
            DEALLOCATE( VARS_PER_FILE, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "VARS_PER_FILE" variable'
                CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( ALLOCATED( VTYPESET ) ) THEN
            DEALLOCATE( VTYPESET, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "VTYPESET" variable'
                CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( ALLOCATED( VNAMESET ) ) THEN
            DEALLOCATE( VNAMESET, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "VNAMESET" variable'
                CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( ALLOCATED( VUNITSET ) ) THEN
            DEALLOCATE( VUNITSET, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "VUNITSET" variable'
                CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        IF( ALLOCATED( VDESCSET ) ) THEN
            DEALLOCATE( VDESCSET, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "VDESCSET" variable'
                CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
            END IF
        END IF

        !......  Check if a single file is requested
        IF( FILENUM /= ALLFILES ) THEN

            !......  Make sure file number is positive
            IF( FILENUM < 1 ) THEN
                MESG = 'File number must be positive or "ALLFILES"'
                CALL M3MSG2( MESG )
                DESCSET = .FALSE.
                RETURN

            !......  Check that number is not greater than total number of files
            ELSE IF( FILENUM > SIZE( FILE_INFO( FILEIDX )%LNAMES ) ) THEN
                WRITE( MESG,92010 ) 'Invalid file number requested; ' //    &
                       'file set "' // TRIM( ROOTNAME ) // '" contains ',   &
                       SIZE( FILE_INFO( FILEIDX )%LNAMES ), ' files '
                CALL M3MSG2( MESG )
                DESCSET = .FALSE.
                RETURN
            ELSE

                !......  Set logical file name
                LNAME = FILE_INFO( FILEIDX )%LNAMES( FILENUM )

                !......  Try to get file description
                IF( .NOT. DESC3( LNAME ) ) THEN
                    DESCSET = .FALSE.
                    RETURN
                END IF

                !......  Store file set information
                NFILESET = 1
                NVARSET = NVARS3D

                !......  Allocate variable information arrays
                ALLOCATE( VARS_PER_FILE( NFILESET ),    &
                          VTYPESET( NVARSET ),          &
                          VNAMESET( NVARSET ),          &
                          VUNITSET( NVARSET ),          &
                          VDESCSET( NVARSET ), STAT=IOS )
                CALL CHECKMEM( IOS, 'VARS_PER_FILE...VDESCSET', FUNCNAME )

                !......  Store info for this file
                VARS_PER_FILE( 1 ) = NVARS3D
                VTYPESET( 1:NVARS3D ) = VTYPE3D( 1:NVARS3D )
                VNAMESET( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
                VUNITSET( 1:NVARS3D ) = UNITS3D( 1:NVARS3D )
                VDESCSET( 1:NVARS3D ) = VDESC3D( 1:NVARS3D )
            END IF
        ELSE

            !......  Open first file of set
            LNAME = FILE_INFO( FILEIDX )%LNAMES( 1 )
            IF( .NOT. DESC3( LNAME ) ) THEN
                DESCSET = .FALSE.
                RETURN
            END IF

            !......  Get total number of files from file header
            NFILESTR = GETCFDSC( FDESC3D, '/NUMBER OF FILES/', .TRUE. )
            NFILEINT = STR2INT( NFILESTR )
            IF( NFILEINT == IMISS3 ) THEN
                MESG = 'Invalid number of files in header of file set "' // TRIM( ROOTNAME ) // '"'
                CALL M3MSG2( MESG )
                DESCSET = .FALSE.
                RETURN
            END IF

            !......  Get total number of variables from file header
            NVARSTR  = GETCFDSC( FDESC3D, '/NUMBER OF VARIABLES/', .TRUE. )
            NVARINT  = STR2INT( NVARSTR )
            IF( NVARINT == IMISS3 ) THEN
                MESG = 'Invalid number of variables in header of file set "' // TRIM( ROOTNAME ) // '"'
                CALL M3MSG2( MESG )
                DESCSET = .FALSE.
                RETURN
            END IF

            !......  Store file set information from header
            NFILESET = NFILEINT
            NVARSET = NVARINT

            !......  Allocate variable information arrays
            ALLOCATE( VARS_PER_FILE( NFILESET ),    &
                      VTYPESET( NVARSET ),          &
                      VNAMESET( NVARSET ),          &
                      VUNITSET( NVARSET ),          &
                      VDESCSET( NVARSET ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VARS_PER_FILE...VDESCSET', FUNCNAME )

            !......  Store info for this file
            VARS_PER_FILE( 1 ) = NVARS3D
            VTYPESET( 1:NVARS3D ) = VTYPE3D( 1:NVARS3D )
            VNAMESET( 1:NVARS3D ) = VNAME3D( 1:NVARS3D )
            VUNITSET( 1:NVARS3D ) = UNITS3D( 1:NVARS3D )
            VDESCSET( 1:NVARS3D ) = VDESC3D( 1:NVARS3D )
            VARPOS = NVARS3D + 1

            !......  Loop through remaining files
            DO I = 2, NFILESET

                !......  Set logical file name
                LNAME = FILE_INFO( FILEIDX )%LNAMES( I )

                !......  Try to get file description
                IF( .NOT. DESC3( LNAME ) ) THEN
                    DESCSET = .FALSE.
                    RETURN
                END IF

                !......  Store current file description
                VARS_PER_FILE( I ) = NVARS3D
                VTYPESET( VARPOS:VARPOS + NVARS3D - 1 ) = VTYPE3D( 1:NVARS3D )
                VNAMESET( VARPOS:VARPOS + NVARS3D - 1 ) = VNAME3D( 1:NVARS3D )
                VUNITSET( VARPOS:VARPOS + NVARS3D - 1 ) = UNITS3D( 1:NVARS3D )
                VDESCSET( VARPOS:VARPOS + NVARS3D - 1 ) = VDESC3D( 1:NVARS3D )
                VARPOS = VARPOS + NVARS3D
            END DO

        END IF

        DESCSET = .TRUE.

        !---------- Format statements --------------

92010   FORMAT ( A, :, I3, :, A )

    END FUNCTION DESCSET


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    CHARACTER(16) FUNCTION PROMPTSET( PROMPT, FMODE, DEFAULT, CALLER )

        !......  Function arguments
        CHARACTER(*), INTENT(IN) :: PROMPT       ! prompt for user
        INTEGER,      INTENT(IN) :: FMODE        ! file opening mode
        CHARACTER(*), INTENT(IN) :: DEFAULT      ! default logical file name
        CHARACTER(*), INTENT(IN) :: CALLER       ! name of calling program

        !......  External functions

        !......  Local parameters
        CHARACTER(16), PARAMETER :: NONE16   = 'NONE'
        CHARACTER(24), PARAMETER :: FUNCNAME = 'MODFILESET/PROMPTSET'      ! function name

        !......  Local variables
        INTEGER            IOS                    ! I/O status
        INTEGER            IDX                    ! index in logical file name

        LOGICAL, SAVE ::   PROMPTON               ! true: prompt for input
        LOGICAL, SAVE ::   INITIAL = .TRUE.       ! true: first time
        LOGICAL            NFLAG                  ! true: "NONE" is in the prompt

        CHARACTER(16)  LNAME                      ! logical file name
        CHARACTER(300) MESG                       ! message buffer
        CHARACTER(512) BUFFER                     ! prompt buffer

        !------------------------------------
        !  Begin body of function PROMPTSET
        !------------------------------------

        !......  On first time, check if prompt should be shown
        IF( INITIAL ) THEN
            PROMPTON = ENVYN( 'PROMPTFLAG', 'Prompt for input flag', .TRUE., IOS )
            INITIAL = .FALSE.
        END IF

        !......  Decide if 'NONE' is a valid response
        NFLAG = ( INDEX( PROMPT, '"NONE"' ) > 0 )

        IF( PROMPTON ) THEN

            !......  Construct actual prompt
            BUFFER = TRIM( PROMPT ) // ' [' // TRIM( DEFAULT ) // '] >>'

            !......  Loop until valid name is given and file set is opened
            DO
                LNAME = ' '
                WRITE( *,95000 ) TRIM( BUFFER ) // ' '
                READ ( *,'(A16)',IOSTAT=IOS ) LNAME

                !......  Check if response was read
                IF( IOS > 0 ) THEN
                    MESG = 'Could not read your response'
                    CALL M3MSG2( MESG )
                    IF( GETYN( 'Try again?', .TRUE. ) ) THEN
                        CYCLE
                    ELSE
                        MESG = 'Could not read logical name for file set'
                        CALL M3EXIT( FUNCNAME, 0, 0, MESG, 2 )
                    END IF
                ELSE IF( IOS < 0 ) THEN
                    MESG = 'Ending program "' // TRIM( CALLER ) // '".'
                    CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                END IF

                !......  Check logical name
                IDX = INDEX( LNAME, '!' )
                IF( IDX > 0 ) THEN
                    LNAME( IDX:LEN( LNAME ) ) = ' '
                END IF

                IF( LNAME == ' ' ) THEN
                    LNAME = DEFAULT
                END IF

                IF( NFLAG .AND. ( LNAME == NONE16 ) ) THEN
                    PROMPTSET = NONE16
                    RETURN
                END IF

                !......  Try to open file set
                IF( .NOT. OPENSET( LNAME, FMODE, CALLER ) ) THEN
                    MESG = 'Could not open file set "' // TRIM( LNAME ) // '".'
                    CALL M3MSG2( MESG )
                    IF( GETYN( 'Try again?', .TRUE. ) ) THEN
                        CYCLE
                    ELSE
                        MESG = 'Ending program "' // TRIM( CALLER ) // '".'
                        CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
                    END IF
                ELSE
                    EXIT
                END IF
            END DO

        !......  Otherwise, don't prompt for input
        ELSE
            LNAME = DEFAULT

            !......  Check if NONE is valid
            IF( NFLAG ) THEN
                IF( LNAME == NONE16 ) THEN
                    PROMPTSET = NONE16
                    RETURN
                END IF

                !......  Check if logical file name is set
                CALL ENVSTR( LNAME, 'Input file name', ' ', BUFFER, IOS )

                !......  If not set (-2) or empty (-1)
                IF( IOS < 0 ) THEN
                    PROMPTSET = NONE16
                    RETURN
                END IF
            END IF

            !......  Try to open file set
            IF( .NOT. OPENSET( LNAME, FMODE, CALLER ) ) THEN
                MESG = 'Could not open file set "' // TRIM( LNAME ) // '".'
                CALL M3MSG2( MESG )
                CALL M3EXIT( CALLER, 0, 0, MESG, 2 )
            END IF
        END IF

        PROMPTSET = LNAME
        RETURN

        !---------- Format statements --------------

    94000 FORMAT( 12( A, : ) )

    95000 FORMAT (/5X, A )              !  generic prompt format.

    END FUNCTION PROMPTSET


    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-
    !!-=-=-=-=-=-=-=-  PRIVATE ROUTINES  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    !!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-


    LOGICAL FUNCTION CREATESET( FIDX, RNAME, PHYSNAME, FSTATUS )

        !......  Subroutine arguments
        INTEGER,        INTENT(IN) :: FIDX          ! file index
        CHARACTER(16),  INTENT(IN) :: RNAME         ! root logical file name
        CHARACTER(256), INTENT(IN) :: PHYSNAME      ! physical file name
        INTEGER,        INTENT(IN) :: FSTATUS       ! file mode (read, read/write)

        !......  Local variables
        INTEGER I, J, K                  ! counters
        INTEGER IOS                      ! I/O status
        INTEGER VARPOS                   ! variable position
        INTEGER BLANKS                   ! number of blank spaces in FDESC3D
        INTEGER FDESCPOS                 ! position in FDESC3D array

        LOGICAL EFLAG                    ! true: error has occurred

        CHARACTER(2)   INTBUF            ! integer string buffer
        CHARACTER(16)  LNAME             ! temporary logical file name
        CHARACTER(256) TEMPNAME          ! temporary physical file name
        CHARACTER(300) MESG              ! message buffer

        CHARACTER(24), PARAMETER :: PROGNAME = 'MODFILESET/CREATESET'      ! program name

        !-------------------------------------
        !  Begin body of function CREATESET
        !-------------------------------------

        !......  Allocate memory for logical file names and variables
        ALLOCATE( FILE_INFO( FIDX )%LNAMES( NFILESET ),     &
                  FILE_INFO( FIDX )%VARS( NVARSET,2 ), STAT=IOS )
        CALL CHECKMEM( IOS, 'FILE_INFO%VARS', PROGNAME )

        !......  Loop through individual files in the set
        VARPOS = 1
        DO I = 1, NFILESET

            !......  Set file description values
            NVARS3D = VARS_PER_FILE( I )
            VTYPE3D( 1:NVARS3D ) = VTYPESET( VARPOS:VARPOS + NVARS3D - 1 )
            VNAME3D( 1:NVARS3D ) = VNAMESET( VARPOS:VARPOS + NVARS3D - 1 )
            UNITS3D( 1:NVARS3D ) = VUNITSET( VARPOS:VARPOS + NVARS3D - 1 )
            VDESC3D( 1:NVARS3D ) = VDESCSET( VARPOS:VARPOS + NVARS3D - 1 )

            !......  On first time through, find space in FDESC3D array
            IF( I == 1 ) THEN
                BLANKS = 0
                DO J = 1, MXDESC3
                    IF( FDESC3D( J ) /= '' ) THEN
                        BLANKS = 0
                    ELSE
                        BLANKS = BLANKS + 1
                    END IF

                    !........  If we have 3 consecutive open spots, add file info
                    IF( BLANKS == 3 ) THEN
                        WRITE( FDESC3D( J - 2 ),93010 ) '/NUMBER OF FILES/', NFILESET
                        WRITE( FDESC3D( J - 1 ),93010 ) '/FILE POSITION/', I
                        WRITE( FDESC3D( J ),    93010 ) '/NUMBER OF VARIABLES/', NVARSET

                        !..........  Save position for writing subsequent info
                        FDESCPOS = J - 1
                        EXIT
                    END IF

                    !........  Couldn't find space
                    IF( J == MXDESC3 ) THEN
                        MESG = 'No spaces available in FDESC3D array'
                        CALL M3WARN( PROGNAME, 0, 0, MESG )
                        CREATESET = .FALSE.
                        RETURN
                    END IF
                END DO
            ELSE
            !......  Add current file position number
                WRITE( FDESC3D( FDESCPOS ),93010 )  '/FILE POSITION/', I
            END IF

            !......  Create new logical and physical file names if needed
            IF( NFILESET > 1 ) THEN
                WRITE( INTBUF, '(I2)' ) I
                INTBUF = ADJUSTL( INTBUF )
                LNAME = TRIM( RNAME ) // TRIM( INTBUF )
                CALL APPENDNAME( PHYSNAME, I, TEMPNAME, EFLAG )
                IF( EFLAG ) THEN
                    CREATESET = .FALSE.
                    RETURN
                END IF

                !......  Set new env. variable
                IF( .NOT. SETENVVAR( LNAME, TEMPNAME ) ) THEN
                    CREATESET = .FALSE.
                    RETURN
                END IF
            ELSE
                LNAME = RNAME
            END IF

            !......  Try to open file
            IF( .NOT. OPEN3( LNAME, FSTATUS, PROGNAME ) ) THEN
                CREATESET = .FALSE.
                RETURN
            END IF

            !......  Store logical file names and variables
            FILE_INFO( FIDX )%LNAMES( I ) = LNAME
            FILE_INFO( FIDX )%VARS( VARPOS:VARPOS+NVARS3D-1,1 ) = VNAME3D( 1:NVARS3D )
            FILE_INFO( FIDX )%VARS( VARPOS:VARPOS+NVARS3D-1,2 ) = LNAME
            VARPOS = VARPOS + NVARS3D

        END DO      ! loop over individual files

        !......  Reset FDESC3D array back to original
        !..        (remove file number and position info)
        FDESC3D( FDESCPOS-1:FDESCPOS+1 ) = ''

        CREATESET = .TRUE.
        RETURN

        !---------- Format statements --------------

93010   FORMAT ( A, I4 )

    END FUNCTION CREATESET


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    LOGICAL FUNCTION CHKFILESET( FIDX )

        !......  Function arguments
        INTEGER, INTENT(IN) :: FIDX      ! file index

        !......  Local variables
        INTEGER TOTALFILES      ! total number of files
        INTEGER TOTALVARS       ! total number of variables

        CHARACTER(16)  FNAME            ! logical file name
        CHARACTER(300) MESG             ! message buffer

        LOGICAL :: EFLAG = .FALSE.      ! true: error found

        !------------------------------------
        !  Begin body of function CHKFILESET
        !------------------------------------

        EFLAG = .FALSE.
        FNAME = RNAMES( FIDX )

        !......  Check that total number of files is consistent
        TOTALFILES = SIZE( FILE_INFO( FIDX )%LNAMES )
        IF( TOTALFILES /= NFILESET ) THEN
            WRITE( MESG,91010 )
            CALL M3MSG2( MESG )
            WRITE( MESG,91020 ) 'Inconsistent number of files ' //  &
                                'for file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            WRITE( MESG,91030 ) 'Value from file:  ', TOTALFILES
            CALL M3MSG2( MESG )
            WRITE( MESG,91030 ) 'Value from caller:', NFILESET
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        !......  Check that total number of vars is consistent
        TOTALVARS = SIZE( FILE_INFO( FIDX )%VARS,1 )
        IF( TOTALVARS /= NVARSET ) THEN
            WRITE( MESG,91010 )
            CALL M3MSG2( MESG )
            WRITE( MESG,91020 ) 'Inconsistent number of variables ' //  &
                                'for file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            WRITE( MESG,91030 ) 'Value from file:  ', TOTALVARS
            CALL M3MSG2( MESG )
            WRITE( MESG,91030 ) 'Value from caller:', NVARSET
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        !......  Set function value and return
        IF( EFLAG ) THEN
            CHKFILESET = .FALSE.
        ELSE
            CHKFILESET = .TRUE.
        END IF

        RETURN

    !---------- Format statements --------------

91010   FORMAT ( 5X, '>>> WARNING in function CHKFILESET <<<' )

91020   FORMAT ( 5X , A )

91030   FORMAT ( 5X , A , I4 )

    END FUNCTION CHKFILESET



    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    LOGICAL FUNCTION CHKSETDESC( FNAME )

        !......  Function arguments
        CHARACTER(16), INTENT(IN) :: FNAME      ! logical file name

        !......  Local variables
        INTEGER              I                            ! counter
        INTEGER              IOS                          ! I/O status
        INTEGER              TOTALVARS                    ! total number of variables

        REAL                 MXVARS                       ! real MXVARS3 from ioapi param

        LOGICAL ::           EFLAG = .FALSE.              ! true: error found

        CHARACTER(300)   MESG                         ! message buffer

        CHARACTER(16) :: FUNCNAME = 'CHKSETDESC'      ! function name

        !------------------------------------
        !  Begin body of function CHKSETDESC
        !------------------------------------

        EFLAG = .FALSE.
        MXVARS = MXVARS3

        !......  Check that total number of variables is valid
        IF( NVARSET <= 0 ) THEN
            MESG = 'Illegal number of variable for file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            CHKSETDESC = .FALSE.
            RETURN
        END IF

        !......  Check that number of variables per file is set
        IF( .NOT. ALLOCATED( VARS_PER_FILE ) ) THEN
            MESG = 'Number of variables per file array is not ' //  &
                   'allocated for file set ' // TRIM( FNAME ) // ';'
            CALL M3MSG2( MESG )
            WRITE( MESG,'( 10(A,I5) )' ) 'using default of ',MXVARS3,' variables per file'
            CALL M3MSG2( MESG )

            NFILESET = CEILING( NVARSET / MXVARS )      ! real division, convert back to int
            ALLOCATE( VARS_PER_FILE( NFILESET ), STAT=IOS )
            CALL CHECKMEM( IOS, 'VARS_PER_FILE', FUNCNAME )

            TOTALVARS = NVARSET
            I = 1
            DO
                IF( TOTALVARS - MXVARS3 <= 0 ) THEN
                    VARS_PER_FILE( I ) = TOTALVARS
                    EXIT
                ELSE
                    VARS_PER_FILE( I ) = MXVARS3
                END IF

                TOTALVARS = TOTALVARS - MXVARS3
                I = I + 1
            END DO

        END IF

        !......  Check that total number of files is consistent
        IF( SIZE( VARS_PER_FILE ) /= NFILESET ) THEN
            MESG = 'Inconsistent number of files in file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        !......  Count total number of variables in VARS_PER_FILE;
        !..        also check that number of variables per file is not more than MXVARS3
        TOTALVARS = 0
        DO I = 1, NFILESET
            TOTALVARS = TOTALVARS + VARS_PER_FILE( I )
            IF( VARS_PER_FILE( I ) > MXVARS3 ) THEN
                WRITE( MESG,'( 10(A,I5) )' )        &
                    'More than',MXVARS3,' variables per file in file set ' // TRIM( FNAME )
                CALL M3MSG2( MESG )
                EFLAG = .TRUE.
            END IF
        END DO

        !......  Check that total number of vars is consistent
        IF( TOTALVARS /= NVARSET ) THEN
            MESG = 'Inconsistent number of variables in file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        !......  Make sure variable information is set
        IF( .NOT. ALLOCATED( VTYPESET ) .OR.    &
            .NOT. ALLOCATED( VNAMESET ) .OR.    &
            .NOT. ALLOCATED( VUNITSET ) .OR.    &
            .NOT. ALLOCATED( VDESCSET ) ) THEN
            MESG = 'One or more variable information arrays are ' //    &
                   'not allocated for file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            CHKSETDESC = .FALSE.
            RETURN
        END IF

        !......  Check size of variable information arrays
        IF( SIZE( VTYPESET ) < NVARSET .OR.     &
            SIZE( VNAMESET ) < NVARSET .OR.     &
            SIZE( VUNITSET ) < NVARSET .OR.     &
            SIZE( VDESCSET ) < NVARSET ) THEN
            MESG = 'One or more variable information arrays are ' //    &
                   'not as large as the number of variables in ' //     &
                   'file set ' // TRIM( FNAME )
            CALL M3MSG2( MESG )
            EFLAG = .TRUE.
        END IF

        !......  Set function value and return
        IF( EFLAG ) THEN
            CHKSETDESC = .FALSE.
        ELSE
            CHKSETDESC = .TRUE.
        END IF

        RETURN

    END FUNCTION CHKSETDESC


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    SUBROUTINE APPENDNAME( ORIGNAME, FILENUM, NEWNAME, EFLAG )

        !......  Subroutine arguments
        CHARACTER(*), INTENT(IN)  :: ORIGNAME      ! original file name
        INTEGER,      INTENT(IN)  :: FILENUM       ! number to append to name
        CHARACTER(*), INTENT(OUT) :: NEWNAME       ! new file name
        LOGICAL,      INTENT(OUT) :: EFLAG         ! true if error occurs

        !......  Local arguments
        INTEGER          IDX         ! string index
        INTEGER          INTLEN      ! length of integer string
        INTEGER          ORGLEN      ! length of original name
        INTEGER          NEWLEN      ! length of new name
        INTEGER          TMPLEN      ! length of temporary name
        CHARACTER(2)     INTBUF      ! buffer for integer value

        !--------------------------------------
        !  Begin body of subroutine APPENDNAME
        !--------------------------------------

        EFLAG = .FALSE.

        !......  Find .ncf extension by looking from the end of the string
        IDX = INDEX( ORIGNAME, '.', .TRUE. )

        !......  Check that a location was found
        IF( IDX == 0 ) THEN
            EFLAG = .TRUE.
            RETURN
        END IF

        !......  Make sure NEWNAME is as long as ORIGNAME
        IF( LEN( NEWNAME ) < LEN( ORIGNAME ) ) THEN
            EFLAG = .TRUE.
            RETURN
        END IF

        !......  Check that FILENUM is two characters or less
        IF( FILENUM < 1 .OR. FILENUM >= 100 ) THEN
            EFLAG = .TRUE.
            RETURN
        END IF

        !......  Write integer to string
        WRITE( INTBUF, '(I2)' ) FILENUM
        INTBUF = ADJUSTL( INTBUF )
        INTLEN = LEN_TRIM( INTBUF )

        !......  Get length of names
        ORGLEN = LEN( ORIGNAME )
        NEWLEN = LEN( NEWNAME )

        !......  Create new name from original name
        NEWNAME( 1:IDX-1 ) = ORIGNAME( 1:IDX-1 )
        NEWNAME( IDX:IDX+INTLEN ) = '.' // TRIM( INTBUF )

        TMPLEN = IDX+INTLEN
        NEWNAME( TMPLEN+1:NEWLEN ) = ORIGNAME( IDX:ORGLEN-INTLEN-1 )

        RETURN

    END SUBROUTINE APPENDNAME


    !!---------------------------------------------------------------
    !!---------------------------------------------------------------

    SUBROUTINE CLEANUP( FIDX )

        !......  Subroutine arguments
        INTEGER, INTENT(IN) :: FIDX      ! file index

        !......  Local variables
        INTEGER  IOS      ! I/O status

        CHARACTER(300) MESG                        ! message buffer

        CHARACTER(24), PARAMETER :: PROGNAME = 'MODFILESET/CLEANUP'      ! program name

        !-----------------------------------
        !  Begin body of subroutine CLEANUP
        !-----------------------------------

        !......  Set root name back to missing
        RNAMES( FIDX ) = CMISS3

        !......  Check if logical names pointer is associated, then deallocate
        IF( ASSOCIATED( FILE_INFO( FIDX )%LNAMES ) ) THEN
            DEALLOCATE( FILE_INFO( FIDX )%LNAMES, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "FILE_INFO%LNAMES" variable'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            NULLIFY( FILE_INFO( FIDX )%LNAMES )
        END IF

        !......  Check if variable pointer is associated, then deallocate
        IF( ASSOCIATED( FILE_INFO( FIDX )%VARS ) ) THEN
            DEALLOCATE( FILE_INFO( FIDX )%VARS, STAT=IOS )
            IF( IOS > 0 ) THEN
                MESG = 'Failure deallocating memory for "FILE_INFO%VARS" variable'
                CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
            END IF
            NULLIFY( FILE_INFO( FIDX )%VARS )
        END IF

        RETURN

    END SUBROUTINE CLEANUP


END MODULE MODFILESET

