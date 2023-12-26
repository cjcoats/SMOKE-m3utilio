
INTEGER FUNCTION IOAPI_GRD_SIZE( NCOLS, NROWS, NLAYS, NVARS, NSTEPS )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !      This function returns the approximate size in megabytes of a
    !      gridded I/O API file based on the number of columns, rows,
    !      layers, variables, and time steps.
    !
    !  PRECONDITIONS REQUIRED:

    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created by C. Seppanen 3/03
    !
    !       Version 10/2016 by C. Coats:
    !            Update for I/O API 3.1 or later;
    !            USE M3UTILIO
    !       Version 11/2023 by CJC:  conversion to ".f90",  and
    !       related changes
    !
    !**************************************************************************
    !
    ! Project Title: EDSS Tools Library
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

    IMPLICIT NONE

    !....  Function arguments
    INTEGER, INTENT (IN) :: NCOLS       ! number of columns
    INTEGER, INTENT (IN) :: NROWS       ! number of rows
    INTEGER, INTENT (IN) :: NLAYS       ! number of layers
    INTEGER, INTENT (IN) :: NVARS       ! number of variables
    INTEGER, INTENT (IN) :: NSTEPS      ! number of time steps

    !....  Other local variables
    INTEGER              NCELLS         ! number of grid cells
    INTEGER              HDRSIZE        ! size of header in bytes
    INTEGER              RECSIZE        ! size of single record in bytes

    !***********************************************************************
    !   begin body of function IOAPI_GRD_SIZE

    !....  Calculate number of grid cells
    NCELLS = NCOLS * NROWS

    !....  Calculate size of header
    HDRSIZE = 9860 + 4 * MXLAYS3 + 116 * MXVARS3 + 160 * MXDESC3

    !....  Calculate size of individual records
    RECSIZE = ( 8 * NVARS + 4 * NLAYS * NVARS * NCELLS )

    !....  Calculate total number of bytes in file
    IOAPI_GRD_SIZE = ( HDRSIZE + RECSIZE * NSTEPS ) / 1000000

    RETURN

END FUNCTION IOAPI_GRD_SIZE
