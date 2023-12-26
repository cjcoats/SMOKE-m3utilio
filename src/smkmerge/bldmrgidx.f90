
SUBROUTINE BLDMRGIDX( MXGRP, MXVARPGP, NGRP )

    !***********************************************************************
    !  subroutine BLDMRGIDX body starts at line
    !
    !  DESCRIPTION:
    !      The purpose of this subroutine is to allocate and populate indicator
    !      arrays that say which pollutants and species are present for each
    !      type of input file (inventory, speciation matrix, multiplicative
    !      control matrix, reactivity matrix, etc.) for each source category
    !      (area, biogenic, mobile, point).  The indicator arrays store the
    !      position in the list of i/o api file variables that the given pollutant
    !      or species matches.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created 2/99 by M. Houyoux
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
    !       related changes
    !***********************************************************************
    !
    ! Project Title: Sparse Matrix Operator Kernel Emissions (SMOKE) Modeling
    !       System
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

    !.....  MODULES for public variables
    !.....  This module contains the major data structure and control flags
    USE MODMERGE, ONLY: AFLAG, BFLAG, MFLAG, PFLAG,                     &
                        AUFLAG, MUFLAG, PUFLAG,                         &
                        ARFLAG, MRFLAG, PRFLAG, SFLAG,                  &
                        ANIPOL, PNIPOL, MNIPPA, NIPPA,                  &
                        ANSMATV, BNSMATV, MNSMATV, PNSMATV, NSMATV,     &
                        ANUMATV, MNUMATV, PNUMATV,                      &
                        ANRMATV, MNRMATV, PNRMATV,                      &
                        ARNMSPC, MRNMSPC, PRNMSPC, NMSPC,               &
                        AEINAM, MEANAM, PEINAM,                         &
                        ASVDESC, BSVDESC, MSVDESC, PSVDESC, TSVDESC,    &
                        ARVDESC, MRVDESC, PRVDESC,                      &
                        AUVNAMS, MUVNAMS, PUVNAMS,                      &
                        SIINDEX, SPINDEX,                               &
                        A_EXIST, M_EXIST, P_EXIST,                      &
                        AU_EXIST, MU_EXIST, PU_EXIST,                   &
                        AR_EXIST, MR_EXIST, PR_EXIST,                   &
                        AS_EXIST, BS_EXIST, MS_EXIST, PS_EXIST,         &
                        EMNAM, EANAM, VGRPCNT, IDVGP, GVNAMES, GVLOUT

    IMPLICIT NONE

    !.....   INCLUDES:

    INCLUDE 'EMCNST3.h90'       !  emissions constant parameters

    !.....  SUBROUTINE ARGUMENTS

    INTEGER, INTENT  (IN) :: MXGRP    ! Maximum no. groups
    INTEGER, INTENT  (IN) :: MXVARPGP     ! Max no. pollutants per group
    INTEGER, INTENT (OUT) :: NGRP     ! Actual number of groups

    !.....   Local allocatable arrays
    INTEGER, ALLOCATABLE :: PPSCNT( : )      ! count of pols per spc or emtype
    INTEGER, ALLOCATABLE :: KAU   ( : )      ! helper array for area mult matrix indices
    INTEGER, ALLOCATABLE :: KMU   ( : )      ! helper array for mobile mult matrix indices
    INTEGER, ALLOCATABLE :: KPU   ( : )      ! helper array for point mult matrix indices

    !.....   Call allocated arrays
    !.....   Group index counter for each source-category-specific list of
    !       pollutants and activities.
    INTEGER  KA( NIPPA )    !  area
    INTEGER  KM( NIPPA )    !  mobile
    INTEGER  KP( NIPPA )    !  point

    !.....   Other local variables
    INTEGER         J, K, L1, L2, N, V    !  counters and indices

    INTEGER         ACNT     ! area src var counter
    INTEGER         BCNT     ! biogenic src var counter
    INTEGER         GCNT     ! group counter
    INTEGER         IOS      ! i/o error status
    INTEGER         JA, JM, JP     ! position in rctvty var desc of pol-spc
    INTEGER         MAXPPS       ! max no. pols per species or per emis type
    INTEGER         MCNT     ! mobile src var counter
    INTEGER         NLOOP    ! tmp loop total
    INTEGER         NSIZE    ! tmp size of PPSCNT
    INTEGER         PCNT     ! point src var counter
    INTEGER         PGCNT    ! previous iteration group count
    INTEGER         PJ       ! previous iteration J index
    INTEGER         PVCNT    ! previous iteration variable count
    INTEGER         VCNT     ! variable counter
    INTEGER         TGRP     ! tmp group number

    LOGICAL      :: EFLAG = .FALSE.      ! true: error found
    LOGICAL         NEXTGRP      ! true: time to increment group number cntr

    CHARACTER(300)     :: MESG     ! message buffer
    CHARACTER(MXDLEN3) :: CBUF     ! tmp pol-to-species buffer
    CHARACTER(NAMLEN3) :: CPOL     ! tmp pollutant buffer
    CHARACTER(NAMLEN3) :: CSPC     ! tmp species buffer
    CHARACTER(NAMLEN3) :: PSPC     ! tmp previous species
    CHARACTER(NAMLEN3) :: PPOL     ! tmp previous pollutant
    CHARACTER(NAMLEN3) :: VBUF     ! tmp variable name buffer
    CHARACTER(PLSLEN3) :: SVBUF     ! tmp speciation name buffer

    CHARACTER(16), PARAMETER :: PROGNAME = 'BLDMRGIDX'     ! program name

    !***********************************************************************
    !   begin body of subroutine BLDMRGIDX

    !.....  Ensure that the max no. of variables per group is large
    !       enough so that all pol-to-species combos for the same species
    !       are in one group...

    !.....  For speciation or not, set number of iterations for upcoming loop
    !       and the dimension for the counter
    IF ( SFLAG ) THEN
        NLOOP = NSMATV
        NSIZE = NMSPC
    ELSE
        NLOOP = NIPPA
        NSIZE = NIPPA
    END IF

    !.....  Allocate tmp pol per species count (or dummy ol-per-pol)
    ALLOCATE( PPSCNT( NSIZE ), STAT=IOS )
    CALL CHECKMEM( IOS, 'PPSCNT', PROGNAME )
    PPSCNT = 0.      ! array

    !.....  Set the number of pollutants per species or emission type
    DO V = 1, NLOOP

        !.....  For speciation, get position of species
        IF( SFLAG ) THEN
            CBUF = TSVDESC( V )
            L1 = INDEX   ( CBUF, SPJOIN )
            L2 = LEN_TRIM( CBUF )
            CSPC = CBUF( L1+1:L2   )
            J = INDEX1( CSPC, NMSPC, EMNAM )           ! get position of species

        !.....  Otherwise, determine position in the list of pols/activities
        ELSE
            CPOL = EANAM( V )
            CBUF = CPOL
            L2 = LEN_TRIM( CBUF )
            J = INDEX1( CPOL, NIPPA, EANAM )           ! get position of pol/act

        END IF

        !.....  Make sure variable was found in its list
        IF( J .LE. 0 ) THEN

            EFLAG = .TRUE.
            MESG = 'ERROR: variable "'// CBUF( 1:L2 ) //    &
                   '" not found in master list.'
            CALL M3MSG2( MESG )

        !.....  Count the number per species or pollutant
        ELSE
            PPSCNT( J ) = PPSCNT( J ) + 1

        END IF

    END DO

    IF( EFLAG ) THEN
        MESG = 'Problem building indices'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF

    !.....  Get the maximum PPSCNT value
    MAXPPS = MAXVAL( PPSCNT )

    !.....  Make sure MXVARPGP is large enough for MAXPPS
    IF( MAXPPS .GT. MXVARPGP ) THEN

        WRITE( MESG,94010 )  'ERROR: the maximum variables '//    &
               'per group of', MXVARPGP, CRLF() // BLANK10 //    &
               'could not support the maximum pollutants ' //    &
               'per species of', MAXPPS, CRLF() // BLANK10 //    &
               'Not enough memory available to run for selected settings.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )

    END IF

    !.....  Allocate memory the group variables...
    ALLOCATE( VGRPCNT( MXGRP ),    &
                IDVGP( MXGRP ),    &
              GVNAMES( MXVARPGP, MXGRP ),    &
              GVLOUT ( MXVARPGP, MXGRP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'GVLOUT', PROGNAME )

    !.....  Initialize group variables
    VGRPCNT = 0        ! array
    IDVGP   = 0        ! array
    GVNAMES = ' '      ! array
    GVLOUT  = .FALSE.      ! array

    !.....  Allocate memory for all *EXIST arrays for a source category if
    !       that source category is present.  This is necessary to simplify
    !       using these arrays, and it will be a small memory penalty.
    !.....  NOTE - all arrays are dimensioned by MXVARGRP, but not all
    !       need that many (esp. the arrays that contain pollutant info when
    !       program is run with speciation)
    !.....  Inventory emissions: e.g., A_EXIST
    !.....  Multiplicative control matrices: e.g., AU_EXIST
    !.....  Reactivity control matrices: e.g., AR_EXIST
    !.....  Speciation matrices: e.g., AS_EXIST

    !.....  Allocate and initialize to zero...
    IF( AFLAG ) THEN       ! area

        ALLOCATE( A_EXIST ( MXVARPGP, MXGRP ),    &
                  AU_EXIST( MXVARPGP, MXGRP ),    &
                  AR_EXIST( MXVARPGP, MXGRP ),    &
                  AS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'AS_EXIST', PROGNAME )
        A_EXIST  = 0      ! array
        AU_EXIST = 0      ! array
        AR_EXIST = 0      ! array
        AS_EXIST = 0      ! array
    ENDIF

    IF( BFLAG ) THEN       ! biogenics
        ALLOCATE( BS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'BS_EXIST', PROGNAME )
        BS_EXIST = 0      ! array
    END IF

    IF( MFLAG ) THEN       ! mobile
        ALLOCATE( M_EXIST( MXVARPGP, MXGRP ),    &
                 MU_EXIST( MXVARPGP, MXGRP ),    &
                 MR_EXIST( MXVARPGP, MXGRP ),    &
                 MS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'MS_EXIST', PROGNAME )
        M_EXIST  = 0      ! array
        MU_EXIST = 0      ! array
        MR_EXIST = 0      ! array
        MS_EXIST = 0      ! array
    ENDIF

    IF( PFLAG ) THEN       ! point
        ALLOCATE(  P_EXIST( MXVARPGP, MXGRP ),    &
                  PU_EXIST( MXVARPGP, MXGRP ),    &
                  PR_EXIST( MXVARPGP, MXGRP ),    &
                  PS_EXIST( MXVARPGP, MXGRP ), STAT=IOS )
        CALL CHECKMEM( IOS, 'PS_EXIST', PROGNAME )
        P_EXIST  = 0      ! array
        PU_EXIST = 0      ! array
        PR_EXIST = 0      ! array
        PS_EXIST = 0      ! array
    ENDIF

    !.....  Allocate memory for variable to pollutant and variable to species
    !       indexes. Variable can be either pollutant or pol-to-species
    !.....  NOTE - SPINDEX will not be used if there is no speciation
    ALLOCATE( SIINDEX( MXVARPGP, MXGRP ),    &
              SPINDEX( MXVARPGP, MXGRP ), STAT=IOS )
    CALL CHECKMEM( IOS, 'SPINDEX', PROGNAME )
    SIINDEX = 0      ! array
    SPINDEX = 0      ! array

    !.....  Create pollutant names array in groups structure, and store counts
    !.....  For speciation, also build group indices for pol-to-species...

    !.....  Loop through pollutants (if no speciation) or pol-to-species names
    VCNT = 0         ! initialize variable count
    GCNT = 1         ! initialize group count
    PJ   = 0
    PSPC = ' '
    DO V = 1, NLOOP

        !.....  Extract pollutant name and species name (for speciation only)
        IF( SFLAG ) THEN
            CBUF = TSVDESC( V )
            L1 = INDEX   ( CBUF, SPJOIN )
            L2 = LEN_TRIM( CBUF )
            CPOL = CBUF(    1:L1-1 )
            CSPC = CBUF( L1+1:L2   )
        ELSE
            CPOL = EANAM( V )
            CBUF = CPOL
        END IF

        !.....  Determine location in the list of pollutants/activities
        J = INDEX1( CPOL, NIPPA, EANAM )
        K = J
        IF( SFLAG ) J = INDEX1( CSPC, NMSPC, EMNAM )

        !.....  Ensure that all pollutants for the same species or emission
        !       type will be in the same group
        IF( J .NE. PJ .AND. VCNT .GT. MXVARPGP ) THEN
            GCNT = GCNT + 1
            VCNT = 1

        !.....  Group count is not up to the maximum yet, increment variable
        !       count
        ELSE IF( VCNT .LT. MXVARPGP ) THEN

            VCNT = VCNT + 1

        !.....  Otherwise increment group count and initialize count
        ELSE
            GCNT = GCNT + 1
            VCNT = 1

        END IF

        !.....  Store variable names and number in group structure
        VGRPCNT(       GCNT ) = VCNT         ! store no. variables in group
        GVNAMES( VCNT, GCNT ) = CBUF         ! store rearranged variable names

        !.....  Store pol/act and species indices
        SIINDEX( VCNT, GCNT ) = K                 ! store pol/act index
        IF( SFLAG ) SPINDEX( VCNT, GCNT ) = J         ! store species index

        !.....  Determine the iterations on which output should occur
        IF( V     .NE. 1    .AND.    &
            SFLAG           .AND.    &
            CSPC  .NE. PSPC       ) THEN

            GVLOUT( PVCNT, PGCNT ) = .TRUE.

        ELSE IF( V    .NE.  1     .AND.    &
                      .NOT. SFLAG .AND.    &
                 CPOL .NE.  PPOL        ) THEN

            GVLOUT( PVCNT, PGCNT ) = .TRUE.

        END IF

        PSPC  = CSPC
        PPOL  = CPOL
        PVCNT = VCNT
        PGCNT = GCNT
        PJ    = J

    END DO

    !.....  Store actual number of groups
    NGRP = GCNT

    !.....  Ensure output on last variable in all groups
    GVLOUT( VGRPCNT( NGRP ), NGRP ) = .TRUE.

    !.....  Create pollutant-group array by looping through groups and checking
    !       if the pollutants in the current group are the same as those in the
    !       previous group...
    !.....  This array will be used in helping to decide whether pollutant
    !       data need to be read for a given group
    IDVGP( 1 ) = 1       ! initialize group number for the first group
    DO N = 2, NGRP

        NEXTGRP = .FALSE.    ! reset indicator for incrementing group
        TGRP = IDVGP( N-1 )      ! set temporary group number

        !.....  If any of the pollutants in this group are different from the
        !       previous group, turn on indicator for incrementing group
        DO V = 1, VGRPCNT( N )
            J = SIINDEX( V,N )
            K = SIINDEX( V,N-1 )
            IF( K .LE. 0 ) THEN
                NEXTGRP = .TRUE.
            ELSE IF( EANAM( J ) .NE. EANAM( K ) ) THEN
                NEXTGRP = .TRUE.
            END IF
        END DO

        !.....  If indicator is true, increment group
        IF( NEXTGRP ) THEN
            TGRP = TGRP + 1
        ENDIF

        !.....  Set group number for current group
        IDVGP( N ) = TGRP

    END DO

    !.....  Precompute the position of the species variables in the reactivity
    !       matrices
    JA   = ANRMATV - ARNMSPC + 1
    JM   = MNRMATV - MRNMSPC + 1
    JP   = PNRMATV - PRNMSPC + 1

    !.....  Create source-categeory-specific indicies for speciation, inventory
    !       emissions, and control matrices...

    !.....  If speciation, loop through the number of variables per
    !       group and store the count for each source category.
    !.....  The counts are kept for each source category because each
    !       source category has only its own speciation factors stored
    !       for only the pollutants in that source category.
    !.....  The counts are used because of the variable groups.  If the
    !       position in the variable list were stored, this would not be
    !       the correct index for a pollutant within a group.  The counts
    !       work because all of the pollutant/species are sorted in the same
    !       order for all source categories and files.
    !.....  If there is reactivity, then also store the position in the
    !       reactivity matrix

    IF ( SFLAG ) THEN

        !.....  Loop through all groups, then number of variables per group
        DO N = 1, NGRP

            ACNT = 0
            BCNT = 0
            MCNT = 0
            PCNT = 0

            DO V = 1, VGRPCNT( N )

                SVBUF = GVNAMES( V,N )

                IF ( AFLAG ) THEN    ! Area sources
                    K = INDEX1( SVBUF, ANSMATV, ASVDESC )
                    IF( K .GT. 0 ) THEN
                        ACNT = ACNT + 1
                        AS_EXIST( V,N ) = ACNT
                    END IF
                END IF

                IF ( ARFLAG ) THEN    ! Area reactivity
                    K = INDEX1( SVBUF, ARNMSPC, ARVDESC( JA ) )
                    AR_EXIST( V,N ) = K
                END IF

                IF ( BFLAG ) THEN    ! Biogenic sources
                    K = INDEX1( SVBUF, BNSMATV, BSVDESC )
                    IF( K .GT. 0 ) THEN
                        BCNT = BCNT + 1
                        BS_EXIST( V,N ) = BCNT
                    END IF
                END IF

                IF ( MFLAG ) THEN    ! Mobile sources
                    K = INDEX1( SVBUF, MNSMATV, MSVDESC )
                    IF( K .GT. 0 ) THEN
                        MCNT = MCNT + 1
                        MS_EXIST( V,N ) = MCNT
                    END IF
                END IF

                IF ( MRFLAG ) THEN    ! Mobile reactivity
                    K = INDEX1( SVBUF, MRNMSPC, MRVDESC( JM ) )
                    MR_EXIST( V,N ) = K
                END IF

                IF ( PFLAG ) THEN    ! Point sources
                    K = INDEX1( SVBUF, PNSMATV, PSVDESC )
                    IF( K .GT. 0 ) THEN
                        PCNT = PCNT + 1
                        PS_EXIST( V,N ) = PCNT
                    END IF
                END IF

                IF ( PRFLAG ) THEN    ! Point reactivity
                    K = INDEX1( SVBUF, PRNMSPC, PRVDESC( JP ) )
                    PR_EXIST( V,N ) = K
                END IF

            END DO      ! end loop on variables per group
        END DO          ! end loop on groups
    END IF              ! For speciation

    !.....  Now build indices for inventory emissions...
    !.....  For each source category, store per-group position for
    !       inventory pollutants...

    !.....  Loop through all groups, then number of variables per group
    DO N = 1, NGRP

        ACNT = 0
        MCNT = 0
        PCNT = 0
        KA   = 0      ! array
        KM   = 0      ! array
        KP   = 0      ! array

        DO V = 1, VGRPCNT( N )

            VBUF = EANAM( SIINDEX( V,N ) )

            IF ( AFLAG ) THEN    ! Area sources
                K = INDEX1( VBUF, ANIPOL, AEINAM )
                IF( K .GT. 0 ) THEN

                    !.....  If index has already been set for the pollutant,
                    !       then reuse it, otherwise update the counter and
                    !       store
                    IF( KA( K ) .NE. 0 ) THEN
                        A_EXIST( V,N ) = KA( K )
                    ELSE
                        ACNT = ACNT + 1
                        A_EXIST( V,N ) = ACNT
                        KA( K ) = ACNT
                    END IF
                END IF
            END IF

            IF ( MFLAG ) THEN    ! Mobile sources
                K = INDEX1( VBUF, MNIPPA, MEANAM )
                IF( K .GT. 0 ) THEN
                    IF( KM( K ) .NE. 0 ) THEN
                        M_EXIST( V,N ) = KM( K )
                    ELSE
                        MCNT = MCNT + 1
                        M_EXIST( V,N ) = MCNT
                        KM( K ) = MCNT
                    END IF
                END IF
            END IF

            IF ( PFLAG ) THEN    ! Point sources
                K = INDEX1( VBUF, PNIPOL, PEINAM )
                IF( K .GT. 0 ) THEN
                    IF( KP( K ) .NE. 0 ) THEN
                        P_EXIST( V,N ) = KP( K )
                    ELSE
                        PCNT = PCNT + 1
                        P_EXIST( V,N ) = PCNT
                        KP( K ) = PCNT
                    END IF
                END IF
            END IF
        END DO      ! end loop on pollutants per group
    END DO          ! end loop on groups

    !.....  Now build indices for multiplicative controls...
    !.....  For each source category, store per-group position for
    !       inventory pollutants...

    !.....  Allocate helper arrays for setup of mulitplicative control matrix
    IF( AUFLAG ) THEN
        ALLOCATE( KAU( ANUMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'KAU', PROGNAME )
    END IF

    IF( MUFLAG ) THEN
        ALLOCATE( KMU( MNUMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'KMU', PROGNAME )
    END IF

    IF( PUFLAG ) THEN
        ALLOCATE( KPU( PNUMATV ), STAT=IOS )
        CALL CHECKMEM( IOS, 'KPU', PROGNAME )
    END IF

    !.....  Loop through all groups, then number of pollutants per group
    DO N = 1, NGRP

        ACNT = 0
        MCNT = 0
        PCNT = 0
        IF( AUFLAG ) KAU  = 0       ! array
        IF( MUFLAG ) KMU  = 0       ! array
        IF( PUFLAG ) KPU  = 0       ! array

        DO V = 1, VGRPCNT( N )

            VBUF = EANAM( SIINDEX( V,N ) )

            IF ( AUFLAG ) THEN    ! Area sources
                K = INDEX1( VBUF, ANUMATV, AUVNAMS )
                IF( K .GT. 0 ) THEN
                    IF( KAU( K ) .NE. 0 ) THEN
                        AU_EXIST( V,N ) = KAU( K )
                    ELSE
                        ACNT = ACNT + 1
                        AU_EXIST( V,N ) = ACNT
                        KAU( K ) = ACNT
                    END IF
                END IF
            END IF

            IF ( MUFLAG ) THEN    ! Mobile sources
                K = INDEX1( VBUF, MNUMATV, MUVNAMS )
                IF( K .GT. 0 ) THEN
                    IF( KMU( K ) .NE. 0 ) THEN
                        MU_EXIST( V,N ) = KMU( K )
                    ELSE
                        MCNT = MCNT + 1
                        MU_EXIST( V,N ) = MCNT
                        KMU( K ) = MCNT
                    END IF
                END IF
            END IF

            IF ( PUFLAG ) THEN    ! Point sources
                K = INDEX1( VBUF, PNUMATV, PUVNAMS )
                IF( K .GT. 0 ) THEN
                    IF( KPU( K ) .NE. 0 ) THEN
                        PU_EXIST( V,N ) = KPU( K )
                    ELSE
                        PCNT = PCNT + 1
                        PU_EXIST( V,N ) = PCNT
                        KPU( K ) = PCNT
                    END IF
                END IF
            END IF

        END DO      ! end loop on pollutants per group
    END DO          ! end loop on groups

    !..... Deallocate local memory
    DEALLOCATE( PPSCNT )
    IF( AUFLAG ) DEALLOCATE( KAU )
    IF( MUFLAG ) DEALLOCATE( KMU )
    IF( PUFLAG ) DEALLOCATE( KPU )

    RETURN

    !******************  FORMAT  STATEMENTS   ******************************

    !.....   Internal buffering formats.....94xxx

94010 FORMAT( 10( A, :, I8, :, 1X ) )

END SUBROUTINE BLDMRGIDX
