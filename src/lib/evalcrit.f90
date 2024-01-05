
LOGICAL FUNCTION EVALCRIT( NV, NORS, MXAND, VALS, REFS, RANK,   &
                           CHRS, COMPARE, COMPCHR, TYPES, STATUS )

    !***********************************************************************
    !  function body starts at line
    !
    !  DESCRIPTION:
    !    This function evaluates values versus criteria for each of the value
    !    and a logical construction of rules with one or more levels of
    !    logic.  STATUS will provide TRUE for the first OR criterion that
    !    is TRUE. And ORs that are partially true (one but not all AND components
    !    of the OR are true) will return all components with a FALSE for STATUS.
    !
    !  PRECONDITIONS REQUIRED:
    !
    !  SUBROUTINES AND FUNCTIONS CALLED:
    !
    !  REVISION  HISTORY:
    !       Created ??/???? by ???
    !       Version 11/2023 by CJC:  USE M3UTILIO, conversion to ".f90",  and
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
    !***********************************************************************

    USE M3UTILIO

    IMPLICIT NONE

    !.......   ARGUMENTS and their descriptions:
    INTEGER     , INTENT (IN) :: NV          ! Number of values
    INTEGER     , INTENT (IN) :: NORS        ! Number of OR conditions
    INTEGER     , INTENT (IN) :: MXAND       ! Max no.  ANDs for single data val
    REAL        , INTENT (IN) :: VALS   ( NV )           ! Data values (real)
    REAL        , INTENT (IN) :: REFS   ( NV )           ! Reference values
    REAL        , INTENT (IN) :: RANK   ( NV )           ! Ranking order
    CHARACTER(*), INTENT (IN) :: CHRS   ( NV )           ! String values
    REAL        , INTENT (IN) :: COMPARE( NORS, MXAND, NV )     ! Formula values
    CHARACTER(*), INTENT (IN) :: COMPCHR( NORS, MXAND, NV )     ! Formula strings
    CHARACTER(6), INTENT (IN) :: TYPES  ( NORS, MXAND, NV )     ! Condition
    LOGICAL     , INTENT(OUT) :: STATUS ( NORS, MXAND, NV )     ! true: condition met

    !.......   OTHER LOCAL VARIABLES and their descriptions:
    INTEGER     L, L2, M, N          ! counters and indices

    REAL        REFMIN           ! tmp minimum for range check
    REAL        REFMAX           ! tmp maximum for range check
    REAL        TMPVAL           ! tmp value

    LOGICAL     ANDSTAT              ! true: all ands apply
    LOGICAL     EFLAG                ! true: error encountered

    CHARACTER(300)  MESG

    CHARACTER(16), PARAMETER :: PROGNAME = 'EVALCRIT'       !  program name

    !***********************************************************************
    !   begin body of function EVALCRIT

    EFLAG     = .FALSE.
    EVALCRIT  = .FALSE.
    STATUS    = .FALSE.        ! array

    !.....  Loop through OR conditions.  If any are true, then loop is done
    DO L = 1, NORS

        !.....  Loop through variables for each OR and check AND conditions
        !               if they are present
        ANDSTAT = .TRUE.            ! Initialize AND status for this OR
        DO M = 1, MXAND

            DO N = 1, NV

                SELECT CASE( TYPES( L,M,N ) )
                  CASE( '=', '==' )
                    IF ( VALS(N) .NE. COMPARE(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( 'IS' )
                    IF ( CHRS(N) .NE. COMPCHR(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '=>', '>=' )
                    IF ( VALS(N) .LT. COMPARE(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '=<', '<=' )
                    IF ( VALS(N) .GT. COMPARE(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '<' )
                    IF ( VALS(N) .GE. COMPARE(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '>' )
                    IF ( VALS(N) .LE. COMPARE(L,M,N) ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '+/-', '-/+' )
                    REFMIN = REFS( N ) - COMPARE( L,M,N )
                    REFMAX = REFS( N ) + COMPARE( L,M,N )
                    IF( VALS(N) .LT. REFMIN .OR.    &
                        VALS(N) .GT. REFMAX      ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( '%' )
                    TMPVAL = REFS( N ) * COMPARE( L,M,N ) * 0.01
                    REFMIN = REFS( N ) - TMPVAL
                    REFMAX = REFS( N ) + TMPVAL
                    IF( VALS(N) .LE. REFMIN .OR.    &
                        VALS(N) .GT. REFMAX      ) THEN
                        ANDSTAT = .FALSE.
                    END IF
                  CASE( 'TOP' )
                    IF ( RANK(N) .GT. COMPARE(L,M,N) .OR.   &
                         VALS(N) .EQ. AMISS3 ) THEN
                        ANDSTAT = .FALSE.

                    END IF
                  CASE( ' ' )      ! Skip if no (blank) case

                  CASE DEFAULT     ! Internal error

                    EFLAG = .TRUE.
                    MESG = 'INTERNAL ERROR: Do not know how to ' //     &
                           'interpret "' // TRIM( TYPES( L,M,N ) ) //   &
                           '" operation.'
                    CALL M3MSG2( MESG )

                END SELECT

                !.....  If any of the AND conditions are false, exit from loop
                !.....  Also reset status for all ANDs on this OR to FALSE
                IF ( .NOT. ANDSTAT ) THEN
                    STATUS( L, 1:MXAND, 1:NV ) = .FALSE.
                    EXIT
                ELSE IF ( TYPES( L,M,N ) .NE. ' ' ) THEN
                    STATUS( L, M, N ) = .TRUE.
                END IF

            END DO      ! End of variables loop

        END DO          ! End of ANDs loop

        !.....  Update OR status
        EVALCRIT = ANDSTAT

        !.....  If any OR status is true, then whole thing is true
        IF( EVALCRIT ) EXIT

    END DO        ! End of ORs loop

    !.....  Abort if error occurred
    IF( EFLAG ) THEN
        MESG = 'Problem interpreting selection criteria.'
        CALL M3EXIT( PROGNAME, 0, 0, MESG, 2 )
    END IF
   
    RETURN

END FUNCTION EVALCRIT
