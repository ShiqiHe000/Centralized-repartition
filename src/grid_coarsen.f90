

MODULE LOAD_COAREN

USE MPI
USE MESH, ONLY: NEL_LOCAL
USE ADAPT, ONLY: H_DEPTH, LPOLY_ORDER
USE SIZE

IMPLICIT NONE

!-----------------------------------------------------------------------
INTEGER, ALLOCATABLE, DIMENSION(:) :: H_COARSEN_FLAG    ! H_COARSEN_FLAG
INTEGER :: GMAX_H_DEPTH  ! MAXIMUM H_DEPTH AMONG LOCAL ELEMENTS

REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: GSPLIT_TRACKING   ! GLOBAL SPLIT TRACKING
!-----------------------------------------------------------------------

CONTAINS

SUBROUTINE P_COARSEN

    !-------------------------------------------------------------------
    ! DO P_COARSEN FIRST. 50% ELEMENT DO P_COARSEN. ONCE POLYNOMIAL ORDER
    ! <= 8 START H_COARSEN
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    INTEGER :: IEL
    
    REAL(KIND=8) :: RAND_NUM
    !-------------------------------------------------------------------
    
    ! DETERMINE WHICH ELEMENT TO COARSE---------------------------------
    DO IEL =1, NEL_LOCAL
        CALL RANDOM_NUMBER(RAND_NUM)
        
        IF (RAND_NUM >= 0.5D0 .AND. LPOLY_ORDER(IEL) > 4) THEN ! P_COARSEN
            LPOLY_ORDER(IEL) = LPOLY_ORDER(IEL)-2
        
        ENDIF
    
    ENDDO
    !-------------------------------------------------------------------

END SUBROUTINE P_COARSEN

SUBROUTINE SET_H_COARSEN_FLAG
    !-------------------------------------------------------------------
    ! IF POLYNOMIAL DEGREE <= 8, 50% CHANCE H_COARSEN
    !-------------------------------------------------------------------
    
    !-------------------------------------------------------------------
    INTEGER :: IEL
    
    REAL(KIND=8) :: RAND_NUM
    !-------------------------------------------------------------------

    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(H_COARSEN_FLAG(NEL_LOCAL))
    !-------------------------------------------------------------------
    
    ! INITIALIZATION----------------------------------------------------
    H_COARSEN_FLAG=0
    !-------------------------------------------------------------------
    
    ! GET THE MAXIMUM H_DEPTH AMONG PROCESSORS--------------------------
    CALL MARK_MAX_H_DEPTH
    
    IF(RANK == 0) THEN
        IF(GMAX_H_DEPTH == 1) THEN
            PRINT *, "ELEMENTS COARSEN COMPLETED."
        ENDIF
    ENDIF
    !-------------------------------------------------------------------
    
    
    ! H_COARSEN DECIDE--------------------------------------------------
    DO IEL=1, NEL_LOCAL
        
        CALL RANDOM_NUMBER(RAND_NUM)
        
        ! ONLY THE ELEMENT WITH THE HIGHEST H_DEPTH WILL BE COARSEN
        IF((RAND_NUM >= 0.5) .AND. (LPOLY_ORDER(IEL) <= 8) &
            .AND. (H_DEPTH(IEL) > 1) .AND. (H_DEPTH(IEL) == GMAX_H_DEPTH)) THEN
            H_COARSEN_FLAG(IEL) = 1 ! IF H_COARSEN, THEN FLAG == 1
        ENDIF
        
    ENDDO
    !-------------------------------------------------------------------
    
    
END SUBROUTINE SET_H_COARSEN_FLAG

SUBROUTINE MARK_MAX_H_DEPTH
    
    !-------------------------------------------------------------------
    INTEGER :: IEL
    INTEGER :: MAX_H_DEPTH 
    !-------------------------------------------------------------------
    
    ! INITIALIZE--------------------------------------------------------
    MAX_H_DEPTH = 0
    !-------------------------------------------------------------------
    
    ! GET MAXIMUM H_DEPTH, LOCAL----------------------------------------
    DO IEL=1, NEL_LOCAL
        IF(H_DEPTH(IEL) > MAX_H_DEPTH) MAX_H_DEPTH=H_DEPTH(IEL)
    ENDDO
    !-------------------------------------------------------------------
    
    ! SYCHRONIZE--------------------------------------------------------
    CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! FIND THE MAXIMUM H_DEPTH AMONG PROCESSORS-------------------------
    CALL MPI_ALLREDUCE(MAX_H_DEPTH, GMAX_H_DEPTH, 1, MPI_INT, MPI_MAX, &
                        MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
!    PRINT *, "MAX_H_DEPTH", MAX_H_DEPTH

END SUBROUTINE MARK_MAX_H_DEPTH



END MODULE LOAD_COAREN
