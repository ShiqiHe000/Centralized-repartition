
MODULE H_COARSENING

USE MPI
USE LOAD_COAREN
USE ADAPT, ONLY: SPLIT_TRACKING
USE FIELDS
USE SIZE

IMPLICIT NONE

INTEGER, ALLOCATABLE, DIMENSION(:) :: GH_COARSEN_FLAG
INTEGER :: NEL_TOTAL_NEW
INTEGER, ALLOCATABLE, DIMENSION(:) :: H_DEPTH_NEW

REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: GSPLIT_TRACKING_NEW
REAL(KIND=8), ALLOCATABLE, DIMENSION(:, :) :: X_GLOBAL_NEW, Y_GLOBAL_NEW

CONTAINS

SUBROUTINE H_COARSEN
    !-------------------------------------------------------------------
    ! IMPLEMENT H_COARSEN ON PROCESSOR 1
    !-------------------------------------------------------------------
    USE MESH, ONLY: NEL_LOCAL, NEL_TOTAL, NELMAX
    
    !-------------------------------------------------------------------
    INTEGER, ALLOCATABLE, DIMENSION(:) :: RECVCOUNT     ! RANK 0 RECV NUMBER OF ELEMENT ON EACH PROC
    INTEGER, ALLOCATABLE, DIMENSION(:) :: DISPLS1        ! THE DISPLACEMANT RELATIVE TO RECVBUF AT WHICH TO PLACE THE INCOMING DATA
    INTEGER :: I, J
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(RECVCOUNT(NUM_PROC))
    ALLOCATE(DISPLS1(NUM_PROC))
    IF (RANK == 0) THEN
        ALLOCATE(GH_COARSEN_FLAG(NEL_TOTAL))
        ALLOCATE(X_GLOBAL_NEW(4, NELMAX), Y_GLOBAL_NEW(4, NELMAX))
        ALLOCATE(GSPLIT_TRACKING_NEW(NELMAX))
        ALLOCATE(H_DEPTH_NEW(NELMAX))
    ENDIF
    !-------------------------------------------------------------------
    
    ! GATHRE NEL_LOCAL ON PROCESS 0-------------------------------------
    CALL MPI_GATHER(NEL_LOCAL, 1, MPI_INT, RECVCOUNT, 1, MPI_INT, 0, &
                    MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------

    ! THE POSITION TO PLACE DATA FROM DIFFERENT PROCESSES---------------
    IF (RANK == 0) THEN
         DISPLS1=0    ! INITIALIZE
         DO I=2, NUM_PROC
!            DISPLS1(I)=DISPLS1(I-1)+RECVCOUNT(I-1)*4
            DISPLS1(I)=DISPLS1(I-1)+RECVCOUNT(I-1)
        ENDDO
    
    ENDIF
    !-------------------------------------------------------------------
    
    ! GATHER LOCAL INFORMATION-------------------------------------------
    CALL MPI_GATHERV(H_COARSEN_FLAG, NEL_LOCAL, MPI_INT, GH_COARSEN_FLAG, &
            RECVCOUNT, DISPLS1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------

    ! H_COARSEN ON PROC 1-----------------------------------------------
    IF (RANK == 0) THEN
        ! INITIALIZE----------------------------------------------------
        NEL_TOTAL_NEW = 1
        !---------------------------------------------------------------
        
        DO I=1, NEL_TOTAL
!            WRITE(*,*) "GSPLIT_TRACKING=", GSPLIT_TRACKING(I)
            
            IF((GH_COARSEN_FLAG(I) == 1) .AND. (GH_COARSEN_FLAG(I+1) == 1) &
                .AND. (NINT(GSPLIT_TRACKING(I)) /= 4)) THEN
                CALL GATHER_ELEMENTS(I)
                
                
            ELSE
                X_GLOBAL_NEW(:, NEL_TOTAL_NEW) = X_GLOBAL(:, I)
                Y_GLOBAL_NEW(:, NEL_TOTAL_NEW) = Y_GLOBAL(:, I)
                
                GSPLIT_TRACKING_NEW(NEL_TOTAL_NEW) = GSPLIT_TRACKING(I)
                
                H_DEPTH_NEW(NEL_TOTAL_NEW) = GH_DEPTH(I)
                
                NEL_TOTAL_NEW = NEL_TOTAL_NEW+1
            ENDIF
        ENDDO
        
        ! UPDATE NEL_TOTAL----------------------------------------------
        NEL_TOTAL_NEW = NEL_TOTAL_NEW-1
        NEL_TOTAL = NEL_TOTAL_NEW
        !---------------------------------------------------------------
        
        ! UPDATE--------------------------------------------------------
        DEALLOCATE(X_GLOBAL, Y_GLOBAL)
        ALLOCATE(X_GLOBAL(4, NEL_TOTAL), Y_GLOBAL(4, NEL_TOTAL))
        X_GLOBAL(:, :) = X_GLOBAL_NEW(:, 1:NEL_TOTAL)
        Y_GLOBAL(:, :) = Y_GLOBAL_NEW(:, 1:NEL_TOTAL)
        
        
        DEALLOCATE(GSPLIT_TRACKING)
        ALLOCATE(GSPLIT_TRACKING(NEL_TOTAL))
        GSPLIT_TRACKING(:) = GSPLIT_TRACKING_NEW(1:NEL_TOTAL)
        
        
        DEALLOCATE(GH_DEPTH)
        ALLOCATE(GH_DEPTH(NEL_TOTAL))
        GH_DEPTH(:) = H_DEPTH_NEW(1:NEL_TOTAL)
        !---------------------------------------------------------------
        
        ! DEALLOCATE----------------------------------------------------
        DEALLOCATE(X_GLOBAL_NEW, Y_GLOBAL_NEW)
        DEALLOCATE(GSPLIT_TRACKING_NEW)
        DEALLOCATE(H_DEPTH_NEW)
        DEALLOCATE(GH_COARSEN_FLAG)
        !---------------------------------------------------------------
        
    ENDIF
    !-------------------------------------------------------------------
    
    ! BROADCAST NEL_TOTAL-----------------------------------------------
    CALL MPI_BCAST(NEL_TOTAL, 1, MPI_INT, 0, MPI_COMM_WORLD, IERROR)
    !-------------------------------------------------------------------
    
    ! DEALLOCATE--------------------------------------------------------
    DEALLOCATE(H_COARSEN_FLAG)
    !-------------------------------------------------------------------

    
END SUBROUTINE H_COARSEN

SUBROUTINE GATHER_ELEMENTS(I)

    !-------------------------------------------------------------------
    INTEGER :: I
    !-------------------------------------------------------------------

    ! 3 COMBIATIONS-----------------------------------------------------
    IF(NINT(GSPLIT_TRACKING(I)) == 1) THEN
        CALL CASE_ONE(I)
    
    ELSEIF(NINT(GSPLIT_TRACKING(I)) == 2) THEN
        CALL CASE_TWO(I)
        
    ELSEIF(NINT(GSPLIT_TRACKING(I)) == 3) THEN
        CALL CASE_THREE(I)
        
    ELSE
        PRINT *, "BUG"
        PRINT *, "GSPLIT_TRACKING=", GSPLIT_TRACKING(I)
        PRINT *, "I=", I
    
    ENDIF

END SUBROUTINE GATHER_ELEMENTS

SUBROUTINE CASE_ONE(I)
    
    !-------------------------------------------------------------------
    INTEGER :: I
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: XX, YY
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(XX(2), YY(2))
    !-------------------------------------------------------------------
    
    ! INITIALIZE--------------------------------------------------------
    XX=0; YY=0
    !-------------------------------------------------------------------
    
    ! CALCULATE NEW COORDINATE------------------------------------------
    XX(1) = X_GLOBAL(1, I)
    XX(2) = X_GLOBAL(2, I+1)
    
    YY(1) = Y_GLOBAL(1, I)
    YY(2) = Y_GLOBAL(3, I+2)
    
    X_GLOBAL_NEW(1, NEL_TOTAL_NEW) = XX(1)
    X_GLOBAL_NEW(2, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(3, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(4, NEL_TOTAL_NEW) = XX(1)
    
    Y_GLOBAL_NEW(1, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(2, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(3, NEL_TOTAL_NEW) = YY(2)
    Y_GLOBAL_NEW(4, NEL_TOTAL_NEW) = YY(2)
    !-------------------------------------------------------------------
    
    ! UPDATE GSPLIT_TRACKING--------------------------------------------
    CALL UPDATE_GSPLIT_TRACKING(I)
    !-------------------------------------------------------------------
    
    ! UPDATE  ----------------------------------------------------------
    H_DEPTH_NEW(NEL_TOTAL_NEW) = GH_DEPTH(I)-1
    I=I+3
    NEL_TOTAL_NEW=NEL_TOTAL_NEW+1
    !-------------------------------------------------------------------
    
END SUBROUTINE CASE_ONE

SUBROUTINE CASE_TWO(I)
    
    !-------------------------------------------------------------------
    INTEGER :: I
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: XX, YY
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(XX(2), YY(2))
    !-------------------------------------------------------------------
    
    ! INITIALIZE--------------------------------------------------------
    XX=0; YY=0
    NEL_TOTAL_NEW = NEL_TOTAL_NEW-1
    !-------------------------------------------------------------------
    
    ! CALCULATE NEW COORDINATE------------------------------------------
    XX(1) = X_GLOBAL(1, I+1)
    XX(2) = X_GLOBAL(2, I)
    
    YY(1) = Y_GLOBAL(1, I)
    YY(2) = Y_GLOBAL(3, I+1)
    
    X_GLOBAL_NEW(1, NEL_TOTAL_NEW) = XX(1)
    X_GLOBAL_NEW(2, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(3, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(4, NEL_TOTAL_NEW) = XX(1)
    
    Y_GLOBAL_NEW(1, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(2, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(3, NEL_TOTAL_NEW) = YY(2)
    Y_GLOBAL_NEW(4, NEL_TOTAL_NEW) = YY(2)
    !-------------------------------------------------------------------
    
    ! UPDATE GSPLIT_TRACKING--------------------------------------------
    CALL UPDATE_GSPLIT_TRACKING(I)
    !-------------------------------------------------------------------
    
    ! UPDATE  ----------------------------------------------------------
    H_DEPTH_NEW(NEL_TOTAL_NEW) = GH_DEPTH(I)-1
    I=I+2
    NEL_TOTAL_NEW=NEL_TOTAL_NEW+1
    !-------------------------------------------------------------------

END SUBROUTINE CASE_TWO

SUBROUTINE CASE_THREE(I)
    
    !-------------------------------------------------------------------
    INTEGER :: I
    
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: XX, YY
    !-------------------------------------------------------------------
    
    ! ALLOCATE----------------------------------------------------------
    ALLOCATE(XX(2), YY(2))
    !-------------------------------------------------------------------
    
    ! INITIALIZE--------------------------------------------------------
    XX=0; YY=0
    NEL_TOTAL_NEW = NEL_TOTAL_NEW-2
    !-------------------------------------------------------------------
    
    ! CALCULATE NEW COORDINATE------------------------------------------
    XX(1) = X_GLOBAL(1, I)
    XX(2) = X_GLOBAL(2, I+1)
    
    YY(1) = Y_GLOBAL(1, I-1)
    YY(2) = Y_GLOBAL(3, I)
    
    X_GLOBAL_NEW(1, NEL_TOTAL_NEW) = XX(1)
    X_GLOBAL_NEW(2, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(3, NEL_TOTAL_NEW) = XX(2)
    X_GLOBAL_NEW(4, NEL_TOTAL_NEW) = XX(1)
    
    Y_GLOBAL_NEW(1, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(2, NEL_TOTAL_NEW) = YY(1)
    Y_GLOBAL_NEW(3, NEL_TOTAL_NEW) = YY(2)
    Y_GLOBAL_NEW(4, NEL_TOTAL_NEW) = YY(2)
    !-------------------------------------------------------------------
    
    ! UPDATE GSPLIT_TRACKING--------------------------------------------
    CALL UPDATE_GSPLIT_TRACKING(I)
    !-------------------------------------------------------------------
    
    ! UPDATE  ----------------------------------------------------------
    H_DEPTH_NEW(NEL_TOTAL_NEW) = GH_DEPTH(I)-1
    I=I+1
    NEL_TOTAL_NEW=NEL_TOTAL_NEW+1
    !-------------------------------------------------------------------

END SUBROUTINE CASE_THREE

SUBROUTINE UPDATE_GSPLIT_TRACKING(I)

    !-------------------------------------------------------------------
    INTEGER :: I
    !-------------------------------------------------------------------

    GSPLIT_TRACKING_NEW(NEL_TOTAL_NEW) = (GSPLIT_TRACKING(I) - &
                                        REAL(NINT(GSPLIT_TRACKING(I)), KIND=8))*10.0D0

END SUBROUTINE UPDATE_GSPLIT_TRACKING



END MODULE H_COARSENING
