
MODULE INPUT
    
    ! MESH FILE
    CHARACTER(LEN=*), PARAMETER :: MESHFILE = "mesh.rea"

    ! ELEMENT NUMBER
!    INTEGER, PARAMETER :: NELMAX = 768  ! MAXIMUM NUMBER OF ELEMENT IN EACH PROC
    INTEGER, PARAMETER :: SPLIT_MAX_NUM = 3

    ! CHOOSE LOAD PARTITION SCHEME
    LOGICAL :: LP_SERIAL = .TRUE.  ! USE SERIAL LP
    LOGICAL :: START_COARSEN = .FALSE. 


END MODULE INPUT
