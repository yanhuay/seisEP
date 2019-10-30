module seismo_parameters
! yanhuay@princeton.edu

use constants, only: IIN, IOUT, MAX_STRING_LEN,MAX_FILENAME_LEN,MAX_KERNEL_NUM, &
    MAX_LINES,MAX_MISFIT_TYPE, SIZE_REAL, SIZE_DOUBLE, CUSTOM_REAL,CUSTOM_COMPLEX,&
    LARGE_VAL, SMALL_VAL,PI

implicit none

!----------------------------------------------------------------------
! constants
! number of Gauss-Lobatto-Legendre (GLL) points (i.e., polynomial degree + 1)
INTEGER, PARAMETER :: NGLLX=5
! number of Gauss-Lobatto-Jacobi (GLJ) points in the axial elements (i.e.,
! polynomial degree + 1)
! the code does NOT work if NGLLZ /= NGLLX because it then cannot handle a
! non-structured mesh
! due to non matching polynomial degrees along common edges
INTEGER, PARAMETER :: NGLLZ=5
INTEGER, PARAMETER :: NGLLY=1

!! solver
CHARACTER (LEN=20) :: solver='specfem2D'
CHARACTER (LEN=MAX_STRING_LEN) :: LOCAL_PATH='OUTPUT_FILES'
CHARACTER (LEN=MAX_STRING_LEN) :: IBOOL_NAME='NSPEC_ibool.bin'

!! FORWARD MODELNG INFO
INTEGER, PARAMETER :: NSTEP=4800
REAL(KIND=CUSTOM_REAL), PARAMETER :: deltat=0.06
REAL(KIND=CUSTOM_REAL), PARAMETER :: t0=0.0
REAL(KIND=CUSTOM_REAL), PARAMETER :: f0=0.084
INTEGER, PARAMETER :: NREC=2
INTEGER, PARAMETER :: NSRC=1
CHARACTER (LEN=20) :: seismotype='displacement'
REAL(KIND=CUSTOM_REAL), PARAMETER :: noise_level=0.0

!! PRE-PROCESSING
! wavelet
INTEGER, PARAMETER :: Wscale=0
!window
LOGICAL :: TIME_WINDOW=.false.
INTEGER, PARAMETER :: window_type=3
REAL(KIND=CUSTOM_REAL), PARAMETER :: taper_percentage=0.2
CHARACTER (LEN=4) :: taper_type='hann'
REAL(KIND=CUSTOM_REAL), PARAMETER :: min_window_len=1.0/f0
REAL(KIND=CUSTOM_REAL), PARAMETER :: VEL_TOP=3900
REAL(KIND=CUSTOM_REAL), PARAMETER :: VEL_BOT=3100
! damping
LOGICAL :: DAMPING=.false.
REAL(KIND=CUSTOM_REAL), PARAMETER :: X_decay=1.0
REAL(KIND=CUSTOM_REAL), PARAMETER :: T_decay=1.0
! mute
LOGICAL :: MUTE_NEAR=.false.
REAL(KIND=CUSTOM_REAL), PARAMETER :: offset_near=0
LOGICAL :: MUTE_FAR=.false.
REAL(KIND=CUSTOM_REAL), PARAMETER :: offset_far=0
! event nomralize
LOGICAL :: EVENT_NORMALIZE=.false.
! trace nomralize
LOGICAL :: TRACE_NORMALIZE=.false.

!! measurement type weight 
INTEGER, PARAMETER :: mtype=MAX_MISFIT_TYPE
REAL(KIND=CUSTOM_REAL), DIMENSION(mtype) :: measurement_weight=1
LOGICAL :: uncertainty=.false.

!! DD par
REAL(KIND=CUSTOM_REAL), PARAMETER :: cc_threshold=0.9
REAL(KIND=CUSTOM_REAL), PARAMETER :: DD_min=SMALL_VAL
REAL(KIND=CUSTOM_REAL), PARAMETER :: DD_max=LARGE_VAL

!! OPTIMIZATION
CHARACTER (LEN=2) :: opt_scheme='QN'
INTEGER, PARAMETER :: CGSTEPMAX=10
CHARACTER (LEN=2) :: CG_scheme='PR'
INTEGER, PARAMETER :: BFGS_STEPMAX=4
REAL(KIND=CUSTOM_REAL), PARAMETER :: initial_step_length=0.04
INTEGER, PARAMETER :: max_step=5
REAL(KIND=CUSTOM_REAL), PARAMETER :: min_step_length=0.01
LOGICAL :: backtracking=.false.

!! CONVERGENCE?
INTEGER, PARAMETER :: iter_start=1
INTEGER, PARAMETER :: iter_end=1
REAL(KIND=CUSTOM_REAL), PARAMETER :: misfit_ratio_initial=0.001
REAL(KIND=CUSTOM_REAL), PARAMETER :: misfit_ratio_previous=0.001

!! POST-PROCESSING
LOGICAL :: smooth=.false.
LOGICAL :: MASK_SOURCE=.false.
LOGICAL :: MASK_STATION=.false.
REAL(KIND=CUSTOM_REAL), PARAMETER :: source_radius=0.0
REAL(KIND=CUSTOM_REAL), PARAMETER :: station_radius=0.0
LOGICAL :: precond=.false.
CHARACTER (LEN=MAX_STRING_LEN) :: precond_name=''
REAL(KIND=CUSTOM_REAL), PARAMETER :: wtr_precond=0.1

!! DISPLAY 
LOGICAL :: DISPLAY_DETAILS=.false.

!!!!!!!!!!!!!!!!! gloabl variables !!!!!!!!!!!!!!!!!!!!!!
INTEGER :: myrank,nproc,iproc

!! ADJOINT?
LOGICAL :: compute_adjoint

!! data
INTEGER :: ndata
INTEGER,PARAMETER :: MAX_DATA_NUM=4
INTEGER, DIMENSION(:), ALLOCATABLE :: which_proc_receiver
INTEGER :: nrec_proc  ! trace from a single proc
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_obs
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_syn
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj_AD
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:), ALLOCATABLE :: seism_adj_DD

REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: st_xval,st_yval,st_zval
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: dis_sr
REAL(KIND=CUSTOM_REAL) :: x_source, y_source, z_source
INTEGER(KIND=4) :: r4head(60)
INTEGER(KIND=2) :: header2(2)
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: time
REAL(KIND=CUSTOM_REAL) :: ratio_data_syn=0.01

!! model
INTEGER,PARAMETER :: MAX_PAR_NUM=3
!specfem intrinsic input model parameters in modeling
CHARACTER (LEN=50), DIMENSION(MAX_PAR_NUM) :: model_list=['rho','vp','vs']

!! measurement
CHARACTER (LEN=MAX_STRING_LEN) :: measurement_list
CHARACTER (LEN=MAX_STRING_LEN) :: misfit_type_list

!! window 
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE  :: win_start, win_end
REAL(KIND=CUSTOM_REAL) ::  tstart, tend
REAL(KIND=CUSTOM_REAL) ::  tstart_ref, tend_ref

! event nomralize
REAL(KIND=CUSTOM_REAL) :: event_norm
! trace nomralize
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: trace_norm

!! damping 
REAL(KIND=CUSTOM_REAL) :: lambda_x=1.0/X_decay
REAL(KIND=CUSTOM_REAL) :: lambda_t=1.0/T_decay

!! source-timefunction
LOGICAL :: conv_stf=.false.
CHARACTER (LEN=MAX_FILENAME_LEN) :: stf_file='source.txt'
REAL(KIND=CUSTOM_REAL) :: tshift_stf, integral_stf
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: stf
INTEGER :: stf_len

!! misfit
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: measurement_AD,measurement_AD_error
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: measurement_DD,measurement_DD_error
REAL(KIND=CUSTOM_REAL) :: mean_AD, var_AD, std_AD
REAL(KIND=CUSTOM_REAL) :: mean_DD, var_DD, std_DD
REAL(KIND=CUSTOM_REAL) :: misfit_AD, misfit_DD
INTEGER :: num_AD, num_DD
INTEGER, DIMENSION(:,:), ALLOCATABLE :: is_pair
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: misfit_proc
REAL(KIND=CUSTOM_REAL) :: misfit

!! kernels
INTEGER :: nspec
INTEGER, DIMENSION(:), ALLOCATABLE :: nspec_proc
INTEGER :: nker
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: g_new
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: p_new

!! models
INTEGER :: nmod
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: m_new
REAL(KIND=CUSTOM_REAL), DIMENSION(:), ALLOCATABLE :: m_try

!! linesearch
INTEGER :: is_done, is_cont, is_brak
REAL(KIND=CUSTOM_REAL) :: step_length, next_step_length,optimal_step_length

!! mask source
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: xstore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: ystore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: zstore
REAL(KIND=CUSTOM_REAL), DIMENSION(:,:,:,:), ALLOCATABLE :: mask
!----------------------------------------------------------------------

end module seismo_parameters
