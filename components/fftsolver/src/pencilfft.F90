!> This Pencil FFT performs 3D forward and backwards FFTs using pencil decomposition. It uses FFTPACK for the actual FFT kernel
!! and this module contains all the data decomposition around this. There is no FFT required in Z, so this performs FFTs in
!! Y and X (in that order forward and reversed backwards.) The data decomposition is the complex aspect, there is the concept of
!! forward and backwards transformations. Forward transformations will go from pencil Z to Y to X and the backwards transformations
!! undo these, so go from X to Y to Z.
!! Note that we use quite a lot of buffer space here, this could be cut down if Y=X dimensions so some optimisation on memory
!! could be done there in that case
module pencil_fft_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX, global_grid_type
  use state_mod, only : model_state_type
  use, intrinsic :: iso_c_binding
  use mpi, only : MPI_DOUBLE_COMPLEX, MPI_INT, MPI_COMM_SELF,MPI_Wtime, mpi_comm_rank
  use fftnorth_mod, only: fftn_r2c, fftn_c2r, fftn_init, fftn_finalise
  use optionsdatabase_mod, only: options_get_integer
  use logging_mod, only: LOG_INFO, log_master_log, log_master_newline
  use conversions_mod, only : conv_to_string
#ifdef GPU
#if defined(CUDA)
  use cudafor
  use cufft
#elif defined(ROCM)
#endif
#endif
#ifdef CUDECOMP
  use cudecomp
#endif
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Describes a specific pencil transposition, from one pencil decomposition to another
  type pencil_transposition
     integer :: my_pencil_size(3), process_decomposition_layout(3), my_process_location(3), dim
     integer, dimension(:), allocatable :: send_sizes, send_offsets, recv_sizes, recv_offsets
     integer, dimension(:,:), allocatable :: recv_dims, send_dims
  end type pencil_transposition

  integer, parameter :: FORWARD=1, BACKWARD=2   !< Transposition directions
  integer :: dim_y_comm, dim_x_comm     !< Communicators for each dimension
  ! Transpositions from one pencil to another
  type(pencil_transposition) :: y_from_z_transposition, x_from_y_transposition, y_from_x_transposition, z_from_y_transposition, &
       y_from_z_2_transposition, x_from_y_2_transposition, y_from_x_2_transposition, z_from_y_2_transposition

  ! Temporary buffers used in transposition
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer :: real_buffer1, real_buffer2, real_buffer3, &
       fft_in_y_buffer , fft_in_x_buffer       
  complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer :: buffer1, buffer2
  
  public initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  
  !counters for number of times the fft routines are called and the time spent in them
  integer :: nforward, nback
  real(kind=DEFAULT_PRECISION) :: tforward, tback

  ! plans for cufft/rocfft
  integer :: plans(4)=0
  logical :: plan_defined(4)=.false.

  integer :: rank

#ifdef CUDECOMP
  type(cudecompGridDescConfig) :: config_cudecomp_T1, config_cudecomp_T2, config_cudecomp_T3
  type(cudecompHandle) :: handle_cudecomp_T1, handle_cudecomp_T2, handle_cudecomp_T3
  type(cudecompGridDesc) :: grid_desc_cudecomp_T1, grid_desc_cudecomp_T2, grid_desc_cudecomp_T3
  type(cudecompPencilInfo) :: pinfo_T1_x, pinfo_T1_y, pinfo_T1_z
  type(cudecompPencilInfo) :: pinfo_T2_x, pinfo_T2_y, pinfo_T2_z
  type(cudecompPencilInfo) :: pinfo_T3_x, pinfo_T3_y, pinfo_T3_z
  real(real64), allocatable :: data_cudecomp(:)
  real(real64), target, allocatable, device :: data_cudecomp_d(:)
  integer(8) :: data_num_elements, data_num_elements_T1, data_num_elements_T2, data_num_elements_T3
  integer(8) :: transpose_work_num_elements
  real(kind=DEFAULT_PRECISION), pointer, device, contiguous :: transpose_work_d(:)
  real(kind=DEFAULT_PRECISION), pointer, device, contiguous, dimension(:,:,:) :: data_cudecomp_T1_x_d, data_cudecomp_T1_y_d, data_cudecomp_T1_z_d
  real(kind=DEFAULT_PRECISION), pointer, device, contiguous, dimension(:,:,:) :: data_cudecomp_T2_x_d, data_cudecomp_T2_y_d, data_cudecomp_T2_z_d
  real(kind=DEFAULT_PRECISION), pointer, device, contiguous, dimension(:,:,:) :: data_cudecomp_T3_x_d, data_cudecomp_T3_y_d, data_cudecomp_T3_z_d
#endif

contains

  !> Initialises the pencil FFT functionality, this will create the transposition structures needed
  !! @param current_state The current model state
  !! @param my_y_start My global start in fourier space for Y
  !! @param my_x_start My global start in fourier space for X
  !! @returns Size of local dimensions in fourier space for this process
  function initialise_pencil_fft(current_state, my_y_start, my_x_start)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(out) :: my_y_start, my_x_start
    integer :: initialise_pencil_fft(3) 
    logical :: print_info

    integer :: ierr, y_distinct_sizes(current_state%parallel%dim_sizes(Y_INDEX)), &
         x_distinct_sizes(current_state%parallel%dim_sizes(X_INDEX))

    my_y_start=deduce_my_global_start(current_state, Y_INDEX)
    my_x_start=deduce_my_global_start(current_state, X_INDEX)

    if (current_state%parallel%dim_sizes(Y_INDEX) .gt. 1 .and. current_state%parallel%dim_sizes(X_INDEX) .gt. 1) then
      call mpi_cart_sub(current_state%parallel%neighbour_comm, (/1,0/), dim_y_comm, ierr)
      call mpi_cart_sub(current_state%parallel%neighbour_comm, (/0,1/), dim_x_comm, ierr)

      call mpi_allgather(current_state%local_grid%size(Y_INDEX), 1, MPI_INT, y_distinct_sizes, 1, MPI_INT, dim_y_comm, ierr)
      call mpi_allgather(current_state%local_grid%size(X_INDEX), 1, MPI_INT, x_distinct_sizes, 1, MPI_INT, dim_x_comm, ierr)
    else if (current_state%parallel%dim_sizes(Y_INDEX) .gt. 1) then
      dim_y_comm=current_state%parallel%monc_communicator
      dim_x_comm=MPI_COMM_SELF
      call mpi_allgather(current_state%local_grid%size(Y_INDEX), 1, MPI_INT, y_distinct_sizes, 1, MPI_INT, dim_y_comm, ierr)
      x_distinct_sizes=current_state%local_grid%size(X_INDEX)
    else if (current_state%parallel%dim_sizes(X_INDEX) .gt. 1) then      
      dim_y_comm=MPI_COMM_SELF
      dim_x_comm=current_state%parallel%monc_communicator
      y_distinct_sizes=current_state%local_grid%size(Y_INDEX)
      call mpi_allgather(current_state%local_grid%size(X_INDEX), 1, MPI_INT, x_distinct_sizes, 1, MPI_INT, dim_x_comm, ierr)
    else
      dim_y_comm=MPI_COMM_SELF
      dim_x_comm=MPI_COMM_SELF
      y_distinct_sizes=current_state%local_grid%size(Y_INDEX)
      x_distinct_sizes=current_state%local_grid%size(X_INDEX)
    end if

    call initialise_transpositions(current_state, y_distinct_sizes, x_distinct_sizes)
#ifdef CUDECOMP
    call initialise_cudecomp(current_state)
#endif
    call mpi_comm_rank(current_state%parallel%monc_communicator, rank, ierr)
    call initialise_buffers()

    initialise_pencil_fft=z_from_y_transposition%my_pencil_size

    call log_master_newline()

    

    print_info=.true.
    if (print_info) then
      write(*,*)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T1 cudecompinfo x-pen (monc z), (x,y,z):", pinfo_T1_x%shape(1), pinfo_T1_x%shape(2), pinfo_T1_x%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T1 cudecompinfo y-pen(y,z,x):", pinfo_T1_y%shape(1), pinfo_T1_y%shape(2), pinfo_T1_y%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T1 cudecompinfo z-pen (monc x), (z,x,y):", pinfo_T1_z%shape(1), pinfo_T1_z%shape(2), pinfo_T1_z%shape(3)

      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T2 cudecompinfo x-pen (monc z), (x,y,z):", pinfo_T2_x%shape(1), pinfo_T2_x%shape(2), pinfo_T2_x%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T2 cudecompinfo y-pen(y,z,x):", pinfo_T2_y%shape(1), pinfo_T2_y%shape(2), pinfo_T2_y%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T2 cudecompinfo z-pen (monc x), (z,x,y):", pinfo_T2_z%shape(1), pinfo_T2_z%shape(2), pinfo_T2_z%shape(3)

      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T3 cudecompinfo x-pen (monc z), (x,y,z):", pinfo_T3_x%shape(1), pinfo_T3_x%shape(2), pinfo_T3_x%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T3 cudecompinfo y-pen(y,z,x):", pinfo_T3_y%shape(1), pinfo_T3_y%shape(2), pinfo_T3_y%shape(3)
      write(*,'(a5,i2,a42,3i4)') "rank:",rank,",T3 cudecompinfo z-pen (monc x), (z,x,y):", pinfo_T3_z%shape(1), pinfo_T3_z%shape(2), pinfo_T3_z%shape(3)

      write(*,'(a5,i2,a40,3i4)') "rank:",rank,", monc decompinfo z-pen(z,y,x):", z_from_y_transposition%my_pencil_size(1),z_from_y_transposition%my_pencil_size(2),z_from_y_transposition%my_pencil_size(3)
      write(*,'(a5,i2,a40,3i4)') "rank:", rank," monc decompinfo y-pen(y,z,x):", size(real_buffer1,1), size(real_buffer1,2), size(real_buffer1,3)
      write(*,'(a5,i2,a40,3i4)') "rank:",rank,", monc decompinfo x-pen(x,y,z):", size(real_buffer2,1), size(real_buffer2,2), size(real_buffer2,3)
    endif
    
    nforward = 0
    nback = 0

    tforward = 0.d0
    tback = 0d0

    
  end function initialise_pencil_fft

  subroutine print_shape_r(array, arrname)
    real(C_DOUBLE_COMPLEX), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: arrname
    integer :: n1, n2, n3
    n1=size(array,1)
    n2=size(array,2)
    n3=size(array,3)
    write(*,'(a5,i2,a1,a20,a11,3i4)') "RANK:",rank, ",",arrname, ", n1,n2,n3:", n1, n2, n3
  end subroutine print_shape_r

  subroutine print_shape_c(array, arrname)
    complex(C_DOUBLE_COMPLEX), intent(in) :: array(:,:,:)
    character(len=*), intent(in) :: arrname
    integer :: n1, n2, n3
    n1=size(array,1)
    n2=size(array,2)
    n3=size(array,3)
    write(*,'(a5,i2,a1,a20,a11,3i4)') "RANK:",rank, ",",arrname, ", n1,n2,n3:", n1, n2, n3
  end subroutine print_shape_c


#ifdef CUDECOMP
  subroutine initialise_cudecomp(current_state)
    type(model_state_type), intent(inout) :: current_state
    integer :: istat, p_rows, p_cols, ierr
    logical :: print_info
    ierr = cudaSetDevice(0)    
    istat = cudecompInit(handle_cudecomp_T1, current_state%parallel%monc_communicator)
    !istat = cudecompInit(handle_cudecomp_T2, current_state%parallel%monc_communicator)
    !istat = cudecompInit(handle_cudecomp_T3, current_state%parallel%monc_communicator)

    istat = cudecompGridDescConfigSetDefaults(config_cudecomp_T1)
    istat = cudecompGridDescConfigSetDefaults(config_cudecomp_T2)
    istat = cudecompGridDescConfigSetDefaults(config_cudecomp_T3)

    config_cudecomp_T1%pdims = [current_state%parallel%dim_sizes(Y_INDEX), current_state%parallel%dim_sizes(X_INDEX)]
    config_cudecomp_T2%pdims = [current_state%parallel%dim_sizes(Y_INDEX), current_state%parallel%dim_sizes(X_INDEX)]
    config_cudecomp_T3%pdims = [current_state%parallel%dim_sizes(Y_INDEX), current_state%parallel%dim_sizes(X_INDEX)]

    config_cudecomp_T1%gdims = [current_state%global_grid%size(Z_INDEX), (current_state%global_grid%size(Y_INDEX)),(current_state%global_grid%size(X_INDEX))]
    config_cudecomp_T2%gdims = [current_state%global_grid%size(Z_INDEX), (current_state%global_grid%size(Y_INDEX)/2+1)*2,current_state%global_grid%size(X_INDEX)] !< change these
    config_cudecomp_T3%gdims = [current_state%global_grid%size(Z_INDEX), (current_state%global_grid%size(Y_INDEX)/2+1)*2,(current_state%global_grid%size(X_INDEX)/2+1)*2] !< change these

    config_cudecomp_T1%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
    config_cudecomp_T2%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P
    config_cudecomp_T3%transpose_comm_backend = CUDECOMP_TRANSPOSE_COMM_MPI_P2P

    config_cudecomp_T1%transpose_axis_contiguous = [.true., .true., .true.]
    config_cudecomp_T2%transpose_axis_contiguous = [.true., .true., .true.]
    config_cudecomp_T3%transpose_axis_contiguous = [.true., .true., .true.]

    istat = cudecompGridDescCreate(handle_cudecomp_T1, grid_desc_cudecomp_T1, config_cudecomp_T1)
    istat = cudecompGridDescCreate(handle_cudecomp_T1, grid_desc_cudecomp_T2, config_cudecomp_T2)
    istat = cudecompGridDescCreate(handle_cudecomp_T1, grid_desc_cudecomp_T3, config_cudecomp_T3)

    ! Get pencil information for T1 (with halo elements)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T1, pinfo_T1_x, 1)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T1, pinfo_T1_y, 2)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T1, pinfo_T1_z, 3)

    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T2, pinfo_T2_x, 1)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T2, pinfo_T2_y, 2)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T2, pinfo_T2_z, 3)

    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T3, pinfo_T3_x, 1)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T3, pinfo_T3_y, 2)
    istat = cudecompGetPencilInfo(handle_cudecomp_T1, grid_desc_cudecomp_T3, pinfo_T3_z, 3)

    data_num_elements_T1 = max(pinfo_T1_x%size, pinfo_T1_y%size, pinfo_T1_z%size)
    data_num_elements_T2 = max(pinfo_T2_x%size, pinfo_T2_y%size, pinfo_T2_z%size)
    data_num_elements_T3 = max(pinfo_T3_x%size, pinfo_T3_y%size, pinfo_T3_z%size)

    data_num_elements = max(data_num_elements_T1, data_num_elements_T2, data_num_elements_T3)

    ! Allocate device buffer
    allocate(data_cudecomp_d(data_num_elements))

    ! Allocate host buffer
    allocate(data_cudecomp(data_num_elements))

    !> Workspace allocation - this might not work properly, so can try to use CUDECOMP_WORKSPACE_AUTO in transpose calls
    !istat = cudecompGetTransposeWorkspaceSize(handle_cudecomp_T1, grid_desc_cudecomp_T1, transpose_work_num_elements)
    !istat = cudecompGetTransposeWorkspaceSize(handle_cudecomp_T2, grid_desc_cudecomp_T2, transpose_work_num_elements)
    !istat = cudecompGetTransposeWorkspaceSize(handle_cudecomp_T3, grid_desc_cudecomp_T3, transpose_work_num_elements)
    !call CHECK_CUDECOMP_EXIT(istat)

    !istat = cudecompMalloc(handle_cudecomp_T3, grid_desc_cudecomp_T3, transpose_work_d, transpose_work_num_elements)
    !call CHECK_CUDECOMP_EXIT(istat)

    data_cudecomp_T1_x_d(1:pinfo_T1_x%shape(1), 1:pinfo_T1_x%shape(2), 1:pinfo_T1_x%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T1_y_d(1:pinfo_T1_y%shape(1), 1:pinfo_T1_y%shape(2), 1:pinfo_T1_y%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T1_z_d(1:pinfo_T1_z%shape(1), 1:pinfo_T1_z%shape(2), 1:pinfo_T1_z%shape(3)) => data_cudecomp_d(:)

    data_cudecomp_T2_x_d(1:pinfo_T2_x%shape(1), 1:pinfo_T2_x%shape(2), 1:pinfo_T2_x%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T2_y_d(1:pinfo_T2_y%shape(1), 1:pinfo_T2_y%shape(2), 1:pinfo_T2_y%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T2_z_d(1:pinfo_T2_z%shape(1), 1:pinfo_T2_z%shape(2), 1:pinfo_T2_z%shape(3)) => data_cudecomp_d(:)

    data_cudecomp_T3_x_d(1:pinfo_T3_x%shape(1), 1:pinfo_T3_x%shape(2), 1:pinfo_T3_x%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T3_y_d(1:pinfo_T3_y%shape(1), 1:pinfo_T3_y%shape(2), 1:pinfo_T3_y%shape(3)) => data_cudecomp_d(:)
    data_cudecomp_T3_z_d(1:pinfo_T3_z%shape(1), 1:pinfo_T3_z%shape(2), 1:pinfo_T3_z%shape(3)) => data_cudecomp_d(:)
    
  end subroutine initialise_cudecomp
#endif

  !> Cleans up allocated buffer memory
  subroutine finalise_pencil_fft(monc_communicator)
    integer, intent(in) :: monc_communicator
    integer :: ierr, i
    real(kind=DEFAULT_PRECISION) :: t_per_forward, t_per_backward

    t_per_forward = tforward/nforward
    t_per_backward= tback/nback
    !call MPI_Allreduce(tforward,tgf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr))
    !call MPI_Allreduce(tbackward,tgb,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr))
    call MPI_Comm_rank(monc_communicator,i,ierr)
#if defined(GPU)
    call log_master_log(LOG_INFO, "Using GPU for FFT")
#else
    call log_master_log(LOG_INFO, "Using CPU for FFT")
#endif
    call log_master_log(LOG_INFO, "Total time in forward FFT ="//trim(conv_to_string(tforward))//" s")
    call log_master_log(LOG_INFO, "Total time in backward FFT ="//trim(conv_to_string(tback))//" s")
    call log_master_log(LOG_INFO, "Time per forward FFT ="//trim(conv_to_string(t_per_forward))//" s")
    call log_master_log(LOG_INFO, "Time per backward FFT ="//trim(conv_to_string(t_per_backward))//" s")

    if (dim_y_comm .ne. MPI_COMM_SELF .and. dim_y_comm .ne. monc_communicator) call mpi_comm_free(dim_y_comm, ierr)
    if (dim_x_comm .ne. MPI_COMM_SELF .and. dim_x_comm .ne. monc_communicator) call mpi_comm_free(dim_x_comm, ierr)
    deallocate(buffer1, buffer2, real_buffer1, real_buffer2, real_buffer3, fft_in_y_buffer , fft_in_x_buffer)
  end subroutine finalise_pencil_fft  

  !> Performs a forward 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
  !! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
  !! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and 
  !! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
  !! @param current_state The current model state
  !! @param source_data The source real data to in the time domain
  !! @param target_data Frequency domain real representation of the time domain source which is allocated here
  subroutine perform_forward_3dfft(current_state, source_data, target_data)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: target_data
   
    call transpose_and_forward_fft_in_y(current_state, source_data, buffer1, real_buffer1)
    call transpose_and_forward_fft_in_x(current_state, real_buffer1, buffer2, real_buffer2)
#ifndef CUDECOMP
    call transpose_to_pencil(y_from_x_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
         real_buffer2, real_buffer3)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_x with monc (Backward)"
#else
    call transpose_x_to_y_cudecomp_T3(real_buffer2, real_buffer3)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_x with cudecomp (Backward), T3"
#endif

#ifndef CUDECOMP
    call transpose_to_pencil(z_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
       real_buffer3, target_data)  
    if (rank==0) write(*,*) "[DEBUG] Transpose z_from_y with monc (Backward)"
#else 
    call transpose_y_to_z_cudecomp_T3(real_buffer3, target_data)
    if (rank==0) write(*,*) "[DEBUG] Transpose z_from_y with cudecomp (Backward), T3"
#endif
  end subroutine perform_forward_3dfft

  !> Performs a backwards 3D FFT and currently results in target data which is the X, Z, Y oriented pencil
  !! Note that the source_data here takes no account for the halo, it is up to caller to exclude this.
  !! This does no FFT in Z, but transposes to Y, does FFT in Y, then transposes to X and 
  !! performs an FFT in that dimension. Pencil decomposition is used which has already been set up.
  !! @param current_state The current model state
  !! @param source_data The source real data to in the frequency domain
  !! @param target_data Time domain complex representation of the frequency domain source
  subroutine perform_backwards_3dfft(current_state, source_data, target_data)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: target_data
#ifndef CUDECOMP
    call transpose_to_pencil(y_from_z_2_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
       source_data, real_buffer3)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_z_2 with monc (Forward)"
#else
    call transpose_z_to_y_cudecomp_T3(source_data, real_buffer3)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_z_2 with cudecomp (Forward), T3"
#endif

#ifndef CUDECOMP
    call transpose_to_pencil(x_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
       real_buffer3, real_buffer2)
    if (rank==0) write(*,*) "[DEBUG] Transpose x_from_y_2 with monc (Forward)"
#else
    call transpose_y_to_x_cudecomp_T3(real_buffer3, real_buffer2)
    if (rank==0) write(*,*) "[DEBUG] Transpose x_from_y_2 with cudecomp (Forward), T3"
#endif

    call transpose_and_backward_fft_in_x(current_state, real_buffer2, buffer2, real_buffer1)
    call transpose_and_backward_fft_in_y(current_state, real_buffer1, buffer1, target_data)
  end subroutine perform_backwards_3dfft

  !> Initialises memory for the buffers used in the FFT

  !  source is current_state%p%data
!perform_forward_3dfft------------------
!----------------------------------------
    !--tx_+_for_in_y
    !----------------------------------------   
        ! source : zyx-rr_rr_rr
        !                T1 (t_z_to_y)
        ! fft_in_y_buffer : yxz-rr_rr_rr
        !                FFT
        ! buffer1 : yxz-sc_rr_rr
        !                c2r
        ! real_buffer1 : yxz-2sr_rr-rr
    !--tx+for_in_x
    !----------------------------------------
        !                T2 (t_y_to_x)
        ! fft_in_x_buffer : xzy-rr_rr_2sr
        !                FFT
        ! buffer2 : xzy-sc_rr_2sr
        !                c2r
        ! real_buffer_2 : xzy-2sr_rr_2sr
    !                T3 (t_x_to_y)
    ! real_buffer3 : yxz-2sr-2sr-rr  ?
    !                T3 (t_y_to_z)
    ! target : zyx-rr_2sr_2sr  ----> Source
  ! target is pt (from fftsolver.F90)
!---------------------------------------------

!< do the pressure solve here (tridiagonal_solver in fftsolver.F90)>

! source is pt (from fftsolver.F90)
! Perform_back_3dfft
!---------------------------------------------
    ! source : zyx-rr_2sr_2sr
    !                T3  (t_z_to_y)
    ! real_buffer3 : yxz-2sr-2sr-rr 
    !                T3 (t_y_to_x)
    ! real_buffer_2 : xzy-2sr_rr_2sr

    ! tx_and_back_fft_in_x
    !----------------------------------
          !                r2c
          ! buffer2 : xzy-2sc_rr_2sr
          !                ifft
          ! fft_in_x_buffer : xzy-rr_rr_2sr
          !                T2 (t_x_to_y)
          ! real_buffer_1 : yxz-2sr_rr_rr
    ! tx_and_back_fft_in_y
    !---------------------------------
        !                r2c
        ! buffer1 : yxz-2sc_rr_rr
        !                ifft
        ! fft_in_y_buffer : yxz-rr_rr_rr
        !                T1 (t_y_to_z)
        ! target : zyx-rr-rr-rr
  ! target is current_state%p%data
!---------------------------------------

  ! 
  subroutine initialise_buffers()
    integer, dimension(3) :: nfft_in_y_buffer, nbuffer1, nreal_buffer1, nfft_in_x_buffer, nbuffer2
    integer, dimension(3) :: nreal_buffer2, nreal_buffer3
    integer :: cy1, cy2, cy3, cx1, cx2, cx3
    ! fft_in_y_buffer
    nfft_in_y_buffer(1) = y_from_z_transposition%my_pencil_size(Y_INDEX)
    nfft_in_y_buffer(2) = y_from_z_transposition%my_pencil_size(X_INDEX)
    nfft_in_y_buffer(3) = y_from_z_transposition%my_pencil_size(Z_INDEX)
    
    nbuffer1(1)         = y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1
    nbuffer1(2)         = y_from_z_transposition%my_pencil_size(X_INDEX)
    nbuffer1(3)         = y_from_z_transposition%my_pencil_size(Z_INDEX)

    nreal_buffer1(1)    = (y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1)*2
    nreal_buffer1(2)    = y_from_z_transposition%my_pencil_size(X_INDEX)
    nreal_buffer1(3)    = y_from_z_transposition%my_pencil_size(Z_INDEX)

    nfft_in_x_buffer(1) = x_from_y_transposition%my_pencil_size(X_INDEX)
    nfft_in_x_buffer(2) = x_from_y_transposition%my_pencil_size(Z_INDEX)
    nfft_in_x_buffer(3) = x_from_y_transposition%my_pencil_size(Y_INDEX)

    nbuffer2(1)         = x_from_y_transposition%my_pencil_size(X_INDEX)/2+1
    nbuffer2(2)         = x_from_y_transposition%my_pencil_size(Z_INDEX)
    nbuffer2(3)         = x_from_y_transposition%my_pencil_size(Y_INDEX)

    nreal_buffer2(1)    = (x_from_y_transposition%my_pencil_size(X_INDEX)/2+1)*2
    nreal_buffer2(2)    = x_from_y_transposition%my_pencil_size(Z_INDEX)
    nreal_buffer2(3)    = x_from_y_transposition%my_pencil_size(Y_INDEX)

    nreal_buffer3(1)    = y_from_x_transposition%my_pencil_size(Y_INDEX)
    nreal_buffer3(2)    = y_from_x_transposition%my_pencil_size(X_INDEX)
    nreal_buffer3(3)    = y_from_x_transposition%my_pencil_size(Z_INDEX)

    write(*,'(a20,3i3)') "MONC:, rank:", rank
    write(*,'(a20,3i3)') "fft_in_y(Y,X,Z):",     nfft_in_y_buffer(1), nfft_in_y_buffer(2), nfft_in_y_buffer(3)
    write(*,'(a20,3i3)') "buffer1(Y,X,Z):",       nbuffer1(1), nbuffer1(2), nbuffer1(3)
    write(*,'(a20,3i3)') "real_buffer1(Y,X,Z):", nreal_buffer1(1), nreal_buffer1(2), nreal_buffer1(3)
    write(*,'(a20,3i3)') "fft_in_x(X,Z,Y) A :",     nfft_in_x_buffer(1),nfft_in_x_buffer(2),nfft_in_x_buffer(3)
    write(*,'(a20,3i3)') "buffer2(x,z,y):",      nbuffer2(1), nbuffer2(2), nbuffer2(3)
    write(*,'(a20,3i3)') "real_buffer2(x,z,y):", nreal_buffer2(1), nreal_buffer2(2), nreal_buffer2(3)
    write(*,'(a20,3i3)') "real_buffer3(y,x,z):", nreal_buffer3(1), nreal_buffer3(2), nreal_buffer3(3)

#ifndef CUDECOMP
    allocate(fft_in_y_buffer(y_from_z_transposition%my_pencil_size(Y_INDEX), y_from_z_transposition%my_pencil_size(X_INDEX), &
                        y_from_z_transposition%my_pencil_size(Z_INDEX)), &
            buffer1(y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1, y_from_z_transposition%my_pencil_size(X_INDEX), &
                        y_from_z_transposition%my_pencil_size(Z_INDEX)), &
            real_buffer1((y_from_z_transposition%my_pencil_size(Y_INDEX)/2+1)*2, y_from_z_transposition%my_pencil_size(X_INDEX), &
                        y_from_z_transposition%my_pencil_size(Z_INDEX)), &
            fft_in_x_buffer(x_from_y_transposition%my_pencil_size(X_INDEX), x_from_y_transposition%my_pencil_size(Z_INDEX), &
                        x_from_y_transposition%my_pencil_size(Y_INDEX)), &
            buffer2(x_from_y_transposition%my_pencil_size(X_INDEX)/2+1, x_from_y_transposition%my_pencil_size(Z_INDEX), &
                        x_from_y_transposition%my_pencil_size(Y_INDEX)), &
            real_buffer2((x_from_y_transposition%my_pencil_size(X_INDEX)/2+1)*2, x_from_y_transposition%my_pencil_size(Z_INDEX), &
                        x_from_y_transposition%my_pencil_size(Y_INDEX)), &
            real_buffer3(y_from_x_transposition%my_pencil_size(Y_INDEX), y_from_x_transposition%my_pencil_size(X_INDEX), &
                        y_from_x_transposition%my_pencil_size(Z_INDEX)))
#else
    

    nfft_in_y_buffer(1) = pinfo_T1_y%shape(1)
    nfft_in_y_buffer(2) = pinfo_T1_y%shape(2)
    nfft_in_y_buffer(3) = pinfo_T1_y%shape(3)

    nbuffer1(1)         = pinfo_T1_y%shape(1)/2+1
    nbuffer1(2)         = pinfo_T1_y%shape(2)
    nbuffer1(3)         = pinfo_T1_y%shape(3)

    nreal_buffer1(1)    = pinfo_T2_y%shape(1)
    nreal_buffer1(2)    = pinfo_T2_y%shape(2)
    nreal_buffer1(3)    = pinfo_T2_y%shape(3)

    nfft_in_x_buffer(1) = pinfo_T2_z%shape(1)
    nfft_in_x_buffer(1) = pinfo_T2_z%shape(2)
    nfft_in_x_buffer(1) = pinfo_T2_z%shape(3)

    nbuffer2(1)         = pinfo_T2_z%shape(1)/2+1
    nbuffer2(2)         = pinfo_T2_z%shape(2)
    nbuffer2(3)         = pinfo_T2_z%shape(3)

    nreal_buffer2(1)    = pinfo_T3_z%shape(1)
    nreal_buffer2(2)    = pinfo_T3_z%shape(2)
    nreal_buffer2(3)    = pinfo_T3_z%shape(3)

    nreal_buffer3(1)    = pinfo_T3_y%shape(1)
    nreal_buffer3(2)    = pinfo_T3_y%shape(2)
    nreal_buffer3(3)    = pinfo_T3_y%shape(3)

    write(*,'(a20,3i3)') "cudecomp:, rank:", rank
    write(*,'(a20,3i3)') "fft_in_y(Y,X,Z):",    pinfo_T1_y%shape(1), pinfo_T1_y%shape(2), pinfo_T1_y%shape(3)
    write(*,'(a20,3i3)') "buffer1(Y,X,Z):",      pinfo_T1_y%shape(1)/2+1, pinfo_T1_y%shape(2), pinfo_T1_y%shape(3)
    write(*,'(a20,3i3)') "real_buffer1(Y,X,Z):",pinfo_T2_y%shape(1), pinfo_T2_y%shape(2), pinfo_T2_y%shape(3)
    write(*,'(a20,3i3)') "fft_in_x(X,Z,Y) B :",    pinfo_T2_z%shape(1), pinfo_T2_z%shape(2), pinfo_T2_z%shape(3)
    write(*,'(a20,3i3)') "buffer2(x,z,y):",     pinfo_T2_z%shape(1)/2+1, pinfo_T2_z%shape(2), pinfo_T2_z%shape(3)
    write(*,'(a20,3i3)') "real_buffer2(x,z,y):",pinfo_T3_z%shape(1), pinfo_T3_z%shape(2), pinfo_T3_z%shape(3)
    write(*,'(a20,3i3)') "real_buffer3(y,x,z):",pinfo_T3_y%shape(1), pinfo_T3_y%shape(2), pinfo_T3_y%shape(3)

    !cy1= pinfo_T1_y%shape(1)
    !cy2= pinfo_T1_y%shape(2)
    !cy3= pinfo_T1_y%shape(3)
    !cx1= pinfo_T1_x%shape(1)
    !cx2= pinfo_T1_x%shape(2)
    !cx3= pinfo_T1_x%shape(3)
    !write(*,'(a3,i4,5(a4,i4))') "cx1:", cx1, ", cx2:", cx2, ", cx3:", cx3, "cy1:", cy1, ", cy2:", cy2, ", cy3:", cy3
    allocate(fft_in_y_buffer(pinfo_T1_y%shape(1), pinfo_T1_y%shape(2), pinfo_T1_y%shape(3)))
    allocate(buffer1(pinfo_T1_y%shape(1)/2+1, pinfo_T1_y%shape(2), pinfo_T1_y%shape(3)))
    allocate(real_buffer1((pinfo_T2_y%shape(1)), pinfo_T2_y%shape(2), pinfo_T2_y%shape(3)))

    allocate(fft_in_x_buffer(pinfo_T2_z%shape(1), pinfo_T2_z%shape(2), pinfo_T2_z%shape(3)))
    allocate(buffer2(pinfo_T2_z%shape(1)/2+1, pinfo_T2_z%shape(2), pinfo_T2_z%shape(3)))
    allocate(real_buffer2((pinfo_T3_z%shape(1)), pinfo_T3_z%shape(2), pinfo_T3_z%shape(3)))
    
    allocate(real_buffer3(pinfo_T3_y%shape(1), pinfo_T3_y%shape(2), pinfo_T3_y%shape(3)))
#endif

    call print_shape_r(fft_in_y_buffer, "fft_in_y_buffer")
    call print_shape_c(buffer1, "buffer1")
    call print_shape_r(real_buffer1, "real_buffer1")

    call print_shape_r(fft_in_x_buffer, "fft_in_x_buffer")
    call print_shape_c(buffer2,"buffer2")
    call print_shape_r(real_buffer2,"real_buffer2")
    
    
    call print_shape_r(real_buffer3, "real_buffer3") ! -> target_data
  end subroutine initialise_buffers

  !> Initialises the pencil transpositions, from a pencil in one dimension to that in another
  !! @param current_state The current model state
  !! @param y_distinct_sizes Y sizes per process 
  !! @param x_distinct_sizes X sizes per process
  subroutine initialise_transpositions(current_state, y_distinct_sizes, x_distinct_sizes)
    type(model_state_type), intent(inout) :: current_state
    integer, dimension(:), intent(in) :: y_distinct_sizes, x_distinct_sizes

    type(pencil_transposition) :: z_pencil

    z_pencil=create_initial_transposition_description(current_state)

    ! Transpositions
    y_from_z_transposition=create_transposition(current_state%global_grid, z_pencil, Y_INDEX, y_distinct_sizes, &
         FORWARD, (/ -1 /))   
    x_from_y_transposition=create_transposition(current_state%global_grid, y_from_z_transposition, X_INDEX, &
         x_distinct_sizes, FORWARD, (/ Y_INDEX /))
    y_from_x_transposition=create_transposition(current_state%global_grid, x_from_y_transposition, Y_INDEX, &
         normal_to_extended_process_dim_sizes(x_distinct_sizes), BACKWARD, (/ Y_INDEX, X_INDEX /))    
    z_from_y_transposition=create_transposition(current_state%global_grid, y_from_x_transposition, Z_INDEX, &
         normal_to_extended_process_dim_sizes(y_distinct_sizes), BACKWARD, (/ Y_INDEX, X_INDEX /))  

    y_from_z_2_transposition=create_transposition(current_state%global_grid, z_from_y_transposition, Y_INDEX, &
          normal_to_extended_process_dim_sizes(y_distinct_sizes), FORWARD, (/ Y_INDEX, X_INDEX /))   
    x_from_y_2_transposition=create_transposition(current_state%global_grid, y_from_z_2_transposition, X_INDEX, &
          normal_to_extended_process_dim_sizes(x_distinct_sizes), FORWARD, (/ Y_INDEX, X_INDEX /))
    y_from_x_2_transposition=create_transposition(current_state%global_grid, x_from_y_2_transposition, Y_INDEX, &
         x_distinct_sizes, BACKWARD, (/ Y_INDEX /))
    z_from_y_2_transposition=create_transposition(current_state%global_grid, y_from_x_2_transposition, Z_INDEX, &
          y_distinct_sizes, BACKWARD, (/ -1 /))   
  end subroutine initialise_transpositions

  !> Creates a specific pencil transposition description. It is maybe more a decomposition description, but the main
  !! complexity comes from the transposition from existing decomposition to new decomposition so therefore it is
  !! called transposition. The new pencil decomposition depends not only on the dimension to split on, but also the existing
  !! pencil decomposition. The new decomposed dimension (i.e. the existing pencil dimension) * other local dimensions is used
  !! as the sending size, receiving though requires knowledge about the data size on the source process so others will send
  !! their pencil dimension size to this process.
  !! @param new_pencil_dim The dimension to use as the new pencil decomposition
  !! @param existing_pencil_dim The dimension used in the current decomposition
  !! @param existing_pencil_process_layout Number of processes per dimension for the current decomposition
  !! @param existing_my_location The current processes block location per dimension for the current decomposition
  !! @param existing_pencil_size Pencil size per dimension for the current decomposition
  !! @param process_dim_sizes Sizes of the pencil dimension from other processes that is used to calculate receive count
  !! @param direction Whether we are transposing forwards or backwards, backwards is just an inverse
  !! @param extended_dimensions The dimensions that this process extends from n to (n/2+1)*2 (i.e. result of fft complex->real)
  type(pencil_transposition) function create_transposition(global_grid, existing_transposition, new_pencil_dim,&
       process_dim_sizes, direction, extended_dimensions)
    type(global_grid_type), intent(inout) :: global_grid
    type(pencil_transposition), intent(in) :: existing_transposition
    integer, dimension(:), intent(in) :: process_dim_sizes
    integer, intent(in) :: new_pencil_dim, direction, extended_dimensions(:)

    create_transposition%process_decomposition_layout=determine_pencil_process_dimensions(&
         new_pencil_dim, existing_transposition%dim, existing_transposition%process_decomposition_layout)

    create_transposition%my_process_location=determine_my_pencil_location(new_pencil_dim, &
         existing_transposition%dim, existing_transposition%my_process_location)

    create_transposition%my_pencil_size=determine_pencil_size(new_pencil_dim, create_transposition%process_decomposition_layout,&
         create_transposition%my_process_location, existing_transposition, global_grid, extended_dimensions)

    allocate(create_transposition%send_dims(3, create_transposition%process_decomposition_layout(existing_transposition%dim)), &
           create_transposition%recv_dims(3, create_transposition%process_decomposition_layout(existing_transposition%dim)))
    if (direction == FORWARD) then            
      call determine_my_process_sizes_per_dim(existing_transposition%dim, &
           existing_transposition%my_pencil_size, create_transposition%process_decomposition_layout, &
           global_grid, extended_dimensions, create_transposition%send_dims)
      call determine_matching_process_dimensions(new_pencil_dim, existing_transposition%dim, process_dim_sizes, &
           create_transposition%my_pencil_size, create_transposition%process_decomposition_layout, create_transposition%recv_dims)
    else
      call determine_my_process_sizes_per_dim(new_pencil_dim, create_transposition%my_pencil_size, &
           existing_transposition%process_decomposition_layout, global_grid, extended_dimensions, create_transposition%recv_dims)
      call determine_matching_process_dimensions(existing_transposition%dim, new_pencil_dim, process_dim_sizes, &
           existing_transposition%my_pencil_size, existing_transposition%process_decomposition_layout, &
           create_transposition%send_dims)
    end if

    allocate(create_transposition%send_sizes(size(create_transposition%send_dims, 2)), &
         create_transposition%send_offsets(size(create_transposition%send_sizes)), &
         create_transposition%recv_sizes(size(create_transposition%recv_dims, 2)), &
         create_transposition%recv_offsets(size(create_transposition%recv_sizes)))

    call concatenate_dimension_sizes(create_transposition%send_dims, create_transposition%send_sizes)
    call determine_offsets_from_size(create_transposition%send_sizes, create_transposition%send_offsets)

    call concatenate_dimension_sizes(create_transposition%recv_dims, create_transposition%recv_sizes)
    call determine_offsets_from_size(create_transposition%recv_sizes, create_transposition%recv_offsets)
    create_transposition%dim=new_pencil_dim
  end function create_transposition

  !> Performs the transposition and forward FFT in the y dimension then converts back to real numbers. The Y size is
  !! (n/2+1)*2 due to the complex to real transformation after the FFT.
  !! @param current_state The current model state
  !! @param source_data Input buffer, Z pencil oriented z,y,x
  !! @param buffer Complex buffer which the FFT writes into
  !! @param real_buffer Output buffer, Y pencil, oriented y,x,z
  subroutine transpose_and_forward_fft_in_y(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(out) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer

    ! Transpose globally from Z pencil to Y pencil
#ifndef CUDECOMP
    call transpose_to_pencil(y_from_z_transposition, (/Z_INDEX, Y_INDEX, X_INDEX/), dim_y_comm, FORWARD, &
       source_data, fft_in_y_buffer)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_z with monc (Forward)"
#else
    call transpose_z_to_y_cudecomp_T1(source_data, fft_in_y_buffer)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_z with cudecomp (Forward), T1"
#endif
    
    call perform_r2c_fft(fft_in_y_buffer, buffer, y_from_z_transposition%my_pencil_size(Y_INDEX), &
         y_from_z_transposition%my_pencil_size(X_INDEX) * y_from_z_transposition%my_pencil_size(Z_INDEX), 1)
    call convert_complex_to_real(buffer, real_buffer)
  end subroutine transpose_and_forward_fft_in_y

#ifdef CUDECOMP
  subroutine transpose_z_to_y_cudecomp(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat, n1, n2, n3, n1_d, n2_d, n3_d
    ! x-pointer to gpu cudecomp data
    n1=size(source_data,1); n2=size(source_data,2); n3=size(source_data,3)
    n1_d=size(data_cudecomp_T1_x_d,1); n2_d=size(data_cudecomp_T1_x_d,2); n3_d=size(data_cudecomp_T1_x_d,3)
    write(*,'(a44,6i4)') "[DEBUG] in - n1, n2, n3, n1_d, n2_d, n3_d:",n1, n2, n3, n1_d, n2_d, n3_d
    data_cudecomp_T1_x_d = source_data
    ! Use x-y instead of z-y because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeXToY(handle_cudecomp_T1, grid_desc_cudecomp_T1, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)

    n1=size(dest_data,1); n2=size(dest_data,2); n3=size(dest_data,3)
    n1_d=size(data_cudecomp_T1_y_d,1); n2_d=size(data_cudecomp_T1_y_d,2); n3_d=size(data_cudecomp_T1_y_d,3)
    write(*,'(a44,6i4)') "[DEBUG] out - n1, n2, n3, n1_d, n2_d, n3_d:",n1, n2, n3, n1_d, n2_d, n3_d
    ! y-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T1_y_d
  end subroutine transpose_z_to_y_cudecomp    

  subroutine transpose_z_to_y_cudecomp_T1(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! x-pointer to gpu cudecomp data
    data_cudecomp_T1_x_d = source_data
    ! Use y-z instead of y-x because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeXToY(handle_cudecomp_T1, grid_desc_cudecomp_T1, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! z-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T1_y_d
  end subroutine transpose_z_to_y_cudecomp_T1

  subroutine transpose_z_to_y_cudecomp_T3(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:) :: dest_data !, intent(inout) :: dest_data
    integer :: istat
    ! x-pointer to gpu cudecomp data
    data_cudecomp_T3_x_d = source_data
    ! Use y-z instead of y-x because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeXToY(handle_cudecomp_T1, grid_desc_cudecomp_T3, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! z-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T3_y_d
  end subroutine transpose_z_to_y_cudecomp_T3
  
  subroutine transpose_y_to_x_cudecomp_T2(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! y-pointer to gpu cudecomp data
    data_cudecomp_T2_y_d = source_data
    ! Use y-z instead of y-x because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeYToZ(handle_cudecomp_T1, grid_desc_cudecomp_T2, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! z-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T2_z_d
  end subroutine transpose_y_to_x_cudecomp_T2    

  subroutine transpose_y_to_x_cudecomp_T3(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! y-pointer to gpu cudecomp data
    data_cudecomp_T3_y_d = source_data
    ! Use y-z instead of y-x because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeYToZ(handle_cudecomp_T1, grid_desc_cudecomp_T3, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! z-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T3_z_d
  end subroutine transpose_y_to_x_cudecomp_T3    

  subroutine transpose_x_to_y_cudecomp_T2(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! z-pointer to gpu cudecomp data
    data_cudecomp_T2_z_d = source_data
    ! Use z-y instead of x-y because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeZToY(handle_cudecomp_T1, grid_desc_cudecomp_T2, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! y-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T2_y_d
  end subroutine transpose_x_to_y_cudecomp_T2

  subroutine transpose_x_to_y_cudecomp_T3(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! z-pointer to gpu cudecomp data
    data_cudecomp_T3_z_d = source_data
    ! Use z-y instead of x-y because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeZToY(handle_cudecomp_T1, grid_desc_cudecomp_T3, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! y-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T3_y_d
  end subroutine transpose_x_to_y_cudecomp_T3

  subroutine transpose_y_to_z_cudecomp_T1(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! y-pointer to gpu cudecomp data
    data_cudecomp_T1_y_d = source_data
    ! Use y-x instead of y-z because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeYToX(handle_cudecomp_T1, grid_desc_cudecomp_T1, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! x-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T1_x_d
  end subroutine transpose_y_to_z_cudecomp_T1   

  subroutine transpose_y_to_z_cudecomp_T3(source_data, dest_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: dest_data
    integer :: istat
    ! y-pointer to gpu cudecomp data
    data_cudecomp_T3_y_d = source_data
    ! Use y-x instead of y-z because monc z direction is equivalent to cudecomp x direction
    istat = cudecompTransposeYToX(handle_cudecomp_T1, grid_desc_cudecomp_T3, data_cudecomp_d, data_cudecomp_d, CUDECOMP_WORKSPACE_AUTO, CUDECOMP_DOUBLE)
    ! x-pointer to gpu cudecomp data
    dest_data = data_cudecomp_T3_x_d
  end subroutine transpose_y_to_z_cudecomp_T3
#endif

  !> Performs the backwards FFT in X and then transposes to Y pencil. The FFT requires complex numbers which are converted to real, 
  !! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
  !! size is n+2, complex number size is n/2+1 and we get n reals out.
  !! @param current_state The current model state
  !! @param source_data Input buffer, X pencil oriented x,z,y
  !! @param buffer Complex buffer which is fed into the FFT
  !! @param real_buffer Output buffer, Y pencil, oriented y,x,z
  subroutine transpose_and_backward_fft_in_x(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(inout) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

    call convert_real_to_complex(source_data, buffer)
    call perform_c2r_fft(buffer, fft_in_x_buffer, x_from_y_2_transposition%my_pencil_size(X_INDEX)-2, &
         x_from_y_2_transposition%my_pencil_size(Y_INDEX) * x_from_y_2_transposition%my_pencil_size(Z_INDEX), 2) 

    ! Transpose globally from X pencil to Y pencil
#ifndef CUDECOMP
    call transpose_to_pencil(y_from_x_2_transposition, (/X_INDEX, Z_INDEX, Y_INDEX/), dim_x_comm, BACKWARD, &
       fft_in_x_buffer, real_buffer)
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_x_2 with monc (Backward)"
#else
    call transpose_x_to_y_cudecomp_T2(fft_in_x_buffer, real_buffer) !< this doesn't work!!!
    if (rank==0) write(*,*) "[DEBUG] Transpose y_from_x_2 with cudecomp (Backward), T2"
#endif
  end subroutine transpose_and_backward_fft_in_x

  !> Performs the transposition and forward FFT in the x dimension. After the FFT the complex space is converted back into
  !! real numbers. The X size is (n/2+1)*2 due to this transformation.
  !! @param current_state The current model state
  !! @param buffer1 Input buffer, Y pencil after the Y dimension FFT oriented y,x,z
  !! @param buffer Complex buffer which results from the FFT
  !! @param buffer2 Output buffer, X pencil after this X FFT, oriented x,z,y
  subroutine transpose_and_forward_fft_in_x(current_state, buffer1, buffer2, buffer3)
    type(model_state_type), target, intent(inout) :: current_state
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:),  contiguous, pointer, intent(out) :: buffer2
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: buffer1, buffer3

    ! Go from global Y pencil to global X pencil
#ifndef CUDECOMP
    call transpose_to_pencil(x_from_y_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_x_comm, FORWARD, &
       buffer1, fft_in_x_buffer)   
    if (rank==0) write(*,*) "[DEBUG] Transpose x_from_y with monc (Forward)"
#else
    call transpose_y_to_x_cudecomp_T2(buffer1, fft_in_x_buffer) !< not working
    if (rank==0) write(*,*) "[DEBUG] Transpose x_from_y with cudecomp (Forward), T2"
#endif
    
    call perform_r2c_fft(fft_in_x_buffer, buffer2, x_from_y_transposition%my_pencil_size(X_INDEX), &
         x_from_y_transposition%my_pencil_size(Y_INDEX) * x_from_y_transposition%my_pencil_size(Z_INDEX), 3)

    call convert_complex_to_real(buffer2, buffer3)
  end subroutine transpose_and_forward_fft_in_x

  !> Performs the backwards FFT in Y and then transposes to Z pencil. The FFT requires complex numbers which are converted to real, 
  !! so the this real to complex operation is performed first. If n is the logical size of the FFT row, then the input
  !! size is n+2, complex number size is n/2+1 and we get n reals out.
  !! @param current_state The current model state
  !! @param source_data Input buffer, Y pencil oriented y,x,z
  !! @param buffer Complex buffer which is fed into the FFT
  !! @param real_buffer Output buffer, Z pencil, oriented z,y,x
  subroutine transpose_and_backward_fft_in_y(current_state, source_data, buffer, real_buffer)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:),  intent(out) :: real_buffer
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: buffer

    call convert_real_to_complex(source_data, buffer)
   
    call perform_c2r_fft(buffer, fft_in_y_buffer,  y_from_x_2_transposition%my_pencil_size(Y_INDEX)-2, &
         y_from_x_2_transposition%my_pencil_size(X_INDEX) * y_from_x_2_transposition%my_pencil_size(Z_INDEX), 4)

    ! Go from global Y pencil to global Z pencil
#ifndef CUDECOMP
    call transpose_to_pencil(z_from_y_2_transposition, (/Y_INDEX, X_INDEX, Z_INDEX/), dim_y_comm, BACKWARD, &
       fft_in_y_buffer, real_buffer)
    if (rank==0) write(*,*) "[DEBUG] Transpose z_from_y_2 with monc (Backward)"
#else
    call transpose_y_to_z_cudecomp_T1(fft_in_y_buffer, real_buffer)
    if (rank==0) write(*,*) "[DEBUG] Transpose z_from_y_2 with cudecomp (Backward), T1"
#endif
  end subroutine transpose_and_backward_fft_in_y

  !> Transposes globally to a new pencil decomposition. This goes from the source dimensions a,b,c to b,c,a (forwards) or c,a,b
  !! (backwards.) It requires multiple steps, first the local data is transposed to c,b,a regardless of direction. 
  !! then it is communicated via alltoall, each process then assembles its own b,c,a or c,a,b data via contiguising 
  !! across blocks as the data layout is nonlinear.
  !! @param transposition_description Description of the transposition
  !! @param source_dims Dimensions of the current pencil that we wish to transpose from, will go from abc to bca
  !! @param communicator The MPI communicator associated with the group of processes who will swap data
  !! @param direction Whether this is going forwards or backwards, it makes a difference to the data arrangement
  !! @param source_data Source data (abc)
  !! @param target_data Target data (bca)
  subroutine transpose_to_pencil(transposition_description, source_dims, communicator, direction, source_data, target_data)
    type(pencil_transposition), intent(in) :: transposition_description
    integer, intent(in) :: source_dims(3), communicator, direction
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_data

    integer :: ierr
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: real_temp
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: real_temp2
    
    
    allocate(real_temp(size(source_data,3), size(source_data,2), size(source_data,1)), &
         real_temp2(product(transposition_description%my_pencil_size)+1))

    call rearrange_data_for_sending(real_source=source_data, real_target=real_temp)    

    call mpi_alltoallv(real_temp, transposition_description%send_sizes, transposition_description%send_offsets, &
         PRECISION_TYPE, real_temp2, transposition_description%recv_sizes, transposition_description%recv_offsets, &
         PRECISION_TYPE, communicator, ierr)
    call contiguise_data(transposition_description, (/source_dims(3), source_dims(2), source_dims(1)/), direction, &
         source_real_buffer=real_temp2, target_real_buffer=target_data)
    deallocate(real_temp, real_temp2)   
  end subroutine transpose_to_pencil  

  !> Contiguises from c,b,a to b,c,a (forwards) or c,a,b (backwards) where these are defined by the source_dims argument. 
  !! It is not as simple as just swapping the required dimensions, as this is after the mpi alltoall and each block lies 
  !! after the previous block running sequentially in a.
  !! @param transposition_description Transposition descriptor
  !! @param source_dims Representation a,b,c of source data, will contiguise to b,a,c
  !! @param direction Whether we wish to contiguise forwards or backwards
  !! @param source_real_buffer Source real data to transform
  !! @param target_real_buffer Target real data which is the result of the operation
  subroutine contiguise_data(transposition_description, source_dims, direction, source_real_buffer, target_real_buffer)
    integer, intent(in) :: source_dims(3), direction
    type(pencil_transposition), intent(in) :: transposition_description
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: source_real_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: target_real_buffer
    
    integer :: number_blocks, i, j, k, n, index_prefix, index_prefix_dim, block_offset, source_index

    number_blocks=size(transposition_description%recv_sizes)
    index_prefix=0
    block_offset=0
    index_prefix_dim=merge(2,1, direction == FORWARD)
    do i=1,number_blocks
      if (i .ge. 2) then
        index_prefix=index_prefix+transposition_description%recv_dims(source_dims(index_prefix_dim), i-1)
        block_offset=block_offset+transposition_description%recv_sizes(i-1)
      end if
      !Transformation is either cba -> bca (forward) or cab (backwards)       
      do j=1, transposition_description%recv_dims(source_dims(3), i) ! a
        do k=1, transposition_description%recv_dims(source_dims(1), i) ! c
          do n=1, transposition_description%recv_dims(source_dims(2), i) ! b
            source_index=block_offset+(j-1)* transposition_description%recv_dims(source_dims(1), i)* &
                 transposition_description%recv_dims(source_dims(2), i)+ (n-1)* &
                 transposition_description%recv_dims(source_dims(1), i)+k            
            if (direction == FORWARD) then
              target_real_buffer(index_prefix+n, k, j)=source_real_buffer(source_index) ! bca
            else
              target_real_buffer(index_prefix+k, j, n)=source_real_buffer(source_index) ! cab
            end if
          end do
        end do
      end do      
    end do    
  end subroutine contiguise_data

#ifdef CUDA
  !> Cuda initialisation forward transform
  !! @param plan fft plan for specific transform
  !! @param n Size of transform
  !! @param istat status flag for cuda function
  !! @param ncol number of batched transforms to carry out
  subroutine fft_r2c_gpu_init(plan, n, istat, ncol)
    integer, intent(inout) :: plan
    integer, intent(in) :: n
    integer, intent(inout) :: istat
    integer, intent(in) :: ncol
    istat = cufftPlan1D(plan,n,CUFFT_D2Z,ncol)
  end subroutine fft_r2c_gpu_init

  !> Cuda initialisation backward transform
  !! @param plan fft plan for specific transform
  !! @param n Size of transform
  !! @param istat status flag for cuda function
  !! @param ncol number of batched transforms to carry out
  subroutine fft_c2r_gpu_init(plan, n, istat,ncol)
    integer, intent(inout) :: plan
    integer, intent(in) :: n
    integer, intent(inout) :: istat
    integer, intent(in) :: ncol
    istat = cufftPlan1D(plan,n,CUFFT_Z2D,ncol)
  end subroutine fft_c2r_gpu_init

  !> Carry out forward fft with cuda
  !! @param in Data to transform (on host)
  !! @param out Transformed data (on host)
  !! @param nt Size of transform
  !! @param plan fft plan for specific transform
  !! @param ncol number of batched transforms to carry out
  subroutine fft_r2c_gpu(in, out, nt, plan)
    real(kind=DEFAULT_PRECISION), intent(inout) :: in(:,:,:)
    complex*16, intent(inout) :: out(:,:,:)
    integer, intent(in) :: nt
    integer, intent(inout) :: plan
    real(kind=DEFAULT_PRECISION), allocatable, device, dimension(:,:,:) :: in_d
    complex*16, allocatable, device, dimension(:,:,:) :: out_d
    integer :: istat, n2, n3
    n2=size(in,2)
    n3=size(in,3)
    allocate(in_d(nt,n2,n3), out_d(nt/2+1,n2,n3))
    in_d = in
    istat = cufftExecD2Z(plan, in_d, out_d)
    out = out_d
  end subroutine fft_r2c_gpu

  !> Carry out backward fft with cuda
  !! @param in Data to transform (on host)
  !! @param out Transformed data (on host)
  !! @param nt Size of transform
  !! @param plan fft plan for specific transform
  !! @param ncol number of batched transforms to carry out
  subroutine fft_c2r_gpu(in, out, nt, plan)
    integer, intent(inout) :: plan
    complex*16, intent(inout) :: in(:,:,:)
    real(kind=DEFAULT_PRECISION), intent(inout) :: out(:,:,:)
    integer, intent(in) :: nt
    complex*16, allocatable, device, dimension(:,:,:) :: in_d
    real(kind=DEFAULT_PRECISION), allocatable, device, dimension(:,:,:) :: out_d
    integer :: istat, n2, n3
    n2=size(in,2)
    n3=size(in,3)
    allocate(in_d(nt/2+1,n2,n3), out_d(nt,n2,n3))
    in_d = in
    istat = cufftExecZ2D(plan, in_d, out_d)
    out = out_d/nt
  end subroutine fft_c2r_gpu
#endif

  !> Actually performs a forward real to complex FFT
  !! @param source_data Source (real) data in the time domain
  !! @param transformed_data Resulting complex data in the frequency domain
  !! @param row_size Number of elements for each FFT
  !! @param num_rows The number of FFTs to perform on the next data elements in the source_data
  !! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
  subroutine perform_r2c_fft(source_data, transformed_data, row_size, num_rows, plan_id)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
    integer, intent(in) :: row_size, num_rows, plan_id
    integer :: i, j, istat, ncols, ncols2, ncols3
    real(kind=DEFAULT_PRECISION) :: tstart, tstop
    ncols2=size(source_data,2)
    ncols3=size(source_data,3)
    ncols=size(source_data,2)* size(source_data,3)

    tstart = MPI_Wtime()
#if !defined(GPU)
    call fftn_init(row_size)

    do i=1,size(source_data,3)
      do j=1,size(source_data,2)
        call fftn_r2c(source_data(:,j,i),transformed_data(:,j,i),row_size)
      enddo
    enddo

    call fftn_finalise()
#else
    ! create plan if not already created
    if (.not. plan_defined(plan_id)) then
      call fft_r2c_gpu_init(plans(plan_id),row_size, istat, ncols)
      plan_defined(plan_id) = .true.
    end if
    call fft_r2c_gpu(source_data(:,1:ncols2,1:ncols3),transformed_data(:,1:ncols2,1:ncols3),row_size, plans(plan_id))
#endif

    tstop = mpi_wtime()

    nforward = nforward +1
    tforward = tforward + (tstop-tstart)
  end subroutine perform_r2c_fft

  !> Performs the complex to real (backwards) FFT
  !! @param source_data Source (complex) data in the frequency domain
  !! @param transformed_data Resulting real data in the time domain
  !! @param row_size Number of elements for each FFT
  !! @param num_rows The number of FFTs to perform on the next data elements in the source_data
  !! @param plan_id Id number of the plan that tracks whether we need to create it or can reuse the existing one
  subroutine perform_c2r_fft(source_data, transformed_data, row_size, num_rows, plan_id)
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(inout) :: source_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), contiguous, pointer, intent(inout) :: transformed_data
    integer, intent(in) :: row_size, num_rows, plan_id
    integer :: i,j, istat, ncols, ncols2,ncols3
    real(kind=DEFAULT_PRECISION) :: tstart, tstop
    ncols2=size(source_data,2)
    ncols3=size(source_data,3)
    ncols=size(source_data,2)* size(source_data,3)

    tstart = MPI_Wtime()
#if !defined(GPU)
    call fftn_init(row_size)

    do i=1,size(source_data,3)
      do j=1,size(source_data,2)
        call fftn_c2r(source_data(:,j,i),transformed_data(:,j,i),row_size)
      enddo
    enddo

    call fftn_finalise()
#else
    ! create plan if not already created
    if (.not. plan_defined(plan_id)) then
      call fft_c2r_gpu_init(plans(plan_id),row_size, istat,ncols)
      plan_defined(plan_id) = .true.
    end if
    call fft_c2r_gpu(source_data(:,1:ncols2,1:ncols3),transformed_data(:,1:ncols2,1:ncols3),row_size, plans(plan_id))
#endif
    tstop = mpi_wtime()

    nback = nback +1
    tback = tback + (tstop-tstart)

  end subroutine perform_c2r_fft

  !> Rearranges data for sending, transposing a,b,c into c,b,a . This is done as alltoall splits on dimension c
  !! so to go from one pencil to another we assume here that a is the existing pencil as it is contiguous
  !! @param real_source Source data to transpose from
  !! @param real_target Target data to transpose to
  subroutine rearrange_data_for_sending(real_source, real_target)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: real_source
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: real_target

    integer :: i

    do i=1, size(real_source,2)
      real_target(:,i,:)=transpose(real_source(:,i,:))
    end do
  end subroutine rearrange_data_for_sending

  !> Determines the number of elements to on my process per dimension which either need to be sent to (forwards transformation) or
  !! received from (backwards) each target process (in the row or column)
  !! This depends on the existing pencil decomposition, as effectively we are breaking that contigulity and
  !! decomposing it into n blocks in that dimension now (provided by new_pencil_procs_per_dim)
  !! @param existing_pencil_dim The pencil dimension that we are transforming from
  !! @param existing_pencil_size Existing pencil decomposition sizes per dimension
  !! @param new_pencil_procs_per_dim For the target decomposition the number of processes per dimension
  !! @param global_grid Description of the global grid which we use for sizing information
  !! @param extended_dimensions List of dimensions where we extend from n to n+2 (i.e. result of FFT complex-> real transformation)
  subroutine determine_my_process_sizes_per_dim(existing_pencil_dim, existing_pencil_size, new_pencil_procs_per_dim, &
       global_grid, extended_dimensions, specific_sizes_per_dim)
    integer, intent(in) :: existing_pencil_dim, existing_pencil_size(:), new_pencil_procs_per_dim(:), extended_dimensions(:)
    type(global_grid_type), intent(inout) :: global_grid
    integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim

    integer :: i, split_size, split_remainder, j, s

    do i=1,3
      if (i == existing_pencil_dim) then
        s=global_grid%size(i)
        if (is_extended_dimension(i, extended_dimensions)) s=s+2
        split_size = s / new_pencil_procs_per_dim(i)
        split_remainder = s - split_size * new_pencil_procs_per_dim(i)
        do j=1,new_pencil_procs_per_dim(existing_pencil_dim)
          specific_sizes_per_dim(i,j)=merge(split_size+1, split_size, j .le. split_remainder)
        end do        
      else
        specific_sizes_per_dim(i,:) = existing_pencil_size(i)
      end if
    end do    
  end subroutine determine_my_process_sizes_per_dim

  !> Simple helper function to deduce send or receive offsets from the sizes
  !! @param source_sizes Sizes that we are using to build the offsets
  subroutine determine_offsets_from_size(source_sizes, determined_offsets)
    integer, intent(in) :: source_sizes(:)
    integer, dimension(:), intent(inout) :: determined_offsets

    integer :: i

    determined_offsets(1)=0
    do i=2,size(source_sizes)
      determined_offsets(i)=determined_offsets(i-1)+source_sizes(i-1)
    end do    
  end subroutine determine_offsets_from_size  

  !> Determines the number of processes in each dimension for the target decomposition. This depends heavily
  !! on the existing decomposition, as we basically contiguise our pencil dimension and decompose the existing
  !! pencil dimension. The third dimension remains unchanged
  !! @param new_pencil_dim New pencil dimension
  !! @param existing_pencil_dim Current decomposition pencil dimension
  !! @param existing_pencil_procs Current decomposition process layout
  function determine_pencil_process_dimensions(new_pencil_dim, existing_pencil_dim, existing_pencil_procs)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_pencil_procs(3)
    integer :: determine_pencil_process_dimensions(3)

    integer :: i

    do i=1,3
      if (i == new_pencil_dim) then
        determine_pencil_process_dimensions(i)=1
      else if (i == existing_pencil_dim) then
        determine_pencil_process_dimensions(i)=existing_pencil_procs(new_pencil_dim)
      else
        determine_pencil_process_dimensions(i)=existing_pencil_procs(i)
      end if
    end do    
  end function determine_pencil_process_dimensions

  !> Determines my location for each dimension in the new pencil decomposition. I.e. which block I am
  !! operating on
  !! @param new_pencil_dim New pencil decomposition dimension
  !! @param existing_pencil_dim Current pencil dimension
  !! @param existing_locations Location for the current decomposition
  function determine_my_pencil_location(new_pencil_dim, existing_pencil_dim, existing_locations)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, existing_locations(3)
    integer :: determine_my_pencil_location(3)

    integer :: i

    do i=1,3
      if (i == new_pencil_dim) then
        determine_my_pencil_location(i)=1
      else if (i == existing_pencil_dim) then
        determine_my_pencil_location(i)=existing_locations(new_pencil_dim)
      else
        determine_my_pencil_location(i)=existing_locations(i)
      end if
    end do    
  end function determine_my_pencil_location  

  !> Concatenates sizes in multiple dimensions for each target process (in a row or column) into a product of
  !! that. This represents all the dimension sizes per process
  !! @param dims The sizes, per dimension and per process that we will fold into target process
  subroutine concatenate_dimension_sizes(dims, concatenated_dim_sizes)
    integer, dimension(:,:), intent(in) :: dims
    integer, dimension(:), intent(inout) :: concatenated_dim_sizes

    integer :: i

    do i=1,size(dims, 2)
      concatenated_dim_sizes(i)=product(dims(:,i))
    end do    
  end subroutine concatenate_dimension_sizes  

  !> Determines the sizes per dimension on the matching process either to receive from (forward transposition) or send to 
  !! (backwards transposition) each source process. Not only does this depend on the
  !! my pencil sizes, but it also depends on the amount of data that the source process has to send over
  !! @param new_pencil_dim The dimension for the new pencil decomposition
  !! @param existing_pencil_dim Dimension for the existing pencil decomposition
  !! @param proc_sizes Size of dimension on the source processes (index in array corresponds to source PID)
  !! @param my_pencil_size My (new) pencil size per dimension
  !! @param pencil_processes_per_dim The process layout per dimension
  subroutine determine_matching_process_dimensions(new_pencil_dim, existing_pencil_dim, proc_sizes, &
       my_pencil_size, pencil_processes_per_dim, specific_sizes_per_dim)
    integer, intent(in) :: new_pencil_dim, existing_pencil_dim, proc_sizes(:), my_pencil_size(:), pencil_processes_per_dim(:)
    integer, dimension(:,:), intent(inout) :: specific_sizes_per_dim

    integer :: i, j

    do i=1,pencil_processes_per_dim(existing_pencil_dim)
      do j=1,3
        if (j==new_pencil_dim) then
          specific_sizes_per_dim(j, i)=proc_sizes(i)
        else
          specific_sizes_per_dim(j, i)=my_pencil_size(j)
        end if
      end do      
    end do    
  end subroutine determine_matching_process_dimensions

  !> Creates an initial transposition representation of the Z pencil that MONC is normally decomposed in. This is then
  !! fed into the create transposition procedure which will generate transpositions to other pencils
  !! @param current_state The current model state
  type(pencil_transposition) function create_initial_transposition_description(current_state)
    type(model_state_type), intent(inout) :: current_state

    create_initial_transposition_description%dim=Z_INDEX
    create_initial_transposition_description%process_decomposition_layout=current_state%parallel%dim_sizes
    create_initial_transposition_description%my_process_location=current_state%parallel%my_coords
    create_initial_transposition_description%my_pencil_size=current_state%local_grid%size
  end function create_initial_transposition_description

  !> Deduces the size of my (local) pencil based upon the new decomposition. This depends heavily on the current
  !! pencil decomposition, the new pencil dimension is the global size, the existing pencil dimension becomes 
  !! decomposed based on the number of processes in that dimension. The third dimension remains unchanged.
  !! @param new_pencil_dim Dimension for the new pencil decomposition
  !! @param pencil_process_layout The processes per dimension layout for the new decomposition
  !! @param my_pencil_location My location in the block layout
  !! @param existing_pencil_dim Current decomposition dimension
  !! @param existing_pencil_size Current decomposition sizes
  !! @param global_grid Description of the global grid which we use for sizing information
  !! @param extended_dimensions List of dimensions where we extend from n to n+2 (i.e. result of FFT complex-> real transformation)
  function determine_pencil_size(new_pencil_dim, pencil_process_layout, my_pencil_location, existing_transposition,&
       global_grid, extended_dimensions)

    type(pencil_transposition), intent(in) :: existing_transposition
    integer, intent(in) :: new_pencil_dim, pencil_process_layout(3), my_pencil_location(3), extended_dimensions(:)
    type(global_grid_type), intent(inout) :: global_grid
    integer :: determine_pencil_size(3)

    integer :: i, split_size, split_remainder, s

    do i=1,3
      if (i == new_pencil_dim) then
        if (is_extended_dimension(i, extended_dimensions)) then
          ! If complex and Y dim then /2+1 for the global size
          determine_pencil_size(i)=(global_grid%size(new_pencil_dim)/2+1)*2
        else
          determine_pencil_size(i)=global_grid%size(new_pencil_dim)
        end if
      else if (i == existing_transposition%dim) then
        s=global_grid%size(i)
        ! If complex and Y dim then use s/2+1 for the size to split
        if (is_extended_dimension(i, extended_dimensions)) s=(s/2+1)*2
        split_size=s/pencil_process_layout(i)
        split_remainder=s - split_size * pencil_process_layout(i)
        determine_pencil_size(i)=merge(split_size+1, split_size, my_pencil_location(i)+1 .le. split_remainder)        
      else        
        determine_pencil_size(i)=existing_transposition%my_pencil_size(i)       
      end if      
    end do    
  end function determine_pencil_size

  !> Determines whether or not the specific dimension is in the list of extended dimensions
  !! @param dimension The dimension to test for
  !! @param extended_dimensions Array of dimensions that will be searched
  !! @returns Whether the dimension is found in the array
  logical function is_extended_dimension(dimension, extended_dimensions)
    integer, intent(in) :: dimension, extended_dimensions(:)

    integer :: i
    do i=1,size(extended_dimensions)
      if (extended_dimensions(i) == dimension) then
        is_extended_dimension=.true.
        return
      end if      
    end do
    is_extended_dimension=.false.
  end function is_extended_dimension
  
  !> Transforms real process dimension sizes into their real after FFT complex->real transformation. The way this works is that
  !! it goes from n to (n/2+1)*2 numbers which is distributed amongst the processes deterministically
  !! @param process_dim_sizes Real process dimension sizes
  !! @returns The extended process dimension sizes
  function normal_to_extended_process_dim_sizes(process_dim_sizes)
    integer, dimension(:), intent(in) :: process_dim_sizes
    integer, dimension(size(process_dim_sizes)) :: normal_to_extended_process_dim_sizes

    integer :: temp_total, split_size, remainder

    temp_total=(sum(process_dim_sizes) /2 + 1) * 2
    split_size=temp_total/size(process_dim_sizes)
    remainder=temp_total - split_size*size(process_dim_sizes)

    normal_to_extended_process_dim_sizes=split_size
    normal_to_extended_process_dim_sizes(1:remainder)=split_size+1    
  end function normal_to_extended_process_dim_sizes

  !> Converts complex representation to its real data counterpart and is called after each forward FFT. 
  !! After a r2c FFT, there are n/2+1 complex numbers - which means that there will be more real numbers in Fourier space
  !! than are provided into the forward FFT call (due to the extra +1). Note that the real size n will always be complex size * 2
  !! This always unpacks the complex dimension in the first dimension
  !! @param complex_data Complex data in Z,Y,X orientation to be unpacked into its real representation
  !! @param real_data The real representation is written into here
  subroutine convert_complex_to_real(complex_data, real_data)
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), intent(in) :: complex_data
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(out) :: real_data

    integer :: i, j, k

    do i=1,size(real_data,3)
      do j=1,size(real_data,2)
        do k=1,size(real_data,1),2
          real_data(k,j,i)=real(real(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
          real_data(k+1,j,i)=real(aimag(complex_data((k+1)/2,j,i)), kind=DEFAULT_PRECISION)
        end do
      end do
    end do
  end subroutine convert_complex_to_real

  !> Converts reals into their complex representation, this is called for backwards FFTs as we need to feed in complex numbers
  !! to force FFTW to do a backwards. It is a relatively simple transformation, as n goes into n/2 complex numbers and as this
  !! is the result of the `convert_complex_to_real` procedure, n always divides evenly.
  !! This is always applied to the first dimension of the real data
  !! @param real_data The source real data to pack into the complex data, it is oriented Z,Y,X
  !! @param complex_data Target complex data which the real data is packaged into
  subroutine convert_real_to_complex(real_data, complex_data)
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), intent(in) :: real_data
    complex(C_DOUBLE_COMPLEX), dimension(:,:,:), contiguous, pointer, intent(out) :: complex_data
    
    integer :: i, j, k

    complex_data=cmplx(0.0d0, 0.0d0, kind=C_DOUBLE_COMPLEX)

    do i=1,size(real_data,3)
      do j=1,size(real_data,2)
        do k=1,size(real_data,1),2
          complex_data((k+1)/2,j,i)=cmplx(real_data(k,j,i), real_data(k+1,j,i), kind=C_DOUBLE_COMPLEX)
        end do
      end do
    end do
  end subroutine convert_real_to_complex

  !> Determines my global start coordinate in Fourier space. 
  !! This is required for cos y and cos x calculation which is fed into the tridiagonal solver. After the forward FFTs, 
  !! each process has ((n/2+1)/p+r) * 2 elements, where p is the number of processes and r
  !! is the uneven process remainder (1 or 0 depending on p). Therefore some processes will have t elements, and some t-2 elements
  !! to feed into the solver
  !! @param current_state The current model state
  !! @param dimension The dimension that we are calculating this for (Y or X)
  !! @returns My global start in Fourier space
  integer function deduce_my_global_start(current_state, dimension)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dimension

    integer complex_size, distributed_size, remainder, larger_nums, smaller_nums

    complex_size=(current_state%global_grid%size(dimension)/2+1)*2
    distributed_size=complex_size / current_state%parallel%dim_sizes(dimension)
    remainder=complex_size - distributed_size * current_state%parallel%dim_sizes(dimension) 
    larger_nums=min(remainder, current_state%parallel%my_coords(dimension))
    smaller_nums=current_state%parallel%my_coords(dimension)-remainder
    deduce_my_global_start=((distributed_size+1)*larger_nums + merge(distributed_size*smaller_nums, 0, smaller_nums .gt. 0)) + 1
  end function deduce_my_global_start
end module pencil_fft_mod    
    
