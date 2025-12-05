!> Fourier space spectral filter. Adapted from fftsolver. It uses the pencil FFT
!! module for 3D FFTs in pencil decomposition. These use FFTW for the actual FFT kernel.
module filter_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use state_mod, only : model_state_type
  use monc_component_mod, only : component_descriptor_type
  use pencil_fft_mod, only : initialise_pencil_fft, finalise_pencil_fft, perform_forward_3dfft, perform_backwards_3dfft
  use communication_types_mod, only : halo_communication_type, neighbour_description_type, field_data_wrapper_type
  use halo_communication_mod, only : copy_buffer_to_field, copy_field_to_buffer, perform_local_data_copy_for_field, &
       init_halo_communication, finalise_halo_communication, blocking_halo_swap, get_single_field_per_halo_cell
  use registry_mod, only : is_component_enabled
  use logging_mod, only : LOG_ERROR, log_master_log
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real
  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: u_filter_enabled     ! activate filter
  logical :: v_filter_enabled     ! activate filter
  logical :: w_filter_enabled     ! activate filter
  logical :: th_filter_enabled     ! activate filter
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: tmp_filt
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: uf
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: vf
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: wf
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: thf
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: filter_coefs
  integer :: fourier_space_sizes(3), kmax
  type(halo_communication_type), save :: halo_swap_state

  public filter_get_descriptor
contains

  !> Descriptor of this component for registration
  !! @returns The filter component descriptor
  type(component_descriptor_type) function filter_get_descriptor()
    filter_get_descriptor%name="filter"
    filter_get_descriptor%version=0.1
    filter_get_descriptor%initialisation=>initialisation_callback
    filter_get_descriptor%timestep=>timestep_callback
    filter_get_descriptor%finalisation=>finalisation_callback
  end function filter_get_descriptor

  !> This initialisation callback sets up the pencil fft module, allocates data for the fourier space pressure term (might be
  !! different size than p) and populates cos x and y
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: i, j, my_y_start, my_x_start, x_size, y_size
    real(kind=DEFAULT_PRECISION) :: tmp1, tmp2, tmp3, fn, ff, cutoff

    u_filter_enabled=options_get_logical(current_state%options_database, "u_filter_enabled")
    v_filter_enabled=options_get_logical(current_state%options_database, "v_filter_enabled")
    w_filter_enabled=options_get_logical(current_state%options_database, "w_filter_enabled")
    th_filter_enabled=options_get_logical(current_state%options_database, "th_filter_enabled")

    if (u_filter_enabled .or. v_filter_enabled .or. w_filter_enabled .or. th_filter_enabled) then
      fourier_space_sizes=initialise_pencil_fft(current_state, my_y_start, my_x_start)
      call init_halo_communication(current_state, get_single_field_per_halo_cell, halo_swap_state, 1, .false.)
      y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
      x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
      allocate(tmp_filt(current_state%local_grid%size(Z_INDEX), y_size, x_size))
      allocate(filter_coefs(fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
      kmax=options_get_integer(current_state%options_database, "f_kmax")
      cutoff=options_get_real(current_state%options_database, "f_cutoff")
      fn=options_get_real(current_state%options_database, "f_order")
      if (kmax .eq. -1) kmax = fourier_space_sizes(Z_INDEX)
      ff = real(current_state%global_grid%size(X_INDEX), kind=DEFAULT_PRECISION)/cutoff
      do i=1, fourier_space_sizes(X_INDEX)
        do j=1, fourier_space_sizes(Y_INDEX)
          tmp1=real(((i-1)+(my_x_start-1)), kind=DEFAULT_PRECISION)
          tmp2=real(((j-1)+(my_y_start-1)), kind=DEFAULT_PRECISION)
          tmp3 = sqrt(tmp1*tmp1 + tmp2*tmp2)
          filter_coefs(j,i) = 1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION + (tmp3 /ff)**(2.0_DEFAULT_PRECISION*fn))
        end do
      end do
    end if


    if (u_filter_enabled) then
      allocate(uf(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    end if ! u_filter_enabled

    if (v_filter_enabled) then
      allocate(vf(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    end if ! v_filter_enabled

    if (w_filter_enabled) then
      allocate(wf(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    end if ! w_filter_enabled

    if (th_filter_enabled) then
      allocate(thf(fourier_space_sizes(Z_INDEX), fourier_space_sizes(Y_INDEX), fourier_space_sizes(X_INDEX)))
    end if ! th_filter_enabled

  end subroutine initialisation_callback  


  !> Timestep call back 
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer ipi, ipj, ipk
    character (len=50):: fname

    integer :: start_loc(3), end_loc(3), i
    if (u_filter_enabled .or. v_filter_enabled .or. w_filter_enabled .or. th_filter_enabled) then

        do i=1,3
          start_loc(i)=current_state%local_grid%local_domain_start_index(i)
          end_loc(i)=current_state%local_grid%local_domain_end_index(i)
        end do
      
      if (u_filter_enabled) then
        call perform_forward_3dfft(current_state, current_state%u%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
             start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), uf)
        do ipi=1,fourier_space_sizes(X_INDEX)
          do ipj=1,fourier_space_sizes(Y_INDEX)
            do ipk=1,kmax
            uf(ipk,ipj,ipi) = uf(ipk,ipj,ipi)*filter_coefs(ipj, ipi)
            end do
          end do
        end do

        if (current_state%immersed%ib_enabled) then ! apply only in proximity of IB
          call perform_backwards_3dfft(current_state, uf, tmp_filt(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
          
          do ipi=start_loc(X_INDEX),end_loc(X_INDEX)
            do ipj=start_loc(Y_INDEX),end_loc(Y_INDEX)
              do ipk=1,kmax
                current_state%u%data(ipk,ipj,ipi) = current_state%u%data(ipk,ipj,ipi)*current_state%immersed%ibprox(ipk,ipj,ipi) +&
                                        tmp_filt(ipk,ipj,ipi)*(1.0_DEFAULT_PRECISION-current_state%immersed%ibprox(ipk,ipj,ipi))
              end do
            end do
          end do
        else
          call perform_backwards_3dfft(current_state, uf, current_state%u%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
        end if
        call blocking_halo_swap(current_state, halo_swap_state, copy_u_to_halo_buffer, &
             perform_local_data_copy_for_u, copy_halo_buffer_to_u)
      end if ! filter_enabled



      if (v_filter_enabled) then
        call perform_forward_3dfft(current_state, current_state%v%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
             start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), vf)
        do ipi=1,fourier_space_sizes(X_INDEX)
          do ipj=1,fourier_space_sizes(Y_INDEX)
            do ipk=1,kmax
            vf(ipk,ipj,ipi) = vf(ipk,ipj,ipi)*filter_coefs(ipj, ipi)
            end do
          end do
        end do

        if (current_state%immersed%ib_enabled) then ! apply only in proximity of IB
          call perform_backwards_3dfft(current_state, vf, tmp_filt(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
          
          do ipi=start_loc(X_INDEX),end_loc(X_INDEX)
            do ipj=start_loc(Y_INDEX),end_loc(Y_INDEX)
              do ipk=1,kmax
                current_state%v%data(ipk,ipj,ipi) = current_state%v%data(ipk,ipj,ipi)*current_state%immersed%ibprox(ipk,ipj,ipi) +&
                                        tmp_filt(ipk,ipj,ipi)*(1.0_DEFAULT_PRECISION-current_state%immersed%ibprox(ipk,ipj,ipi))
              end do
            end do
          end do
        else
          call perform_backwards_3dfft(current_state, vf, current_state%v%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
        end if
    
        call blocking_halo_swap(current_state, halo_swap_state, copy_v_to_halo_buffer, &
             perform_local_data_copy_for_v, copy_halo_buffer_to_v)
      end if ! filter_enabled
           
      if (w_filter_enabled) then
        call perform_forward_3dfft(current_state, current_state%w%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
             start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), wf)
        do ipi=1,fourier_space_sizes(X_INDEX)
          do ipj=1,fourier_space_sizes(Y_INDEX)
            do ipk=1,kmax
            wf(ipk,ipj,ipi) = wf(ipk,ipj,ipi)*filter_coefs(ipj, ipi)
            end do
          end do
        end do
    
        if (current_state%immersed%ib_enabled) then ! apply only in proximity of IB
          call perform_backwards_3dfft(current_state, wf, tmp_filt(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
          
          do ipi=start_loc(X_INDEX),end_loc(X_INDEX)
            do ipj=start_loc(Y_INDEX),end_loc(Y_INDEX)
              do ipk=1,kmax
                current_state%w%data(ipk,ipj,ipi) = current_state%w%data(ipk,ipj,ipi)*current_state%immersed%ibprox(ipk,ipj,ipi) +&
                                        tmp_filt(ipk,ipj,ipi)*(1.0_DEFAULT_PRECISION-current_state%immersed%ibprox(ipk,ipj,ipi))
              end do
            end do
          end do
        else
          call perform_backwards_3dfft(current_state, wf, current_state%w%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
        end if
    
        call blocking_halo_swap(current_state, halo_swap_state, copy_w_to_halo_buffer, &
             perform_local_data_copy_for_w, copy_halo_buffer_to_w)
      end if ! filter_enabled
           
      if (th_filter_enabled) then
        call perform_forward_3dfft(current_state, current_state%th%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
             start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)), thf)
        do ipi=1,fourier_space_sizes(X_INDEX)
          do ipj=1,fourier_space_sizes(Y_INDEX)
            do ipk=1,kmax
            thf(ipk,ipj,ipi) = thf(ipk,ipj,ipi)*filter_coefs(ipj, ipi)
            end do
          end do
        end do

        if (current_state%immersed%ib_enabled) then ! apply only in proximity of IB
          call perform_backwards_3dfft(current_state, thf, tmp_filt(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
          
          do ipi=start_loc(X_INDEX),end_loc(X_INDEX)
            do ipj=start_loc(Y_INDEX),end_loc(Y_INDEX)
              do ipk=1,kmax
                current_state%th%data(ipk,ipj,ipi) = &
                  current_state%th%data(ipk,ipj,ipi)*current_state%immersed%ibprox(ipk,ipj,ipi) +&
                  tmp_filt(ipk,ipj,ipi)*(1.0_DEFAULT_PRECISION-current_state%immersed%ibprox(ipk,ipj,ipi))
              end do
            end do
          end do
        else
          call perform_backwards_3dfft(current_state, thf, current_state%th%data(start_loc(Z_INDEX):end_loc(Z_INDEX), &
               start_loc(Y_INDEX):end_loc(Y_INDEX), start_loc(X_INDEX):end_loc(X_INDEX)))   
        end if
!    
    
        call blocking_halo_swap(current_state, halo_swap_state, copy_th_to_halo_buffer, &
             perform_local_data_copy_for_th, copy_halo_buffer_to_th)
      end if ! filter_enabled
           
           
    end if ! end timestep modulo
           
         

  end subroutine timestep_callback  

  !> Called at MONC finalisation, will call to the pencil fft module to clean itself up and free the pressure term
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

!     call finalise_pencil_fft(current_state%parallel%monc_communicator)
    if (allocated(uf)) deallocate(uf)
    if (allocated(vf)) deallocate(vf)
    if (allocated(wf)) deallocate(wf)
    if (allocated(thf)) deallocate(thf)
  end subroutine finalisation_callback

  
  
  !> Copies the u field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_u_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%u%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_u_to_halo_buffer

  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_u(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%u%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_u

  !> Does local data copying for u variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_u(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%u%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_u





  
  
  !> Copies the v field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_v_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%v%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_v_to_halo_buffer

  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_v(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%v%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_v

  !> Does local data copying for v variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_v(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%v%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_v





  
  
  !> Copies the w field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_w_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%w%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_w_to_halo_buffer

  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_w(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%w%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_w

  !> Does local data copying for w variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_w(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%w%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_w





  
  
  !> Copies the th field data to halo buffers for a specific process in a dimension and halo cell
  !! @param current_state The current model state
  !! @param neighbour_descriptions Description of the neighbour halo swapping status
  !! @param dim Dimension to copy from
  !! @param source_index The source index of the dimension we are reading from in the prognostic field
  !! @param pid_location Location of the neighbouring process in the local stored data structures
  !! @param current_page The current (next) buffer page to copy into
  !! @param source_data Optional source data which is read from
  subroutine copy_th_to_halo_buffer(current_state, neighbour_description, dim, source_index, &
       pid_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, pid_location, source_index
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_field_to_buffer(current_state%local_grid, neighbour_description%send_halo_buffer, current_state%th%data, &
         dim, source_index, current_page(pid_location))

    current_page(pid_location)=current_page(pid_location)+1
  end subroutine copy_th_to_halo_buffer

  !! @param dim The dimension we receive for
  !! @param target_index The target index for the dimension we are receiving for
  !! @param neighbour_location The location in the local neighbour data stores of this neighbour
  !! @param current_page The current, next, halo swap page to read from (all previous have been read and copied already)
  !! @param source_data Optional source data which is written into
  subroutine copy_halo_buffer_to_th(current_state, neighbour_description, dim, target_index, &
       neighbour_location, current_page, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: dim, target_index, neighbour_location
    integer, intent(inout) :: current_page(:)
    type(neighbour_description_type), intent(inout) :: neighbour_description
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call copy_buffer_to_field(current_state%local_grid, neighbour_description%recv_halo_buffer, current_state%th%data, &
         dim, target_index, current_page(neighbour_location))

    current_page(neighbour_location)=current_page(neighbour_location)+1
  end subroutine copy_halo_buffer_to_th

  !> Does local data copying for th variable halo swap
  !! @param current_state The current model state_mod
  !! @param source_data Optional source data which is written into
  subroutine perform_local_data_copy_for_th(current_state, halo_depth, involve_corners, source_data)
    type(model_state_type), intent(inout) :: current_state
    integer, intent(in) :: halo_depth
    logical, intent(in) :: involve_corners
    type(field_data_wrapper_type), dimension(:), intent(in), optional :: source_data

    call perform_local_data_copy_for_field(current_state%th%data, current_state%local_grid, &
         current_state%parallel%my_rank, halo_depth, involve_corners)
  end subroutine perform_local_data_copy_for_th






end module filter_mod
