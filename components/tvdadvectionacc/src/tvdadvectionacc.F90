!> Implements TVD advection for prognostic fields
module tvdadvectionacc_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use grids_mod, only : local_grid_type, global_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_ERROR, LOG_INFO, log_master_log
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_information_type, component_field_value_type
  use mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE
  use optionsdatabase_mod, only : options_get_string
  use prognostics_mod, only : prognostic_field_type, prognostic_field_ptr_type
  use state_mod, only : model_state_type, parallel_state_type, FORWARD_STEPPING
  use stencil_mod, only : grid_stencil_type, interpolate_to_dual, create_stencil, free_stencil
  use ultimatefluxacc_mod, only : ultflx
  implicit none

  private
  logical :: advect_flow, advect_th, advect_q
  logical :: th_flux_pending = .false.

  type(grid_stencil_type), save :: star_stencil
  integer, save :: u_index=0, v_index=0, w_index=0
  type(prognostic_field_type), dimension(:), allocatable :: interpolated_fields

  integer, dimension(3) :: local_domain_start_index, local_domain_end_index, local_grid_size, &
      parallel_my_coords, parallel_dim_sizes
  integer :: kmax
  real(kind=DEFAULT_PRECISION) :: cx, cy

  ! Dynamically allocated arrays for persistent device residency
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: rdzn, dzn, rdz
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: u_field, v_field, w_field, d_flux_x, d_flux_y, d_flux_z
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: flux_y_buffer

  ! Dedicated buffer for accelerated theta tendency so current_state%sth is not double-summed
  type(prognostic_field_type), target :: sth_acc_field

  public tvdadvectionacc_get_descriptor
contains
  type(component_descriptor_type) function tvdadvectionacc_get_descriptor()
    tvdadvectionacc_get_descriptor%name="tvd_advection_acc"
    tvdadvectionacc_get_descriptor%version=0.1
    tvdadvectionacc_get_descriptor%initialisation=>initialisation_callback
    tvdadvectionacc_get_descriptor%timestep=>timestep_callback
    tvdadvectionacc_get_descriptor%finalisation=>finalisation_callback
    tvdadvectionacc_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    tvdadvectionacc_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(tvdadvectionacc_get_descriptor%published_fields(1))
    tvdadvectionacc_get_descriptor%published_fields(1)="tvd_advection_acc_sth"
  end function tvdadvectionacc_get_descriptor

  subroutine field_information_retrieval_callback(current_state, name, field_sizing)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_sizing

    field_sizing%field_type = COMPONENT_ARRAY_FIELD_TYPE
    field_sizing%data_type  = COMPONENT_DOUBLE_DATA_TYPE

    if (name .eq. "tvd_advection_acc_sth") then
      field_sizing%enabled = advect_th
      if (advect_th) then
        !$acc wait(11)
        field_sizing%number_dimensions = 3
        field_sizing%dimension_sizes   = [size(sth_acc_field%data, 1), &
                                          size(sth_acc_field%data, 2), &
                                          size(sth_acc_field%data, 3), 0]
      else
        field_sizing%number_dimensions = 0
      end if
    else
      call log_master_log(LOG_ERROR, "Mismatch between published_fields and field_information_retrieval_callback in tvdadvectionacc_mod.")
    end if
  end subroutine field_information_retrieval_callback

  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    integer :: send_handle, recv_handle, x, y

    if (name .eq. "tvd_advection_acc_sth") then
      if (advect_th .and. th_flux_pending) then
        !$acc wait(11)
        send_handle = register_y_flux_wrap_send_if_required( &
            d_flux_y(:, current_state%local_grid%local_domain_start_index(Y_INDEX)-1, :), &
            flux_y_buffer, current_state%parallel, current_state%local_grid)
        call complete_y_flux_wrap_send_if_required(current_state%parallel, current_state%local_grid, send_handle)

        recv_handle = register_y_flux_wrap_recv_if_required(flux_y_buffer, current_state%parallel, current_state%local_grid)
        call complete_y_flux_wrap_recv_if_required( &
            d_flux_y(:, current_state%local_grid%local_domain_end_index(Y_INDEX), :), flux_y_buffer, &
            current_state%parallel, current_state%local_grid, recv_handle)

        do x = current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do y = current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            call differentiate_face_values(y, x, current_state%u, current_state%v, current_state%w, &
                                           y, x, sth_acc_field, current_state%local_grid, current_state%global_grid, &
                                           d_flux_y(:, y-1, x), d_flux_x(:, y, x+1), &
                                           current_state%global_grid%configuration%vertical%tzc1, &
                                           current_state%global_grid%configuration%vertical%tzc2, &
                                           .true., d_flux_x(:, y, x), d_flux_y(:, y, x), d_flux_z(:, y, x))
          end do
        end do
        th_flux_pending = .false.
      end if

      field_value%real_3d_array_ptr => sth_acc_field%data
    else
      call log_master_log(LOG_ERROR, "Mismatch between published_fields and field_value_retrieval_callback in tvdadvectionacc_mod.")
    end if
  end subroutine field_value_retrieval_callback

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(prognostic_field_ptr_type), dimension(3) :: fields
    integer, dimension(3, 2) :: sizes
    integer :: num_fields
    logical :: xdim, ydim
    integer :: sz, sy, sx

    xdim=.false.
    ydim=.false.
    num_fields=0
#ifdef U_ACTIVE
    xdim=.true.
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%u
    sizes(num_fields,:) = (/ 2, 2 /)
    u_index = num_fields
#endif

#ifdef V_ACTIVE
    ydim=.true.
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%v
    sizes(num_fields,:) = (/ 1, 1 /)
    v_index=num_fields
#endif

#ifdef W_ACTIVE
    num_fields = num_fields + 1
    fields(num_fields)%ptr => current_state%w
    sizes(num_fields,:) = (/ 1, 1 /)
    w_index=num_fields
#endif
    allocate(interpolated_fields(0:num_fields))
#ifdef U_ACTIVE
    allocate(interpolated_fields(u_index)%data(current_state%global_grid%size(Z_INDEX), -1:3, -1:3))
    interpolated_fields(u_index)%active=.true.
#endif
#ifdef V_ACTIVE
    allocate(interpolated_fields(v_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(v_index)%active=.true.
#endif
#ifdef W_ACTIVE
    allocate(interpolated_fields(w_index)%data(current_state%global_grid%size(Z_INDEX), 0:2, 0:2))
    interpolated_fields(w_index)%active=.true.
#endif

    star_stencil = create_stencil(current_state%local_grid, fields, num_fields, 3, sizes, xdim, ydim)

    advect_flow = determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))
    advect_th   = determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q    = determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))

    advect_th = advect_th .and. current_state%th%active
    advect_q  = advect_q .and. (current_state%number_q_fields .gt. 0)

    if (advect_flow) call log_master_log(LOG_INFO, "Accelerating TVD advection: flow fields")
    if (advect_th)   call log_master_log(LOG_INFO, "Accelerating TVD advection: theta field")
    if (advect_q)    call log_master_log(LOG_INFO, "Accelerating TVD advection: q fields")

    if (advect_flow .or. advect_th .or. advect_q) then
      local_domain_start_index = current_state%local_grid%local_domain_start_index
      local_domain_end_index   = current_state%local_grid%local_domain_end_index
      local_grid_size          = current_state%local_grid%size
      parallel_my_coords       = current_state%parallel%my_coords
      parallel_dim_sizes       = current_state%parallel%dim_sizes
      cx   = current_state%global_grid%configuration%horizontal%cx
      cy   = current_state%global_grid%configuration%horizontal%cy
      kmax = current_state%local_grid%size(Z_INDEX)

      allocate(rdz,  source=current_state%global_grid%configuration%vertical%rdz)
      allocate(rdzn, source=current_state%global_grid%configuration%vertical%rdzn)
      allocate(dzn,  source=current_state%global_grid%configuration%vertical%dzn)

      ! Full halo-padded dimensions
      sz = size(current_state%u%data, 1)
      sy = size(current_state%u%data, 2)
      sx = size(current_state%u%data, 3)

      allocate(u_field(sz, sy, sx))
      allocate(v_field(sz, sy, sx))
      allocate(w_field(sz, sy, sx))
      allocate(d_flux_x(sz, sy, sx))
      allocate(d_flux_y(sz, sy, sx))
      allocate(d_flux_z(sz, sy, sx))
      allocate(flux_y_buffer(sz, sx))

      d_flux_x = 0.0_DEFAULT_PRECISION
      d_flux_y = 0.0_DEFAULT_PRECISION
      d_flux_z = 0.0_DEFAULT_PRECISION

      !$acc enter data create(u_field, v_field, w_field) &
      !$acc copyin(local_domain_start_index, local_domain_end_index, local_grid_size, &
      !$acc        parallel_my_coords, parallel_dim_sizes, cx, cy, kmax, &
      !$acc        rdz, rdzn, dzn, d_flux_x, d_flux_y, d_flux_z)
    end if

    if (advect_th) then
      allocate(sth_acc_field%data(size(current_state%sth%data, 1), &
                                  size(current_state%sth%data, 2), &
                                  size(current_state%sth%data, 3)))
      sth_acc_field%data = 0.0_DEFAULT_PRECISION
      sth_acc_field%active = .true.
    end if
  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    call free_stencil(star_stencil)
    if (allocated(interpolated_fields)) deallocate(interpolated_fields)

    if (advect_th) then
      if (allocated(sth_acc_field%data)) deallocate(sth_acc_field%data)
    end if

    if (advect_flow .or. advect_th .or. advect_q) then
      !$acc exit data delete(u_field, v_field, w_field, d_flux_x, d_flux_y, d_flux_z, &
      !$acc                  rdz, rdzn, dzn, local_domain_start_index, local_domain_end_index, &
      !$acc                  local_grid_size, parallel_my_coords, parallel_dim_sizes, cx, cy, kmax)
      if (allocated(u_field))       deallocate(u_field)
      if (allocated(v_field))       deallocate(v_field)
      if (allocated(w_field))       deallocate(w_field)
      if (allocated(d_flux_x))      deallocate(d_flux_x)
      if (allocated(d_flux_y))      deallocate(d_flux_y)
      if (allocated(d_flux_z))      deallocate(d_flux_z)
      if (allocated(flux_y_buffer)) deallocate(flux_y_buffer)
      if (allocated(rdz))           deallocate(rdz)
      if (allocated(rdzn))          deallocate(rdzn)
      if (allocated(dzn))           deallocate(dzn)
    end if
  end subroutine finalisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (advect_flow .or. advect_th .or. advect_q) then
      u_field = current_state%u%data
      v_field = current_state%v%data
      w_field = current_state%w%data
      !$acc update device(u_field, v_field, w_field) async(10)
    end if

    if (advect_flow) call advect_flow_fields(current_state)
    if (advect_th)   call advect_theta(current_state)
    if (advect_q)    call advect_q_fields(current_state)
  end subroutine timestep_callback

  subroutine advect_flow_fields(current_state)
    type(model_state_type), intent(inout) :: current_state

#ifdef U_ACTIVE
    call advect_u(current_state)
#endif
#ifdef V_ACTIVE
    call advect_v(current_state)
#endif
#ifdef W_ACTIVE
    call advect_w(current_state)
#endif
  end subroutine advect_flow_fields

  subroutine advect_q_fields(current_state)
    type(model_state_type), intent(inout) :: current_state
    integer :: i

    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        if (current_state%field_stepping .eq. FORWARD_STEPPING) then
          call advect_scalar_field(current_state, current_state%q(i), current_state%sq(i), .false.)
        else
          call advect_scalar_field(current_state, current_state%zq(i), current_state%sq(i), .false.)
        end if
      end if
    end do
  end subroutine advect_q_fields

  subroutine advect_theta(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (current_state%th%active) then
      current_state%sth%data = 0.0_DEFAULT_PRECISION
      if (current_state%field_stepping .eq. FORWARD_STEPPING) then
        call advect_scalar_field(current_state, current_state%th, sth_acc_field, .true.)
      else
        call advect_scalar_field(current_state, current_state%zth, sth_acc_field, .true.)
      end if
    end if
  end subroutine advect_theta

  ! Explicit-shape dummy arguments (sz, sy, sx) prevent nvfortran fort2 from segfaulting
  subroutine ultflx_scalar_accelerated(sx, sy, sz, local_domain_start_idx, u, v, w, zf, &
       local_domain_end_idx, local_grid_sz, cx_val, cy_val, parallel_coords, parallel_dims, &
       kdof, dt, flux_y, flux_z, flux_x, rdz_arr, rdzn_arr, dzn_arr, kmin, kmax_val)

    integer, intent(in) :: sx, sy, sz, kdof, kmin, kmax_val
    integer, dimension(3), intent(in) :: local_domain_start_idx, local_domain_end_idx, local_grid_sz
    integer, dimension(3), intent(in) :: parallel_coords, parallel_dims
    real(kind=DEFAULT_PRECISION), intent(in) :: dt, cx_val, cy_val
    real(kind=DEFAULT_PRECISION), dimension(sz,sy,sx), intent(in) :: u, v, w, zf
    real(kind=DEFAULT_PRECISION), dimension(sz), intent(in) :: rdz_arr, rdzn_arr, dzn_arr
    real(kind=DEFAULT_PRECISION), dimension(sz,sy,sx), intent(inout) :: flux_y, flux_z, flux_x

    integer :: x, y

    !$acc parallel loop collapse(2) wait(10) async(10) default(none) &
    !$acc copyin(zf, dt, kdof, kmin) &
    !$acc present(u, v, w, flux_x, flux_y, flux_z, rdz_arr, rdzn_arr, dzn_arr, &
    !$acc         local_domain_start_idx, local_domain_end_idx, local_grid_sz, &
    !$acc         parallel_coords, parallel_dims, cx_val, cy_val, kmax_val)
    do x = local_domain_start_idx(X_INDEX)-1, local_domain_end_idx(X_INDEX)
      do y = local_domain_start_idx(Y_INDEX)-1, local_domain_end_idx(Y_INDEX)
        call ultflx(y, x, u, v, w, &
                    y, x, zf, local_domain_end_idx, local_grid_sz, &
                    cx_val, cy_val, &
                    parallel_coords, parallel_dims, kdof, dt, &
                    flux_y(:, y, x), flux_z(:, y, x), flux_x(:, y, x+1), &
                    rdz_arr, rdzn_arr, dzn_arr, kmin, kmax_val)
      end do
    end do
    !$acc update host(flux_x, flux_y, flux_z) wait(10) async(11)
  end subroutine ultflx_scalar_accelerated

  subroutine advect_scalar_field(current_state, field, source_field, defer_differentiation)
    type(model_state_type), intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: field
    type(prognostic_field_type), intent(inout) :: source_field
    logical, intent(in) :: defer_differentiation

    real(kind=DEFAULT_PRECISION) :: dt
    integer :: x, y, recv_handle, send_handle

    dt = current_state%dtm
    if (current_state%scalar_stepping .ne. FORWARD_STEPPING) then
      dt = dt * 2.0_DEFAULT_PRECISION
    end if

    ! Pass full halo-padded sizes (size(u_field, ...)) and distinct kmax scalar
    call ultflx_scalar_accelerated(size(u_field, 3), size(u_field, 2), size(u_field, 1), &
                                   local_domain_start_index, u_field, v_field, w_field, field%data, &
                                   local_domain_end_index, local_grid_size, cx, cy, &
                                   parallel_my_coords, parallel_dim_sizes, &
                                   0, dt, d_flux_y, d_flux_z, d_flux_x, &
                                   rdz, rdzn, dzn, 2, kmax)

    if (defer_differentiation) then
      th_flux_pending = .true.
      return
    end if

    !$acc wait(11)
    send_handle = register_y_flux_wrap_send_if_required( &
        d_flux_y(:, current_state%local_grid%local_domain_start_index(Y_INDEX)-1, :), &
        flux_y_buffer, current_state%parallel, current_state%local_grid)
    call complete_y_flux_wrap_send_if_required(current_state%parallel, current_state%local_grid, send_handle)

    recv_handle = register_y_flux_wrap_recv_if_required(flux_y_buffer, current_state%parallel, current_state%local_grid)
    call complete_y_flux_wrap_recv_if_required( &
        d_flux_y(:, current_state%local_grid%local_domain_end_index(Y_INDEX), :), &
        flux_y_buffer, current_state%parallel, current_state%local_grid, recv_handle)

    do x = current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do y = current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        call differentiate_face_values(y, x, current_state%u, current_state%v, current_state%w, &
                                       y, x, source_field, current_state%local_grid, current_state%global_grid, &
                                       d_flux_y(:, y-1, x), d_flux_x(:, y, x+1), &
                                       current_state%global_grid%configuration%vertical%tzc1, &
                                       current_state%global_grid%configuration%vertical%tzc2, &
                                       .true., d_flux_x(:, y, x), d_flux_y(:, y, x), d_flux_z(:, y, x))
      end do
    end do
  end subroutine advect_scalar_field

  subroutine advect_flow_field(current_state, field_index, interpolate_field, advect_field, source_field, kdof, rdz1, rdz2, dz, tz1, tz2, kmin, kmax_arg)
    type(model_state_type), intent(inout) :: current_state
    type(prognostic_field_type), intent(inout) :: interpolate_field
    type(prognostic_field_type), intent(inout) :: advect_field
    type(prognostic_field_type), intent(inout) :: source_field
    integer, intent(in) :: field_index, kdof, kmin, kmax_arg
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: rdz1, rdz2, dz, tz1, tz2

    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: local_flux_y_buffer
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: flux_x, flux_y, flux_z
    real(kind=DEFAULT_PRECISION) :: dt
    integer :: x, y, recv_handle, send_handle

    allocate(flux_x(current_state%local_grid%size(Z_INDEX), &
                    current_state%local_grid%size(Y_INDEX)+2*current_state%local_grid%halo_size(Y_INDEX), &
                    current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)))
    allocate(flux_y(current_state%local_grid%size(Z_INDEX), &
                    current_state%local_grid%size(Y_INDEX)+2*current_state%local_grid%halo_size(Y_INDEX), &
                    current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)))
    allocate(flux_z(current_state%local_grid%size(Z_INDEX), &
                    current_state%local_grid%size(Y_INDEX)+2*current_state%local_grid%halo_size(Y_INDEX), &
                    current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)))

    allocate(local_flux_y_buffer(current_state%local_grid%size(Z_INDEX), &
                                 current_state%local_grid%size(X_INDEX)+2*current_state%local_grid%halo_size(X_INDEX)))

    dt = current_state%dtm
    if (current_state%momentum_stepping .ne. FORWARD_STEPPING) then
      dt = dt * 2.0_DEFAULT_PRECISION
    end if

    do x = current_state%local_grid%local_domain_start_index(X_INDEX)-1, current_state%local_grid%local_domain_end_index(X_INDEX)
      do y = current_state%local_grid%local_domain_start_index(Y_INDEX)-1, current_state%local_grid%local_domain_end_index(Y_INDEX)
        call interpolate_to_dual(current_state%local_grid, interpolate_field, star_stencil, x, y, interpolated_fields, field_index)
        call ultflx(1, 1, interpolated_fields(u_index)%data, interpolated_fields(v_index)%data, interpolated_fields(w_index)%data, &
                    y, x, advect_field%data, current_state%local_grid%local_domain_end_index, current_state%local_grid%size, &
                    current_state%global_grid%configuration%horizontal%cx, current_state%global_grid%configuration%horizontal%cy, &
                    current_state%parallel%my_coords, current_state%parallel%dim_sizes, &
                    kdof, dt, flux_y(:, y, x), flux_z(:, y, x), flux_x(:, y, x+1), &
                    rdz1, rdz2, dz, kmin, kmax_arg)
      end do
    end do

    send_handle = register_y_flux_wrap_send_if_required(flux_y(:, current_state%local_grid%local_domain_start_index(Y_INDEX)-1, :), local_flux_y_buffer, current_state%parallel, current_state%local_grid)
    call complete_y_flux_wrap_send_if_required(current_state%parallel, current_state%local_grid, send_handle)

    recv_handle = register_y_flux_wrap_recv_if_required(local_flux_y_buffer, current_state%parallel, current_state%local_grid)
    call complete_y_flux_wrap_recv_if_required(flux_y(:, current_state%local_grid%local_domain_end_index(Y_INDEX), :), local_flux_y_buffer, current_state%parallel, current_state%local_grid, recv_handle)

    do x = current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do y = current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        call differentiate_face_values(1, 1, interpolated_fields(u_index), interpolated_fields(v_index), interpolated_fields(w_index), &
                                       y, x, source_field, current_state%local_grid, current_state%global_grid, &
                                       flux_y(:, y-1, x), flux_x(:, y, x+1), &
                                       tz1, tz2, .true., flux_x(:, y, x), flux_y(:, y, x), flux_z(:, y, x))
      end do
    end do

    deallocate(flux_x)
    deallocate(flux_y)
    deallocate(flux_z)
    deallocate(local_flux_y_buffer)
  end subroutine advect_flow_field

  subroutine advect_u(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (current_state%field_stepping .eq. FORWARD_STEPPING) then
      call advect_flow_field(current_state, u_index, current_state%u, current_state%u, current_state%su, 0, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%dzn, &
                             current_state%global_grid%configuration%vertical%tzc1, &
                             current_state%global_grid%configuration%vertical%tzc2, &
                             2, current_state%local_grid%size(Z_INDEX))
    else
      call advect_flow_field(current_state, u_index, current_state%u, current_state%zu, current_state%su, 0, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%dzn, &
                             current_state%global_grid%configuration%vertical%tzc1, &
                             current_state%global_grid%configuration%vertical%tzc2, &
                             2, current_state%local_grid%size(Z_INDEX))
    end if
  end subroutine advect_u

  subroutine advect_v(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (current_state%field_stepping .eq. FORWARD_STEPPING) then
      call advect_flow_field(current_state, v_index, current_state%v, current_state%v, current_state%sv, 0, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%dzn, &
                             current_state%global_grid%configuration%vertical%tzc1, &
                             current_state%global_grid%configuration%vertical%tzc2, &
                             2, current_state%local_grid%size(Z_INDEX))
    else
      call advect_flow_field(current_state, v_index, current_state%v, current_state%zv, current_state%sv, 0, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%dzn, &
                             current_state%global_grid%configuration%vertical%tzc1, &
                             current_state%global_grid%configuration%vertical%tzc2, &
                             2, current_state%local_grid%size(Z_INDEX))
    end if
  end subroutine advect_v

  subroutine advect_w(current_state)
    type(model_state_type), intent(inout) :: current_state

    if (current_state%field_stepping .eq. FORWARD_STEPPING) then
      call advect_flow_field(current_state, w_index, current_state%w, current_state%w, current_state%sw, 1, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%dz, &
                             current_state%global_grid%configuration%vertical%tzd1, &
                             current_state%global_grid%configuration%vertical%tzd2, &
                             1, current_state%local_grid%size(Z_INDEX)-1)
    else
      call advect_flow_field(current_state, w_index, current_state%w, current_state%zw, current_state%sw, 1, &
                             current_state%global_grid%configuration%vertical%rdzn, &
                             current_state%global_grid%configuration%vertical%rdz, &
                             current_state%global_grid%configuration%vertical%dz, &
                             current_state%global_grid%configuration%vertical%tzd1, &
                             current_state%global_grid%configuration%vertical%tzd2, &
                             1, current_state%local_grid%size(Z_INDEX)-1)
    end if
  end subroutine advect_w

  subroutine differentiate_face_values(y_flow_index, x_flow_index, u, v, w, y_source_index, x_source_index, source_field, &
       local_grid, global_grid, flux_y_previous, flux_x_previous, tzc1, tzc2, differentiate_top, flux_x, flux_y, flux_z)

    integer, intent(in) :: y_flow_index, x_flow_index, y_source_index, x_source_index
    logical, intent(in) :: differentiate_top
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: tzc1, tzc2
    type(prognostic_field_type), intent(in) :: u, w, v
    type(prognostic_field_type), intent(inout) :: source_field
    type(global_grid_type), intent(in) :: global_grid
    type(local_grid_type), intent(in) :: local_grid
    real(kind=DEFAULT_PRECISION), intent(in), dimension(:) :: flux_y_previous, flux_x_previous, flux_x, flux_y, flux_z

    integer :: k

    do k=2,local_grid%size(Z_INDEX)-1
      source_field%data(k, y_source_index, x_source_index)=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*(w%data(k-1, y_flow_index, x_flow_index)* flux_z(k)*tzc1(k) - &
           w%data(k, y_flow_index, x_flow_index)*flux_z(k+1)*tzc2(k))
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) - u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end do
    if (differentiate_top) then
      k=local_grid%size(Z_INDEX)
      source_field%data(k, y_source_index, x_source_index)=0.0_DEFAULT_PRECISION
#ifdef V_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (v%data(k, y_flow_index-1, x_flow_index)* flux_y_previous(k) - v%data(k, y_flow_index, x_flow_index)*flux_y(k))*&
           global_grid%configuration%horizontal%cy
#endif
#ifdef W_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           4.0_DEFAULT_PRECISION*tzc1(k)* w%data(k-1, y_flow_index, x_flow_index)*flux_z(k)
#endif
#ifdef U_ACTIVE
      source_field%data(k, y_source_index, x_source_index)=source_field%data(k, y_source_index, x_source_index)+&
           (u%data(k, y_flow_index, x_flow_index-1)* flux_x(k) -u%data(k, y_flow_index, x_flow_index)*flux_x_previous(k))*&
           global_grid%configuration%horizontal%cx
#endif
    end if
  end subroutine differentiate_face_values

  subroutine complete_y_flux_wrap_send_if_required(parallel, local_grid, handle)
    type(parallel_state_type), intent(in) :: parallel
    type(local_grid_type), intent(in) :: local_grid
    integer, intent(in) :: handle
    integer :: ierr

    if (parallel%my_coords(Y_INDEX) == 0 .and. handle .ne. MPI_REQUEST_NULL) then
      call mpi_wait(handle, MPI_STATUS_IGNORE, ierr)
    end if
  end subroutine complete_y_flux_wrap_send_if_required

  integer function register_y_flux_wrap_send_if_required(flux_y, flux_y_buffer, parallel, local_grid)
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: flux_y
    real(kind=DEFAULT_PRECISION), dimension(:,:), contiguous, intent(out) :: flux_y_buffer
    type(parallel_state_type), intent(in) :: parallel
    type(local_grid_type), intent(in) :: local_grid
    integer :: ierr

    register_y_flux_wrap_send_if_required = MPI_REQUEST_NULL
    if (parallel%my_coords(Y_INDEX) == 0) then
      flux_y_buffer = flux_y
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,1)) then
        call mpi_isend(flux_y_buffer, size(flux_y_buffer), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,1), 0, &
             parallel%neighbour_comm, register_y_flux_wrap_send_if_required, ierr)
      end if
    end if
  end function register_y_flux_wrap_send_if_required

  subroutine complete_y_flux_wrap_recv_if_required(flux_y, flux_y_buffer, parallel, local_grid, handle)
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(out) :: flux_y
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: flux_y_buffer
    type(parallel_state_type), intent(in) :: parallel
    type(local_grid_type), intent(in) :: local_grid
    integer, intent(in) :: handle
    integer :: ierr

    if (parallel%my_coords(Y_INDEX) == parallel%dim_sizes(Y_INDEX)-1) then
      if (handle .ne. MPI_REQUEST_NULL) then
        call mpi_wait(handle, MPI_STATUS_IGNORE, ierr)
      end if
      flux_y = flux_y_buffer
    end if
  end subroutine complete_y_flux_wrap_recv_if_required

  integer function register_y_flux_wrap_recv_if_required(flux_y_buffer, parallel, local_grid)
    real(kind=DEFAULT_PRECISION), dimension(:,:), contiguous, intent(out) :: flux_y_buffer
    type(parallel_state_type), intent(in) :: parallel
    type(local_grid_type), intent(in) :: local_grid
    integer :: ierr

    register_y_flux_wrap_recv_if_required = MPI_REQUEST_NULL
    if (parallel%my_coords(Y_INDEX) == parallel%dim_sizes(Y_INDEX)-1) then
      if (parallel%my_rank .ne. local_grid%neighbours(Y_INDEX,3)) then
        call mpi_irecv(flux_y_buffer, size(flux_y_buffer), PRECISION_TYPE, local_grid%neighbours(Y_INDEX,3), 0, &
             parallel%neighbour_comm, register_y_flux_wrap_recv_if_required, ierr)
      end if
    end if
  end function register_y_flux_wrap_recv_if_required

  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "tvdacc" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here
end module tvdadvectionacc_mod
