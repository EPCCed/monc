!> Does the field stepping
!! Stepping is called at thedynamics-join stage over the entire local domain
module stepfields_mod
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_integer
  use state_mod, only : FORWARD_STEPPING, CENTRED_STEPPING, model_state_type
  use prognostics_mod, only : prognostic_field_type
  use grids_mod, only : local_grid_type, X_INDEX, Y_INDEX, Z_INDEX
  use datadefn_mod, only : DEFAULT_PRECISION
  use registry_mod, only : is_component_enabled, get_component_field_information, get_component_field_value
  use science_constants_mod, only : rlargep
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: determine_flow_minmax=.false., cfl_is_enabled

  real(kind=DEFAULT_PRECISION), allocatable :: resetq_min(:)
  logical :: l_nonconservative_positive_q=.true.

  ! Local tendency diagnostic variables for this component
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable ::     &
       tend_3d_th,tend_3d_qv,       &
       tend_3d_ql,tend_3d_qi,tend_3d_qr,tend_3d_qs,tend_3d_qg,       &
       tend_3d_tabs
  logical :: l_tend_3d_th,l_tend_3d_qv,       &
             l_tend_3d_ql,l_tend_3d_qi,l_tend_3d_qr,l_tend_3d_qs,l_tend_3d_qg,       &
             l_tend_3d_tabs
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable ::         &
       tend_pr_tot_th,tend_pr_tot_qv,       &
       tend_pr_tot_ql,tend_pr_tot_qi,tend_pr_tot_qr,tend_pr_tot_qs,tend_pr_tot_qg,       &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_th,l_tend_pr_tot_qv,       &
             l_tend_pr_tot_ql,l_tend_pr_tot_qi,l_tend_pr_tot_qr,l_tend_pr_tot_qs,l_tend_pr_tot_qg,       &
             l_tend_pr_tot_tabs
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  public stepfields_get_descriptor

contains

  type(component_descriptor_type) function stepfields_get_descriptor()
    stepfields_get_descriptor%name="stepfields"
    stepfields_get_descriptor%version=0.1
    stepfields_get_descriptor%initialisation=>initialisation_callback
    stepfields_get_descriptor%timestep=>timestep_callback
    stepfields_get_descriptor%finalisation=>finalisation_callback

    stepfields_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    stepfields_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(stepfields_get_descriptor%published_fields(8+8))

    stepfields_get_descriptor%published_fields(1)="tend_th_total_3d_local"
    stepfields_get_descriptor%published_fields(2)="tend_qv_total_3d_local"
    stepfields_get_descriptor%published_fields(3)="tend_ql_total_3d_local"
    stepfields_get_descriptor%published_fields(4)="tend_qi_total_3d_local"
    stepfields_get_descriptor%published_fields(5)="tend_qr_total_3d_local"
    stepfields_get_descriptor%published_fields(6)="tend_qs_total_3d_local"
    stepfields_get_descriptor%published_fields(7)="tend_qg_total_3d_local"
    stepfields_get_descriptor%published_fields(8)="tend_tabs_total_3d_local"

    stepfields_get_descriptor%published_fields(8+1)="tend_th_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+2)="tend_qv_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+3)="tend_ql_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+4)="tend_qi_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+5)="tend_qr_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+6)="tend_qs_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+7)="tend_qg_total_profile_total_local"
    stepfields_get_descriptor%published_fields(8+8)="tend_tabs_total_profile_total_local"
  end function stepfields_get_descriptor

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag

#ifdef U_ACTIVE
    allocate(current_state%global_grid%configuration%vertical%savolubar(current_state%local_grid%size(Z_INDEX)))
#endif
#ifdef V_ACTIVE
    allocate(current_state%global_grid%configuration%vertical%savolvbar(current_state%local_grid%size(Z_INDEX)))
#endif

    allocate(resetq_min(current_state%number_q_fields))
    cfl_is_enabled=is_component_enabled(current_state%options_database, "cfltest")
    if (cfl_is_enabled) call reset_local_minmax_values(current_state)

    l_qdiag = (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0)

    l_tend_pr_tot_th   = current_state%th%active
    l_tend_pr_tot_qv   = l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0
    l_tend_pr_tot_ql   = l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0
    l_tend_pr_tot_qi   = l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0
    l_tend_pr_tot_qr   = l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0
    l_tend_pr_tot_qs   = l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0
    l_tend_pr_tot_qg   = l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0
    l_tend_pr_tot_tabs = l_tend_pr_tot_th

    l_tend_3d_th   = current_state%th%active .or. l_tend_pr_tot_th
    l_tend_3d_qv   = (l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0) .or. l_tend_pr_tot_qv
    l_tend_3d_ql   = (l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_ql
    l_tend_3d_qi   = (l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qi
    l_tend_3d_qr   = (l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qr
    l_tend_3d_qs   = (l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qs
    l_tend_3d_qg   = (l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    if (l_tend_3d_th) then
      allocate(tend_3d_th(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qv) then
      iqv=get_q_index(standard_q_names%VAPOUR, 'stepfields')
      allocate(tend_3d_qv(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'stepfields')
      allocate(tend_3d_ql(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'stepfields')
      allocate(tend_3d_qi(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'stepfields')
      allocate(tend_3d_qr(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'stepfields')
      allocate(tend_3d_qs(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'stepfields')
      allocate(tend_3d_qg(current_state%local_grid%size(Z_INDEX), &
                          current_state%local_grid%size(Y_INDEX), &
                          current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_tabs) then
      allocate(tend_3d_tabs(current_state%local_grid%size(Z_INDEX), &
                            current_state%local_grid%size(Y_INDEX), &
                            current_state%local_grid%size(X_INDEX)))
    endif

    if (l_tend_pr_tot_th) allocate(tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qv) allocate(tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_ql) allocate(tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qi) allocate(tend_pr_tot_qi(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qr) allocate(tend_pr_tot_qr(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qs) allocate(tend_pr_tot_qs(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qg) allocate(tend_pr_tot_qg(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_tabs) allocate(tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)))
  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    deallocate(resetq_min)

    if (allocated(tend_3d_th)) deallocate(tend_3d_th)
    if (allocated(tend_3d_qv)) deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql)) deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi)) deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr)) deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs)) deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg)) deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_th)) deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv)) deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql)) deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi)) deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr)) deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs)) deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg)) deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)
  end subroutine finalisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    type(component_field_information_type) :: source_field_accelerated_info
    type(component_field_value_type) :: source_field_accelerated

    integer :: n, x, y, iq, target_x_index, target_y_index
    integer :: compare_timestep

    compare_timestep = current_state%timestep + current_state%reconfig_timestep_offset

    if (cfl_is_enabled) then
      if ((mod(compare_timestep, current_state%cfl_frequency) == 1 &
          .or. compare_timestep - current_state%start_timestep .le. current_state%cfl_frequency) &
          .or. compare_timestep .ge. (current_state%last_cfl_timestep + current_state%cfl_frequency)) then
        determine_flow_minmax=.true.
        call reset_local_minmax_values(current_state)
      else
        determine_flow_minmax=.false.
      end if
    end if

    if (is_component_enabled(current_state%options_database, "tvd_advection_acc")) then
      source_field_accelerated_info = get_component_field_information(current_state, "tvd_advection_acc_sth")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "tvd_advection_acc_sth")
        current_state%sth%data = current_state%sth%data + source_field_accelerated%real_3d_array_ptr
      end if
    end if

    if (is_component_enabled(current_state%options_database, "pw_advection_acc")) then
      source_field_accelerated_info = get_component_field_information(current_state, "pw_advection_acc_su")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "pw_advection_acc_su")
        current_state%su%data = current_state%su%data + source_field_accelerated%real_3d_array_ptr
      end if

      source_field_accelerated_info = get_component_field_information(current_state, "pw_advection_acc_sv")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "pw_advection_acc_sv")
        current_state%sv%data = current_state%sv%data + source_field_accelerated%real_3d_array_ptr
      end if

      source_field_accelerated_info = get_component_field_information(current_state, "pw_advection_acc_sw")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "pw_advection_acc_sw")
        current_state%sw%data = current_state%sw%data + source_field_accelerated%real_3d_array_ptr
      end if

      source_field_accelerated_info = get_component_field_information(current_state, "pw_advection_acc_sth")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "pw_advection_acc_sth")
        current_state%sth%data = current_state%sth%data + source_field_accelerated%real_3d_array_ptr
      end if

      source_field_accelerated_info = get_component_field_information(current_state, "pw_advection_acc_sq")
      if (source_field_accelerated_info%enabled) then
        source_field_accelerated = get_component_field_value(current_state, "pw_advection_acc_sq")
        do n = 1, current_state%number_q_fields
          current_state%sq(n)%data = current_state%sq(n)%data + source_field_accelerated%real_4d_array_ptr(:, :, :, n)
        end do
      end if
    end if

    ! Save vertical profiles once per timestep for swapsmooth Robert filter
#ifdef U_ACTIVE
    current_state%global_grid%configuration%vertical%savolubar = &
        current_state%global_grid%configuration%vertical%olzubar
#endif
#ifdef V_ACTIVE
    current_state%global_grid%configuration%vertical%savolvbar = &
        current_state%global_grid%configuration%vertical%olzvbar
#endif

    ! Zero profile tendency totals before looping over domain
    if (l_tend_pr_tot_th)   tend_pr_tot_th(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_qv)   tend_pr_tot_qv(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_ql)   tend_pr_tot_ql(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_qi)   tend_pr_tot_qi(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_qr)   tend_pr_tot_qr(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_qs)   tend_pr_tot_qs(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_qg)   tend_pr_tot_qg(:)   = 0.0_DEFAULT_PRECISION
    if (l_tend_pr_tot_tabs) tend_pr_tot_tabs(:) = 0.0_DEFAULT_PRECISION

    do x=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
      do y=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
        if (determine_flow_minmax .and. cfl_is_enabled) call determine_local_flow_minmax(current_state, y, x)
        call step_all_fields(current_state, y, x)

        if (current_state%diagnostic_sample_timestep) then
          target_y_index = y - current_state%local_grid%halo_size(Y_INDEX)
          target_x_index = x - current_state%local_grid%halo_size(X_INDEX)
          call compute_component_tendencies(current_state, x, y, target_x_index, target_y_index)
        end if
      end do
    end do

    ! Remove negative rounding errors
    if (l_nonconservative_positive_q) then
      do iq=1,current_state%number_q_fields
        resetq_min(iq)=minval(current_state%zq(iq)%data)
        current_state%zq(iq)%data = max(current_state%zq(iq)%data, 0.0_DEFAULT_PRECISION)
      end do
    end if

  end subroutine timestep_callback

  subroutine step_all_fields(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_x, local_y

    integer :: x_prev, y_prev, i
    real(kind=DEFAULT_PRECISION) :: c1, c2

    x_prev = local_x-2
    y_prev = local_y-1

    c1 = 1.0_DEFAULT_PRECISION - 2.0_DEFAULT_PRECISION*current_state%tsmth
    c2 = current_state%tsmth

#ifdef U_ACTIVE
    call step_single_field(local_x, local_y, &
         x_prev, y_prev, current_state%u, current_state%zu, current_state%su, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%ugal, c1, c2, .false., current_state%savu)
#endif
#ifdef V_ACTIVE
    call step_single_field(local_x, local_y, &
         x_prev, y_prev, current_state%v, current_state%zv, current_state%sv, current_state%local_grid, .true., &
         current_state%field_stepping, current_state%dtm, current_state%vgal, c1, c2, .false., current_state%savv)
#endif
#ifdef W_ACTIVE
    call step_single_field(local_x, local_y, &
         x_prev, y_prev, current_state%w, current_state%zw, current_state%sw, current_state%local_grid, .false., &
         current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, .false., current_state%savw)
#endif
    if (current_state%th%active) then
      call step_single_field(local_x, local_y, &
           x_prev, y_prev, current_state%th, current_state%zth, current_state%sth, current_state%local_grid, .false., &
           current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, &
           current_state%field_stepping == CENTRED_STEPPING)
    end if
    do i=1,current_state%number_q_fields
      if (current_state%q(i)%active) then
        call step_single_field(local_x, local_y, x_prev, y_prev, &
             current_state%q(i), current_state%zq(i), current_state%sq(i), current_state%local_grid, .false., &
             current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, &
             current_state%field_stepping == CENTRED_STEPPING)
      end if
    end do
    do i=1,current_state%n_tracers
      call step_single_field(local_x, local_y, x_prev, y_prev, &
           current_state%tracer(i), current_state%ztracer(i), current_state%stracer(i), current_state%local_grid, .false., &
           current_state%field_stepping, current_state%dtm, real(0., kind=DEFAULT_PRECISION), c1, c2, &
           current_state%field_stepping == CENTRED_STEPPING)
    end do
  end subroutine step_all_fields

  subroutine determine_local_flow_minmax(current_state, local_y, local_x)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: local_x, local_y

    integer :: k

    do k=2, current_state%local_grid%local_domain_end_index(Z_INDEX)
#ifdef U_ACTIVE
      current_state%local_zumax = max(current_state%local_zumax, current_state%zu%data(k,local_y,local_x))
      current_state%local_zumin = min(current_state%local_zumin, current_state%zu%data(k,local_y,local_x))
#endif
#ifdef V_ACTIVE
      current_state%local_zvmax = max(current_state%local_zvmax, current_state%zv%data(k,local_y,local_x))
      current_state%local_zvmin = min(current_state%local_zvmin, current_state%zv%data(k,local_y,local_x))
#endif
#ifdef W_ACTIVE
      if (k .lt. current_state%local_grid%local_domain_end_index(Z_INDEX)) then
        current_state%abswmax(k) = max(current_state%abswmax(k), abs(current_state%zw%data(k,local_y,local_x)))
      end if
#endif
    end do
  end subroutine determine_local_flow_minmax

  subroutine reset_local_minmax_values(current_state)
    type(model_state_type), intent(inout), target :: current_state

    current_state%local_zumin=rlargep
    current_state%local_zumax=-rlargep
    current_state%local_zvmin=rlargep
    current_state%local_zvmax=-rlargep
    current_state%abswmax=-rlargep
  end subroutine reset_local_minmax_values

  subroutine step_single_field(x_local_index, y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
    integer, intent(in) :: x_local_index, y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field, do_timesmoothing
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (x_prev .ge. local_grid%local_domain_start_index(X_INDEX)) then
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
      else
        call step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)
      end if
    end if

    if (x_local_index == local_grid%local_domain_end_index(X_INDEX)) then
      if (x_local_index .gt. 1) then
        if (present(sav)) then
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
        else
          call step_column_in_slice(y_local_index, x_local_index-1, y_prev, field, zfield, sfield, local_grid, &
               flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)
        end if
      end if
      if (present(sav)) then
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
      else
        call step_column_in_slice(y_local_index, x_local_index, y_prev, field, zfield, sfield, local_grid, &
             flow_field, direction, dtm, gal, c1, c2, do_timesmoothing)
      end if
    end if
  end subroutine step_single_field

  subroutine perform_timesmooth_for_field(field, zfield, local_grid, x_index, y_index, c1, c2)
    type(prognostic_field_type), intent(inout) :: field, zfield
    type(local_grid_type), intent(inout) :: local_grid
    integer, intent(in) :: x_index, y_index
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2

    integer :: k

    do k=1,local_grid%size(Z_INDEX)
      field%data(k, y_index, x_index)=c1*field%data(k, y_index, x_index)+c2*zfield%data(k, y_index, x_index)
    end do
  end subroutine perform_timesmooth_for_field

  subroutine step_column_in_slice(y_local_index, x_prev, y_prev, field, zfield, sfield, local_grid,&
       flow_field, direction, dtm, gal, c1, c2, do_timesmoothing, sav)
    integer, intent(in) :: y_local_index, x_prev, y_prev, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field, do_timesmoothing
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    real(kind=DEFAULT_PRECISION), intent(in) :: c1, c2
    type(prognostic_field_type), optional, intent(inout) :: sav

    if (y_prev .ge. local_grid%local_domain_start_index(Y_INDEX)) then
      if (do_timesmoothing) then
        call perform_timesmooth_for_field(field, zfield, local_grid, x_prev, y_prev, c1, c2)
      end if
      if (present(sav)) then
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_prev, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if

    if (y_local_index == local_grid%local_domain_end_index(Y_INDEX)) then
      if (do_timesmoothing) then
        call perform_timesmooth_for_field(field, zfield, local_grid, x_prev, y_local_index, c1, c2)
      end if
      if (present(sav)) then
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
      else
        call step_field(x_prev, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal)
      end if
    end if
  end subroutine step_column_in_slice

  subroutine step_field(x_local_index, y_local_index, field, zfield, sfield, local_grid, flow_field, direction, dtm, gal, sav)
    integer, intent(in) :: x_local_index, y_local_index, direction
    real(kind=DEFAULT_PRECISION), intent(in) :: dtm, gal
    logical, intent(in) :: flow_field
    type(local_grid_type), intent(inout) :: local_grid
    type(prognostic_field_type), intent(inout) :: field, zfield, sfield
    type(prognostic_field_type), optional, intent(inout) :: sav

    integer :: k
    real(kind=DEFAULT_PRECISION) :: actual_gal, dtm_x2

    dtm_x2 = 2.0_DEFAULT_PRECISION * dtm
    actual_gal = merge(gal, real(0.0_DEFAULT_PRECISION, kind=DEFAULT_PRECISION), flow_field)

    sfield%data(1,y_local_index, x_local_index)=0.0_DEFAULT_PRECISION

    do k=1,local_grid%size(Z_INDEX)
      if (present(sav) .and. direction .eq. CENTRED_STEPPING) &
           sav%data(k,y_local_index, x_local_index) = zfield%data(k, y_local_index, x_local_index) + actual_gal
      if (flow_field) field%data(k, y_local_index, x_local_index) = actual_gal + field%data(k, y_local_index, x_local_index)
      if (direction == FORWARD_STEPPING) then
        zfield%data(k, y_local_index, x_local_index) = field%data(k, y_local_index, x_local_index) + dtm * &
             sfield%data(k, y_local_index, x_local_index)
      else
        zfield%data(k, y_local_index, x_local_index) = actual_gal+zfield%data(k, y_local_index, x_local_index)+dtm_x2*&
             sfield%data(k, y_local_index, x_local_index)
      end if
    end do
  end subroutine step_field

  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: cxn, cyn, txn, tyn

    if (l_tend_3d_th) tend_3d_th(:,tyn,txn)=current_state%sth%data(:,cyn,cxn)
    if (l_tend_3d_qv) tend_3d_qv(:,tyn,txn)=current_state%sq(iqv)%data(:,cyn,cxn)
    if (l_tend_3d_ql) tend_3d_ql(:,tyn,txn)=current_state%sq(iql)%data(:,cyn,cxn)
    if (l_tend_3d_qi) tend_3d_qi(:,tyn,txn)=current_state%sq(iqi)%data(:,cyn,cxn)
    if (l_tend_3d_qr) tend_3d_qr(:,tyn,txn)=current_state%sq(iqr)%data(:,cyn,cxn)
    if (l_tend_3d_qs) tend_3d_qs(:,tyn,txn)=current_state%sq(iqs)%data(:,cyn,cxn)
    if (l_tend_3d_qg) tend_3d_qg(:,tyn,txn)=current_state%sq(iqg)%data(:,cyn,cxn)
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn)= &
          current_state%sth%data(:,cyn,cxn) * current_state%global_grid%configuration%vertical%rprefrcp(:)
    endif

    if (l_tend_pr_tot_th)   tend_pr_tot_th(:)   = tend_pr_tot_th(:)   + tend_3d_th(:,tyn,txn)
    if (l_tend_pr_tot_qv)   tend_pr_tot_qv(:)   = tend_pr_tot_qv(:)   + tend_3d_qv(:,tyn,txn)
    if (l_tend_pr_tot_ql)   tend_pr_tot_ql(:)   = tend_pr_tot_ql(:)   + tend_3d_ql(:,tyn,txn)
    if (l_tend_pr_tot_qi)   tend_pr_tot_qi(:)   = tend_pr_tot_qi(:)   + tend_3d_qi(:,tyn,txn)
    if (l_tend_pr_tot_qr)   tend_pr_tot_qr(:)   = tend_pr_tot_qr(:)   + tend_3d_qr(:,tyn,txn)
    if (l_tend_pr_tot_qs)   tend_pr_tot_qs(:)   = tend_pr_tot_qs(:)   + tend_3d_qs(:,tyn,txn)
    if (l_tend_pr_tot_qg)   tend_pr_tot_qg(:)   = tend_pr_tot_qg(:)   + tend_3d_qg(:,tyn,txn)
    if (l_tend_pr_tot_tabs) tend_pr_tot_tabs(:) = tend_pr_tot_tabs(:) + tend_3d_tabs(:,tyn,txn)
  end subroutine compute_component_tendencies

  subroutine field_information_retrieval_callback(current_state, name, field_information)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_information
    integer :: strcomp

    strcomp=INDEX(name, "_total_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=3
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2)=current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3)=current_state%local_grid%size(X_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_total_3d_local") then
        field_information%enabled=l_tend_3d_th
      else if (name .eq. "tend_qv_total_3d_local") then
        field_information%enabled=l_tend_3d_qv
      else if (name .eq. "tend_ql_total_3d_local") then
        field_information%enabled=l_tend_3d_ql
      else if (name .eq. "tend_qi_total_3d_local") then
        field_information%enabled=l_tend_3d_qi
      else if (name .eq. "tend_qr_total_3d_local") then
        field_information%enabled=l_tend_3d_qr
      else if (name .eq. "tend_qs_total_3d_local") then
        field_information%enabled=l_tend_3d_qs
      else if (name .eq. "tend_qg_total_3d_local") then
        field_information%enabled=l_tend_3d_qg
      else if (name .eq. "tend_tabs_total_3d_local") then
        field_information%enabled=l_tend_3d_tabs
      else
        field_information%enabled=.true.
      end if
    end if

    strcomp=INDEX(name, "_total_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type=COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions=1
      field_information%dimension_sizes(1)=current_state%local_grid%size(Z_INDEX)
      field_information%data_type=COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_th_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_th
      else if (name .eq. "tend_qv_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_total_profile_total_local") then
        field_information%enabled=l_tend_pr_tot_tabs
      else
        field_information%enabled=.true.
      end if
    end if
  end subroutine field_information_retrieval_callback

  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    if      (name .eq. "tend_th_total_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_total_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_total_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_total_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_total_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_total_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_total_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_total_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)
    else if (name .eq. "tend_th_total_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_total_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_total_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_total_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_total_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_total_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_total_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_total_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_tabs)
    end if
  end subroutine field_value_retrieval_callback

  subroutine set_published_field_value(field_value, real_1d_field, real_2d_field, real_3d_field)
    type(component_field_value_type), intent(inout) :: field_value
    real(kind=DEFAULT_PRECISION), dimension(:), optional :: real_1d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:), optional :: real_2d_field
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), optional :: real_3d_field

    if (present(real_1d_field)) then
      allocate(field_value%real_1d_array(size(real_1d_field)), source=real_1d_field)
    else if (present(real_2d_field)) then
      allocate(field_value%real_2d_array(size(real_2d_field, 1), size(real_2d_field, 2)), source=real_2d_field)
    else if (present(real_3d_field)) then
      allocate(field_value%real_3d_array(size(real_3d_field, 1), size(real_3d_field, 2), size(real_3d_field, 3)), &
               source=real_3d_field)
    end if
  end subroutine set_published_field_value

end module stepfields_mod
