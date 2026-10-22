!> Piacsek-Williams advection scheme with synchronous OpenMP Target offload
module pwadvectionompt_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, &
      COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, &
      component_field_value_type, component_field_information_type
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use collections_mod, only : map_type
  use state_mod, only : model_state_type
  use grids_mod, only : Z_INDEX, Y_INDEX, X_INDEX
  use q_indices_mod, only : get_q_index, standard_q_names
  use logging_mod, only : LOG_INFO, log_master_log

  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: advect_flow, advect_th, advect_q
  logical :: l_toplevel = .true.

  ! Persistent host/device buffers for OpenMP Target offload
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: u, v, w
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: su, sv, sw
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: th, sth
  real(kind=DEFAULT_PRECISION), dimension(:,:,:,:), allocatable :: q, sq

  integer :: number_q_fields
  integer, dimension(3) :: start_index, end_index
  real(kind=DEFAULT_PRECISION) :: tcx, tcy, cx, cy
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tzc1, tzc2, tzd1, tzd2

  ! Local tendency diagnostic variables for this component
  real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: &
       tend_3d_u, tend_3d_v, tend_3d_w, tend_3d_th, tend_3d_qv, &
       tend_3d_ql, tend_3d_qi, tend_3d_qr, tend_3d_qs, tend_3d_qg, &
       tend_3d_tabs
  logical :: l_tend_3d_u, l_tend_3d_v, l_tend_3d_w, l_tend_3d_th, l_tend_3d_qv, &
             l_tend_3d_ql, l_tend_3d_qi, l_tend_3d_qr, l_tend_3d_qs, l_tend_3d_qg, &
             l_tend_3d_tabs
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: &
       tend_pr_tot_u, tend_pr_tot_v, tend_pr_tot_w, tend_pr_tot_th, tend_pr_tot_qv, &
       tend_pr_tot_ql, tend_pr_tot_qi, tend_pr_tot_qr, tend_pr_tot_qs, tend_pr_tot_qg, &
       tend_pr_tot_tabs
  logical :: l_tend_pr_tot_u, l_tend_pr_tot_v, l_tend_pr_tot_w, l_tend_pr_tot_th, l_tend_pr_tot_qv, &
             l_tend_pr_tot_ql, l_tend_pr_tot_qi, l_tend_pr_tot_qr, l_tend_pr_tot_qs, l_tend_pr_tot_qg, &
             l_tend_pr_tot_tabs
  integer :: iqv=0, iql=0, iqr=0, iqi=0, iqs=0, iqg=0

  public pwadvectionompt_get_descriptor

contains

  !> Provides a description of this component for the core to register
  type(component_descriptor_type) function pwadvectionompt_get_descriptor()
    pwadvectionompt_get_descriptor%name="pw_advection_ompt"
    pwadvectionompt_get_descriptor%version=0.1
    pwadvectionompt_get_descriptor%initialisation=>initialisation_callback
    pwadvectionompt_get_descriptor%timestep=>timestep_callback
    pwadvectionompt_get_descriptor%finalisation=>finalisation_callback

    pwadvectionompt_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    pwadvectionompt_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(pwadvectionompt_get_descriptor%published_fields(11+11))

    pwadvectionompt_get_descriptor%published_fields(1)= "tend_u_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(2)= "tend_v_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(3)= "tend_w_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(4)= "tend_th_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(5)= "tend_qv_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(6)= "tend_ql_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(7)= "tend_qi_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(8)= "tend_qr_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(9)= "tend_qs_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(10)="tend_qg_pwadvection_3d_local"
    pwadvectionompt_get_descriptor%published_fields(11)="tend_tabs_pwadvection_3d_local"

    pwadvectionompt_get_descriptor%published_fields(11+1)= "tend_u_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+2)= "tend_v_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+3)= "tend_w_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+4)= "tend_th_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+5)= "tend_qv_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+6)= "tend_ql_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+7)= "tend_qi_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+8)= "tend_qr_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+9)= "tend_qs_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+10)="tend_qg_pwadvection_profile_total_local"
    pwadvectionompt_get_descriptor%published_fields(11+11)="tend_tabs_pwadvection_profile_total_local"
  end function pwadvectionompt_get_descriptor

  !> Initialisation callback, sets up configuration and persistent OpenMP device buffers
  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    logical :: l_qdiag
    integer :: sz, sy, sx

    number_q_fields = current_state%number_q_fields

    advect_flow = determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))
    advect_th   = determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q    = determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))

    advect_th = advect_th .and. current_state%th%active
    advect_q  = advect_q .and. (number_q_fields .gt. 0)

    start_index = current_state%local_grid%local_domain_start_index
    end_index   = current_state%local_grid%local_domain_end_index

    if (advect_flow .or. advect_th .or. advect_q) then
      sz = size(current_state%u%data, 1)
      sy = size(current_state%u%data, 2)
      sx = size(current_state%u%data, 3)

      allocate(u(sz, sy, sx), v(sz, sy, sx), w(sz, sy, sx))

      tcx = current_state%global_grid%configuration%horizontal%tcx
      tcy = current_state%global_grid%configuration%horizontal%tcy
      cx  = current_state%global_grid%configuration%horizontal%cx
      cy  = current_state%global_grid%configuration%horizontal%cy

      allocate(tzc1, source=current_state%global_grid%configuration%vertical%tzc1)
      allocate(tzc2, source=current_state%global_grid%configuration%vertical%tzc2)
      allocate(tzd1, source=current_state%global_grid%configuration%vertical%tzd1)
      allocate(tzd2, source=current_state%global_grid%configuration%vertical%tzd2)

      !$omp target enter data map(alloc: u, v, w) map(to: tzc1, tzc2, tzd1, tzd2)
    end if

    if (advect_flow) then
      call log_master_log(LOG_INFO, "Accelerating OpenMP PW advection: flow fields")
      allocate(su(size(current_state%su%data, 1), size(current_state%su%data, 2), size(current_state%su%data, 3)))
      allocate(sv(size(current_state%sv%data, 1), size(current_state%sv%data, 2), size(current_state%sv%data, 3)))
      allocate(sw(size(current_state%sw%data, 1), size(current_state%sw%data, 2), size(current_state%sw%data, 3)))
      su = 0.0_DEFAULT_PRECISION
      sv = 0.0_DEFAULT_PRECISION
      sw = 0.0_DEFAULT_PRECISION
      !$omp target enter data map(to: su, sv, sw)
    end if

    if (advect_th) then
      call log_master_log(LOG_INFO, "Accelerating OpenMP PW advection: theta field")
      allocate(th(size(current_state%th%data, 1), size(current_state%th%data, 2), size(current_state%th%data, 3)))
      allocate(sth(size(current_state%sth%data, 1), size(current_state%sth%data, 2), size(current_state%sth%data, 3)))
      th  = 0.0_DEFAULT_PRECISION
      sth = 0.0_DEFAULT_PRECISION
      !$omp target enter data map(alloc: th) map(to: sth)
    end if

    if (advect_q) then
      call log_master_log(LOG_INFO, "Accelerating OpenMP PW advection: q fields")
      allocate(q(size(current_state%q(1)%data, 1), size(current_state%q(1)%data, 2), &
                 size(current_state%q(1)%data, 3), number_q_fields))
      allocate(sq(size(current_state%sq(1)%data, 1), size(current_state%sq(1)%data, 2), &
                  size(current_state%sq(1)%data, 3), number_q_fields))
      q  = 0.0_DEFAULT_PRECISION
      sq = 0.0_DEFAULT_PRECISION
      !$omp target enter data map(alloc: q) map(to: sq)
    end if

    ! Set tendency diagnostic logicals based on availability
    l_qdiag = (.not. current_state%passive_q .and. current_state%number_q_fields .gt. 0) .and. advect_q

    l_tend_pr_tot_u    = current_state%u%active .and. advect_flow
    l_tend_pr_tot_v    = current_state%v%active .and. advect_flow
    l_tend_pr_tot_w    = current_state%w%active .and. advect_flow
    l_tend_pr_tot_th   = current_state%th%active .and. advect_th
    l_tend_pr_tot_qv   = l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0
    l_tend_pr_tot_ql   = l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0
    l_tend_pr_tot_qi   = l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0
    l_tend_pr_tot_qr   = l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0
    l_tend_pr_tot_qs   = l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0
    l_tend_pr_tot_qg   = l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0
    l_tend_pr_tot_tabs = l_tend_pr_tot_th

    l_tend_3d_u    = (current_state%u%active .and. advect_flow) .or. l_tend_pr_tot_u
    l_tend_3d_v    = (current_state%v%active .and. advect_flow) .or. l_tend_pr_tot_v
    l_tend_3d_w    = (current_state%w%active .and. advect_flow) .or. l_tend_pr_tot_w
    l_tend_3d_th   = (current_state%th%active .and. advect_th)  .or. l_tend_pr_tot_th
    l_tend_3d_qv   = (l_qdiag .and. current_state%water_vapour_mixing_ratio_index > 0) .or. l_tend_pr_tot_qv
    l_tend_3d_ql   = (l_qdiag .and. current_state%liquid_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_ql
    l_tend_3d_qi   = (l_qdiag .and. current_state%ice_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qi
    l_tend_3d_qr   = (l_qdiag .and. current_state%rain_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qr
    l_tend_3d_qs   = (l_qdiag .and. current_state%snow_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qs
    l_tend_3d_qg   = (l_qdiag .and. current_state%graupel_water_mixing_ratio_index > 0) .or. l_tend_pr_tot_qg
    l_tend_3d_tabs = l_tend_3d_th

    if (l_tend_3d_u)  allocate(tend_3d_u(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    if (l_tend_3d_v)  allocate(tend_3d_v(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    if (l_tend_3d_w)  allocate(tend_3d_w(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    if (l_tend_3d_th) allocate(tend_3d_th(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    if (l_tend_3d_qv) then
      iqv=get_q_index(standard_q_names%VAPOUR, 'pw_advection_ompt')
      allocate(tend_3d_qv(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_ql) then
      iql=get_q_index(standard_q_names%CLOUD_LIQUID_MASS, 'pw_advection_ompt')
      allocate(tend_3d_ql(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qi) then
      iqi=get_q_index(standard_q_names%ICE_MASS, 'pw_advection_ompt')
      allocate(tend_3d_qi(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qr) then
      iqr=get_q_index(standard_q_names%RAIN_MASS, 'pw_advection_ompt')
      allocate(tend_3d_qr(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qs) then
      iqs=get_q_index(standard_q_names%SNOW_MASS, 'pw_advection_ompt')
      allocate(tend_3d_qs(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_qg) then
      iqg=get_q_index(standard_q_names%GRAUPEL_MASS, 'pw_advection_ompt')
      allocate(tend_3d_qg(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))
    endif
    if (l_tend_3d_tabs) allocate(tend_3d_tabs(current_state%local_grid%size(Z_INDEX), current_state%local_grid%size(Y_INDEX), current_state%local_grid%size(X_INDEX)))

    if (l_tend_pr_tot_u)    allocate(tend_pr_tot_u(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_v)    allocate(tend_pr_tot_v(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_w)    allocate(tend_pr_tot_w(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_th)   allocate(tend_pr_tot_th(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qv)   allocate(tend_pr_tot_qv(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_ql)   allocate(tend_pr_tot_ql(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qi)   allocate(tend_pr_tot_qi(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qr)   allocate(tend_pr_tot_qr(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qs)   allocate(tend_pr_tot_qs(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_qg)   allocate(tend_pr_tot_qg(current_state%local_grid%size(Z_INDEX)))
    if (l_tend_pr_tot_tabs) allocate(tend_pr_tot_tabs(current_state%local_grid%size(Z_INDEX)))
  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (advect_q) then
      !$omp target exit data map(delete: q, sq)
      if (allocated(q))  deallocate(q)
      if (allocated(sq)) deallocate(sq)
    end if

    if (advect_th) then
      !$omp target exit data map(delete: th, sth)
      if (allocated(th))  deallocate(th)
      if (allocated(sth)) deallocate(sth)
    end if

    if (advect_flow) then
      !$omp target exit data map(delete: su, sv, sw)
      if (allocated(su)) deallocate(su)
      if (allocated(sv)) deallocate(sv)
      if (allocated(sw)) deallocate(sw)
    end if

    if (advect_flow .or. advect_th .or. advect_q) then
      !$omp target exit data map(delete: u, v, w, tzc1, tzc2, tzd1, tzd2)
      if (allocated(u))    deallocate(u)
      if (allocated(v))    deallocate(v)
      if (allocated(w))    deallocate(w)
      if (allocated(tzc1)) deallocate(tzc1)
      if (allocated(tzc2)) deallocate(tzc2)
      if (allocated(tzd1)) deallocate(tzd1)
      if (allocated(tzd2)) deallocate(tzd2)
    end if

    if (allocated(tend_3d_u))    deallocate(tend_3d_u)
    if (allocated(tend_3d_v))    deallocate(tend_3d_v)
    if (allocated(tend_3d_w))    deallocate(tend_3d_w)
    if (allocated(tend_3d_th))   deallocate(tend_3d_th)
    if (allocated(tend_3d_qv))   deallocate(tend_3d_qv)
    if (allocated(tend_3d_ql))   deallocate(tend_3d_ql)
    if (allocated(tend_3d_qi))   deallocate(tend_3d_qi)
    if (allocated(tend_3d_qr))   deallocate(tend_3d_qr)
    if (allocated(tend_3d_qs))   deallocate(tend_3d_qs)
    if (allocated(tend_3d_qg))   deallocate(tend_3d_qg)
    if (allocated(tend_3d_tabs)) deallocate(tend_3d_tabs)

    if (allocated(tend_pr_tot_u))    deallocate(tend_pr_tot_u)
    if (allocated(tend_pr_tot_v))    deallocate(tend_pr_tot_v)
    if (allocated(tend_pr_tot_w))    deallocate(tend_pr_tot_w)
    if (allocated(tend_pr_tot_th))   deallocate(tend_pr_tot_th)
    if (allocated(tend_pr_tot_qv))   deallocate(tend_pr_tot_qv)
    if (allocated(tend_pr_tot_ql))   deallocate(tend_pr_tot_ql)
    if (allocated(tend_pr_tot_qi))   deallocate(tend_pr_tot_qi)
    if (allocated(tend_pr_tot_qr))   deallocate(tend_pr_tot_qr)
    if (allocated(tend_pr_tot_qs))   deallocate(tend_pr_tot_qs)
    if (allocated(tend_pr_tot_qg))   deallocate(tend_pr_tot_qg)
    if (allocated(tend_pr_tot_tabs)) deallocate(tend_pr_tot_tabs)
  end subroutine finalisation_callback

  !> Called once per timestep in the 'dynamics-acc' (entire domain) group
  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: x, y, n, target_x_index, target_y_index
    integer :: sx_idx, ex_idx, sy_idx, ey_idx, k_top
    logical :: calculate_diagnostics

    calculate_diagnostics = current_state%diagnostic_sample_timestep

    sx_idx = start_index(X_INDEX)
    ex_idx = end_index(X_INDEX)
    sy_idx = start_index(Y_INDEX)
    ey_idx = end_index(Y_INDEX)
    k_top  = merge(end_index(Z_INDEX), end_index(Z_INDEX) - 1, l_toplevel)

    if (calculate_diagnostics) then
      if (l_tend_pr_tot_u)    tend_pr_tot_u(:)    = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_v)    tend_pr_tot_v(:)    = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_w)    tend_pr_tot_w(:)    = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_th)   tend_pr_tot_th(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_qv)   tend_pr_tot_qv(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_ql)   tend_pr_tot_ql(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_qi)   tend_pr_tot_qi(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_qr)   tend_pr_tot_qr(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_qs)   tend_pr_tot_qs(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_qg)   tend_pr_tot_qg(:)   = 0.0_DEFAULT_PRECISION
      if (l_tend_pr_tot_tabs) tend_pr_tot_tabs(:) = 0.0_DEFAULT_PRECISION

      do x = sx_idx, ex_idx
        do y = sy_idx, ey_idx
          target_y_index = y - current_state%local_grid%halo_size(Y_INDEX)
          target_x_index = x - current_state%local_grid%halo_size(X_INDEX)
          call save_precomponent_tendencies(current_state, x, y, target_x_index, target_y_index)
        end do
      end do
    end if

    if (advect_flow .or. advect_th .or. advect_q) then
      u = current_state%u%data
      v = current_state%v%data
      w = current_state%w%data
      !$omp target update to(u, v, w)
    end if

    if (advect_flow) then
      call advect_flow_fields_ompt(size(u, 1), size(u, 2), size(u, 3), &
                                   sx_idx, ex_idx, sy_idx, ey_idx, end_index(Z_INDEX), &
                                   u, v, w, su, sv, sw, tcx, tcy, tzc1, tzc2, tzd1, tzd2, l_toplevel)
#ifdef U_ACTIVE
      current_state%su%data(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx) = &
           su(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx)
#endif
#ifdef V_ACTIVE
      current_state%sv%data(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx) = &
           sv(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx)
#endif
#ifdef W_ACTIVE
      current_state%sw%data(2:end_index(Z_INDEX)-1, sy_idx:ey_idx, sx_idx:ex_idx) = &
           sw(2:end_index(Z_INDEX)-1, sy_idx:ey_idx, sx_idx:ex_idx)
#endif
    end if

    if (advect_th) then
      th = current_state%th%data
      !$omp target update to(th)
      call advect_th_field_ompt(size(u, 1), size(u, 2), size(u, 3), &
                                sx_idx, ex_idx, sy_idx, ey_idx, end_index(Z_INDEX), &
                                u, v, w, th, sth, cx, cy, tzc1, tzc2, l_toplevel)
      current_state%sth%data(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx) = &
           sth(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx)
    end if

    if (advect_q) then
      do n = 1, number_q_fields
        q(:, :, :, n) = current_state%q(n)%data
      end do
      !$omp target update to(q)
      call advect_q_fields_ompt(size(u, 1), size(u, 2), size(u, 3), number_q_fields, &
                                sx_idx, ex_idx, sy_idx, ey_idx, end_index(Z_INDEX), &
                                u, v, w, q, sq, cx, cy, tzc1, tzc2, l_toplevel)
      do n = 1, number_q_fields
        current_state%sq(n)%data(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx) = &
             current_state%sq(n)%data(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx) + &
             sq(2:k_top, sy_idx:ey_idx, sx_idx:ex_idx, n)
      end do
    end if

    if (calculate_diagnostics) then
      do x = sx_idx, ex_idx
        do y = sy_idx, ey_idx
          target_y_index = y - current_state%local_grid%halo_size(Y_INDEX)
          target_x_index = x - current_state%local_grid%halo_size(X_INDEX)
          call compute_component_tendencies(current_state, x, y, target_x_index, target_y_index)
        end do
      end do
    end if
  end subroutine timestep_callback

  !> Synchronous OpenMP Target kernel for flow fields (u, v, w)
  subroutine advect_flow_fields_ompt(sz, sy, sx, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, &
       u_arr, v_arr, w_arr, su_arr, sv_arr, sw_arr, tcx_val, tcy_val, &
       tzc1_arr, tzc2_arr, tzd1_arr, tzd2_arr, l_top)
    integer, intent(in) :: sz, sy, sx, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx
    real(kind=DEFAULT_PRECISION), intent(in) :: tcx_val, tcy_val
    logical, intent(in) :: l_top
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx), intent(in) :: u_arr, v_arr, w_arr
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx), intent(inout) :: su_arr, sv_arr, sw_arr
    real(kind=DEFAULT_PRECISION), dimension(sz), intent(in) :: tzc1_arr, tzc2_arr, tzd1_arr, tzd2_arr

    integer :: x, y, k

    !$omp target map(alloc: u_arr, v_arr, w_arr, su_arr, sv_arr, sw_arr, tzc1_arr, tzc2_arr, tzd1_arr, tzd2_arr) &
    !$omp map(to: sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, tcx_val, tcy_val)
    !$omp teams distribute parallel do collapse(3)
    do x = sx_idx, ex_idx
      do y = sy_idx, ey_idx
        do k = 2, ez_idx - 1
#ifdef U_ACTIVE
          su_arr(k, y, x) = tcx_val * (u_arr(k, y, x-1) * (u_arr(k, y, x) + u_arr(k, y, x-1)) - &
                                       u_arr(k, y, x+1) * (u_arr(k, y, x) + u_arr(k, y, x+1)))
#ifdef V_ACTIVE
          su_arr(k, y, x) = su_arr(k, y, x) + &
                            tcy_val * (u_arr(k, y-1, x) * (v_arr(k, y-1, x) + v_arr(k, y-1, x+1)) - &
                                       u_arr(k, y+1, x) * (v_arr(k, y, x)   + v_arr(k, y, x+1)))
#endif
#ifdef W_ACTIVE
          su_arr(k, y, x) = su_arr(k, y, x) + &
                            tzc1_arr(k) * u_arr(k-1, y, x) * (w_arr(k-1, y, x) + w_arr(k-1, y, x+1)) - &
                            tzc2_arr(k) * u_arr(k+1, y, x) * (w_arr(k, y, x)   + w_arr(k, y, x+1))
#endif
#endif

#ifdef V_ACTIVE
          sv_arr(k, y, x) = tcy_val * (v_arr(k, y-1, x) * (v_arr(k, y, x) + v_arr(k, y-1, x)) - &
                                       v_arr(k, y+1, x) * (v_arr(k, y, x) + v_arr(k, y+1, x)))
#ifdef U_ACTIVE
          sv_arr(k, y, x) = sv_arr(k, y, x) + &
                            tcx_val * (v_arr(k, y, x-1) * (u_arr(k, y, x-1) + u_arr(k, y+1, x-1)) - &
                                       v_arr(k, y, x+1) * (u_arr(k, y, x)   + u_arr(k, y+1, x)))
#endif
#ifdef W_ACTIVE
          sv_arr(k, y, x) = sv_arr(k, y, x) + &
                            tzc1_arr(k) * v_arr(k-1, y, x) * (w_arr(k-1, y, x) + w_arr(k-1, y+1, x)) - &
                            tzc2_arr(k) * v_arr(k+1, y, x) * (w_arr(k, y, x)   + w_arr(k, y+1, x))
#endif
#endif

#ifdef W_ACTIVE
          sw_arr(k, y, x) = tzd1_arr(k) * w_arr(k-1, y, x) * (w_arr(k, y, x) + w_arr(k-1, y, x)) - &
                            tzd2_arr(k) * w_arr(k+1, y, x) * (w_arr(k, y, x) + w_arr(k+1, y, x))
#ifdef U_ACTIVE
          sw_arr(k, y, x) = sw_arr(k, y, x) + &
                            tcx_val * (w_arr(k, y, x-1) * (u_arr(k, y, x-1) + u_arr(k+1, y, x-1)) - &
                                       w_arr(k, y, x+1) * (u_arr(k, y, x)   + u_arr(k+1, y, x)))
#endif
#ifdef V_ACTIVE
          sw_arr(k, y, x) = sw_arr(k, y, x) + &
                            tcy_val * (w_arr(k, y-1, x) * (v_arr(k, y-1, x) + v_arr(k+1, y-1, x)) - &
                                       w_arr(k, y+1, x) * (v_arr(k, y, x)   + v_arr(k+1, y, x)))
#endif
#endif
        end do
      end do
    end do
    !$omp end target

    if (l_top) then
      k = ez_idx
      !$omp target map(alloc: u_arr, v_arr, w_arr, su_arr, sv_arr, tzc1_arr) &
      !$omp map(to: sx_idx, ex_idx, sy_idx, ey_idx, k, tcx_val, tcy_val)
      !$omp teams distribute parallel do collapse(2)
      do x = sx_idx, ex_idx
        do y = sy_idx, ey_idx
#ifdef U_ACTIVE
          su_arr(k, y, x) = tcx_val * (u_arr(k, y, x-1) * (u_arr(k, y, x) + u_arr(k, y, x-1)) - &
                                       u_arr(k, y, x+1) * (u_arr(k, y, x) + u_arr(k, y, x+1)))
#ifdef V_ACTIVE
          su_arr(k, y, x) = su_arr(k, y, x) + &
                            tcy_val * (u_arr(k, y-1, x) * (v_arr(k, y-1, x) + v_arr(k, y-1, x+1)) - &
                                       u_arr(k, y+1, x) * (v_arr(k, y, x)   + v_arr(k, y, x+1)))
#endif
#ifdef W_ACTIVE
          su_arr(k, y, x) = su_arr(k, y, x) + &
                            tzc1_arr(k) * u_arr(k-1, y, x) * (w_arr(k-1, y, x) + w_arr(k-1, y, x+1))
#endif
#endif

#ifdef V_ACTIVE
          sv_arr(k, y, x) = tcy_val * (v_arr(k, y-1, x) * (v_arr(k, y, x) + v_arr(k, y-1, x)) - &
                                       v_arr(k, y+1, x) * (v_arr(k, y, x) + v_arr(k, y+1, x)))
#ifdef U_ACTIVE
          sv_arr(k, y, x) = sv_arr(k, y, x) + &
                            tcx_val * (v_arr(k, y, x-1) * (u_arr(k, y, x-1) + u_arr(k, y+1, x-1)) - &
                                       v_arr(k, y, x+1) * (u_arr(k, y, x)   + u_arr(k, y+1, x)))
#endif
#ifdef W_ACTIVE
          sv_arr(k, y, x) = sv_arr(k, y, x) + &
                            tzc1_arr(k) * v_arr(k-1, y, x) * (w_arr(k-1, y, x) + w_arr(k-1, y+1, x))
#endif
#endif
        end do
      end do
      !$omp end target
    end if

    !$omp target update from(su_arr, sv_arr, sw_arr)
  end subroutine advect_flow_fields_ompt

  !> Synchronous OpenMP Target kernel for theta field
  subroutine advect_th_field_ompt(sz, sy, sx, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, &
       u_arr, v_arr, w_arr, th_arr, sth_arr, cx_val, cy_val, tzc1_arr, tzc2_arr, l_top)
    integer, intent(in) :: sz, sy, sx, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx
    real(kind=DEFAULT_PRECISION), intent(in) :: cx_val, cy_val
    logical, intent(in) :: l_top
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx), intent(in) :: u_arr, v_arr, w_arr, th_arr
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx), intent(inout) :: sth_arr
    real(kind=DEFAULT_PRECISION), dimension(sz), intent(in) :: tzc1_arr, tzc2_arr

    integer :: x, y, k

    !$omp target map(alloc: u_arr, v_arr, w_arr, th_arr, sth_arr, tzc1_arr, tzc2_arr) &
    !$omp map(to: sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, cx_val, cy_val)
    !$omp teams distribute parallel do collapse(3)
    do x = sx_idx, ex_idx
      do y = sy_idx, ey_idx
        do k = 2, ez_idx - 1
          sth_arr(k, y, x) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 0.5_DEFAULT_PRECISION * cx_val * &
                             (u_arr(k, y, x-1) * th_arr(k, y, x-1) - u_arr(k, y, x) * th_arr(k, y, x+1))
#endif
#ifdef V_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 0.5_DEFAULT_PRECISION * cy_val * &
                             (v_arr(k, y-1, x) * th_arr(k, y-1, x) - v_arr(k, y, x) * th_arr(k, y+1, x))
#endif
#ifdef W_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 2.0_DEFAULT_PRECISION * &
                             (tzc1_arr(k) * w_arr(k-1, y, x) * th_arr(k-1, y, x) - &
                              tzc2_arr(k) * w_arr(k, y, x)   * th_arr(k+1, y, x))
#endif
        end do
      end do
    end do
    !$omp end target

    if (l_top) then
      k = ez_idx
      !$omp target map(alloc: u_arr, v_arr, w_arr, th_arr, sth_arr, tzc1_arr) &
      !$omp map(to: sx_idx, ex_idx, sy_idx, ey_idx, k, cx_val, cy_val)
      !$omp teams distribute parallel do collapse(2)
      do x = sx_idx, ex_idx
        do y = sy_idx, ey_idx
          sth_arr(k, y, x) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 0.5_DEFAULT_PRECISION * cx_val * &
                             (u_arr(k, y, x-1) * th_arr(k, y, x-1) - u_arr(k, y, x) * th_arr(k, y, x+1))
#endif
#ifdef V_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 0.5_DEFAULT_PRECISION * cy_val * &
                             (v_arr(k, y-1, x) * th_arr(k, y-1, x) - v_arr(k, y, x) * th_arr(k, y+1, x))
#endif
#ifdef W_ACTIVE
          sth_arr(k, y, x) = sth_arr(k, y, x) + 2.0_DEFAULT_PRECISION * &
                             tzc1_arr(k) * w_arr(k-1, y, x) * th_arr(k-1, y, x)
#endif
        end do
      end do
      !$omp end target
    end if

    !$omp target update from(sth_arr)
  end subroutine advect_th_field_ompt

  !> Synchronous OpenMP Target kernel for packed q moisture fields
  subroutine advect_q_fields_ompt(sz, sy, sx, nq, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, &
       u_arr, v_arr, w_arr, q_arr, sq_arr, cx_val, cy_val, tzc1_arr, tzc2_arr, l_top)
    integer, intent(in) :: sz, sy, sx, nq, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx
    real(kind=DEFAULT_PRECISION), intent(in) :: cx_val, cy_val
    logical, intent(in) :: l_top
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx), intent(in) :: u_arr, v_arr, w_arr
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx, nq), intent(in) :: q_arr
    real(kind=DEFAULT_PRECISION), dimension(sz, sy, sx, nq), intent(inout) :: sq_arr
    real(kind=DEFAULT_PRECISION), dimension(sz), intent(in) :: tzc1_arr, tzc2_arr

    integer :: n, x, y, k

    !$omp target map(alloc: u_arr, v_arr, w_arr, q_arr, sq_arr, tzc1_arr, tzc2_arr) &
    !$omp map(to: nq, sx_idx, ex_idx, sy_idx, ey_idx, ez_idx, cx_val, cy_val)
    !$omp teams distribute parallel do collapse(4)
    do n = 1, nq
      do x = sx_idx, ex_idx
        do y = sy_idx, ey_idx
          do k = 2, ez_idx - 1
            sq_arr(k, y, x, n) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 0.5_DEFAULT_PRECISION * cx_val * &
                                 (u_arr(k, y, x-1) * q_arr(k, y, x-1, n) - u_arr(k, y, x) * q_arr(k, y, x+1, n))
#endif
#ifdef V_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 0.5_DEFAULT_PRECISION * cy_val * &
                                 (v_arr(k, y-1, x) * q_arr(k, y-1, x, n) - v_arr(k, y, x) * q_arr(k, y+1, x, n))
#endif
#ifdef W_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 2.0_DEFAULT_PRECISION * &
                                 (tzc1_arr(k) * w_arr(k-1, y, x) * q_arr(k-1, y, x, n) - &
                                  tzc2_arr(k) * w_arr(k, y, x)   * q_arr(k+1, y, x, n))
#endif
          end do
        end do
      end do
    end do
    !$omp end target

    if (l_top) then
      k = ez_idx
      !$omp target map(alloc: u_arr, v_arr, w_arr, q_arr, sq_arr, tzc1_arr) &
      !$omp map(to: nq, sx_idx, ex_idx, sy_idx, ey_idx, k, cx_val, cy_val)
      !$omp teams distribute parallel do collapse(3)
      do n = 1, nq
        do x = sx_idx, ex_idx
          do y = sy_idx, ey_idx
            sq_arr(k, y, x, n) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 0.5_DEFAULT_PRECISION * cx_val * &
                                 (u_arr(k, y, x-1) * q_arr(k, y, x-1, n) - u_arr(k, y, x) * q_arr(k, y, x+1, n))
#endif
#ifdef V_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 0.5_DEFAULT_PRECISION * cy_val * &
                                 (v_arr(k, y-1, x) * q_arr(k, y-1, x, n) - v_arr(k, y, x) * q_arr(k, y+1, x, n))
#endif
#ifdef W_ACTIVE
            sq_arr(k, y, x, n) = sq_arr(k, y, x, n) + 2.0_DEFAULT_PRECISION * &
                                 tzc1_arr(k) * w_arr(k-1, y, x) * q_arr(k-1, y, x, n)
#endif
          end do
        end do
      end do
      !$omp end target
    end if

    !$omp target update from(sq_arr)
  end subroutine advect_q_fields_ompt

  !> Save the 3d tendencies coming into the component
  subroutine save_precomponent_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(in) :: current_state
    integer, intent(in) :: cxn, cyn, txn, tyn

    if (l_tend_3d_u)    tend_3d_u(:,tyn,txn)    = current_state%su%data(:,cyn,cxn)
    if (l_tend_3d_v)    tend_3d_v(:,tyn,txn)    = current_state%sv%data(:,cyn,cxn)
    if (l_tend_3d_w)    tend_3d_w(:,tyn,txn)    = current_state%sw%data(:,cyn,cxn)
    if (l_tend_3d_th)   tend_3d_th(:,tyn,txn)   = current_state%sth%data(:,cyn,cxn)
    if (l_tend_3d_qv)   tend_3d_qv(:,tyn,txn)   = current_state%sq(iqv)%data(:,cyn,cxn)
    if (l_tend_3d_ql)   tend_3d_ql(:,tyn,txn)   = current_state%sq(iql)%data(:,cyn,cxn)
    if (l_tend_3d_qi)   tend_3d_qi(:,tyn,txn)   = current_state%sq(iqi)%data(:,cyn,cxn)
    if (l_tend_3d_qr)   tend_3d_qr(:,tyn,txn)   = current_state%sq(iqr)%data(:,cyn,cxn)
    if (l_tend_3d_qs)   tend_3d_qs(:,tyn,txn)   = current_state%sq(iqs)%data(:,cyn,cxn)
    if (l_tend_3d_qg)   tend_3d_qg(:,tyn,txn)   = current_state%sq(iqg)%data(:,cyn,cxn)
    if (l_tend_3d_tabs) tend_3d_tabs(:,tyn,txn) = current_state%sth%data(:,cyn,cxn) * &
                                                  current_state%global_grid%configuration%vertical%rprefrcp(:)
  end subroutine save_precomponent_tendencies

  !> Computation of component tendencies
  subroutine compute_component_tendencies(current_state, cxn, cyn, txn, tyn)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: cxn, cyn, txn, tyn

    if (l_tend_3d_u)    tend_3d_u(:,tyn,txn)  = current_state%su%data(:,cyn,cxn)      - tend_3d_u(:,tyn,txn)
    if (l_tend_3d_v)    tend_3d_v(:,tyn,txn)  = current_state%sv%data(:,cyn,cxn)      - tend_3d_v(:,tyn,txn)
    if (l_tend_3d_w)    tend_3d_w(:,tyn,txn)  = current_state%sw%data(:,cyn,cxn)      - tend_3d_w(:,tyn,txn)
    if (l_tend_3d_th)   tend_3d_th(:,tyn,txn) = current_state%sth%data(:,cyn,cxn)     - tend_3d_th(:,tyn,txn)
    if (l_tend_3d_qv)   tend_3d_qv(:,tyn,txn) = current_state%sq(iqv)%data(:,cyn,cxn) - tend_3d_qv(:,tyn,txn)
    if (l_tend_3d_ql)   tend_3d_ql(:,tyn,txn) = current_state%sq(iql)%data(:,cyn,cxn) - tend_3d_ql(:,tyn,txn)
    if (l_tend_3d_qi)   tend_3d_qi(:,tyn,txn) = current_state%sq(iqi)%data(:,cyn,cxn) - tend_3d_qi(:,tyn,txn)
    if (l_tend_3d_qr)   tend_3d_qr(:,tyn,txn) = current_state%sq(iqr)%data(:,cyn,cxn) - tend_3d_qr(:,tyn,txn)
    if (l_tend_3d_qs)   tend_3d_qs(:,tyn,txn) = current_state%sq(iqs)%data(:,cyn,cxn) - tend_3d_qs(:,tyn,txn)
    if (l_tend_3d_qg)   tend_3d_qg(:,tyn,txn) = current_state%sq(iqg)%data(:,cyn,cxn) - tend_3d_qg(:,tyn,txn)
    if (l_tend_3d_tabs) then
      tend_3d_tabs(:,tyn,txn) = current_state%sth%data(:,cyn,cxn) * &
           current_state%global_grid%configuration%vertical%rprefrcp(:) - tend_3d_tabs(:,tyn,txn)
    endif

    if (l_tend_pr_tot_u)    tend_pr_tot_u(:)    = tend_pr_tot_u(:)    + tend_3d_u(:,tyn,txn)
    if (l_tend_pr_tot_v)    tend_pr_tot_v(:)    = tend_pr_tot_v(:)    + tend_3d_v(:,tyn,txn)
    if (l_tend_pr_tot_w)    tend_pr_tot_w(:)    = tend_pr_tot_w(:)    + tend_3d_w(:,tyn,txn)
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

    strcomp = INDEX(name, "pwadvection_3d_local")
    if (strcomp .ne. 0) then
      field_information%field_type       = COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions= 3
      field_information%dimension_sizes(1) = current_state%local_grid%size(Z_INDEX)
      field_information%dimension_sizes(2) = current_state%local_grid%size(Y_INDEX)
      field_information%dimension_sizes(3) = current_state%local_grid%size(X_INDEX)
      field_information%data_type        = COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_u
      else if (name .eq. "tend_v_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_v
      else if (name .eq. "tend_w_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_w
      else if (name .eq. "tend_th_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_th
      else if (name .eq. "tend_qv_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_qv
      else if (name .eq. "tend_ql_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_ql
      else if (name .eq. "tend_qi_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_qi
      else if (name .eq. "tend_qr_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_qr
      else if (name .eq. "tend_qs_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_qs
      else if (name .eq. "tend_qg_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_qg
      else if (name .eq. "tend_tabs_pwadvection_3d_local") then
        field_information%enabled = l_tend_3d_tabs
      else
        field_information%enabled = .true.
      end if
    end if

    strcomp = INDEX(name, "pwadvection_profile_total_local")
    if (strcomp .ne. 0) then
      field_information%field_type       = COMPONENT_ARRAY_FIELD_TYPE
      field_information%number_dimensions= 1
      field_information%dimension_sizes(1) = current_state%local_grid%size(Z_INDEX)
      field_information%data_type        = COMPONENT_DOUBLE_DATA_TYPE

      if      (name .eq. "tend_u_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_u
      else if (name .eq. "tend_v_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_v
      else if (name .eq. "tend_w_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_w
      else if (name .eq. "tend_th_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_th
      else if (name .eq. "tend_qv_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_qv
      else if (name .eq. "tend_ql_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_ql
      else if (name .eq. "tend_qi_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_qi
      else if (name .eq. "tend_qr_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_qr
      else if (name .eq. "tend_qs_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_qs
      else if (name .eq. "tend_qg_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_qg
      else if (name .eq. "tend_tabs_pwadvection_profile_total_local") then
        field_information%enabled = l_tend_pr_tot_tabs
      else
        field_information%enabled = .true.
      end if
    end if
  end subroutine field_information_retrieval_callback

  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

    if      (name .eq. "tend_u_pwadvection_3d_local" .and. allocated(tend_3d_u)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_u)
    else if (name .eq. "tend_v_pwadvection_3d_local" .and. allocated(tend_3d_v)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_v)
    else if (name .eq. "tend_w_pwadvection_3d_local" .and. allocated(tend_3d_w)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_w)
    else if (name .eq. "tend_th_pwadvection_3d_local" .and. allocated(tend_3d_th)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_th)
    else if (name .eq. "tend_qv_pwadvection_3d_local" .and. allocated(tend_3d_qv)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qv)
    else if (name .eq. "tend_ql_pwadvection_3d_local" .and. allocated(tend_3d_ql)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_ql)
    else if (name .eq. "tend_qi_pwadvection_3d_local" .and. allocated(tend_3d_qi)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qi)
    else if (name .eq. "tend_qr_pwadvection_3d_local" .and. allocated(tend_3d_qr)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qr)
    else if (name .eq. "tend_qs_pwadvection_3d_local" .and. allocated(tend_3d_qs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qs)
    else if (name .eq. "tend_qg_pwadvection_3d_local" .and. allocated(tend_3d_qg)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_qg)
    else if (name .eq. "tend_tabs_pwadvection_3d_local" .and. allocated(tend_3d_tabs)) then
      call set_published_field_value(field_value, real_3d_field=tend_3d_tabs)
    else if (name .eq. "tend_u_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_u)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_u)
    else if (name .eq. "tend_v_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_v)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_v)
    else if (name .eq. "tend_w_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_w)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_w)
    else if (name .eq. "tend_th_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_th)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_th)
    else if (name .eq. "tend_qv_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qv)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qv)
    else if (name .eq. "tend_ql_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_ql)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_ql)
    else if (name .eq. "tend_qi_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qi)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qi)
    else if (name .eq. "tend_qr_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qr)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qr)
    else if (name .eq. "tend_qs_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qs)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qs)
    else if (name .eq. "tend_qg_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_qg)) then
      call set_published_field_value(field_value, real_1d_field=tend_pr_tot_qg)
    else if (name .eq. "tend_tabs_pwadvection_profile_total_local" .and. allocated(tend_pr_tot_tabs)) then
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

  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "pwompt" .or. trim(field) .eq. "any") then
        determine_if_advection_here = .true.
      else
        determine_if_advection_here = .false.
      end if
    else
      determine_if_advection_here = .true.
    end if
  end function determine_if_advection_here

end module pwadvectionompt_mod
