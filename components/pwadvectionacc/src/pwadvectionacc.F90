module pwadvectionacc_mod
  use datadefn_mod, only : DEFAULT_PRECISION
  use grids_mod, only : X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log
  use monc_component_mod, only : COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, component_descriptor_type, component_field_information_type, component_field_value_type
  use optionsdatabase_mod, only : options_get_string, options_get_integer
  use state_mod, only : model_state_type
  implicit none

  private
  logical :: advect_flow, advect_th, advect_q

  integer :: work_scaling

#ifdef USE_SINGLE_PRECISION
  real, dimension(:, :, :), allocatable, target :: su, sv, sw, sth
  real, dimension(:, :, :, :), allocatable, target :: sq
  real, dimension(:, :, :, :), allocatable :: q
  real, dimension(:, :, :), allocatable :: u, v, w, th
  real(kind=DEFAULT_PRECISION), dimension(:, :, :), allocatable, target :: dsu, dsv, dsw, dsth
  real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), allocatable, target :: dsq
#else
  real(kind=DEFAULT_PRECISION), dimension(:, :, :), allocatable, target :: su, sv, sw, sth
  real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), allocatable, target :: sq

  ! We pack the number_q_fields rank-3 fields into one rank-4 array for the accelerator
  real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), allocatable :: q

  ! These need updated each timestep, but are needed for all fields.
  real(kind=DEFAULT_PRECISION), dimension(:, :, :), allocatable :: u, v, w

#endif

  ! These are constant for the entire simulation.
  integer :: number_q_fields ! q
  integer, dimension(3) :: start_index, end_index ! flow, th, q
  real(kind=DEFAULT_PRECISION) :: tcx, tcy        ! flow
  real(kind=DEFAULT_PRECISION) :: cx, cy          ! th, q
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tzc1, tzc2 ! flow, th, q
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tzd1, tzd2 ! flow

  public pwadvectionacc_get_descriptor
contains
  type(component_descriptor_type) function pwadvectionacc_get_descriptor()
    pwadvectionacc_get_descriptor%name="pw_advection_acc"
    pwadvectionacc_get_descriptor%version=0.1
    pwadvectionacc_get_descriptor%initialisation=>initialisation_callback
    pwadvectionacc_get_descriptor%timestep=>timestep_callback
    pwadvectionacc_get_descriptor%finalisation=>finalisation_callback
    pwadvectionacc_get_descriptor%field_value_retrieval=>field_value_retrieval_callback
    pwadvectionacc_get_descriptor%field_information_retrieval=>field_information_retrieval_callback
    allocate(pwadvectionacc_get_descriptor%published_fields(5))
    pwadvectionacc_get_descriptor%published_fields(1)="pw_advection_acc_su"
    pwadvectionacc_get_descriptor%published_fields(2)="pw_advection_acc_sv"
    pwadvectionacc_get_descriptor%published_fields(3)="pw_advection_acc_sw"
    pwadvectionacc_get_descriptor%published_fields(4)="pw_advection_acc_sth"
    pwadvectionacc_get_descriptor%published_fields(5)="pw_advection_acc_sq"
  end function pwadvectionacc_get_descriptor

  subroutine convert_double_to_real_3d(source_a, target_a)
    real(kind=DEFAULT_PRECISION), dimension(:, :, :), intent(inout) :: source_a
    real, dimension(:, :, :), intent(inout) :: target_a

    integer :: i, j, k

    do i=1, size(source_a, 3)
      do j=1, size(source_a, 2)
        do k=1, size(source_a, 1)
          target_a(k, j, i)=real(source_a(k, j, i), kind=4)
        end do
      end do
    end do
  end subroutine convert_double_to_real_3d

  subroutine convert_double_to_real_4d(source_a, target_a)
    real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), intent(inout) :: source_a
    real, dimension(:, :, :, :), intent(inout) :: target_a

    integer :: i, j, k, n

    do n=1, size(source_a, 4)
    do i=1, size(source_a, 3)
      do j=1, size(source_a, 2)
        do k=1, size(source_a, 1)
          target_a(k, j, i, n)=real(source_a(k, j, i, n), kind=4)
        end do
      end do
    end do
    end do
  end subroutine convert_double_to_real_4d

  subroutine convert_real_to_double_3d(source_a, target_a)
    real, dimension(:, :, :), intent(inout) :: source_a
    real(kind=DEFAULT_PRECISION), dimension(:, :, :), intent(inout) :: target_a

    integer :: i, j, k

    do i=1, size(source_a, 3)
      do j=1, size(source_a, 2)
        do k=1, size(source_a, 1)
          target_a(k, j, i)=real(source_a(k, j, i), kind=DEFAULT_PRECISION)
        end do
      end do
    end do
  end subroutine convert_real_to_double_3d

  subroutine convert_real_to_double_4d(source_a, target_a)
    real, dimension(:, :, :, :), intent(inout) :: source_a
    real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), intent(inout) :: target_a
    
    integer :: i, j, k, n
    
    do n=1, size(source_a, 4)
    do i=1, size(source_a, 3)
      do j=1, size(source_a, 2)
        do k=1, size(source_a, 1)
          target_a(k, j, i, n)=real(source_a(k, j, i, n), kind=DEFAULT_PRECISION)
        end do
      end do
    end do
    end do
  end subroutine convert_real_to_double_4d

  subroutine field_information_retrieval_callback(current_state, name, field_sizing)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_information_type), intent(out) :: field_sizing

    field_sizing%field_type = COMPONENT_ARRAY_FIELD_TYPE
    field_sizing%data_type  = COMPONENT_DOUBLE_DATA_TYPE

    if (name .eq. "pw_advection_acc_su") then
      !$acc wait(31)
      !$acc update host(su)
      field_sizing%number_dimensions = 3
      field_sizing%dimension_sizes   = [size(su, 1), size(su, 2), size(su, 3), 0]
      field_sizing%enabled           = advect_flow
    else if (name .eq. "pw_advection_acc_sv") then
      !$acc wait(31)
      !$acc update host(sv)
      field_sizing%number_dimensions = 3
      field_sizing%dimension_sizes   = [size(sv, 1), size(sv, 2), size(sv, 3), 0]
      field_sizing%enabled           = advect_flow
    else if (name .eq. "pw_advection_acc_sw") then
      !$acc wait(31)
      !$acc update host(sw)
      field_sizing%number_dimensions = 3
      field_sizing%dimension_sizes   = [size(sw, 1), size(sw, 2), size(sw, 3), 0]
      field_sizing%enabled           = advect_flow
    else if (name .eq. "pw_advection_acc_sth") then
      !$acc wait(21)
      !$acc update host(sth)
      field_sizing%number_dimensions = 3
      field_sizing%dimension_sizes   = [size(sth, 1), size(sth, 2), size(sth, 3), 0]
      field_sizing%enabled           = advect_th
    else if (name .eq. "pw_advection_acc_sq") then
      !$acc wait(11)
      !$acc update host(sq)
      field_sizing%number_dimensions = 4
      field_sizing%dimension_sizes   = [size(sq, 1), size(sq, 2), size(sq, 3), size(sq, 4)]
      field_sizing%enabled           = advect_q
    else
      call log_master_log(LOG_ERROR, "Mismatch between published_fields and field_information_retrieval_callback in pwadvectionacc_mod.")
    end if
  end subroutine field_information_retrieval_callback

  subroutine field_value_retrieval_callback(current_state, name, field_value)
    type(model_state_type), target, intent(inout) :: current_state
    character(len=*), intent(in) :: name
    type(component_field_value_type), intent(out) :: field_value

#ifdef USE_SINGLE_PRECISION
   if (name .eq. "pw_advection_acc_su") then
        call convert_real_to_double_3d(su, dsu)
        field_value%real_3d_array_ptr => dsu
   else if (name .eq. "pw_advection_acc_sv") then
        call convert_real_to_double_3d(sv, dsv)
        field_value%real_3d_array_ptr => dsv
   else if (name .eq. "pw_advection_acc_sw") then
        call convert_real_to_double_3d(sw, dsw)
        field_value%real_3d_array_ptr => dsw
   else if (name .eq. "pw_advection_acc_sth") then
        call convert_real_to_double_3d(sth, dsth)
        field_value%real_3d_array_ptr => dsth
   else if (name .eq. "pw_advection_acc_sq") then
        call convert_real_to_double_4d(sq, dsq)
        field_value%real_4d_array_ptr => dsq
   else
        call log_master_log(LOG_ERROR, "Mismatch between published_fields and field_value_retrieval_callback in pwadvectionacc_mod.")
   end if
#else
   if (name .eq. "pw_advection_acc_su") then
      field_value%real_3d_array_ptr => su
   else if (name .eq. "pw_advection_acc_sv") then
      field_value%real_3d_array_ptr => sv
   else if (name .eq. "pw_advection_acc_sw") then
      field_value%real_3d_array_ptr => sw
   else if (name .eq. "pw_advection_acc_sth") then
      field_value%real_3d_array_ptr => sth
   else if (name .eq. "pw_advection_acc_sq") then
      field_value%real_4d_array_ptr => sq
   else
      call log_master_log(LOG_ERROR, "Mismatch between published_fields and field_value_retrieval_callback in pwadvectionacc_mod.")
    end if
#endif
  end subroutine field_value_retrieval_callback      

  subroutine initialisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    advect_flow = determine_if_advection_here(options_get_string(current_state%options_database, "advection_flow_fields"))
    advect_th   = determine_if_advection_here(options_get_string(current_state%options_database, "advection_theta_field"))
    advect_q    = determine_if_advection_here(options_get_string(current_state%options_database, "advection_q_fields"))

#ifdef SCALING_WORKLOAD
    work_scaling=options_get_integer(current_state%options_database, "work_scaling")
#endif

    advect_th = advect_th .and. current_state%th%active
    advect_q = advect_q .and. number_q_fields .gt. 0

    if (advect_flow .or. advect_th .or. advect_q) then
      allocate(u(size(current_state%u%data, 1), size(current_state%u%data, 2), size(current_state%u%data, 3)))
      allocate(v(size(current_state%v%data, 1), size(current_state%v%data, 2), size(current_state%v%data, 3)))
      allocate(w(size(current_state%w%data, 1), size(current_state%w%data, 2), size(current_state%w%data, 3)))

      ! TODO: A few of these are only required if e.g. advect_flow, rather than
      !       for all fields. These could be checked and copied individually,
      !       but we would need to take care to call 'enter' and 'exit data' the
      !       right number of times.
      number_q_fields = current_state%number_q_fields

      start_index = current_state%local_grid%local_domain_start_index
      end_index   = current_state%local_grid%local_domain_end_index

      tcx = current_state%global_grid%configuration%horizontal%tcx
      tcy = current_state%global_grid%configuration%horizontal%tcy

      cx = current_state%global_grid%configuration%horizontal%cx
      cy = current_state%global_grid%configuration%horizontal%cy

      allocate(tzc1, source=current_state%global_grid%configuration%vertical%tzc1)
      allocate(tzc2, source=current_state%global_grid%configuration%vertical%tzc2)

      allocate(tzd1, source=current_state%global_grid%configuration%vertical%tzd1)
      allocate(tzd2, source=current_state%global_grid%configuration%vertical%tzd2)

      !$acc enter data create(u, v, w) copyin(number_q_fields, work_scaling, start_index, end_index, tcx, tcy, cx, cy, tzc1, tzc2, tzd1, tzd2)
    end if

#ifdef USE_SINGLE_PRECISION
 call log_master_log(LOG_INFO, "Working with single precision")
#endif

    if (advect_flow) then
      call log_master_log(LOG_INFO, "Accelerating advection: flow fields")

      allocate(su(size(current_state%su%data, 1), size(current_state%su%data, 2), size(current_state%su%data, 3)))
      allocate(sv(size(current_state%sv%data, 1), size(current_state%sv%data, 2), size(current_state%sv%data, 3)))
      allocate(sw(size(current_state%sw%data, 1), size(current_state%sw%data, 2), size(current_state%sw%data, 3)))

#ifdef USE_SINGLE_PRECISION
      allocate(dsu(size(current_state%su%data, 1), size(current_state%su%data, 2), size(current_state%su%data, 3)))
      allocate(dsv(size(current_state%sv%data, 1), size(current_state%sv%data, 2), size(current_state%sv%data, 3)))
      allocate(dsw(size(current_state%sw%data, 1), size(current_state%sw%data, 2), size(current_state%sw%data, 3)))
#endif

      su = 0
      sv = 0
      sw = 0

      !$acc enter data copyin(su, sv, sw)
    end if

    if (advect_th) then
      call log_master_log(LOG_INFO, "Accelerating advection: theta field")

      allocate(sth(size(current_state%sth%data, 1), size(current_state%sth%data, 2), size(current_state%sth%data, 3)))
#ifdef USE_SINGLE_PRECISION
      allocate(dsth(size(current_state%sth%data, 1), size(current_state%sth%data, 2), size(current_state%sth%data, 3)))
      allocate(th(size(current_state%sth%data, 1), size(current_state%sth%data, 2), size(current_state%sth%data, 3)))
#endif

      sth = 0

      !$acc enter data copyin(sth)
    end if

    if (advect_q) then
      call log_master_log(LOG_INFO, "Accelerating advection: q fields")

      allocate(q(size(current_state%q(1)%data, 1), size(current_state%q(1)%data, 2), size(current_state%q(1)%data, 3), number_q_fields))
      allocate(sq(size(current_state%sq(1)%data, 1), size(current_state%sq(1)%data, 2), size(current_state%sq(1)%data, 3), number_q_fields))

#ifdef USE_SINGLE_PRECISION
      allocate(dsq(size(current_state%sq(1)%data, 1), size(current_state%sq(1)%data, 2), size(current_state%sq(1)%data, 3), number_q_fields))
#endif

      sq = 0

      !$acc enter data create(q) copyin(sq)
    end if
  end subroutine initialisation_callback

  subroutine finalisation_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    if (advect_q) then
      ! q, sq
      !!!$acc exit data

      deallocate(q)
      deallocate(sq)
    end if

    if (advect_th) then
      ! sth
      !!!$acc exit data

      deallocate(sth)
    end if

    if (advect_flow) then
      ! su, sv, sw
      !!!$acc exit data

      deallocate(su)
      deallocate(sv)
      deallocate(sw)
    end if

    if (advect_flow .or. advect_th .or. advect_q) then
      ! u, v, w, start_index, end_index, tcx, tcy, cx, cy, tzc1, tzc2, tzd1, tzd2
      !!!$acc exit data

      deallocate(u)
      deallocate(v)
      deallocate(w)

      deallocate(tzc1)
      deallocate(tzc2)
      deallocate(tzd1)
      deallocate(tzd2)
    end if
  end subroutine finalisation_callback

  subroutine timestep_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state
    integer :: n

    if (advect_flow .or. advect_th .or. advect_q) then
#ifdef USE_SINGLE_PRECISION
     call convert_double_to_real_3d(current_state%u%data, u)
     call convert_double_to_real_3d(current_state%v%data, v)
     call convert_double_to_real_3d(current_state%w%data, w)
#else
      u = current_state%u%data
      v = current_state%v%data
      w = current_state%w%data
#endif

      !$acc update device(u, v, w) async(1)
    end if

    if (advect_flow) then
      ! stepfields sums our published field with other source terms, but we must
      ! zero the current_state arrays as there is no reset between steps other
      ! than that provided by the advection component i.e. if we don't do this,
      ! we add on top of the previous steps' source terms.
      current_state%su%data = 0
      current_state%sv%data = 0
      current_state%sw%data = 0

      call advect_flow_fields_accelerated()
    end if
    
    if (advect_th) then
      ! There doesn't appear to be any dependency on the value of th (or sq)
      ! between here and stepfields, but zero the arrays just to be safe.
      !current_state%th%data = 0
#ifdef USE_SINGLE_PRECISION
      call convert_double_to_real_3d(current_state%th%data, th)
      call advect_th_field_accelerated(th)
#else
      call advect_th_field_accelerated(current_state%th%data)
#endif
    end if

    if (advect_q) then
      ! advect_q_fields needs to iterate at a higher level than the above two i.e. current_state%q(n)%data(k, y, x) vs. current_state%q%data(k, y, x)
      do n = 1, number_q_fields
        current_state%q(n)%data = 0
      end do

      call advect_q_fields_accelerated(current_state)
    end if
  end subroutine timestep_callback

  ! TODO: Verify result!
  subroutine advect_q_fields_accelerated(current_state)
    type(model_state_type), intent(inout) :: current_state
    integer :: n

    do n = 1, current_state%number_q_fields
#ifdef USE_SINGLE_PRECISION
      call convert_double_to_real_3d(current_state%q(n)%data, q(:, :, :, n))
#else
      q(:, :, :, n) = current_state%q(n)%data
#endif
    end do

    !$acc update device(q) async(2)
    call advect_q_field_accelerated(q)
  end subroutine advect_q_fields_accelerated

  ! TODO: Verify result!
  subroutine advect_q_field_accelerated(q)
#ifdef USE_SINGLE_PRECISION
    real, dimension(:, :, :, :), intent(in) :: q
#else
    real(kind=DEFAULT_PRECISION), dimension(:, :, :, :), intent(in) :: q
#endif
    integer :: n, x, y, k, r

#ifdef SCALING_WORKLOAD
   !$acc parallel loop collapse(5) wait(1, 2) async(10) default(none) present(q, sq, u, v, w, number_q_fields, start_index, end_index, cx, cy, tzc1, tzc2, work_scaling)
   do r=1, work_scaling
#else
    !$acc parallel loop collapse(4) wait(1, 2) async(10) default(none) present(q, sq, u, v, w, number_q_fields, start_index, end_index, cx, cy, tzc1, tzc2)
#endif
    do n = 1, number_q_fields
      do x = start_index(X_INDEX), end_index(X_INDEX)
        do y = start_index(Y_INDEX), end_index(Y_INDEX)
          do k = start_index(Z_INDEX), end_index(Z_INDEX)
            sq(k, y, x, n) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
            sq(k, y, x, n) = 0.5_DEFAULT_PRECISION&
                             * cx * (u(k, y, x-1) * q(k, y, x-1, n)&
                                     - u(k, y, x) * q(k, y, x+1, n))
#endif
#ifdef V_ACTIVE
            sq(k, y, x, n) = sq(k, y, x, n) +&
                             0.5_DEFAULT_PRECISION&
                             * cy * (v(k, y-1, x) * q(k, y-1, x, n)&
                                     - v(k, y, x) * q(k, y+1, x, n))
#endif
#ifdef W_ACTIVE
            sq(k, y, x, n) = sq(k, y, x, n) +&
                             2.0_DEFAULT_PRECISION&
                             * (tzc1(k) * w(k-1, y, x) * q(k-1, y, x, n)&
                                - tzc2(k) * w(k, y, x) * q(k+1, y, x, n))
#endif
          end do
        end do
      end do
    end do

#ifdef SCALING_WORKLOAD
    end do
#endif

    !$acc update host(sq) wait(10) async(11)

    ! TODO: if (l_toplevel) ...
  end subroutine advect_q_field_accelerated

  ! TODO: Verify result!
  subroutine advect_th_field_accelerated(th)
#ifdef USE_SINGLE_PRECISION
    real, dimension(:, :, :), intent(in) :: th
#else
    real(kind=DEFAULT_PRECISION), dimension(:, :, :), intent(in) :: th
#endif
    integer :: x, y, k, n

#ifdef SCALING_WORKLOAD
    !$acc parallel loop collapse(4) wait(1) async(20) default(none) copyin(th) present(sth, u, v, w, start_index, end_index, cx, cy, tzc1, tzc2, work_scaling)
    do n=1, work_scaling
#else
    !$acc parallel loop collapse(3) wait(1) async(20) default(none) copyin(th) present(sth, u, v, w, start_index, end_index, cx, cy, tzc1, tzc2)
#endif
    do x = start_index(X_INDEX), end_index(X_INDEX)
      do y = start_index(Y_INDEX), end_index(Y_INDEX)
        do k = start_index(Z_INDEX), end_index(Z_INDEX)
          sth(k, y, x) = 0.0_DEFAULT_PRECISION
#ifdef U_ACTIVE
          sth(k, y, x) = 0.5_DEFAULT_PRECISION&
                         * cx * (u(k, y, x-1) * th(k, y, x-1)&
                                 - u(k, y, x) * th(k, y, x+1))
#endif
#ifdef V_ACTIVE
          sth(k, y, x) = sth(k, y, x) +&
                         0.5_DEFAULT_PRECISION&
                         * cy * (v(k, y-1, x) * th(k, y-1, x)&
                                 - v(k, y, x) * th(k, y+1, x))
#endif
#ifdef W_ACTIVE
          sth(k, y, x) = sth(k, y, x) +&
                         2.0_DEFAULT_PRECISION&
                         * (tzc1(k) * w(k-1, y, x) * th(k-1, y, x)&
                            - tzc2(k) * w(k, y, x) * th(k+1, y, x))
#endif
        end do
      end do
    end do
#ifdef SCALING_WORKLOAD
    end do
#endif

    !$acc update host(sth) wait(20) async(21)

    ! TODO: if (l_toplevel) ...
  end subroutine advect_th_field_accelerated

  subroutine advect_flow_fields_accelerated()
    integer :: x, y, k, n

#ifdef SCALING_WORKLOAD
    !$acc parallel loop collapse(4) wait(1) async(30) default(none) present(u, v, w, su, sv, sw, start_index, end_index, tcx, tcy, tzc1, tzc2, tzd1, tzd2, work_scaling)
    do n=1, work_scaling
#else
    !$acc parallel loop collapse(3) wait(1) async(30) default(none) present(u, v, w, su, sv, sw, start_index, end_index, tcx, tcy, tzc1, tzc2, tzd1, tzd2)
#endif
    do x = start_index(X_INDEX), end_index(X_INDEX)
      do y = start_index(Y_INDEX), end_index(Y_INDEX)
        do k = start_index(Z_INDEX), end_index(Z_INDEX)
#ifdef U_ACTIVE
          su(k, y, x) = tcx * (u(k, y, x-1)&
                               * (u(k, y, x) + u(k, y, x-1))&
                               - u(k, y, x+1)&
                                 * (u(k, y, x) + u(k, y, x+1)))
#ifdef V_ACTIVE
          su(k, y, x) = su(k, y, x) +&
                        tcy * (u(k, y-1, x)&
                               * (v(k, y-1, x) + v(k, y-1, x+1))&
                               - u(k, y+1, x)&
                                 * (v(k, y, x) + v(k, y, x+1)))
#endif
#ifdef W_ACTIVE
          su(k, y, x) = su(k, y, x) +&
                        tzc1(k) * u(k-1, y, x)&
                                * (w(k-1, y, x) + w(k-1, y, x+1))&
                                - tzc2(k) * u(k+1, y, x)&
                                          * (w(k, y, x) + w(k, y, x+1))
#endif
#endif

#ifdef V_ACTIVE
          sv(k, y, x) = tcy * (v(k, y-1, x)&
                               * (v(k, y, x) + v(k, y-1, x))&
                               - v(k, y+1, x) &
                                 * (v(k, y, x) + v(k, y+1, x)))
#ifdef U_ACTIVE
          sv(k, y, x) = sv(k, y, x) +&
                        tcx * (v(k, y, x-1)&
                               * (u(k, y, x-1) + u(k, y+1, x-1))&
                               - v(k, y, x+1)&
                                 * (u(k, y, x) + u(k, y+1, x)))
#endif
#ifdef W_ACTIVE
          sv(k, y, x) = sv(k, y, x) +&
                        tzc1(k) * v(k-1, y, x)&
                                * (w(k-1, y, x) + w(k-1, y+1, x))&
                                - tzc2(k) * v(k+1, y, x)&
                                          * (w(k, y, x) + w(k, y+1, x))
#endif
#endif

#ifdef W_ACTIVE
          sw(k, y, x) = tzd1(k) * w(k-1, y, x)&
                                * (w(k, y, x) + w(k-1, y, x))&
                                - tzd2(k) * w(k+1, y, x)&
                                          * (w(k, y, x) + w(k+1, y, x))
#ifdef U_ACTIVE
          sw(k, y, x) = sw(k, y, x) +&
                        tcx * (w(k, y, x-1) &
                               * (u(k, y, x-1) + u(k+1, y, x-1))&
                               - w(k, y, x+1)&
                                 * (u(k, y, x) + u(k+1, y, x)))
#endif
#ifdef V_ACTIVE
          sw(k, y, x) = sw(k, y, x) +&
                        tcy * (w(k, y-1, x)&
                               * (v(k, y-1, x) + v(k+1, y-1, x))&
                               - w(k, y+1, x)&
                                 * (v(k, y, x) + v(k+1, y, x)))
#endif
#endif
        end do
      end do
    end do

#ifdef SCALING_WORKLOAD
    end do
#endif

    !$acc update host(su, sv, sw) wait(30) async(31)

    ! TODO: if (l_toplevel) ...
  end subroutine advect_flow_fields_accelerated

  !> Parses a field string (read in from the configuration file) and determines whether this algorithm should be used
  !! for advecting that field
  !! @param field The string configuration of field advection
  !! @returns Whether or not the field is advected here
  logical function determine_if_advection_here(field)
    character(len=*), intent(in) :: field

    if (len_trim(field) .ne. 0) then
      if (trim(field) .eq. "pwacc" .or. trim(field) .eq. "any") then
        determine_if_advection_here=.true.
      else
        determine_if_advection_here=.false.
      end if
    else
      determine_if_advection_here=.true.
    end if
  end function determine_if_advection_here  
end module pwadvectionacc_mod
