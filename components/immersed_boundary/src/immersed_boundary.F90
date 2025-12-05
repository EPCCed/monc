! Calculates immersed boundary node values to be applied in ib_finalise.
! Includes IB versions of MOSTBC and u* lookup table subroutines
! Used ghost point (GP) or Velocity reconstruction (VR) methods (e.g. DeLeon et al. Boundary Layer
! Meteorology 2018, 167)

module immersed_boundary_mod

  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type

  use state_mod, only : FORWARD_STEPPING, PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, &
   model_state_type
  use grids_mod
  use prognostics_mod, only : prognostic_field_type
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_array_size, &
     options_get_real_array, options_get_string, options_get_logical
  use science_constants_mod
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_strerror, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_master_log, log_log, log_get_logging_level
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, PRECISION_TYPE
  use conversions_mod, only : conv_to_string
  use q_indices_mod, only: get_q_index, standard_q_names

  implicit none

#ifndef TEST_MODE
  private
#endif

  real(kind=DEFAULT_PRECISION) :: dump_freq
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: lograt_u, lograt_v, lograt_w ! VR coefficients
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: hrat_u, hrat_v, hrat_w, hrat_s, hrat_diff_s ! VR coefficients
  integer :: dump_number, didi, didj, nghost_u, nghost_v, nghost_w, nghost_s
  logical :: enable_theta=.false.
  ! arrays for MOST solver on IB points
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: zlogth_ib, zlogm_ib, tstrcona_ib, rhmbc_ib
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ddbc_ib, cmbc_ib, eecon_ib, delz_ib
  real(kind=DEFAULT_PRECISION) :: xx0con, yy0con, tstrconb
  integer, parameter :: CONVERGENCE_SUCCESS=1, CONVERGENCE_RICHARDSON_TOO_LARGE=2, CONVERGENCE_FAILURE=3
  real(kind=DEFAULT_PRECISION),parameter :: prandtl=0.7 ! neutral turbulent Prandtl number from Mason & Brown 1999
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: w2s_store_diff ! array for storing diffusivity on w BI points 
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: fbuoy_ib
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: lookup_ib_ustr, lookup_ib_vel
  integer(kind=DEFAULT_PRECISION), dimension(:), allocatable :: lookup_map_g2t, lookup_map_t2g
  real(kind=DEFAULT_PRECISION) :: velmax, velmin, aloginv
  integer(kind=DEFAULT_PRECISION) :: nunst, iqv
 
  public immersed_boundary_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function immersed_boundary_get_descriptor()
    immersed_boundary_get_descriptor%name="immersed_boundary"
    immersed_boundary_get_descriptor%version=0.1
    immersed_boundary_get_descriptor%timestep=>timestep_callback
    immersed_boundary_get_descriptor%initialisation=>initialisation_callback
    immersed_boundary_get_descriptor%finalisation=>finalisation_callback

  end function immersed_boundary_get_descriptor
  


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    ! IB interpolation arrays for scalar surface points
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid
    real(kind=DEFAULT_PRECISION) :: tmp1, tmp2, tmp3, zzd, delz
    real(kind=DEFAULT_PRECISION) :: bhbc
    integer :: ighost, k,kk, myproc
    character (len=50):: fname

    enable_theta=options_get_logical(current_state%options_database, "enable_theta")
    myproc = current_state%parallel%my_rank

    nghost_u = size(current_state%immersed%gp_ijk_u,1)
    nghost_v = size(current_state%immersed%gp_ijk_v,1)
    nghost_w = size(current_state%immersed%gp_ijk_w,1)
    nghost_s = size(current_state%immersed%gp_ijk_s,1)
    

    if(enable_theta)then
      current_state%immersed%theta_surf_ib_s(:)=current_state%theta_virtual_surf
      current_state%immersed%theta_surf_ib_w(:)=current_state%theta_virtual_surf
    end if

    
    if(current_state%immersed%ib_type==1) then ! allocate VR arrays
      allocate(lograt_u(nghost_u))
      allocate(lograt_v(nghost_v))
      allocate(lograt_w(nghost_w))
      allocate(hrat_u(nghost_u))
      allocate(hrat_v(nghost_v))
      allocate(hrat_w(nghost_w))
      allocate(hrat_s(nghost_s))
      if(current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES)then
        allocate(hrat_diff_s(nghost_s))
      end if


      do ighost =1,nghost_u
        lograt_u(ighost) = log((current_state%immersed%gp_len_u(ighost)+z0)/z0)/&
                           log((current_state%immersed%ip_len_u(ighost)+z0)/z0)
        hrat_u(ighost) = current_state%immersed%gp_len_u(ighost)/&
                         current_state%immersed%ip_len_u(ighost)
      end do
      do ighost =1,nghost_v
        lograt_v(ighost) = log((current_state%immersed%gp_len_v(ighost)+z0)/z0)/&
                           log((current_state%immersed%ip_len_v(ighost)+z0)/z0)
        hrat_v(ighost) = current_state%immersed%gp_len_v(ighost)/&
                         current_state%immersed%ip_len_v(ighost)
      end do
      do ighost =1,nghost_w
        lograt_w(ighost) = log((current_state%immersed%gp_len_w(ighost)+z0)/z0)/&
                           log((current_state%immersed%ip_len_w(ighost)+z0)/z0)
        hrat_w(ighost) = current_state%immersed%gp_len_w(ighost)/&
                         current_state%immersed%ip_len_w(ighost)
      end do
      do ighost =1,nghost_s
        hrat_s(ighost) = current_state%immersed%gp_len_s(ighost)/&
                         current_state%immersed%ip_len_s(ighost)
      end do
      if(current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES)then
        do ighost =1,nghost_s
          hrat_diff_s(ighost) = (current_state%immersed%gp_len_s(ighost)+&
                                 current_state%immersed%ip_len_s(ighost))/2.0_DEFAULT_PRECISION
        end do
      end if

    end if ! VR only


    if(current_state%use_surface_boundary_conditions)then
      ! set MOST variables for IB points
      allocate(zlogth_ib(nghost_w))
      allocate(zlogm_ib(nghost_w))
      allocate(delz_ib(nghost_w))
      do ighost=1,nghost_w
        if(current_state%immersed%ib_type==0)delz_ib(ighost)=current_state%immersed%k2e_len(ighost)
        if(current_state%immersed%ib_type==1)delz_ib(ighost)=current_state%immersed%ip_len_w(ighost)
        zlogm_ib(ighost)=log(1.0_DEFAULT_PRECISION+delz_ib(ighost)/z0)
        zlogth_ib(ighost)=log((delz_ib(ighost)+z0)/z0th)
      end do
      if(current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
        allocate(tstrcona_ib(nghost_w))
        allocate(rhmbc_ib(nghost_w))
        allocate(ddbc_ib(nghost_w))
        allocate(cmbc_ib(nghost_w))
        allocate(eecon_ib(nghost_w))
        do ighost=1,nghost_w
          tstrcona_ib(ighost)=von_karman_constant/alphah*zlogth_ib(ighost)
          bhbc=alphah*zlogth_ib(ighost)
          rhmbc_ib(ighost)=betah*(delz_ib(ighost)+z0-z0th)/&
                           (betam*delz_ib(ighost))
          ddbc_ib(ighost)=zlogm_ib(ighost)*(bhbc-&
              rhmbc_ib(ighost)*zlogm_ib(ighost))
          eecon_ib(ighost)=2.0_DEFAULT_PRECISION*rhmbc_ib(ighost)*zlogm_ib(ighost)-bhbc
          cmbc_ib(ighost)=betam*delz_ib(ighost)*G*von_karman_constant/&
                          current_state%immersed%theta_surf_ib_w(ighost)
          xx0con=gammam*z0
          yy0con=gammah*z0th
          tstrconb=von_karman_constant/alphah
        end do
      end if
    end if

    if ( current_state%use_surface_boundary_conditions .and.&
         current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
      allocate(w2s_store_diff(nghost_w))
      allocate(fbuoy_ib(nghost_w))

      ! initialise to something more sensible !
      w2s_store_diff=0.0_DEFAULT_PRECISION

      if(current_state%number_q_fields > 0)then
        iqv = current_state%water_vapour_mixing_ratio_index
      end if

      fbuoy_ib=0.
      nunst = 0 ! count number of unstable points
      if(.not. current_state%passive_th)then
        zgrid = current_state%global_grid%configuration%vertical%z(:)
        do ighost=1,nghost_w   
          tmp3 = current_state%immersed%bi_xyz_w(ighost,3)
          do kk=1, size(zgrid)-1
            if (tmp3 .ge. zgrid(kk) .and. tmp3 .lt. zgrid(kk+1))k = kk
          end do
          tmp1 = current_state%global_grid%configuration%vertical%buoy_co(k) 
          tmp2 = current_state%global_grid%configuration%vertical%buoy_co(k+1)
          zzd = (tmp3-zgrid(k))/(zgrid(k+1)-zgrid(k))
          fbuoy_ib(ighost)=(((1.-zzd)*tmp1)+(zzd*tmp2))*current_state%surface_temperature_flux
          if(.not. current_state%passive_q .and. current_state%number_q_fields > 0)then
            fbuoy_ib(ighost)=fbuoy_ib(ighost)+current_state%cq(iqv)*current_state%surface_vapour_flux*G
          end if ! q active
          if(fbuoy_ib(ighost).gt.0.0_DEFAULT_PRECISION)nunst=nunst+1
        end do ! ighost

        if(nunst.gt.0)then
          allocate(lookup_ib_ustr(nunst,current_state%lookup_table_entries))
          allocate(lookup_ib_vel(nunst,current_state%lookup_table_entries))
          allocate(lookup_map_g2t(nghost_w))
          allocate(lookup_map_t2g(nunst))
          nunst=0
          lookup_map_g2t=-1
          lookup_map_t2g=-1
          do ighost=1,nghost_w
            if(fbuoy_ib(ighost).gt.0.0_DEFAULT_PRECISION)then
              nunst=nunst+1
              lookup_map_g2t(ighost)=nunst
              lookup_map_t2g(nunst)=ighost
            end if
          end do
          call set_look_ib(current_state)


        end if


      end if ! theta active
    end if ! use BCs & fixed fluxes

    current_state%immersed%gp_store_u = 0.0_DEFAULT_PRECISION
    current_state%immersed%gp_store_v = 0.0_DEFAULT_PRECISION
    current_state%immersed%gp_store_w = 0.0_DEFAULT_PRECISION
    if(enable_theta)then
      current_state%immersed%gp_store_th = 0.0_DEFAULT_PRECISION
    end if
    if (current_state%use_viscosity_and_diffusion)then
      current_state%immersed%gp_store_vis = 0.0_DEFAULT_PRECISION
      current_state%immersed%gp_store_diff = 0.0_DEFAULT_PRECISION
    end if
    if(current_state%number_q_fields > 0)then
      current_state%immersed%gp_store_q = 0.0_DEFAULT_PRECISION
    end if


  end subroutine initialisation_callback




  !> Timestep callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: i,j,k, ighost, jghost, myproc, q
    integer ipi, ipj, ipk
    real(kind=DEFAULT_PRECISION) :: xxd, yyd, zzd, cc
    real(kind=DEFAULT_PRECISION) :: ut_ip, vt_ip, wt_ip, un_ip, vn_ip, wn_ip
    real(kind=DEFAULT_PRECISION) :: pvel_at_k2e, uvel, vvel, wvel
    real(kind=DEFAULT_PRECISION) :: nhat_x, nhat_y, nhat_z, that_x, that_y, that_z
    real(kind=DEFAULT_PRECISION) :: tvel_gp, nvel_gp, tvel_ip
    real(kind=DEFAULT_PRECISION) :: vis_ip, diff_ip, thref_ip, vistmp, difftmp
    real(kind=DEFAULT_PRECISION) :: tmp1, tmp2, tmp3 
    real(kind=DEFAULT_PRECISION) :: visneut, diffneut
    real(kind=DEFAULT_PRECISION) :: dthv_surf, ustr, thvstr, flag, ema_gamma
    integer :: convergence_status, n, x_size, y_size, dcheck, ierr
    logical :: docontrol

    myproc = current_state%parallel%my_rank

!    call mpi_barrier(current_state%parallel%monc_communicator, ierr)

    flag = -1

    if (current_state%immersed%ib_enabled) then 
      current_state%immersed%dodiags = .false.
      if ((current_state%time .gt. 0.0_DEFAULT_PRECISION) .and.&
           current_state%immersed%diags_enabled .and. &
          (abs(mod(current_state%time, current_state%immersed%dump_freq)) .lt. current_state%dtm)) then
        current_state%immersed%dodiags = .true.
      end if
    end if


    if (current_state%halo_column) return

    if (current_state%immersed%ib_enabled) then 
      j=current_state%column_local_y
      i=current_state%column_local_x
      if (current_state%immersed%ib_col(j,i)) then

  
        ! no-slip 
        if (current_state%immersed%ib_type==0)then
!=============================================================================
          do k=1,current_state%immersed%kmax_ji(j,i)
    
            
            if (current_state%immersed%ghost_u(k,j,i) .gt. 0) then ! u points
              ighost = current_state%immersed%ghost_u(k,j,i)
  
              ! diagnostics
              if(current_state%immersed%dump_u.and.current_state%immersed%dodiags)then
                ipi = current_state%immersed%bi_ijk_u(ighost,1)
                ipj = current_state%immersed%bi_ijk_u(ighost,2)
                ipk = current_state%immersed%bi_ijk_u(ighost,3)
                xxd = current_state%immersed%biint_xyz_u(ighost,1)
                yyd = current_state%immersed%biint_xyz_u(ighost,2)
                zzd = current_state%immersed%biint_xyz_u(ighost,3)
                call interp_field(current_state%zu, ipi, ipj, ipk, xxd, yyd, zzd, cc)
                current_state%immersed%ib_u(ighost) = cc
              end if
  
              ipi = current_state%immersed%ip_ijk_uu(ighost,1)
              ipj = current_state%immersed%ip_ijk_uu(ighost,2)
              ipk = current_state%immersed%ip_ijk_uu(ighost,3)
              xxd = current_state%immersed%int_xyz_uu(ighost,1)
              yyd = current_state%immersed%int_xyz_uu(ighost,2)
              zzd = current_state%immersed%int_xyz_uu(ighost,3)
              call interp_field(current_state%u, ipi, ipj, ipk, xxd, yyd, zzd, cc)
  
              current_state%immersed%gp_store_u(ighost) =  -cc 


            end if
  
  !=============================================================================
    
            if (current_state%immersed%ghost_v(k,j,i) .gt. 0) then ! v points
              ighost = current_state%immersed%ghost_v(k,j,i)
  
              if(current_state%immersed%dump_v.and.current_state%immersed%dodiags)then
                ! diagnostics 
                ipi = current_state%immersed%bi_ijk_v(ighost,1)
                ipj = current_state%immersed%bi_ijk_v(ighost,2)
                ipk = current_state%immersed%bi_ijk_v(ighost,3)
                xxd = current_state%immersed%biint_xyz_v(ighost,1)
                yyd = current_state%immersed%biint_xyz_v(ighost,2)
                zzd = current_state%immersed%biint_xyz_v(ighost,3)
                call interp_field(current_state%zv, ipi, ipj, ipk, xxd, yyd, zzd, cc)
                current_state%immersed%ib_v(ighost) = cc
              end if
  
  
              ipi = current_state%immersed%ip_ijk_vv(ighost,1)
              ipj = current_state%immersed%ip_ijk_vv(ighost,2)
              ipk = current_state%immersed%ip_ijk_vv(ighost,3)
              xxd = current_state%immersed%int_xyz_vv(ighost,1)
              yyd = current_state%immersed%int_xyz_vv(ighost,2)
              zzd = current_state%immersed%int_xyz_vv(ighost,3)
              call interp_field(current_state%v, ipi, ipj, ipk, xxd, yyd, zzd, cc)
  
              current_state%immersed%gp_store_v(ighost) =  -cc
  
            end if
  
    
  !=============================================================================
    
            if (current_state%immersed%ghost_w(k,j,i) .gt. 0) then ! w points
              ighost = current_state%immersed%ghost_w(k,j,i)
  
              ! diagnostics 
              if (current_state%immersed%dump_w.and.current_state%immersed%dodiags)then
                ipi = current_state%immersed%bi_ijk_w(ighost,1)
                ipj = current_state%immersed%bi_ijk_w(ighost,2)
                ipk = current_state%immersed%bi_ijk_w(ighost,3)
                xxd = current_state%immersed%biint_xyz_w(ighost,1)
                yyd = current_state%immersed%biint_xyz_w(ighost,2)
                zzd = current_state%immersed%biint_xyz_w(ighost,3)
                call interp_field(current_state%zw, ipi, ipj, ipk, xxd, yyd, zzd, cc)
                current_state%immersed%ib_w(ighost) = cc
              end if
  
              ipi = current_state%immersed%ip_ijk_ww(ighost,1)
              ipj = current_state%immersed%ip_ijk_ww(ighost,2)
              ipk = current_state%immersed%ip_ijk_ww(ighost,3)
              xxd = current_state%immersed%int_xyz_ww(ighost,1)
              yyd = current_state%immersed%int_xyz_ww(ighost,2)
              zzd = current_state%immersed%int_xyz_ww(ighost,3)
              call interp_field(current_state%w, ipi, ipj, ipk, xxd, yyd, zzd, cc)
  
              current_state%immersed%gp_store_w(ighost) =  -cc 
    
  !            ! viscosity and diffusion
    
              if (current_state%use_viscosity_and_diffusion)then
                ! parallel velocity at k2 equivalent
                ipi = current_state%immersed%ip_ijk_k2e_wu(ighost,1)
                ipj = current_state%immersed%ip_ijk_k2e_wu(ighost,2)
                ipk = current_state%immersed%ip_ijk_k2e_wu(ighost,3)
                xxd = current_state%immersed%int_xyz_k2e_wu(ighost,1)
                yyd = current_state%immersed%int_xyz_k2e_wu(ighost,2)
                zzd = current_state%immersed%int_xyz_k2e_wu(ighost,3)
                call interp_field(current_state%u, ipi, ipj, ipk, xxd, yyd, zzd, uvel)
      
                ipi = current_state%immersed%ip_ijk_k2e_wv(ighost,1)
                ipj = current_state%immersed%ip_ijk_k2e_wv(ighost,2)
                ipk = current_state%immersed%ip_ijk_k2e_wv(ighost,3)
                xxd = current_state%immersed%int_xyz_k2e_wv(ighost,1)
                yyd = current_state%immersed%int_xyz_k2e_wv(ighost,2)
                zzd = current_state%immersed%int_xyz_k2e_wv(ighost,3)
                call interp_field(current_state%v, ipi, ipj, ipk, xxd, yyd, zzd, vvel)
      
                ipi = current_state%immersed%ip_ijk_k2e_w(ighost,1)
                ipj = current_state%immersed%ip_ijk_k2e_w(ighost,2)
                ipk = current_state%immersed%ip_ijk_k2e_w(ighost,3)
                xxd = current_state%immersed%int_xyz_k2e_w(ighost,1)
                yyd = current_state%immersed%int_xyz_k2e_w(ighost,2)
                zzd = current_state%immersed%int_xyz_k2e_w(ighost,3)
                call interp_field(current_state%w, ipi, ipj, ipk, xxd, yyd, zzd, wvel)
      
    
                pvel_at_k2e = ((uvel+current_state%ugal)*current_state%immersed%norms_w(ighost,1))+&
                              ((vvel+current_state%vgal)*current_state%immersed%norms_w(ighost,2))+&
                              (wvel*current_state%immersed%norms_w(ighost,3))
                tmp1 = (uvel+current_state%ugal) - pvel_at_k2e*current_state%immersed%norms_w(ighost,1)
                tmp2 = (vvel+current_state%vgal) - pvel_at_k2e*current_state%immersed%norms_w(ighost,2)
                tmp3 = wvel - pvel_at_k2e*current_state%immersed%norms_w(ighost,3)
                pvel_at_k2e = sqrt(tmp1**2.0 + tmp2**2.0 + tmp3**2.0) + smallp

                ipi = current_state%immersed%ip_ijk_ww(ighost,1)
                ipj = current_state%immersed%ip_ijk_ww(ighost,2)
                ipk = current_state%immersed%ip_ijk_ww(ighost,3)
                xxd = current_state%immersed%int_xyz_ww(ighost,1)
                yyd = current_state%immersed%int_xyz_ww(ighost,2)
                zzd = current_state%immersed%int_xyz_ww(ighost,3)
                call interp_field(current_state%vis_coefficient, ipi, ipj, ipk, xxd, yyd, zzd, vis_ip)
                call interp_field(current_state%diff_coefficient, ipi, ipj, ipk, xxd, yyd, zzd, diff_ip)
  
                ipi = current_state%immersed%ip_ijk_k2e_ws(ighost,1)
                ipj = current_state%immersed%ip_ijk_k2e_ws(ighost,2)
                ipk = current_state%immersed%ip_ijk_k2e_ws(ighost,3)
                xxd = current_state%immersed%int_xyz_k2e_ws(ighost,1)
                yyd = current_state%immersed%int_xyz_k2e_ws(ighost,2)
                zzd = current_state%immersed%int_xyz_k2e_ws(ighost,3)
                call interp_field(current_state%th, ipi, ipj, ipk, xxd, yyd, zzd, cc)
                thref_ip = current_state%global_grid%configuration%vertical%thref(ipk)*(1.0-zzd)+&
                           current_state%global_grid%configuration%vertical%thref(ipk+1)*zzd
                if(.not. current_state%passive_q)then ! i.e. q is active 
                  call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, tmp1)
                  thref_ip = thref_ip*(1.0_DEFAULT_PRECISION + current_state%cq(iqv)*tmp1)
                end if


                if(current_state%use_surface_boundary_conditions)then
                  if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
  
                    dthv_surf = cc + thref_ip - current_state%immersed%theta_surf_ib_w(ighost)
          
                    convergence_status = mostbc_ib(current_state, pvel_at_k2e, dthv_surf,&
                                         delz_ib(ighost), ighost, ustr, thvstr)
          
                    if (current_state%immersed%dump_ustar .and. current_state%immersed%dodiags) then
                      current_state%immersed%ib_ustar(ighost) = ustr
                    end if
          
                    vistmp = delz_ib(ighost)*ustr**2/pvel_at_k2e
                    difftmp = delz_ib(ighost)*ustr*thvstr/(dthv_surf+smallp)
          
                    ! neutral Prandtl number
                    ustr=von_karman_constant*pvel_at_k2e/zlogm_ib(ighost)
                    visneut=delz_ib(ighost)*ustr**2/(pvel_at_k2e+smallp)
                    diffneut=visneut/prandtl
  
                    tmp1 = current_state%immersed%norms_w(ighost,3)
                    difftmp = difftmp*tmp1 + diffneut*(1.0 - tmp1) 
                    vistmp = vistmp*tmp1 + visneut*(1.0 - tmp1) 
        
                  end if
                

                  if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then
                    tmp1 = current_state%immersed%norms_w(ighost,3)
                    if((tmp1.lt.1.0_DEFAULT_PRECISION).or.(fbuoy_ib(ighost).eq.0.0_DEFAULT_PRECISION))then
                      call handle_neutral_fluxes_ib&
                      (current_state, ighost, pvel_at_k2e, delz_ib(ighost),&
                       zlogth_ib(ighost),zlogm_ib(ighost), visneut, diffneut)
                    end if

                    if(fbuoy_ib(ighost).gt.0.0_DEFAULT_PRECISION)then
                      call handle_convective_fluxes_ib&
                      (current_state, ighost, pvel_at_k2e, delz_ib(ighost),zlogth_ib(ighost), vistmp, difftmp)
               
                    else if(fbuoy_ib(ighost).eq.0.0_DEFAULT_PRECISION)then
                      vistmp = visneut
                      difftmp = diffneut

                    else
                      call handle_stable_fluxes_ib&
                      (current_state, ighost, pvel_at_k2e, delz_ib(ighost),&
                       zlogth_ib(ighost),zlogm_ib(ighost), vistmp, difftmp)
                    end if
                    if(tmp1.lt.1.0_DEFAULT_PRECISION)then
                      difftmp = difftmp*tmp1 + diffneut*(1.0 - tmp1) 
                      vistmp = vistmp*tmp1 + visneut*(1.0 - tmp1) 
                    end if
                    if(current_state%th%active)w2s_store_diff(ighost)=difftmp


                  end if ! fixed fluxes



                else ! use_surface_bcs = false
                  ! neutral
                  ustr=current_state%global_grid%configuration%vertical%vk_on_zlogm*pvel_at_k2e
                  thvstr=0.0_DEFAULT_PRECISION 
                  vistmp=current_state%global_grid%configuration%vertical%zn(2)*ustr**2/pvel_at_k2e
                  difftmp=0.0_DEFAULT_PRECISION

                end if 
                current_state%immersed%gp_store_vis(ighost) = 2.0_DEFAULT_PRECISION*vistmp - vis_ip
                current_state%immersed%gp_store_diff(ighost)= 2.0_DEFAULT_PRECISION*difftmp - diff_ip                
 
              end if ! viscosity & diffusion
      
            end if ! w point
    
  !=============================================================================
    
            ! scalars
  
            ! diagnostics
            if(current_state%immersed%dodiags.and.current_state%immersed%ghost_s(k,j,i).gt.0)then
             ighost = current_state%immersed%ghost_s(k,j,i)
             ipi = current_state%immersed%bi_ijk_s(ighost,1)
             ipj = current_state%immersed%bi_ijk_s(ighost,2)
             ipk = current_state%immersed%bi_ijk_s(ighost,3)
             xxd = current_state%immersed%biint_xyz_s(ighost,1)
             yyd = current_state%immersed%biint_xyz_s(ighost,2)
             zzd = current_state%immersed%biint_xyz_s(ighost,3)
             if (enable_theta .and. current_state%immersed%dump_th) then
               call interp_field(current_state%zth, ipi, ipj, ipk, xxd, yyd, zzd, cc)
               current_state%immersed%ib_th(ighost) = cc
             end if
             if (current_state%immersed%dump_p) then
               call interp_field(current_state%p, ipi, ipj, ipk, xxd, yyd, zzd, cc)
               current_state%immersed%ib_p(ighost) = cc
             end if
            end if
  
  
            if ((enable_theta.or.current_state%number_q_fields .gt. 0)&
                .and.current_state%immersed%ghost_s(k,j,i).gt.0)then ! theta points

                ighost = current_state%immersed%ghost_s(k,j,i)
                
                ipi = current_state%immersed%ip_ijk_ss(ighost,1)
                ipj = current_state%immersed%ip_ijk_ss(ighost,2)
                ipk = current_state%immersed%ip_ijk_ss(ighost,3)
                xxd = current_state%immersed%int_xyz_ss(ighost,1)
                yyd = current_state%immersed%int_xyz_ss(ighost,2)
                zzd = current_state%immersed%int_xyz_ss(ighost,3)
                call interp_field(current_state%th, ipi, ipj, ipk, xxd, yyd, zzd, cc)
  
                if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES)then
                  call interp_w2s(current_state,w2s_store_diff, ighost, tmp1)
                  thref_ip = current_state%global_grid%configuration%vertical%thref(ipk)*(1.0-zzd)+&
                             current_state%global_grid%configuration%vertical%thref(ipk+1)*zzd
                  if(.not. current_state%passive_q)then ! i.e. q is active 
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, tmp2)
                    thref_ip = thref_ip*(1.0_DEFAULT_PRECISION + current_state%cq(iqv)*tmp2)
                  end if
                  current_state%immersed%gp_store_th(ighost)=&
                     (current_state%surface_temperature_flux*2.0*current_state%immersed%gp_len_s(ighost)/&
                     tmp1)+cc-current_state%immersed%thref_surf_ib_w(ighost) +thref_ip

                  ! Surface Flux of vapour
                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, cc)
                    current_state%immersed%gp_store_q(ighost)=cc+&
                    (current_state%surface_vapour_flux*2.0*current_state%immersed%gp_len_s(ighost)/tmp1)
                  endif

                else if(current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES)then
                  tmp1 = current_state%immersed%theta_surf_ib_s(ighost)-&
                         current_state%immersed%thref_surf_ib_s(ighost)
                  current_state%immersed%gp_store_th(ighost) =  2.0_DEFAULT_PRECISION*tmp1 -cc

                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, cc)
                    current_state%immersed%gp_store_q(ighost)=&
                    (current_state%surface_vapour_flux*2.0) - cc
                  endif


                else
                  current_state%immersed%gp_store_th(ighost) = cc 
                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, cc)
                    current_state%immersed%gp_store_q(ighost) = cc
                  end if
                end if
  
            end if
  
  
          end do
        end if ! no slip ib_type






        ! VR
        if (current_state%immersed%ib_type==1)then
!=============================================================================
          do k=1,current_state%immersed%kmax_ji(j,i)
            
            if (current_state%immersed%ghost_u(k,j,i) .gt. 0) then ! u points
              ighost = current_state%immersed%ghost_u(k,j,i)
  
              ipi = current_state%immersed%ip_ijk_uu(ighost,1)
              ipj = current_state%immersed%ip_ijk_uu(ighost,2)
              ipk = current_state%immersed%ip_ijk_uu(ighost,3)
              xxd = current_state%immersed%int_xyz_uu(ighost,1)
              yyd = current_state%immersed%int_xyz_uu(ighost,2)
              zzd = current_state%immersed%int_xyz_uu(ighost,3)
              call interp_field(current_state%u, ipi, ipj, ipk, xxd, yyd, zzd, uvel)

              ipi = current_state%immersed%ip_ijk_uv(ighost,1)
              ipj = current_state%immersed%ip_ijk_uv(ighost,2)
              ipk = current_state%immersed%ip_ijk_uv(ighost,3)
              xxd = current_state%immersed%int_xyz_uv(ighost,1)
              yyd = current_state%immersed%int_xyz_uv(ighost,2)
              zzd = current_state%immersed%int_xyz_uv(ighost,3)
              call interp_field(current_state%v, ipi, ipj, ipk, xxd, yyd, zzd, vvel)

              ipi = current_state%immersed%ip_ijk_uw(ighost,1)
              ipj = current_state%immersed%ip_ijk_uw(ighost,2)
              ipk = current_state%immersed%ip_ijk_uw(ighost,3)
              xxd = current_state%immersed%int_xyz_uw(ighost,1)
              yyd = current_state%immersed%int_xyz_uw(ighost,2)
              zzd = current_state%immersed%int_xyz_uw(ighost,3)
              call interp_field(current_state%w, ipi, ipj, ipk, xxd, yyd, zzd, wvel)

              nhat_x = current_state%immersed%norms_u(ighost,1)
              nhat_y = current_state%immersed%norms_u(ighost,2)
              nhat_z = current_state%immersed%norms_u(ighost,3)
              tmp1 = (uvel*nhat_x)+(vvel*nhat_y)+(wvel*nhat_z) ! projection of vel at IP on to normal
              un_ip = tmp1*nhat_x
              vn_ip = tmp1*nhat_y
              wn_ip = tmp1*nhat_z
              ut_ip = uvel - un_ip
              vt_ip = vvel - vn_ip
              wt_ip = wvel - wn_ip
              tmp2 = sqrt(ut_ip**2.0 + vt_ip**2.0 + wt_ip**2.0) ! tangential speed at IP
              if(tmp2.gt.0.0)that_x = ut_ip/tmp2
              tvel_gp = tmp2*lograt_u(ighost)
              nvel_gp = tmp1*hrat_u(ighost)
              current_state%immersed%gp_store_u(ighost) = tvel_gp*that_x + nvel_gp*nhat_x

                
            end if ! u points -----------------------------------------------
!

            if (current_state%immersed%ghost_v(k,j,i) .gt. 0) then ! v points
              ighost = current_state%immersed%ghost_v(k,j,i)
  
              ipi = current_state%immersed%ip_ijk_vu(ighost,1)
              ipj = current_state%immersed%ip_ijk_vu(ighost,2)
              ipk = current_state%immersed%ip_ijk_vu(ighost,3)
              xxd = current_state%immersed%int_xyz_vu(ighost,1)
              yyd = current_state%immersed%int_xyz_vu(ighost,2)
              zzd = current_state%immersed%int_xyz_vu(ighost,3)
              call interp_field(current_state%u, ipi, ipj, ipk, xxd, yyd, zzd, uvel)

              ipi = current_state%immersed%ip_ijk_vv(ighost,1)
              ipj = current_state%immersed%ip_ijk_vv(ighost,2)
              ipk = current_state%immersed%ip_ijk_vv(ighost,3)
              xxd = current_state%immersed%int_xyz_vv(ighost,1)
              yyd = current_state%immersed%int_xyz_vv(ighost,2)
              zzd = current_state%immersed%int_xyz_vv(ighost,3)
              call interp_field(current_state%v, ipi, ipj, ipk, xxd, yyd, zzd, vvel)

              ipi = current_state%immersed%ip_ijk_vw(ighost,1)
              ipj = current_state%immersed%ip_ijk_vw(ighost,2)
              ipk = current_state%immersed%ip_ijk_vw(ighost,3)
              xxd = current_state%immersed%int_xyz_vw(ighost,1)
              yyd = current_state%immersed%int_xyz_vw(ighost,2)
              zzd = current_state%immersed%int_xyz_vw(ighost,3)
              call interp_field(current_state%w, ipi, ipj, ipk, xxd, yyd, zzd, wvel)


              nhat_x = current_state%immersed%norms_v(ighost,1)
              nhat_y = current_state%immersed%norms_v(ighost,2)
              nhat_z = current_state%immersed%norms_v(ighost,3)
              tmp1 = (uvel*nhat_x)+(vvel*nhat_y)+(wvel*nhat_z) ! normal speed at IP
              un_ip = tmp1*nhat_x
              vn_ip = tmp1*nhat_y
              wn_ip = tmp1*nhat_z
              ut_ip = uvel - un_ip
              vt_ip = vvel - vn_ip
              wt_ip = wvel - wn_ip
              tmp2 = sqrt(ut_ip**2.0 + vt_ip**2.0 + wt_ip**2.0) ! tangential speed at IP
              if(tmp2.gt.0.0)that_y = vt_ip/tmp2
              tvel_gp = tmp2*lograt_v(ighost)
              nvel_gp = tmp1*hrat_v(ighost)
              current_state%immersed%gp_store_v(ighost) = tvel_gp*that_y + nvel_gp*nhat_y
              

            end if ! v points -----------------------------------------------


            if (current_state%immersed%ghost_w(k,j,i) .gt. 0) then ! w points
              ighost = current_state%immersed%ghost_w(k,j,i)
  
              ipi = current_state%immersed%ip_ijk_wu(ighost,1)
              ipj = current_state%immersed%ip_ijk_wu(ighost,2)
              ipk = current_state%immersed%ip_ijk_wu(ighost,3)
              xxd = current_state%immersed%int_xyz_wu(ighost,1)
              yyd = current_state%immersed%int_xyz_wu(ighost,2)
              zzd = current_state%immersed%int_xyz_wu(ighost,3)
              call interp_field(current_state%u, ipi, ipj, ipk, xxd, yyd, zzd, uvel)

              ipi = current_state%immersed%ip_ijk_wv(ighost,1)
              ipj = current_state%immersed%ip_ijk_wv(ighost,2)
              ipk = current_state%immersed%ip_ijk_wv(ighost,3)
              xxd = current_state%immersed%int_xyz_wv(ighost,1)
              yyd = current_state%immersed%int_xyz_wv(ighost,2)
              zzd = current_state%immersed%int_xyz_wv(ighost,3)
              call interp_field(current_state%v, ipi, ipj, ipk, xxd, yyd, zzd, vvel)

              ipi = current_state%immersed%ip_ijk_ww(ighost,1)
              ipj = current_state%immersed%ip_ijk_ww(ighost,2)
              ipk = current_state%immersed%ip_ijk_ww(ighost,3)
              xxd = current_state%immersed%int_xyz_ww(ighost,1)
              yyd = current_state%immersed%int_xyz_ww(ighost,2)
              zzd = current_state%immersed%int_xyz_ww(ighost,3)
              call interp_field(current_state%w, ipi, ipj, ipk, xxd, yyd, zzd, wvel)


              nhat_x = current_state%immersed%norms_w(ighost,1)
              nhat_y = current_state%immersed%norms_w(ighost,2)
              nhat_z = current_state%immersed%norms_w(ighost,3)
              tmp1 = (uvel*nhat_x)+(vvel*nhat_y)+(wvel*nhat_z) ! normal speed at IP
              un_ip = tmp1*nhat_x
              vn_ip = tmp1*nhat_y
              wn_ip = tmp1*nhat_z
              ut_ip = uvel - un_ip
              vt_ip = vvel - vn_ip
              wt_ip = wvel - wn_ip
              tmp2 = sqrt(ut_ip**2.0 + vt_ip**2.0 + wt_ip**2.0) ! tangential speed at IP
              if(tmp2.gt.0.0)that_z = wt_ip/tmp2
              tvel_gp = tmp2*lograt_w(ighost)
              nvel_gp = tmp1*hrat_w(ighost)
              current_state%immersed%gp_store_w(ighost) = tvel_gp*that_z + nvel_gp*nhat_z
              tvel_ip=tmp2


              ! viscosity and diffusivity
              if (current_state%use_viscosity_and_diffusion)then
                if (enable_theta) then
                  ! get thetav at IP

                  ipi = current_state%immersed%ip_ijk_ws(ighost,1)
                  ipj = current_state%immersed%ip_ijk_ws(ighost,2)
                  ipk = current_state%immersed%ip_ijk_ws(ighost,3)
                  xxd = current_state%immersed%int_xyz_ws(ighost,1)
                  yyd = current_state%immersed%int_xyz_ws(ighost,2)
                  zzd = current_state%immersed%int_xyz_ws(ighost,3)
                  call interp_field(current_state%th, ipi, ipj, ipk, xxd, yyd, zzd, tmp1)

                  if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES) then
                    tmp2 = current_state%immersed%theta_surf_ib_w(ighost)-&
                           current_state%immersed%thref_surf_ib_w(ighost)
                    dthv_surf = tmp1-tmp2
                    convergence_status = mostbc_ib(current_state, tvel_ip, dthv_surf,&
                      current_state%immersed%ip_len_w(ighost), ighost, ustr, thvstr)
    
                    vistmp =(current_state%immersed%ip_len_w(ighost)*ustr**2)/(tvel_ip+smallp)
                    difftmp = current_state%immersed%ip_len_w(ighost)*ustr*thvstr/(dthv_surf+smallp)
          
                    ! neutral Prandtl number
                    ustr=von_karman_constant*tvel_ip/zlogm_ib(ighost)
                    visneut=current_state%immersed%ip_len_w(ighost)*ustr**2/(tvel_ip+smallp)
                    diffneut=visneut/prandtl
         
                    tmp1 = current_state%immersed%norms_w(ighost,3)
                    difftmp = difftmp*tmp1 + diffneut*(1.0 - tmp1) 
                    vistmp = vistmp*tmp1 + visneut*(1.0 - tmp1) 
                  end if !surface_values
                  
                  if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES) then

                    tmp1 = current_state%immersed%norms_w(ighost,3)
                    if((tmp1.lt.1.0_DEFAULT_PRECISION).or.(fbuoy_ib(ighost).eq.0.0_DEFAULT_PRECISION))then
                      call handle_neutral_fluxes_ib&
                      (current_state, ighost, tvel_ip, delz_ib(ighost),&
                       zlogth_ib(ighost),zlogm_ib(ighost), visneut, diffneut)
                    end if

                    if(fbuoy_ib(ighost).gt.0.0_DEFAULT_PRECISION)then
                      call handle_convective_fluxes_ib&
                      (current_state, ighost, tvel_ip, delz_ib(ighost),&
                      zlogth_ib(ighost), vistmp, difftmp)
                    else if(fbuoy_ib(ighost).eq.0.0_DEFAULT_PRECISION)then
                      vistmp = visneut
                      difftmp = diffneut
                    else
                      call handle_stable_fluxes_ib&
                      (current_state, ighost, tvel_ip, delz_ib(ighost),&
                       zlogth_ib(ighost),zlogm_ib(ighost), vistmp, difftmp)
                    end if

                    if(tmp1.lt.1.0_DEFAULT_PRECISION)then
                      difftmp = difftmp*tmp1 + diffneut*(1.0 - tmp1) 
                      vistmp = vistmp*tmp1 + visneut*(1.0 - tmp1) 
                    end if
                    w2s_store_diff(ighost)=difftmp
                  end if ! surface_fluxes

                else ! no theta
                  ! neutral
                  ustr=von_karman_constant*tvel_ip/zlogm_ib(ighost)
                  thvstr=0.0_DEFAULT_PRECISION
                  vistmp=current_state%immersed%ip_len_w(ighost)*ustr**2/tvel_ip
                  difftmp=0.0_DEFAULT_PRECISION
                end if
  
                ! get visc and diff at IP
                ipi = current_state%immersed%ip_ijk_ww(ighost,1)
                ipj = current_state%immersed%ip_ijk_ww(ighost,2)
                ipk = current_state%immersed%ip_ijk_ww(ighost,3)
                xxd = current_state%immersed%int_xyz_ww(ighost,1)
                yyd = current_state%immersed%int_xyz_ww(ighost,2)
                zzd = current_state%immersed%int_xyz_ww(ighost,3)
                call interp_field(current_state%vis_coefficient, ipi, ipj, ipk, xxd, yyd, zzd, tmp1)
                call interp_field(current_state%diff_coefficient, ipi, ipj, ipk, xxd, yyd, zzd, tmp2)

                current_state%immersed%gp_store_vis(ighost) =  &
                hrat_w(ighost)*(tmp1-vistmp) + vistmp
                current_state%immersed%gp_store_diff(ighost) =  &
                hrat_w(ighost)*(tmp2-difftmp) + difftmp

              end if ! if use viscosity and diffusion
  
            end if ! w points -----------------------------------------------



            if (enable_theta) then
              if (current_state%immersed%ghost_s(k,j,i) .gt. 0) then ! s points
                ighost = current_state%immersed%ghost_s(k,j,i)
    
                ipi = current_state%immersed%ip_ijk_ss(ighost,1)
                ipj = current_state%immersed%ip_ijk_ss(ighost,2)
                ipk = current_state%immersed%ip_ijk_ss(ighost,3)
                xxd = current_state%immersed%int_xyz_ss(ighost,1)
                yyd = current_state%immersed%int_xyz_ss(ighost,2)
                zzd = current_state%immersed%int_xyz_ss(ighost,3)
                call interp_field(current_state%th, ipi, ipj, ipk, xxd, yyd, zzd, cc)

                if (current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_FLUXES)then
                  thref_ip = current_state%global_grid%configuration%vertical%thref(ipk)*(1.0-zzd)+&
                             current_state%global_grid%configuration%vertical%thref(ipk+1)*zzd

                  if(.not. current_state%passive_q)then ! i.e. q is active 
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, tmp2)
                    thref_ip = thref_ip*(1.0_DEFAULT_PRECISION + current_state%cq(iqv)*tmp2)
                  end if

                  ipi = current_state%immersed%ip_ijk_sw(ighost,1)
                  ipj = current_state%immersed%ip_ijk_sw(ighost,2)
                  ipk = current_state%immersed%ip_ijk_sw(ighost,3)
                  xxd = current_state%immersed%int_xyz_sw(ighost,1)
                  yyd = current_state%immersed%int_xyz_sw(ighost,2)
                  zzd = current_state%immersed%int_xyz_sw(ighost,3)
                  call interp_field(current_state%diff_coefficient, ipi, ipj, ipk, xxd, yyd, zzd, tmp1)
                  call interp_w2s(current_state,w2s_store_diff, ighost, tmp2)

                  ! diffusivity at VR-IP mid point
                  tmp3 = hrat_diff_s(ighost)*(tmp1-tmp2) + tmp2
                  
                  current_state%immersed%gp_store_th(ighost)=&
                  (current_state%surface_temperature_flux*(current_state%immersed%ip_len_s(ighost)-&
                   current_state%immersed%gp_len_s(ighost))/tmp3)+&
                   cc-current_state%immersed%thref_surf_ib_s(ighost) +thref_ip

                  ! Surface Flux of vapour
                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, tmp1)
                    current_state%immersed%gp_store_q(ighost)=&
                  (current_state%surface_vapour_flux*(current_state%immersed%ip_len_s(ighost)-&
                   current_state%immersed%gp_len_s(ighost))/tmp3)+tmp1
                  end if


                else if(current_state%type_of_surface_boundary_conditions == PRESCRIBED_SURFACE_VALUES)then
                  tmp2 = current_state%immersed%theta_surf_ib_s(ighost)-&
                         current_state%immersed%thref_surf_ib_s(ighost)
                  current_state%immersed%gp_store_th(ighost) =  &
                  hrat_s(ighost)*(cc-tmp2) + tmp2

                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, cc)
                    current_state%immersed%gp_store_q(ighost)=&
                    hrat_s(ighost)*(cc-current_state%surface_vapour_flux)+&
                    current_state%surface_vapour_flux
                  end if

                else
                  current_state%immersed%gp_store_th(ighost) = cc
                  if (current_state%number_q_fields .gt. 0) then
                    call interp_field(current_state%q(iqv), ipi, ipj, ipk, xxd, yyd, zzd, cc)
                    current_state%immersed%gp_store_q(ighost) = cc
                  end if
                end if
  
              end if ! s points -----------------------------------------------
            end if ! if theta enabled

          end do

        end if ! VR ib_type

    

      end if ! if ib_col
    end if ! if ib_enabled


  end subroutine timestep_callback

  !> Finalisation callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
  end subroutine finalisation_callback

  subroutine interp_field(xf, ipi, ipj, ipk, xxd, yyd, zzd, val)
    type(prognostic_field_type), intent(in) :: xf
    integer, intent(in) :: ipi, ipj, ipk
    real(kind=DEFAULT_PRECISION), intent(in) :: xxd, yyd, zzd
    real(kind=DEFAULT_PRECISION), intent(out) :: val
    real(kind=DEFAULT_PRECISION) :: c00, c01, c10, c11, c0, c1, cc
    c00 = xf%data(ipk,ipj,ipi)*(1.0-xxd)+xf%data(ipk,ipj,ipi+1)*(xxd)
    c01 = xf%data(ipk+1,ipj,ipi)*(1.0-xxd)+xf%data(ipk+1,ipj,ipi+1)*(xxd)
    c10 = xf%data(ipk,ipj+1,ipi)*(1.0-xxd)+xf%data(ipk,ipj+1,ipi+1)*(xxd)
    c11 = xf%data(ipk+1,ipj+1,ipi)*(1.0-xxd)+xf%data(ipk+1,ipj+1,ipi+1)*(xxd)
    c0 = c00*(1.0-yyd) + c10*yyd
    c1 = c01*(1.0-yyd) + c11*yyd
    val = c0*(1.0-zzd) + c1*zzd
  end subroutine interp_field

  subroutine interp_w2s(current_state, inarr, ighost, val)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), intent(in) :: inarr(:)
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(out) :: val
    val=&
    (current_state%immersed%w2s_idw_wgt(ighost,1)*inarr(current_state%immersed%w2s_idw_idx(ighost,1)))+&
    (current_state%immersed%w2s_idw_wgt(ighost,2)*inarr(current_state%immersed%w2s_idw_idx(ighost,2)))+&
    (current_state%immersed%w2s_idw_wgt(ighost,3)*inarr(current_state%immersed%w2s_idw_idx(ighost,3)))+&
    (current_state%immersed%w2s_idw_wgt(ighost,4)*inarr(current_state%immersed%w2s_idw_idx(ighost,4)))

  end subroutine interp_w2s


  ! MOST solver taken from lowerbc.F90, adapted for IB
  integer function mostbc_ib(current_state, delu, delt, delz, ighost, ustrdg, tstrdg)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), intent(in)  :: delu, delt, delz
    integer, intent(in)  :: ighost
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg
    real(kind=DEFAULT_PRECISION) :: ellmocon

    if (delt .lt. 0.0_DEFAULT_PRECISION) then
      if (delu .le. smallp) then 
        ustrdg=0.0_DEFAULT_PRECISION
        tstrdg=tstrcona_ib(ighost)*delt
        mostbc_ib=CONVERGENCE_SUCCESS
      else
         ! The unstable case
         ellmocon=current_state%immersed%theta_surf_ib_w(ighost)/(G*von_karman_constant)
         mostbc_ib=solve_monin_obukhov_unstable_case_ib(delu, delt, ellmocon, delz, ustrdg, tstrdg, ighost)
      end if

    else if (delt .gt. 0.0_DEFAULT_PRECISION) then
      ! The stable case
      mostbc_ib=solve_monin_obukhov_stable_case_ib(delu, delt, zlogm_ib(ighost), &
           cmbc_ib(ighost), ighost, ustrdg, tstrdg)

    else
      ! Trivial neutral case
        ustrdg=delu*von_karman_constant/zlogm_ib(ighost)
        tstrdg=0.0_DEFAULT_PRECISION 
        mostbc_ib=CONVERGENCE_SUCCESS
    end if

    if (mostbc_ib .ne. CONVERGENCE_SUCCESS) then
      if (mostbc_ib .eq. CONVERGENCE_RICHARDSON_TOO_LARGE) then
        !call log_log(LOG_WARN, "Richardson number greater than critical value")
      else if(mostbc_ib .eq. CONVERGENCE_FAILURE) then
        write(*,*)ighost,delu, delt, delz
        call log_log(LOG_ERROR, "Convergence failure (MOSTBC_IB) after 200 iterations")
      end if
    end if   
  end function mostbc_ib

  integer function solve_monin_obukhov_unstable_case_ib(delu, delt, ellmocon, delz, ustrdg, tstrdg, ighost)
    real(kind=DEFAULT_PRECISION), intent(in) :: delu, delt, ellmocon, delz
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg
    real(kind=DEFAULT_PRECISION), parameter :: smth = 0.05_DEFAULT_PRECISION,& ! Smoothing between iterations
        tolm=1.0E-4_DEFAULT_PRECISION,  tolt=1.0E-4_DEFAULT_PRECISION ! Convergence tollerance for u and t star
    integer :: i
    real(kind=DEFAULT_PRECISION) :: ellmo, psim, psih, x4, xx, xx0, y2, yy, yy0, err_ustr, err_tstr, &
         ustrl, tstrl, & ! U and T star at start of iteration
         ustrpl, tstrpl,& ! U and T star at end of iteration
         x4con, y2con
    
    ! First set initial values
    ustrl=von_karman_constant*delu/zlogm_ib(ighost)
    tstrl=tstrcona_ib(ighost)*delt
    
    x4con=gammam*(delz+z0)
    y2con=gammah*(delz+z0)

    ! Now start iteration
    do i=1, 200
      ellmo=ustrl*ustrl*ellmocon/tstrl

      ! Test for possible square root of negative quantity
      x4=1.0_DEFAULT_PRECISION-(x4con)/ellmo      
      if (x4 .lt. 0.0_DEFAULT_PRECISION) call log_log(LOG_ERROR, "Negative square root in x4")        

      xx=sqrt(sqrt(x4))
      xx0=sqrt(sqrt(1.0_DEFAULT_PRECISION-xx0con / ellmo))        
      psim=2.*( log((xx+1.0_DEFAULT_PRECISION)/(xx0+1.0_DEFAULT_PRECISION))-atan(xx)+atan(xx0) )+&
           log((xx*xx+1.0_DEFAULT_PRECISION)/(xx0*xx0+1.0_DEFAULT_PRECISION))
      ustrpl=von_karman_constant*delu/(zlogm_ib(ighost)-psim)

      ! Test for possible square root of negative quantity
      y2=1.-y2con/ellmo
      if (y2 .lt. 0.0_DEFAULT_PRECISION) call log_log(LOG_ERROR, "Negative square root in y2")
      yy=sqrt(y2)
      yy0=sqrt(1.0_DEFAULT_PRECISION-yy0con/ellmo)
      psih=2.*log((1.0_DEFAULT_PRECISION+yy)/(1.0_DEFAULT_PRECISION+yy0))
      tstrpl=tstrconb*delt/(zlogth_ib(ighost)-psih)
      err_ustr=abs((ustrpl-ustrl)/ ustrl)
      err_tstr=abs((tstrpl-tstrl)/ tstrl)                                       
      if ((err_tstr .lt. tolt) .and. (err_ustr .lt. tolm))  then                                                 
        ustrdg=ustrpl
        tstrdg=tstrpl
        solve_monin_obukhov_unstable_case_ib=CONVERGENCE_SUCCESS
        return
      else                                                                 
        ustrl=(1.0_DEFAULT_PRECISION-smth)*ustrpl+smth*ustrl
        tstrl=(1.0_DEFAULT_PRECISION-smth)*tstrpl+smth*tstrl
      end if
    end do
    solve_monin_obukhov_unstable_case_ib=CONVERGENCE_FAILURE
  end function solve_monin_obukhov_unstable_case_ib

  integer function solve_monin_obukhov_stable_case_ib(delu, delt, zlogm, cmbc, ighost, ustrdg, tstrdg)
    real(kind=DEFAULT_PRECISION), intent(in) :: delu, delt, zlogm, cmbc
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(out) :: ustrdg, tstrdg

    real(kind=DEFAULT_PRECISION) :: am, ah, ee, ff, det

    am=von_karman_constant*delu
    ah=von_karman_constant*delt
    ee=am*eecon_ib(ighost)
    ff=ah*cmbc-rhmbc_ib(ighost)*am*am  !
    det=ee*ee-(4.*ddbc_ib(ighost))*ff
    solve_monin_obukhov_stable_case_ib=CONVERGENCE_SUCCESS
    ! Test for laminar flow
    if (ff .gt. 0.0_DEFAULT_PRECISION) then
      if ((ee .lt. 0.0_DEFAULT_PRECISION).and.(det .gt. 0.0_DEFAULT_PRECISION)) then
        ustrdg=(-ee+sqrt(det))*0.5_DEFAULT_PRECISION/ddbc_ib(ighost)
        tstrdg=ustrdg*(am-zlogm*ustrdg)*1.0_DEFAULT_PRECISION/cmbc_ib(ighost)
      else
        solve_monin_obukhov_stable_case_ib=CONVERGENCE_RICHARDSON_TOO_LARGE        
        ustrdg=0.0_DEFAULT_PRECISION
        tstrdg=0.0_DEFAULT_PRECISION
      end if      
    else if (ddbc_ib(ighost) .eq. 0.0_DEFAULT_PRECISION) then
      ! Degenerate case
      ustrdg=-ff/ee
      tstrdg=delt*ustrdg/delu
    else
      ! Solve quadratic for USTRDG
      ustrdg=(-ee+sqrt(det))*0.5_DEFAULT_PRECISION/ddbc_ib(ighost)
      tstrdg=ustrdg*(am-zlogm*ustrdg)*1.0_DEFAULT_PRECISION/cmbc_ib(ighost)
    end if
  end function solve_monin_obukhov_stable_case_ib


! PRESCRIBED FLUXES subroutines ***********************************************
  subroutine handle_convective_fluxes_ib(current_state, ighost, pvel, dz, zlogth, vis, diff)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(in) :: pvel, dz, zlogth
    real(kind=DEFAULT_PRECISION), intent(out) :: vis, diff

    integer :: n
    real(kind=DEFAULT_PRECISION) :: ustr

    ustr=look_ib(current_state, ighost, pvel)

    vis=dz*ustr**2/pvel
    diff=(von_karman_constant*dz*ustr/alphah)/&
         (zlogth- 2.*log((1.+sqrt(1.+gammah*von_karman_constant*fbuoy_ib(ighost)*(dz+z0)&
         /ustr**3))/ (1.+sqrt(1.+gammah*von_karman_constant*fbuoy_ib(ighost)*z0th/ustr**3))))


  end subroutine handle_convective_fluxes_ib

  subroutine handle_neutral_fluxes_ib(current_state, ighost, pvel, dz, zlogth, zlogm, vis, diff)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(in) :: pvel, dz, zlogth, zlogm
    real(kind=DEFAULT_PRECISION), intent(out) :: vis, diff

    integer :: n
    real(kind=DEFAULT_PRECISION) :: ustr


    ustr=pvel*von_karman_constant/zlogm
    vis=dz*ustr**2/pvel
    diff=vis*zlogm/(alphah*zlogth)

  end subroutine handle_neutral_fluxes_ib

  subroutine handle_stable_fluxes_ib(current_state, ighost, pvel, dz, zlogth, zlogm, vis, diff)
    type(model_state_type), target, intent(inout) :: current_state
    integer, intent(in) :: ighost
    real(kind=DEFAULT_PRECISION), intent(in) :: pvel, dz, zlogth, zlogm
    real(kind=DEFAULT_PRECISION), intent(out) :: vis, diff

    integer :: n
    real(kind=DEFAULT_PRECISION) :: ustr


    ustr=pvel*von_karman_constant/zlogm
    if((fbuoy_ib(ighost) - 1.E-9_DEFAULT_PRECISION) .lt. -4.0_DEFAULT_PRECISION*&
         von_karman_constant**2*pvel**3/ (27.0_DEFAULT_PRECISION*betam*dz*zlogm**2)) then
      ! Too stable for turbulence
      vis=0.0_DEFAULT_PRECISION
      diff=0.0_DEFAULT_PRECISION
    else
      ustr=ustr/3.0_DEFAULT_PRECISION*(1.0_DEFAULT_PRECISION-2.0_DEFAULT_PRECISION*&
           cos((acos(-27.0_DEFAULT_PRECISION*betam*von_karman_constant*dz*fbuoy_ib(ighost)/(zlogm*&
           2.0_DEFAULT_PRECISION*ustr**3)-1.0_DEFAULT_PRECISION)+ 2.0_DEFAULT_PRECISION*pi)/3.0_DEFAULT_PRECISION))
      vis=dz*ustr**2/pvel
      diff=vis*(zlogm-betam*dz*von_karman_constant*fbuoy_ib(ighost)/ustr**3)/(alphah*zlogth-betah*&
           von_karman_constant*fbuoy_ib(ighost)*(dz+z0-z0th)/ustr**3)
    end if


  end subroutine handle_stable_fluxes_ib

  subroutine set_look_ib(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer I, ik, ilook, ighost    ! Loop counters
    real(kind=DEFAULT_PRECISION) :: smth, &       ! Relaxation parameter for unstable case
         ob, x1, x0, zlogm, fbuoy, dz

    velmax=100.0_DEFAULT_PRECISION
    velmin=0.1_DEFAULT_PRECISION
    aloginv=1.0_DEFAULT_PRECISION/log(velmax/velmin)
    smth=0.1_DEFAULT_PRECISION  ! _relaxation parameter for unstable case


    do ilook=1, nunst
      ighost=lookup_map_t2g(ilook)

      do ik=1, current_state%lookup_table_entries
        lookup_ib_vel(ilook,ik)=velmin*(velmax/velmin)**&
             (real(ik-1)/real(current_state%lookup_table_entries-1))            
        lookup_ib_ustr(ilook,ik)=lookup_ib_vel(ilook,ik)*von_karman_constant/zlogm_ib(ighost)
        do i=1, 30 
          ob=-lookup_ib_ustr(ilook,ik)**3/(von_karman_constant*fbuoy_ib(ighost))
          x1=sqrt(sqrt(1.-gammam*(delz_ib(ighost)+z0)/ob))
          x0=sqrt(sqrt(1.-gammam*z0/ob))
          lookup_ib_ustr(ilook,ik)=(1.-smth)*lookup_ib_ustr(ilook,ik) + smth*&
               (von_karman_constant*lookup_ib_vel(ilook,ik)) / (zlogm_ib(ighost)-(2.0_DEFAULT_PRECISION*&
               log((x1+1.0_DEFAULT_PRECISION)/(x0+1.0_DEFAULT_PRECISION)) + log((x1*x1+1.0_DEFAULT_PRECISION)/&
               (x0*x0+1.0_DEFAULT_PRECISION)) + 2.0_DEFAULT_PRECISION*atan(x0) -2.0_DEFAULT_PRECISION*atan(x1)))
        end do
      end do
    end do
    ! apply out-of-bounds conditions (cneut, cfc) during lookup
  end subroutine set_look_ib
 
  real(kind=DEFAULT_PRECISION) function look_ib(current_state, ighost, vel)
    type(model_state_type), target, intent(inout) :: current_state
    real(kind=DEFAULT_PRECISION), intent(in) :: vel     ! Horizontal speed at lowest model level
    integer, intent(in) :: ighost   

    real(kind=DEFAULT_PRECISION) :: lookup_real_posn
    integer :: lookup_int_posn, ilook

    lookup_real_posn=1.0_DEFAULT_PRECISION+real(current_state%lookup_table_entries-1)*&
         log(vel/velmin)*aloginv
    lookup_int_posn=int(lookup_real_posn)                                                         

    ilook=lookup_map_g2t(ighost)
    if (lookup_int_posn .ge. 1) then                                                       
      if (lookup_int_posn .lt. current_state%lookup_table_entries) then      ! Linear interpolation
        look_ib=lookup_ib_ustr(ilook,lookup_int_posn)+ (lookup_real_posn-real(lookup_int_posn))*&
             (lookup_ib_ustr(ilook,lookup_int_posn+1)-lookup_ib_ustr(ilook,lookup_int_posn))
      else     ! Near neutral
        look_ib=vel*lookup_ib_ustr(ilook,current_state%lookup_table_entries)/&
             lookup_ib_vel(ilook,current_state%lookup_table_entries)
      end if
    else       ! Nearly free convection                                     
      look_ib=vel**(-convective_limit)*lookup_ib_ustr(ilook,1)*lookup_ib_vel(ilook,1)**convective_limit  ! _Businger-Dyer
    end if
  end function look_ib


end module immersed_boundary_mod
