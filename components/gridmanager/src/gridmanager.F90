!> Manages the grid based upon the model state_mod
module gridmanager_mod
  use datadefn_mod, only : STRING_LENGTH
  use monc_component_mod, only : component_descriptor_type
  use state_mod, only : model_state_type, PRESCRIBED_SURFACE_FLUXES
  use optionsdatabase_mod, only : options_get_integer, options_get_logical, options_get_real, &
     options_get_logical_array, options_get_real_array, options_get_string_array, options_get_array_size, options_get_string
  use grids_mod, only : vertical_grid_configuration_type, X_INDEX, Y_INDEX, Z_INDEX
  use logging_mod, only : LOG_INFO, LOG_ERROR, log_master_log, log_log
  use conversions_mod, only : conv_to_string
  use datadefn_mod, only : DEFAULT_PRECISION, PRECISION_TYPE
  use science_constants_mod, only : von_karman_constant, z0, z0th, cp, r_over_cp, r, g, rlvap_over_cp
  use saturation_mod, only : qsaturation, dqwsatdt
  use q_indices_mod, only: get_q_index, standard_q_names
  use interpolation_mod, only: piecewise_linear_1d
  use rcemip_mod, only: rcemip_init
  use tracers_mod, only: reinitialise_trajectories
  use mpi, only: MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
  use netcdf, only : nf90_noerr, nf90_global, nf90_nowrite, nf90_inquire_attribute, nf90_open, nf90_strerror, &
       nf90_inq_dimid, nf90_inquire_dimension, nf90_inq_varid, nf90_get_var, nf90_inquire, nf90_close, nf90_get_att
  use logging_mod, only : LOG_INFO, LOG_WARN, LOG_ERROR, LOG_DEBUG, log_master_log, log_log, log_get_logging_level


  implicit none

#ifndef TEST_MODE
  private
#endif
  
  ! 1 = No adjustment
  ! 2 = ensure P0=PSF by adjusting THREF profile by constant factor
  ! 3 = ensure P0=PSF by adjusting PSF (not advised)
  ! 4 = ensure P0=PSF by adjusting PTOP
  integer, parameter :: ANELASTIC_PROFILE_MODE=4
  real, parameter :: DEFAULT_SPACING = 1.E9  !< The default spacing used if no grid is active in a specific dimension
  real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: qinit
  integer :: iqv  ! index for vapour

  public gridmanager_get_descriptor

contains

  !> Provides the descriptor back to the caller and is used in component registration
  !! @returns The GridManager component descriptor
  type(component_descriptor_type) function gridmanager_get_descriptor()
    gridmanager_get_descriptor%name="grid_manager"
    gridmanager_get_descriptor%version=0.1
    gridmanager_get_descriptor%initialisation=>initialise_callback
    gridmanager_get_descriptor%finalisation=>finalise_callback
  end function gridmanager_get_descriptor

  !> Called during initialisation and will initialise the horizontal and vertical grid configurations
  !! Note that the model state_mod (from a checkpoint or external file) must have been initialised already
  !! @param current_state The current model state_mod
  subroutine initialise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    integer :: dimensions, k

    if (.not. current_state%initialised) then
      call log_log(LOG_ERROR, "Must initialise the model state_mod before constructing the grid properties")
    end if

    call initialise_horizontalgrid_configuration_types(current_state)
    call initialise_verticalgrid_configuration_type(current_state)
    current_state%immersed%ib_enabled= & 
       options_get_logical(current_state%options_database, "immersed_boundary_enabled")
    if (current_state%immersed%ib_enabled) then
      call initialise_immersed_boundary(current_state)
    end if
    dimensions=1
    if (current_state%global_grid%active(X_INDEX)) dimensions = dimensions+1
    if (current_state%global_grid%active(Y_INDEX)) dimensions = dimensions+1   
    call log_master_log(LOG_INFO, trim(conv_to_string(dimensions))//"D system; z="//&
         trim(conv_to_string(current_state%global_grid%size(Z_INDEX)))//", y="//&
         trim(conv_to_string(current_state%global_grid%size(Y_INDEX)))//", x="//&
         trim(conv_to_string(current_state%global_grid%size(X_INDEX))))
  end subroutine initialise_callback

  !> Called as MONC exits, for good practice this will deallocate the memory used for the grids
  !! @param current_state The current model state_mod
  subroutine finalise_callback(current_state)
    type(model_state_type), target, intent(inout) :: current_state

    type(vertical_grid_configuration_type) :: vertical_grid

    vertical_grid=current_state%global_grid%configuration%vertical

    deallocate(vertical_grid%z, vertical_grid%zn, vertical_grid%dz, vertical_grid%dzn, vertical_grid%czb, vertical_grid%cza, &
         vertical_grid%czg, vertical_grid%czh, vertical_grid%rdz, vertical_grid%rdzn, vertical_grid%tzc1, vertical_grid%tzc2,&
         vertical_grid%tzd1, vertical_grid%tzd2, vertical_grid%thref, vertical_grid%theta_init,  vertical_grid%temp_init, &
         vertical_grid%rh_init, vertical_grid%tref, &
         vertical_grid%prefn, vertical_grid%pdiff, vertical_grid%prefrcp, vertical_grid%rprefrcp, vertical_grid%rho, &
         vertical_grid%rhon, vertical_grid%tstarpr, vertical_grid%qsat, vertical_grid%dqsatdt, vertical_grid%qsatfac, &
         vertical_grid%dthref, vertical_grid%rneutml, vertical_grid%rneutml_sq, vertical_grid%buoy_co, &
         vertical_grid%u_init, vertical_grid%v_init, vertical_grid%q_init,                              &
         vertical_grid%q_rand, vertical_grid%theta_rand, vertical_grid%w_subs, vertical_grid%w_rand, &
         vertical_grid%q_force, vertical_grid%theta_force, vertical_grid%u_force, vertical_grid%v_force &
         )
  end subroutine finalise_callback  

  !> Will initialise the immersed boundary ghost points, image points, surface
  !  exchange points etc.
  !! @param current_state The current model state_mod
  subroutine initialise_immersed_boundary(current_state)
    type(model_state_type), intent(inout) :: current_state
    integer :: nni, nnj, x_size, y_size, idiag, i_diag
    integer(DEFAULT_PRECISION) :: nnodes
    integer :: i_nrst, j_nrst, i_surf, j_surf, ip0, jp0, kp0, ighost, jghost
    integer :: nelem, nsurf_u, nsurf_v, nsurf_w, nsurf_s
    integer :: nghosts_u, nghosts_v, nghosts_w, nghosts_s
    integer :: i0, j0, xhalo, yhalo, myproc, ncid, ierr, nxg, nyg, nfloats
    integer :: ip1, jp1, i_min, j_min, i_min2, j_min2, kmax, ip_dir
    integer :: i,j,k, ii, jj, kk,  i_ib, j_ib , iid, jjd, ielem, p, q, r ! loop counters
    integer, parameter :: MAX_FILE_LEN=200       !< Maximum length of surface input filename
    integer, dimension(:), allocatable :: elem_lst ! element search indices
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: elem_cent  ! element centre positions (x,y,z)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: elem_xlims  ! element x limits
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: elem_ylims  ! element y limits
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: elem_zlims  ! element z limits
    real(DEFAULT_PRECISION) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6
    real(DEFAULT_PRECISION) :: tmp1u, tmp1v, tmp1w, tmp1s, proj
    real(DEFAULT_PRECISION) :: nrmi, nrmj, nrmk, maxh, llen
    real(DEFAULT_PRECISION) :: xdd, ydd, zdd, xl, xr, yl, yr
    real(DEFAULT_PRECISION) :: x_ug, y_ug, z_ug, dx, dy, dz
    real(DEFAULT_PRECISION) :: x_vg, y_vg, z_vg
    real(DEFAULT_PRECISION) :: x_wg, y_wg, z_wg
    real(DEFAULT_PRECISION) :: x_sg, y_sg, z_sg
    real(DEFAULT_PRECISION) :: x_g, y_g, z_g
    real(DEFAULT_PRECISION) :: x_bi, y_bi, z_bi
    real(DEFAULT_PRECISION) :: x_ip, y_ip, z_ip
    real(DEFAULT_PRECISION) :: prime_wgt
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid
    character (len=50):: fname
    character(MAX_FILE_LEN) :: input_file
    logical :: lwrite
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: nodes, norms  ! node locations (x,y,z), and unit normal components (i,j,k)
    integer(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: elems  ! node indices associated with element (x4, clockwise)
    character(len=*), parameter :: NNODES_KEY =  "nnodes"  
    character(len=*), parameter :: NODES_KEY =  "nodes"  
    character(len=*), parameter :: ELEMS_KEY =  "elems"  
    character(len=*), parameter :: NORMS_KEY =  "norms"  
    integer, dimension(:,:), allocatable :: tmp_ijk_u, tmp_ijk_v, tmp_ijk_w, tmp_ijk_s
    integer :: i00, i01, i10, i11, inode
 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: diag_xyz ! real diagnostic array
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: diag2_xyz ! real diagnostic array
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: diag3_xyz ! real diagnostic array
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: gp_xyz_u, ip_xyz_u, bi_xyz_u 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: gp_xyz_v, ip_xyz_v, bi_xyz_v 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: gp_xyz_w, ip_xyz_w, bi_xyz_w, diag4_xyz_w 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: gp_xyz_s, ip_xyz_s, bi_xyz_s 
    real(kind=DEFAULT_PRECISION), parameter :: projtol=1.e-10

    logical :: u_filter_enabled, v_filter_enabled,  w_filter_enabled, th_filter_enabled, vert_filter_enabled ! filter switches for ib proximity array
    integer, parameter :: idw_nmax=4  ! maximum stencil size for IDW (inverse distance weighting)
    integer :: gsearchrad, jstart, jstop, surface_boundary_type !IDW variables
    real(kind=DEFAULT_PRECISION) :: idw_hmax, idw_vmax, dmax

    real(kind=DEFAULT_PRECISION) :: ib_filter_prox
    logical :: enable_theta=.false., use_surface_bcs=.false.


    ! Determine vapour index
    if (.not. current_state%passive_q) then 
       iqv = get_q_index(standard_q_names%VAPOUR, 'gridmanager (IB)')
    endif

    use_surface_bcs=options_get_logical(current_state%options_database,&
                    "use_surface_boundary_conditions")
    if(use_surface_bcs)&
    surface_boundary_type=options_get_integer(current_state%options_database,&
                          "type_of_surface_boundary_conditions")
    
    enable_theta=options_get_logical(current_state%options_database, "enable_theta")
    input_file=options_get_string(current_state%options_database, "ib_input_file")
    myproc = current_state%parallel%my_rank
    if (myproc == 0) then 
      call check_status(nf90_open(path = trim(input_file), mode = nf90_nowrite, ncid = ncid))
      call read_dims(trim(input_file), ncid, nnodes, NNODES_KEY)
      call check_status(nf90_close(ncid))
    end if ! proc 0 only

    call mpi_bcast(nnodes, 1, PRECISION_TYPE, 0, &
                   current_state%parallel%monc_communicator, ierr)

    i0 = current_state%local_grid%start(X_INDEX)
    j0 = current_state%local_grid%start(Y_INDEX)
    xhalo = current_state%local_grid%halo_size(X_INDEX)
    yhalo = current_state%local_grid%halo_size(Y_INDEX)
    nxg = current_state%global_grid%size(X_INDEX)
    nyg = current_state%global_grid%size(Y_INDEX)
    nnj = current_state%local_grid%size(Y_INDEX)
    nni = current_state%local_grid%size(X_INDEX)
    y_size=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX) * 2
    x_size=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX) * 2
    dx = current_state%global_grid%configuration%horizontal%dx
    dy = current_state%global_grid%configuration%horizontal%dy
    lwrite = .false. ! diagnostic write switch

    u_filter_enabled = .false.
    v_filter_enabled = .false.
    w_filter_enabled = .false.
    th_filter_enabled = .false.
    vert_filter_enabled = .false.
    u_filter_enabled=options_get_logical(current_state%options_database, "u_filter_enabled")
    v_filter_enabled=options_get_logical(current_state%options_database, "v_filter_enabled")
    w_filter_enabled=options_get_logical(current_state%options_database, "w_filter_enabled")
    th_filter_enabled=options_get_logical(current_state%options_database, "th_filter_enabled")
    if(.not.enable_theta)th_filter_enabled=.false.
    vert_filter_enabled=options_get_logical(current_state%options_database, "vert_filter_enabled")
    if (u_filter_enabled .or. v_filter_enabled .or. w_filter_enabled .or. th_filter_enabled .or. vert_filter_enabled) then
      allocate(current_state%immersed%ibprox(current_state%local_grid%size(Z_INDEX), y_size, x_size))
      ib_filter_prox=options_get_real(current_state%options_database, "f_prox")
    end if

    current_state%immersed%ib_type=options_get_integer(current_state%options_database, "ib_type")

    ! diagnostic options
    current_state%immersed%diags_enabled = options_get_logical(current_state%options_database, "ib_diags_enabled")
    if (current_state%immersed%diags_enabled) then
      current_state%immersed%dump_u = options_get_logical(current_state%options_database, "ib_dump_u")
      current_state%immersed%dump_v=options_get_logical(current_state%options_database, "ib_dump_v")
      current_state%immersed%dump_w=options_get_logical(current_state%options_database, "ib_dump_w")
      current_state%immersed%dump_th=options_get_logical(current_state%options_database, "ib_dump_th")
      if(.not.enable_theta)current_state%immersed%dump_th=.false.
      current_state%immersed%dump_p=options_get_logical(current_state%options_database, "ib_dump_p")
      current_state%immersed%dump_ustar=options_get_logical(current_state%options_database, "ib_dump_ustar")
      current_state%immersed%dump_freq=options_get_real(current_state%options_database, "diag_write_freq")
    end if


    allocate(nodes(nnodes,3))
    allocate(elems(nnodes,4))
    allocate(norms(nnodes,3))
    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))
    allocate(elem_lst(5000))
    allocate(elem_cent(nnodes,3))
    allocate(elem_xlims(nnodes,2))
    allocate(elem_ylims(nnodes,2))
    allocate(elem_zlims(nnodes,2))
    allocate(tmp_ijk_u(nni*nnj*5,3)) ! this may need to be increased
    allocate(tmp_ijk_v(nni*nnj*5,3)) ! this may need to be increased
    allocate(tmp_ijk_w(nni*nnj*5,3)) ! this may need to be increased
    allocate(tmp_ijk_s(nni*nnj*5,3)) ! this may need to be increased
    allocate(diag_xyz(nni*nnj*5,3))
    allocate(diag2_xyz(nni*nnj*5,3))
    allocate(diag3_xyz(nni*nnj*5,3))
    allocate(current_state%immersed%ib_col(y_size, x_size))
    allocate(current_state%immersed%ib_col_s(y_size, x_size))
    allocate(current_state%immersed%ib_col_u(y_size, x_size))
    allocate(current_state%immersed%ib_col_v(y_size, x_size))
    allocate(current_state%immersed%ib_col_w(y_size, x_size))
    allocate(current_state%immersed%indic_u(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%indic_v(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%indic_w(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%indic_s(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%ghost_u(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%ghost_v(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%ghost_w(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%ghost_s(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%rneutml_sq(current_state%local_grid%size(Z_INDEX), y_size, x_size))
    allocate(current_state%immersed%kmax_ji(y_size, x_size))

    if (myproc == 0) then 
      
      call check_status(nf90_open(path = trim(input_file), mode = nf90_nowrite, ncid = ncid))
      if (log_get_logging_level() .ge. LOG_DEBUG) then
        call log_master_log(LOG_DEBUG, "Reading in IB mesh data from:"//trim(input_file) )
      end if
      call read_variables(trim(input_file), ncid, nodes, elems, norms, &
                          NODES_KEY, ELEMS_KEY, NORMS_KEY)
      call check_status(nf90_close(ncid))

    end if ! proc 0 only

    nfloats = nnodes*3
    call mpi_bcast(nodes, nfloats, PRECISION_TYPE, 0, &
                   current_state%parallel%monc_communicator, ierr)
    call mpi_bcast(norms, nfloats, PRECISION_TYPE, 0, &
                   current_state%parallel%monc_communicator, ierr)
    nfloats = nnodes*4
    call mpi_bcast(elems, nfloats, PRECISION_TYPE, 0, &
                   current_state%parallel%monc_communicator, ierr)
    

    
    do ielem=1,nnodes
      elem_xlims(ielem,1) = min(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
                                nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      elem_xlims(ielem,2) = max(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
                                nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      elem_ylims(ielem,1) = min(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
                                nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      elem_ylims(ielem,2) = max(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
                                nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      elem_zlims(ielem,1) = min(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
                                nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      elem_zlims(ielem,2) = max(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
                                nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      tmp2 = elem_xlims(ielem,2)-elem_xlims(ielem,1)
      tmp3 = elem_ylims(ielem,2)-elem_ylims(ielem,1)
      if (tmp2 .gt. nxg*dx/2.0) tmp2 = tmp2 - (nxg*dx) 
      if (tmp2 .lt. -nxg*dx/2.0) tmp2 = tmp2 + (nxg*dx) 
      if (tmp3 .gt. nyg*dy/2.0) tmp3 = tmp3 - (nyg*dy) 
      if (tmp3 .lt. -nyg*dy/2.0) tmp3 = tmp3 + (nyg*dy) 
      elem_cent(ielem,1) = nodes(elems(ielem,1),1) + tmp2/2.
      elem_cent(ielem,2) = nodes(elems(ielem,1),2) + tmp3/2.
      elem_cent(ielem,3) = (nodes(elems(ielem,1),3) + nodes(elems(ielem,2),3) + &
                            nodes(elems(ielem,3),3) + nodes(elems(ielem,4),3) )/4.
    end do
    maxh = maxval(nodes(:,3))
    do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX) !run through local column
      if (current_state%global_grid%configuration%vertical%zn(k) .gt. maxh) then
        kmax=k
        exit
      end if
    end do

 
     current_state%immersed%indic_u(:,:,:) = 0 ! (0=fluid, 1=solid) 
     current_state%immersed%indic_v(:,:,:) = 0 ! (0=fluid, 1=solid) 
     current_state%immersed%indic_w(:,:,:) = 0 ! (0=fluid, 1=solid) 
     current_state%immersed%indic_s(:,:,:) = 0 ! (0=fluid, 1=solid) 
     current_state%immersed%rneutml_sq(:,:,:) = 0.0_DEFAULT_PRECISION 
     if (allocated(current_state%immersed%ibprox)) then
       current_state%immersed%ibprox(:,:,:) = 0.0_DEFAULT_PRECISION 
     end if
     nsurf_u = 0
     nsurf_v = 0
     nsurf_w = 0
     nsurf_s = 0
     current_state%immersed%kmax_ji(:,:) = 0


     do i=1, x_size
       do j=1, y_size
         jj = j+j0-1-yhalo
         ii = i+i0-1-xhalo
         if (jj .lt. 1) jj = jj + current_state%global_grid%size(Y_INDEX)
         if (ii .lt. 1) ii = ii + current_state%global_grid%size(X_INDEX)
         if (jj .gt. current_state%global_grid%size(Y_INDEX)) jj = jj - current_state%global_grid%size(Y_INDEX)
         if (ii .gt. current_state%global_grid%size(X_INDEX)) ii = ii - current_state%global_grid%size(X_INDEX)
 
         x_ug = (ii-1)*dx + (0.5_DEFAULT_PRECISION*dx)  ! x position of u 
         y_ug = (jj-1)*dy  ! y position of u
         x_vg = (ii-1)*dx  ! x position of v 
         y_vg = (jj-1)*dy + (0.5_DEFAULT_PRECISION*dy)  ! y position of v
         x_wg = (ii-1)*dx  ! x position of w 
         y_wg = (jj-1)*dy  ! y position of w
         x_sg = (ii-1)*dx  ! x position of scalars
         y_sg = (jj-1)*dy  ! y position of scalars


 
 ! limit max k value based on max elem_cent(ielem,3)
 ! localise element search based on horizontal minimum
         ! find list of all surface elements within 2.0dx of nearest grid point
         tmp1 = 5.0*max(dx,dy)
         nelem = 0
         do ielem=1,nnodes
           tmp2 = elem_cent(ielem,1) - x_sg
           tmp3 = elem_cent(ielem,2) - y_sg
           tmp4 = sqrt(tmp2*tmp2 + tmp3*tmp3)
           if ((tmp4 .le. tmp1)) then
             nelem = nelem+1
             elem_lst(nelem) = ielem
           end if
         end do
 
         do k=1, current_state%local_grid%local_domain_end_index(Z_INDEX)
             z_wg = current_state%global_grid%configuration%vertical%z(k)
           if (k .gt. kmax+1)then
             current_state%immersed%indic_u(k,j,i) = 0
             current_state%immersed%indic_v(k,j,i) = 0
             current_state%immersed%indic_w(k,j,i) = 0
             current_state%immersed%indic_s(k,j,i) = 0
           else
             !kmax !run through local column
             z_ug = current_state%global_grid%configuration%vertical%zn(k)
             z_vg = current_state%global_grid%configuration%vertical%zn(k)
             z_sg = current_state%global_grid%configuration%vertical%zn(k)
             current_state%immersed%indic_u(k,j,i) = -1
             current_state%immersed%indic_v(k,j,i) = -1
             current_state%immersed%indic_w(k,j,i) = -1
             current_state%immersed%indic_s(k,j,i) = -1
  
             ! run through surface elements and find nearest
             do ii=1,nelem
               ielem=elem_lst(ii)
  
               ! u points
               if (current_state%immersed%indic_u(k,j,i) .eq. -1)then 
                 x_g = x_ug
                 y_g = y_ug
                 z_g = z_ug
                 xl = elem_xlims(ielem,1)
                 xr = elem_xlims(ielem,2)
                 yl = elem_ylims(ielem,1)
                 yr = elem_ylims(ielem,2)
                 if ((x_g - xl) .gt. nxg*dx/2.0) xl = xl + (nxg*dx) 
                 if ((x_g - xl) .lt. -nxg*dx/2.0) xl = xl - (nxg*dx) 
                 if ((x_g - xr) .gt. nxg*dx/2.0) xr = xr + (nxg*dx) 
                 if ((x_g - xr) .lt. -nxg*dx/2.0) xr = xr - (nxg*dx) 
                 if ((y_g - yl) .gt. nyg*dy/2.0) yl = yl + (nyg*dy) 
                 if ((y_g - yl) .lt. -nyg*dy/2.0) yl = yl - (nyg*dy) 
                 if ((y_g - yr) .gt. nyg*dy/2.0) yr = yr + (nyg*dy) 
                 if ((y_g - yr) .lt. -nyg*dy/2.0) yr = yr - (nyg*dy) 
  
                 if (norms(ielem,3) .gt. 0.0 .and. x_g .ge. xl .and. x_g .le. xr .and. y_g .ge. yl .and. y_g .le. yr ) then
                   i00=0
                   i01=0
                   i10=0
                   i11=0
                   do inode=1,4 
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i00 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i01 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i10 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i11 = inode
                   end do
                   tmp1 = (x_g - xl)/(xr-xl)
                   tmp2 = (y_g - yl)/(yr-yl)
                   tmp3 = nodes(elems(ielem,i00),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i10),3)
                   tmp4 = nodes(elems(ielem,i01),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i11),3)
                   tmp5 = tmp3*(1.0-tmp2) + tmp4*tmp2
                   if (tmp5 .lt. z_g)current_state%immersed%indic_u(k,j,i) = 0
                   if (tmp5 .ge. z_g)then
                     current_state%immersed%indic_u(k,j,i) = 1
                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
                   end if
!                   if (abs(tmp5 - z_g) .lt. projtol)then
!                     current_state%immersed%indic_u(k,j,i) = 1
!                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
!                   end if
                 end if
               end if
  
               ! v points
               if (current_state%immersed%indic_v(k,j,i) .eq. -1)then 
                 x_g = x_vg
                 y_g = y_vg
                 z_g = z_vg
                 xl = elem_xlims(ielem,1)
                 xr = elem_xlims(ielem,2)
                 yl = elem_ylims(ielem,1)
                 yr = elem_ylims(ielem,2)
                 if ((x_g - xl) .gt. nxg*dx/2.0) xl = xl + (nxg*dx) 
                 if ((x_g - xl) .lt. -nxg*dx/2.0) xl = xl - (nxg*dx) 
                 if ((x_g - xr) .gt. nxg*dx/2.0) xr = xr + (nxg*dx) 
                 if ((x_g - xr) .lt. -nxg*dx/2.0) xr = xr - (nxg*dx) 
                 if ((y_g - yl) .gt. nyg*dy/2.0) yl = yl + (nyg*dy) 
                 if ((y_g - yl) .lt. -nyg*dy/2.0) yl = yl - (nyg*dy) 
                 if ((y_g - yr) .gt. nyg*dy/2.0) yr = yr + (nyg*dy) 
                 if ((y_g - yr) .lt. -nyg*dy/2.0) yr = yr - (nyg*dy) 
  
  
                 if (norms(ielem,3) .gt. 0.0 .and. x_g .ge. xl .and. x_g .le. xr .and. y_g .ge. yl .and. y_g .le. yr ) then
                   i00=0
                   i01=0
                   i10=0
                   i11=0
                   do inode=1,4 
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i00 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i01 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i10 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i11 = inode
                   end do
  
                   tmp1 = (x_g - xl)/(xr-xl)
                   tmp2 = (y_g - yl)/(yr-yl)
                   tmp3 = nodes(elems(ielem,i00),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i10),3)
                   tmp4 = nodes(elems(ielem,i01),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i11),3)
                   tmp5 = tmp3*(1.0-tmp2) + tmp4*tmp2
                   if (tmp5 .lt. z_g)current_state%immersed%indic_v(k,j,i) = 0
                   if (tmp5 .ge. z_g)then
                     current_state%immersed%indic_v(k,j,i) = 1
                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
                   end if
!                   if (abs(tmp5 - z_g) .lt. projtol)then
!                     current_state%immersed%indic_v(k,j,i) = 1
!                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
!                   end if
                 end if
               end if
  
               ! w points
               if (current_state%immersed%indic_w(k,j,i) .eq. -1)then 
                 x_g = x_wg
                 y_g = y_wg
                 z_g = z_wg
                 xl = elem_xlims(ielem,1)
                 xr = elem_xlims(ielem,2)
                 yl = elem_ylims(ielem,1)
                 yr = elem_ylims(ielem,2)
                 if ((x_g - xl) .gt. nxg*dx/2.0) xl = xl + (nxg*dx) 
                 if ((x_g - xl) .lt. -nxg*dx/2.0) xl = xl - (nxg*dx) 
                 if ((x_g - xr) .gt. nxg*dx/2.0) xr = xr + (nxg*dx) 
                 if ((x_g - xr) .lt. -nxg*dx/2.0) xr = xr - (nxg*dx) 
                 if ((y_g - yl) .gt. nyg*dy/2.0) yl = yl + (nyg*dy) 
                 if ((y_g - yl) .lt. -nyg*dy/2.0) yl = yl - (nyg*dy) 
                 if ((y_g - yr) .gt. nyg*dy/2.0) yr = yr + (nyg*dy) 
                 if ((y_g - yr) .lt. -nyg*dy/2.0) yr = yr - (nyg*dy) 
  
                 if (norms(ielem,3) .gt. 0.0 .and. x_g .ge. xl .and. x_g .le. xr .and. y_g .ge. yl .and. y_g .le. yr ) then
                   i00=0
                   i01=0
                   i10=0
                   i11=0
                   do inode=1,4 
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i00 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i01 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i10 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i11 = inode
                   end do
                   tmp1 = (x_g - xl)/(xr-xl)
                   tmp2 = (y_g - yl)/(yr-yl)
                   tmp3 = nodes(elems(ielem,i00),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i10),3)
                   tmp4 = nodes(elems(ielem,i01),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i11),3)
                   tmp5 = tmp3*(1.0-tmp2) + tmp4*tmp2
                   if (tmp5 .lt. z_g)current_state%immersed%indic_w(k,j,i) = 0
                   if (tmp5 .ge. z_g)then
                     current_state%immersed%indic_w(k,j,i) = 1
                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
                   end if
!                   if (abs(tmp5 - z_g) .lt. projtol)then
!                     current_state%immersed%indic_w(k,j,i) = 1
!                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
!                   end if
                 end if
               end if
  
               ! s points
               if (current_state%immersed%indic_s(k,j,i) .eq. -1)then 
                 x_g = x_sg
                 y_g = y_sg
                 z_g = z_sg
                 xl = elem_xlims(ielem,1)
                 xr = elem_xlims(ielem,2)
                 yl = elem_ylims(ielem,1)
                 yr = elem_ylims(ielem,2)
                 if ((x_g - xl) .gt. nxg*dx/2.0) xl = xl + (nxg*dx) 
                 if ((x_g - xl) .lt. -nxg*dx/2.0) xl = xl - (nxg*dx) 
                 if ((x_g - xr) .gt. nxg*dx/2.0) xr = xr + (nxg*dx) 
                 if ((x_g - xr) .lt. -nxg*dx/2.0) xr = xr - (nxg*dx) 
                 if ((y_g - yl) .gt. nyg*dy/2.0) yl = yl + (nyg*dy) 
                 if ((y_g - yl) .lt. -nyg*dy/2.0) yl = yl - (nyg*dy) 
                 if ((y_g - yr) .gt. nyg*dy/2.0) yr = yr + (nyg*dy) 
                 if ((y_g - yr) .lt. -nyg*dy/2.0) yr = yr - (nyg*dy) 
  
                 if (norms(ielem,3) .gt. 0.0 .and. x_g .ge. xl .and. x_g .le. xr .and. y_g .ge. yl .and. y_g .le. yr ) then
                   i00=0
                   i01=0
                   i10=0
                   i11=0
                   do inode=1,4 
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i00 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,1)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i01 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,1)) .lt. 1.e-10) )i10 = inode
                     if ((abs(nodes(elems(ielem,inode),1) - elem_xlims(ielem,2)) .lt. 1.e-10) .and.&
                         (abs(nodes(elems(ielem,inode),2) - elem_ylims(ielem,2)) .lt. 1.e-10) )i11 = inode
                   end do
                   tmp1 = (x_g - xl)/(xr-xl)
                   tmp2 = (y_g - yl)/(yr-yl)
                   tmp3 = nodes(elems(ielem,i00),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i10),3)
                   tmp4 = nodes(elems(ielem,i01),3)*(1.0-tmp1) + tmp1*nodes(elems(ielem,i11),3)
                   tmp5 = tmp3*(1.0-tmp2) + tmp4*tmp2
                   if (tmp5 .lt. z_g)current_state%immersed%indic_s(k,j,i) = 0
                   if (tmp5 .ge. z_g)then
                     current_state%immersed%indic_s(k,j,i) = 1
                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
                   end if
!                   if (abs(tmp5 - z_g) .lt. projtol)then
!                     current_state%immersed%indic_s(k,j,i) = 1
!                     if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
!                   end if
                 end if
               end if
  
    
             end do ! ielem loop
           end if ! if k .lt. kmax


           ! calculate rneut_ml for all points in column
           current_state%immersed%rneutml_sq(k,j,i) = 0.0_DEFAULT_PRECISION
           if (current_state%immersed%indic_w(k,j,i) .eq. 0) then
             tmp1w = 1.e10
             do ii=1,nelem
               ielem=elem_lst(ii)
               tmp3 = elem_cent(ielem,1) - x_wg
               tmp4 = elem_cent(ielem,2) - y_wg
               tmp5 = elem_cent(ielem,3) - z_wg
               tmp2 = sqrt(tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5)
               if (tmp2 .lt. tmp1w) then
                 tmp1w = tmp2
               end if
             end do

             if (abs(tmp1w) .lt. 4.0*dx)then ! refine search if point is close to surface
               zgrid=current_state%global_grid%configuration%vertical%z(:)
               dz = zgrid(k+1)-zgrid(k)
!               idiag=-1
               call find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_wg, y_wg, z_wg, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen, ip_dir)
               tmp2 = sqrt((x_bi-x_wg)*(x_bi-x_wg) + (y_bi-y_wg)*(y_bi-y_wg) + (z_bi-z_wg)*(z_bi-z_wg))
               if (tmp2 .lt. tmp1w)tmp1w=tmp2
             end if
             current_state%immersed%rneutml_sq(k,j,i) = tmp1w ! store distance for ml calculation

           end if ! if fluid w point


         end do ! k loop
       end do ! j loop
     end do ! i loop
 

 
    ! find ghost points or VR points
     current_state%immersed%ghost_u = 0
     current_state%immersed%ghost_v = 0
     current_state%immersed%ghost_w = 0
     current_state%immersed%ghost_s = 0
     nghosts_u = 0
     nghosts_v = 0
     nghosts_w = 0
     nghosts_s = 0
     current_state%immersed%ib_col = .false.
     current_state%immersed%ib_col_s = .false.
     current_state%immersed%ib_col_u = .false.
     current_state%immersed%ib_col_v = .false.
     current_state%immersed%ib_col_w = .false.
     do i=1+xhalo,current_state%local_grid%size(X_INDEX)+xhalo
       do j=1+yhalo,current_state%local_grid%size(Y_INDEX)+yhalo
         do k=2, current_state%local_grid%size(Z_INDEX)-1 !run through local column
           
           ! No slip -----------------------
           if (current_state%immersed%ib_type==0)then
             ! u points
             if (    current_state%immersed%indic_u(k,j,i)==1 )then
               if (sum(current_state%immersed%indic_u(k-1:k+1,j-1:j+1,i-1:i+1)) .lt. 27) then
                 nghosts_u = nghosts_u + 1
                 current_state%immersed%ghost_u(k,j,i) = nghosts_u ! ghost point
                 tmp_ijk_u(nghosts_u,1) = i
                 tmp_ijk_u(nghosts_u,2) = j
                 tmp_ijk_u(nghosts_u,3) = k
                 current_state%immersed%ib_col_u(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! v points
             if (    current_state%immersed%indic_v(k,j,i)==1 )then
               if (sum(current_state%immersed%indic_v(k-1:k+1,j-1:j+1,i-1:i+1)) .lt. 27) then
                 nghosts_v = nghosts_v + 1
                 current_state%immersed%ghost_v(k,j,i) = nghosts_v ! ghost point
                 tmp_ijk_v(nghosts_v,1) = i
                 tmp_ijk_v(nghosts_v,2) = j
                 tmp_ijk_v(nghosts_v,3) = k
                 current_state%immersed%ib_col_v(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! w points
             if (    current_state%immersed%indic_w(k,j,i)==1 )then
               if (sum(current_state%immersed%indic_w(k-1:k+1,j-1:j+1,i-1:i+1)) .lt. 27) then
                 nghosts_w = nghosts_w + 1
                 current_state%immersed%ghost_w(k,j,i) = nghosts_w ! ghost point
                 tmp_ijk_w(nghosts_w,1) = i
                 tmp_ijk_w(nghosts_w,2) = j
                 tmp_ijk_w(nghosts_w,3) = k
                 current_state%immersed%ib_col_w(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! scalar points
             if (    current_state%immersed%indic_s(k,j,i)==1 )then
               if (sum(current_state%immersed%indic_s(k-1:k+1,j-1:j+1,i-1:i+1)) .lt. 27) then
                 nghosts_s = nghosts_s + 1
                 current_state%immersed%ghost_s(k,j,i) = nghosts_s ! ghost point
                 tmp_ijk_s(nghosts_s,1) = i
                 tmp_ijk_s(nghosts_s,2) = j
                 tmp_ijk_s(nghosts_s,3) = k
                 current_state%immersed%ib_col_s(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
           end if


           ! VR -----------------------
           if (current_state%immersed%ib_type==1)then
             ! u points
             if (    current_state%immersed%indic_u(k,j,i)==0 )then
               if (sum(current_state%immersed%indic_u(k-1:k+1,j-1:j+1,i-1:i+1)) .gt. 0) then
                 nghosts_u = nghosts_u + 1
                 current_state%immersed%ghost_u(k,j,i) = nghosts_u ! ghost point
                 tmp_ijk_u(nghosts_u,1) = i
                 tmp_ijk_u(nghosts_u,2) = j
                 tmp_ijk_u(nghosts_u,3) = k
                 current_state%immersed%ib_col_u(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! v points
             if (    current_state%immersed%indic_v(k,j,i)==0 )then
               if (sum(current_state%immersed%indic_v(k-1:k+1,j-1:j+1,i-1:i+1)) .gt. 0) then
                 nghosts_v = nghosts_v + 1
                 current_state%immersed%ghost_v(k,j,i) = nghosts_v ! ghost point
                 tmp_ijk_v(nghosts_v,1) = i
                 tmp_ijk_v(nghosts_v,2) = j
                 tmp_ijk_v(nghosts_v,3) = k
                 current_state%immersed%ib_col_v(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! w points
             if (    current_state%immersed%indic_w(k,j,i)==0 )then
               if (sum(current_state%immersed%indic_w(k-1:k+1,j-1:j+1,i-1:i+1)) .gt. 0) then
                 nghosts_w = nghosts_w + 1
                 current_state%immersed%ghost_w(k,j,i) = nghosts_w ! ghost point
                 tmp_ijk_w(nghosts_w,1) = i
                 tmp_ijk_w(nghosts_w,2) = j
                 tmp_ijk_w(nghosts_w,3) = k
                 current_state%immersed%ib_col_w(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
  
             ! scalar points
             if (    current_state%immersed%indic_s(k,j,i)==0 )then
               if (sum(current_state%immersed%indic_s(k-1:k+1,j-1:j+1,i-1:i+1)) .gt. 0) then
                 nghosts_s = nghosts_s + 1
                 current_state%immersed%ghost_s(k,j,i) = nghosts_s ! ghost point
                 tmp_ijk_s(nghosts_s,1) = i
                 tmp_ijk_s(nghosts_s,2) = j
                 tmp_ijk_s(nghosts_s,3) = k
                 current_state%immersed%ib_col_s(j,i) = .true.
                 current_state%immersed%ib_col(j,i) = .true.
               end if
             end if
           end if


         end do
       end do
     end do


    ! allocate common arrays (NS and VR)
    allocate(current_state%immersed%gp_ijk_u(nghosts_u,3))
    allocate(current_state%immersed%gp_ijk_v(nghosts_v,3))
    allocate(current_state%immersed%gp_ijk_w(nghosts_w,3))
    allocate(current_state%immersed%gp_ijk_s(nghosts_s,3))
    allocate(current_state%immersed%norms_u(nghosts_u,3))
    allocate(current_state%immersed%norms_v(nghosts_v,3))
    allocate(current_state%immersed%norms_w(nghosts_w,3))
    allocate(current_state%immersed%norms_s(nghosts_s,3))
    allocate(current_state%immersed%idiag_u(nghosts_u))
    allocate(current_state%immersed%idiag_v(nghosts_v))
    allocate(current_state%immersed%idiag_w(nghosts_w))
    allocate(current_state%immersed%idiag_s(nghosts_s))
    allocate(current_state%immersed%gp_store_u(nghosts_u))
    allocate(current_state%immersed%gp_store_v(nghosts_v))
    allocate(current_state%immersed%gp_store_w(nghosts_w))
    if (current_state%use_viscosity_and_diffusion)then
      allocate(current_state%immersed%gp_store_vis(nghosts_w))
      allocate(current_state%immersed%gp_store_diff(nghosts_w))
    end if
    if(enable_theta)then
      allocate(current_state%immersed%gp_store_th(nghosts_s))
    end if
    if(current_state%number_q_fields .gt. 0)then
      allocate(current_state%immersed%gp_store_q(nghosts_s))
    end if
 
    ! local arrays
    allocate(gp_xyz_u(nghosts_u,3))
    allocate(ip_xyz_u(nghosts_u,3))
    allocate(bi_xyz_u(nghosts_u,3))
    allocate(gp_xyz_v(nghosts_v,3))
    allocate(ip_xyz_v(nghosts_v,3))
    allocate(bi_xyz_v(nghosts_v,3))
    allocate(gp_xyz_w(nghosts_w,3))
    allocate(ip_xyz_w(nghosts_w,3))
    allocate(bi_xyz_w(nghosts_w,3))
    allocate(diag4_xyz_w(nghosts_w,3))
    allocate(gp_xyz_s(nghosts_s,3))
    allocate(ip_xyz_s(nghosts_s,3))
    allocate(bi_xyz_s(nghosts_s,3))

    ! diagnostic/control arrays
    if (current_state%immersed%diags_enabled) then
      if(current_state%immersed%dump_u.and.(.not.allocated(current_state%immersed%ib_u)))then
        allocate(current_state%immersed%ib_u(nghosts_u))
      end if
      if(current_state%immersed%dump_v.and.(.not.allocated(current_state%immersed%ib_v)))then
        allocate(current_state%immersed%ib_v(nghosts_v))
      end if
      if(current_state%immersed%dump_w.and.(.not.allocated(current_state%immersed%ib_w)))then
        allocate(current_state%immersed%ib_w(nghosts_w))
      end if
      if(current_state%immersed%dump_th.and.(.not.allocated(current_state%immersed%ib_th)))then
        allocate(current_state%immersed%ib_th(nghosts_s))
      end if
      if(current_state%immersed%dump_p.and.(.not.allocated(current_state%immersed%ib_p)))then
        allocate(current_state%immersed%ib_p(nghosts_s))
      end if
      if(current_state%immersed%dump_ustar.and.(.not.allocated(current_state%immersed%ib_ustar)))then
        allocate(current_state%immersed%ib_ustar(nghosts_w))
      end if
    end if

    if(enable_theta)then
      allocate(current_state%immersed%theta_surf_ib_s(nghosts_s))
      allocate(current_state%immersed%thref_surf_ib_s(nghosts_s))
      allocate(current_state%immersed%theta_surf_ib_w(nghosts_w))
      allocate(current_state%immersed%thref_surf_ib_w(nghosts_w))
    end if

    if (use_surface_bcs.and.surface_boundary_type.eq. PRESCRIBED_SURFACE_FLUXES)then
      allocate(current_state%immersed%w2s_idw_idx(nghosts_s,idw_nmax))
      allocate(current_state%immersed%w2s_idw_wgt(nghosts_s,idw_nmax))
    end if

    ! allocate NS arrays
    if (current_state%immersed%ib_type==0)then
      allocate(current_state%immersed%ip_ijk_uu(nghosts_u,3))
      allocate(current_state%immersed%ip_ijk_vv(nghosts_v,3))
      allocate(current_state%immersed%ip_ijk_ww(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_k2e_w(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_k2e_ws(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_k2e_wu(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_k2e_wv(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_ss(nghosts_s,3))
      allocate(current_state%immersed%int_xyz_uu(nghosts_u,3))
      allocate(current_state%immersed%int_xyz_vv(nghosts_v,3))
      allocate(current_state%immersed%int_xyz_ww(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_k2e_w(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_k2e_ws(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_k2e_wu(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_k2e_wv(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_ss(nghosts_s,3))
      
      allocate(current_state%immersed%bi_ijk_u(nghosts_u,3))
      allocate(current_state%immersed%bi_ijk_v(nghosts_v,3))
      allocate(current_state%immersed%bi_ijk_w(nghosts_w,3))
      allocate(current_state%immersed%bi_ijk_s(nghosts_s,3))
      allocate(current_state%immersed%bi_xyz_u(nghosts_u,3))
      allocate(current_state%immersed%bi_xyz_v(nghosts_v,3))
      allocate(current_state%immersed%bi_xyz_w(nghosts_w,3))
      allocate(current_state%immersed%bi_xyz_s(nghosts_s,3))
      allocate(current_state%immersed%k2e_xyz_w(nghosts_w,3))
      allocate(current_state%immersed%biint_xyz_u(nghosts_u,3))
      allocate(current_state%immersed%biint_xyz_v(nghosts_v,3))
      allocate(current_state%immersed%biint_xyz_w(nghosts_w,3))
      allocate(current_state%immersed%biint_xyz_s(nghosts_s,3))
      allocate(current_state%immersed%gp_len_w(nghosts_w))
      allocate(current_state%immersed%gp_len_s(nghosts_s))
      allocate(current_state%immersed%k2e_len(nghosts_w))
    end if !NS arrays

    
    ! allocate VR arrays
    if (current_state%immersed%ib_type==1)then
      allocate(current_state%immersed%ip_ijk_uu(nghosts_u,3))
      allocate(current_state%immersed%ip_ijk_uv(nghosts_u,3))
      allocate(current_state%immersed%ip_ijk_uw(nghosts_u,3))
      allocate(current_state%immersed%ip_ijk_vu(nghosts_v,3))
      allocate(current_state%immersed%ip_ijk_vv(nghosts_v,3))
      allocate(current_state%immersed%ip_ijk_vw(nghosts_v,3))
      allocate(current_state%immersed%ip_ijk_wu(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_wv(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_ww(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_ws(nghosts_w,3))
      allocate(current_state%immersed%ip_ijk_ss(nghosts_s,3))
      
      allocate(current_state%immersed%int_xyz_uu(nghosts_u,3))
      allocate(current_state%immersed%int_xyz_uv(nghosts_u,3))
      allocate(current_state%immersed%int_xyz_uw(nghosts_u,3))
      allocate(current_state%immersed%int_xyz_vu(nghosts_v,3))
      allocate(current_state%immersed%int_xyz_vv(nghosts_v,3))
      allocate(current_state%immersed%int_xyz_vw(nghosts_v,3))
      allocate(current_state%immersed%int_xyz_wu(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_wv(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_ww(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_ws(nghosts_w,3))
      allocate(current_state%immersed%int_xyz_ss(nghosts_s,3))
      if (surface_boundary_type.eq.PRESCRIBED_SURFACE_FLUXES)then
        allocate(current_state%immersed%int_xyz_sw(nghosts_s,3))
        allocate(current_state%immersed%ip_ijk_sw(nghosts_s,3))
      end if
      
      allocate(current_state%immersed%gp_len_u(nghosts_u))
      allocate(current_state%immersed%gp_len_v(nghosts_v))
      allocate(current_state%immersed%gp_len_w(nghosts_w))
      allocate(current_state%immersed%gp_len_s(nghosts_s))
      allocate(current_state%immersed%ip_len_u(nghosts_u))
      allocate(current_state%immersed%ip_len_v(nghosts_v))
      allocate(current_state%immersed%ip_len_w(nghosts_w))
      allocate(current_state%immersed%ip_len_s(nghosts_s))
      
      ! needed for specified_flux BC
      allocate(current_state%immersed%bi_xyz_w(nghosts_w,3))

    end if !VR arrays





    do ighost = 1,nghosts_u ! u ghosts
      current_state%immersed%gp_ijk_u(ighost,1) = tmp_ijk_u(ighost,1)
      current_state%immersed%gp_ijk_u(ighost,2) = tmp_ijk_u(ighost,2)
      current_state%immersed%gp_ijk_u(ighost,3) = tmp_ijk_u(ighost,3)
      i = current_state%immersed%gp_ijk_u(ighost,1)
      j = current_state%immersed%gp_ijk_u(ighost,2)
      k = current_state%immersed%gp_ijk_u(ighost,3)
      if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
      current_state%immersed%ghost_u(k,j,i) = ighost ! reset ghost index
      x_g = (i0 + i - 2 - xhalo)*dx + (0.5_DEFAULT_PRECISION*dx)  ! x position of u 
      y_g = (j0 + j - 2 - yhalo)*dy  ! y position of u
      z_g = current_state%global_grid%configuration%vertical%zn(k)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      ! find list of all surface elements within 2dx of nearest
      tmp1 = 2.*max(dx,dy)
      nelem = 0
      do ielem=1,nnodes
        tmp3 = elem_cent(ielem,1) - x_g
        tmp4 = elem_cent(ielem,2) - y_g
        tmp5 = elem_cent(ielem,3) - z_g
        if (tmp3 .gt. nxg*dx/2.0) tmp3 = tmp3 - (nxg*dx) 
        if (tmp3 .lt. -nxg*dx/2.0) tmp3 = tmp3 + (nxg*dx) 
        if (tmp4 .gt. nyg*dy/2.0) tmp4 = tmp4 - (nyg*dy) 
        if (tmp4 .lt. -nyg*dy/2.0) tmp4 = tmp4 + (nyg*dy) 
        tmp2 = sqrt(tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5)
        if (tmp2 .le. tmp1) then
          nelem = nelem+1
          elem_lst(nelem) = ielem
        end if
      end do

      ! no-slip ---------------------
      if (current_state%immersed%ib_type==0)then
  !      if (myproc .eq. 0 .and. (ighost .eq. 2817)) idiag=-1
        call find_image_point(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen)
        if(z_ip .lt. z_g .or. idiag .eq. 0)then
          write(*,*)'WARNING (gridmanager): BI not found for u proc, ijk', myproc, i, j, k 
          write(*,*)'WARNING (gridmanager): idiag, pos', idiag, x_g, y_g, z_g 
          write(*,*)'WARNING (gridmanager): ighost, image point location', ighost, x_ip, y_ip, z_ip
          current_state%immersed%ghost_u(k,j,i) = -1
        end if
  
        gp_xyz_u(ighost,1) = x_g
        gp_xyz_u(ighost,2) = y_g
        gp_xyz_u(ighost,3) = z_g
        ip_xyz_u(ighost,1) = x_ip
        ip_xyz_u(ighost,2) = y_ip
        ip_xyz_u(ighost,3) = z_ip
        bi_xyz_u(ighost,1) = x_bi
        bi_xyz_u(ighost,2) = y_bi
        bi_xyz_u(ighost,3) = z_bi
        current_state%immersed%idiag_u(ighost)=idiag
  
  
        ! find indices for image point interpolation
        ip0 = int((x_ip/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_ip/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%ip_ijk_uu(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_uu(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_uu(ighost,3) = kp0
        current_state%immersed%int_xyz_uu(ighost,1) = xdd
        current_state%immersed%int_xyz_uu(ighost,2) = ydd
        current_state%immersed%int_xyz_uu(ighost,3) = zdd
    
        ! same for BI run-time diagnostics
        ip0 = int((x_bi/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(y_bi/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_bi/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = y_bi/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_bi - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%bi_ijk_u(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%bi_ijk_u(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%bi_ijk_u(ighost,3) = kp0
        current_state%immersed%biint_xyz_u(ighost,1) = xdd
        current_state%immersed%biint_xyz_u(ighost,2) = ydd
        current_state%immersed%biint_xyz_u(ighost,3) = zdd
        current_state%immersed%bi_xyz_u(ighost,1) = x_bi
        current_state%immersed%bi_xyz_u(ighost,2) = y_bi
        current_state%immersed%bi_xyz_u(ighost,3) = z_bi
      end if ! if ib_type==0


      ! VR ---------------------
      if (current_state%immersed%ib_type==1)then
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        dz = zgrid(k+1)-zgrid(k)
!        if (myproc .eq. 0 .and. (ighost .eq. 1)) idiag=-1
        call find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen, ip_dir)
        if(z_ip .lt. z_g .or. idiag .eq. 0)then
          write(*,*)'WARNING (gridmanager): BI not found for u proc, ijk', myproc, i, j, k 
          write(*,*)'WARNING (gridmanager): idiag, pos', idiag, x_g, y_g, z_g 
          write(*,*)'WARNING (gridmanager): ighost, image point location', ighost, x_ip, y_ip, z_ip
          current_state%immersed%ghost_u(k,j,i) = -1
        end if
  
        gp_xyz_u(ighost,1) = x_g
        gp_xyz_u(ighost,2) = y_g
        gp_xyz_u(ighost,3) = z_g
        ip_xyz_u(ighost,1) = x_ip
        ip_xyz_u(ighost,2) = y_ip
        ip_xyz_u(ighost,3) = z_ip
        bi_xyz_u(ighost,1) = x_bi
        bi_xyz_u(ighost,2) = y_bi
        bi_xyz_u(ighost,3) = z_bi
        current_state%immersed%idiag_u(ighost)=idiag
  
        current_state%immersed%norms_u(ighost,1) = nrmi
        current_state%immersed%norms_u(ighost,2) = nrmj
        current_state%immersed%norms_u(ighost,3) = nrmk
        current_state%immersed%gp_len_u(ighost) = sqrt((x_g - x_bi)**2.0 +&
                                                       (y_g - y_bi)**2.0 +&
                                                       (z_g - z_bi)**2.0)
        current_state%immersed%ip_len_u(ighost) = sqrt((x_ip - x_bi)**2.0 +&
                                                       (y_ip - y_bi)**2.0 +&
                                                       (z_ip - z_bi)**2.0)

                                                     
        ! find indices for image point interpolation
        ! u grid (u ghosts) --------------------------------------
        ip0 = int((x_ip/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_uu(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_uu(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_uu(ighost,3) = kp0
        current_state%immersed%int_xyz_uu(ighost,1) = xdd
        current_state%immersed%int_xyz_uu(ighost,2) = ydd
        current_state%immersed%int_xyz_uu(ighost,3) = zdd
                                                     
        ! v grid (u ghosts) --------------------------------------
        ip0 = int(x_ip/dx) + 1
        jp0 = int((y_ip/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_uv(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_uv(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_uv(ighost,3) = kp0
        current_state%immersed%int_xyz_uv(ighost,1) = xdd
        current_state%immersed%int_xyz_uv(ighost,2) = ydd
        current_state%immersed%int_xyz_uv(ighost,3) = zdd
                                                     
        ! w grid (u ghosts) --------------------------------------
        zgrid=current_state%global_grid%configuration%vertical%z(:)
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_uw(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_uw(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_uw(ighost,3) = kp0
        current_state%immersed%int_xyz_uw(ighost,1) = xdd
        current_state%immersed%int_xyz_uw(ighost,2) = ydd
        current_state%immersed%int_xyz_uw(ighost,3) = zdd
                                                     
      end if ! if ib_type==1

    end do ! u ghosts

    !--------------------------------------------------------------------------------------


    do ighost = 1,nghosts_v ! v ghosts
      current_state%immersed%gp_ijk_v(ighost,1) = tmp_ijk_v(ighost,1)
      current_state%immersed%gp_ijk_v(ighost,2) = tmp_ijk_v(ighost,2)
      current_state%immersed%gp_ijk_v(ighost,3) = tmp_ijk_v(ighost,3)
      i = current_state%immersed%gp_ijk_v(ighost,1)
      j = current_state%immersed%gp_ijk_v(ighost,2)
      k = current_state%immersed%gp_ijk_v(ighost,3)
      if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
      current_state%immersed%ghost_v(k,j,i) = ighost ! reset ghost index
      x_g = (i0 + i - 2 - xhalo)*dx   ! x position of v
      y_g = (j0 + j - 2 - yhalo)*dy + (0.5_DEFAULT_PRECISION*dy) ! y position of v
      z_g = current_state%global_grid%configuration%vertical%zn(k)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)

      ! find list of all surface elements within 2dx of nearest
      tmp1 = 2.*max(dx,dy)
      nelem = 0
      do ielem=1,nnodes
        tmp3 = elem_cent(ielem,1) - x_g
        tmp4 = elem_cent(ielem,2) - y_g
        tmp5 = elem_cent(ielem,3) - z_g
        if (tmp3 .gt. nxg*dx/2.0) tmp3 = tmp3 - (nxg*dx) 
        if (tmp3 .lt. -nxg*dx/2.0) tmp3 = tmp3 + (nxg*dx) 
        if (tmp4 .gt. nyg*dy/2.0) tmp4 = tmp4 - (nyg*dy) 
        if (tmp4 .lt. -nyg*dy/2.0) tmp4 = tmp4 + (nyg*dy) 
        tmp2 = sqrt(tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5)
        if (tmp2 .le. tmp1) then
          nelem = nelem+1
          elem_lst(nelem) = ielem
        end if
      end do

      ! no-slip ---------------------
      if (current_state%immersed%ib_type==0)then
        call find_image_point(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen)
  
        gp_xyz_v(ighost,1) = x_g
        gp_xyz_v(ighost,2) = y_g
        gp_xyz_v(ighost,3) = z_g
        ip_xyz_v(ighost,1) = x_ip
        ip_xyz_v(ighost,2) = y_ip
        ip_xyz_v(ighost,3) = z_ip
        bi_xyz_v(ighost,1) = x_bi
        bi_xyz_v(ighost,2) = y_bi
        bi_xyz_v(ighost,3) = z_bi
        current_state%immersed%idiag_v(ighost)=idiag
  
        ! find indices for image point interpolation
        ip0 = int(x_ip/dx) + 1
        jp0 = int((y_ip/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%ip_ijk_vv(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_vv(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_vv(ighost,3) = kp0
        current_state%immersed%int_xyz_vv(ighost,1) = xdd
        current_state%immersed%int_xyz_vv(ighost,2) = ydd
        current_state%immersed%int_xyz_vv(ighost,3) = zdd
    
        ! same for BI run-time diagnostics
        ip0 = int(x_bi/dx) + 1
        jp0 = int((y_bi/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_bi/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_bi/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (z_bi - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%bi_ijk_v(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%bi_ijk_v(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%bi_ijk_v(ighost,3) = kp0
        current_state%immersed%biint_xyz_v(ighost,1) = xdd
        current_state%immersed%biint_xyz_v(ighost,2) = ydd
        current_state%immersed%biint_xyz_v(ighost,3) = zdd
        current_state%immersed%bi_xyz_v(ighost,1) = x_bi
        current_state%immersed%bi_xyz_v(ighost,2) = y_bi
        current_state%immersed%bi_xyz_v(ighost,3) = z_bi
      end if ! if ib_type==0
    

      ! VR ---------------------
      if (current_state%immersed%ib_type==1)then
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        dz = zgrid(k+1)-zgrid(k)
        call find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen, ip_dir)
        if(z_ip .lt. z_g .or. idiag .eq. 0)then
          write(*,*)'WARNING (gridmanager): BI not found for u proc, ijk', myproc, i, j, k 
          write(*,*)'WARNING (gridmanager): idiag, pos', idiag, x_g, y_g, z_g 
          write(*,*)'WARNING (gridmanager): ighost, image point location', ighost, x_ip, y_ip, z_ip
          current_state%immersed%ghost_v(k,j,i) = -1
        end if
  
        gp_xyz_v(ighost,1) = x_g
        gp_xyz_v(ighost,2) = y_g
        gp_xyz_v(ighost,3) = z_g
        ip_xyz_v(ighost,1) = x_ip
        ip_xyz_v(ighost,2) = y_ip
        ip_xyz_v(ighost,3) = z_ip
        bi_xyz_v(ighost,1) = x_bi
        bi_xyz_v(ighost,2) = y_bi
        bi_xyz_v(ighost,3) = z_bi
        current_state%immersed%idiag_v(ighost)=idiag
  
        current_state%immersed%norms_v(ighost,1) = nrmi
        current_state%immersed%norms_v(ighost,2) = nrmj
        current_state%immersed%norms_v(ighost,3) = nrmk
        current_state%immersed%gp_len_v(ighost) = sqrt((x_g - x_bi)**2.0 +&
                                                       (y_g - y_bi)**2.0 +&
                                                       (z_g - z_bi)**2.0)
        current_state%immersed%ip_len_v(ighost) = sqrt((x_ip - x_bi)**2.0 +&
                                                       (y_ip - y_bi)**2.0 +&
                                                       (z_ip - z_bi)**2.0)
  

        ! find indices for image point interpolation
        ! u grid (v ghosts) --------------------------------------
        ip0 = int((x_ip/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_vu(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_vu(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_vu(ighost,3) = kp0
        current_state%immersed%int_xyz_vu(ighost,1) = xdd
        current_state%immersed%int_xyz_vu(ighost,2) = ydd
        current_state%immersed%int_xyz_vu(ighost,3) = zdd
                                                     
        ! v grid (v ghosts) --------------------------------------
        ip0 = int(x_ip/dx) + 1
        jp0 = int((y_ip/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_vv(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_vv(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_vv(ighost,3) = kp0
        current_state%immersed%int_xyz_vv(ighost,1) = xdd
        current_state%immersed%int_xyz_vv(ighost,2) = ydd
        current_state%immersed%int_xyz_vv(ighost,3) = zdd
                                                     
        ! w grid (v ghosts) --------------------------------------
        zgrid=current_state%global_grid%configuration%vertical%z(:)
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_vw(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_vw(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_vw(ighost,3) = kp0
        current_state%immersed%int_xyz_vw(ighost,1) = xdd
        current_state%immersed%int_xyz_vw(ighost,2) = ydd
        current_state%immersed%int_xyz_vw(ighost,3) = zdd
                                                     
      end if ! if ib_type==1

    end do ! v ghosts

    !--------------------------------------------------------------------------------------




    do ighost = 1,nghosts_w ! w ghosts
      current_state%immersed%gp_ijk_w(ighost,1) = tmp_ijk_w(ighost,1)
      current_state%immersed%gp_ijk_w(ighost,2) = tmp_ijk_w(ighost,2)
      current_state%immersed%gp_ijk_w(ighost,3) = tmp_ijk_w(ighost,3)
      i = current_state%immersed%gp_ijk_w(ighost,1)
      j = current_state%immersed%gp_ijk_w(ighost,2)
      k = current_state%immersed%gp_ijk_w(ighost,3)
      if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
      current_state%immersed%ghost_w(k,j,i) = ighost ! reset ghost index
      x_g = (i0 + i - 2 - xhalo)*dx  ! x position of w
      y_g = (j0 + j - 2 - yhalo)*dy  ! y position of w
      z_g = current_state%global_grid%configuration%vertical%z(k)
      zgrid=current_state%global_grid%configuration%vertical%z(:)

      ! find list of all surface elements within 2dx of nearest
      tmp1 = 2.*max(dx,dy)
      nelem = 0
      do ielem=1,nnodes
        tmp3 = elem_cent(ielem,1) - x_g
        tmp4 = elem_cent(ielem,2) - y_g
        tmp5 = elem_cent(ielem,3) - z_g
        if (tmp3 .gt. nxg*dx/2.0) tmp3 = tmp3 - (nxg*dx) 
        if (tmp3 .lt. -nxg*dx/2.0) tmp3 = tmp3 + (nxg*dx) 
        if (tmp4 .gt. nyg*dy/2.0) tmp4 = tmp4 - (nyg*dy) 
        if (tmp4 .lt. -nyg*dy/2.0) tmp4 = tmp4 + (nyg*dy) 
        tmp2 = sqrt(tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5)
        if (tmp2 .le. tmp1) then
          nelem = nelem+1
          elem_lst(nelem) = ielem
        end if
      end do

      ! no-slip ---------------------
      if (current_state%immersed%ib_type==0)then
        call find_image_point(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen)
  
  
        gp_xyz_w(ighost,1) = x_g
        gp_xyz_w(ighost,2) = y_g
        gp_xyz_w(ighost,3) = z_g
        ip_xyz_w(ighost,1) = x_ip
        ip_xyz_w(ighost,2) = y_ip
        ip_xyz_w(ighost,3) = z_ip
        bi_xyz_w(ighost,1) = x_bi
        bi_xyz_w(ighost,2) = y_bi
        bi_xyz_w(ighost,3) = z_bi
        current_state%immersed%idiag_w(ighost)=idiag
  
        ! find indices for image point interpolation
        ! for w and scalar points this needs to be done for u,v, and w points as well
        current_state%immersed%norms_w(ighost,1) = nrmi
        current_state%immersed%norms_w(ighost,2) = nrmj
        current_state%immersed%norms_w(ighost,3) = nrmk
  
        ! w
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%ip_ijk_ww(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_ww(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_ww(ighost,3) = kp0
        current_state%immersed%int_xyz_ww(ighost,1) = xdd
        current_state%immersed%int_xyz_ww(ighost,2) = ydd
        current_state%immersed%int_xyz_ww(ighost,3) = zdd
  
        ! same for BI run-time diagnostics
        ip0 = int(x_bi/dx) + 1
        jp0 = int(y_bi/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_bi/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_bi/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_bi - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%bi_ijk_w(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%bi_ijk_w(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%bi_ijk_w(ighost,3) = kp0
        current_state%immersed%biint_xyz_w(ighost,1) = xdd
        current_state%immersed%biint_xyz_w(ighost,2) = ydd
        current_state%immersed%biint_xyz_w(ighost,3) = zdd
        current_state%immersed%bi_xyz_w(ighost,1) = x_bi
        current_state%immersed%bi_xyz_w(ighost,2) = y_bi
        current_state%immersed%bi_xyz_w(ighost,3) = z_bi
        current_state%immersed%gp_len_w(ighost) = llen
  
        ! k2 equivalents for u, v and w grids
        ! not less than 50z0 from surface to avoid roughness sublayer (Basu, BLM 163, 2017)
!        llen = max(current_state%global_grid%configuration%vertical%zn(2), z0*50.0)
        llen = current_state%global_grid%configuration%vertical%zn(2)
        current_state%immersed%k2e_len(ighost) = llen
        tmp1 = x_bi + nrmi*llen
        tmp2 = y_bi + nrmj*llen
        tmp3 = z_bi + nrmk*llen
        ip0 = int(tmp1/dx) + 1
        jp0 = int(tmp2/dy) + 1
        do kk=1, size(zgrid)-1
          if (tmp3 .ge. zgrid(kk) .and. tmp3 .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = tmp1/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = tmp2/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (tmp3 - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_k2e_w(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_k2e_w(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_k2e_w(ighost,3) = kp0
        current_state%immersed%int_xyz_k2e_w(ighost,1) = xdd
        current_state%immersed%int_xyz_k2e_w(ighost,2) = ydd
        current_state%immersed%int_xyz_k2e_w(ighost,3) = zdd
        diag4_xyz_w(ighost,1) = tmp1
        diag4_xyz_w(ighost,2) = tmp2
        diag4_xyz_w(ighost,3) = tmp3
        current_state%immersed%k2e_xyz_w(ighost,1) = tmp1
        current_state%immersed%k2e_xyz_w(ighost,2) = tmp2
        current_state%immersed%k2e_xyz_w(ighost,3) = tmp3
  
        ! s,u,v  vertical grid
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        ! s  grid (x and y are same)
        do kk=1, size(zgrid)-1
          if (tmp3 .ge. zgrid(kk) .and. tmp3 .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        zdd = (tmp3 - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_k2e_ws(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_k2e_ws(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_k2e_ws(ighost,3) = kp0
        current_state%immersed%int_xyz_k2e_ws(ighost,1) = xdd
        current_state%immersed%int_xyz_k2e_ws(ighost,2) = ydd
        current_state%immersed%int_xyz_k2e_ws(ighost,3) = zdd
        
        ! u  grid
        ip0 = int((tmp1/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(tmp2/dy) + 1
        do kk=1, size(zgrid)-1
          if (tmp3 .ge. zgrid(kk) .and. tmp3 .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = tmp1/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = tmp2/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (tmp3 - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_k2e_wu(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_k2e_wu(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_k2e_wu(ighost,3) = kp0
        current_state%immersed%int_xyz_k2e_wu(ighost,1) = xdd
        current_state%immersed%int_xyz_k2e_wu(ighost,2) = ydd
        current_state%immersed%int_xyz_k2e_wu(ighost,3) = zdd
  
        ! v  grid
        ip0 = int(tmp1/dx) + 1 
        jp0 = int((tmp2/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (tmp3 .ge. zgrid(kk) .and. tmp3 .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        ! u  grid
        xdd = tmp1/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = tmp2/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (tmp3 - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_k2e_wv(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_k2e_wv(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_k2e_wv(ighost,3) = kp0
        current_state%immersed%int_xyz_k2e_wv(ighost,1) = xdd
        current_state%immersed%int_xyz_k2e_wv(ighost,2) = ydd
        current_state%immersed%int_xyz_k2e_wv(ighost,3) = zdd
      
      end if

      ! VR ---------------------
      if (current_state%immersed%ib_type==1)then
        zgrid=current_state%global_grid%configuration%vertical%z(:)
        dz = zgrid(k+1)-zgrid(k)
        call find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen, ip_dir)
        if(z_ip .lt. z_g .or. idiag .eq. 0)then
          write(*,*)'WARNING (gridmanager): BI not found for u proc, ijk', myproc, i, j, k 
          write(*,*)'WARNING (gridmanager): idiag, pos', idiag, x_g, y_g, z_g 
          write(*,*)'WARNING (gridmanager): ighost, image point location', ighost, x_ip, y_ip, z_ip
          current_state%immersed%ghost_w(k,j,i) = -1
        end if
  
        gp_xyz_w(ighost,1) = x_g
        gp_xyz_w(ighost,2) = y_g
        gp_xyz_w(ighost,3) = z_g
        ip_xyz_w(ighost,1) = x_ip
        ip_xyz_w(ighost,2) = y_ip
        ip_xyz_w(ighost,3) = z_ip
        bi_xyz_w(ighost,1) = x_bi
        bi_xyz_w(ighost,2) = y_bi
        bi_xyz_w(ighost,3) = z_bi
        current_state%immersed%idiag_w(ighost)=idiag
  
        current_state%immersed%norms_w(ighost,1) = nrmi
        current_state%immersed%norms_w(ighost,2) = nrmj
        current_state%immersed%norms_w(ighost,3) = nrmk
   
        current_state%immersed%gp_len_w(ighost) = sqrt((x_g - x_bi)**2.0 +&
                                                       (y_g - y_bi)**2.0 +&
                                                       (z_g - z_bi)**2.0)
        current_state%immersed%ip_len_w(ighost) = sqrt((x_ip - x_bi)**2.0 +&
                                                       (y_ip - y_bi)**2.0 +&
                                                       (z_ip - z_bi)**2.0)
 
        current_state%immersed%bi_xyz_w(ighost,1) = x_bi
        current_state%immersed%bi_xyz_w(ighost,2) = y_bi
        current_state%immersed%bi_xyz_w(ighost,3) = z_bi


        ! find indices for image point interpolation
        ! u grid (w ghosts) --------------------------------------
        ip0 = int((x_ip/dx) + 0.5_DEFAULT_PRECISION)
        jp0 = int(y_ip/dy) + 1
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - (real(ip0)- 0.5_DEFAULT_PRECISION)
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_wu(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_wu(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_wu(ighost,3) = kp0
        current_state%immersed%int_xyz_wu(ighost,1) = xdd
        current_state%immersed%int_xyz_wu(ighost,2) = ydd
        current_state%immersed%int_xyz_wu(ighost,3) = zdd
                                                     
        ! v grid (w ghosts) --------------------------------------
        ip0 = int(x_ip/dx) + 1
        jp0 = int((y_ip/dy) + 0.5_DEFAULT_PRECISION)
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - (real(jp0)- 0.5_DEFAULT_PRECISION)
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_wv(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_wv(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_wv(ighost,3) = kp0
        current_state%immersed%int_xyz_wv(ighost,1) = xdd
        current_state%immersed%int_xyz_wv(ighost,2) = ydd
        current_state%immersed%int_xyz_wv(ighost,3) = zdd
                                                     
        ! s grid (w ghosts) --------------------------------------
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_ws(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_ws(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_ws(ighost,3) = kp0
        current_state%immersed%int_xyz_ws(ighost,1) = xdd
        current_state%immersed%int_xyz_ws(ighost,2) = ydd
        current_state%immersed%int_xyz_ws(ighost,3) = zdd

        ! w grid (w ghosts) --------------------------------------
        zgrid=current_state%global_grid%configuration%vertical%z(:)
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_ww(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_ww(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_ww(ighost,3) = kp0
        current_state%immersed%int_xyz_ww(ighost,1) = xdd
        current_state%immersed%int_xyz_ww(ighost,2) = ydd
        current_state%immersed%int_xyz_ww(ighost,3) = zdd


      end if ! if ib_type==1

    end do ! w ghosts

    !--------------------------------------------------------------------------------------


    do ighost = 1,nghosts_s ! scalar ghosts
      current_state%immersed%gp_ijk_s(ighost,1) = tmp_ijk_s(ighost,1)
      current_state%immersed%gp_ijk_s(ighost,2) = tmp_ijk_s(ighost,2)
      current_state%immersed%gp_ijk_s(ighost,3) = tmp_ijk_s(ighost,3)
      i = current_state%immersed%gp_ijk_s(ighost,1)
      j = current_state%immersed%gp_ijk_s(ighost,2)
      k = current_state%immersed%gp_ijk_s(ighost,3)
      if (k .gt. current_state%immersed%kmax_ji(j,i))current_state%immersed%kmax_ji(j,i) = k
      current_state%immersed%ghost_s(k,j,i) = ighost ! reset ghost index
      x_g = (i0 + i - 2 - xhalo)*dx  ! x position of s 
      y_g = (j0 + j - 2 - yhalo)*dy  ! y position of s
      z_g = current_state%global_grid%configuration%vertical%zn(k)
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      ! find list of all surface elements within 2dx of nearest
      tmp1 = 2.*max(dx,dy)
      nelem = 0
      do ielem=1,nnodes
        tmp3 = elem_cent(ielem,1) - x_g
        tmp4 = elem_cent(ielem,2) - y_g
        tmp5 = elem_cent(ielem,3) - z_g
        if (tmp3 .gt. nxg*dx/2.0) tmp3 = tmp3 - (nxg*dx) 
        if (tmp3 .lt. -nxg*dx/2.0) tmp3 = tmp3 + (nxg*dx) 
        if (tmp4 .gt. nyg*dy/2.0) tmp4 = tmp4 - (nyg*dy) 
        if (tmp4 .lt. -nyg*dy/2.0) tmp4 = tmp4 + (nyg*dy) 
        tmp2 = sqrt(tmp3*tmp3 + tmp4*tmp4 + tmp5*tmp5)
        if (tmp2 .le. tmp1) then
          nelem = nelem+1
          elem_lst(nelem) = ielem
        end if
      end do

      ! no-slip ---------------------
      if (current_state%immersed%ib_type==0)then
        call find_image_point(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen)
  
        gp_xyz_s(ighost,1) = x_g
        gp_xyz_s(ighost,2) = y_g
        gp_xyz_s(ighost,3) = z_g
        ip_xyz_s(ighost,1) = x_ip
        ip_xyz_s(ighost,2) = y_ip
        ip_xyz_s(ighost,3) = z_ip
        bi_xyz_s(ighost,1) = x_bi
        bi_xyz_s(ighost,2) = y_bi
        bi_xyz_s(ighost,3) = z_bi
        current_state%immersed%idiag_s(ighost)=idiag
        current_state%immersed%gp_len_s(ighost) = llen
  
  
        current_state%immersed%norms_s(ighost,1) = nrmi
        current_state%immersed%norms_s(ighost,2) = nrmj
        current_state%immersed%norms_s(ighost,3) = nrmk
  
        ! find indices for image point interpolation
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
    
        ! indices and interpolation coefficients for timestep callback
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%ip_ijk_ss(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_ss(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_ss(ighost,3) = kp0
        current_state%immersed%int_xyz_ss(ighost,1) = xdd
        current_state%immersed%int_xyz_ss(ighost,2) = ydd
        current_state%immersed%int_xyz_ss(ighost,3) = zdd
  
        ! same for BI run-time diagnostics
        ip0 = int(x_bi/dx) + 1
        jp0 = int(y_bi/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
  
        xdd = x_bi/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_bi/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_bi - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
    
        current_state%immersed%bi_ijk_s(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%bi_ijk_s(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%bi_ijk_s(ighost,3) = kp0
        current_state%immersed%biint_xyz_s(ighost,1) = xdd
        current_state%immersed%biint_xyz_s(ighost,2) = ydd
        current_state%immersed%biint_xyz_s(ighost,3) = zdd
        current_state%immersed%bi_xyz_s(ighost,1) = x_bi
        current_state%immersed%bi_xyz_s(ighost,2) = y_bi
        current_state%immersed%bi_xyz_s(ighost,3) = z_bi
 

      end if ! if ib_type=0


      ! VR ---------------------
      if (current_state%immersed%ib_type==1)then
        zgrid=current_state%global_grid%configuration%vertical%zn(:)
        dz = zgrid(k+1)-zgrid(k)
        call find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo,  myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi, y_bi, z_bi, nrmi, nrmj, nrmk, llen, ip_dir)
        if(z_ip .lt. z_g .or. idiag .eq. 0)then
          write(*,*)'WARNING (gridmanager): BI not found for u proc, ijk', myproc, i, j, k 
          write(*,*)'WARNING (gridmanager): idiag, pos', idiag, x_g, y_g, z_g 
          write(*,*)'WARNING (gridmanager): ighost, image point location', ighost, x_ip, y_ip, z_ip
          current_state%immersed%ghost_s(k,j,i) = -1
        end if
  
        gp_xyz_s(ighost,1) = x_g
        gp_xyz_s(ighost,2) = y_g
        gp_xyz_s(ighost,3) = z_g
        ip_xyz_s(ighost,1) = x_ip
        ip_xyz_s(ighost,2) = y_ip
        ip_xyz_s(ighost,3) = z_ip
        bi_xyz_s(ighost,1) = x_bi
        bi_xyz_s(ighost,2) = y_bi
        bi_xyz_s(ighost,3) = z_bi
        current_state%immersed%idiag_s(ighost)=idiag
  
        current_state%immersed%norms_s(ighost,1) = nrmi
        current_state%immersed%norms_s(ighost,2) = nrmj
        current_state%immersed%norms_s(ighost,3) = nrmk
   
        current_state%immersed%gp_len_s(ighost) = sqrt((x_g - x_bi)**2.0 +&
                                                       (y_g - y_bi)**2.0 +&
                                                       (z_g - z_bi)**2.0)
        current_state%immersed%ip_len_s(ighost) = sqrt((x_ip - x_bi)**2.0 +&
                                                       (y_ip - y_bi)**2.0 +&
                                                       (z_ip - z_bi)**2.0)
  




        ! s grid (s ghosts) --------------------------------------
        ip0 = int(x_ip/dx) + 1
        jp0 = int(y_ip/dy) + 1
        do kk=1, size(zgrid)-1
          if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
            kp0 = kk
          end if
        end do
        xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
        ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
        zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
        current_state%immersed%ip_ijk_ss(ighost,1) = ip0-i0+1+xhalo
        current_state%immersed%ip_ijk_ss(ighost,2) = jp0-j0+1+yhalo
        current_state%immersed%ip_ijk_ss(ighost,3) = kp0
        current_state%immersed%int_xyz_ss(ighost,1) = xdd
        current_state%immersed%int_xyz_ss(ighost,2) = ydd
        current_state%immersed%int_xyz_ss(ighost,3) = zdd


        if(surface_boundary_type == PRESCRIBED_SURFACE_FLUXES)then
          ! w grid (s ghosts) --------------------------------------
          zgrid=current_state%global_grid%configuration%vertical%z(:)
          ip0 = int(x_ip/dx) + 1
          jp0 = int(y_ip/dy) + 1
          do kk=1, size(zgrid)-1
            if (z_ip .ge. zgrid(kk) .and. z_ip .lt. zgrid(kk+1)) then
              kp0 = kk
            end if
          end do
          xdd = x_ip/dx - real(ip0) + 1.0_DEFAULT_PRECISION
          ydd = y_ip/dy - real(jp0) + 1.0_DEFAULT_PRECISION
          zdd = (z_ip - zgrid(kp0))/(zgrid(kp0+1) - zgrid(kp0))
          current_state%immersed%ip_ijk_sw(ighost,1) = ip0-i0+1+xhalo
          current_state%immersed%ip_ijk_sw(ighost,2) = jp0-j0+1+yhalo
          current_state%immersed%ip_ijk_sw(ighost,3) = kp0
          current_state%immersed%int_xyz_sw(ighost,1) = xdd
          current_state%immersed%int_xyz_sw(ighost,2) = ydd
          current_state%immersed%int_xyz_sw(ighost,3) = zdd
        end if


      end if ! if ib_type==1

    end do ! s ghosts

    !--------------------------------------------------------------------------------------

    ! set theta and theta_ref for s and w BIs
    ! s ghosts
    if(enable_theta)then
      zgrid = current_state%global_grid%configuration%vertical%zn(:)
      do ighost = 1, nghosts_s
        z_bi = bi_xyz_s(ighost,3)
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            k = kk
          end if
        end do
        tmp1 = current_state%global_grid%configuration%vertical%thref(k) 
        tmp2 = current_state%global_grid%configuration%vertical%thref(k+1)
        zdd = (z_bi-zgrid(k))/(zgrid(k+1)-zgrid(k))
        current_state%immersed%thref_surf_ib_s(ighost)=((1.-zdd)*tmp1)+(zdd*tmp2)
      end do

      ! w ghosts
      zgrid = current_state%global_grid%configuration%vertical%z(:)
      do ighost = 1, nghosts_w
        z_bi = bi_xyz_w(ighost,3)
        do kk=1, size(zgrid)-1
          if (z_bi .ge. zgrid(kk) .and. z_bi .lt. zgrid(kk+1)) then
            k = kk
          end if
        end do
        tmp1 = current_state%global_grid%configuration%vertical%thref(k) 
        tmp2 = current_state%global_grid%configuration%vertical%thref(k+1)
        zdd = (z_bi-zgrid(k))/(zgrid(k+1)-zgrid(k))
        current_state%immersed%thref_surf_ib_w(ighost)=((1.-zdd)*tmp1)+(zdd*tmp2)
      end do
    end if

    ! set neutral mixing length for ib columns

    tmp1 = ib_filter_prox*max(dx,dy)
    do i = 1, x_size
      do j = 1, y_size
        do k=1, current_state%global_grid%size(Z_INDEX)
          if (allocated(current_state%immersed%ibprox)) then
            current_state%immersed%ibprox(k,j,i) = (tanh((current_state%immersed%rneutml_sq(k,j,i)-&
                                             (tmp1/2.0_DEFAULT_PRECISION))/tmp1)+1.0_DEFAULT_PRECISION)/2.0_DEFAULT_PRECISION
          end if

          current_state%immersed%rneutml_sq(k,j,i) = 1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION/(von_karman_constant*&
           (current_state%immersed%rneutml_sq(k,j,i)+z0))**2+1.0_DEFAULT_PRECISION/current_state%rmlmax**2)
        end do
      end do
    end do
       

    ! find nearest w point BIs for IDW interpolation
    if (use_surface_bcs.and.surface_boundary_type.eq. PRESCRIBED_SURFACE_FLUXES)then
      gsearchrad = (2*current_state%local_grid%size(Y_INDEX))+10
      current_state%immersed%w2s_idw_wgt=1.e10 ! initalise to large number
      current_state%immersed%w2s_idw_idx=-1 ! initalise as flag
      k=maxval(current_state%immersed%kmax_ji)
      dz = current_state%global_grid%configuration%vertical%dzn(k-1)
      dmax=sqrt(dx*dx + dy*dy + dz*dz)*1.1
      do ighost=1,nghosts_s
        jstart = max(ighost-gsearchrad,1)
        jstop  = min(ighost+gsearchrad,nghosts_w)
        do jghost=jstart,jstop
          tmp3=sqrt((bi_xyz_w(jghost,1)-bi_xyz_s(ighost,1))**2.0+&
                    (bi_xyz_w(jghost,2)-bi_xyz_s(ighost,2))**2.0+&
                    (bi_xyz_w(jghost,3)-bi_xyz_s(ighost,3))**2.0 )
          if(tmp3.lt.maxval(current_state%immersed%w2s_idw_wgt(ighost,:),dim=1))then
            k=maxloc(current_state%immersed%w2s_idw_wgt(ighost,:),dim=1)
            current_state%immersed%w2s_idw_wgt(ighost,k)=tmp3
            current_state%immersed%w2s_idw_idx(ighost,k)=jghost
          end if
        end do

        tmp1=0.0
        do k=1,4 ! clean up and normalise
          if(current_state%immersed%w2s_idw_idx(ighost,k).eq.-1)then
            current_state%immersed%w2s_idw_wgt(ighost,k)=0.0
            current_state%immersed%w2s_idw_idx(ighost,k)=1
          else
            current_state%immersed%w2s_idw_wgt(ighost,k)=&
              max((1.0-current_state%immersed%w2s_idw_wgt(ighost,k)/dmax),0.0)
          end if
          tmp1=tmp1+current_state%immersed%w2s_idw_wgt(ighost,k)
        end do

        if(tmp1.gt.0.0)then
          current_state%immersed%w2s_idw_wgt(ighost,:)=&
          current_state%immersed%w2s_idw_wgt(ighost,:)/tmp1
        else
          write(*,*)'WARNING: no points found for IDW on ighost',ighost
        end if
      end do
    end if

    ! Diagnostics
    ! ====================================================================

 
    if (myproc == 0) then
      zgrid=current_state%global_grid%configuration%vertical%z(:)
      write(fname,'(a)')"./ib_diagnostic_files/zn_z.dat"
      open(20,file=trim(fname), form='formatted', status='replace')
      write(20,'(a)')'zn,z,dzn,dz'
      do k=1, size(zgrid)
          write(20,'(4f15.5)') current_state%global_grid%configuration%vertical%zn(k), zgrid(k),&
                               current_state%global_grid%configuration%vertical%dzn(k), &
                               current_state%global_grid%configuration%vertical%dz(k)
      enddo
      close(20)
    end if

    write(fname,'(a,i3.3,a)') "./ib_diagnostic_files/gridm_xyz_u_",myproc,".dat"
    open(20,file=trim(fname), form='formatted', status='replace')
    do ighost=1,nghosts_u
      write(20,'(12f15.5)') gp_xyz_u(ighost,1), gp_xyz_u(ighost,2), gp_xyz_u(ighost,3),&
                           ip_xyz_u(ighost,1), ip_xyz_u(ighost,2), ip_xyz_u(ighost,3),&
                           bi_xyz_u(ighost,1), bi_xyz_u(ighost,2), bi_xyz_u(ighost,3),&
                           current_state%immersed%norms_u(ighost,1),&
                           current_state%immersed%norms_u(ighost,2),&
                           current_state%immersed%norms_u(ighost,3)
    enddo
    close(20)

    write(fname,'(a,i3.3,a)') "./ib_diagnostic_files/gridm_xyz_v_",myproc,".dat"
    open(20,file=trim(fname), form='formatted', status='replace')
    do ighost=1,nghosts_v
      write(20,'(12f15.5)') gp_xyz_v(ighost,1), gp_xyz_v(ighost,2), gp_xyz_v(ighost,3),&
                           ip_xyz_v(ighost,1), ip_xyz_v(ighost,2), ip_xyz_v(ighost,3),&
                           bi_xyz_v(ighost,1), bi_xyz_v(ighost,2), bi_xyz_v(ighost,3),&
                           current_state%immersed%norms_v(ighost,1),&
                           current_state%immersed%norms_v(ighost,2),&
                           current_state%immersed%norms_v(ighost,3)
    enddo
    close(20)

    write(fname,'(a,i3.3,a)') "./ib_diagnostic_files/gridm_xyz_w_",myproc,".dat"
    open(20,file=trim(fname), form='formatted', status='replace')
    do ighost=1,nghosts_w
      write(20,'(12f15.5)') gp_xyz_w(ighost,1), gp_xyz_w(ighost,2), gp_xyz_w(ighost,3),&
                           ip_xyz_w(ighost,1), ip_xyz_w(ighost,2), ip_xyz_w(ighost,3),&
                           bi_xyz_w(ighost,1), bi_xyz_w(ighost,2), bi_xyz_w(ighost,3),&
                           current_state%immersed%norms_w(ighost,1),&
                           current_state%immersed%norms_w(ighost,2),&
                           current_state%immersed%norms_w(ighost,3)
    enddo
    close(20)

    write(fname,'(a,i3.3,a)') "./ib_diagnostic_files/gridm_xyz_s_",myproc,".dat"
    open(20,file=trim(fname), form='formatted', status='replace')
    do ighost=1,nghosts_s
      write(20,'(12f15.5)') gp_xyz_s(ighost,1), gp_xyz_s(ighost,2), gp_xyz_s(ighost,3),&
                           ip_xyz_s(ighost,1), ip_xyz_s(ighost,2), ip_xyz_s(ighost,3),&
                           bi_xyz_s(ighost,1), bi_xyz_s(ighost,2), bi_xyz_s(ighost,3),&
                           current_state%immersed%norms_s(ighost,1),&
                           current_state%immersed%norms_s(ighost,2),&
                           current_state%immersed%norms_s(ighost,3)
    enddo
    close(20)




   deallocate(elem_cent,elem_xlims,elem_ylims,elem_zlims, &
   zgrid, nodes, norms, elems, tmp_ijk_u, tmp_ijk_v, tmp_ijk_w,tmp_ijk_s, &
   diag_xyz, diag2_xyz, diag3_xyz, &
   gp_xyz_u, ip_xyz_u, bi_xyz_u, gp_xyz_v, ip_xyz_v, bi_xyz_v, &
   gp_xyz_w, ip_xyz_w, bi_xyz_w, gp_xyz_s, ip_xyz_s, bi_xyz_s )


  end subroutine initialise_immersed_boundary



  !> Will initialise the vertical grid configuration
  !! @param current_state The current model state_mod
  subroutine initialise_verticalgrid_configuration_type(current_state)
    type(model_state_type), intent(inout) :: current_state

    call allocate_vertical_grid_data(current_state%global_grid%configuration%vertical, &
       current_state%global_grid%size(Z_INDEX), current_state%number_q_fields )
    call set_up_and_smooth_grid(current_state%global_grid%configuration%vertical, &
         current_state%global_grid%configuration%vertical%kgd, current_state%global_grid%configuration%vertical%hgd, &
         size(current_state%global_grid%configuration%vertical%kgd), current_state%global_grid%size(Z_INDEX), &
         current_state%global_grid%top(Z_INDEX), options_get_integer(current_state%options_database, "nsmth"), &
         current_state%origional_vertical_grid_setup, current_state%continuation_run)
    call set_vertical_reference_profile(current_state, current_state%global_grid%configuration%vertical, &
         current_state%global_grid%size(Z_INDEX))
    
  end subroutine initialise_verticalgrid_configuration_type

  !> Sets up the vertical grid reference profile at each point
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine set_vertical_reference_profile(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k

    call calculate_initial_profiles(current_state, vertical_grid)
    call set_up_vertical_reference_properties(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call set_anelastic_pressure(current_state)
    ! 
    call set_qv_init_from_rh(current_state)

    do k=2,kkp-1
      ! for diffusion onto p-level from below
      vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
      ! for diffusion onto p-level from above
      vertical_grid%cza(k)=(vertical_grid%rho(k)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k+1))
      vertical_grid%czg(k)=-vertical_grid%czb(k)-vertical_grid%cza(k)
      if (k .gt. 2) vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    end do
    do k=2,kkp-1
      ! advection onto p-level from below
      vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k) 
      ! advection onto p-level from above
      vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k) 
    end do
    do k=2,kkp-1
      ! advection onto w-level (K) from below
      vertical_grid%tzd1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k)/vertical_grid%rho(k)
      ! advection onto w-level (K) from above
      vertical_grid%tzd2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdzn(k+1)*vertical_grid%rhon(k+1)/vertical_grid%rho(k)
    end do
    k=kkp
    vertical_grid%czb(k)=(vertical_grid%rho(k-1)/vertical_grid%rhon(k))/(vertical_grid%dz(k)*vertical_grid%dzn(k))
    vertical_grid%cza(k)=0.0_DEFAULT_PRECISION
    vertical_grid%czg(k)=-vertical_grid%czb(k)
    vertical_grid%czh(k)=vertical_grid%czb(k)*vertical_grid%cza(k-1)
    vertical_grid%tzc2(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k)/vertical_grid%rhon(k)
    vertical_grid%tzc1(k)=0.25_DEFAULT_PRECISION*vertical_grid%rdz(k)*vertical_grid%rho(k-1)/vertical_grid%rhon(k) 
    vertical_grid%czn=vertical_grid%dzn(2)*0.5_DEFAULT_PRECISION
    vertical_grid%zlogm=log(1.0_DEFAULT_PRECISION+vertical_grid%zn(2)/z0)
    vertical_grid%zlogth=log((vertical_grid%zn(2)+z0)/z0th)
    vertical_grid%vk_on_zlogm=von_karman_constant/vertical_grid%zlogm
    call setup_reference_state_liquid_water_temperature_and_saturation(&
         current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call calculate_mixing_length_for_neutral_case(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
    call set_buoyancy_coefficient(current_state, vertical_grid, current_state%global_grid%size(Z_INDEX))
  end subroutine set_vertical_reference_profile

  !> Calculates the initial profiles for U, V, TH & Q if required
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  subroutine calculate_initial_profiles(current_state, vertical_grid)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid

    integer :: nq_init ! The number of q fields to initialize
    integer :: nzq     ! The number of input levels for q_init
    integer :: i,j,n, k,  i_tracer ! loop counters
    integer :: iq  ! temporary q varible index

    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: f_init_pl_q       ! Initial node values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_q      ! Initial node height values for q variables
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_theta  ! Initial node values for potential temperature variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_theta  ! Initial node height values for potential temperature variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_u      ! Initial node values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_u      ! Initial node height values for u variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_v      ! Initial node values for v variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_v      ! Initial node height values for v variable

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_thref   ! Initial node values for thref
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_thref   ! Initial node height values for thref

    logical :: l_init_pl_u     ! if .true. then initialize u field
    logical :: l_init_pl_v     ! if .true. then initialize v field
    logical :: l_init_pl_theta ! if .true. then initialize potential temperature field
    logical :: l_init_pl_rh    ! if .true. then initialize relative humidity field
    logical :: l_init_pl_q     ! if .true. then initialize q fields
    logical :: l_thref         ! if .true. then initialize thref profile (overrides thref0)
    logical :: l_matchthref    ! if .true. then initialize thref to be the same as theta_init

    character(len=STRING_LENGTH), dimension(:), allocatable :: names_init_pl_q ! names of q variables to initialize
    
    real(kind=DEFAULT_PRECISION), allocatable :: f_init_pl_q_tmp(:) !temporary 1D storage of initial q field
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation

    real(kind=DEFAULT_PRECISION) :: zztop ! top of the domain
    real(kind=DEFAULT_PRECISION) :: qsat

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))
    
    zztop = current_state%global_grid%top(Z_INDEX)

    ! Initialize everything to zero.  This won't make sense for theta.
    vertical_grid%q_init = 0.0_DEFAULT_PRECISION
    vertical_grid%u_init = 0.0_DEFAULT_PRECISION
    vertical_grid%v_init = 0.0_DEFAULT_PRECISION
    vertical_grid%theta_init = 0.0_DEFAULT_PRECISION

    l_init_pl_theta=options_get_logical(current_state%options_database, "l_init_pl_theta")
    l_init_pl_rh=options_get_logical(current_state%options_database, "l_init_pl_rh") 
    l_init_pl_q=options_get_logical(current_state%options_database, "l_init_pl_q")
    if (l_init_pl_q) then
      allocate(names_init_pl_q(options_get_array_size(current_state%options_database, "names_init_pl_q")))
      call options_get_string_array(current_state%options_database, "names_init_pl_q", names_init_pl_q)
      if (size(names_init_pl_q) .eq. 0) then
        call log_master_log(LOG_ERROR, "Model configured with l_init_pl_q=.true., but "//&
                                       "no names_init_pl_q have been specified")
      end if
      do n = 1,size(names_init_pl_q)
         if (trim(names_init_pl_q(n)) .eq. 'vapour' .and. l_init_pl_rh) then 
            call log_master_log(LOG_ERROR, "Initialisation of vapour and RH - STOP")
         endif
      enddo
    end if
    l_init_pl_u=options_get_logical(current_state%options_database, "l_init_pl_u")
    l_init_pl_v=options_get_logical(current_state%options_database, "l_init_pl_v")

    l_thref=options_get_logical(current_state%options_database, "l_thref")
    l_matchthref=options_get_logical(current_state%options_database, "l_matchthref")

    if (.not. current_state%continuation_run) then  ! For continuations, ensure thref is that from checkpoint
      if (l_thref)then
        if (.not. l_matchthref)then
          allocate(z_thref(options_get_array_size(current_state%options_database, "z_thref")), &
               f_thref(options_get_array_size(current_state%options_database, "f_thref")))
          call options_get_real_array(current_state%options_database, "z_thref", z_thref)
          call options_get_real_array(current_state%options_database, "f_thref", f_thref)
          call check_top(zztop, z_thref(size(z_thref)), 'z_thref')
          call check_input_levels(size(z_thref), size(f_thref), "f_thref")
          zgrid=current_state%global_grid%configuration%vertical%zn(:)
          call piecewise_linear_1d(z_thref(1:size(z_thref)), f_thref(1:size(f_thref)), zgrid, &
             current_state%global_grid%configuration%vertical%thref)
          deallocate(z_thref, f_thref)
        end if
      else
        current_state%global_grid%configuration%vertical%thref(:)=current_state%thref0
      end if
    end if 

    if (l_init_pl_theta)then
      allocate(z_init_pl_theta(options_get_array_size(current_state%options_database, "z_init_pl_theta")), &
             f_init_pl_theta(options_get_array_size(current_state%options_database, "f_init_pl_theta")))
      call options_get_real_array(current_state%options_database, "z_init_pl_theta", z_init_pl_theta)
      call options_get_real_array(current_state%options_database, "f_init_pl_theta", f_init_pl_theta)
      call check_top(zztop, z_init_pl_theta(size(z_init_pl_theta)), 'z_init_pl_theta')
      call check_input_levels(size(z_init_pl_theta), size(f_init_pl_theta), "f_init_pl_theta")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_theta(1:size(z_init_pl_theta)), f_init_pl_theta(1:size(f_init_pl_theta)), zgrid, &
         current_state%global_grid%configuration%vertical%theta_init)
      if (l_matchthref .and. .not. current_state%continuation_run) then ! For continuations, ensure thref is that from checkpoint
         if(.not. current_state%use_anelastic_equations) then
           call log_master_log(LOG_ERROR, "Non-anelastic equation set and l_matchthref are incompatible")
         end if
         current_state%global_grid%configuration%vertical%thref = current_state%global_grid%configuration%vertical%theta_init
      end if
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%th%data(:,j,i) = current_state%global_grid%configuration%vertical%theta_init(:) - &
                 current_state%global_grid%configuration%vertical%thref(:) 
          end do
        end do
      end if
      deallocate(z_init_pl_theta, f_init_pl_theta)
    end if

    if (l_init_pl_u)then
      allocate(z_init_pl_u(options_get_array_size(current_state%options_database, "z_init_pl_u")), &
             f_init_pl_u(options_get_array_size(current_state%options_database, "f_init_pl_u")))
      call options_get_real_array(current_state%options_database, "z_init_pl_u", z_init_pl_u)
      call options_get_real_array(current_state%options_database, "f_init_pl_u", f_init_pl_u)
      call check_top(zztop, z_init_pl_u(size(z_init_pl_u)), 'z_init_pl_u')
      call check_input_levels(size(z_init_pl_u), size(f_init_pl_u), "f_init_pl_u")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_u(1:size(z_init_pl_u)), f_init_pl_u(1:size(f_init_pl_u)), &
         zgrid, current_state%global_grid%configuration%vertical%u_init)
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%u%data(:,j,i) = current_state%global_grid%configuration%vertical%u_init(:)
          end do
        end do
      end if
      deallocate(z_init_pl_u, f_init_pl_u)
    end if

    if (l_init_pl_v)then
      allocate(z_init_pl_v(options_get_array_size(current_state%options_database, "z_init_pl_v")), &
             f_init_pl_v(options_get_array_size(current_state%options_database, "f_init_pl_v")))
      call options_get_real_array(current_state%options_database, "z_init_pl_v", z_init_pl_v)
      call options_get_real_array(current_state%options_database, "f_init_pl_v", f_init_pl_v)
      call check_top(zztop, z_init_pl_v(size(z_init_pl_v)), 'z_init_pl_v')
      call check_input_levels(size(z_init_pl_v), size(f_init_pl_v), "f_init_pl_v")
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      call piecewise_linear_1d(z_init_pl_v(1:size(z_init_pl_v)), f_init_pl_v(1:size(f_init_pl_v)), &
         zgrid, current_state%global_grid%configuration%vertical%v_init)
      if (.not. current_state%continuation_run) then
        do i=current_state%local_grid%local_domain_start_index(X_INDEX), current_state%local_grid%local_domain_end_index(X_INDEX)
          do j=current_state%local_grid%local_domain_start_index(Y_INDEX), current_state%local_grid%local_domain_end_index(Y_INDEX)
            current_state%v%data(:,j,i) = current_state%global_grid%configuration%vertical%v_init(:)
          end do
        end do
      end if
      deallocate(z_init_pl_v, f_init_pl_v)
    end if

    if (l_init_pl_q)then
      nq_init=size(names_init_pl_q)
      allocate(z_init_pl_q(options_get_array_size(current_state%options_database, "z_init_pl_q")))
      call options_get_real_array(current_state%options_database, "z_init_pl_q", z_init_pl_q)
      nzq=size(z_init_pl_q)
      call check_top(zztop, z_init_pl_q(nzq), 'z_init_pl_q')
      zgrid=current_state%global_grid%configuration%vertical%zn(:)
      allocate(f_init_pl_q_tmp(options_get_array_size(current_state%options_database, "f_init_pl_q")))
      if (nq_init*nzq .ne. size(f_init_pl_q_tmp)) then
        call log_master_log(LOG_ERROR, "There is a mismatch between the number of initial moisture heights, "// &
                                       "size(z_init_pl_q)="//trim(conv_to_string(nzq))//                        &
                                       ", and the initial moisture values, "//                                  &
                                       "size(f_init_pl_q)="//trim(conv_to_string(size(f_init_pl_q_tmp)))//      &
                                       ".  The length of f_init_pl_q should equal the length of z_init_pl_q "// &
                                       "multiplied by the number of names_init_pl_q.")
      end if
      call options_get_real_array(current_state%options_database, "f_init_pl_q", f_init_pl_q_tmp)
      allocate(f_init_pl_q(nzq, nq_init))
      f_init_pl_q(1:nzq, 1:nq_init)=reshape(f_init_pl_q_tmp, (/nzq, nq_init/))
      do n=1, nq_init
         iq=get_q_index(trim(names_init_pl_q(n)), 'piecewise_initialization')
         call check_input_levels(size(z_init_pl_q), size(f_init_pl_q(1:nzq,n)), "f_init_pl_q")
        call piecewise_linear_1d(z_init_pl_q(1:nzq), f_init_pl_q(1:nzq,n), zgrid, &
           current_state%global_grid%configuration%vertical%q_init(:,iq))
        if (.not. current_state%continuation_run) then
          do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
               current_state%local_grid%local_domain_end_index(X_INDEX)
            do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                 current_state%local_grid%local_domain_end_index(Y_INDEX)
              current_state%q(iq)%data(:,j,i) = current_state%global_grid%configuration%vertical%q_init(:, iq)
            end do
          end do
        end if
      end do
      deallocate(f_init_pl_q_tmp, z_init_pl_q, f_init_pl_q, names_init_pl_q)
    end if

    ! Override with RCEMIP initial conditions, if logical is set.
    if (options_get_logical(current_state%options_database, "l_rcemip_initial")) &
        call rcemip_init(current_state)
   
    if (current_state%n_tracers .gt. 0) then
      if (.not. current_state%continuation_run) then
        do i_tracer = 1,current_state%n_tracers
          current_state%tracer(i_tracer)%data(:,:,:) = 0.0_DEFAULT_PRECISION
          current_state%ztracer(i_tracer)%data(:,:,:) = 0.0_DEFAULT_PRECISION
        end do
      end if
    end if
    
    if (current_state%traj_tracer_index .gt. 0) then
      if (.not. current_state%continuation_run) then
        call reinitialise_trajectories(current_state)
      end if
    end if
    
    deallocate(zgrid)      
  end subroutine calculate_initial_profiles

  !> Calculates the mixing length for the neutral case
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine calculate_mixing_length_for_neutral_case(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k

    do k=2, kkp-1
      vertical_grid%rneutml(k)=sqrt(1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION/(von_karman_constant*&
           (vertical_grid%z(k)+z0))**2+1.0_DEFAULT_PRECISION/current_state%rmlmax**2) )
      vertical_grid%rneutml_sq(k)=vertical_grid%rneutml(k)*vertical_grid%rneutml(k)
    end do    
  end subroutine calculate_mixing_length_for_neutral_case 

  !> Sets the buoyancy coefficient from the grid configuration and configuration
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine set_buoyancy_coefficient(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout), target :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    integer :: k    

    if (.not. current_state%passive_th) then
      if(current_state%use_anelastic_equations)then                                                      
        do k=1, kkp-1          
          vertical_grid%buoy_co(k)=cp*(vertical_grid%prefn(k)**r_over_cp-vertical_grid%prefn(k+1)**r_over_cp)/&
               ((current_state%surface_reference_pressure**r_over_cp)*vertical_grid%dzn(k+1))
        end do
      else                                                                     
        vertical_grid%buoy_co(1:kkp-1)=G/current_state%thref0        ! _Boussinesq
      end if
      ! Dummy value at top level
      vertical_grid%buoy_co(kkp)=0.
    else
      vertical_grid%buoy_co(:)=0.
    end if
  end subroutine set_buoyancy_coefficient

  !> Setting up reference state liquid water temperature and saturation mixing ratio on main levels.
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp Number of grid points in a vertical column
  subroutine setup_reference_state_liquid_water_temperature_and_saturation(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp

    real(kind=DEFAULT_PRECISION) :: delta_t=1.0_DEFAULT_PRECISION, qlinit, tinit, qsatin, dqsatdtin, dsatfacin
    integer :: iter, k
    
    do k=1,kkp
       vertical_grid%tref(k)=vertical_grid%thref(k)*(vertical_grid%prefn(k)/current_state%surface_reference_pressure)**r_over_cp
      vertical_grid%tstarpr(k)=0.0_DEFAULT_PRECISION
   end do
   if (current_state%th%active) then
      ! PREFRCP is used and hence calculated if theta is active
      do k = 1,kkp   
         vertical_grid%prefrcp(k)=(current_state%surface_reference_pressure/vertical_grid%prefn(k))**r_over_cp
         vertical_grid%rprefrcp(k)=1.0_DEFAULT_PRECISION/vertical_grid%prefrcp(k)
         ! Denotion between setup run and chain run in LEM - need to consider here too
         vertical_grid%qsat(k)=qsaturation(vertical_grid%tref(k), 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k))
         vertical_grid%dqsatdt(k)=(qsaturation(vertical_grid%tref(k)+delta_t, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k)) -&
              qsaturation(vertical_grid%tref(k)-delta_t, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k)))/&
              (2.0_DEFAULT_PRECISION*delta_t)
         vertical_grid%qsatfac(k)=1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION+rlvap_over_cp*vertical_grid%dqsatdt(k))
      end do
      if (current_state%calculate_th_and_q_init) then
         do k=1,kkp
            !       !Note that at this point THETA_INIT and QINIT(IQ=1) are still
            !       !theta_l and q_t, as read in from configuration.
            !       ! start from input QL profile
            qlinit=qinit(k, current_state%liquid_water_mixing_ratio_index)
            do iter=1,5
               !         ! calculate T and thence new q_l from Taylor expansion
               !         ! keeping theta_l and q_t fixed
               !         ! Note theta_l = theta - (L/c_p)*q_l here
               tinit     = vertical_grid%theta_init(k)*vertical_grid%rprefrcp(k) + rlvap_over_cp*qlinit
               qsatin    = qsaturation(tinit, 0.01_DEFAULT_PRECISION*vertical_grid%prefn(k))
               dqsatdtin = dqwsatdt(qsatin, tinit)
               dsatfacin=( 1.0_DEFAULT_PRECISION/(1.0_DEFAULT_PRECISION + rlvap_over_cp*dqsatdtin*vertical_grid%rprefrcp(k)))
               qlinit=max(0.0_DEFAULT_PRECISION, (qinit(k, current_state%water_vapour_mixing_ratio_index)-&
                    (qsatin+dqsatdtin*(vertical_grid%theta_init(k)*vertical_grid%rprefrcp(k)-tinit) ))*dsatfacin)
            end do
            qinit(k, current_state%liquid_water_mixing_ratio_index)=qlinit
            qinit(k, current_state%water_vapour_mixing_ratio_index)=qinit(k,current_state%water_vapour_mixing_ratio_index)-qlinit
            vertical_grid%theta_init(k)=vertical_grid%theta_init(k)+rlvap_over_cp*qlinit

            ! Denotion between setup run and chain run in LEM - need to consider here too
            vertical_grid%tstarpr(k)= tinit-vertical_grid%tref(k)
            vertical_grid%qsat(k)=qsatin
            vertical_grid%dqsatdt(k)=dqsatdtin        
            vertical_grid%qsatfac(k)= ( 1.0_DEFAULT_PRECISION/ ( 1.0_DEFAULT_PRECISION + rlvap_over_cp*dqsatdtin ) )
         end do
      endif
   end if
  end subroutine setup_reference_state_liquid_water_temperature_and_saturation  

  !> Sets up the reference properties for the vertical grid at each point
  !! @param current_state The current model state_mod
  !! @param vertical_grid The vertical grid that we are working on
  !! @param kkp The number of grid points in a vertical column
  subroutine set_up_vertical_reference_properties(current_state, vertical_grid, kkp)
    type(model_state_type), intent(inout) :: current_state
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: kkp
    
    integer :: k

    do k=1,kkp
!       vertical_grid%thref(k)=current_state%thref0
!       vertical_grid%theta_init(k)=vertical_grid%thref(k) ! In LEM this can also be set from configuration (TODO)       
       vertical_grid%prefn(k)=0.0_DEFAULT_PRECISION
       vertical_grid%pdiff(k)=0.0_DEFAULT_PRECISION
       vertical_grid%rho(k)=current_state%rhobous
       vertical_grid%rhon(k)=current_state%rhobous
    end do
    do k=1,kkp-1
      vertical_grid%dthref(k)=vertical_grid%thref(k+1)-vertical_grid%thref(k)
    end do
    vertical_grid%dthref(kkp)=0.0_DEFAULT_PRECISION
  end subroutine set_up_vertical_reference_properties  

  !> Sets up and smooths the vertical grid. This is based upon the grid configuration already read in
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  !! @param origional_setup To use the origional vertical grid setup routine or the new one
  subroutine set_up_and_smooth_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth, origional_setup, continuation_run)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp, nsmth
    real(kind=DEFAULT_PRECISION),intent(in) :: zztop
    logical, intent(in) :: origional_setup, continuation_run

    integer :: k

    if (.not. continuation_run) then
      if (origional_setup) then
        call original_vertical_grid_setup(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth)
      else
        call new_vertical_grid_setup(vertical_grid, kgd, kkp, zztop)
      end if
    end if
    
    ! Regardless of the vertical grid computation method, set the level deltas
    do k=2,kkp
       vertical_grid%dz(k)=vertical_grid%z(k)-vertical_grid%z(k-1)
       vertical_grid%dzn(k)= vertical_grid%zn(k)-vertical_grid%zn(k-1)                                                     
       vertical_grid%rdz(k)=1./vertical_grid%dz(k)                                                          
       vertical_grid%rdzn(k)=1./vertical_grid%dzn(k)                                                        
    end do
    vertical_grid%dzn(1)=0.d0
  end subroutine set_up_and_smooth_grid

  !> The original vertical grid setup and smoothing as per the LEM
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  subroutine original_vertical_grid_setup(vertical_grid, kgd, hgd, ninitp, kkp, zztop, nsmth)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp, nsmth
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    integer :: n, k
    
    call create_linear_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop)
    ! Smooth grid
    vertical_grid%z(1)=0.0_DEFAULT_PRECISION
    vertical_grid%z(kkp)=zztop
    do n=1,nsmth
      do k=2,kkp
        vertical_grid%zn(k)=0.5_DEFAULT_PRECISION*(vertical_grid%z(k)+vertical_grid%z(k-1))
      end do
      do k=2,kkp-1
        vertical_grid%z(k)=0.5_DEFAULT_PRECISION*(vertical_grid%zn(k)+vertical_grid%zn(k+1))
      end do
    end do
    ! Fourth order interpolation
    do k=3,kkp-1
      vertical_grid%zn(k)=0.0625_DEFAULT_PRECISION*(9.0_DEFAULT_PRECISION*&
           (vertical_grid%z(k-1)+vertical_grid%z(k))-vertical_grid%z(k+1)-vertical_grid%z(k-2))
    end do
    vertical_grid%zn(2)=0.5_DEFAULT_PRECISION*(vertical_grid%z(1)+vertical_grid%z(2))
    vertical_grid%zn(1)=-vertical_grid%zn(2)
    vertical_grid%zn(kkp)=0.5_DEFAULT_PRECISION*(vertical_grid%z(kkp-1)+vertical_grid%z(kkp))
  end subroutine original_vertical_grid_setup

  !> The newer vertical grid setup routine
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  subroutine new_vertical_grid_setup(vertical_grid, kgd, kkp, zztop)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    integer, intent(in) :: kkp
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    real(kind=DEFAULT_PRECISION) :: a(2*kkp), r1, d1, dd, d0, a0
    logical :: first_gt=.true.
    integer :: k, k0

    r1=1.10_DEFAULT_PRECISION
    d1=10.0_DEFAULT_PRECISION

    dd=0.5_DEFAULT_PRECISION*d1
    d0=dd
    a(1)=-dd
    a(2)=0.0_DEFAULT_PRECISION
    do k=3, kkp*2
      if (d0 .gt. dd .or. k==3) then
        a(k)=a(k-1)+dd
        if (.not. (dd .gt. 25.0_DEFAULT_PRECISION .and. a(k) .lt. 2000.0_DEFAULT_PRECISION)) then
          if (a(k) .lt. 2000.0_DEFAULT_PRECISION) then
            dd=dd*r1
          else
            dd=dd*(1.0_DEFAULT_PRECISION+(r1-1.0_DEFAULT_PRECISION)/1.5_DEFAULT_PRECISION)
          end if
        end if
        d0=(zztop-a(k))/real(kkp*2-k, kind=DEFAULT_PRECISION)
      else
        if (first_gt) then
          k0=k
          a0=a(k-1)+d0
          first_gt=.false.
        end if
        a(k)=a0+real(k-k0, kind=DEFAULT_PRECISION)*d0
      end if
    end do    
    
    do k=1, kkp
      vertical_grid%z(k)=a(k*2)
      vertical_grid%zn(k)=a(2*k-1)
    end do
  end subroutine new_vertical_grid_setup

  !> Creates the linear vertical grid based upon the configuration properties.
  !!
  !! This will correspond the vertical k grid points to the configuration which might be only
  !! a few reference properties.
  !! @param vertical_grid The vertical grid
  !! @param kgd The grid heights per division
  !! @param hgd The real world (m) heights per division
  !! @param kkp Number of points in the vertical domain
  !! @param zztop The real world (m) height of the top of the column
  !! @param nsmth Number of smoothing iterations to run on the grid
  subroutine create_linear_grid(vertical_grid, kgd, hgd, ninitp, kkp, zztop)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, dimension(:), intent(in) :: kgd
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: hgd
    integer, intent(in) :: ninitp, kkp
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop

    integer :: kmax, k, i
    real(kind=DEFAULT_PRECISION) :: zmax

    vertical_grid%z(1) = 0.0_DEFAULT_PRECISION
    kmax=kgd(1)
    zmax=0.0_DEFAULT_PRECISION
    if (kgd(1) .gt. 1) then
      do k=1,kgd(1)
        ! Loop up to first division point
        vertical_grid%z(k)=real(k-1, kind=DEFAULT_PRECISION)*hgd(1)/real(kgd(1)-1, kind=DEFAULT_PRECISION)
      end do
      do i=2,ninitp
        if(kgd(i) .gt. 0) then
          kmax=kgd(i)
          zmax=hgd(i)
          
          do k=kgd(i-1)+1,kgd(i)
            vertical_grid%z(k)=hgd(i-1)+(hgd(i)-hgd(i-1))*real(k-kgd(i-1), kind=DEFAULT_PRECISION)&
                 /real(kgd(i)-kgd(i-1), kind=DEFAULT_PRECISION)
          end do
        end if
      end do
   end if
    if(kmax .lt. kkp)then
      do k=kmax,kkp
        ! Handle any points above the kth max division
        vertical_grid%z(k)=zmax+(zztop-zmax)*real(k-kmax, kind=DEFAULT_PRECISION)/real(kkp-kmax, kind=DEFAULT_PRECISION)
     end do
    end if
  end subroutine create_linear_grid

  !> Allocates the data required for the vertical grid configuration
  !! @param vertical_grid The vertical grid that we are working with
  !! @param n The number of grid points in the vertical
  !! @param nq The number of q fields
  subroutine allocate_vertical_grid_data(vertical_grid, n, nq)
    type(vertical_grid_configuration_type), intent(inout) :: vertical_grid
    integer, intent(in) :: n
    integer, intent(in) :: nq

    allocate(vertical_grid%dz(n), vertical_grid%dzn(n),&
         vertical_grid%czb(n), vertical_grid%cza(n), vertical_grid%czg(n), vertical_grid%czh(n),&
         vertical_grid%rdz(n), vertical_grid%rdzn(n), vertical_grid%tzc1(n), vertical_grid%tzc2(n),&
         vertical_grid%tzd1(n), vertical_grid%tzd2(n), vertical_grid%theta_init(n), vertical_grid%temp_init(n), &
         vertical_grid%rh_init(n), &
         vertical_grid%tref(n), vertical_grid%prefn(n), vertical_grid%pdiff(n), vertical_grid%prefrcp(n), &
         vertical_grid%rprefrcp(n), vertical_grid%rho(n), vertical_grid%rhon(n), vertical_grid%tstarpr(n), &
         vertical_grid%qsat(n), vertical_grid%dqsatdt(n), vertical_grid%qsatfac(n), vertical_grid%dthref(n), &
         vertical_grid%rneutml(n), vertical_grid%rneutml_sq(n), vertical_grid%buoy_co(n), &
         vertical_grid%u_init(n), vertical_grid%v_init(n), vertical_grid%theta_rand(n), vertical_grid%w_rand(n), &
         vertical_grid%w_subs(n), vertical_grid%u_force(n), vertical_grid%v_force(n), vertical_grid%theta_force(n))

    if (.not. allocated(vertical_grid%thref)) allocate(vertical_grid%thref(n))
    if (.not. allocated(vertical_grid%z)) allocate(vertical_grid%z(n))
    if (.not. allocated(vertical_grid%zn)) allocate(vertical_grid%zn(n))

    allocate(vertical_grid%q_rand(n,nq), vertical_grid%q_init(n,nq), vertical_grid%q_force(n,nq))
  end subroutine allocate_vertical_grid_data  

  !> Initialises the horizontal grid configurations
  !! @param current_state The current model state_mod
  subroutine initialise_horizontalgrid_configuration_types(current_state)
    type(model_state_type), intent(inout) :: current_state

    current_state%global_grid%configuration%horizontal%dx = merge(real(current_state%global_grid%resolution(X_INDEX)), &
         DEFAULT_SPACING, current_state%global_grid%active(X_INDEX))
    current_state%global_grid%configuration%horizontal%dy = merge(real(current_state%global_grid%resolution(Y_INDEX)), &
         DEFAULT_SPACING, current_state%global_grid%active(Y_INDEX))

    current_state%global_grid%configuration%horizontal%cx=1./current_state%global_grid%configuration%horizontal%dx
    current_state%global_grid%configuration%horizontal%cy=1./current_state%global_grid%configuration%horizontal%dy
    current_state%global_grid%configuration%horizontal%cx2=current_state%global_grid%configuration%horizontal%cx ** 2
    current_state%global_grid%configuration%horizontal%cy2=current_state%global_grid%configuration%horizontal%cy ** 2
    current_state%global_grid%configuration%horizontal%cxy=current_state%global_grid%configuration%horizontal%cx * &
         current_state%global_grid%configuration%horizontal%cy
    current_state%global_grid%configuration%horizontal%tcx=&
         0.25_DEFAULT_PRECISION/current_state%global_grid%configuration%horizontal%dx
    current_state%global_grid%configuration%horizontal%tcy=&
         0.25_DEFAULT_PRECISION/current_state%global_grid%configuration%horizontal%dy
  end subroutine initialise_horizontalgrid_configuration_types  

  !> Set reference profile of potential temperature for the Boussinesq/Anelastic approximation
  !! Note that this is not in general the same as the profile defining the initial vertical distribution of potential
  !! temperature for the integration. In particular, while the later may contain sharp changes in gradient representing
  !! a capping inversion or the tropopause, for example, the reference profile should be smooth
  !! @param current_state The current model state
  subroutine set_anelastic_pressure(current_state)
    type(model_state_type), intent(inout) :: current_state
    
    if (current_state%use_anelastic_equations) then
      call compute_anelastic_pressure_profile_and_density(current_state)
    else
      if (current_state%passive_th) then
        current_state%global_grid%configuration%vertical%prefn=0.0_DEFAULT_PRECISION
      else
        call compute_anelastic_pressure_profile_only(current_state)
      end if
        current_state%global_grid%configuration%vertical%rho=current_state%rhobous
        current_state%global_grid%configuration%vertical%rhon=current_state%rhobous
        current_state%global_grid%configuration%vertical%pdiff=0.0_DEFAULT_PRECISION
    end if
  end subroutine set_anelastic_pressure  

  !> Computes the anelastic pressure only - if we are using Boussinesq approximation
  !! @param current_state The current model state
  subroutine compute_anelastic_pressure_profile_only(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: ipass, k
    real(kind=DEFAULT_PRECISION) ::    p0    &!pressure at z=0 adjustments made after 1st iteration so P0=PSF after 2nd iteration
        ,   ptop  &!pressure at z=ZN(KKP)
        , thprof(current_state%local_grid%size(Z_INDEX))

    
    ! TODO: NOTE - we are mocking in thprof at the moment, this should be read from a configuration and used here instead
    thprof=0.0_DEFAULT_PRECISION
    ptop=0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%pdiff(current_state%local_grid%size(Z_INDEX))=0.0_DEFAULT_PRECISION

    do ipass=1,2 ! _after first pass adjust PTOP
      current_state%global_grid%configuration%vertical%prefn(current_state%local_grid%size(Z_INDEX))=&
           (ptop/current_state%surface_reference_pressure)**r_over_cp
      do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
        current_state%global_grid%configuration%vertical%pdiff(k)=G*&
                     current_state%global_grid%configuration%vertical%dzn(k+1)/(0.5_DEFAULT_PRECISION*cp*&
                     (current_state%global_grid%configuration%vertical%thref(k)+&
                     current_state%global_grid%configuration%vertical%thref(k+1)))
      end do
      do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
        current_state%global_grid%configuration%vertical%prefn(k)=&
             current_state%global_grid%configuration%vertical%prefn(k+1)+&
             current_state%global_grid%configuration%vertical%pdiff(k)
      end do
      do k=current_state%local_grid%size(Z_INDEX),1,-1
        current_state%global_grid%configuration%vertical%prefn(k)=current_state%surface_reference_pressure*&
             current_state%global_grid%configuration%vertical%prefn(k)**(1.0_DEFAULT_PRECISION/r_over_cp)
      end do
      p0=0.5_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%prefn(1)+&
             current_state%global_grid%configuration%vertical%prefn(2))
      if (ipass .eq. 1) then                                                       
        ptop=current_state%surface_pressure**r_over_cp+ptop**r_over_cp-p0**r_over_cp
        if (ptop .le. 0.0_DEFAULT_PRECISION .and. current_state%parallel%my_rank==0) then
          call log_log(LOG_ERROR, "Negative ptop in setup of anelastic. Need a warmer THREF or different setup options")
        end if
        ptop=ptop**(1.0_DEFAULT_PRECISION/r_over_cp)
      end if
    end do
  end subroutine compute_anelastic_pressure_profile_only  

  !> Computes the anelastic pressure and density - if we are using Anelastic approximation
  !! @param current_state The current model state
  subroutine compute_anelastic_pressure_profile_and_density(current_state)
    type(model_state_type), intent(inout) :: current_state

    integer :: ipass, k
    real(kind=DEFAULT_PRECISION) ::    p0    &!pressure at z=0 adjustments made after 1st iteration so P0=PSF after 2nd iteration
        ,   ptop  &!pressure at z=ZN(KKP)
        ,   thfactor !factor for multiplying TH profile (if IADJANELP=2)

    ptop=0.0_DEFAULT_PRECISION
    current_state%global_grid%configuration%vertical%pdiff(current_state%local_grid%size(Z_INDEX))=0.0_DEFAULT_PRECISION
    do ipass=1,2 ! _after first pass, may adjust basic states
        if (ipass .eq. 1 .or. ANELASTIC_PROFILE_MODE .gt. 1) then                                     
            current_state%global_grid%configuration%vertical%prefn(current_state%local_grid%size(Z_INDEX))=&
                 (ptop/current_state%surface_reference_pressure)**r_over_cp
            do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
                current_state%global_grid%configuration%vertical%pdiff(k)=G*&
                     current_state%global_grid%configuration%vertical%dzn(k+1)/(0.5_DEFAULT_PRECISION*cp*&
                     (current_state%global_grid%configuration%vertical%thref(k)+&
                     current_state%global_grid%configuration%vertical%thref(k+1)))                
            end do
            do k=current_state%local_grid%size(Z_INDEX)-1,1,-1
                current_state%global_grid%configuration%vertical%prefn(k)=&
                     current_state%global_grid%configuration%vertical%prefn(k+1)+&
                     current_state%global_grid%configuration%vertical%pdiff(k)
            end do
            do k=current_state%local_grid%size(Z_INDEX),1,-1
                current_state%global_grid%configuration%vertical%prefn(k)=current_state%surface_reference_pressure*&
                     current_state%global_grid%configuration%vertical%prefn(k)**(1.0_DEFAULT_PRECISION/r_over_cp)
            end do           
        end if                                                                    
        p0=0.5_DEFAULT_PRECISION*(current_state%global_grid%configuration%vertical%prefn(1)+&
             current_state%global_grid%configuration%vertical%prefn(2))
        !       !-------------------------------------------------------------
        !       ! _If IADJANELP>1 we adjust the basic states to ensure P0=PSF,
        !       !                                as follows:
        !       !    IADJANELP=2     adjust THREF profile by constant factor
        !       !    IADJANELP=3     adjust PSF
        !       !    IADJANELP=4     adjust PTOP
        !       ! _Option 3 tends to give rather large changes in PSF, so
        !       !   I prefer 2 or 4 for most purposes
        !       !-------------------------------------------------------------
        if (ipass .eq. 1 .and. ANELASTIC_PROFILE_MODE .eq. 2) then                                    
            !             ! _adjust THREF profile by constant factor to enforce
            !             !    P0 = (fixed) PSF
            thfactor=((p0/current_state%surface_reference_pressure)**r_over_cp-(ptop/current_state%surface_reference_pressure)**&
                 r_over_cp)/((current_state%surface_pressure/current_state%surface_reference_pressure)**r_over_cp-&
                (ptop/current_state%surface_reference_pressure)**r_over_cp)
            do k=1,current_state%local_grid%size(Z_INDEX)
                current_state%global_grid%configuration%vertical%thref(k)=&
                     current_state%global_grid%configuration%vertical%thref(k)*thfactor
            end do
        end if
        if (ipass .eq. 1 .and. ANELASTIC_PROFILE_MODE .eq. 4) then                                    
            !             ! _adjust PTOP so that P0 = (fixed) PSF
            ptop=current_state%surface_pressure**r_over_cp+ptop**r_over_cp-p0**r_over_cp
            if (ptop .le. 0.0_DEFAULT_PRECISION .and. current_state%parallel%my_rank==0) then
                call log_log(LOG_ERROR, "Negative ptop in setup of anelastic. Need a warmer THREF or different setup options")
            end if
            ptop=ptop**(1.0_DEFAULT_PRECISION/r_over_cp)
        end if                                                                    
    end do
    !     !---------------------------------------
    !     ! _Finally compute density from pressure
    !     !---------------------------------------
    do k=1,current_state%local_grid%size(Z_INDEX)
        current_state%global_grid%configuration%vertical%rhon(k)=current_state%global_grid%configuration%vertical%prefn(k)&
            /(r*current_state%global_grid%configuration%vertical%thref(k)*&
            (current_state%global_grid%configuration%vertical%prefn(k)/current_state%surface_reference_pressure)**r_over_cp)
    end do
    do k=1,current_state%local_grid%size(Z_INDEX)-1
        current_state%global_grid%configuration%vertical%rho(k)=sqrt(current_state%global_grid%configuration%vertical%rhon(k)*&
             current_state%global_grid%configuration%vertical%rhon(k+1))                                           
    end do
    current_state%global_grid%configuration%vertical%rho(current_state%local_grid%size(Z_INDEX))=&
         current_state%global_grid%configuration%vertical%rhon(current_state%local_grid%size(Z_INDEX))**2/&
         current_state%global_grid%configuration%vertical%rhon(current_state%local_grid%size(Z_INDEX)-1)    
  end subroutine compute_anelastic_pressure_profile_and_density  


  subroutine check_top(zztop, z, info)
    real(kind=DEFAULT_PRECISION), intent(in) :: zztop
    real(kind=DEFAULT_PRECISION), intent(in) :: z
    character(*), intent(in) :: info

    if (z<zztop)then
      call log_master_log(LOG_ERROR, "Top of input profile is below the top of the domain: "//trim(info))
    end if

  end subroutine check_top

  subroutine check_input_levels(z_levels, field_levels, field)
    integer, intent(in) :: z_levels
    integer, intent(in) :: field_levels
    character(*), intent(in) :: field

    if (z_levels /= field_levels)then
       call log_master_log(LOG_ERROR, "Input levels not equal for "//trim(field)//", z_levels = "// &
            trim(conv_to_string(z_levels))//" field_levels = "//trim(conv_to_string(field_levels)))
    end if

  end subroutine check_input_levels
  
  subroutine set_qv_init_from_rh(current_state)

    type(model_state_type), intent(inout) :: current_state

    logical :: l_init_pl_rh    ! if .true. then initialize relative humidity field
    real(kind=DEFAULT_PRECISION) :: zztop ! top of the domain
    real(kind=DEFAULT_PRECISION) :: qsat
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: f_init_pl_rh     ! Initial node values for relative humidity variable
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: z_init_pl_rh     ! Initial node height values for relative humidity variable
    real(kind=DEFAULT_PRECISION), allocatable :: zgrid(:)  ! z grid to use in interpolation
    real(kind=DEFAULT_PRECISION), allocatable :: TdegK(:)  ! temperature in Kelvin
    integer :: i,j,n, k ! loop counters
    integer :: iq  ! temporary q varible index

    type(vertical_grid_configuration_type) :: vertical_grid

    vertical_grid=current_state%global_grid%configuration%vertical

    allocate(zgrid(current_state%local_grid%local_domain_end_index(Z_INDEX)))
    
    zztop = current_state%global_grid%top(Z_INDEX)

    l_init_pl_rh=options_get_logical(current_state%options_database, "l_init_pl_rh") 

    if (l_init_pl_rh)then
       allocate(z_init_pl_rh(options_get_array_size(current_state%options_database, "z_init_pl_rh")), &
             f_init_pl_rh(options_get_array_size(current_state%options_database, "f_init_pl_rh")))
       call options_get_real_array(current_state%options_database, "z_init_pl_rh", z_init_pl_rh)
       call options_get_real_array(current_state%options_database, "f_init_pl_rh", f_init_pl_rh)
       call check_top(zztop, z_init_pl_rh(size(z_init_pl_rh)), 'z_init_pl_rh')
       call check_input_levels(size(z_init_pl_rh), size(f_init_pl_rh), "f_init_pl_rh")
       zgrid=current_state%global_grid%configuration%vertical%zn(:)
       call piecewise_linear_1d(z_init_pl_rh(1:size(z_init_pl_rh)), f_init_pl_rh(1:size(f_init_pl_rh)), zgrid, &
            current_state%global_grid%configuration%vertical%rh_init)
      
      if (.not. current_state%passive_q .and. current_state%th%active) then
         iq=get_q_index('vapour', 'piecewise_initialization')
         allocate(TdegK(current_state%local_grid%local_domain_end_index(Z_INDEX)))
         TdegK(:) = current_state%global_grid%configuration%vertical%theta_init(:)* &
              (vertical_grid%prefn(:)/current_state%surface_reference_pressure)**r_over_cp
         do k = current_state%local_grid%local_domain_start_index(Z_INDEX), &
              current_state%local_grid%local_domain_end_index(Z_INDEX)
            qsat=qsaturation(TdegK(k), current_state%global_grid%configuration%vertical%prefn(k)/100.)    
            current_state%global_grid%configuration%vertical%q_init(k, iq) = & 
                 (current_state%global_grid%configuration%vertical%rh_init(k)/100.0)*qsat
         enddo
         if (.not. current_state%continuation_run) then
            do i=current_state%local_grid%local_domain_start_index(X_INDEX), &
                 current_state%local_grid%local_domain_end_index(X_INDEX)
               do j=current_state%local_grid%local_domain_start_index(Y_INDEX), &
                    current_state%local_grid%local_domain_end_index(Y_INDEX)
                  current_state%q(iq)%data(:,j,i) = current_state%global_grid%configuration%vertical%q_init(:, iq)
               end do
            end do
         end if

         deallocate(TdegK)
      else
         call log_master_log(LOG_ERROR, "Initialising with RH but q and/or theta passive")
      end if
      
      deallocate(z_init_pl_rh, f_init_pl_rh)
    end if

  end subroutine set_qv_init_from_rh


  
  subroutine read_dims(filename, ncid, nnodes, NNODES_KEY)
    character(*), intent(in) :: filename, NNODES_KEY
    integer, intent(in) :: ncid
    integer :: status, variable_id
    integer(kind=DEFAULT_PRECISION), intent(inout) :: nnodes

    status=nf90_inq_varid(ncid, NNODES_KEY, variable_id)
    if (status==nf90_noerr)then
      call check_status(nf90_inq_varid(ncid, NNODES_KEY, variable_id))
      call check_status(nf90_get_var(ncid, variable_id, nnodes))
    end if
  end subroutine read_dims



  subroutine read_variables(filename, ncid, nodes, elems, norms, NODES_KEY,&
    ELEMS_KEY, NORMS_KEY)
    character(*), intent(in) :: filename, NODES_KEY, ELEMS_KEY, NORMS_KEY
    integer, intent(in) :: ncid
    integer :: status, variable_id
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout) :: nodes
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout) :: norms
    integer(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout) :: elems
    status=nf90_inq_varid(ncid, NODES_KEY, variable_id)
    if (status==nf90_noerr)then
      call read_single_variable(ncid, NODES_KEY, data2d=nodes)
    end if

    status=nf90_inq_varid(ncid, ELEMS_KEY, variable_id)
    if (status==nf90_noerr)then
      call read_single_variable_int(ncid, ELEMS_KEY, data2d=elems)
    end if

    status=nf90_inq_varid(ncid, NORMS_KEY, variable_id)
    if (status==nf90_noerr)then
      call read_single_variable(ncid, NORMS_KEY, data2d=norms)
    end if

  end subroutine read_variables


  subroutine read_single_variable(ncid, key, data2d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout) :: data2d
    integer :: variable_id
    call check_status(nf90_inq_varid(ncid, key, variable_id))
    call check_status(nf90_get_var(ncid, variable_id, data2d))
  end subroutine read_single_variable


  subroutine read_single_variable_int(ncid, key, data2d)
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: key
    integer(kind=DEFAULT_PRECISION), dimension(:,:), intent(inout) :: data2d
    integer :: variable_id
    call check_status(nf90_inq_varid(ncid, key, variable_id))
    call check_status(nf90_get_var(ncid, variable_id, data2d))
  end subroutine read_single_variable_int


  !> Will check a NetCDF status and write to log_log error any decoded statuses
  !! @param status The NetCDF status flag
  subroutine check_status(status)
    integer, intent(in) :: status

    if (status /= nf90_noerr) then
      call log_log(LOG_ERROR, "NetCDF returned error code of "//trim(nf90_strerror(status)))
    end if
  end subroutine check_status



  subroutine find_image_point(idiag, nxg, nyg, i0, j0, xhalo, yhalo, myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi2, y_bi2, z_bi2, nrmi2, &
                              nrmj2, nrmk2, llen2)

    integer, intent(in) :: nxg, nyg
    integer, intent(in) :: i0, j0, xhalo, yhalo, myproc, nelem
    real(kind=DEFAULT_PRECISION), intent(in) :: x_g, y_g, z_g, dx, dy
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: nodes, norms
    integer(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: elems
    integer, dimension(:), intent(in) :: elem_lst
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: zgrid
    real(kind=DEFAULT_PRECISION), intent(out) :: x_ip, y_ip, z_ip
    real(kind=DEFAULT_PRECISION), intent(out) :: nrmi2, nrmj2, nrmk2, llen2
    integer, intent(inout) :: idiag

    integer :: i, iid, ielem, inode! loop counters
    integer :: i_surf, ip0, jp0, kp0
    integer :: ip1, ielem0, inode0, nfelem
    real(DEFAULT_PRECISION) :: diffx, diffy, diffz, z0
    real(DEFAULT_PRECISION) :: diffx2, diffy2, diffz2
    real(DEFAULT_PRECISION) :: x_bi, y_bi, z_bi, llen
    real(DEFAULT_PRECISION) :: x_bi2, y_bi2, z_bi2
    real(DEFAULT_PRECISION) :: tmp1, tmp2, tmp3, tmp4, tmp5
    real(DEFAULT_PRECISION) :: nrmi, nrmj, nrmk
    real(DEFAULT_PRECISION) :: xl, xr, yl, yr, zl, zr
    logical :: dwrite, l_found

    dwrite=.false.
    if (idiag .lt. 0)dwrite=.true.

    if(dwrite)then
      write(*,*)nelem
    end if

    llen = 1.e10
    i_surf = -999
    idiag = 0
    do i = 1, nelem
      ielem = elem_lst(i)
      diffx = nodes(elems(ielem,1),1) - x_g
      diffy = nodes(elems(ielem,1),2) - y_g
      if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
      if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
      if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
      if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 
      z0 = nodes(elems(ielem,1),3)
      nrmi = norms(ielem,1)
      nrmj = norms(ielem,2)
      nrmk = norms(ielem,3)
      tmp2 = diffx*nrmi + diffy*nrmj + (z0-z_g)*nrmk ! distance from grid point to BI
      x_bi = x_g + tmp2*nrmi ! boundary intercept coords
      y_bi = y_g + tmp2*nrmj
      z_bi = z_g + tmp2*nrmk

      ! check if intercept is within surface element limits
      xl = min(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
               nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      xr = max(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
               nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      yl = min(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
               nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      yr = max(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
               nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      zl = min(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
               nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      zr = max(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
               nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      if (abs(xl - x_bi) .gt. nxg*dx/2.0) xl = xl - (nxg*dx)
      if (abs(xr - x_bi) .gt. nxg*dx/2.0) xr = xr + (nxg*dx)
      if (abs(yl - y_bi) .gt. nyg*dy/2.0) yl = yl - (nyg*dy)
      if (abs(yr - y_bi) .gt. nyg*dy/2.0) yr = yr + (nyg*dy)
      if(dwrite.and.ielem.eq.17424)then
        write(*,*)i,ielem
      end if
      if (      x_bi .ge. xl .and. x_bi .le. xr &
          .and. y_bi .ge. yl .and. y_bi .le. yr &
          .and. z_bi .ge. zl .and. z_bi .le. zr) then
!        if (abs(tmp2) .lt. llen) then
        if (tmp2 .gt. 0.0 .and. tmp2 .lt. llen) then
          llen = tmp2
          i_surf = ielem
          x_ip = x_g + 2.*llen*nrmi ! image point coords
          y_ip = y_g + 2.*llen*nrmj
          z_ip = z_g + 2.*llen*nrmk
          x_bi2 = x_bi
          y_bi2 = y_bi
          z_bi2 = z_bi
          nrmi2 = nrmi
          nrmj2 = nrmj
          nrmk2 = nrmk
          llen2 = llen
          idiag=1
          if(dwrite)then
            write(*,*)'find_IP 1'
            write(*,*)ielem, tmp2, llen
            write(*,*)x_g, y_g, z_g
            write(*,*)diffx, diffy, (z0-z_g)
            write(*,*)nrmi, nrmj, nrmk
            write(*,*)x_bi, y_bi, z_bi
            write(*,*)xl, yl, zl
            write(*,*)xr, yr, zr
          end if ! diags

        end if
      end if
!      if(dwrite)then
!        write(*,*)'find_IP 1'
!        write(*,*)ielem, tmp2, llen
!        write(*,*)x_g, y_g, z_g
!        write(*,*)diffx, diffy, (z0-z_g)
!        write(*,*)nrmi, nrmj, nrmk
!        write(*,*)x_bi, y_bi, z_bi
!        write(*,*)xl, yl, zl
!        write(*,*)xr, yr, zr
!      end if

    end do ! ielem loop

    if (i_surf .eq. -999) then ! normal surface intercept not found. Find nearest edge normal:
      llen = 1.e10
      do i = 1, nelem
        ielem = elem_lst(i)
        ! check all edges associated with each element
        do inode = 1,4
          ip1 = inode+1
          if (ip1 .gt. 4)ip1=1
          diffx = nodes(elems(ielem,inode),1) - x_g
          diffy = nodes(elems(ielem,inode),2) - y_g
          diffz = nodes(elems(ielem,inode),3) - z_g
          if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
          if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
          if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
          if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 

          diffx2 = nodes(elems(ielem,ip1),1) - nodes(elems(ielem,inode),1)
          diffy2 = nodes(elems(ielem,ip1),2) - nodes(elems(ielem,inode),2)
          diffz2 = nodes(elems(ielem,ip1),3) - nodes(elems(ielem,inode),3)
          if (diffx2 .gt. nxg*dx/2.0) diffx2 = diffx2 - (nxg*dx) 
          if (diffx2 .lt. -nxg*dx/2.0) diffx2 = diffx2 + (nxg*dx) 
          if (diffy2 .gt. nyg*dy/2.0) diffy2 = diffy2 - (nyg*dy) 
          if (diffy2 .lt. -nyg*dy/2.0) diffy2 = diffy2 + (nyg*dy)
          tmp2 = (diffx2*diffx2 + diffy2*diffy2 + diffz2*diffz2)
          if(tmp2.gt.0.0)then
            tmp1 = -(diffx*diffx2 + diffy*diffy2 + diffz*diffz2)/tmp2
          else
            tmp1=-1.
          end if

          if (tmp1 .ge. 0.0 .and. tmp1 .lt. 1.0 ) then
            x_bi = nodes(elems(ielem,inode),1) + tmp1*diffx2
            y_bi = nodes(elems(ielem,inode),2) + tmp1*diffy2
            z_bi = nodes(elems(ielem,inode),3) + tmp1*diffz2
  
            tmp2 = sqrt((x_bi-x_g)*(x_bi-x_g) + &
                        (y_bi-y_g)*(y_bi-y_g) + &
                        (z_bi-z_g)*(z_bi-z_g)   )
            if (tmp2 .lt. llen .and. (z_bi-z_g) .gt. 0.0) then
              llen = tmp2
              i_surf = elems(ielem,inode)
              x_ip = x_g + 2.*(x_bi-x_g)! image point coords
              y_ip = y_g + 2.*(y_bi-y_g)
              z_ip = z_g + 2.*(z_bi-z_g)
              x_bi2 = x_bi
              y_bi2 = y_bi
              z_bi2 = z_bi
              nrmi2 = (x_bi-x_g)/llen
              nrmj2 = (y_bi-y_g)/llen
              nrmk2 = (z_bi-z_g)/llen
              llen2 = llen
              idiag=2
            end if
          end if

  
        end do ! inode loop
      end do ! ielem loop
    end if


    if (i_surf .eq. -999) then ! normal edge intercept not found. Use nearest node
      llen = 1.e10
      do i = 1, nelem
        ielem = elem_lst(i)
        if(dwrite)then
          write(*,*)'t3',i,ielem
        end if
        do inode = 1,4
          diffx = nodes(elems(ielem,inode),1) - x_g
          diffy = nodes(elems(ielem,inode),2) - y_g
          diffz = nodes(elems(ielem,inode),3) - z_g
          if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
          if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
          if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
          if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 

          tmp2 = sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
          if (tmp2 .lt. llen .and. diffz .gt. 0.0) then
            llen = tmp2
            i_surf = elems(ielem,inode)
            x_bi = nodes(elems(ielem,inode),1)
            y_bi = nodes(elems(ielem,inode),2)
            z_bi = nodes(elems(ielem,inode),3)
            x_ip = x_g + 2.*(x_bi-x_g)! image point coords
            y_ip = y_g + 2.*(y_bi-y_g)
            z_ip = z_g + 2.*(z_bi-z_g)
            x_bi2 = x_bi
            y_bi2 = y_bi
            z_bi2 = z_bi
            nrmi2 = (x_bi-x_g)/llen
            nrmj2 = (y_bi-y_g)/llen
            nrmk2 = (z_bi-z_g)/llen
            llen2 = llen
            idiag=3
          end if

        end do ! inode loop
      end do ! ielem loop
!      if (nrmk2 .lt. 0.0) then
!        write(*,*)'WARNING (gridmanager): negative k norm for image point at', myproc, x_g, y_g, z_g 
!      end if

    end if

  end subroutine find_image_point



  ! find_image_point_vr finds the face intercept point for the VR scheme
  subroutine find_image_point_vr(idiag, nxg, nyg, i0, j0, xhalo, yhalo, myproc, &
                              x_g, y_g, z_g, zgrid, dx, dy, dz, nodes, norms, elems, elem_lst, nelem, &
                              x_ip, y_ip, z_ip, x_bi2, y_bi2, z_bi2, nrmi2, &
                              nrmj2, nrmk2, llen2, ip_dir)

    integer, intent(in) :: nxg, nyg
    integer, intent(in) :: i0, j0, xhalo, yhalo, myproc, nelem
    real(kind=DEFAULT_PRECISION), intent(in) :: x_g, y_g, z_g, dx, dy, dz
    real(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: nodes, norms
    integer(kind=DEFAULT_PRECISION), dimension(:,:), intent(in) :: elems
    integer, dimension(:), intent(in) :: elem_lst
    real(kind=DEFAULT_PRECISION), dimension(:), intent(in) :: zgrid
    real(kind=DEFAULT_PRECISION), intent(out) :: x_ip, y_ip, z_ip
    real(kind=DEFAULT_PRECISION), intent(out) :: nrmi2, nrmj2, nrmk2, llen2
    integer, intent(inout) :: idiag
    integer, intent(out) :: ip_dir

    integer :: i, iid, ielem, inode! loop counters
    integer :: i_surf, ip0, jp0, kp0
    integer :: ip1, ielem0, inode0, nfelem
    real(DEFAULT_PRECISION) :: diffx, diffy, diffz, z0
    real(DEFAULT_PRECISION) :: diffx2, diffy2, diffz2
    real(DEFAULT_PRECISION) :: x_bi, y_bi, z_bi, llen
    real(DEFAULT_PRECISION) :: x_bi2, y_bi2, z_bi2
    real(DEFAULT_PRECISION) :: tmp1, tmp2, tmp3, tmp4, tmp5
    real(DEFAULT_PRECISION) :: nrmi, nrmj, nrmk
    real(DEFAULT_PRECISION) :: xl, xr, yl, yr, zl, zr
    logical :: dwrite, l_found

    dwrite=.false.
    if (idiag .lt. 0)dwrite=.true.


    llen = 1.e10
    i_surf = -999
    idiag = 0
    do i = 1, nelem
      ielem = elem_lst(i)

      diffx = nodes(elems(ielem,1),1) - x_g
      diffy = nodes(elems(ielem,1),2) - y_g
      if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
      if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
      if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
      if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 
      z0 = nodes(elems(ielem,1),3)
      nrmi = norms(ielem,1)
      nrmj = norms(ielem,2)
      nrmk = norms(ielem,3)
      tmp2 = diffx*nrmi + diffy*nrmj + (z0-z_g)*nrmk ! distance from grid point to BI
      x_bi = x_g + tmp2*nrmi ! boundary intercept coords
      y_bi = y_g + tmp2*nrmj
      z_bi = z_g + tmp2*nrmk

      ! check if intercept is within surface element limits
      xl = min(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
               nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      xr = max(nodes(elems(ielem,1),1),nodes(elems(ielem,2),1), &
               nodes(elems(ielem,4),1),nodes(elems(ielem,3),1) )
      yl = min(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
               nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      yr = max(nodes(elems(ielem,1),2),nodes(elems(ielem,2),2), &
               nodes(elems(ielem,4),2),nodes(elems(ielem,3),2) )
      zl = min(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
               nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      zr = max(nodes(elems(ielem,1),3),nodes(elems(ielem,2),3), &
               nodes(elems(ielem,4),3),nodes(elems(ielem,3),3) )
      if (abs(xl - x_bi) .gt. nxg*dx/2.0) xl = xl - (nxg*dx)
      if (abs(xr - x_bi) .gt. nxg*dx/2.0) xr = xr + (nxg*dx)
      if (abs(yl - y_bi) .gt. nyg*dy/2.0) yl = yl - (nyg*dy)
      if (abs(yr - y_bi) .gt. nyg*dy/2.0) yr = yr + (nyg*dy)
      if (      x_bi .ge. xl .and. x_bi .le. xr &
          .and. y_bi .ge. yl .and. y_bi .le. yr &
          .and. z_bi .ge. zl .and. z_bi .le. zr) then
        if (tmp2 .lt. 0.0 .and. abs(tmp2) .lt. llen) then
          llen = abs(tmp2)
          i_surf = ielem
          x_bi2 = x_bi
          y_bi2 = y_bi
          z_bi2 = z_bi
          nrmi2 = nrmi
          nrmj2 = nrmj
          nrmk2 = nrmk
          llen2 = llen
          idiag=1
!          if(dwrite)then
!            write(*,*)'find_IP'
!            write(*,*)ielem, tmp2, llen
!            write(*,*)x_g, y_g, z_g
!            write(*,*)diffx, diffy, (z0-z_g)
!            write(*,*)nrmi, nrmj, nrmk
!            write(*,*)x_bi, y_bi, z_bi
!            write(*,*)xl, yl, zl
!            write(*,*)xr, yr, zr
!          end if ! diags

        end if
      end if
      if(dwrite)then
        write(*,*)'find_IP'
        write(*,*)ielem, tmp2, llen
        write(*,*)x_g, y_g, z_g
        write(*,*)diffx, diffy, (z0-z_g)
        write(*,*)nrmi, nrmj, nrmk
        write(*,*)x_bi, y_bi, z_bi
        write(*,*)xl, yl, zl
        write(*,*)xr, yr, zr
      end if

    end do ! ielem loop

    if (i_surf .eq. -999) then ! normal surface intercept not found. Find nearest edge normal:
      llen = 1.e10
      do i = 1, nelem
        ielem = elem_lst(i)
        ! check all edges associated with each element
        do inode = 1,4
          ip1 = inode+1
          if (ip1 .gt. 4)ip1=1
          diffx = nodes(elems(ielem,inode),1) - x_g
          diffy = nodes(elems(ielem,inode),2) - y_g
          diffz = nodes(elems(ielem,inode),3) - z_g
          if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
          if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
          if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
          if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 

          diffx2 = nodes(elems(ielem,ip1),1) - nodes(elems(ielem,inode),1)
          diffy2 = nodes(elems(ielem,ip1),2) - nodes(elems(ielem,inode),2)
          diffz2 = nodes(elems(ielem,ip1),3) - nodes(elems(ielem,inode),3)
          if (diffx2 .gt. nxg*dx/2.0) diffx2 = diffx2 - (nxg*dx) 
          if (diffx2 .lt. -nxg*dx/2.0) diffx2 = diffx2 + (nxg*dx) 
          if (diffy2 .gt. nyg*dy/2.0) diffy2 = diffy2 - (nyg*dy) 
          if (diffy2 .lt. -nyg*dy/2.0) diffy2 = diffy2 + (nyg*dy)

          tmp2 = (diffx2*diffx2 + diffy2*diffy2 + diffz2*diffz2)
          if(tmp2.gt.0.0)then
            tmp1 = -(diffx*diffx2 + diffy*diffy2 + diffz*diffz2)/tmp2
          else
            tmp1=-1.
          end if

          if (tmp1 .ge. 0.0 .and. tmp1 .lt. 1.0 ) then
            x_bi = nodes(elems(ielem,inode),1) + tmp1*diffx2
            y_bi = nodes(elems(ielem,inode),2) + tmp1*diffy2
            z_bi = nodes(elems(ielem,inode),3) + tmp1*diffz2
  
            tmp2 = sqrt((x_bi-x_g)*(x_bi-x_g) + &
                        (y_bi-y_g)*(y_bi-y_g) + &
                        (z_bi-z_g)*(z_bi-z_g)   )
            if (tmp2 .lt. llen .and. (z_g-z_bi) .gt. 1.e-9) then
              llen = tmp2
              i_surf = elems(ielem,inode)
              x_bi2 = x_bi
              y_bi2 = y_bi
              z_bi2 = z_bi
              nrmi2 = -(x_bi-x_g)/llen
              nrmj2 = -(y_bi-y_g)/llen
              nrmk2 = -(z_bi-z_g)/llen
              llen2 = llen
              idiag=2
              if(dwrite)then
                write(*,*)ielem, inode, ip1
                write(*,*)x_g, y_g, z_g
                write(*,*)diffx, diffy, diffz
                write(*,*)diffx2, diffy2, diffz2
                write(*,*)x_bi, y_bi, z_bi
                write(*,*)tmp1, tmp2, llen
              end if
!
            end if
          end if
  
        end do ! inode loop
      end do ! ielem loop
    end if


    if (i_surf .eq. -999) then ! normal edge intercept not found. Use nearest node
      llen = 1.e10
      do i = 1, nelem
        ielem = elem_lst(i)
        do inode = 1,4
          diffx = nodes(elems(ielem,inode),1) - x_g
          diffy = nodes(elems(ielem,inode),2) - y_g
          diffz = nodes(elems(ielem,inode),3) - z_g
          if (diffx .gt. nxg*dx/2.0) diffx = diffx - (nxg*dx) 
          if (diffx .lt. -nxg*dx/2.0) diffx = diffx + (nxg*dx) 
          if (diffy .gt. nyg*dy/2.0) diffy = diffy - (nyg*dy) 
          if (diffy .lt. -nyg*dy/2.0) diffy = diffy + (nyg*dy) 

          tmp2 = sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
          if (tmp2 .lt. llen .and. diffz .lt. 1.e-9) then
            llen = tmp2
            i_surf = elems(ielem,inode)
            x_bi = nodes(elems(ielem,inode),1)
            y_bi = nodes(elems(ielem,inode),2)
            z_bi = nodes(elems(ielem,inode),3)
            x_bi2 = x_bi
            y_bi2 = y_bi
            z_bi2 = z_bi
            nrmi2 = -(x_bi-x_g)/llen
            nrmj2 = -(y_bi-y_g)/llen
            nrmk2 = -(z_bi-z_g)/llen
            llen2 = llen
            idiag=3
          end if

        end do ! inode loop
      end do ! ielem loop

    end if

    x_ip = 0. 
    y_ip = 0.
    z_ip = 0.
    if (idiag .gt. 0) then ! calc image point
      tmp1 = 1.e10
      ip_dir=0
      if (abs(nrmi2).gt.0.0.and.abs(dx/nrmi2).lt.tmp1) then
        tmp1=abs(dx/nrmi2)
        ip_dir=1
      end if
      if (abs(nrmj2).gt.0.0.and.abs(dy/nrmj2).lt.tmp1) then
        tmp1=abs(dy/nrmj2)
        ip_dir=2
      end if
      if (abs(nrmk2).gt.0.0.and.abs(dz/nrmk2).lt.tmp1) then
        tmp1=abs(dz/nrmk2)
        ip_dir=3
      end if

      ! move IP into next grid cell for interpolation
      if(ip_dir.eq.1)tmp1=tmp1+(dx/2.0)
      if(ip_dir.eq.2)tmp1=tmp1+(dy/2.0)
      if(ip_dir.eq.3)tmp1=tmp1+(dz/2.0)

      x_ip = x_g + nrmi2*tmp1
      y_ip = y_g + nrmj2*tmp1
      z_ip = z_g + nrmk2*tmp1

    end if
    if(dwrite)then
      write(*,*)idiag
      write(*,*)x_ip,y_ip,z_ip
      write(*,*)nrmi2, nrmj2, nrmk2
      write(*,*)ip_dir, tmp1
    end if ! diags

  end subroutine find_image_point_vr




end module gridmanager_mod
