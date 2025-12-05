!> Functionality to support the different types of grid and abstraction between global grids
!! and local ones after decomposition
!!
!! Currently MONC supports the Arakawa C grid
module grids_mod
  use datadefn_mod, only : DEFAULT_PRECISION, SINGLE_PRECISION
  implicit none

#ifndef TEST_MODE
  private
#endif

  !> Grid index parameters
  integer, parameter, public :: Z_INDEX = 1, Y_INDEX = 2, X_INDEX = 3
  !> Grid type parameters (usually applied to each dimension of a prognostic)
  integer, parameter, public :: PRIMAL_GRID=1, DUAL_GRID=2

  !> The configuration of the grid horizontally
  type, public :: horizontal_grid_configuration_type
     real(kind=DEFAULT_PRECISION) dx, dy     !< Grid spacings in x and y directions (m)
     real(kind=DEFAULT_PRECISION) cx, cy     !< Reciprocal of the horizontal grid spacings
     real(kind=DEFAULT_PRECISION) cx2, cy2   !< Reciprocal of the square of the grid spacings
     real(kind=DEFAULT_PRECISION) cxy        !< 1/(DX*DY)
     real(kind=DEFAULT_PRECISION) tcx, tcy
  end type horizontal_grid_configuration_type

  type, public :: immersed_boundary_type
    ! Immersed boundary arrays
    logical :: ib_enabled

    ! diagnostics
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_u  ! u on boundary points
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_v  ! v on boundary points
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_w  ! w on boundary points
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_th ! th on boundary points
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_p ! p on boundary points
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ib_ustar ! ustar on k2 equivalent
    logical :: dump_u, dump_v, dump_w, dump_th, dump_p, dump_ustar, dodiags
    logical :: diags_enabled
    real(kind=DEFAULT_PRECISION) :: dump_freq
 
    integer, dimension(:,:), allocatable :: kmax_ji
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: surface_temperature_flux_ib_s, surface_temperature_flux_ib_w
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: surface_vapour_flux_ib_s, surface_vapour_flux_ib_w
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: theta_surf_ib_s, theta_surf_ib_w, thref_surf_ib_s, thref_surf_ib_w
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: theta_virtual_surf_ib_s, theta_virtual_surf_ib_w
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: surface_vapour_mixing_ratio_ib_s, surface_vapour_mixing_ratio_ib_w
    integer, dimension(:,:,:), allocatable :: indic_u, indic_v, indic_w, indic_s ! fluid / solid indicator
    logical, dimension(:,:), allocatable :: ib_col_s, ib_col_u, ib_col_v, ib_col_w, ib_col
    integer, dimension(:,:,:), allocatable :: ghost_u, ghost_v, ghost_w, ghost_s ! ghost point indicator
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: rneutml_sq ! neutral mixing length wrt immersed boundary
    real(kind=DEFAULT_PRECISION), dimension(:,:,:), allocatable :: ibprox ! immersed boundary proximity function for spectral filter
    integer, dimension(:,:), allocatable :: gp_ijk_u, gp_ijk_v, gp_ijk_w, gp_ijk_s ! ghost point indices
    integer, dimension(:,:), allocatable :: ip_ijk_u, ip_ijk_v, ip_ijk_w, ip_ijk_s ! image point base indices
    integer, dimension(:,:), allocatable :: bi_ijk_u, bi_ijk_v, bi_ijk_w, bi_ijk_s ! boundary intercept base indices
    integer, dimension(:,:), allocatable :: ip_ijk_k2e_s, ip_ijk_k2e_su, ip_ijk_k2e_sv, ip_ijk_k2e_sw ! image point base indices for k2 equivalent height from IB
    integer, dimension(:,:), allocatable :: ip_ijk_k2e_w, ip_ijk_k2e_ws, ip_ijk_k2e_wu, ip_ijk_k2e_wv! image point base indices for k2 equivalent height from IB
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_u, int_xyz_v, int_xyz_w, int_xyz_s ! trilinear interpolation fractional distances (xyz)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: biint_xyz_u, biint_xyz_v, biint_xyz_w, biint_xyz_s ! trilinear interpolation fractional distances for boundary intercept (xyz)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_k2e_w, int_xyz_k2e_ws, int_xyz_k2e_wu, int_xyz_k2e_wv! trilinear interpolation fractional distances (xyz)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: bi_xyz_u, bi_xyz_v, bi_xyz_w, bi_xyz_s, k2e_xyz_w ! boundary intercept locations
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: norms_s, norms_w, norms_u, norms_v
    integer, dimension(:), allocatable :: idiag_u, idiag_v, idiag_w, idiag_s
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: gp_store_u, gp_store_v, gp_store_w, gp_store_th
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: gp_store_vis, gp_store_diff
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: gp_store_q

    integer, dimension(:,:), allocatable :: ip_ijk_uu, ip_ijk_uv, ip_ijk_uw ! VR image point base indices
    integer, dimension(:,:), allocatable :: ip_ijk_vu, ip_ijk_vv, ip_ijk_vw 
    integer, dimension(:,:), allocatable :: ip_ijk_wu, ip_ijk_wv, ip_ijk_ww,  ip_ijk_ws
    integer, dimension(:,:), allocatable :: ip_ijk_ss,  ip_ijk_sw
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_uu, int_xyz_uv, int_xyz_uw ! VR trilinear interpolation fractional distances (xyz)
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_vu, int_xyz_vv, int_xyz_vw 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_wu, int_xyz_wv, int_xyz_ww, int_xyz_ws 
    real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: int_xyz_ss, int_xyz_sw

    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: gp_len_u, gp_len_v, gp_len_w, gp_len_s, k2e_len
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: ip_len_u, ip_len_v, ip_len_w, ip_len_s

    ! For specified flux BCs, scalars are interpolated from nearby BIs on the w grid
    integer, dimension(:,:), allocatable :: w2s_idw_idx ! wghost indices
    real, dimension(:,:), allocatable :: w2s_idw_wgt ! w point weights 


    integer :: ib_type ! no slip, velocity reconstruction etc.
  
  end type immersed_boundary_type


  !> The configuration of the grid vertically
  type, public :: vertical_grid_configuration_type
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: &
          z,&      !< Heights at w levels (m)
          zn,&     !< Heights at pressure levels (m)
          dz,&     !< Vertical spacing between K and K-1 w levels
          dzn,&    !< Vertical spacing between K and K-1 p levels
          rdz,&    !< Reciprocal of DZ
          rdzn,&   !< Reciprocal of DZN
          rho,&    !< Density at w levels (kg/m3)
          rhon,&   !< Density at p levels (kg/m3)
          thref,&  !< Reference potential temperature (K)
          dthref,& !< Gradient of thref (K)
          tref,&   !< Reference temperature (K)
          theta_init,& !< Initial profile of potential temperature
          temp_init,&  !< Initial profile of absolute temperature
          rh_init, &   !< Initial profile of relative humidity
          u_init,&     !< Initial profile of u
          v_init,&     !< Initial profile of v
          prefn,& !< Reference pressure (Pa)
          pdiff,& !< Difference between pressure levels (Pa)
          prefrcp,&
          rprefrcp,&
          czb,&   !< (rho(k-1)/rhon(k))/(dz(k)*dzn(k)) use for diffusion onto p level from below
          cza,&   !< (rho(k)/rhon(k))/(dz(k)*dzn(k+1)) use for diffusion onto p level from above
          czg,&   !< CZB-CZA for tridiagonal solver in POISSON
          czh,&   !< CZB*CZA for tridiagonal solver in POISSON      
          tzc1,&  !< 0.25*rdz(k)*rho(k-1)/rhon(k) for advection onto p-level from below
          tzc2,&  !< 0.25*rdz(k)*rho(k)/rhon(k) for advection onto p-level from above
          tzd1,&  !< 0.25*rdzn(k+1)*rhon(k)/rho(k) for advection onto w-level from below
          tzd2,&  !< 0.25*rdzn(k+1)*rhon(k+1)/rho(k) for advection onto w-level from above 
          w_subs,& !< Subsidence velocity
          olubar,& !< Current U mean
          savolubar,&          
          olvbar,& !< Current V mean
          savolvbar,&
          olthbar,& !< Current theta mean
          olzubar,& !< Previous timestep U mean
          olzvbar,& !< Previous timestep V mean
          olzthbar,& !< Previous theta mean
          dmpco,& !< Damping coefficient on pressure levels
          dmpcoz,& !< Damping coefficient on w-levels
          tstarpr,& ! Temperature about which Taylor Expansion
          qsat,&
          dqsatdt,&
          qsatfac,&
          rneutml,&
          rneutml_sq,&
          buoy_co, &
          theta_rand, &!< profile of amplitude of theta perturbations
          theta_force, & !<profile of forcing term for theta
          u_force, & !<profile of forcing term for u
          v_force, & !<profile of forcing term for v
          w_up,    & !<profile of updraft threshold
          w_dwn,   & !<profile of downdraft threshold
          w_rand     !<profile of amplitude of w perturbations



     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_init !< Initial profile of q variables
     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_rand !< Initial profile of amplitude of q perturbations
     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: q_force !< Profiles of forcing terms for q variables

     real(kind=DEFAULT_PRECISION), dimension(:,:), allocatable :: olqbar,olzqbar
     
     real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: hgd
     real(kind=DEFAULT_PRECISION) :: czn, zlogm, zlogth, vk_on_zlogm
     integer, dimension(:), allocatable :: kgd
     integer :: kdmpmin
  end type vertical_grid_configuration_type

  !> Wraps the dimensional configuration types
  type, public :: grid_configuration_type
     type(horizontal_grid_configuration_type) :: horizontal !< Configuration horizontally
     type(vertical_grid_configuration_type) :: vertical !< Configuration vertically
  end type grid_configuration_type

  !> Defines the global grid
  type, public :: global_grid_type
     type(grid_configuration_type) :: configuration !< Configuration of the grid
     real(kind=DEFAULT_PRECISION), dimension(3) :: bottom,&     !< Bottom (lowest) bounds in each dimension
          top,&        !< Top (highest) bounds in each dimension
          resolution   !< The resolution of the grid in each dimension
     integer, dimension(3) ::    size         !< Number of grid points in each dimension
     logical, dimension(3) :: active = (/ .false., .false., .false. /) !< Whether a specific dimension is active
     integer :: dimensions = 0 !< Number of active dimensions
  end type global_grid_type

  !> Defined the local grid, i.e. the grid held on this process after decomposition
  type, public :: local_grid_type
     integer, dimension(3) :: start,& !< Global start coordinates of local grid
          end,& !< Global end coordinates of local grid
          size,& !< Grid points held locally
          halo_size,& !< Grid point size of the halo (halo_depth)
          local_domain_start_index,& !< The start index of computation data (local data is halo->data->halo so this precomputes the data start)
          local_domain_end_index !< The end index of computation data (local data is halo->data->halo so this precomputes the data end)
     logical, dimension(3) :: active = (/ .false., .false., .false. /) !< Dimensions which are active
     integer, dimension(:,:), allocatable :: neighbours , corner_neighbours !< Neighbouring process Id per dimension
     integer :: dimensions = 0 !< Number of active dimensions
  end type local_grid_type

end module grids_mod
