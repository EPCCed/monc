! Apply Gaussian filter in vertical (part of immersed boundary code)

module vert_filter_mod

  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type
  use grids_mod
  use optionsdatabase_mod, only : options_get_logical, options_get_integer, options_get_real
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, PRECISION_TYPE
  use state_mod, only : model_state_type

  implicit none

#ifndef TEST_MODE
  private
#endif

  logical :: l_dofilt
  logical :: u_filter_enabled     ! activate filter
  logical :: v_filter_enabled     ! activate filter
  logical :: w_filter_enabled     ! activate filter
  logical :: th_filter_enabled     ! activate filter
  real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: work1  ! shifted column for conv
  real(kind=DEFAULT_PRECISION), dimension(3),parameter :: filt=(/0.25, 0.5, 0.25/)  ! filter kernel
  integer :: nvert  ! number of points above IB to filter
  integer :: nk, kstart, kstop, nfilt
  

  
  public vert_filter_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function vert_filter_get_descriptor()
    vert_filter_get_descriptor%name="vert_filter"
    vert_filter_get_descriptor%version=0.1
    vert_filter_get_descriptor%timestep=>timestep_callback
    vert_filter_get_descriptor%initialisation=>initialisation_callback
    vert_filter_get_descriptor%finalisation=>finalisation_callback

  end function vert_filter_get_descriptor
  


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    real(kind=DEFAULT_PRECISION) :: tmp2
    real(kind=DEFAULT_PRECISION), dimension(:), allocatable :: tmp
    integer :: i

    l_dofilt=options_get_logical(current_state%options_database, "vert_filter_enabled")
    if (.not.current_state%immersed%ib_enabled) then
      l_dofilt=.false.
    end if

    if (l_dofilt) then
      u_filter_enabled=options_get_logical(current_state%options_database, "u_filter_enabled")
      v_filter_enabled=options_get_logical(current_state%options_database, "v_filter_enabled")
      w_filter_enabled=options_get_logical(current_state%options_database, "w_filter_enabled")
      th_filter_enabled=options_get_logical(current_state%options_database, "th_filter_enabled")
      nfilt=options_get_integer(current_state%options_database, "nfilt")
      nvert=options_get_integer(current_state%options_database, "nvert")

      ! column points above GP to process
      allocate(work1(nvert))
    end if

  end subroutine initialisation_callback



  !> Timestep callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    integer :: i,j,k,p,q, myproc

    if (l_dofilt) then
!      myproc = current_state%parallel%my_rank

      if (current_state%halo_column) return
      j=current_state%column_local_y
      i=current_state%column_local_x

      do k=1,current_state%immersed%kmax_ji(j,i)

        ! u points
        if (u_filter_enabled.and.(current_state%immersed%ghost_u(k,j,i).gt.0)) then
          do q=1,nfilt 
            do p=1,nvert
              work1(p)=current_state%u%data(k+p-1,j,i)*filt(1)+&
                       current_state%u%data(k+p,j,i)*filt(2)+&
                       current_state%u%data(k+p+1,j,i)*filt(3)
            end do
            do p=1,nvert
              current_state%u%data(k+p,j,i) = work1(p)
            end do
          end do
        end if
  
        ! v points
        if (v_filter_enabled.and.(current_state%immersed%ghost_v(k,j,i).gt.0)) then
          do q=1,nfilt 
            do p=1,nvert
              work1(p)=current_state%v%data(k+p-1,j,i)*filt(1)+&
                       current_state%v%data(k+p,j,i)*filt(2)+&
                       current_state%v%data(k+p+1,j,i)*filt(3)
            end do
            do p=1,nvert
              current_state%v%data(k+p,j,i) = work1(p)              
            end do
          end do
        end if
  
        ! w points
        if (w_filter_enabled.and.(current_state%immersed%ghost_w(k,j,i).gt.0)) then
          do q=1,nfilt 
            do p=1,nvert
              work1(p)=current_state%w%data(k+p-1,j,i)*filt(1)+&
                       current_state%w%data(k+p,j,i)*filt(2)+&
                       current_state%w%data(k+p+1,j,i)*filt(3)
            end do
            do p=1,nvert
              current_state%w%data(k+p,j,i) = work1(p)
            end do
          end do
        end if
  
        ! theta points
        if (th_filter_enabled.and.(current_state%immersed%ghost_s(k,j,i).gt.0)) then
          do q=1,nfilt 
            do p=1,nvert
              work1(p)=current_state%th%data(k+p-1,j,i)*filt(1)+&
                       current_state%th%data(k+p,j,i)*filt(2)+&
                       current_state%th%data(k+p+1,j,i)*filt(3)
            end do
            do p=1,nvert
              current_state%th%data(k+p,j,i) = work1(p)
            end do
          end do
        end if
  
      end do



    end if ! end l_dofilt


  end subroutine timestep_callback

  !> Finalisation callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
  end subroutine finalisation_callback



end module vert_filter_mod
