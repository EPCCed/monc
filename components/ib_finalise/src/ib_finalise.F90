! Finalise immersed boundary variables, get diagnostics, and zero sources on
! solid points

module ib_finalise_mod

  use monc_component_mod, only : COMPONENT_SCALAR_FIELD_TYPE, COMPONENT_DOUBLE_DATA_TYPE, &
       COMPONENT_ARRAY_FIELD_TYPE, COMPONENT_INTEGER_DATA_TYPE, component_descriptor_type, &
       component_field_value_type, component_field_information_type

  use state_mod, only : FORWARD_STEPPING, PRESCRIBED_SURFACE_FLUXES, PRESCRIBED_SURFACE_VALUES, &
   model_state_type
  use grids_mod
  use optionsdatabase_mod, only : options_get_real, options_get_integer, options_get_array_size, &
     options_get_real_array, options_get_string, options_get_logical
  use datadefn_mod, only : DEFAULT_PRECISION, STRING_LENGTH, PRECISION_TYPE
  use conversions_mod, only : conv_to_string

  implicit none

#ifndef TEST_MODE
  private
#endif

  integer :: dump_number, didi, didj, iqv
  character (len=50):: dumpdir
  logical :: l_doset
!  real(kind=DEFAULT_PRECISION), parameter :: t_rlx=1./100. ! relaxation time for subterrainian points
  real(kind=DEFAULT_PRECISION) :: t_rlx ! relaxation time for subterrainian points
  logical :: enable_theta=.false.
  real(kind=DEFAULT_PRECISION) :: viscous_courant_coefficient
  public ib_finalise_get_descriptor

contains

  !> Provides a description of this component for the core to register
  !! @returns The descriptor containing registration information for this component
  type(component_descriptor_type) function ib_finalise_get_descriptor()
    ib_finalise_get_descriptor%name="ib_finalise"
    ib_finalise_get_descriptor%version=0.1
    ib_finalise_get_descriptor%timestep=>timestep_callback
    ib_finalise_get_descriptor%initialisation=>initialisation_callback
    ib_finalise_get_descriptor%finalisation=>finalisation_callback

  end function ib_finalise_get_descriptor
  


  subroutine initialisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
    integer ::  myproc,i,j,k
    character (len=50):: fname

    if(.not.current_state%immersed%ib_enabled)return

    enable_theta=options_get_logical(current_state%options_database, "enable_theta")
    
    t_rlx=1./50.
!    if (current_state%immersed%ib_type==1)t_rlx = 1.

    if(current_state%number_q_fields > 0)then
      iqv = current_state%water_vapour_mixing_ratio_index
    end if

    if (current_state%immersed%diags_enabled) then 
      dumpdir=options_get_string(current_state%options_database, "ib_dump_dir")
      dump_number = 0
  
      ! column index for dump switch
      didj=current_state%local_grid%size(Y_INDEX) + current_state%local_grid%halo_size(Y_INDEX)
      didi=current_state%local_grid%size(X_INDEX) + current_state%local_grid%halo_size(X_INDEX)
  
      myproc = current_state%parallel%my_rank
  
      if (current_state%immersed%dump_u) then
        write(fname,'(a,a,i4.4,a)') trim(dumpdir),"/ib_u_",myproc,"_xyz"
        open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
        write(20)current_state%immersed%bi_xyz_u
        close(20)
      end if
      if (current_state%immersed%dump_v) then
        write(fname,'(a,a,i4.4,a)') trim(dumpdir),"/ib_v_",myproc,"_xyz"
        open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
        write(20)current_state%immersed%bi_xyz_v
        close(20)
      end if
      if (current_state%immersed%dump_w) then
        write(fname,'(a,a,i4.4,a)') trim(dumpdir),"/ib_w_",myproc,"_xyz"
        open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
        write(20)current_state%immersed%bi_xyz_w
        close(20)
      end if
      if (current_state%immersed%dump_th .or. current_state%immersed%dump_p) then
        write(fname,'(a,a,i4.4,a)') trim(dumpdir),"/ib_s_",myproc,"_xyz"
        open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
        write(20)current_state%immersed%bi_xyz_s
        close(20)
      end if
      if (current_state%immersed%dump_ustar) then
        write(fname,'(a,a,i4.4,a)') trim(dumpdir),"/ib_ustar_",myproc,"_xyz"
        open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
        write(20)current_state%immersed%k2e_xyz_w
        close(20)
      end if


    end if

    ! extra diagnostics...
    if (.true.) then
      dumpdir=options_get_string(current_state%options_database, "ib_dump_dir")
      myproc = current_state%parallel%my_rank
      write(fname,'(a,a,i4.4)') trim(dumpdir),"/indic_u_",myproc
      open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
      write(20)current_state%immersed%indic_u
      close(20)
      write(fname,'(a,a,i4.4)') trim(dumpdir),"/indic_v_",myproc
      open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
      write(20)current_state%immersed%indic_v
      close(20)
      write(fname,'(a,a,i4.4)') trim(dumpdir),"/indic_w_",myproc
      open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
      write(20)current_state%immersed%indic_w
      close(20)
      write(fname,'(a,a,i4.4)') trim(dumpdir),"/indic_s_",myproc
      open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
      write(20)current_state%immersed%indic_s
      close(20)
    end if

    do i=1,current_state%local_grid%size(X_INDEX)
      do j=1,current_state%local_grid%size(Y_INDEX)
        do k=1,current_state%local_grid%size(Z_INDEX)
          if (current_state%immersed%indic_u(k,j,i) .eq. 1) then
              current_state%u%data(k,j,i) = 0.0_DEFAULT_PRECISION
          end if
          if (current_state%immersed%indic_v(k,j,i) .eq. 1) then
              current_state%v%data(k,j,i) = 0.0_DEFAULT_PRECISION
          end if
          if (current_state%immersed%indic_w(k,j,i) .eq. 1) then
              current_state%w%data(k,j,i) = 0.0_DEFAULT_PRECISION
          end if
          if (current_state%immersed%indic_s(k,j,i) .eq. 1) then
              current_state%th%data(k,j,i) = 0.0_DEFAULT_PRECISION
          end if
        end do
      end do
    end do

    ! Initialise arrays for solid points

    l_doset=current_state%immersed%ib_enabled

    ! use non-IB geometry for now
    viscous_courant_coefficient=1.0_DEFAULT_PRECISION/current_state%global_grid%configuration%vertical%dzn(2)**2+&
         1.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dx*&
         current_state%global_grid%configuration%horizontal%dx)+&
         1.0_DEFAULT_PRECISION/(current_state%global_grid%configuration%horizontal%dy*&
         current_state%global_grid%configuration%horizontal%dy)

  end subroutine initialisation_callback




  !> Timestep callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine timestep_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state

    integer :: i,j,k, ighost, myproc
    integer :: dcheck, ierr
    character (len=50):: fname
    logical :: dodump 
    integer :: ipi, ipj, ipk, nramp=10
    real(kind=DEFAULT_PRECISION) :: dtm_scale, tmp
          
    if(.not.current_state%immersed%ib_enabled)return
    myproc = current_state%parallel%my_rank
    
    if (current_state%immersed%ib_enabled .and.&
        current_state%immersed%diags_enabled) then 
  
      j=current_state%column_local_y
      i=current_state%column_local_x
 

      dodump  = .false.
      if (current_state%immersed%dodiags .and. i .eq. didi .and. j .eq. didj) then
        dcheck = nint(current_state%time/current_state%immersed%dump_freq)
        if (dcheck .ne. dump_number) then
          dodump = .true.
          dump_number = dump_number + 1

 
          if (current_state%immersed%dump_u) then
            write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_u_",myproc,"_",trim(conv_to_string(int(current_state%time)))
            open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
              write(20) current_state%immersed%ib_u
            close(20)
          end if
  
          if (current_state%immersed%dump_v) then
            write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_v_",myproc,"_",trim(conv_to_string(int(current_state%time)))
            open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
              write(20) current_state%immersed%ib_v
            close(20)
          end if
  
          if (current_state%immersed%dump_w) then
            write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_w_",myproc,"_",trim(conv_to_string(int(current_state%time)))
            open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
              write(20) current_state%immersed%ib_w
            close(20)
          end if
  
          if (current_state%immersed%dump_p) then
            write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_p_",myproc,"_",trim(conv_to_string(int(current_state%time)))
            open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
              write(20) current_state%immersed%ib_p
            close(20)
          end if
  
          if(enable_theta)then
            if (current_state%immersed%dump_th) then
              write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_th_",myproc,"_",trim(conv_to_string(int(current_state%time)))
              open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
                write(20) current_state%immersed%ib_th
              close(20)
            end if
          end if
    
          if (current_state%immersed%dump_ustar) then
            write(fname,'(a,a,i4.4,a,a)') trim(dumpdir),"/ib_ustar_",myproc,"_",trim(conv_to_string(int(current_state%time)))
            open(20,file=trim(fname), form='unformatted', status='replace', access='stream')
              write(20) current_state%immersed%ib_ustar
            close(20)
          end if

        end if ! end dcheck
      end if ! end dodiags and i,j check
    end if ! end ib and diags enabled
 



    !***********************************************************************
    ! set variables for GPs, zero source terms and relax non GPs to default.

    if (current_state%halo_column) return
     
    if (l_doset) then 
      j=current_state%column_local_y
      i=current_state%column_local_x

      if (current_state%immersed%ib_col(j,i)) then

        if (current_state%immersed%ib_type==0)then

          do k=1,current_state%immersed%kmax_ji(j,i)
            if (current_state%immersed%ghost_u(k,j,i) .gt. 0) then ! u points
              ighost = current_state%immersed%ghost_u(k,j,i)
              current_state%u%data(k,j,i) = current_state%immersed%gp_store_u(ighost)
              current_state%su%data(k,j,i) = 0.0_DEFAULT_PRECISION
            end if
            if (current_state%immersed%ghost_v(k,j,i) .gt. 0) then ! v points
              ighost = current_state%immersed%ghost_v(k,j,i)
              current_state%v%data(k,j,i) = current_state%immersed%gp_store_v(ighost)
              current_state%sv%data(k,j,i) = 0.0_DEFAULT_PRECISION
            end if
            if (current_state%immersed%ghost_w(k,j,i) .gt. 0) then ! w points
              ighost = current_state%immersed%ghost_w(k,j,i)
              current_state%w%data(k,j,i) = current_state%immersed%gp_store_w(ighost)
              current_state%sw%data(k,j,i) = 0.0_DEFAULT_PRECISION
              if(current_state%use_viscosity_and_diffusion)then
                current_state%vis_coefficient%data(k,j,i) = current_state%immersed%gp_store_vis(ighost)
                current_state%diff_coefficient%data(k,j,i) = current_state%immersed%gp_store_diff(ighost)
              end if
            end if
            if(enable_theta)then
              if (current_state%immersed%ghost_s(k,j,i) .gt. 0) then ! s points
                ighost = current_state%immersed%ghost_s(k,j,i)
                current_state%th%data(k,j,i) = current_state%immersed%gp_store_th(ighost)
                current_state%sth%data(k,j,i) = 0.0_DEFAULT_PRECISION
              end if
            end if
            if(current_state%number_q_fields .gt. 0)then
              if (current_state%immersed%ghost_s(k,j,i) .gt. 0) then ! s points
                ighost = current_state%immersed%ghost_s(k,j,i)
                current_state%q(iqv)%data(k,j,i) = current_state%immersed%gp_store_q(ighost)
                current_state%sq(iqv)%data(k,j,i) = 0.0_DEFAULT_PRECISION
              end if
            end if
  
            if (current_state%immersed%indic_u(k,j,i) .eq. 1) then
              current_state%su%data(k,j,i) = 0.0_DEFAULT_PRECISION
              if (current_state%immersed%ghost_u(k,j,i) .lt. 1) then
                current_state%su%data(k,j,i) = t_rlx*&
                (0.0_DEFAULT_PRECISION - current_state%u%data(k,j,i))
              end if
            end if
            if (current_state%immersed%indic_v(k,j,i) .eq. 1) then
              current_state%sv%data(k,j,i) = 0.0_DEFAULT_PRECISION
              if (current_state%immersed%ghost_v(k,j,i) .lt. 1) then
                current_state%sv%data(k,j,i) = t_rlx*&
                (0.0_DEFAULT_PRECISION - current_state%v%data(k,j,i))
              end if
            end if
            if (current_state%immersed%indic_w(k,j,i) .eq. 1) then
              current_state%sw%data(k,j,i) = 0.0_DEFAULT_PRECISION
              if (current_state%immersed%ghost_w(k,j,i) .lt. 1) then
                current_state%sw%data(k,j,i) = t_rlx*&
                (0.0_DEFAULT_PRECISION - current_state%w%data(k,j,i))
                current_state%vis_coefficient%data(k,j,i) = current_state%vis_coefficient%data(k+1,j,i)
                current_state%diff_coefficient%data(k,j,i) = current_state%diff_coefficient%data(k+1,j,i)
              end if
            end if
  
            if(enable_theta)then
              if (current_state%immersed%indic_s(k,j,i) .eq. 1) then
                current_state%sth%data(k,j,i) = 0.0_DEFAULT_PRECISION
                if (current_state%immersed%ghost_s(k,j,i) .lt. 1) then
                  current_state%sth%data(k,j,i) = t_rlx*&
                  (0.0_DEFAULT_PRECISION- current_state%th%data(k,j,i))
                end if
              end if
            end if

            if(current_state%number_q_fields .gt. 0)then
              if (current_state%immersed%indic_s(k,j,i) .eq. 1) then
                current_state%sq(iqv)%data(k,j,i) = 0.0_DEFAULT_PRECISION
                if (current_state%immersed%ghost_s(k,j,i) .lt. 1) then
                  current_state%sq(iqv)%data(k,j,i) = t_rlx*&
                  (0.0_DEFAULT_PRECISION- current_state%q(iqv)%data(k,j,i))
                end if
              end if
            end if
          end do
        end if ! ib_type=0
  
        if (current_state%immersed%ib_type==1)then

          do k=1,current_state%immersed%kmax_ji(j,i)
            if (current_state%immersed%ghost_u(k,j,i) .gt. 0) then ! u points
              ighost = current_state%immersed%ghost_u(k,j,i)
              current_state%u%data(k,j,i) = current_state%immersed%gp_store_u(ighost)
            end if
            if (current_state%immersed%ghost_v(k,j,i) .gt. 0) then ! v points
              ighost = current_state%immersed%ghost_v(k,j,i)
              current_state%v%data(k,j,i) = current_state%immersed%gp_store_v(ighost)
            end if
            if (current_state%immersed%ghost_w(k,j,i) .gt. 0) then ! w points
              ighost = current_state%immersed%ghost_w(k,j,i)
              current_state%w%data(k,j,i) = current_state%immersed%gp_store_w(ighost)
              if(current_state%use_viscosity_and_diffusion)then
                current_state%vis_coefficient%data(k,j,i) = current_state%immersed%gp_store_vis(ighost)
                current_state%diff_coefficient%data(k,j,i) = current_state%immersed%gp_store_diff(ighost)
              end if
            end if
            if(enable_theta)then
              if (current_state%immersed%ghost_s(k,j,i) .gt. 0) then ! s points
                ighost = current_state%immersed%ghost_s(k,j,i)
                current_state%th%data(k,j,i) = current_state%immersed%gp_store_th(ighost)
              end if
            end if
            if(current_state%number_q_fields .gt. 0)then
              if (current_state%immersed%ghost_s(k,j,i) .gt. 0) then ! s points
                ighost = current_state%immersed%ghost_s(k,j,i)
                current_state%q(iqv)%data(k,j,i) = current_state%immersed%gp_store_q(ighost)
              end if
            end if
!  
            if (current_state%immersed%indic_u(k,j,i) .eq. 1) then
              current_state%su%data(k,j,i) = t_rlx*&
              (0.0_DEFAULT_PRECISION - current_state%u%data(k,j,i))
            end if
            if (current_state%immersed%indic_v(k,j,i) .eq. 1) then
              current_state%sv%data(k,j,i) = t_rlx*&
              (0.0_DEFAULT_PRECISION - current_state%v%data(k,j,i))
            end if
            if (current_state%immersed%indic_w(k,j,i) .eq. 1) then
              current_state%sw%data(k,j,i) = t_rlx*&
              (0.0_DEFAULT_PRECISION - current_state%w%data(k,j,i))
              current_state%vis_coefficient%data(k,j,i) = 0.0_DEFAULT_PRECISION
              current_state%diff_coefficient%data(k,j,i) = 0.0_DEFAULT_PRECISION
            end if
  
            if(enable_theta)then
              if (current_state%immersed%indic_s(k,j,i) .eq. 1) then
                current_state%sth%data(k,j,i) = t_rlx*&
                (0.0_DEFAULT_PRECISION- current_state%th%data(k,j,i))
              end if
            end if

            if(current_state%number_q_fields .gt. 0)then
              if (current_state%immersed%indic_s(k,j,i) .eq. 1) then
                current_state%sq(iqv)%data(k,j,i) = t_rlx*&
                (0.0_DEFAULT_PRECISION- current_state%q(iqv)%data(k,j,i))
              end if
            end if

          end do ! k loop

        end if! ib_type=1

        if (current_state%use_viscosity_and_diffusion)then
          current_state%cvis=max(current_state%cvis,maxval(current_state%immersed%gp_store_vis),&
                                 maxval(current_state%immersed%gp_store_diff)*viscous_courant_coefficient)
        end if
        

      end if ! if ib_col
    end if ! if l_doset

  end subroutine timestep_callback

  !> Finalisation callback hook which will <DO WHAT?>
  !! @param current_state The current model state
  subroutine finalisation_callback(current_state)
    type(model_state_type), intent(inout), target :: current_state
  end subroutine finalisation_callback


end module ib_finalise_mod
