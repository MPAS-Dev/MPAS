module atm_time_integration

   use mpas_grid_types
   use mpas_kind_types
   use mpas_configure
   use mpas_constants
   use mpas_dmpar
   use mpas_vector_reconstruction
   ! Added only clause to keep xlf90 from getting confused from the overloaded abs intrinsic in mpas_timekeeping
   use mpas_timekeeping, only: MPAS_Time_type, MPAS_TimeInterval_type, &
                               mpas_set_time, mpas_set_timeInterval, mpas_get_time, operator(+), add_t_ti

#ifdef DO_PHYSICS
   use mpas_atmphys_driver_microphysics
   use mpas_atmphys_todynamics
   use mpas_atmphys_utilities
#endif

   contains


   subroutine atm_timestep(domain, dt, timeStamp, itimestep)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      character(len=*), intent(in) :: timeStamp
      integer, intent(in) :: itimestep

      type (block_type), pointer :: block
      type (MPAS_Time_type) :: currTime
      type (MPAS_TimeInterval_type) :: dtInterval
      character (len=StrKIND) :: xtime

      if (trim(config_time_integration) == 'SRK3') then
         call atm_srk3(domain, dt, itimestep)
      else
         write(0,*) 'Unknown time integration option '//trim(config_time_integration)
         write(0,*) 'Currently, only ''SRK3'' is supported.'
         stop
      end if

      call mpas_set_time(currTime, dateTimeString=timeStamp)
      call mpas_set_timeInterval(dtInterval, dt=dt)
      currTime = currTime + dtInterval
      call mpas_get_time(currTime, dateTimeString=xtime)

      block => domain % blocklist
      do while (associated(block))
         block % state % time_levs(2) % state % xtime % scalar = xtime
         block => block % next
      end do

   end subroutine atm_timestep


   subroutine atm_srk3(domain, dt, itimestep)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Advance model state forward in time by the specified time step using 
   !   time-split RK3 scheme
   !
   ! Hydrostatic (primitive eqns.) solver
   !
   ! Input: domain - current model state in time level 1 (e.g., time_levs(1)state%h(:,:)) 
   !                 plus grid meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:)) contains 
   !                  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (domain_type), intent(inout) :: domain
      real (kind=RKIND), intent(in) :: dt
      integer, intent(in) :: itimestep

      integer :: iCell, k, iEdge
      type (block_type), pointer :: block

      integer :: rk_step, number_of_sub_steps
      integer :: iScalar

      real (kind=RKIND), dimension(3) :: rk_timestep, rk_sub_timestep
      integer, dimension(3) :: number_sub_steps
      integer :: small_step
      logical, parameter :: debug = .false.
!      logical, parameter :: debug = .true.
      logical, parameter :: debug_mass_conservation = .true.

      real (kind=RKIND) :: scalar_min, scalar_max
      real (kind=RKIND) :: global_scalar_min, global_scalar_max


      !
      ! Initialize RK weights
      !

      number_of_sub_steps = config_number_of_sub_steps
      rk_timestep(1) = dt/3.
      rk_timestep(2) = dt/2.
      rk_timestep(3) = dt

      rk_sub_timestep(1) = dt/3.
      rk_sub_timestep(2) = dt/real(number_of_sub_steps)
      rk_sub_timestep(3) = dt/real(number_of_sub_steps)

      number_sub_steps(1) = 1
      number_sub_steps(2) = number_of_sub_steps/2
      number_sub_steps(3) = number_of_sub_steps

      if(debug) write(0,*) ' copy step in rk solver '

! WCS-parallel: it appears we have chosen to update all edges of nCellsSolve (to cut down on communications on acoustic steps).
! Do our communications patterns and loops (specifically in compute_dyn_tend) reflect this?  Or do they assume we are only updating 
! the so-called owned edges?



! WCS-parallel: first three and rtheta_p arise from scalar transport and microphysics update (OK).  Others come from where?

! theta_m
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(1) % state % theta_m)
 
! scalars
      call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(1) % state % scalars)

! pressure_p
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % pressure_p)

! rtheta_p
      call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rtheta_p)


      block => domain % blocklist
      do while (associated(block))
         ! We are setting values in the halo here, so no communications are needed.
         ! Alternatively, we could just set owned cells and edge values and communicate after this block loop.
         call atm_rk_integration_setup( block % state % time_levs(2) % state, block % state % time_levs(1) % state, block % diag )
         block => block % next
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! BEGIN RK loop 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      do rk_step = 1, 3  ! Runge-Kutta loop

         if(debug) write(0,*) ' rk substep ', rk_step

        block => domain % blocklist
        do while (associated(block))
           ! The coefficients are set for owned cells (cqw) and for all edges of owned cells, 
           ! thus no communications should be needed after this call.  
           ! We could consider combining this and the next block loop.
           call atm_compute_moist_coefficients( block % state % time_levs(2) % state, block % diag, block % mesh )
           block => block % next
        end do


        if (debug) write(0,*) ' compute_dyn_tend '
        block => domain % blocklist
        do while (associated(block))
           call atm_compute_dyn_tend( block % tend, block % state % time_levs(2) % state, block % diag, block % mesh, rk_step, dt )
           block => block % next
        end do
        if (debug) write(0,*) ' finished compute_dyn_tend '


#ifdef DO_PHYSICS
        if (debug) write(0,*) ' add physics tendencies '
        block => domain % blocklist
        do while (associated(block))
           call physics_addtend( block % mesh, &
                        block % state % time_levs(1) % state, &
                        block % diag, &
                        block % tend, &
                        block % tend_physics, &
                        block % state % time_levs(2) % state % rho_zz % array(:,:), &
                        block % diag % rho_edge % array(:,:), & 
                        rk_step )
           block => block % next
        end do
        if (debug) write(0,*) ' finished add physics tendencies '
#endif

!***********************************
!  we will need to communicate the momentum tendencies here - we want tendencies for all edges of owned cells
!  because we are solving for all edges of owned cells
!***********************************

! tend_u
         call mpas_dmpar_exch_halo_field(domain % blocklist % tend % u, (/ 1 /))

         block => domain % blocklist
            do while (associated(block))
!               call atm_set_smlstep_pert_variables( block % state % time_levs(1) % state, block % state % time_levs(2) % state,  &
               call atm_set_smlstep_pert_variables( block % tend, block % diag, block % mesh )
               call atm_compute_vert_imp_coefs( block % state % time_levs(2) % state, block % mesh, block % diag, rk_sub_timestep(rk_step) )
            block => block % next
         end do

         do small_step = 1, number_sub_steps(rk_step)

            if(debug) write(0,*) ' acoustic step ',small_step
      
            block => domain % blocklist
            do while (associated(block))
               call atm_advance_acoustic_step( block % state % time_levs(2) % state, block % diag, block % tend, &
                                           block % mesh, rk_sub_timestep(rk_step) )
               block => block % next
            end do

            if(debug) write(0,*) ' acoustic step complete '
  
            !  will need communications here for rtheta_pp

!  WCS-parallel: is this a candidate for a smaller stencil?  we need only communicate cells that share edges with owned cells.

! rtheta_pp
            call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rtheta_pp, (/ 1 /))
 
         end do  ! end of small stimestep loop

         !  will need communications here for rho_pp

! WCS-parallel: is communication of rw_p and rho_pp because of limiter (pd or mono scheme?),
!  or is it needed for the large-step variable recovery (to get decoupled variables)?
!  seems like only rho_pp needed...
!
!  or, do we need ru and u in the halo for diagnostics that are computed later in compute_solve_diagnostics?
!
!  rho_pp might be candidate for smaller stencil (same stencil as rtheta_pp above).

! MGD seems necessary
! rw_p
         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % diag % rw_p, (/ 1 /))
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rw_p)

! MGD seems necessary
! ru_p
         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % diag % ru_p, (/ 2 /))
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % ru_p)

! rho_pp
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rho_pp)

         ! the second layer of halo cells must be exchanged before calling atm_recover_large_step_variables
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rtheta_pp, (/ 2 /))


         block => domain % blocklist
         do while (associated(block))
            call atm_recover_large_step_variables( block % state % time_levs(2) % state,                 &
                                               block % diag, block % tend, block % mesh,                 &
                                               rk_timestep(rk_step), number_sub_steps(rk_step), rk_step  )
            block => block % next
         end do

!  ************  advection of moist variables here...

! u
         !CR: SMALLER STENCIL?: call mpas_dmpar_exch_halo_field(block % state % time_levs(2) % state % u, (/ 3 /))
         call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % u)


         if (config_scalar_advection) then

         block => domain % blocklist
         do while (associated(block))
            !
            ! Note: The advance_scalars_mono routine can be used without limiting, and thus, encompasses 
            !       the functionality of the advance_scalars routine; however, it is noticeably slower, 
            !       so we keep the advance_scalars routine as well
            !
            if (rk_step < 3 .or. (.not. config_monotonic .and. .not. config_positive_definite)) then
               call atm_advance_scalars( block % tend, &
                                     block % state % time_levs(1) % state, block % state % time_levs(2) % state, &
                                     block % diag, &
                                     block % mesh, rk_timestep(rk_step) )
            else
               block % domain = domain 
               call atm_advance_scalars_mono( block % tend, &
                                              block % state % time_levs(1) % state, block % state % time_levs(2) % state, &
                                              block % diag, block % mesh, &
                                              rk_timestep(rk_step))
            end if
            block => block % next
         end do

! For now, we do scalar halo updates later on...
!         block => domain % blocklist
!         do while (associated(block))
!            call mpas_dmpar_exch_halo_field3d_real(domain % dminfo, block % tend % scalars % array(:,:,:), &
!                                             block % tend % num_scalars, block % mesh % nVertLevels, block % mesh % nCells, &
!                                             block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
!            call mpas_dmpar_exch_halo_field3d_real(domain % dminfo, block % state % time_levs(2) % state % scalars % array(:,:,:), &
!                                             block % state % time_levs(2) % state % num_scalars, block % mesh % nVertLevels, block % mesh % nCells, &
!                                             block % parinfo % cellsToSend, block % parinfo % cellsToRecv)
!            block => block % next
!         end do

         else
 
            write(0,*) ' no scalar advection '

         end if

! WCS-parallel: seems like we already use u and w (and other state variables) as if they were already correctly set in the halo,
! but they are not (at least to the outer edges - see communications below? or are those communications redundent?).  
! Perhaps we should communicate u, w, theta_m, rho_zz, etc after recover_large_step_variables),
! cover it with the scalar communications, and then compute solve_diagnostics.  I do not think we need to communicate the stuff we compute
! in compute_solve_diagnostics if we compute it out in the halo (and I think we do - the halos should be large enough).

         block => domain % blocklist
         do while (associated(block))
            call atm_compute_solve_diagnostics( dt, block % state % time_levs(2) % state, block % diag, block % mesh )
            block => block % next
         end do

         if(debug) write(0,*) ' diagnostics complete '


      ! need communications here to fill out u, w, theta_m, p, and pp, scalars, etc  
      ! so that they are available for next RK step or the first rk substep of the next timestep

!MGD seems necessary
! w
         call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % w)

! pv_edge
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % pv_edge)

! rho_edge
         call mpas_dmpar_exch_halo_field(domain % blocklist % diag % rho_edge)

!  ****  this will always be needed - perhaps we can cover this with compute_solve_diagnostics

! scalars
         if(rk_step < 3) then
            call mpas_dmpar_exch_halo_field(domain % blocklist % state % time_levs(2) % state % scalars)
         end if

       end do ! rk_step loop

!...  compute full velocity vectors at cell centers:
      block => domain % blocklist
        do while (associated(block))
           call mpas_reconstruct(block % mesh, block % state % time_levs(2) % state % u % array, &
                                 block % diag % uReconstructX % array,                           &
                                 block % diag % uReconstructY % array,                           &
                                 block % diag % uReconstructZ % array,                           &
                                 block % diag % uReconstructZonal % array,                       &
                                 block % diag % uReconstructMeridional % array                   &
                                )
           block => block % next
        end do

!... call to parameterizations of cloud microphysics. calculation of the tendency of water vapor to horizontal and
!... vertical advection needed for the Tiedtke parameterization of convection.
#ifdef DO_PHYSICS
      block => domain % blocklist
      do while(associated(block))

         !NOTE: The calculation of the tendency due to horizontal and vertical advection for the water vapor mixing ratio
         !requires that the subroutine atm_advance_scalars_mono was called on the third Runge Kutta step, so that a halo
         !update for the scalars at time_levs(1) is applied. A halo update for the scalars at time_levs(2) is done above. 
         if(config_monotonic) then
            block % tend_physics % rqvdynten % array(:,:) = &
                 ( block % state % time_levs(2) % state % scalars % array(block % state % time_levs(2) % state % index_qv,:,:)   &
                 - block % state % time_levs(1) % state % scalars % array(block % state % time_levs(1) % state % index_qv,:,:) ) &
                 / config_dt
         else
            block % tend_physics % rqvdynten % array(:,:) = 0._RKIND
         endif

         !simply set to zero negative mixing ratios of different water species (for now):
         where ( block % state % time_levs(2) % state % scalars % array(:,:,:) .lt. 0.) &
            block % state % time_levs(2) % state % scalars % array(:,:,:) = 0.

         !call microphysics schemes:
         if(config_microp_scheme .ne. 'off') &
            call microphysics_driver ( block % state % time_levs(2) % state, block % diag, block % diag_physics, &
                                       block % tend, block % mesh, itimestep )

         block => block % next
      end do
#endif


!      if(debug) then
!        101 format(' local  min, max scalar',i4,2(1x,e17.10))
        102 format(' global min, max scalar',i4,2(1x,e17.10))
        write(0,*)
        block => domain % blocklist
          do while (associated(block))
             scalar_min = 0.
             scalar_max = 0.
             do iCell = 1, block % mesh % nCellsSolve
             do k = 1, block % mesh % nVertLevels
               scalar_min = min(scalar_min, block % state % time_levs(2) % state % w % array(k,iCell))
               scalar_max = max(scalar_max, block % state % time_levs(2) % state % w % array(k,iCell))
             enddo
             enddo
             call mpas_dmpar_min_real(domain%dminfo, scalar_min, global_scalar_min)
             call mpas_dmpar_max_real(domain%dminfo, scalar_max, global_scalar_max)
!             write(0,*) 'local  min, max w ',scalar_min, scalar_max
             write(0,*) 'global min, max w ',global_scalar_min, global_scalar_max

             scalar_min = 0.
             scalar_max = 0.
             do iEdge = 1, block % mesh % nEdgesSolve
             do k = 1, block % mesh % nVertLevels
               scalar_min = min(scalar_min, block % state % time_levs(2) % state % u % array(k,iEdge))
               scalar_max = max(scalar_max, block % state % time_levs(2) % state % u % array(k,iEdge))
             enddo
             enddo
             call mpas_dmpar_min_real(domain%dminfo, scalar_min, global_scalar_min)
             call mpas_dmpar_max_real(domain%dminfo, scalar_max, global_scalar_max)
!             write(0,*) 'local  min, max u ',scalar_min, scalar_max
             write(0,*) 'global min, max u ',global_scalar_min, global_scalar_max

             do iScalar = 1, block % state % time_levs(2) % state % num_scalars
                scalar_min = 0.
                scalar_max = 0.
                do iCell = 1, block % mesh % nCellsSolve
                do k = 1, block % mesh % nVertLevels
                  scalar_min = min(scalar_min, block % state % time_levs(2) % state % scalars % array(iScalar,k,iCell))
                  scalar_max = max(scalar_max, block % state % time_levs(2) % state % scalars % array(iScalar,k,iCell))
                enddo
                enddo
                call mpas_dmpar_min_real(domain%dminfo, scalar_min, global_scalar_min)
                call mpas_dmpar_max_real(domain%dminfo, scalar_max, global_scalar_max)
!                write(0,101) iScalar,scalar_min,scalar_max
                write(0,102) iScalar,global_scalar_min,global_scalar_max
             end do
             block => block % next

          end do
!      end if


   end subroutine atm_srk3

!---

   subroutine atm_rk_integration_setup( s_old, s_new, diag )

     implicit none
     type (state_type) :: s_new, s_old
     type (diag_type) :: diag

     diag % ru_save % array = diag % ru % array
     diag % rw_save % array = diag % rw % array
     diag % rtheta_p_save % array = diag % rtheta_p % array
     diag % rho_p_save % array = diag % rho_p % array

     s_old % u % array = s_new % u % array
     s_old % w % array = s_new % w % array
     s_old % theta_m % array = s_new % theta_m % array
     s_old % rho_zz % array = s_new % rho_zz % array
     s_old % scalars % array = s_new % scalars % array

   end subroutine atm_rk_integration_setup

!-----

   subroutine atm_compute_moist_coefficients( state, diag, grid )

      implicit none
      type (state_type) :: state
      type (diag_type) :: diag
      type (mesh_type) :: grid

      integer :: iEdge, iCell, k, cell1, cell2, iq
      integer :: nCells, nEdges, nVertLevels, nCellsSolve
      real (kind=RKIND) :: qtot

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels
      nCellsSolve = grid % nCellsSolve

        do iCell = 1, nCellsSolve
          do k = 2, nVertLevels
            qtot = 0.
            do iq = state % moist_start, state % moist_end
              qtot = qtot + 0.5 * (state % scalars % array (iq, k, iCell) + state % scalars % array (iq, k-1, iCell))
            end do
            diag % cqw % array(k,iCell) = 1./(1.+qtot)
          end do
        end do

        do iEdge = 1, nEdges
          cell1 = grid % cellsOnEdge % array(1,iEdge)
          cell2 = grid % cellsOnEdge % array(2,iEdge)
          if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
            do k = 1, nVertLevels
              qtot = 0.
              do iq = state % moist_start, state % moist_end
                 qtot = qtot + 0.5 * ( state % scalars % array (iq, k, cell1) + state % scalars % array (iq, k, cell2) )
              end do
              diag % cqu % array(k,iEdge) = 1./( 1. + qtot)
            end do
          end if
        end do

   end subroutine atm_compute_moist_coefficients

!---

   subroutine atm_compute_vert_imp_coefs(s, grid, diag, dts)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Compute coefficients for vertically implicit gravity-wave/acoustic computations
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: diag - cofrz, cofwr, cofwz, coftz, cofwt, a, alpha and gamma
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      type (state_type), intent(in) :: s
      type (mesh_type), intent(in) :: grid
      type (diag_type), intent(inout) :: diag
      real (kind=RKIND), intent(in) :: dts

      integer :: i, k, iq

      integer :: nCells, nVertLevels, nCellsSolve
      real (kind=RKIND), dimension(:,:), pointer :: zz, cqw, p, t, rb, rtb, pb, rt
      real (kind=RKIND), dimension(:,:), pointer :: cofwr, cofwz, coftz, cofwt, a_tri, alpha_tri, gamma_tri
      real (kind=RKIND), dimension(:), pointer :: cofrz, rdzw, fzm, fzp, rdzu

      real (kind=RKIND), dimension( grid % nVertLevels ) :: b_tri,c_tri
      real (kind=RKIND) :: epssm, dtseps, c2, qtot, rcv

!  set coefficients

      nCells      = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels
!      epssm = grid % epssm  !  this should come in through the namelist  ******************
!      epssm = 0.1
      epssm = config_epssm

      rdzu => grid % rdzu % array
      rdzw => grid % rdzw % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      zz => grid % zz % array
      cqw => diag % cqw % array

      p => diag % exner % array
      pb => diag % exner_base % array
      rt => diag % rtheta_p % array
      rtb => diag % rtheta_base % array
      rb => diag % rho_base % array

      alpha_tri => diag % alpha_tri % array
      gamma_tri => diag % gamma_tri % array
      a_tri => diag % a_tri % array
      cofwr => diag % cofwr % array      
      cofwz => diag % cofwz % array      
      coftz => diag % coftz % array      
      cofwt => diag % cofwt % array      
      cofrz => diag % cofrz % array      

      t => s % theta_m % array

      dtseps = .5*dts*(1.+epssm)
      rcv = rgas/(cp-rgas)
      c2 = cp*rcv

      do k=1,nVertLevels
         cofrz(k) = dtseps*rdzw(k)
      end do

      do i = 1, nCellsSolve  !  we only need to do cells we are solving for, not halo cells

        do k=2,nVertLevels
          cofwr(k,i) =.5*dtseps*gravity*(fzm(k)*zz(k,i)+fzp(k)*zz(k-1,i))
        end do
        coftz(1,i) = 0.0
        do k=2,nVertLevels
           cofwz(k,i) = dtseps*c2*(fzm(k)*zz(k,i)+fzp(k)*zz(k-1,i))  &
                *rdzu(k)*cqw(k,i)*(fzm(k)*p (k,i)+fzp(k)*p (k-1,i))
           coftz(k,i) = dtseps*   (fzm(k)*t (k,i)+fzp(k)*t (k-1,i))
        end do
        coftz(nVertLevels+1,i) = 0.0
        do k=1,nVertLevels

          qtot = 0.
          do iq = s % moist_start, s % moist_end
            qtot = qtot + s % scalars % array (iq, k, i)
          end do

          cofwt(k,i) = .5*dtseps*rcv*zz(k,i)*gravity*rb(k,i)/(1.+qtot)  &
                              *p(k,i)/((rtb(k,i)+rt(k,i))*pb(k,i))
        end do

        a_tri(1,i) = 0.  ! note, this value is never used
        b_tri(1) = 1.    ! note, this value is never used
        c_tri(1) = 0.    ! note, this value is never used
        gamma_tri(1,i) = 0.
        alpha_tri(1,i) = 0.  ! note, this value is never used

        do k=2,nVertLevels
          a_tri(k,i) = -cofwz(k  ,i)* coftz(k-1,i)*rdzw(k-1)*zz(k-1,i)   &
                       +cofwr(k  ,i)* cofrz(k-1  )                       &
                       -cofwt(k-1,i)* coftz(k-1,i)*rdzw(k-1)
          b_tri(k) = 1.                                                  &
                       +cofwz(k  ,i)*(coftz(k  ,i)*rdzw(k  )*zz(k  ,i)   &
                                    +coftz(k  ,i)*rdzw(k-1)*zz(k-1,i))   &
                       -coftz(k  ,i)*(cofwt(k  ,i)*rdzw(k  )             &
                                     -cofwt(k-1,i)*rdzw(k-1))            &
                       +cofwr(k,  i)*(cofrz(k    )-cofrz(k-1))
          c_tri(k) =   -cofwz(k  ,i)* coftz(k+1,i)*rdzw(k  )*zz(k  ,i)   &
                       -cofwr(k  ,i)* cofrz(k    )                       &
                       +cofwt(k  ,i)* coftz(k+1,i)*rdzw(k  )
        end do
        do k=2,nVertLevels
          alpha_tri(k,i) = 1./(b_tri(k)-a_tri(k,i)*gamma_tri(k-1,i))
          gamma_tri(k,i) = c_tri(k)*alpha_tri(k,i)
        end do

      end do ! loop over cells

   end subroutine atm_compute_vert_imp_coefs

!------------------------

   subroutine atm_set_smlstep_pert_variables( tend, diag, grid )

      implicit none
      type (tend_type) :: tend
      type (diag_type) :: diag
      type (mesh_type) :: grid
      !SHP-w
      integer :: iCell, iEdge, k, cell1, cell2, coef_3rd_order
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, dvEdge, areaCell
      real (kind=RKIND) :: flux
      !SHP-w
      coef_3rd_order = config_coef_3rd_order
      if(config_theta_adv_order /=3) coef_3rd_order = 0

      !SHP-w
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      dvEdge => grid % dvEdge % array
      areaCell => grid % areaCell % array
      cellsOnEdge => grid % cellsOnEdge % array

      diag % rho_pp % array = diag % rho_p_save % array - diag % rho_p % array

      diag % ru_p % array = diag % ru_save % array - diag % ru % array
      diag % rtheta_pp % array = diag % rtheta_p_save % array - diag % rtheta_p % array
      diag % rtheta_pp_old % array = diag % rtheta_pp % array
      diag % rw_p % array = diag % rw_save % array - diag % rw % array

      do iCell = 1, grid % nCellsSolve
      do k = 2, grid % nVertLevels
        tend % w % array(k,iCell) = ( fzm(k) * grid % zz % array(k  ,iCell) +   &
                                      fzp(k) * grid % zz % array(k-1,iCell)   ) &
                                     * tend % w % array(k,iCell)
      end do
      end do

      do iEdge = 1,grid % nEdges

         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)

         !SHP-w
         do k = 2, grid%nVertLevels
            flux = fzm(k) * tend % u % array(k,iEdge) + fzp(k) * tend % u % array(k-1,iEdge)
            tend % w % array(k,cell2) = tend % w % array(k,cell2)   &
                     + (grid % zb % array(k,2,iEdge) + coef_3rd_order*sign(1.0_RKIND,tend % u % array(k,iEdge))*grid %zb3 % array(k,2,iEdge))*flux   &
                     * (fzm(k) * grid % zz % array(k,cell2) + fzp(k) * grid % zz % array(k-1,cell2)) 
            tend % w % array(k,cell1) = tend % w % array(k,cell1)   &
                     - (grid % zb % array(k,1,iEdge) + coef_3rd_order*sign(1.0_RKIND,tend % u % array(k,iEdge))*grid %zb3 % array(k,1,iEdge))*flux   &
                     * (fzm(k) * grid % zz % array(k,cell1) + fzp(k) * grid % zz % array(k-1,cell1)) 
         end do

      end do

      diag % ruAvg % array = 0.
      diag % wwAvg % array = 0.

   end subroutine atm_set_smlstep_pert_variables

!-------------------------------

   subroutine atm_advance_acoustic_step( s, diag, tend, grid, dts )

      implicit none

      type (state_type) :: s
      type (diag_type) :: diag
      type (tend_type) :: tend
      type (mesh_type) :: grid
      real (kind=RKIND), intent(in) :: dts

      real (kind=RKIND), dimension(:,:), pointer :: rho_zz, theta_m, ru_p, rw_p, rtheta_pp,    &
                                                    rtheta_pp_old, zz, exner, cqu, ruAvg, &
                                                    wwAvg, rho_pp, cofwt, coftz, zx,      &
                                                    a_tri, alpha_tri, gamma_tri, dss,     &
                                                    tend_ru, tend_rho, tend_rt, tend_rw,  &
                                                    zgrid, cofwr, cofwz, w, h_divergence
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, rdzw, dcEdge, AreaCell, cofrz, dvEdge

      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, pzp, pzm

      real (kind=RKIND) :: smdiv, c2, rcv
      real (kind=RKIND), dimension( grid % nVertLevels ) :: du
      real (kind=RKIND), dimension( grid % nVertLevels + 1 ) :: dpzx
      real (kind=RKIND), dimension( grid % nVertLevels, grid % nCells+1 ) :: ts, rs

      integer :: cell1, cell2, iEdge, iCell, k
      real (kind=RKIND) :: pgrad, flux, resm, epssm

      real (kind=RKIND) :: cf1, cf2, cf3, pr, pl
      integer :: kr, kl

      integer :: nEdges, nCells, nCellsSolve, nVertLevels

      logical, parameter :: debug = .false.
!      logical, parameter :: debug = .true.
      logical, parameter :: debug1 = .false.
      logical :: newpx

!--

      rho_zz => s % rho_zz % array
      theta_m => s % theta_m % array
      w => s % w % array

      rtheta_pp => diag % rtheta_pp % array
      rtheta_pp_old => diag % rtheta_pp_old % array
      h_divergence => diag % h_divergence % array
      ru_p => diag % ru_p % array
      rw_p => diag % rw_p % array
      exner => diag % exner % array
      cqu => diag % cqu % array
      ruAvg => diag % ruAvg % array
      wwAvg => diag % wwAvg % array
      rho_pp => diag % rho_pp % array
      cofwt => diag % cofwt % array
      coftz => diag % coftz % array
      cofrz => diag % cofrz % array
      cofwr => diag % cofwr % array
      cofwz => diag % cofwz % array
      a_tri => diag % a_tri % array
      alpha_tri => diag % alpha_tri % array
      gamma_tri => diag % gamma_tri % array
      dss => grid % dss % array

      pzp  => grid % pzp % array
      pzm  => grid % pzm % array

      tend_ru => tend % u % array
      tend_rho => tend % rho_zz % array
      tend_rt => tend % theta_m % array
      tend_rw => tend % w % array

      zz => grid % zz % array
      zx => grid % zx % array
      zgrid => grid % zgrid % array
      fzm => grid % fzm % array
      fzp => grid % fzp % array
      rdzw => grid % rdzw % array
      dcEdge => grid % dcEdge % array
      dvEdge => grid % dvEdge % array
      AreaCell => grid % AreaCell % array

!  might these be pointers instead? **************************

      nEdges = grid % nEdges
      nCells = grid % nCells
      nCellsSolve = grid % nCellsSolve
      nVertLevels = grid % nVertLevels

      cf1 = grid % cf1 % scalar
      cf2 = grid % cf2 % scalar
      cf3 = grid % cf3 % scalar

      cpr         => grid % cpr % array
      cpl         => grid % cpl % array
      newpx = config_newpx

      epssm = config_epssm
      smdiv = config_smdiv

      rcv = rgas/(cp-rgas)
      c2 = cp*rcv
      resm   = (1.-epssm)/(1.+epssm)

      ts = 0.
      rs = 0.

      ! acoustic step divergence damping - forward weight rtheta_pp
      rtheta_pp_old = rtheta_pp + smdiv*(rtheta_pp - rtheta_pp_old)

      if(debug) write(0,*) ' updating ru_p '

      do iEdge = 1, nEdges
 
         cell1 = grid % cellsOnEdge % array (1,iEdge)
         cell2 = grid % cellsOnEdge % array (2,iEdge)
         ! update edge for block-owned cells
         if (cell1 <= grid % nCellsSolve .or. cell2 <= grid % nCellsSolve ) then

            if (newpx) then

               k = 1
               pr  =   cpr(k  ,iEdge)*zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)   &
                     + cpr(k+1,iEdge)*zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)   &
                     + cpr(k+2,iEdge)*zz(k+2,cell2)*rtheta_pp_old(k+2,cell2)

               pl  =   cpl(k  ,iEdge)*zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)   &
                     + cpl(k+1,iEdge)*zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)   &
                     + cpl(k+2,iEdge)*zz(k+2,cell1)*rtheta_pp_old(k+2,cell1)
               pgrad = 2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)
               pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
               du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 

               do k=2,nVertLevels

                  kr = min(nVertLevels,k+ nint(.5-sign(0.5_RKIND,zx(k,iEdge)+zx(k+1,iEdge))))
                  kl = min(nVertLevels,2*k+1-kr)
                  pr = zz(k,cell2)*rtheta_pp_old(k ,cell2)+.5*(zgrid(k   ,cell1) +zgrid(k +1,cell1)   &
                                                              -zgrid(k   ,cell2) -zgrid(k +1,cell2))  &
                                                             /(zgrid(kr+1,cell2) -zgrid(kr-1,cell2))  &
                    *(zz(kr,cell2)*rtheta_pp_old(kr,cell2)-zz(kr-1,cell2)*rtheta_pp_old(kr-1,cell2))
                  pl = zz(k,cell1)*rtheta_pp_old(k ,cell1)+.5*(zgrid(k   ,cell2) +zgrid(k +1,cell2)   &
                                                              -zgrid(k   ,cell1) -zgrid(k +1,cell1))  &
                                                             /(zgrid(kl+1,cell1) -zgrid(kl-1,cell1))  &
                    *(zz(kl,cell1)*rtheta_pp_old(kl,cell1)-zz(kl-1,cell1)*rtheta_pp_old(kl-1,cell1))
                  pgrad = 2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)
                  pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
                  du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 
               end do

            else

               k = 1
!               dpzx(k) = .5*zx(k,iEdge)*(cf1*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)    &
!                                             +zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))   &
!                                        +cf2*(zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)    &
!                                             +zz(k+1,cell1)*rtheta_pp_old(k+1,cell1))   &
!                                        +cf3*(zz(k+2,cell2)*rtheta_pp_old(k+2,cell2)    &
!                                             +zz(k+2,cell1)*rtheta_pp_old(k+2,cell1)))

               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                               &
                         *(pzm(k,cell2)*(zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)        &
                                        -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))       &
                          +pzm(k,cell1)*(zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)        &
                                        -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))       &
                          +pzp(k,cell2)*(zz(k+2,cell2)*rtheta_pp_old(k+2,cell2)        &
                                        -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))       &
                          +pzp(k,cell1)*(zz(k+2,cell1)*rtheta_pp_old(k+2,cell1)        &
                                        -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)))

               do k=2,grid % nVertLevels-1
!                  dpzx(k)=.5*zx(k,iEdge)*(fzm(k)*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)   &
!                                                 +zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))  &
!                                         +fzp(k)*(zz(k-1,cell2)*rtheta_pp_old(k-1,cell2)   &
!                                                 +zz(k-1,cell1)*rtheta_pp_old(k-1,cell1)))
                  dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                                   &
                                   *(pzp(k,cell2)*(zz(k+1,cell2)*rtheta_pp_old(k+1,cell2)     &
                                                  -zz(k  ,cell2)*rtheta_pp_old(k  ,cell2))    &
                                    +pzm(k,cell2)*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)     &
                                                  -zz(k-1,cell2)*rtheta_pp_old(k-1,cell2))    &
                                    +pzp(k,cell1)*(zz(k+1,cell1)*rtheta_pp_old(k+1,cell1)     &
                                                  -zz(k  ,cell1)*rtheta_pp_old(k  ,cell1))    &
                                    +pzm(k,cell1)*(zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)     &
                                                  -zz(k-1,cell1)*rtheta_pp_old(k-1,cell1)))
               end do

               k=grid % nVertLevels
               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                                   &
                                *(pzm(k,cell2)*(zz(k  ,cell2)*rtheta_pp_old(k  ,cell2)     &
                                               -zz(k-1,cell2)*rtheta_pp_old(k-1,cell2))    &
                                 +pzm(k,cell1)*(zz(k  ,cell1)*rtheta_pp_old(k  ,cell1)     &
                                               -zz(k-1,cell1)*rtheta_pp_old(k-1,cell1)))

!               dpzx(nVertLevels + 1) = 0.

               do k=1,nVertLevels
!                  pgrad =  (rtheta_pp_old(k,cell2)-rtheta_pp_old(k,cell1))/dcEdge(iEdge)  &
!                               - rdzw(k)*(dpzx(k+1)-dpzx(k))
                  pgrad =     ((rtheta_pp_old(k,cell2)*zz(k,cell2)                    &
                               -rtheta_pp_old(k,cell1)*zz(k,cell1))/dcEdge(iEdge)     &
                            -dpzx(k))/(.5*(zz(k,cell2)+zz(k,cell1)))
                  pgrad = 0.5*c2*(exner(k,cell1)+exner(k,cell2))*pgrad
                  du(k) = dts*(tend_ru(k,iEdge) - cqu(k,iEdge) * pgrad) 
!                          + (0.05/6.)*dcEdge(iEdge)*(h_divergence(k,cell2)-h_divergence(k,cell1))
               end do
            end if

            do k=1,nVertLevels
               ru_p(k,iEdge) = ru_p(k,iEdge) + du(k)

               if(debug) then
                 if(iEdge == 3750) then
                   write(0,*) ' k, pgrad, tend_ru ',k,pgrad,tend_ru(k,3750)
                 end if
               end if

!  need to add horizontal fluxes into density update, rtheta update and w update

               flux = dts*dvEdge(iEdge)*ru_p(k,iEdge)
               rs(k,cell1) = rs(k,cell1)-flux/AreaCell(cell1)
               rs(k,cell2) = rs(k,cell2)+flux/AreaCell(cell2)
   
               flux = flux*0.5*(theta_m(k,cell2)+theta_m(k,cell1))
               ts(k,cell1) = ts(k,cell1)-flux/AreaCell(cell1)
               ts(k,cell2) = ts(k,cell2)+flux/AreaCell(cell2)
   
               ruAvg(k,iEdge) = ruAvg(k,iEdge) + ru_p(k,iEdge)

            end do

        end if ! end test for block-owned cells

      end do ! end loop over edges

      ! saving rtheta_pp before update for use in divergence damping in next acoustic step
      rtheta_pp_old(:,:) = rtheta_pp(:,:)

      do iCell = 1, nCellsSolve

        do k=1, nVertLevels
          rs(k,iCell) = rho_pp(k,iCell) + dts*tend_rho(k,iCell) + rs(k,iCell)      &
                          - cofrz(k)*resm*(rw_p(k+1,iCell)-rw_p(k,iCell))
          ts(k,iCell) = rtheta_pp(k,iCell) + dts*tend_rt(k,iCell) + ts(k,iCell)    &
                             - resm*rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)      &
                             -coftz(k,iCell)*rw_p(k,iCell))
        enddo

        do k=2, nVertLevels

          wwavg(k,iCell) = wwavg(k,iCell) + 0.5*(1.-epssm)*rw_p(k,iCell)

          rw_p(k,iCell) = rw_p(k,iCell) +  dts*tend_rw(k,iCell)                       &
                     - cofwz(k,iCell)*((zz(k  ,iCell)*ts (k  ,iCell)                  &
                                   -zz(k-1,iCell)*ts (k-1,iCell))                     &
                             +resm*(zz(k  ,iCell)*rtheta_pp(k  ,iCell)                &
                                   -zz(k-1,iCell)*rtheta_pp(k-1,iCell)))              &
                     - cofwr(k,iCell)*((rs (k,iCell)+rs (k-1,iCell))                  &
                             +resm*(rho_pp(k,iCell)+rho_pp(k-1,iCell)))               &
                     + cofwt(k  ,iCell)*(ts (k  ,iCell)+resm*rtheta_pp(k  ,iCell))    &
                     + cofwt(k-1,iCell)*(ts (k-1,iCell)+resm*rtheta_pp(k-1,iCell))
        enddo

        do k=2,nVertLevels
          rw_p(k,iCell) = (rw_p(k,iCell)-a_tri(k,iCell)*rw_p(k-1,iCell))*alpha_tri(k,iCell)
        end do

        do k=nVertLevels,1,-1
          rw_p(k,iCell) = rw_p(k,iCell) - gamma_tri(k,iCell)*rw_p(k+1,iCell)     
        end do

        do k=2,nVertLevels
           rw_p(k,iCell) = (rw_p(k,iCell)-dts*dss(k,iCell)*               &
                       (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))        &
                       *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell))       &
                                *w(k,iCell)    )/(1.+dts*dss(k,iCell))

           wwAvg(k,iCell) = wwAvg(k,iCell) + 0.5*(1.+epssm)*rw_p(k,iCell)

        end do

        do k=1,nVertLevels
          rho_pp(k,iCell) = rs(k,iCell) - cofrz(k) *(rw_p(k+1,iCell)-rw_p(k  ,iCell))
          rtheta_pp(k,iCell) = ts(k,iCell) - rdzw(k)*(coftz(k+1,iCell)*rw_p(k+1,iCell)  &
                             -coftz(k  ,iCell)*rw_p(k  ,iCell))
        end do

      end do !  end of loop over cells

   end subroutine atm_advance_acoustic_step

!------------------------

   subroutine atm_recover_large_step_variables( s, diag, tend, grid, dt, ns, rk_step )

      implicit none
      type (state_type) :: s
      type (diag_type) :: diag
      type (tend_type) :: tend
      type (mesh_type) :: grid
      integer, intent(in) :: ns, rk_step
      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(:,:), pointer :: wwAvg, rw_save, w, rw, rw_p, rtheta_p, rtheta_pp,   &
                                                    rtheta_p_save, rt_diabatic_tend, rho_p, rho_p_save, &
                                                    rho_pp, rho_zz, rho_base, ruAvg, ru_save, ru_p, u, ru, &
                                                    exner, exner_base, rtheta_base, pressure_p,         &
                                                    zz, theta_m, pressure_b, qvapor
      real (kind=RKIND), dimension(:), pointer :: fzm, fzp, dvEdge, AreaCell
      real (kind=RKIND), dimension(:,:,:), pointer :: zb, zb3 
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer :: iCell, iEdge, k, cell1, cell2
      integer :: nVertLevels, nCells, nCellsSolve, nEdges, nEdgesSolve
      real (kind=RKIND) :: rcv, p0, cf1, cf2, cf3, flux, coef_3rd_order

!      logical, parameter :: debug=.true.
      logical, parameter :: debug=.false.

!---
       wwAvg => diag % wwAvg % array
       rw_save => diag % rw_save % array
       rw => diag % rw % array
       rw_p => diag % rw_p % array
       w => s % w % array

       rtheta_p => diag % rtheta_p % array
       rtheta_p_save => diag % rtheta_p_save % array
       rtheta_pp => diag % rtheta_pp % array
       rtheta_base => diag % rtheta_base % array
       rt_diabatic_tend => tend % rt_diabatic_tend % array
       theta_m => s % theta_m % array
       qvapor => s % scalars % array(s%index_qv,:,:)

       rho_zz => s % rho_zz % array
       rho_p => diag % rho_p % array
       rho_p_save => diag % rho_p_save % array
       rho_pp => diag % rho_pp % array
       rho_base => diag % rho_base % array

       ruAvg => diag % ruAvg % array
       ru_save => diag % ru_save % array
       ru_p => diag % ru_p % array
       ru => diag % ru % array
       u => s % u % array

       exner => diag % exner % array
       exner_base => diag % exner_base % array

       pressure_p => diag % pressure_p % array
       pressure_b => diag % pressure_base % array

       zz => grid % zz % array
       zb => grid % zb % array
       zb3 => grid % zb3 % array
       fzm => grid % fzm % array
       fzp => grid % fzp % array
       dvEdge => grid % dvEdge % array
       AreaCell => grid % AreaCell % array
       CellsOnEdge => grid % CellsOnEdge % array

       nVertLevels = grid % nVertLevels
       nCells = grid % nCells
       nCellsSolve = grid % nCellsSolve
       nEdges = grid % nEdges
       nEdgesSolve = grid % nEdgesSolve

       rcv = rgas/(cp-rgas)
       p0 = 1.e+05  ! this should come from somewhere else...

       cf1 = grid % cf1 % scalar
       cf2 = grid % cf2 % scalar
       cf3 = grid % cf3 % scalar
       coef_3rd_order = config_coef_3rd_order
       if(config_theta_adv_order /=3) coef_3rd_order = 0

      ! compute new density everywhere so we can compute u from ru.
      ! we will also need it to compute theta_m below

      do iCell = 1, nCells

        do k = 1, nVertLevels

          rho_p(k,iCell) = rho_p(k,iCell) + rho_pp(k,iCell)

          rho_zz(k,iCell) = rho_p(k,iCell) + rho_base(k,iCell)
        end do

      !  recover owned-cell values in block

!        if( iCell <= nCellsSolve ) then    ! If using this test, more halo exchanges will be needed
!  WCS-parallel: OK


          w(1,iCell) = 0.
          do k = 2, nVertLevels
            wwAvg(k,iCell) = rw(k,iCell) + (wwAvg(k,iCell) / float(ns))

            rw(k,iCell) = rw(k,iCell) + rw_p(k,iCell)


          ! pick up part of diagnosed w from omega
            w(k,iCell) = rw(k,iCell)/( (fzm(k)*zz (k,iCell)+fzp(k)*zz (k-1,iCell))   &
                                      *(fzm(k)*rho_zz(k,iCell)+fzp(k)*rho_zz(k-1,iCell)) )
          end do
          w(nVertLevels+1,iCell) = 0.

          if(rk_step == 3) then
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p(k,iCell) + rtheta_pp(k,iCell) &
                                 - dt * rho_zz(k,iCell) * rt_diabatic_tend(k,iCell)
            end do
          else
            do k = 1, nVertLevels
               rtheta_p(k,iCell) = rtheta_p(k,iCell) + rtheta_pp(k,iCell)
            end do
          end if

          do k = 1, nVertLevels
            theta_m(k,iCell) = (rtheta_p(k,iCell) + rtheta_base(k,iCell))/rho_zz(k,iCell)
            exner(k,iCell) = (zz(k,iCell)*(rgas/p0)*(rtheta_p(k,iCell)+rtheta_base(k,iCell)))**rcv
             ! pressure below is perturbation pressure
            pressure_p(k,iCell) = zz(k,iCell) * rgas * (exner(k,iCell)*rtheta_p(k,iCell)+rtheta_base(k,iCell)  &
                                                          * (exner(k,iCell)-exner_base(k,iCell)))
          end do

      end do

      ! recover time-averaged ruAvg on all edges of owned cells (for upcoming scalar transport).  
      ! we solved for these in the acoustic-step loop.  
      ! we will compute ru and u here also, given we are here, even though we only need them on nEdgesSolve

      do iEdge = 1, nEdges

        cell1 = CellsOnEdge(1,iEdge)
        cell2 = CellsOnEdge(2,iEdge)

!        if( cell1 <= nCellsSolve .or. cell2 <= nCellsSolve ) then   ! If using this test, more halo exchanges will be needed,
                                                                     ! though we can avoid division by zero, e.g., by rho_zz(:,cell2)
          do k = 1, nVertLevels


!  WCS-parallel: we could pick this up in the last acoustic step (ruAvg) because we solve for all edges of owned cells

            ruAvg(k,iEdge) = ru(k,iEdge) + (ruAvg(k,iEdge) / float(ns))

            ru(k,iEdge) = ru(k,iEdge) + ru_p(k,iEdge)

            u(k,iEdge) = 2.*ru(k,iEdge)/(rho_zz(k,cell1)+rho_zz(k,cell2))
          enddo


! WCS-parallel: we likely only need this for owned cells

          !SHP-mtn
          flux = cf1*ru(1,iEdge) + cf2*ru(2,iEdge) + cf3*ru(3,iEdge)
          w(1,cell2) = w(1,cell2) - (zb(1,2,iEdge) + sign(1.0_RKIND,flux)*coef_3rd_order*zb3(1,2,iEdge))  &
                                 *flux/(cf1*rho_zz(1,cell2)+cf2*rho_zz(2,cell2)+cf3*rho_zz(3,cell2))
          w(1,cell1) = w(1,cell1) + (zb(1,1,iEdge) + sign(1.0_RKIND,flux)*coef_3rd_order*zb3(1,1,iEdge))  &
                                 *flux/(cf1*rho_zz(1,cell1)+cf2*rho_zz(2,cell1)+cf3*rho_zz(3,cell1))

          do k = 2, nVertLevels
            flux = (fzm(k)*ru(k,iEdge)+fzp(k)*ru(k-1,iEdge))
            w(k,cell2) = w(k,cell2) - (zb(k,2,iEdge)+sign(1.0_RKIND,flux)*coef_3rd_order*zb3(k,2,iEdge)) &
                                 *flux/(fzm(k)*rho_zz(k,cell2)+fzp(k)*rho_zz(k-1,cell2))
            w(k,cell1) = w(k,cell1) + (zb(k,1,iEdge)+sign(1.0_RKIND,flux)*coef_3rd_order*zb3(k,1,iEdge)) &
                                 *flux/(fzm(k)*rho_zz(k,cell1)+fzp(k)*rho_zz(k-1,cell1))
          enddo

!        end if

      enddo

   end subroutine atm_recover_large_step_variables

!---------------------------------------------------------------------------------------

   subroutine atm_advance_scalars( tend, s_old, s_new, diag, grid, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed scalar tendencies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(in) :: tend
      type (state_type), intent(in) :: s_old
      type (state_type), intent(inout) :: s_new
      type (diag_type), intent(in) :: diag
      type (mesh_type), intent(in) :: grid
      real (kind=RKIND) :: dt

      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_old, scalar_new, scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg, rho_edge, rho_zz, zgrid, kdiff
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, qv_init
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND), dimension( s_old % num_scalars, grid % nVertLevels ) :: flux_arr

      real (kind=RKIND), dimension( s_old % num_scalars, grid % nVertLevels + 1 ) :: wdtn
      integer :: nVertLevels

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND) :: coef_3rd_order

      real (kind=RKIND) :: h_theta_eddy_visc2, v_theta_eddy_visc2

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      integer, parameter :: hadv_opt = 2

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      coef_3rd_order = config_coef_3rd_order

      scalar_old  => s_old % scalars % array
      scalar_new  => s_new % scalars % array
      kdiff       => diag % kdiff % array
      deriv_two   => grid % deriv_two % array
      uhAvg       => diag % ruAvg % array
      dvEdge      => grid % dvEdge % array
      dcEdge      => grid % dcEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      scalar_tend => tend % scalars % array
      h_old       => s_old % rho_zz % array
      h_new       => s_new % rho_zz % array
      wwAvg       => diag % wwAvg % array
      areaCell    => grid % areaCell % array

      fnm         => grid % fzm % array
      fnp         => grid % fzp % array
      rdnw        => grid % rdzw % array
      meshScalingDel2 => grid % meshScalingDel2 % array
      meshScalingDel4 => grid % meshScalingDel4 % array

      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array

      nVertLevels = grid % nVertLevels

      h_theta_eddy_visc2 = config_h_theta_eddy_visc2
      v_theta_eddy_visc2 = config_v_theta_eddy_visc2
      rho_edge     => diag % rho_edge % array
      rho_zz       => s_new % rho_zz % array
      qv_init      => grid % qv_init % array
      zgrid        => grid % zgrid % array

#ifndef DO_PHYSICS
      scalar_tend = 0.  !  testing purposes - we have no sources or sinks
#endif

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts form scalar_old
      !
      !
      !  horizontal flux divergence, accumulate in scalar_tend

         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then  ! only for owned cells
  
               flux_arr(:,:) = 0.
               do i=1,nAdvCellsForEdge(iEdge)
                 iCell = advCellsForEdge(i,iEdge)
                 do k=1,grid % nVertLevels
                   scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge)
                   do iScalar=1,s_old % num_scalars
                     flux_arr(iScalar,k) = flux_arr(iScalar,k) + scalar_weight* scalar_new(iScalar,k,iCell)
                   end do
                 end do
               end do

               do k=1,grid % nVertLevels
               do iScalar=1,s_old % num_scalars
                     scalar_tend(iScalar,k,cell1) = scalar_tend(iScalar,k,cell1) &
                            - uhAvg(k,iEdge)*flux_arr(iScalar,k)/areaCell(cell1)
                     scalar_tend(iScalar,k,cell2) = scalar_tend(iScalar,k,cell2) &
                            + uhAvg(k,iEdge)*flux_arr(iScalar,k)/areaCell(cell2)
               end do
               end do

            end if
          end do

      !
      !  vertical flux divergence
      !

      ! zero fluxes at top and bottom

         wdtn(:,1) = 0.
         wdtn(:,grid % nVertLevels+1) = 0.

         if (config_scalar_vadv_order == 2) then

           do iCell=1,grid % nCellsSolve
             do k = 2, nVertLevels
                do iScalar=1,s_old % num_scalars
                  wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
                end do
             end do
             do k=1,grid % nVertLevels   ! Could be nVertLevelsSolve?
                do iScalar=1,s_old % num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
                end do
             end do
           end do

         else if ( config_scalar_vadv_order == 3 ) then

           do iCell=1,grid % nCellsSolve

             k = 2
             do iScalar=1,s_old % num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
             enddo
             
             do k=3,nVertLevels-1
               do iScalar=1,s_old % num_scalars
                 wdtn(iScalar,k) = flux3( scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  &
                                          scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell),  &
                                          wwAvg(k,iCell), coef_3rd_order )
               end do
             end do
             k = nVertLevels
             do iScalar=1,s_old % num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
             enddo

             do k=1,grid % nVertLevels   ! Could be nVertLevelsSolve?
                do iScalar=1,s_old % num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
                end do
             end do

           end do

         else if ( config_scalar_vadv_order == 4 ) then

           do iCell=1,grid % nCellsSolve

             k = 2
             do iScalar=1,s_old % num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
             enddo
             do k=3,nVertLevels-1
               do iScalar=1,s_old % num_scalars
                 wdtn(iScalar,k) = flux4( scalar_new(iScalar,k-2,iCell),scalar_new(iScalar,k-1,iCell),  &
                                          scalar_new(iScalar,k  ,iCell),scalar_new(iScalar,k+1,iCell), wwAvg(k,iCell) )
               end do
             end do
             k = nVertLevels
             do iScalar=1,s_old % num_scalars
               wdtn(iScalar,k) = wwAvg(k,iCell)*(fnm(k)*scalar_new(iScalar,k,iCell)+fnp(k)*scalar_new(iScalar,k-1,iCell))
             enddo

             do k=1,grid % nVertLevels   ! Could be nVertLevelsSolve?
                do iScalar=1,s_old % num_scalars
                  scalar_new(iScalar,k,iCell) = (   scalar_old(iScalar,k,iCell)*h_old(k,iCell) &
                        + dt*( scalar_tend(iScalar,k,iCell) -rdnw(k)*(wdtn(iScalar,k+1)-wdtn(iScalar,k)) ) )/h_new(k,iCell)
                end do
             end do

           end do
                                                                                        
         else 

           write(0,*) ' bad value for config_scalar_vadv_order - ', config_scalar_vadv_order

         end if

   end subroutine atm_advance_scalars

!---------------------------

   subroutine atm_advance_scalars_mono(tend, s_old, s_new, diag, grid, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
      implicit none

      type (tend_type),intent(in)     :: tend
      type (state_type),intent(inout) :: s_old
      type (state_type),intent(inout) :: s_new
      type (diag_type),intent(in)     :: diag
      type (mesh_type),intent(in)     :: grid
      real (kind=RKIND),intent(in)    :: dt

      type (block_type), pointer :: block
      integer :: i, iCell, iEdge, k, iScalar, cell1, cell2
      real (kind=RKIND) :: flux, scalar_weight

      real (kind=RKIND), dimension(:,:,:), pointer :: scalar_tend
      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      real (kind=RKIND), dimension(:,:), pointer :: uhAvg, h_old, h_new, wwAvg, rho_edge, rho_zz, zgrid, kdiff
      real (kind=RKIND), dimension(:), pointer :: dvEdge, dcEdge, areaCell, qv_init
      integer, dimension(:,:), pointer :: cellsOnEdge

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd

      type (field3DReal), pointer :: tempField
      type (field3DReal), target :: tempFieldTarget

      real (kind=RKIND), dimension( grid % nVertLevels, grid % nCells ) :: scalar_old, scalar_new
      real (kind=RKIND), dimension( grid % nVertLevels, grid % nCells ) :: s_max, s_min
      real (kind=RKIND), dimension( 2, grid % nVertLevels, grid % nCells ), target :: scale_arr

      integer, parameter :: SCALE_IN = 1, SCALE_OUT = 2

      real (kind=RKIND), dimension( grid % nVertLevels, grid % nEdges ) :: flux_arr
      real (kind=RKIND), dimension( grid % nVertLevels + 1, grid % nCells ) :: wdtn

      integer :: nVertLevels, num_scalars, icellmax, kmax

      real (kind=RKIND), dimension(:), pointer :: fnm, fnp, rdnw, meshScalingDel2, meshScalingDel4
      real (kind=RKIND) :: coef_3rd_order

      real (kind=RKIND) :: h_theta_eddy_visc2, v_theta_eddy_visc2

      real (kind=RKIND) :: flux3, flux4, flux_upwind
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3, scmin,scmax
      real (kind=RKIND) :: s_min_update, s_max_update, s_upwind, scale_factor

      integer, parameter :: hadv_opt = 2
      real (kind=RKIND), parameter :: eps=1.e-20
      logical, parameter :: debug_print = .false.

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
          ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

      block => grid % block

      coef_3rd_order = config_coef_3rd_order

      kdiff       => diag % kdiff % array
      deriv_two   => grid % deriv_two % array
      uhAvg       => diag % ruAvg % array
      dvEdge      => grid % dvEdge % array
      dcEdge      => grid % dcEdge % array
      cellsOnEdge => grid % cellsOnEdge % array
      scalar_tend => tend % scalars % array
      h_old       => s_old % rho_zz % array
      h_new       => s_new % rho_zz % array
      wwAvg       => diag % wwAvg % array
      areaCell    => grid % areaCell % array

      fnm         => grid % fzm % array
      fnp         => grid % fzp % array
      rdnw        => grid % rdzw % array
      meshScalingDel2 => grid % meshScalingDel2 % array
      meshScalingDel4 => grid % meshScalingDel4 % array

      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array

      nVertLevels = grid % nVertLevels

      h_theta_eddy_visc2 = config_h_theta_eddy_visc2
      v_theta_eddy_visc2 = config_v_theta_eddy_visc2
      rho_edge     => diag % rho_edge % array
      rho_zz       => s_new % rho_zz % array
      qv_init      => grid % qv_init % array
      zgrid        => grid % zgrid % array

#ifndef DO_PHYSICS
      scalar_tend = 0.  !  testing purposes - we have no sources or sinks
#endif

      !
      ! Update scalars for physics (i.e., what is in scalar_tend)
      !   we should probably move this to another routine called before mono advection ****
      !

      do iCell = 1,grid%nCellsSolve
      do k = 1, grid%nVertLevels
      do iScalar = 1,s_old%num_scalars
         s_old % scalars % array(iScalar,k,iCell) = s_old % scalars % array(iScalar,k,iCell)+dt*scalar_tend(iScalar,k,iCell) / h_old(k,iCell)
         scalar_tend(iScalar,k,iCell) = 0.
      end do
      end do
      end do

      !  halo exchange

      call mpas_dmpar_exch_halo_field(s_old % scalars)

      !
      ! Runge Kutta integration, so we compute fluxes from scalar_new values, update starts from scalar_old
      !

      ! do one scalar at a time

      num_scalars = 1

      do iScalar = 1, s_old % num_scalars
        write(0,*) ' mono transport for scalar ',iScalar

        do iCell = 1, grid%nCells
        do k=1, grid%nVertLevels
          scalar_old(k,iCell) = s_old % scalars % array(iScalar,k,iCell)
          scalar_new(k,iCell) = s_new % scalars % array(iScalar,k,iCell)
        end do
        end do

         if(debug_print) then
            scmin = scalar_old(1,1)
            scmax = scalar_old(1,1)
            do iCell = 1, grid%nCells
            do k=1, grid%nVertLevels
               scmin = min(scmin,scalar_old(k,iCell))
               scmax = max(scmax,scalar_old(k,iCell))
            enddo
            enddo
            write(0,*) ' scmin, scmin old in ',scmin,scmax

            scmin = scalar_new(1,1)
            scmax = scalar_new(1,1)
            do iCell = 1, grid%nCells
            do k=1, grid%nVertLevels
               scmin = min(scmin,scalar_new(k,iCell))
               scmax = max(scmax,scalar_new(k,iCell))
            enddo
            enddo
            write(0,*) ' scmin, scmin new in ',scmin,scmax
         end if


      !
      !  vertical flux divergence, and min and max bounds for flux limiter
      !

 
         do iCell=1,grid % nCellsSolve

          ! zero flux at top and bottom
           wdtn(1,iCell) = 0.
           wdtn(grid % nVertLevels+1,iCell) = 0.

           k = 1
           s_max(k,iCell) = max(scalar_old(1,iCell),scalar_old(2,iCell))
           s_min(k,iCell) = min(scalar_old(1,iCell),scalar_old(2,iCell))

           k = 2
           wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(k,iCell)+fnp(k)*scalar_new(k-1,iCell))
           s_max(k,iCell) = max(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
           s_min(k,iCell) = min(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
             
           do k=3,nVertLevels-1
              wdtn(k,iCell) = flux3( scalar_new(k-2,iCell),scalar_new(k-1,iCell),  &
                                     scalar_new(k  ,iCell),scalar_new(k+1,iCell),  &
                                     wwAvg(k,iCell), coef_3rd_order )
              s_max(k,iCell) = max(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
              s_min(k,iCell) = min(scalar_old(k-1,iCell),scalar_old(k,iCell),scalar_old(k+1,iCell))
           end do
 
           k = nVertLevels
           wdtn(k,iCell) = wwAvg(k,iCell)*(fnm(k)*scalar_new(k,iCell)+fnp(k)*scalar_new(k-1,iCell))
           s_max(k,iCell) = max(scalar_old(k,iCell),scalar_old(k-1,iCell))
           s_min(k,iCell) = min(scalar_old(k,iCell),scalar_old(k-1,iCell))

      ! pull s_min and s_max from the (horizontal) surrounding cells

           do i=1, grid % nEdgesOnCell % array(iCell)
             do k=1, grid % nVertLevels
               s_max(k,iCell) = max(s_max(k,iCell),scalar_old(k, grid % CellsOnCell % array(i,iCell)))
               s_min(k,iCell) = min(s_min(k,iCell),scalar_old(k, grid % CellsOnCell % array(i,iCell)))
             end do
           end do

         end do

      !
      !  horizontal flux divergence

         flux_arr(:,:) = 0.
         do iEdge=1,grid%nEdges

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then  ! only for owned cells
  
               do i=1,nAdvCellsForEdge(iEdge)
                 iCell = advCellsForEdge(i,iEdge)
                 do k=1,grid % nVertLevels
                   scalar_weight = uhAvg(k,iEdge)*(adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,uhAvg(k,iEdge))*adv_coefs_3rd(i,iEdge))
                   flux_arr(k,iEdge) = flux_arr(k,iEdge) + scalar_weight* scalar_new(k,iCell)
                 end do
               end do

            end if

          end do

!  vertical flux divergence for upwind update, we will put upwind update into scalar_new, and put factor of dt in fluxes

          do iCell = 1, grid % nCellsSolve

            k = 1
            scalar_new(k,iCell) = scalar_old(k,iCell)*h_old(k,iCell)

            do k = 2, nVertLevels
              scalar_new(k,iCell) = scalar_old(k,iCell)*h_old(k,iCell)
              flux_upwind = dt*(max(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k-1,iCell) + min(0.0_RKIND,wwAvg(k,iCell))*scalar_old(k,iCell))
              scalar_new(k-1,iCell) = scalar_new(k-1,iCell) - flux_upwind*rdnw(k-1)
              scalar_new(k  ,iCell) = scalar_new(k  ,iCell) + flux_upwind*rdnw(k)
              wdtn(k,iCell) = dt*wdtn(k,iCell) - flux_upwind
            end do

! scale_arr(SCALE_IN,:,:) and scale_arr(SCALE_OUT:,:) are used here to store the incoming and outgoing perturbation flux 
! contributions to the update:  first the vertical flux component, then the horizontal

            do k=1,nVertLevels
              scale_arr(SCALE_IN, k,iCell) = - rdnw(k)*(min(0.0_RKIND,wdtn(k+1,iCell))-max(0.0_RKIND,wdtn(k,iCell)))
              scale_arr(SCALE_OUT,k,iCell) = - rdnw(k)*(max(0.0_RKIND,wdtn(k+1,iCell))-min(0.0_RKIND,wdtn(k,iCell)))
            end do

          end do

!  horizontal flux divergence for upwind update

         !  upwind flux computation

         do iEdge=1,grid%nEdges
           cell1 = cellsOnEdge(1,iEdge)
           cell2 = cellsOnEdge(2,iEdge)
             if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then  ! only for owned cells
               do k=1,grid % nVertLevels
                 flux_upwind = grid % dvEdge % array(iEdge) * dt *   &
                        (max(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell1) + min(0.0_RKIND,uhAvg(k,iEdge))*scalar_old(k,cell2))
                 flux_arr(k,iEdge) = dt*flux_arr(k,iEdge) - flux_upwind
                 scalar_new(k,cell1) = scalar_new(k,cell1) - flux_upwind / areaCell(cell1)
                 scalar_new(k,cell2) = scalar_new(k,cell2) + flux_upwind / areaCell(cell2)

                 scale_arr(SCALE_OUT,k,cell1) = scale_arr(SCALE_OUT,k,cell1) - max(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell1)
                 scale_arr(SCALE_IN, k,cell1) = scale_arr(SCALE_IN, k,cell1) - min(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell1)
                 scale_arr(SCALE_OUT,k,cell2) = scale_arr(SCALE_OUT,k,cell2) + min(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell2)
                 scale_arr(SCALE_IN, k,cell2) = scale_arr(SCALE_IN, k,cell2) + max(0.0_RKIND,flux_arr(k,iEdge)) / areaCell(cell2)

               end do
             end if
          end do

!  next, the limiter

          do iCell = 1, grid % nCellsSolve
            do k = 1, nVertLevels
               s_min_update = (scalar_new(k,iCell)+scale_arr(SCALE_OUT,k,iCell))/h_new(k,iCell)
               s_max_update = (scalar_new(k,iCell)+scale_arr(SCALE_IN,k,iCell))/h_new(k,iCell)
               s_upwind = scalar_new(k,iCell)/h_new(k,iCell)
           
               scale_factor = (s_max(k,iCell)-s_upwind)/(s_max_update-s_upwind+eps)
               scale_arr(SCALE_IN,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

               scale_factor = (s_upwind-s_min(k,iCell))/(s_upwind-s_min_update+eps)
               scale_arr(SCALE_OUT,k,iCell) = min( 1.0_RKIND, max( 0.0_RKIND, scale_factor) )

            end do
          end do
!
!  communicate scale factors here
!
!
!  WCS_halo_opt_2 - communicate only first halo row in these next two exchanges
!
      tempField => tempFieldTarget

      tempField % block => block
      tempField % dimSizes(1) = 2
      tempField % dimSizes(2) = grid % nVertLevels
      tempField % dimSizes(3) = grid % nCells
      tempField % sendList => block % parinfo % cellsToSend
      tempField % recvList => block % parinfo % cellsToRecv
      tempField % copyList => block % parinfo % cellsToCopy
      tempField % prev => null()
      tempField % next => null()

      tempField % array => scale_arr
      call mpas_dmpar_exch_halo_field(tempField, (/ 1 /))

!
!  rescale the fluxes
!
            do iEdge = 1, grid % nEdges
               cell1 = grid % cellsOnEdge % array(1,iEdge)
               cell2 = grid % cellsOnEdge % array(2,iEdge)
               if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then
                  do k = 1, nVertLevels
                     flux = flux_arr(k,iEdge)
                     flux = max(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k,cell1), scale_arr(SCALE_IN, k,cell2)) &
                          + min(0.0_RKIND,flux) * min(scale_arr(SCALE_IN, k,cell1), scale_arr(SCALE_OUT,k,cell2))
                     flux_arr(k,iEdge) = flux
                  end do
               end if
            end do
 
       ! rescale the vertical flux
 
            do iCell=1,grid % nCells
               do k = 2, nVertLevels
                  flux =  wdtn(k,iCell)
                  flux = max(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k-1,iCell), scale_arr(SCALE_IN,k  ,iCell)) &
                       + min(0.0_RKIND,flux) * min(scale_arr(SCALE_OUT,k  ,iCell), scale_arr(SCALE_IN,k-1,iCell))
                  wdtn(k,iCell) = flux
               end do
            end do
!
!  do the scalar update now that we have the fluxes
!
         do iEdge=1,grid%nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCellsSolve .or. cell2 <= grid%nCellsSolve) then  ! only for owned cells
               do k=1,grid % nVertLevels
                  scalar_new(k,cell1) = scalar_new(k,cell1) - flux_arr(k,iEdge)/areaCell(cell1)
                  scalar_new(k,cell2) = scalar_new(k,cell2) + flux_arr(k,iEdge)/areaCell(cell2)
               end do
            end if
         end do

         do iCell=1,grid % nCellsSolve
            do k=1,grid % nVertLevels
              scalar_new(k,iCell) = (   scalar_new(k,iCell)  &
                  + (-rdnw(k)*(wdtn(k+1,iCell)-wdtn(k,iCell)) ) )/h_new(k,iCell)
            end do
         end do

         if(debug_print) then

            scmin = scalar_new(1,1)
            scmax = scalar_new(1,1)
            do iCell = 1, grid%nCellsSolve
            do k=1, grid%nVertLevels
               scmax = max(scmax,scalar_new(k,iCell))
               scmin = min(scmin,scalar_new(k,iCell))
               if(s_max(k,iCell) < scalar_new(k,iCell)) then
                  write(32,*) ' over - k,iCell,s_min,s_max,scalar_new ',k,iCell,s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)
               end if
               if(s_min(k,iCell) > scalar_new(k,iCell)) then
                  write(32,*) ' under - k,iCell,s_min,s_max,scalar_new ',k,iCell,s_min(k,iCell),s_max(k,iCell),scalar_new(k,iCell)
               end if
            enddo
            enddo
            write(0,*) ' scmin, scmax new out ',scmin,scmax
            write(0,*) ' icell_min, k_min ',icellmax, kmax

         end if

         do iCell = 1, grid%nCells
         do k=1, grid%nVertLevels
            s_new % scalars % array(iScalar,k,iCell) = max(0.0_RKIND,scalar_new(k,iCell))
         end do
         end do

      end do !  loop over scalars

   end subroutine atm_advance_scalars_mono

!----

   subroutine atm_compute_dyn_tend(tend, s, diag, grid, rk_step, dt)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute height and normal wind tendencies, as well as diagnostic variables
   !
   ! Input: s - current model state
   !        grid - grid metadata
   !
   ! Output: tend - computed diagnostics (parallel velocities, v; mass fluxes, rv; 
   !                circulation; vorticity; and kinetic energy, ke) and the 
   !                tendencies for height (h) and u (u)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      type (tend_type), intent(inout) :: tend
      type (state_type), intent(in) :: s
      type (diag_type), intent(in) :: diag
      type (mesh_type), intent(in) :: grid
      integer, intent(in) :: rk_step
      real (kind=RKIND), intent(in) :: dt

      logical, parameter :: rk_diffusion = .false.

      integer :: iEdge, iCell, iVertex, k, cell1, cell2, vertex1, vertex2, eoe, i, j, iq
      real (kind=RKIND) :: flux, workpv

      integer :: nCells, nEdges, nVertices, nVertLevels, nCellsSolve
      real (kind=RKIND) :: h_mom_eddy_visc2,   v_mom_eddy_visc2,   h_mom_eddy_visc4
      real (kind=RKIND) :: h_theta_eddy_visc2, v_theta_eddy_visc2, h_theta_eddy_visc4
      real (kind=RKIND) :: u_diffusion
      real (kind=RKIND), dimension(:), pointer ::  fEdge, dvEdge, dcEdge, areaCell, areaTriangle, meshScalingDel2, meshScalingDel4
      real (kind=RKIND), dimension(:,:), pointer :: weightsOnEdge, zgrid, rho_edge, rho_zz, ru, u, v, tend_u, &
                                                    circulation, divergence, vorticity, ke, pv_edge, theta_m, rw, tend_rho, &
                                                    rt_diabatic_tend, tend_theta, tend_w, w, cqw, rb, rr, pp, pressure_b, zz, zx, cqu, & 
                                                    h_divergence, kdiff

      real (kind=RKIND), dimension(:,:), pointer :: tend_u_euler, tend_w_euler, tend_theta_euler

      real (kind=RKIND), dimension(:,:,:), pointer :: deriv_two
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, edgesOnCell, edgesOnEdge
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      real (kind=RKIND), dimension( grid % nVertLevels + 1 ) :: wduz, wdwz, wdtz, dpzx
      real (kind=RKIND), dimension( grid % nVertLevels ) :: u_mix, ru_edge_w, q
      real (kind=RKIND) :: theta_turb_flux, z1, z2, z3, z4, zm, z0, zp, r
      real (kind=RKIND) :: d2fdx2_cell1, d2fdx2_cell2

      integer, dimension(:,:), pointer :: advCellsForEdge
      integer, dimension(:), pointer :: nAdvCellsForEdge
      real (kind=RKIND), dimension(:,:), pointer :: adv_coefs, adv_coefs_3rd
      real (kind=RKIND) :: scalar_weight

      real (kind=RKIND), dimension(:), pointer :: rdzu, rdzw, fzm, fzp, qv_init
      real (kind=RKIND), dimension(:,:), pointer :: t_init 

      real (kind=RKIND), dimension(:,:), pointer :: cpr, cpl, pzp, pzm
      integer :: kr, kl

      real (kind=RKIND), allocatable, dimension(:,:) :: divergence_ru, qtot 
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_theta, delsq_divergence
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_u
      real (kind=RKIND), allocatable, dimension(:,:) :: delsq_circulation, delsq_vorticity
      real (kind=RKIND) :: cf1, cf2, cf3, pr, pl
      real (kind=RKIND) :: prandtl_inv

!      logical, parameter :: debug = .true.
      logical, parameter :: debug = .false.

      !SHP-curvature
      logical, parameter :: curvature = .true.
      real (kind=RKIND) :: r_earth
      real (kind=RKIND), dimension(:,:), pointer :: ur_cell, vr_cell

      real (kind=RKIND), parameter :: c_s = 0.125
!      real (kind=RKIND), parameter :: c_s = 0.25
      real (kind=RKIND), dimension( grid % nVertLevels ) :: d_diag, d_off_diag, flux_arr
      real (kind=RKIND), dimension(:,:), pointer :: defc_a, defc_b
      logical :: delsq_horiz_mixing, newpx


      real (kind=RKIND) :: coef_3rd_order

      real (kind=RKIND) :: flux3, flux4
      real (kind=RKIND) :: q_im2, q_im1, q_i, q_ip1, ua, coef3

      flux4(q_im2, q_im1, q_i, q_ip1, ua) =                     &
                ua*( 7.*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0

      flux3(q_im2, q_im1, q_i, q_ip1, ua, coef3) =              &
                flux4(q_im2, q_im1, q_i, q_ip1, ua) +           &
                coef3*abs(ua)*((q_ip1 - q_im2)-3.*(q_i-q_im1))/12.0

!-----------

      !SHP-curvature
      r_earth = grid % sphere_radius
      ur_cell => diag % uReconstructZonal % array
      vr_cell => diag % uReconstructMeridional % array

      coef_3rd_order = config_coef_3rd_order

      rho_zz       => s % rho_zz % array
      rho_edge     => diag % rho_edge % array
      rb           => diag % rho_base % array
      rr           => diag % rho_p % array
      u            => s % u % array
      v            => diag % v % array
      kdiff        => diag % kdiff % array
      ru           => diag % ru % array
      w            => s % w % array
      rw           => diag % rw % array
      theta_m      => s % theta_m % array
      circulation  => diag % circulation % array
      divergence   => diag % divergence % array
      vorticity    => diag % vorticity % array
      ke           => diag % ke % array
      pv_edge      => diag % pv_edge % array
      pp           => diag % pressure_p % array
      pressure_b   => diag % pressure_base % array
      h_divergence => diag % h_divergence % array

      pzp          => grid % pzp % array
      pzm          => grid % pzm % array


      weightsOnEdge     => grid % weightsOnEdge % array
      cellsOnEdge       => grid % cellsOnEdge % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnCell       => grid % edgesOnCell % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      fEdge             => grid % fEdge % array
      deriv_two         => grid % deriv_two % array
      zz                => grid % zz % array
      zx                => grid % zx % array

      defc_a             => grid % defc_a % array
      defc_b             => grid % defc_b % array

      meshScalingDel2    => grid % meshScalingDel2 % array
      meshScalingDel4    => grid % meshScalingDel4 % array


      tend_u      => tend % u % array
      tend_theta  => tend % theta_m % array
      tend_w      => tend % w % array
      tend_rho    => tend % rho_zz % array
      rt_diabatic_tend  => tend % rt_diabatic_tend % array

      tend_u_euler      => tend % u_euler % array
      tend_theta_euler  => tend % theta_euler % array
      tend_w_euler      => tend % w_euler % array

      t_init      => grid % t_init % array
      qv_init     => grid % qv_init % array

      rdzu        => grid % rdzu % array
      rdzw        => grid % rdzw % array
      fzm         => grid % fzm % array
      fzp         => grid % fzp % array
      zgrid       => grid % zgrid % array
      cqw         => diag % cqw % array
      cqu         => diag % cqu % array

      cpr         => grid % cpr % array
      cpl         => grid % cpl % array
      newpx       = config_newpx

      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertLevels = grid % nVertLevels
      nVertices   = grid % nVertices
      nCellsSolve = grid % nCellsSolve

      h_mom_eddy_visc2 = config_h_mom_eddy_visc2
!     h_mom_eddy_visc4 = config_h_mom_eddy_visc4
      v_mom_eddy_visc2 = config_v_mom_eddy_visc2
      h_theta_eddy_visc2 = config_h_theta_eddy_visc2
!     h_theta_eddy_visc4 = config_h_theta_eddy_visc4
      v_theta_eddy_visc2 = config_v_theta_eddy_visc2

      nAdvCellsForEdge => grid % nAdvCellsForEdge % array
      advCellsForEdge => grid % advCellsForEdge % array
      adv_coefs => grid % adv_coefs % array
      adv_coefs_3rd => grid % adv_coefs_3rd % array

      prandtl_inv = 1.0_RKIND/prandtl

      !
      ! Compute u (normal) velocity tendency for each edge (cell face)
      !

      write(0,*) ' rk_step in compute_dyn_tend ',rk_step


      delsq_horiz_mixing = .false.
      if (config_horiz_mixing == "2d_smagorinsky" .and. (rk_step == 1 .or. rk_diffusion)) then
         do iCell = 1, nCells
            d_diag(:) = 0.
            d_off_diag(:) = 0.
            do iEdge = 1, grid % nEdgesOnCell % array (iCell)
               do k=1, nVertLevels
                  d_diag(k)     = d_diag(k)     + defc_a(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                - defc_b(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
                  d_off_diag(k) = d_off_diag(k) + defc_b(iEdge,iCell)*u(k,EdgesOnCell(iEdge,iCell))  &
                                                + defc_a(iEdge,iCell)*v(k,EdgesOnCell(iEdge,iCell))
               end do
            end do
            do k=1, nVertLevels
               kdiff(k,iCell) = (c_s * config_len_disp)**2 * sqrt(d_diag(k)**2 + d_off_diag(k)**2)
               kdiff(k,iCell) = min(kdiff(k,iCell),(0.01*config_len_disp**2)/dt)
            end do
         end do
!ldf (2012-10-10):
         h_mom_eddy_visc4   = config_visc4_2dsmag * config_len_disp**3
         h_theta_eddy_visc4 = h_mom_eddy_visc4
         delsq_horiz_mixing = .true.
         write(0,*) '... config_visc4_2dsmag = ', config_visc4_2dsmag
         write(0,*) '... h_mom_eddy_visc4    = ', h_mom_eddy_visc4
         write(0,*) '... h_theta_eddy_visc4  = ', h_theta_eddy_visc4
      else if ( config_horiz_mixing == "2d_fixed") then
         h_mom_eddy_visc4   = config_h_mom_eddy_visc4
         h_theta_eddy_visc4 = config_h_theta_eddy_visc4
         delsq_horiz_mixing = .true.
!ldf (2012-10-10):
      end if

      tend_u(:,:) = 0.0

      cf1 = grid % cf1 % scalar
      cf2 = grid % cf2 % scalar
      cf3 = grid % cf3 % scalar

      !  tendency for density
      !  divergence_ru may calculated in the diagnostic subroutine - it is temporary
      allocate(divergence_ru(nVertLevels, nCells+1))
      allocate(qtot(nVertLevels, nCells+1))

      divergence_ru(:,:) = 0.0
      h_divergence(:,:) = 0.
      do iEdge=1,grid % nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
           flux = ru(k,iEdge)*dvEdge(iEdge)
           divergence_ru(k,cell1) = divergence_ru(k,cell1) + flux
           divergence_ru(k,cell2) = divergence_ru(k,cell2) - flux
         end do
      end do

      qtot(:,:)=0.
      do iCell = 1,nCells
        r = 1.0 / areaCell(iCell)
        do k = 1,nVertLevels
           divergence_ru(k,iCell) = divergence_ru(k,iCell) * r
           h_divergence(k,iCell) = divergence_ru(k,iCell)
           tend_rho(k,iCell) = -divergence_ru(k,iCell)-rdzw(k)*(rw(k+1,iCell)-rw(k,iCell))

           do iq = s % moist_start, s % moist_end
              qtot(k,iCell) = qtot(k,iCell) + s % scalars % array (iq, k, iCell)
           end do

        end do
      end do    

!****       u dyn tend

         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

           if(newpx)  then

               k = 1
               pr  = cpr(k,iEdge)*pp(k,cell2)+cpr(k+1,iEdge)*pp(k+1,cell2)+cpr(k+2,iEdge)*pp(k+2,cell2)
               pl  = cpl(k,iEdge)*pp(k,cell1)+cpl(k+1,iEdge)*pp(k+1,cell1)+cpl(k+2,iEdge)*pp(k+2,cell1)
               tend_u(k,iEdge) =  - cqu(k,iEdge)*2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)

               do k=2,nVertLevels

                  kr = min(nVertLevels,k+ nint(.5-sign(0.5_RKIND,zx(k,iEdge)+zx(k+1,iEdge))))
                  kl = min(nVertLevels,2*k+1-kr)

                  pr = pp(k,cell2)+.5*(zgrid(k   ,cell1)+zgrid(k +1,cell1)-zgrid(k ,cell2)-zgrid(k +1,cell2))   &
                                     /(zgrid(kr+1,cell2)-zgrid(kr-1,cell2))*( pp(kr,cell2)-pp   (kr-1,cell2))
                  pl = pp(k,cell1)+.5*(zgrid(k   ,cell2)+zgrid(k +1,cell2)-zgrid(k ,cell1)-zgrid(k +1,cell1))   &
                                     /(zgrid(kl+1,cell1)-zgrid(kl-1,cell1))*( pp(kl,cell1)-pp   (kl-1,cell1))
                  tend_u(k,iEdge) =  - cqu(k,iEdge)*2./(zz(k,cell1)+zz(k,cell2))*(pr-pl)/dcEdge(iEdge)

               end do

            else
               k = 1
!!              dpzx(k) = .5*zx(k,iEdge)*(cf1*(pp(k  ,cell2)+pp(k  ,cell1))   &
!!                                       +cf2*(pp(k+1,cell2)+pp(k+1,cell1))   &
!!                                       +cf3*(pp(k+2,cell2)+pp(k+2,cell1)))

               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                            *(pzm(k,cell2)*(pp(k+1,cell2)-pp(k,cell2))    &
                             +pzm(k,cell1)*(pp(k+1,cell1)-pp(k,cell1))    &
                             +pzp(k,cell2)*(pp(k+2,cell2)-pp(k,cell2))    &
                             +pzp(k,cell1)*(pp(k+2,cell1)-pp(k,cell1))) 
  
               do k = 2, nVertLevels-1

!!              dpzx(k) = .5*zx(k,iEdge)*(fzm(k)*(pp(k  ,cell2)+pp(k  ,cell1))  &
!!                                   +fzp(k)*(pp(k-1,cell2)+pp(k-1,cell1)))

                  dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                             *(pzp(k,cell2)*(pp(k+1,cell2)-pp(k  ,cell2))    &
                              +pzm(k,cell2)*(pp(k  ,cell2)-pp(k-1,cell2))    &
                              +pzp(k,cell1)*(pp(k+1,cell1)-pp(k  ,cell1))    &
                              +pzm(k,cell1)*(pp(k  ,cell1)-pp(k-1,cell1)))   

               end do

               k = nVertLevels
               dpzx(k) = .25*(zx(k,iEdge)+zx(k+1,iEdge))                  &
                          *(pzm(k,cell2)*(pp(k  ,cell2)-pp(k-1,cell2))    &
                           +pzm(k,cell1)*(pp(k  ,cell1)-pp(k-1,cell1)))   

!!               dpzx(nVertLevels+1) = 0.

               do k=1,nVertLevels

!!                  tend_u(k,iEdge) =  - cqu(k,iEdge)*( (pp(k,cell2)/zz(k,cell2) - pp(k,cell1)/zz(k,cell1))  &
!!                                                   /  dcEdge(iEdge) - rdzw(k)*(dpzx(k+1)-dpzx(k)) )

                  tend_u(k,iEdge) =  - cqu(k,iEdge)*((pp(k,cell2)-pp(k,cell1))/dcEdge(iEdge)   &
                                          - dpzx(k) ) / (.5*(zz(k,cell2)+zz(k,cell1)))
               end do

            end if

            wduz(1) = 0.
            k = 2
            wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2) )*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))  
            do k=3,nVertLevels-1
               wduz(k) = flux3( u(k-2,iEdge),u(k-1,iEdge),u(k,iEdge),u(k+1,iEdge),0.5*(rw(k,cell1)+rw(k,cell2)), 1.0_RKIND )
            end do
            k = nVertLevels
            wduz(k) =  0.5*( rw(k,cell1)+rw(k,cell2) )*(fzm(k)*u(k,iEdge)+fzp(k)*u(k-1,iEdge))  

            wduz(nVertLevels+1) = 0.

            do k=1,nVertLevels
!               tend_u(k,iEdge) =  - cqu(k,iEdge)*( (pp(k,cell2)/zz(k,cell2) - pp(k,cell1)/zz(k,cell1))  &
!                                                /  dcEdge(iEdge) - rdzw(k)*(dpzx(k+1)-dpzx(k)) )
               tend_u(k,iEdge) = tend_u(k,iEdge) - rdzw(k)*(wduz(k+1)-wduz(k)) 
            end do

               q(:) = 0.0
               do j = 1,nEdgesOnEdge(iEdge)
                 eoe = edgesOnEdge(j,iEdge)
                 do k=1,nVertLevels
                   workpv = 0.5 * (pv_edge(k,iEdge) + pv_edge(k,eoe))
                   q(k) = q(k) + weightsOnEdge(j,iEdge) * u(k,eoe) * workpv * rho_edge(k,eoe)
                 end do
               end do

            do k=1,nVertLevels
               tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge)* (q(k) - (ke(k,cell2) - ke(k,cell1))       &
                                                                    / dcEdge(iEdge))                            &
                                                - u(k,iEdge)*0.5*(divergence_ru(k,cell1)+divergence_ru(k,cell2)) 
               !SHP-curvature
               if (curvature) then
               tend_u(k,iEdge) = tend_u(k,iEdge) &
                                - 2.*omega*cos(grid % angleEdge % array(iEdge))*cos(grid % latEdge % array(iEdge))  &
                                  *rho_edge(k,iEdge)*.25*(w(k,cell1)+w(k+1,cell1)+w(k,cell2)+w(k+1,cell2))          & 
                                - u(k,iEdge)*.25*(w(k+1,cell1)+w(k,cell1)+w(k,cell2)+w(k+1,cell2))                  &
                                  *rho_edge(k,iEdge)/r_earth
               !old-err.
               !tend_u(k,iEdge) = tend_u(k,iEdge) &
               !                 - 2.*omega_e*cos(grid % angleEdge % array(iEdge))*cos(grid % latEdge % array(iEdge))  &
               !                   *.25*(rw(k,cell1)+rw(k+1,cell1)+rw(k,cell2)+rw(k+1,cell2))                          &
               !                 - u(k,iEdge)*.25*(rw(k+1,cell1)+rw(k,cell1)+rw(k,cell2)+rw(k+1,cell2))/r_earth
               end if
            end do
         end do

      deallocate(divergence_ru)

      !
      !  horizontal mixing for u
      !

      if (rk_step == 1 .or. rk_diffusion) then

      tend_u_euler = 0.

      if (delsq_horiz_mixing) then

        if ((h_mom_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then
           do iEdge=1,grid % nEdgesSolve
              cell1 = cellsOnEdge(1,iEdge)
              cell2 = cellsOnEdge(2,iEdge)
              vertex1 = verticesOnEdge(1,iEdge)
              vertex2 = verticesOnEdge(2,iEdge)
  
              do k=1,nVertLevels
  
                 !
                 ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                 !                    only valid for h_mom_eddy_visc2 == constant
                 !
                 u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                                -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge),0.25*dcEdge(iEdge))
                 u_diffusion = rho_edge(k,iEdge)*h_mom_eddy_visc2 * u_diffusion
                 u_diffusion = u_diffusion * meshScalingDel2(iEdge)
  
!                 tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
                 tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + u_diffusion
              end do
           end do

        else if ( config_horiz_mixing == "2d_smagorinsky") then

           do iEdge=1,grid % nEdgesSolve
              cell1 = cellsOnEdge(1,iEdge)
              cell2 = cellsOnEdge(2,iEdge)
              vertex1 = verticesOnEdge(1,iEdge)
              vertex2 = verticesOnEdge(2,iEdge)

              do k=1,nVertLevels
                 !
                 ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
                 !                    only valid for h_mom_eddy_visc2 == constant
                 !
                 u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                                -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge),0.25*dcEdge(iEdge))
                 u_diffusion = rho_edge(k,iEdge)* 0.5*(kdiff(k,cell1)+kdiff(k,cell2)) * u_diffusion
                 u_diffusion = u_diffusion * meshScalingDel2(iEdge)

!                 tend_u(k,iEdge) = tend_u(k,iEdge) + u_diffusion
                 tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + u_diffusion
              end do
           end do
        end if

      end if ! delsq_horiz_mixing for u

!ldf (2012-10-10): Modified loop below to allow hyper-diffusion when 2d_smagorinsky is set to true.
!     if ((h_mom_eddy_visc4 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then
      if ((h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
          (h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then

         allocate(delsq_divergence(nVertLevels, nCells+1))
         allocate(delsq_u(nVertLevels, nEdges+1))
         allocate(delsq_circulation(nVertLevels, nVertices+1))
         allocate(delsq_vorticity(nVertLevels, nVertices+1))

         delsq_u(:,:) = 0.0

         do iEdge=1,grid % nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               u_diffusion =   ( divergence(k,cell2)  - divergence(k,cell1) ) / dcEdge(iEdge)  &
                              -( vorticity(k,vertex2) - vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge))

               delsq_u(k,iEdge) = delsq_u(k,iEdge) + u_diffusion
            end do
         end do

         delsq_circulation(:,:) = 0.0
         do iEdge=1,nEdges
               do k=1,nVertLevels
                  delsq_circulation(k,verticesOnEdge(1,iEdge)) = delsq_circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * delsq_u(k,iEdge)
                  delsq_circulation(k,verticesOnEdge(2,iEdge)) = delsq_circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * delsq_u(k,iEdge)
               end do
         end do
         do iVertex=1,nVertices
            r = 1.0 / areaTriangle(iVertex)
            do k=1,nVertLevels
               delsq_vorticity(k,iVertex) = delsq_circulation(k,iVertex) * r
            end do
         end do

         delsq_divergence(:,:) = 0.0
         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
               do k=1,nVertLevels
                 delsq_divergence(k,cell1) = delsq_divergence(k,cell1) + delsq_u(k,iEdge)*dvEdge(iEdge)
                 delsq_divergence(k,cell2) = delsq_divergence(k,cell2) - delsq_u(k,iEdge)*dvEdge(iEdge)
               end do
         end do
         do iCell = 1,nCells
            r = 1.0 / areaCell(iCell)
            do k = 1,nVertLevels
               delsq_divergence(k,iCell) = delsq_divergence(k,iCell) * r
            end do
         end do

         do iEdge=1,grid % nEdgesSolve
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            vertex1 = verticesOnEdge(1,iEdge)
            vertex2 = verticesOnEdge(2,iEdge)

            do k=1,nVertLevels

               !
               ! Compute diffusion, computed as \nabla divergence - k \times \nabla vorticity
               !                    only valid for h_mom_eddy_visc4 == constant
               !
               u_diffusion =  rho_edge(k,iEdge) * (  ( delsq_divergence(k,cell2)  - delsq_divergence(k,cell1) ) / dcEdge(iEdge)  &
                           -( delsq_vorticity(k,vertex2) - delsq_vorticity(k,vertex1) ) / max(dvEdge(iEdge), 0.25*dcEdge(iEdge)) &
                                                  )

!               tend_u(k,iEdge) = tend_u(k,iEdge) - h_mom_eddy_visc4 * u_diffusion
               u_diffusion = u_diffusion * meshScalingDel4(iEdge)
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) - h_mom_eddy_visc4 * u_diffusion
            end do
         end do

         deallocate(delsq_divergence)
         deallocate(delsq_u)
         deallocate(delsq_circulation)
         deallocate(delsq_vorticity)

      end if

      !
      !  vertical mixing for u - 2nd order 
      !
      if ( v_mom_eddy_visc2 > 0.0 ) then

         if (config_mix_full) then

         do iEdge=1,grid % nEdgesSolve

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=2,nVertLevels-1

               z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
               z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
               z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
               z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)

!               tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
!                                  (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      &
!                                 -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                  (u(k+1,iEdge)-u(k  ,iEdge))/(zp-z0)                      &
                                 -(u(k  ,iEdge)-u(k-1,iEdge))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

         do iEdge=1,grid % nEdgesSolve

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            do k=1,nVertLevels
#ifdef ROTATED_GRID
              u_mix(k) = u(k,iEdge) - grid % u_init % array(k) * sin( grid % angleEdge % array(iEdge) )
#else
              u_mix(k) = u(k,iEdge) - grid % u_init % array(k) * cos( grid % angleEdge % array(iEdge) )
#endif
            end do

            do k=2,nVertLevels-1

               z1 = 0.5*(zgrid(k-1,cell1)+zgrid(k-1,cell2))
               z2 = 0.5*(zgrid(k  ,cell1)+zgrid(k  ,cell2))
               z3 = 0.5*(zgrid(k+1,cell1)+zgrid(k+1,cell2))
               z4 = 0.5*(zgrid(k+2,cell1)+zgrid(k+2,cell2))

               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)

               tend_u_euler(k,iEdge) = tend_u_euler(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
                                  (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      &
                                 -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm))
!               tend_u(k,iEdge) = tend_u(k,iEdge) + rho_edge(k,iEdge) * v_mom_eddy_visc2*(  &
!                                  (u_mix(k+1)-u_mix(k  ))/(zp-z0)                      &
!                                 -(u_mix(k  )-u_mix(k-1))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

         end if

      end if

      end if ! (rk_step 1 test for computing mixing terms)

!  add in mixing for u

      do iEdge=1,grid % nEdgesSolve
         do k=1,nVertLevels
            tend_u(k,iEdge) = tend_u(k,iEdge) + tend_u_euler(k,iEdge)
         end do
      end do

!----------- rhs for w

      tend_w(:,:) = 0.

      !
      !  horizontal advection for w
      !

      if (config_w_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               do k=2,grid % nVertLevels
                  flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge) ) &
                                        *(w(k,cell1) + w(k,cell2))*0.5 
                  tend_w(k,cell1) = tend_w(k,cell1) - flux
                  tend_w(k,cell2) = tend_w(k,cell2) + flux
               end do
            end if
         end do

      else if (config_w_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,grid % nVertLevels
                 ru_edge_w(k) = fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)
               end do

               flux_arr(:) = 0.
               do i=1,nAdvCellsForEdge(iEdge)
                 iCell = advCellsForEdge(i,iEdge)
                 do k=2,grid % nVertLevels
                   scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,ru_edge_w(k))*adv_coefs_3rd(i,iEdge)
                   flux_arr(k) = flux_arr(k) + scalar_weight* w(k,iCell)
                 end do
               end do

               do k=1,grid % nVertLevels
                  tend_w(k,cell1) = tend_w(k,cell1) - ru_edge_w(k)*flux_arr(k)
                  tend_w(k,cell2) = tend_w(k,cell2) + ru_edge_w(k)*flux_arr(k)
               end do

!               do k=2,grid % nVertLevels
!
!                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * w(k,cell1)
!                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * w(k,cell2)
!                  do i=1, grid % nEdgesOnCell % array (cell1)
!                     if ( grid % CellsOnCell % array (i,cell1) <= grid%nCells) &
!                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * w(k,grid % CellsOnCell % array (i,cell1))
!                  end do
!                  do i=1, grid % nEdgesOnCell % array (cell2)
!                     if ( grid % CellsOnCell % array (i,cell2) <= grid%nCells) &
!                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * w(k,grid % CellsOnCell % array (i,cell2))
!                  end do
!
!  3rd order stencil
!                  if( u(k,iEdge)+u(k-1,iEdge) > 0) then
!                     flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge))*(  &
!                                             0.5*(w(k,cell1) + w(k,cell2))                 &
!                                            -(dcEdge(iEdge) **2) * (d2fdx2_cell1) / 6. )
!                  else
!                     flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge))*(  &
!                                             0.5*(w(k,cell1) + w(k,cell2))                 &
!                                            -(dcEdge(iEdge) **2) * (d2fdx2_cell2) / 6. )
!                  end if
!
!                  tend_w(k,cell1) = tend_w(k,cell1) - flux
!                  tend_w(k,cell2) = tend_w(k,cell2) + flux
!
!               end do

            end if
         end do

      else  if (config_w_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,grid % nVertLevels

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * w(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * w(k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) <= grid%nCells) &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * w(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) <= grid%nCells) &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * w(k,grid % CellsOnCell % array (i,cell2))
                  end do

                  flux = dvEdge(iEdge) * (fzm(k)*ru(k,iEdge) + fzp(k)*ru(k-1,iEdge)) * (  &
                                          0.5*(w(k,cell1) + w(k,cell2))                   &
                                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                  tend_w(k,cell1) = tend_w(k,cell1) - flux
                  tend_w(k,cell2) = tend_w(k,cell2) + flux
               end do

            end if

         end do
      end if

      !SHP-curvature
      if (curvature) then

         do iCell = 1, grid % nCellsSolve
            do k=2,nVertLevels
               tend_w(k,iCell) = tend_w(k,iCell) + (rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))*          &
                                         ( (fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))**2.             &
                                          +(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell))**2. )/r_earth   &
                                   + 2.*omega*cos(grid % latCell % array(iCell))                               &
                                          *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))                 &
                                          *(rho_zz(k,iCell)*fzm(k)+rho_zz(k-1,iCell)*fzp(k))

               !old_err.
               !tend_w(k,iCell) = tend_w(k,iCell) &
               !                 + rho_zz(k,iCell)*( (fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))**2.          &
               !                                 +(fzm(k)*vr_cell(k,iCell)+fzp(k)*vr_cell(k-1,iCell))**2. )/r_earth   &
               !                 + 2.*omega_e*cos(grid % latCell % array(iCell))*rho_zz(k,iCell)      &
               !                     *(fzm(k)*ur_cell(k,iCell)+fzp(k)*ur_cell(k-1,iCell))
            end do
         end do

      end if

      !
      !  horizontal mixing for w - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      !  Note: we are using quite a bit of the theta code here - could be combined later???

      if (rk_step == 1 .or. rk_diffusion) then

      tend_w_euler = 0.

      if (delsq_horiz_mixing) then

        if ((h_mom_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then

          do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,grid % nVertLevels
                  theta_turb_flux = h_mom_eddy_visc2*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                  flux = 0.5*dvEdge (iEdge) * (rho_edge(k,iEdge)+rho_edge(k-1,iEdge)) * theta_turb_flux
!                  tend_w(k,cell1) = tend_w(k,cell1) + flux
!                  tend_w(k,cell2) = tend_w(k,cell2) - flux
                  tend_w_euler(k,cell1) = tend_w_euler(k,cell1) + flux/areaCell(cell1)
                  tend_w_euler(k,cell2) = tend_w_euler(k,cell2) - flux/areaCell(cell2)
               end do

             end if
           end do
        else if (config_horiz_mixing == "2d_smagorinsky") then

          do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,grid % nVertLevels
!                  theta_turb_flux = h_mom_eddy_visc2*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = 0.25*(kdiff(k,cell1)+kdiff(k,cell2)+kdiff(k-1,cell1)+kdiff(k,cell2))  &
                                        *(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                  flux = 0.5*dvEdge (iEdge) * (rho_edge(k,iEdge)+rho_edge(k-1,iEdge)) * theta_turb_flux
!                  tend_w(k,cell1) = tend_w(k,cell1) + flux
!                  tend_w(k,cell2) = tend_w(k,cell2) - flux
                  tend_w_euler(k,cell1) = tend_w_euler(k,cell1) + flux/areaCell(cell1)
                  tend_w_euler(k,cell2) = tend_w_euler(k,cell2) - flux/areaCell(cell2)
               end do

             end if
           end do
        end if
 
      end if

!ldf (2010-10-10):
!     if ( (h_mom_eddy_visc4 > 0.0) .and. (config_horiz_mixing == "2d_fixed")) then
      if ((h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
          (h_mom_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then

         allocate(delsq_theta(nVertLevels, nCells+1))

         delsq_theta(:,:) = 0.

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            do k=2,grid % nVertLevels
               delsq_theta(k,cell1) = delsq_theta(k,cell1) + dvEdge(iEdge)*0.5*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
               delsq_theta(k,cell2) = delsq_theta(k,cell2) - dvEdge(iEdge)*0.5*(rho_edge(k,iEdge)+rho_edge(k-1,iEdge))*(w(k,cell2) - w(k,cell1))/dcEdge(iEdge)
            end do
         end do

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k=2,nVertLevels
               delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
            end do
         end do

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=2,grid % nVertLevels
                  theta_turb_flux = h_mom_eddy_visc4*(delsq_theta(k,cell2) - delsq_theta(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = theta_turb_flux * meshScalingDel4(iEdge)
                  flux = dvEdge (iEdge) * theta_turb_flux

!                  tend_w(k,cell1) = tend_w(k,cell1) - flux
!                  tend_w(k,cell2) = tend_w(k,cell2) + flux
                  tend_w_euler(k,cell1) = tend_w_euler(k,cell1) - flux/areaCell(cell1)
                  tend_w_euler(k,cell2) = tend_w_euler(k,cell2) + flux/areaCell(cell2)
               end do

            end if
         end do

         deallocate(delsq_theta)

      end if

      end if ! horizontal mixing for w computed in first rk_step

      !
      !  vertical advection, pressure gradient and buoyancy for w
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !

      do iCell = 1, nCells

         wdwz(1) = 0.
         if (config_w_vadv_order == 2) then

            do k=2,nVertLevels
               wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            end do

         else if (config_w_vadv_order == 3) then

            k = 2
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            do k=3,nVertLevels-1
               wdwz(k) = flux3( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)), 1.0_RKIND )
            end do
            k = nVertLevels
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         else if (config_w_vadv_order == 4) then

            k = 2
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))
            do k=3,nVertLevels-1
               wdwz(k) = flux4( w(k-2,iCell),w(k-1,iCell),w(k,iCell),w(k+1,iCell),0.5*(rw(k,iCell)+rw(k-1,iCell)) )
            end do
            k = nVertLevels
            wdwz(k) =  0.25*(rw(k,icell)+rw(k-1,iCell))*(w(k,iCell)+w(k-1,iCell))

         end if

         wdwz(nVertLevels+1) = 0.

         do k=2,nVertLevels

            tend_w(k,iCell) = tend_w(k,iCell)/areaCell(iCell) -rdzu(k)*(wdwz(k+1)-wdwz(k))    &
!SHP-buoy
                                  - cqw(k,iCell)*( rdzu(k)*(pp(k,iCell)-pp(k-1,iCell))        &
                                  + gravity*  &
                                   ( fzm(k)*(rb(k,iCell)*(qtot(k,iCell)) +         &
                                             rr(k,iCell)*(1.+qtot(k,iCell)))                  &
                                    +fzp(k)*(rb(k-1,iCell)*(qtot(k-1,iCell))  +  &
                                             rr(k-1,iCell)*(1.+qtot(k-1,iCell))) ))
         end do
      end do

      !
      !  vertical mixing for w - 2nd order 
      !

      if (rk_step == 1 .or. rk_diffusion) then

      if ( v_mom_eddy_visc2 > 0.0 ) then

         do iCell = 1, grid % nCellsSolve
            do k=2,nVertLevels
               tend_w_euler(k,iCell) = tend_w_euler(k,iCell) + v_mom_eddy_visc2*0.5*(rho_zz(k,iCell)+rho_zz(k-1,iCell))*(  &
                                        (w(k+1,iCell)-w(k  ,iCell))*rdzw(k)                              &
                                       -(w(k  ,iCell)-w(k-1,iCell))*rdzw(k-1) )*rdzu(k)
            end do
         end do

      end if

      end if ! mixing term computed first rk_step


! add in mixing terms for w

      do iCell = 1, grid % nCellsSolve
         do k=2,nVertLevels
            tend_w(k,iCell) = tend_w(k,iCell) + tend_w_euler(k,iCell)
         end do
      end do

      deallocate(qtot)

!----------- rhs for theta

      tend_theta(:,:) = 0.

      !
      !  horizontal advection for theta
      !

      if (config_theta_adv_order == 2) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
               do k=1,grid % nVertLevels
                  flux = dvEdge(iEdge) *  ru(k,iEdge) * ( 0.5*(theta_m(k,cell1) + theta_m(k,cell2)) )
                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
               end do
            end if
         end do

      else if (config_theta_adv_order == 3) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               flux_arr(:) = 0.
               do i=1,nAdvCellsForEdge(iEdge)
                 iCell = advCellsForEdge(i,iEdge)
                 do k=1,grid % nVertLevels
                   scalar_weight = adv_coefs(i,iEdge) + coef_3rd_order*sign(1.0_RKIND,ru(k,iEdge))*adv_coefs_3rd(i,iEdge)
                   flux_arr(k) = flux_arr(k) + scalar_weight* theta_m(k,iCell)
                 end do
               end do

               do k=1,grid % nVertLevels
                  tend_theta(k,cell1) = tend_theta(k,cell1) - ru(k,iEdge)*flux_arr(k)
                  tend_theta(k,cell2) = tend_theta(k,cell2) + ru(k,iEdge)*flux_arr(k)
               end do


!               do k=1,grid % nVertLevels
!
!                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * theta_m(k,cell1)
!                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * theta_m(k,cell2)
!                  do i=1, grid % nEdgesOnCell % array (cell1)
!                     if ( grid % CellsOnCell % array (i,cell1) <= grid%nCells) &
!                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * theta_m(k,grid % CellsOnCell % array (i,cell1))
!                  end do
!                  do i=1, grid % nEdgesOnCell % array (cell2)
!                     if ( grid % CellsOnCell % array (i,cell2) <= grid%nCells) &
!                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * theta_m(k,grid % CellsOnCell % array (i,cell2))
!                  end do
!
!  3rd order stencil
!
!                  if( u(k,iEdge) > 0) then
!                        flux = dvEdge(iEdge) * ru(k,iEdge) * (          &
!                                               0.5*(theta_m(k,cell1) + theta_m(k,cell2))      &
!                                                -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
!                                                -(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
!                     else
!                        flux = dvEdge(iEdge) *  ru(k,iEdge) * (          &
!                                               0.5*(theta_m(k,cell1) + theta_m(k,cell2))      &
!                                                -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12.          &
!                                                +(dcEdge(iEdge) **2) * coef_3rd_order*(d2fdx2_cell1 - d2fdx2_cell2) / 12. )
!                  end if
!
!                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
!                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
!
!               end do

            end if
         end do

      else  if (config_theta_adv_order == 4) then

         do iEdge=1,nEdges
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 <= grid%nCells .and. cell2 <= grid%nCells) then

               do k=1,grid % nVertLevels

                  d2fdx2_cell1 = deriv_two(1,1,iEdge) * theta_m(k,cell1)
                  d2fdx2_cell2 = deriv_two(1,2,iEdge) * theta_m(k,cell2)
                  do i=1, grid % nEdgesOnCell % array (cell1)
                     if ( grid % CellsOnCell % array (i,cell1) <= grid%nCells) &
                     d2fdx2_cell1 = d2fdx2_cell1 + deriv_two(i+1,1,iEdge) * theta_m(k,grid % CellsOnCell % array (i,cell1))
                  end do
                  do i=1, grid % nEdgesOnCell % array (cell2)
                     if ( grid % CellsOnCell % array (i,cell2) <= grid%nCells) &
                     d2fdx2_cell2 = d2fdx2_cell2 + deriv_two(i+1,2,iEdge) * theta_m(k,grid % CellsOnCell % array (i,cell2))
                  end do

                  flux = dvEdge(iEdge) *  ru(k,iEdge) * (                                               &
                                         0.5*(theta_m(k,cell1) + theta_m(k,cell2))                      &
                                          -(dcEdge(iEdge) **2) * (d2fdx2_cell1 + d2fdx2_cell2) / 12. )

                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
               end do

            end if

         end do
      end if

      !
      !  horizontal mixing for theta_m - we could combine this with advection directly (i.e. as a turbulent flux),
      !  but here we can also code in hyperdiffusion if we wish (2nd order at present)
      !

      if (rk_step == 1 .or. rk_diffusion) then

      tend_theta_euler = 0.

      if (delsq_horiz_mixing) then
         if ( (h_theta_eddy_visc2 > 0.0) .and. (config_horiz_mixing == "2d_fixed") ) then

            do iEdge=1,grid % nEdges
               cell1 = grid % cellsOnEdge % array(1,iEdge)
               cell2 = grid % cellsOnEdge % array(2,iEdge)
               if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
  
                  do k=1,grid % nVertLevels
                     theta_turb_flux = h_theta_eddy_visc2*prandtl_inv*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
                     theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                     flux = dvEdge (iEdge) * rho_edge(k,iEdge) * theta_turb_flux
!                     tend_theta(k,cell1) = tend_theta(k,cell1) + flux
!                     tend_theta(k,cell2) = tend_theta(k,cell2) - flux
                     tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) + flux/areaCell(cell1)
                     tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) - flux/areaCell(cell2)
                  end do
  
               end if
            end do

         else if (  ( config_horiz_mixing == "2d_smagorinsky") ) then

            do iEdge=1,grid % nEdges
               cell1 = grid % cellsOnEdge % array(1,iEdge)
               cell2 = grid % cellsOnEdge % array(2,iEdge)
               if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then
 
                  do k=1,grid % nVertLevels
                     theta_turb_flux = 0.5*(kdiff(k,cell1)+kdiff(k,cell2))*prandtl_inv  &
                                           *(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
                     theta_turb_flux = theta_turb_flux * meshScalingDel2(iEdge)
                     flux = dvEdge (iEdge) * rho_edge(k,iEdge) * theta_turb_flux
!                     tend_theta(k,cell1) = tend_theta(k,cell1) + flux
!                     tend_theta(k,cell2) = tend_theta(k,cell2) - flux
                     tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) + flux/areaCell(cell1)
                     tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) - flux/areaCell(cell2)
                  end do
  
               end if
            end do
         end if

      end if

!ldf (2010-10-10):
!     if ( (h_theta_eddy_visc4 > 0.0) .and. (config_horiz_mixing == "2d_fixed") ) then
      if ((h_theta_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_fixed") .or. &
          (h_theta_eddy_visc4 > 0.0 .and. config_horiz_mixing == "2d_smagorinsky")) then

         allocate(delsq_theta(nVertLevels, nCells+1))

         delsq_theta(:,:) = 0.

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            do k=1,grid % nVertLevels
               delsq_theta(k,cell1) = delsq_theta(k,cell1) + dvEdge(iEdge)*rho_edge(k,iEdge)*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
               delsq_theta(k,cell2) = delsq_theta(k,cell2) - dvEdge(iEdge)*rho_edge(k,iEdge)*(theta_m(k,cell2) - theta_m(k,cell1))/dcEdge(iEdge)
            end do
         end do

         do iCell = 1, nCells
            r = 1.0 / areaCell(iCell)
            do k=1,nVertLevels
               delsq_theta(k,iCell) = delsq_theta(k,iCell) * r
            end do
         end do

         do iEdge=1,grid % nEdges
            cell1 = grid % cellsOnEdge % array(1,iEdge)
            cell2 = grid % cellsOnEdge % array(2,iEdge)
            if (cell1 <= nCellsSolve .or. cell2 <= nCellsSolve) then

               do k=1,grid % nVertLevels
                  theta_turb_flux = h_theta_eddy_visc4*prandtl_inv*(delsq_theta(k,cell2) - delsq_theta(k,cell1))/dcEdge(iEdge)
                  theta_turb_flux = theta_turb_flux * meshScalingDel4(iEdge)
                  flux = dvEdge (iEdge) * theta_turb_flux

!                  tend_theta(k,cell1) = tend_theta(k,cell1) - flux
!                  tend_theta(k,cell2) = tend_theta(k,cell2) + flux
                  tend_theta_euler(k,cell1) = tend_theta_euler(k,cell1) - flux/areaCell(cell1)
                  tend_theta_euler(k,cell2) = tend_theta_euler(k,cell2) + flux/areaCell(cell2)
               end do

            end if
         end do

         deallocate(delsq_theta)

      end if

      end if ! theta mixing calculated first rk_step

      !
      !  vertical advection plus diabatic term
      !  Note: we are also dividing through by the cell area after the horizontal flux divergence
      !
      do iCell = 1, nCells

         wdtz(1) = 0.
         if (config_theta_vadv_order == 2) then

            do k=2,nVertLevels
               wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            end do

         else if (config_theta_vadv_order == 3) then

            k = 2
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            do k=3,nVertLevels-1
               wdtz(k) = flux3( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell), coef_3rd_order )
            end do
            k = nVertLevels
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))

         else if (config_theta_vadv_order == 4) then

            k = 2
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))
            do k=3,nVertLevels-1
               wdtz(k) = flux4( theta_m(k-2,iCell),theta_m(k-1,iCell),theta_m(k,iCell),theta_m(k+1,iCell), rw(k,iCell) )
            end do
            k = nVertLevels
            wdtz(k) =  rw(k,icell)*(fzm(k)*theta_m(k,iCell)+fzp(k)*theta_m(k-1,iCell))

         end if


         wdtz(nVertLevels+1) = 0.

         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell)/areaCell(iCell) -rdzw(k)*(wdtz(k+1)-wdtz(k))
            tend_theta(k,iCell) = tend_theta(k,iCell) + rho_zz(k,iCell)*rt_diabatic_tend(k,iCell)
         end do
      end do

      !
      !  vertical mixing for theta - 2nd order 
      !

      if (rk_step == 1 .or. rk_diffusion) then

      if ( v_theta_eddy_visc2 > 0.0 ) then

         if (config_mix_full) then

         do iCell = 1, grid % nCellsSolve
            do k=2,nVertLevels-1
               z1 = zgrid(k-1,iCell)
               z2 = zgrid(k  ,iCell)
               z3 = zgrid(k+1,iCell)
               z4 = zgrid(k+2,iCell)

               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)

!               tend_theta(k,iCell) = tend_theta(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
!                                        (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
!                                       -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
               tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                        (theta_m(k+1,iCell)-theta_m(k  ,iCell))/(zp-z0)                 &
                                       -(theta_m(k  ,iCell)-theta_m(k-1,iCell))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

         else  ! idealized cases where we mix on the perturbation from the initial 1-D state

         do iCell = 1, grid % nCellsSolve
            do k=2,nVertLevels-1
               z1 = zgrid(k-1,iCell)
               z2 = zgrid(k  ,iCell)
               z3 = zgrid(k+1,iCell)
               z4 = zgrid(k+2,iCell)

               zm = 0.5*(z1+z2)
               z0 = 0.5*(z2+z3)
               zp = 0.5*(z3+z4)

!               tend_theta(k,iCell) = tend_theta(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
!                                        ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
!                                       -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
               tend_theta_euler(k,iCell) = tend_theta_euler(k,iCell) + v_theta_eddy_visc2*prandtl_inv*rho_zz(k,iCell)*(&
                                        ((theta_m(k+1,iCell)-t_init(k+1,iCell))-(theta_m(k  ,iCell)-t_init(k,iCell)))/(zp-z0)      &
                                       -((theta_m(k  ,iCell)-t_init(k,iCell))-(theta_m(k-1,iCell)-t_init(k-1,iCell)))/(z0-zm) )/(0.5*(zp-zm))
            end do
         end do

         end if

      end if

      end if ! compute theta mixing on first rk_step

      do iCell = 1, grid % nCellsSolve
         do k=1,nVertLevels
            tend_theta(k,iCell) = tend_theta(k,iCell) + tend_theta_euler(k,iCell)
         end do
      end do

   end subroutine atm_compute_dyn_tend

!-------

   subroutine atm_compute_solve_diagnostics(dt, s, diag, grid)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! Compute diagnostic fields used in the tendency computations
   !
   ! Input: grid - grid metadata
   !
   ! Output: s - computed diagnostics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      implicit none

      real (kind=RKIND), intent(in) :: dt
      type (state_type), intent(inout) :: s
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(in) :: grid


      integer :: iEdge, iCell, iVertex, k, cell1, cell2, eoe, i
      real (kind=RKIND) :: h_vertex, r

      integer :: nCells, nEdges, nVertices, nVertLevels
      real (kind=RKIND), dimension(:), pointer :: fVertex, fEdge, dvEdge, dcEdge, areaCell, areaTriangle
      real (kind=RKIND), dimension(:,:), pointer :: vh, weightsOnEdge, kiteAreasOnVertex, h_edge, h, u, v, &
                                                    circulation, vorticity, ke, pv_edge, pv_vertex, pv_cell, gradPVn, gradPVt, &
                                                    divergence
      integer, dimension(:,:), pointer :: cellsOnEdge, cellsOnVertex, verticesOnEdge, edgesOnCell, edgesOnEdge, edgesOnVertex
      integer, dimension(:), pointer :: nEdgesOnCell, nEdgesOnEdge

      !WCS-instability
      logical, parameter :: hollingsworth=.true.
      real (kind=RKIND), allocatable, dimension(:,:) :: ke_vertex
      real (kind=RKIND)  :: ke_fact


      h           => s % rho_zz % array
      u           => s % u % array
      v           => diag % v % array
      vh          => diag % rv % array
      h_edge      => diag % rho_edge % array
      circulation => diag % circulation % array
      vorticity   => diag % vorticity % array
      divergence  => diag % divergence % array
      ke          => diag % ke % array
      pv_edge     => diag % pv_edge % array
      pv_vertex   => diag % pv_vertex % array
      pv_cell     => diag % pv_cell % array
      gradPVn     => diag % gradPVn % array
      gradPVt     => diag % gradPVt % array

      weightsOnEdge     => grid % weightsOnEdge % array
      kiteAreasOnVertex => grid % kiteAreasOnVertex % array
      cellsOnEdge       => grid % cellsOnEdge % array
      cellsOnVertex     => grid % cellsOnVertex % array
      verticesOnEdge    => grid % verticesOnEdge % array
      nEdgesOnCell      => grid % nEdgesOnCell % array
      edgesOnCell       => grid % edgesOnCell % array
      nEdgesOnEdge      => grid % nEdgesOnEdge % array
      edgesOnEdge       => grid % edgesOnEdge % array
      edgesOnVertex     => grid % edgesOnVertex % array
      dcEdge            => grid % dcEdge % array
      dvEdge            => grid % dvEdge % array
      areaCell          => grid % areaCell % array
      areaTriangle      => grid % areaTriangle % array
      fVertex           => grid % fVertex % array
      fEdge             => grid % fEdge % array
                  
      nCells      = grid % nCells
      nEdges      = grid % nEdges
      nVertices   = grid % nVertices
      nVertLevels = grid % nVertLevels

      !
      ! Compute height on cell edges at velocity locations
      !
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
            h_edge(k,iEdge) = 0.5 * (h(k,cell1) + h(k,cell2))
         end do
      end do


      !
      ! Compute circulation and relative vorticity at each vertex
      !
      circulation(:,:) = 0.0
      do iEdge=1,nEdges
         do k=1,nVertLevels
            circulation(k,verticesOnEdge(1,iEdge)) = circulation(k,verticesOnEdge(1,iEdge)) - dcEdge(iEdge) * u(k,iEdge)
            circulation(k,verticesOnEdge(2,iEdge)) = circulation(k,verticesOnEdge(2,iEdge)) + dcEdge(iEdge) * u(k,iEdge)
         end do
      end do
      do iVertex=1,nVertices
         do k=1,nVertLevels
            vorticity(k,iVertex) = circulation(k,iVertex) / areaTriangle(iVertex)
         end do
      end do


      !
      ! Compute the divergence at each cell center
      !
      divergence(:,:) = 0.0
      do iEdge=1,nEdges
         cell1 = cellsOnEdge(1,iEdge)
         cell2 = cellsOnEdge(2,iEdge)
         do k=1,nVertLevels
           divergence(k,cell1) = divergence(k,cell1) + u(k,iEdge)*dvEdge(iEdge)
           divergence(k,cell2) = divergence(k,cell2) - u(k,iEdge)*dvEdge(iEdge)
         end do
      end do
      do iCell = 1,nCells
         r = 1.0 / areaCell(iCell)
         do k = 1,nVertLevels
            divergence(k,iCell) = divergence(k,iCell) * r
         end do
      end do


      !
      ! Compute kinetic energy in each cell
      !
      ke(:,:) = 0.0
      do iCell=1,nCells
         do i=1,nEdgesOnCell(iCell)
            iEdge = edgesOnCell(i,iCell)
            do k=1,nVertLevels
               ke(k,iCell) = ke(k,iCell) + 0.25 * dcEdge(iEdge) * dvEdge(iEdge) * u(k,iEdge)**2.0
            end do
         end do
         do k=1,nVertLevels
            ke(k,iCell) = ke(k,iCell) / areaCell(iCell)
         end do
      end do

      !WCS-instability
      ! Compute ke at cell vertices - AG's new KE construction, part 1
      ! *** approximation here because we don't have inner triangle areas
      !
      if (hollingsworth) then
      allocate (ke_vertex(nVertLevels,nVertices))
      do iVertex=1,nVertices
         do k=1,nVertLevels
!            ke_vertex(k,iVertex) = (  subTriangleAreasOnVertex(1,iVertex)*u(k,EdgesOnVertex(1,iVertex))**2.0  &
!                                    + subTriangleAreasOnVertex(2,iVertex)*u(k,EdgesOnVertex(2,iVertex))**2.0  &
!                                    + subTriangleAreasOnVertex(3,iVertex)*u(k,EdgesOnVertex(3,iVertex))**2.0  &
!                                                 ) / AreaTriangle(iVertex)

            ke_vertex(k,iVertex) = (  dcEdge(EdgesOnVertex(1,iVertex))*dvEdge(EdgesOnVertex(1,iVertex))*u(k,EdgesOnVertex(1,iVertex))**2.0  &
                                     +dcEdge(EdgesOnVertex(2,iVertex))*dvEdge(EdgesOnVertex(2,iVertex))*u(k,EdgesOnVertex(2,iVertex))**2.0  &
                                     +dcEdge(EdgesOnVertex(3,iVertex))*dvEdge(EdgesOnVertex(3,iVertex))*u(k,EdgesOnVertex(3,iVertex))**2.0  &
                                                 ) * 0.25 / AreaTriangle(iVertex)

         end do
      end do

      ! adjust ke at cell vertices - AG's new KE construction, part 2
      !

      ke_fact = 1.0 - .375

      do iCell=1,nCells
         do k=1,nVertLevels
            ke(k,iCell) = ke_fact*ke(k,iCell)
         end do
      end do

      do iVertex = 1, nVertices
       do i=1,grid % vertexDegree
          iCell = cellsOnVertex(i,iVertex)
          do k = 1,nVertLevels
             ke(k,iCell) = ke(k,iCell) + (1.-ke_fact)*kiteAreasOnVertex(i, iVertex) * ke_vertex(k, iVertex) / areaCell(iCell)
          end do
       end do
      end do
      deallocate (ke_vertex)
      end if
      !END of WCS-instability

      !
      ! Compute v (tangential) velocities
      !
      v(:,:) = 0.0
      do iEdge = 1,nEdges
         do i=1,nEdgesOnEdge(iEdge)
            eoe = edgesOnEdge(i,iEdge)
            do k = 1,nVertLevels
               v(k,iEdge) = v(k,iEdge) + weightsOnEdge(i,iEdge) * u(k, eoe)
            end do
         end do
      end do


      !
      ! Compute height at vertices, pv at vertices, and average pv to edge locations
      !  ( this computes pv_vertex at all vertices bounding real cells )
      !
      do iVertex = 1,nVertices
         do k=1,nVertLevels
            h_vertex = 0.0
            do i=1,grid % vertexDegree
               h_vertex = h_vertex + h(k,cellsOnVertex(i,iVertex)) * kiteAreasOnVertex(i,iVertex)
            end do
            h_vertex = h_vertex / areaTriangle(iVertex)

            pv_vertex(k,iVertex) = (fVertex(iVertex) + vorticity(k,iVertex)) / h_vertex
         end do
      end do


      !
      ! Compute pv at the edges
      !   ( this computes pv_edge at all edges bounding real cells )
      !
      pv_edge(:,:) = 0.0
      do iVertex = 1,nVertices
         do i=1,grid % vertexDegree
            iEdge = edgesOnVertex(i,iVertex)
            do k=1,nVertLevels
               pv_edge(k,iEdge) =  pv_edge(k,iEdge)  + 0.5 * pv_vertex(k,iVertex)
            end do
         end do
      end do

      !
      ! Compute pv at cell centers
      !    ( this computes pv_cell for all real cells )
      !
      pv_cell(:,:) = 0.0
      do iVertex = 1, nVertices
         do i=1,grid % vertexDegree
            iCell = cellsOnVertex(i,iVertex)
            do k = 1,nVertLevels
               pv_cell(k,iCell) = pv_cell(k,iCell) + kiteAreasOnVertex(i, iVertex) * pv_vertex(k, iVertex) / areaCell(iCell)
            end do
         end do
      end do


      if (config_apvm_upwinding > 0.0) then

      !
      ! Modify PV edge with upstream bias. 
      !
      ! Compute gradient of PV in the tangent direction
      !   ( this computes gradPVt at all edges bounding real cells )
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
            gradPVt(k,iEdge) = (pv_vertex(k,verticesOnEdge(2,iEdge)) - pv_vertex(k,verticesOnEdge(1,iEdge))) / &
                               dvEdge(iEdge)
         end do
      end do

      !
      ! Compute gradient of PV in normal direction
      !   (tdr: 2009-10-02: this is not correct because the pv_cell in the halo is not correct)
      !
      gradPVn(:,:) = 0.0
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
            gradPVn(k,iEdge) = (pv_cell(k,cellsOnEdge(2,iEdge)) - pv_cell(k,cellsOnEdge(1,iEdge))) / &
                                 dcEdge(iEdge)
         end do
      end do

      do iEdge = 1,nEdges
         do k = 1,nVertLevels
            pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * v(k,iEdge) * dt * gradPVt(k,iEdge)
         end do
      end do

      ! Modify PV edge with upstream bias.
      !
      do iEdge = 1,nEdges
         do k = 1,nVertLevels
            pv_edge(k,iEdge) = pv_edge(k,iEdge) - config_apvm_upwinding * u(k,iEdge) *dt * gradPVn(k,iEdge)
         end do
      end do

      end if  ! apvm upwinding


   end subroutine atm_compute_solve_diagnostics

!----------

   subroutine atm_init_coupled_diagnostics( state, diag, grid )

      implicit none
   
      type (state_type), intent(inout) :: state
      type (diag_type), intent(inout) :: diag
      type (mesh_type), intent(inout) :: grid

      !SHP-w
      integer :: k,iCell,iEdge,iCell1,iCell2, cell1, cell2, coef_3rd_order
      real (kind=RKIND) :: p0, rcv, flux

      !SHP-w
      coef_3rd_order = config_coef_3rd_order
      if(config_theta_adv_order /=3) coef_3rd_order = 0

      rcv = rgas / (cp-rgas)
      p0 = 1.e5  ! this should come from somewhere else...

      do iCell=1,grid%nCells
         do k=1,grid%nVertLevels
            state % theta_m % array(k,iCell) = diag % theta % array(k,iCell) * (1._RKIND + rvord * state % scalars % array(state % index_qv,k,iCell))
            state % rho_zz % array(k,iCell) = diag % rho % array(k,iCell) / grid % zz % array(k,iCell)
         end do
      end do

      do iEdge = 1, grid % nEdges
         iCell1 = grid % cellsOnEdge % array(1,iEdge)
         iCell2 = grid % cellsOnEdge % array(2,iEdge)
         do k=1,grid % nVertLevels
            diag % ru % array(k,iEdge) = 0.5 * state % u % array(k,iEdge) * (state % rho_zz % array(k,iCell1) + state % rho_zz % array(k,iCell2))
         end do
      end do

!      ! Compute w from rho_zz and rw
!      do iCell=1,grid%nCells
!         diag % rw % array(1,iCell) = 0.
!         diag % rw % array(grid%nVertLevels+1,iCell) = 0.
!         do k=2,grid%nVertLevels
!            diag % rw % array(k,iCell) = state % w % array(k,iCell)     &
!                          * (grid % fzp % array(k) * state % rho_zz % array(k-1,iCell) + grid % fzm % array(k) * state % rho_zz % array(k,iCell))
!         end do
!      end do


! WCS bug fix 20110916

      ! Compute rw (i.e. rho_zz * omega) from rho_zz, w, and ru.
      ! We are reversing the procedure we use in subroutine atm_recover_large_step_variables.
      ! first, the piece that depends on w.
      do iCell=1,grid%nCells
         diag % rw % array(1,iCell) = 0.
         diag % rw % array(grid%nVertLevels+1,iCell) = 0.
         do k=2,grid%nVertLevels
            diag % rw % array(k,iCell) = state % w % array(k,iCell)     &
                          * (grid % fzp % array(k) * state % rho_zz % array(k-1,iCell) + grid % fzm % array(k) * state % rho_zz % array(k,iCell)) &
                          * (grid % fzp % array(k) * grid % zz % array(k-1,iCell) + grid % fzm % array(k) * grid % zz % array(k,iCell))
         end do
      end do
  
      !SHP-w
      ! next, the piece that depends on ru
      do iEdge=1,grid%nEdges
        cell1 = grid % CellsOnEdge % array(1,iEdge)
        cell2 = grid % CellsOnEdge % array(2,iEdge)
          do k = 2, grid % nVertLevels
            flux = (grid % fzm % array(k) * diag % ru % array(k,iEdge)+grid % fzp % array(k) * diag % ru % array(k-1,iEdge))
            diag % rw % array(k,cell2) = diag % rw % array(k,cell2)   &
                          + (grid % zb % array(k,2,iEdge) + coef_3rd_order * sign(1.0_RKIND,flux) * grid % zb3 % array(k,2,iEdge))*flux   &
                          * (grid % fzp % array(k) * grid % zz % array(k-1,cell2) + grid % fzm % array(k) * grid % zz % array(k,cell2))
            diag % rw % array(k,cell1) = diag % rw % array(k,cell1)   &
                          - (grid % zb % array(k,1,iEdge) + coef_3rd_order * sign(1.0_RKIND,flux) * grid % zb3 % array(k,1,iEdge))*flux   &
                          * (grid % fzp % array(k) * grid % zz % array(k-1,cell1) + grid % fzm % array(k) * grid % zz % array(k,cell1))
          end do
      end do

!  end WCS bug fix

      do iCell = 1, grid % nCells
         do k=1,grid % nVertLevels
            diag % rho_p % array(k,iCell) = state % rho_zz % array(k,iCell) - diag % rho_base % array(k,iCell)
         end do
      end do

      do iCell = 1, grid % nCells
         do k=1,grid % nVertLevels
            diag % rtheta_base % array(k,iCell) = diag % theta_base % array(k,iCell) * diag % rho_base % array(k,iCell)
         end do
      end do

      do iCell = 1, grid % nCells
         do k=1,grid % nVertLevels
            diag % rtheta_p % array(k,iCell) = state % theta_m % array(k,iCell) * diag % rho_p % array(k,iCell)  &
                                             + diag % rho_base % array(k,iCell)   * (state % theta_m % array(k,iCell) - diag % theta_base % array(k,iCell))
         end do
      end do

      do iCell=1,grid % nCells
         do k=1,grid % nVertLevels
            diag % exner % array(k,iCell) = (grid % zz % array(k,iCell) * (rgas/p0) * (diag % rtheta_p % array(k,iCell) + diag % rtheta_base % array(k,iCell)))**rcv
         end do
      end do

      do iCell=1,grid % nCells
         do k=1,grid % nVertLevels
            diag % pressure_p % array(k,iCell) = grid % zz % array(k,iCell) * rgas &
                                               * (  diag % exner % array(k,iCell) * diag % rtheta_p % array(k,iCell) &
                                                  + diag % rtheta_base % array(k,iCell) * (diag % exner % array(k,iCell) - diag % exner_base % array(k,iCell)) &
                                                 )
         end do
      end do

   end subroutine atm_init_coupled_diagnostics


end module atm_time_integration
