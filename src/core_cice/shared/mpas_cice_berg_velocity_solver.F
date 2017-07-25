!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_velocity_solver
!
!> \brief   Calculate iceberg velocities
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This module calculates iceberg velocities by solving the
!> iceberg momentum equation,
!>  m * du/dt = F_a + F_o + F_i + F_c + F_t,
!> where m is iceberg mass, u is iceberg velocity, and the forcing terms
!> are atmosphere, ocean, sea ice, Coriolis, and sea surface tilt.
!
!-----------------------------------------------------------------------

module cice_berg_velocity_solver

  use mpas_derived_types
  use mpas_pool_routines
  use mpas_timekeeping
  use mpas_dmpar
  use mpas_timer
  use mpas_log, only: mpas_log_write

  implicit none

  private
  save

  public :: &
       cice_init_berg_velocity_solver, &
       cice_run_berg_velocity_solver

  ! berg velocity solver constants
  real(kind=RKIND), parameter, private :: &
       bergMassMinimum = 0.01_RKIND, &
       strengthThresh = 1.0e4_RKIND

contains

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_init_berg_velocity_solver
!
!> \brief   Initialize berg velocity solver
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine initializes iceberg velocities as 0 and sets
!> the iceberg subcycling time step.
!
!-----------------------------------------------------------------------

  subroutine cice_init_berg_velocity_solver(domain)

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergVelocitySolverPool

    real(kind=RKIND), pointer :: &
         bergDynamicsTimeStep, &
         config_dt

    integer, pointer :: &
         config_berg_dynamics_subcycle_number

    ! set timesteps and initialize velocities to 0
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_config(block % configs, "config_dt", config_dt)
       call MPAS_pool_get_config(block % configs, "config_berg_dynamics_subcycle_number", config_berg_dynamics_subcycle_number)

       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergDynamicsTimeStep", bergDynamicsTimeStep)

       bergDynamicsTimeStep = config_dt / real(config_berg_dynamics_subcycle_number,RKIND)

       call init_berg_velocity(block)

       block => block % next
    enddo

  end subroutine cice_init_berg_velocity_solver

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_velocity
!
!> \brief   Initialize icebergs velocity
!> \author  Darin Comeau, LANL
!> \date    6 June 2017
!> \details This routine is intended to initialize iceberg velocities.
!
!-----------------------------------------------------------------------

  subroutine init_berg_velocity(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity

    integer :: &
         iVertex, &
         iCategory

    call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
    call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)

    call MPAS_pool_get_dimension(meshPool, "nVerticesSolve", nVerticesSolve)
    call MPAS_pool_get_dimension(meshPool, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
    call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)

    do iVertex = 1, nVerticesSolve

       do iCategory = 1, nBergCategories

          uBergVelocity(iCategory,iVertex) = 0.0_RKIND
          vBergVelocity(iCategory,iVertex) = 0.0_RKIND

       enddo

    enddo

  end subroutine init_berg_velocity

!-----------------------------------------------------------------------
! Time stepping
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_run_berg_velocity_solver
!
!> \brief   Calculate iceberg velocities by solving momentum equation
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine is a wrapper for the iceberg velocity solver
!> over a sea ice time step.
!
!-----------------------------------------------------------------------

  subroutine cice_run_berg_velocity_solver(domain, clock)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    logical, pointer :: &
         config_use_berg_velocity_solver

    ! determine if velocity solver switched on
    call MPAS_pool_get_config(domain % configs, "config_use_berg_velocity_solver", config_use_berg_velocity_solver)

    if (config_use_berg_velocity_solver) then

       ! pre subcycle
       call mpas_timer_start("Berg velocity solver pre-cycle")
       call berg_velocity_solver_pre_subcycle(domain)
       call mpas_timer_stop("Berg velocity solver pre-cycle")

       ! subcycle the dynamics
       call mpas_timer_start("Berg velocity solver sub-cycle")
       call berg_velocity_solver_subcycle(domain, clock)
       call mpas_timer_stop("Berg velocity solver sub-cycle")

    endif ! config_use_berg_velocity_solver

  end subroutine cice_run_berg_velocity_solver!}}}

!-----------------------------------------------------------------------
! Pre sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_velocity_solver_pre_subcycle
!
!> \brief   Prepare for iceberg velocity solver subcycling.
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine initializes and prepares variables & masks
!> needed for iceberg velocity subcycling.
!
!-----------------------------------------------------------------------

  subroutine berg_velocity_solver_pre_subcycle(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! initialize subcycle variables
    call mpas_timer_start("init subcycle var")
    call init_berg_subcycle_variables(domain)
    call mpas_timer_stop("init subcycle var")

    ! aggregate categories for area and volume into total mass
    call mpas_timer_start("agg berg mass and area")
    call aggregate_berg_mass(domain)
    call mpas_timer_stop("agg berg mass and area")

    ! calculate velocity calculation masks
    call mpas_timer_start("calc berg masks")
    call berg_calculation_masks(domain)
    call mpas_timer_stop("calc berg masks")

    ! calculate horizontal and vertical surface areas
    call mpas_timer_start("berg vertical area")
    call berg_vertical_areas(domain)
    call mpas_timer_stop("berg vertical area")

  end subroutine berg_velocity_solver_pre_subcycle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  init_berg_subcycle_variables
!
!> \brief   Initialize berg subcycle variables
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine initializes the forcing terms and iceberg
!> velocities before subcycling.
!
!-----------------------------------------------------------------------

  subroutine init_berg_subcycle_variables(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergVelocitySolverPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         uBergVelocityInitial, &
         vBergVelocityInitial, &
         bergDraftVertex, &
         bergMassVertex, &
         oneBergMassVertex, &
         oneBergAreaVertex, &
         bergVertAirVertex, &
         bergVertOcnVertex, &
         bergVertIceVertex, &
         bergCoefAtmU, &
         bergCoefAtmV, &
         bergCoefOcnU, &
         bergCoefOcnV, &
         bergCoefIceU, &
         bergCoefIceV, &
         bergCoefCorU, &
         bergCoefCorV, &
         bergCoefTltU, &
         bergCoefTltV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocityVertex, &
         vAirVelocityVertex, &
         airDensityVertex

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    integer :: &
         iVertex, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocityInitial", uBergVelocityInitial)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocityInitial", vBergVelocityInitial)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergDraftVertex", bergDraftVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergMassVertex", bergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergMassVertex", oneBergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergAreaVertex", oneBergAreaVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertAirVertex", bergVertAirVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertOcnVertex", bergVertOcnVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertIceVertex", bergVertIceVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmU", bergCoefAtmU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmV", bergCoefAtmV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnU", bergCoefOcnU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnV", bergCoefOcnV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceU", bergCoefIceU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceV", bergCoefIceV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorU", bergCoefCorU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorV", bergCoefCorV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltU", bergCoefTltU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltV", bergCoefTltV)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             uBergVelocityInitial(iCategory,iVertex) = uBergVelocity(iCategory,iVertex)
             vBergVelocityInitial(iCategory,iVertex) = vBergVelocity(iCategory,iVertex)
             bergVelocityMask(iCategory,iVertex) = 0
             bergDraftVertex(iCategory,iVertex) = 0.0_RKIND
             bergMassVertex(iCategory,iVertex)  = 0.0_RKIND
             oneBergMassVertex(iCategory,iVertex) = 0.0_RKIND
             oneBergAreaVertex(iCategory,iVertex) = 0.0_RKIND
             bergVertAirVertex(iCategory,iVertex) = 0.0_RKIND
             bergVertOcnVertex(iCategory,iVertex) = 0.0_RKIND
             bergVertIceVertex(iCategory,iVertex) = 0.0_RKIND
             bergCoefAtmU(iCategory,iVertex) = 0.0_RKIND
             bergCoefAtmV(iCategory,iVertex) = 0.0_RKIND
             bergCoefOcnU(iCategory,iVertex) = 0.0_RKIND
             bergCoefOcnV(iCategory,iVertex) = 0.0_RKIND
             bergCoefIceU(iCategory,iVertex) = 0.0_RKIND
             bergCoefIceV(iCategory,iVertex) = 0.0_RKIND
             bergCoefCorU(iCategory,iVertex) = 0.0_RKIND
             bergCoefCorV(iCategory,iVertex) = 0.0_RKIND
             bergCoefTltU(iCategory,iVertex) = 0.0_RKIND
             bergCoefTltV(iCategory,iVertex) = 0.0_RKIND

          enddo ! iCategory

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine init_berg_subcycle_variables!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  aggregate_berg_mass
!
!> \brief   Aggregate berg masses
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine is intended to update iceberg masses for the
!> velocity solver.
!
!-----------------------------------------------------------------------

  subroutine aggregate_berg_mass(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         tracersPool, &
         tracersAggregatePool

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)

       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

       call MPAS_pool_get_array(tracersPool, "bergMassCategory", bergMassCategory, 1)

       call MPAS_pool_get_array(tracersAggregatePool, "bergMassCell", bergMassCell)

       do iCell = 1, nCellsSolve

          bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

       enddo ! iCell

       block => block % next
    enddo

    ! halo exchange berg mass
    call mpas_timer_start("berg mass halo")
    call MPAS_dmpar_field_halo_exch(domain, 'bergMassCell')
    call mpas_timer_stop("berg mass halo")

  end subroutine aggregate_berg_mass!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_calculation_masks
!
!> \brief   Berg calculation masks
!> \author  Darin Comeau, LANL
!> \date    July 5, 2017
!> \details This subroutine creates masks needed for the velocity solver,
!> which at present is only one computational mask.
!
!-----------------------------------------------------------------------

  subroutine berg_calculation_masks(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         tracersPool, &
         bergVelocitySolverPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassVertex

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    integer, pointer :: &
         nBergCategories

    integer :: &
         iCategory

    ! interpolate area and mass from cells to vertices
    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)

       call MPAS_pool_get_array(tracersPool, "bergMassCategory", bergMassCategory, 1)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergMassVertex", bergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)

       do iCategory = 1, nBergCategories

          call cice_interpolate_cell_to_vertex(&
               meshPool, &
               bergMassVertex(iCategory,:), &
               bergMassCategory(iCategory,:))

       enddo

       block => block % next
    end do

    ! calculate computational masks
    call berg_velocity_calculation_mask(domain)

    ! halo exchange velocity mask
    call mpas_timer_start("velocity mask halo")
    call MPAS_dmpar_field_halo_exch(domain, 'bergVelocityMask')
    call mpas_timer_stop("velocity mask halo")

  end subroutine berg_calculation_masks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_velocity_calculation_mask
!
!> \brief   Create mask for berg velocity calculation
!> \author  Darin Comeau, LANL
!> \date    July 17, 2017
!> \details This subroutine creates a mask for calculating iceberg
!> velocities and indicates when the iceberg is captured in the sea ice
!> pack.
!
!-----------------------------------------------------------------------

  subroutine berg_velocity_calculation_mask(domain)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool, &
         boundaryPool, &
         oceanCouplingPool, &
         tracersAggregatePool, &
         velocitySolverPool

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassVertex

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         areacell, &
         icePressure

    integer, dimension(:), pointer :: &
         interiorVertex, &
         landIceMaskVertex

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    integer :: &
         iVertex, &
         iCategory, &
         iVertexDegree, &
         iCell

    integer, pointer :: &
         nVerticesSolve, &
         nVertices, &
         nBergCategories, &
         vertexDegree

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nVertices", nVertices)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)
       call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "boundary", boundaryPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergMassVertex", bergMassVertex)

       call MPAS_pool_get_array(boundaryPool, "interiorVertex", interiorVertex)

       call MPAS_pool_get_array(oceanCouplingPool, "landIceMaskVertex", landIceMaskVertex)

       call MPAS_pool_get_array(velocitySolverPool, "icePressure", icePressure)
       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             bergVelocityMask(iCategory,iVertex) = 0

             if (interiorVertex(iVertex) == 1 .and. &
                landIceMaskVertex(iVertex) == 0 .and. &
                bergMassVertex(iCategory,iVertex) > bergMassMinimum) then

                ! this vertex has sufficient ice
                bergVelocityMask(iCategory,iVertex) = 1

                do iVertexDegree = 1, vertexDegree

                   iCell = cellsOnVertex(iVertexDegree,iVertex)

                   ! this berg is captured by sea ice
                   if (iceAreaCell(iCell)/areaCell(iCell) > 0.9_RKIND .and. icePressure(iCell) > strengthThresh) &
                      bergVelocityMask(iCategory,iVertex) = 2

                enddo

             endif

          enddo

       enddo ! iVertex

       do iVertex = nVerticesSolve+1, nVertices

          bergVelocityMask(:,iVertex) = 0

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine berg_velocity_calculation_mask!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_vertical_areas
!
!> \brief   Calculate berg vertical areas
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details Thi subroutine calculates the vertical areas of the iceberg
!> in contact with air, ocean, and sea ice.
!
!-----------------------------------------------------------------------

  subroutine berg_vertical_areas(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_constants, only: &
         cicePuny, &
         ciceDensityIce, &
         ciceDensitySnow, &
         ciceDensitySeaWater, &
         pii

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         tracersPool, &
         tracersAggregatePool, &
         bergStatePool, &
         bergVelocityPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, &
         bergLength, &
         bergHeight, &
         bergDraft, &
         oneBergArea, &
         oneBergMass, &
         bergVertAir, &
         bergVertOcn, &
         bergVertIce

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergDraftVertex, &
         oneBergAreaVertex, &
         oneBergMassVertex, &
         bergVertAirVertex, &
         bergVertOcnVertex, &
         bergVertIceVertex

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell, &
         iceAreaCell, &
         iceVolumeCell, &
         snowVolumeCell

    character(len=strKIND), pointer :: &
         config_berg_horizontal_shape

    real(kind=RKIND), pointer :: &
         bergDensity

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    integer :: &
         iCell, &
         iCategory

    real(kind=RKIND) :: &
         bergSnow, &
         seaIceThick, &
         seaIceSnow, &
         seaIceDraft, &
         zAir, &
         zOcn, &
         zIce, &
         tab2cyl, &
         p66, &
         p83

    call MPAS_pool_get_config(domain % configs, "config_berg_horizontal_shape", config_berg_horizontal_shape)
    call MPAS_pool_get_config(domain % configs, "config_berg_density", bergDensity)

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nCellsSolve", nCellsSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "tracers", tracersPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)
       call MPAS_pool_get_subpool(block % structs, "berg_state", bergStatePool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocityPool)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)
       call MPAS_pool_get_array(tracersAggregatePool, "iceVolumeCell", iceVolumeCell)
       call MPAS_pool_get_array(tracersAggregatePool, "snowVolumeCell", snowVolumeCell)

       call MPAS_pool_get_array(tracersPool, "bergMassCategory", bergMassCategory,1)
       call MPAS_pool_get_array(tracersPool, "bergLength", bergLength,1)
       call MPAS_pool_get_array(tracersPool, "bergHeight", bergHeight,1)

       call MPAS_pool_get_array(bergStatePool, "bergDraft", bergDraft)
       call MPAS_pool_get_array(bergStatePool, "oneBergArea", oneBergArea)
       call MPAS_pool_get_array(bergStatePool, "oneBergMass", oneBergMass)
       call MPAS_pool_get_array(bergStatePool, "bergVertAir", bergVertAir)
       call MPAS_pool_get_array(bergStatePool, "bergVertOcn", bergVertOcn)
       call MPAS_pool_get_array(bergStatePool, "bergVertIce", bergVertIce)

       call MPAS_pool_get_array(bergVelocityPool, "bergDraftVertex", bergDraftVertex)
       call MPAS_pool_get_array(bergVelocityPool, "oneBergAreaVertex", oneBergAreaVertex)
       call MPAS_pool_get_array(bergVelocityPool, "oneBergMassVertex", oneBergMassVertex)
       call MPAS_pool_get_array(bergVelocityPool, "bergVertAirVertex", bergVertAirVertex)
       call MPAS_pool_get_array(bergVelocityPool, "bergVertOcnVertex", bergVertOcnVertex)
       call MPAS_pool_get_array(bergVelocityPool, "bergVertIceVertex", bergVertIceVertex)

       bergSnow = 0.0_RKIND ! no snow on bergs for now
       p66 = 2.0_RKIND/3.0_RKIND
       p83 = 5.0_RKIND/6.0_RKIND
       tab2cyl = 2.0_RKIND*sqrt(2.0_RKIND/(3.0_RKIND*pii))

       do iCell = 1, nCellsSolve

          seaIceThick = 0.0_RKIND
          seaIceSnow = 0.0_RKIND

          if (iceAreaCell(iCell) > cicePuny) then

             seaIceThick = iceVolumeCell(iCell)/iceAreaCell(iCell)
             seaIceSnow = snowVolumeCell(iCell)/iceAreaCell(iCell)

          endif

          do iCategory = 1, nBergCategories

          if (bergMassCategory(iCategory,iCell) > bergMassMinimum) then

             seaIceDraft = (ciceDensitySnow*seaIceSnow + ciceDensityIce*seaIceThick)/ciceDensitySeaWater
             bergDraft(iCategory,iCell) = (ciceDensitySnow*bergSnow + bergDensity*bergHeight(iCategory,iCell))/ciceDensitySeaWater

             ! vertical (m) areas in contact with each field
             zOcn = bergDraft(iCategory,iCell) - seaIceDraft
             zIce = seaIceDraft
             zAir = bergHeight(iCategory,iCell) - zOcn - zIce

             ! horizontal area of one iceberg
             oneBergArea(iCategory,iCell) = p66*bergLength(iCategory,iCell)**2

             ! mass of one iceberg
             oneBergMass(iCategory,iCell) = oneBergArea(iCategory,iCell)*bergHeight(iCategory,iCell)*bergDensity

             if (trim(config_berg_horizontal_shape) == "cylindrical") then

                oneBergArea(iCategory,iCell) = pii*(tab2cyl*bergLength(iCategory,iCell))**2
                bergVertAir(iCategory,iCell) = 2.0_RKIND*tab2cyl*bergLength(iCategory,iCell)*zAir
                bergVertOcn(iCategory,iCell) = 2.0_RKIND*tab2cyl*bergLength(iCategory,iCell)*zOcn
                bergVertIce(iCategory,iCell) = 2.0_RKIND*tab2cyl*bergLength(iCategory,iCell)*zIce

             elseif (trim(config_berg_horizontal_shape) == "tabular1_1p5") then

                bergVertAir(iCategory,iCell) = p83*bergLength(iCategory,iCell)*zAir
                bergVertOcn(iCategory,iCell) = p83*bergLength(iCategory,iCell)*zOcn
                bergVertIce(iCategory,iCell) = p83*bergLength(iCategory,iCell)*zIce

             else

                call mpas_log_write(&
                  "berg_geometry: config_berg_horizontal_shape unknown:"//trim(config_berg_horizontal_shape), &
                  MPAS_LOG_CRIT)

             endif
          endif
          enddo ! iCategory

       enddo ! iCell

       do iCategory = 1, nBergCategories

          ! interpolate from cells to vertices
          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             bergDraftVertex(iCategory,:), &
             bergDraft(iCategory,:))

          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             oneBergAreaVertex(iCategory,:), &
             oneBergArea(iCategory,:))

          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             oneBergMassVertex(iCategory,:), &
             oneBergMass(iCategory,:))

          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             bergVertAirVertex(iCategory,:), &
             bergVertAir(iCategory,:))

          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             bergVertOcnVertex(iCategory,:), &
             bergVertOcn(iCategory,:))

          call cice_interpolate_cell_to_vertex(&
             meshPool, &
             bergVertIceVertex(iCategory,:), &
             bergVertIce(iCategory,:))

          enddo

       block => block % next
    end do ! block

  end subroutine berg_vertical_areas

!-----------------------------------------------------------------------
! Sub-cycle
!-----------------------------------------------------------------------

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_velocity_solver_subcycle
!> \brief   Solve for berg velocities
!> \author  Darin Comeau, LANL
!> \date    July 17, 2017
!> \details This subroutine solves for iceberg velocities on a subcycle.
!
!-----------------------------------------------------------------------

  subroutine berg_velocity_solver_subcycle(&
       domain, &
       clock)!{{{

    type(domain_type), intent(inout) :: &
         domain !< Input/Output:

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    integer, pointer :: &
         config_berg_dynamics_subcycle_number

    integer :: &
         iBergSubcycle

    call MPAS_pool_get_config(domain % configs, "config_berg_dynamics_subcycle_number", config_berg_dynamics_subcycle_number)

    do iBergSubcycle = 1, config_berg_dynamics_subcycle_number

       call single_subcycle_berg_velocity_solver(&
            domain, &
            clock)

    enddo

  end subroutine berg_velocity_solver_subcycle!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  single_subcycle_berg_velocity_solver
!
!> \brief   Solve for berg velocities
!> \author  Darin Comeau, LANL
!> \date    July 17, 2017
!> \details This subroutine calculates the forcing terms and solves for
!> iceberg velocities within the subcycle timestep.
!
!-----------------------------------------------------------------------

  subroutine single_subcycle_berg_velocity_solver(&
       domain, &
       clock)!{{{

    type(domain_type), intent(inout) :: &
         domain

    type(MPAS_Clock_type), intent(in) :: &
         clock !< Input:

    character(len=strKIND), pointer :: &
         config_berg_velocity_solver_type

    call MPAS_pool_get_config(domain % configs, "config_berg_velocity_solver_type", config_berg_velocity_solver_type)

    ! solve for velocity
    if (trim(config_berg_velocity_solver_type) == "forward_euler") then

       call mpas_timer_start("Calculate berg forcing terms")
       call calculate_berg_forcing(domain)
       call mpas_timer_stop("Calculate berg forcing terms")

       call mpas_timer_start("Berg velocity solver compute")
       call solve_berg_velocity_FE(domain)
       call mpas_timer_stop("Berg velocity solver compute")

    else

       call mpas_log_write(&
            "berg_velocity_solver: config_berg_velocity_solver_type unknown:"//trim(config_berg_velocity_solver_type), &
            MPAS_LOG_CRIT)

    endif

    ! halo exchange
    call mpas_timer_start("Berg velocity solver halo")
    call MPAS_dmpar_field_halo_exch(domain, 'uBergVelocity')
    call MPAS_dmpar_field_halo_exch(domain, 'vBergVelocity')
    call mpas_timer_stop("Berg velocity solver halo")

  end subroutine single_subcycle_berg_velocity_solver!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  calculate_berg_forcing
!
!> \brief   Calculate forcing terms in berg momentum equation
!> \author  Darin Comeau, LANL
!> \date    July 17, 2017
!> \details This subroutine calculates the forcing terms on the right
!> side of the iceberg momentum equation. Note only the sea surface tilt
!> forcing term does not depend on iceberg velocity, and could be pulled
!> out of the subcycle loop.
!
!-----------------------------------------------------------------------

  subroutine calculate_berg_forcing(domain)

    type(domain_type), intent(inout) :: &
         domain

    ! calculate atm forcing
    call mpas_timer_start("berg air forcing")
    call berg_atm_forcing(domain)
    call mpas_timer_stop("berg air forcing")

    ! calculate ocean forcing
    call mpas_timer_start("berg ocean forcing")
    call berg_ocn_forcing(domain)
    call mpas_timer_stop("berg ocean forcing")

    ! calculate sea ice forcing
    call mpas_timer_start("berg ice forcing")
    call berg_ice_forcing(domain)
    call mpas_timer_stop("berg ice forcing")

    ! calculate coriolis forcing
    call mpas_timer_start("berg coriolis forcing")
    call berg_cor_forcing(domain)
    call mpas_timer_stop("berg coriolis forcing")

    ! calculate surface tilt forcing
    call mpas_timer_start("berg surface tilt forcing")
    call berg_tlt_forcing(domain)
    call mpas_timer_stop("berg surface tilt forcing")

  end subroutine calculate_berg_forcing!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_atm_forcing
!
!> \brief   Calculate atmospheric forcing term
!> \author  Darin Comeau, LANL
!> \date    June 28, 2017
!> \details This subroutine calculates the atmospheric forcing term
!> in the iceberg momentum equation. The outputs bergCoefAtmU, bergCoefAtmV,
!> are in units of force (kg*m/s^2).
!
!-----------------------------------------------------------------------

  subroutine berg_atm_forcing(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool, &
         atmosCouplingPool

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocity, &
         vAirVelocity, &
         airDensity

    real(kind=RKIND), dimension(:), pointer :: &
         uAirVelocityVertex, &
         vAirVelocityVertex, &
         airDensityVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         oneBergAreaVertex, &
         bergVertAirVertex, &
         bergCoefAtmU, &
         bergCoefAtmV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    real(kind=RKIND) :: &
         relSpeed, &
         atmCoef

    integer :: &
         iVertex, &
         iCategory

    real(kind=RKIND), parameter :: &
         airVertDragCoef = 0.2_RKIND, &
         airHorizDragCoef = 0.00025_RKIND

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "atmos_coupling", atmosCouplingPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergAreaVertex", oneBergAreaVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertAirVertex", bergVertAirVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmU", bergCoefAtmU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmV", bergCoefAtmV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "uAirVelocityVertex", uAirVelocityVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vAirVelocityVertex", vAirVelocityVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "airDensityVertex", airDensityVertex)

       call MPAS_pool_get_array(atmosCouplingPool, "uAirVelocity", uAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "vAirVelocity", vAirVelocity)
       call MPAS_pool_get_array(atmosCouplingPool, "airDensity", airDensity)

       ! interpolate ait velocities from cells to vertices
       call cice_interpolate_cell_to_vertex(&
          meshPool, &
          uAirVelocityVertex, &
          uAirVelocity)

       call cice_interpolate_cell_to_vertex(&
          meshPool, &
          vAirVelocityVertex, &
          vAirVelocity)

       call cice_interpolate_cell_to_vertex(&
          meshPool, &
          airDensityVertex, &
          airDensity)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             if (bergVelocityMask(iCategory,iVertex) > 0) then

             relSpeed = sqrt((uAirVelocityVertex(iVertex) - uBergVelocity(iCategory,iVertex))**2 + (vAirVelocityVertex(iVertex) - vBergVelocity(iCategory,iVertex))**2)
             atmCoef = airDensityVertex(iVertex)*relSpeed*(airVertDragCoef*bergVertAirVertex(iCategory,iVertex) + airHorizDragCoef*oneBergAreaVertex(iCategory,iVertex))

             bergCoefAtmU(iCategory,iVertex) = atmCoef*(uAirVelocityVertex(iVertex) - uBergVelocity(iCategory,iVertex))
             bergCoefAtmV(iCategory,iVertex) = atmCoef*(vAirVelocityVertex(iVertex) - vBergVelocity(iCategory,iVertex))

             endif

          enddo

       enddo ! iVertex

       block => block % next
    end do

  end subroutine berg_atm_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_ocn_forcing
!
!> \brief   Calculate oceanic forcing term
!> \author  Darin Comeau, LANL
!> \date    July 11, 2017
!> \details This subroutine calculates the oceanic forcing term
!> in the iceberg momentum equation. The outputs bergCoefOcnU, bergCoefOcnV,
!> are in units of force (kg*m/s^2).
!
!-----------------------------------------------------------------------

  subroutine berg_ocn_forcing(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_constants, only: &
         ciceDensitySeaWater

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool, &
         oceanCouplingPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         oneBergAreaVertex, &
         bergVertOcnVertex, &
         bergCoefOcnU, &
         bergCoefOcnV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocity, &
         vOceanVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         uOceanVelocityVertex, &
         vOceanVelocityVertex

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    real(kind=RKIND) :: &
         relSpeed, &
         ocnCoef

    integer :: &
         iVertex, &
         iCategory

    real(kind=RKIND), parameter :: &
         ocnVertDragCoef = 0.425_RKIND, &
         ocnHorizDragCoef = 0.0005_RKIND

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergAreaVertex", oneBergAreaVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertOcnVertex", bergVertOcnVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnU", bergCoefOcnU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnV", bergCoefOcnV)

       call MPAS_pool_get_array(oceanCouplingPool, "uOceanVelocity", uOceanVelocity)
       call MPAS_pool_get_array(oceanCouplingPool, "vOceanVelocity", vOceanVelocity)

       call MPAS_pool_get_array(velocitySolverPool, "uOceanVelocityVertex", uOceanVelocityVertex)
       call MPAS_pool_get_array(velocitySolverPool, "vOceanVelocityVertex", vOceanVelocityVertex)

       ! interpolate cell ocean velocity to vertex
       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            uOceanVelocityVertex, &
            uOceanVelocity)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            vOceanVelocityVertex, &
            vOceanVelocity)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             if (bergVelocityMask(iCategory,iVertex) > 0) then

             relSpeed = sqrt((uOceanVelocityVertex(iVertex) - uBergVelocity(iCategory,iVertex))**2 + (vOceanVelocityVertex(iVertex) - vBergVelocity(iCategory,iVertex))**2)
             ocnCoef = ciceDensitySeaWater*relSpeed*(ocnVertDragCoef*bergVertOcnVertex(iCategory,iVertex) + ocnHorizDragCoef*oneBergAreaVertex(iCategory,iVertex))

             bergCoefOcnU(iCategory,iVertex) = ocnCoef*(uOceanVelocityVertex(iVertex) - uBergVelocity(iCategory,iVertex))
             bergCoefOcnV(iCategory,iVertex) = ocnCoef*(vOceanVelocityVertex(iVertex) - vBergVelocity(iCategory,iVertex))

             endif

          enddo

       enddo ! iVertex

       block => block % next
    end do

  end subroutine berg_ocn_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_ice_forcing
!
!> \brief   Calculate sea ice forcing term
!> \author  Darin Comeau, LANL
!> \date    July 11, 2017
!> \details This subroutine calculates the sea ice forcing term
!> in the iceberg momentum equation. The outputs bergCoefIceU, bergCoefIceV,
!> are in units of force (kg*m/s^2).
!
!-----------------------------------------------------------------------

  subroutine berg_ice_forcing(domain)

    use cice_constants, only: &
         ciceDensityIce

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool, &
         velocitySolverPool, &
         tracersAggregatePool

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    real(kind=RKIND), dimension(:), pointer :: &
         iceAreaCell, &
         areaCell

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         oneBergAreaVertex, &
         bergVertIceVertex, &
         bergCoefIceU, &
         bergCoefIceV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories, &
         vertexDegree

    integer, dimension(:,:), pointer :: &
         cellsOnVertex

    real(kind=RKIND) :: &
         relSpeed, &
         iceCoef

    integer :: &
         iCell, &
         iVertex, &
         iCategory, &
         iVertexDegree

    real(kind=RKIND), parameter :: &
         iceDragCoef = 0.5_RKIND

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)
       call MPAS_pool_get_dimension(block % dimensions, "vertexDegree", vertexDegree)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracersAggregatePool)

       call MPAS_pool_get_array(meshPool, "cellsOnVertex", cellsOnVertex)
       call MPAS_pool_get_array(meshPool, "areaCell", areaCell)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergAreaVertex", oneBergAreaVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVertIceVertex", bergVertIceVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceU", bergCoefIceU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceV", bergCoefIceV)

       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

       call MPAS_pool_get_array(tracersAggregatePool, "iceAreaCell", iceAreaCell)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             if (bergVelocityMask(iCategory,iVertex) > 0) then

             relSpeed = sqrt((uVelocity(iVertex) - uBergVelocity(iCategory,iVertex))**2 + (vVelocity(iVertex) - vBergVelocity(iCategory,iVertex))**2)
             iceCoef = ciceDensityIce*relSpeed*iceDragCoef*bergVertIceVertex(iCategory,iVertex)

             do iVertexDegree = 1, vertexDegree

                iCell = cellsOnVertex(iVertexDegree,iVertex)

                if (iceAreaCell(iCell)/areaCell(iCell) < 0.15_RKIND) iceCoef = 0.0_RKIND

             enddo

             bergCoefIceU(iCategory,iVertex) = iceCoef*(uVelocity(iVertex) - uBergVelocity(iCategory,iVertex))
             bergCoefIceV(iCategory,iVertex) = iceCoef*(vVelocity(iVertex) - vBergVelocity(iCategory,iVertex))

             endif

          enddo

       enddo ! iVertex

       block => block % next
    end do

  end subroutine berg_ice_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_cor_forcing
!
!> \brief   Calculate Coriolis forcing term
!> \author  Darin Comeau, LANL
!> \date    July 11, 2017
!> \details This subroutine calculates the Coriolis forcing term
!> in the iceberg momentum equation. The outputs bergCoefCorU, bergCoefCorV,
!> are in units of force (kg*m/s^2).
!
!-----------------------------------------------------------------------

  subroutine berg_cor_forcing(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_constants, only: &
         omega

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         latVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         oneBergMassVertex, &
         bergCoefCorU, &
         bergCoefCorV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    integer :: &
         iVertex, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergMassVertex", oneBergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorU", bergCoefCorU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorV", bergCoefCorV)

       call MPAS_pool_get_array(meshPool, "latVertex", latVertex)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             if (bergVelocityMask(iCategory,iVertex) > 0) then

             bergCoefCorU(iCategory,iVertex) =  vBergVelocity(iCategory,iVertex)*oneBergMassVertex(iCategory,iVertex)*2.0_RKIND*omega*sin(latVertex(iVertex))
             bergCoefCorV(iCategory,iVertex) = -uBergVelocity(iCategory,iVertex)*oneBergMassVertex(iCategory,iVertex)*2.0_RKIND*omega*sin(latVertex(iVertex))

             endif

          enddo

       enddo ! iVertex

       block => block % next
    end do

  end subroutine berg_cor_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_tlt_forcing
!
!> \brief   Calculate sea surface tilt forcing term
!> \author  Darin Comeau, LANL
!> \date    July 11, 2017
!> \details This subroutine calculates the sea surface tilt forcing term
!> in the iceberg momentum equation. The outputs bergCoefTltU, bergCoefTltV,
!> are in units of force (kg*m/s^2).
!
!-----------------------------------------------------------------------

  subroutine berg_tlt_forcing(domain)

    use cice_mesh, only: &
         cice_interpolate_cell_to_vertex

    use cice_constants, only: &
         ciceGravity

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         meshPool, &
         bergVelocitySolverPool, &
         oceanCouplingPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:), pointer :: &
         latVertex

    real(kind=RKIND), dimension(:,:), pointer :: &
         oneBergMassVertex, &
         bergCoefTltU, &
         bergCoefTltV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTiltU, &
         seaSurfaceTiltV

    real(kind=RKIND), dimension(:), pointer :: &
         seaSurfaceTiltVertexU, &
         seaSurfaceTiltVertexV

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    integer :: &
         iVertex, &
         iCategory

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "mesh", meshPool)
       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "ocean_coupling", oceanCouplingPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(meshPool, "latVertex", latVertex)

       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergMassVertex", oneBergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltU", bergCoefTltU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltV", bergCoefTltV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)

       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltU", seaSurfaceTiltU)
       call MPAS_pool_get_array(oceanCouplingPool, "seaSurfaceTiltV", seaSurfaceTiltV)

       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexU", seaSurfaceTiltVertexU)
       call MPAS_pool_get_array(velocitySolverPool, "seaSurfaceTiltVertexV", seaSurfaceTiltVertexV)

       ! interpolate sea surface tilt from cells to vertices
       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexU, &
            seaSurfaceTiltU)

       call cice_interpolate_cell_to_vertex(&
            meshPool, &
            seaSurfaceTiltVertexV, &
            seaSurfaceTiltV)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             if (bergVelocityMask(iCategory,iVertex) > 0) then

             bergCoefTltU(iCategory,iVertex) = -seaSurfaceTiltVertexU(iVertex)*oneBergMassVertex(iCategory,iVertex)*ciceGravity
             bergCoefTltV(iCategory,iVertex) = -seaSurfaceTiltVertexV(iVertex)*oneBergMassVertex(iCategory,iVertex)*ciceGravity

             endif

          enddo

       enddo ! iVertex

       block => block % next
    end do

  end subroutine berg_tlt_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  solve_berg_velocity_FE
!
!> \brief   Solve for berg velocities by forward Euler
!> \author  Darin Comeau, LANL
!> \date    July 17, 2017
!> \details This subroutine solves the iceberg momentum equation by
!> forward Euler.
!
!-----------------------------------------------------------------------

  subroutine solve_berg_velocity_FE(domain)

    use cice_constants, only: &
         cicePuny

    type(domain_type), intent(inout) :: &
         domain

    type(block_type), pointer :: &
         block

    type(MPAS_pool_type), pointer :: &
         bergVelocitySolverPool, &
         velocitySolverPool

    real(kind=RKIND), dimension(:,:), pointer :: &
         uBergVelocity, &
         vBergVelocity, &
         uBergVelocityInitial, &
         vBergVelocityInitial, &
         oneBergMassVertex, &
         bergCoefAtmU, &
         bergCoefAtmV, &
         bergCoefOcnU, &
         bergCoefOcnV, &
         bergCoefIceU, &
         bergCoefIceV, &
         bergCoefCorU, &
         bergCoefCorV, &
         bergCoefTltU, &
         bergCoefTltV

    integer, dimension(:,:), pointer :: &
         bergVelocityMask

    real(kind=RKIND), pointer :: &
         bergDynamicsTimeStep

    real(kind=RKIND), dimension(:), pointer :: &
         uVelocity, &
         vVelocity

    integer, pointer :: &
         nVerticesSolve, &
         nBergCategories

    integer :: &
         iVertex, &
         iCategory

    real(kind=RKIND) :: &
         rightSideU, &
         rightSideV

    block => domain % blocklist
    do while (associated(block))

       call MPAS_pool_get_dimension(block % dimensions, "nVerticesSolve", nVerticesSolve)
       call MPAS_pool_get_dimension(block % dimensions, "nBergCategories", nBergCategories)

       call MPAS_pool_get_subpool(block % structs, "berg_velocity_solver", bergVelocitySolverPool)
       call MPAS_pool_get_subpool(block % structs, "velocity_solver", velocitySolverPool)

       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocity", uBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocity", vBergVelocity)
       call MPAS_pool_get_array(bergVelocitySolverPool, "uBergVelocityInitial", uBergVelocityInitial)
       call MPAS_pool_get_array(bergVelocitySolverPool, "vBergVelocityInitial", vBergVelocityInitial)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergVelocityMask", bergVelocityMask)
       call MPAS_pool_get_array(bergVelocitySolverPool, "oneBergMassVertex", oneBergMassVertex)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmU", bergCoefAtmU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefAtmV", bergCoefAtmV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnU", bergCoefOcnU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefOcnV", bergCoefOcnV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceU", bergCoefIceU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefIceV", bergCoefIceV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorU", bergCoefCorU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefCorV", bergCoefCorV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltU", bergCoefTltU)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergCoefTltV", bergCoefTltV)
       call MPAS_pool_get_array(bergVelocitySolverPool, "bergDynamicsTimeStep", bergDynamicsTimeStep)

       call MPAS_pool_get_array(velocitySolverPool, "uVelocity", uVelocity)
       call MPAS_pool_get_array(velocitySolverPool, "vVelocity", vVelocity)

       do iVertex = 1, nVerticesSolve

          do iCategory = 1, nBergCategories

             rightSideU = 0.0_RKIND
             rightSideV = 0.0_RKIND

             rightSideU = bergCoefAtmU(iCategory,iVertex) + &
                          bergCoefOcnU(iCategory,iVertex) + &
                          bergCoefIceU(iCategory,iVertex) + &
                          bergCoefCorU(iCategory,iVertex) + &
                          bergCoefTltU(iCategory,iVertex)
             rightSideV = bergCoefAtmV(iCategory,iVertex) + &
                          bergCoefOcnV(iCategory,iVertex) + &
                          bergCoefIceV(iCategory,iVertex) + &
                          bergCoefCorV(iCategory,iVertex) + &
                          bergCoefTltV(iCategory,iVertex)

             if (bergVelocityMask(iCategory,iVertex) == 1 .and. oneBergMassVertex(iCategory,iVertex) > cicePuny) then

                uBergVelocity(iCategory,iVertex) = uBergVelocityInitial(iCategory,iVertex) + &
                                                   rightSideU*bergDynamicsTimeStep/oneBergMassVertex(iCategory,iVertex)

                vBergVelocity(iCategory,iVertex) = vBergVelocityInitial(iCategory,iVertex) + &
                                                   rightSideV*bergDynamicsTimeStep/oneBergMassVertex(iCategory,iVertex)

                ! write to log if iceberg velocities get too large, indicating instability (consider changing to abort)
                if (abs(uBergVelocity(iCategory,iVertex))>1.0e2_RKIND) &
                call mpas_log_write(" Berg velocity large, uBergVelocity $r ", realArgs=(/uBergVelocity(iCategory,iVertex)/))

                if (abs(vBergVelocity(iCategory,iVertex))>1.0e2_RKIND) &
                call mpas_log_write(" Berg velocity large, vBergVelocity $r ", realArgs=(/vBergVelocity(iCategory,iVertex)/))

             elseif (bergVelocityMask(iCategory,iVertex) == 2) then

                uBergVelocity(iCategory,iVertex) = uVelocity(iVertex)
                vBergVelocity(iCategory,iVertex) = vVelocity(iVertex)

             else

                uBergVelocity(iCategory,iVertex) = 0.0_RKIND
                vBergVelocity(iCategory,iVertex) = 0.0_RKIND

             endif

             uBergVelocityInitial(iCategory,iVertex) = uBergVelocity(iCategory,iVertex)
             vBergVelocityInitial(iCategory,iVertex) = vBergVelocity(iCategory,iVertex)

          enddo ! iCategory

       enddo ! iVertex

       block => block % next
    enddo

  end subroutine solve_berg_velocity_FE!}}}

!-----------------------------------------------------------------------

end module cice_berg_velocity_solver
