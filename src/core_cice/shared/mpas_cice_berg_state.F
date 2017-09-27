!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_berg_state
!
!> \brief   Update iceberg model
!> \author  Darin Comeau, LANL
!> \date    26 Sept 2017
!> \details This module is intended to update the iceberg state variables.
!
!-----------------------------------------------------------------------

module cice_berg_state

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_dmpar
   use mpas_timer
   use mpas_log, only: mpas_log_write

   use mpas_stream_manager
   use mpas_io_units

   implicit none

   private
   save

   public :: &
      berg_calving_event, &
      berg_dimensions_update, &
      berg_area_update

contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_calving_event
!
!> \brief   Initialize icebergs size dimensions
!> \author  Darin Comeau, LANL
!> \date    26 Sept 2017
!> \details This routine is intended to add iceberg mass through calving.
!> Note this will be moved out of this module, as will be called outside 
!> of initialization.
!
!-----------------------------------------------------------------------

  subroutine berg_calving_event(block)

    type(block_type), intent(inout) :: &
         block

    type(MPAS_pool_type), pointer :: &
         mesh, &
         tracers, &
         tracers_aggregate, &
         berg_calving

    integer, pointer :: &
         nCellsSolve, &
         nBergCategories

    real(kind=RKIND), dimension(:,:), pointer :: &
         bergMassCategory, &
         bergAreaCategory, &
         bergLength, &
         bergHeight, &
         bergCalvingMass

    real(kind=RKIND), dimension(:), pointer :: &
         bergMassCell, &
         bergAreaCell, &
         areaCell

    real(kind=RKIND), dimension(:), pointer :: &
         bergCalvingLength, &
         bergCalvingHeight

    integer, dimension(:), pointer :: &
         bergCalvingMask

    real(kind=RKIND), pointer :: &
         bergDensity

    integer :: &
         iCell, &
         iCategory

    call MPAS_pool_get_config(block % configs, "config_berg_density", bergDensity)

    call MPAS_pool_get_subpool(block % structs, "mesh", mesh)
    call MPAS_pool_get_subpool(block % structs, "tracers", tracers)
    call MPAS_pool_get_subpool(block % structs, "tracers_aggregate", tracers_aggregate)
    call MPAS_pool_get_subpool(block % structs, "berg_calving", berg_calving)

    call MPAS_pool_get_dimension(mesh, "nCellsSolve", nCellsSolve)
    call MPAS_pool_get_dimension(mesh, "nBergCategories", nBergCategories)

    call MPAS_pool_get_array(mesh, "areaCell", areaCell)
    call MPAS_pool_get_array(tracers, "bergMassCategory", bergMassCategory,1)
    call MPAS_pool_get_array(tracers, "bergAreaCategory", bergAreaCategory,1)
    call MPAS_pool_get_array(tracers, "bergLength", bergLength,1)
    call MPAS_pool_get_array(tracers, "bergHeight", bergHeight,1)
    call MPAS_pool_get_array(tracers_aggregate, "bergMassCell", bergMassCell)
    call MPAS_pool_get_array(tracers_aggregate, "bergAreaCell", bergAreaCell)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMask", bergCalvingMask)
    call MPAS_pool_get_array(berg_calving, "bergCalvingMass", bergCalvingMass)
    call MPAS_pool_get_array(berg_calving, "bergCalvingLength", bergCalvingLength)
    call MPAS_pool_get_array(berg_calving, "bergCalvingHeight", bergCalvingHeight)

    do iCell = 1, nCellsSolve

       if (bergCalvingMask(iCell) == 1) then

          do iCategory = 1, nBergCategories

             ! average dimensions with existing bergs
             call berg_dimensions_update(bergMassCategory(iCategory,iCell), &
                                         bergCalvingMass(iCategory,iCell), &
                                         bergLength(iCategory,iCell), &
                                         bergCalvingLength(iCategory), &
                                         bergHeight(iCategory,iCell), &
                                         bergCalvingHeight(iCategory), &
                                         bergLength(iCategory,iCell), &
                                         bergHeight(iCategory,iCell))

             ! add calving mass
             bergMassCategory(iCategory,iCell) = bergMassCategory(iCategory,iCell) &
                                               + bergCalvingMass(iCategory,iCell)

          enddo

          ! integrated quantities
          bergMassCell(iCell) = sum(bergMassCategory(:,iCell))

          ! update iceberg area
          call berg_area_update(nBergCategories, &
                                bergMassCategory(:,iCell), &
                                bergHeight(:,iCell), &
                                bergDensity, &
                                areaCell(iCell), &
                                bergAreaCategory(:,iCell), &
                                bergAreaCell(iCell))

       endif

    enddo

  end subroutine berg_calving_event

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_dimensions_update
!
!> \brief   Update iceberg dimensions
!> \author  Darin Comeau, LANL
!> \date    Sep 26, 2017
!> \details This subroutine updates iceberg dimensions as the weighted
!> average of two iceberg masses in the same cell.
!
!-----------------------------------------------------------------------

  subroutine berg_dimensions_update(&
       bergMass1, &
       bergMass2, &
       bergLength1, &
       bergLength2, &
       bergHeight1, &
       bergHeight2, &
       newBergLength, &
       newBergHeight)

    real(kind=RKIND), intent(in) :: &
         bergMass1, &
         bergMass2, &
         bergLength1, &
         bergLength2, &
         bergHeight1, &
         bergHeight2

    real(kind=RKIND), intent(out) :: &
         newBergLength, &
         newBergHeight

    real(kind=RKIND) :: &
         weight

    weight = bergMass1 / (bergMass1 + bergMass2)

    newBergLength = weight*bergLength1 + (1.0_RKIND-weight)*bergLength2
    newBergHeight = weight*bergHeight1 + (1.0_RKIND-weight)*bergHeight2

  end subroutine berg_dimensions_update 

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  berg_area_update
!
!> \brief   Update iceberg area
!> \author  Darin Comeau, LANL
!> \date    Sep 26, 2017
!> \details This subroutine updates iceberg area state variable from the
!> iceberg mass state variable.
!
!-----------------------------------------------------------------------

  subroutine berg_area_update(&
       nBergCategories, &
       bergMassCategory, &
       bergHeight, &
       bergDensity, &
       areaCell, &
       bergAreaCategory, &
       bergAreaCell)

    integer, intent(in) :: &
         nBergCategories

    real(kind=RKIND), dimension(:), intent(in) :: &
         bergMassCategory, &
         bergHeight

    real(kind=RKIND), intent(in) :: &
         bergDensity, &
         areaCell

    real(kind=RKIND), dimension(:), intent(out) :: &
         bergAreaCategory

    real(kind=RKIND), intent(out) :: &
         bergAreaCell

    integer :: &
         iCategory

    do iCategory = 1, nBergCategories

       bergAreaCategory(iCategory) = bergMassCategory(iCategory) / (bergHeight(iCategory)*bergDensity)
       bergAreaCategory(iCategory) = bergAreaCategory(iCategory) / areaCell

    enddo

    bergAreaCell = sum(bergAreaCategory(:))

    if (bergAreaCell > 1.0_RKIND) then
       call mpas_log_write(&
            "berg area greater than one: $r ", realArgs=(/bergAreaCell/), &
            messageType = MPAS_LOG_CRIT)

    endif

  end subroutine berg_area_update

!-----------------------------------------------------------------------

end module cice_berg_state
