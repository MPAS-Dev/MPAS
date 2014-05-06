! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_paraview_catalyst
!
!> \brief   Paraview catalyst in-situ adapter
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  These routines set up and pass data structures to the 
!>  paraview catalyst for in-situ visualization and analysis.
!
!-----------------------------------------------------------------------

module mpas_paraview_catalyst

   use mpas_grid_types
   use mpas_dmpar
   use mpas_dmpar_types
   use mpas_sort
   use mpas_configure
   use mpas_io_streams
   use mpas_io_output

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: mpas_insitu_initialize, &
             mpas_insitu_create_geometry, &
             mpas_insitu_load_data, &
             mpas_insitu_coprocess, &
             mpas_insitu_finalize

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

   contains

!***********************************************************************
!
!  routine mpas_insitu_initialize
!
!> \brief   Initialize paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!> Initialize paraview catalyst
!
!-----------------------------------------------------------------------
   subroutine mpas_insitu_initialize!{{{

      call mpas_initialize

   end subroutine mpas_insitu_initialize!}}}

!***********************************************************************
!
!  routine mpas_insitu_create_geometry
!
!> \brief   Create mesh variables for paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  Create all the mesh information needed by the paraview catalyst
!>  for in-situ visualization and analysis.
!
!-----------------------------------------------------------------------
   subroutine mpas_insitu_create_geometry(domain)!{{{
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
!mrp change inout to in when everything is working
      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! Config options
      integer, pointer :: config_num_halos

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertices, maxEdges, vertexDegree, nVertLevels, &
         indexTemperature, indexSalinity, num_tracers
      integer :: timeLevel
      integer :: nGhostCell, nGhostVertex
      integer :: iHalo, iCell, indx
      integer, dimension(:), allocatable :: cellGhost, cellHalo
      integer, dimension(:), allocatable :: vertexGhost, vertexHalo
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnVertex, vertexMask, verticesOnCell, cellMask

      real (kind=RKIND), dimension(:), pointer ::  &
         xCell, yCell, zCell, xVertex, yVertex, zVertex, lonCell, latCell, lonVertex, latVertex

      type (block_type), pointer :: block
      type (field2DReal), pointer :: densityField, pressureField, circulationField, relativeVorticityField
      type (field3DReal), pointer :: tracersField

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool

      type (mpas_exchange_list), pointer :: exchListPtr

      timeLevel=1  ! mrp change to input argument later

      call mpas_pool_get_config(domain % configs,'config_num_halos',config_num_halos)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         ! Collect the number of halo cells in all levels
         nGhostCell = 0
         do iHalo = 1, config_num_halos
            exchListPtr => block % parinfo % cellsToRecv % halos(iHalo) % exchList
            do while(associated(exchListPtr))
               nGhostCell = nGhostCell + exchListPtr % nList
               exchListPtr => exchListPtr % next
            end do
         end do

         nGhostVertex = 0
         do iHalo = 1, config_num_halos
            exchListPtr => block % parinfo % verticesToRecv % halos(iHalo) % exchList
            do while(associated(exchListPtr))
               nGhostVertex = nGhostVertex + exchListPtr % nList
               exchListPtr => exchListPtr % next
            end do
         end do

         ! Allocate arrays to hold halo cell index and halo number
         allocate(cellGhost(nGhostCell))
         allocate(cellHalo(nGhostCell))
         allocate(vertexGhost(nGhostVertex))
         allocate(vertexHalo(nGhostVertex))

         ! Fill in halo cells indices and halo number
         indx = 1
         do iHalo = 1, config_num_halos
            exchListPtr => block %parinfo % cellsToRecv % halos(iHalo) % exchList
            do while(associated(exchListPtr))
               do iCell = 1, exchListPtr % nList
                  cellGhost(indx) = exchListPtr % destList(iCell)
                  cellHalo(indx) = iHalo
                  indx = indx + 1
               end do
               exchListPtr => exchListPtr % next
            end do
         end do

         indx = 1
         do iHalo = 1, config_num_halos
            exchListPtr => block %parinfo % verticesToRecv % halos(iHalo) % exchList
            do while(associated(exchListPtr))
               do iCell = 1, exchListPtr % nList
                  vertexGhost(indx) = exchListPtr % destList(iCell)
                  vertexHalo(indx) = iHalo
                  indx = indx + 1
               end do
               exchListPtr => exchListPtr % next
            end do
         end do

         call mpas_pool_get_dimension(block % dimensions, 'nCells', nCells)
         call mpas_pool_get_dimension(block % dimensions, 'nEdges', nEdges)
         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nVertices', nVertices)
         call mpas_pool_get_dimension(block % dimensions, 'maxEdges', maxEdges)
         call mpas_pool_get_dimension(block % dimensions, 'vertexDegree', vertexDegree)

         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'zCell', zCell)
         call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
         call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
         call mpas_pool_get_array(meshPool, 'latCell', latCell)
         call mpas_pool_get_array(meshPool, 'lonVertex', lonVertex)
         call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
         call mpas_pool_get_array(meshPool, 'vertexMask', vertexMask)
         call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
         call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      print *,'COPROCESS create geometry',nCells,nVertices,maxEdges,vertexDegree,nVertLevels
      call coprocessor_create_grid( \
                       nCells, maxEdges, nGhostCell, cellGhost, cellHalo, \
                       nVertices, vertexDegree, nGhostVertex, vertexGhost, vertexHalo, \
                       nVertLevels, \
                       xCell, \
                       yCell, \
                       zCell, \
                       xVertex, \
                       yVertex, \
                       zVertex, \
                       lonCell, \
                       latCell, \
                       lonVertex, \
                       latVertex, \
                       nEdgesOnCell, \
                       cellsOnVertex, \
                       vertexMask, \
                       verticesOnCell, \
                       cellMask)

         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)
         call mpas_pool_get_field(statePool, 'tracers', tracersField, timeLevel)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
         call coprocessor_register_tracer_data( \
                       indexTemperature, \
                       tracersField % constituentNames(indexTemperature), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, &   !!!!!!!!!!mrp -1 is temp
                       tracersField % array)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)
         call coprocessor_register_tracer_data( \
                       indexSalinity, \
                       tracersField % constituentNames(indexSalinity), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, &   !!!!!!!!!!mrp -1 is temp
                       tracersField % array)

         call mpas_pool_get_field(diagnosticsPool, 'density', densityField)
         call coprocessor_register_data( \
                       densityField % fieldName, \
                       densityField % dimSizes(1), \
                       densityField % dimSizes(2)-1, &   !!!! mrp +1 is temp
                       densityField % array)


         call mpas_pool_get_field(diagnosticsPool, 'pressure', pressureField)
         call coprocessor_register_data( \
                       pressureField % fieldName, \
                       pressureField % dimSizes(1), \
                       pressureField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       pressureField % array)

         call mpas_pool_get_field(diagnosticsPool, 'circulation', circulationField)
         call coprocessor_register_data( \
                       circulationField % fieldName, \
                       circulationField % dimSizes(1), \
                       circulationField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       circulationField % array)

         call mpas_pool_get_field(diagnosticsPool, 'relativeVorticity', relativeVorticityField)
         call coprocessor_register_data( \
                       relativeVorticityField % fieldName, \
                       relativeVorticityField % dimSizes(1), \
                       relativeVorticityField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       relativeVorticityField % array)

         deallocate(cellGhost)
         deallocate(cellHalo)
         deallocate(vertexGhost)
         deallocate(vertexHalo)

       block => block % next
    end do

   end subroutine mpas_insitu_create_geometry!}}}

!***********************************************************************
!
!  routine mpas_insitu_load_data
!
!> \brief   Pass data to paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  Pass data to paraview catalyst
!
!-----------------------------------------------------------------------
   subroutine mpas_insitu_load_data(domain, itime)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
!mrp change inout to in when everything is working
      type (domain_type), intent(inout) :: domain
      integer, intent(inout) :: itime
!      integer, intent(in) :: timeLevel  ! mrp add later

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! Dimensions
      integer, pointer :: indexTemperature, indexSalinity, num_tracers

      type (block_type), pointer :: block
      integer :: timeLevel
      type (field2DReal), pointer :: densityField, pressureField, circulationField, relativeVorticityField
      type (field3DReal), pointer :: tracersField

      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: diagnosticsPool

timeLevel=1  ! mrp change to input argument later

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)

         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)
         call mpas_pool_get_field(statePool, 'tracers', tracersField, timeLevel)
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
         call coprocessor_add_tracer_data( \
                       itime, \
                       indexTemperature, \
                       tracersField % constituentNames(indexTemperature), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, &   !!!!!!!!!!mrp -1 is temp
                       tracersField % array)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)
         call coprocessor_add_tracer_data( \
                       itime, \
                       indexSalinity, \
                       tracersField % constituentNames(indexSalinity), \
                       num_tracers, \
                       tracersField % dimSizes(2), \
                       tracersField % dimSizes(3)-1, &   !!!!!!!!!!mrp -1 is temp
                       tracersField % array)

         call mpas_pool_get_field(diagnosticsPool, 'density', densityField)
         call coprocessor_add_data( \
                       itime, \
                       densityField % fieldName, \
                       densityField % dimSizes(1), \
                       densityField % dimSizes(2)-1, &  !!!!! mrp -1 is temp
                       densityField % array)

         call mpas_pool_get_field(diagnosticsPool, 'pressure', pressureField)
         call coprocessor_add_data( \
                       itime, \
                       pressureField % fieldName, \
                       pressureField % dimSizes(1), \
                       pressureField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       pressureField % array)

         call mpas_pool_get_field(diagnosticsPool, 'circulation', circulationField)
         call coprocessor_add_data( \
                       itime, \
                       circulationField % fieldName, \
                       circulationField % dimSizes(1), \
                       circulationField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       circulationField % array)

         call mpas_pool_get_field(diagnosticsPool, 'relativeVorticity', relativeVorticityField)
         call coprocessor_add_data( \
                       itime, \
                       relativeVorticityField % fieldName, \
                       relativeVorticityField % dimSizes(1), \
                       relativeVorticityField % dimSizes(2)-1, &   !!!!!!!!!!mrp -1 is temp
                       relativeVorticityField % array)

     block => block % next
   end do

   end subroutine mpas_insitu_load_data!}}}

!***********************************************************************
!
!  routine mpas_insitu_coprocess
!
!> \brief   Deterimine if it is time to load data
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!> Deterimine if it is time to load data
!
!-----------------------------------------------------------------------
   subroutine mpas_insitu_coprocess(domain, itime)!{{{

      type (domain_type), intent(inout) :: domain
      integer, intent(inout) :: itime
      integer :: doWork

      doWork = 0
      print *,'checking for COPROCESS ',itime
      call mpas_check_coprocess(itime, doWork)
      if(doWork .eq. 0) then
         return
      endif
      ! we actually need to do coprocessing
      ! the grids should already be loaded.
      call mpas_insitu_load_data(domain, itime)
      print *,'COPROCESS ',itime
      call mpas_coprocess

   end subroutine mpas_insitu_coprocess!}}}

!***********************************************************************
!
!  routine mpas_insitu_finalize
!
!> \brief   finalize paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  finalize paraview catalyst
!
!-----------------------------------------------------------------------
   subroutine mpas_insitu_finalize()!{{{
 
      print *,'COPROCESS finalize'
      call mpas_finalize
 
   end subroutine mpas_insitu_finalize!}}}
 
end module mpas_paraview_catalyst

! vim: foldmethod=marker
