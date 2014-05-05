! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module mpas_paraview_catalyst

   use mpas_grid_types
   use mpas_dmpar
   use mpas_dmpar_types
   use mpas_sort
   use mpas_configure
   use mpas_io_streams
   use mpas_io_output

   contains

   subroutine mpas_insitu_initialize!{{{

      call mpas_initialize

   end subroutine mpas_insitu_initialize!}}}

   subroutine mpas_insitu_create_geometry(domain)!{{{
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block
      type (mpas_exchange_list), pointer :: exchListPtr
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      integer :: nCells, nEdges, nVertices, maxEdges, vertexDegree, nVertLevels
      integer :: nGhostCell, nGhostVertex
      integer :: nHaloLayers, iHalo, iCell, indx
      integer :: indexT, indexS, numTracers
      integer, dimension(:), allocatable :: cellGhost, cellHalo
      integer, dimension(:), allocatable :: vertexGhost, vertexHalo


    block => domain % blocklist
    do while (associated(block))

      nCells = block % mesh % nCells
      nVertices = block % mesh % nVertices
      nVertLevels = block % mesh % nVertLevels
      nEdges = block % mesh % nEdges
      maxEdges = block % mesh % maxEdges
      vertexDegree = block % mesh % vertexDegree

      ! Collect the number of halo cells in all levels
      nHaloLayers = config_num_halos
      nGhostCell = 0
      do iHalo = 1, nHaloLayers
        exchListPtr => block %parinfo % cellsToRecv % halos(iHalo) % exchList
        do while(associated(exchListPtr))
          nGhostCell = nGhostCell + exchListPtr % nList
          exchListPtr => exchListPtr % next
        end do
      end do

      nGhostVertex = 0
      do iHalo = 1, nHaloLayers
        exchListPtr => block %parinfo % verticesToRecv % halos(iHalo) % exchList
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
      do iHalo = 1, nHaloLayers
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
      do iHalo = 1, nHaloLayers
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

      print *,'COPROCESS create geometry',nCells,nVertices,maxEdges,vertexDegree,nVertLevels
      call coprocessor_create_grid( \
                       nCells, maxEdges, nGhostCell, cellGhost, cellHalo, \
                       nVertices, vertexDegree, nGhostVertex, vertexGhost, vertexHalo, \
                       nVertLevels, \
                       block % mesh % xCell % array, \
                       block % mesh % yCell % array, \
                       block % mesh % zCell % array, \
                       block % mesh % xVertex % array, \
                       block % mesh % yVertex % array, \
                       block % mesh % zVertex % array, \
                       block % mesh % lonCell % array, \
                       block % mesh % latCell % array, \
                       block % mesh % lonVertex % array, \
                       block % mesh % latVertex % array, \
                       block % mesh % nEdgesOnCell % array, \
                       block % mesh % cellsOnVertex % array, \
                       block % mesh % vertexMask % array, \
                       block % mesh % verticesOnCell % array, \
                       block % mesh % cellMask % array)

      tracers => block % state % time_levs(2) % state % tracers % array
      numTracers = size(tracers, dim=1)
      indexT = block % state % time_levs(1) % state % index_temperature
      indexS = block % state % time_levs(1) % state % index_salinity
      call coprocessor_register_tracer_data( \
                       indexT, \
                       block % state % time_levs(1) % state % tracers % constituentNames(indexT), \
                       numTracers, \
                       block % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block % state % time_levs(1) % state % tracers % array)
      call coprocessor_register_tracer_data( \
                       indexS, \
                       block % state % time_levs(1) % state % tracers % constituentNames(indexS), \
                       numTracers, \
                       block % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block % state % time_levs(1) % state % tracers % array)

      call coprocessor_register_data( \
                       block % diagnostics % density % fieldName, \
                       block % diagnostics % density % dimSizes(1), \
                       block % diagnostics % density % dimSizes(2), \
                       block % diagnostics % density % array)
      call coprocessor_register_data( \
                       block % diagnostics % pressure % fieldName, \
                       block % diagnostics % pressure % dimSizes(1), \
                       block % diagnostics % pressure % dimSizes(2), \
                       block % diagnostics % pressure % array)
      call coprocessor_register_data( \
                       block % diagnostics % circulation % fieldName, \
                       block % diagnostics % circulation % dimSizes(1), \
                       block % diagnostics % circulation % dimSizes(2), \
                       block % diagnostics % circulation % array)
      call coprocessor_register_data( \
                       block % diagnostics % relativeVorticity % fieldName, \
                       block % diagnostics % relativeVorticity % dimSizes(1), \
                       block % diagnostics % relativeVorticity % dimSizes(2), \
                       block % diagnostics % relativeVorticity % array)

      deallocate(cellGhost)
      deallocate(cellHalo)
      deallocate(vertexGhost)
      deallocate(vertexHalo)

       block => block % next
    end do

   end subroutine mpas_insitu_create_geometry!}}}

   subroutine mpas_insitu_load_data(domain, itime)!{{{

      implicit none
 
      type (domain_type), intent(inout) :: domain
      type (block_type), pointer :: block
      integer, intent(inout) :: itime
      integer :: indexT, indexS, numTracers
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

    block => domain % blocklist
    do while (associated(block))

      tracers => block % state % time_levs(2) % state % tracers % array
      numTracers = size(tracers, dim=1)
      indexT = block % state % time_levs(1) % state % index_temperature
      indexS = block % state % time_levs(1) % state % index_salinity

      call coprocessor_add_tracer_data( \
                       itime, \
                       indexT, \
                       block % state % time_levs(1) % state % tracers % constituentNames(indexT), \
                       numTracers, \
                       block % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block % state % time_levs(1) % state % tracers % array)
      call coprocessor_add_tracer_data( \
                       itime, \
                       indexS, \
                       block % state % time_levs(1) % state % tracers % constituentNames(indexS), \
                       numTracers, \
                       block % state % time_levs(1) % state % tracers % dimSizes(2), \
                       block % state % time_levs(1) % state % tracers % dimSizes(3), \
                       block % state % time_levs(1) % state % tracers % array)

      call coprocessor_add_data( \
                       itime, \
                       block % diagnostics % density % fieldName, \
                       block % diagnostics % density % dimSizes(1), \
                       block % diagnostics % density % dimSizes(2), \
                       block % diagnostics % density % array)
      call coprocessor_add_data( \
                       itime, \
                       block % diagnostics % pressure % fieldName, \
                       block % diagnostics % pressure % dimSizes(1), \
                       block % diagnostics % pressure % dimSizes(2), \
                       block % diagnostics % pressure % array)
      call coprocessor_add_data( \
                       itime, \
                       block % diagnostics % circulation % fieldName, \
                       block % diagnostics % circulation % dimSizes(1), \
                       block % diagnostics % circulation % dimSizes(2), \
                       block % diagnostics % circulation % array)
      call coprocessor_add_data( \
                       itime, \
                       block % diagnostics % relativeVorticity % fieldName, \
                       block % diagnostics % relativeVorticity % dimSizes(1), \
                       block % diagnostics % relativeVorticity % dimSizes(2), \
                       block % diagnostics % relativeVorticity % array)

     block => block % next
   end do

   end subroutine mpas_insitu_load_data!}}}

   subroutine mpas_insitu_coprocess(domain, itime)!{{{

      implicit none

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


   subroutine mpas_insitu_finalize()!{{{
 
      implicit none
 
      print *,'COPROCESS finalize'
      call mpas_finalize
 
   end subroutine mpas_insitu_finalize!}}}
 
end module mpas_paraview_catalyst
