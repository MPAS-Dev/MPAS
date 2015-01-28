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
   use mpas_constants

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

   public :: mpas_init_mesh_paraview_catalyst

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

   contains

!***********************************************************************
!
!  routine mpas_init_mesh_paraview_catalyst
!
!> \brief   Create mesh variables for paraview catalyst
!> \author  Patricia Fasel, Mark Petersen
!> \date    April 2014
!> \details 
!>  Create all the mesh information needed by the paraview catalyst
!>  for in-situ visualization and analysis.
!
!-----------------------------------------------------------------------
   subroutine mpas_init_mesh_paraview_catalyst(domain,longitude_periodic_split)!{{{
   
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: domain
      real (kind=RKIND), intent(in) :: longitude_periodic_split

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      ! Config options
      integer, pointer :: config_num_halos

      ! Dimensions
      integer, pointer :: nCells, nEdges, nVertices, maxEdges, vertexDegree, nVertLevels
      integer :: timeLevel
      integer :: nGhostCell, nGhostVertex
      integer :: iHalo, iCell, indx
      integer, dimension(:), allocatable :: cellGhost, cellHalo
      integer, dimension(:), allocatable :: vertexGhost, vertexHalo
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnVertex, vertexMask, verticesOnCell, cellMask

      real (kind=RKIND) :: split
      real (kind=RKIND), dimension(:), pointer ::  &
         xCell, yCell, zCell, xVertex, yVertex, zVertex, lonCell, latCell, lonVertex, latVertex
      real (kind=RKIND), dimension(:), pointer ::  &
         lonCellDegrees, latCellDegrees, lonVertexDegrees, latVertexDegrees

      type (block_type), pointer :: block

      type (mpas_pool_type), pointer :: meshPool

      type (mpas_exchange_list), pointer :: exchListPtr

      timeLevel=1  ! mrp change to input argument later

      call mpas_pool_get_config(domain % configs,'config_num_halos',config_num_halos)

      block => domain % blocklist
      do while (associated(block))

         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)

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

         ! Allocate arrays to hold lon/lat dimensions to send to catalyst
         allocate(lonCellDegrees(nCells+1))
         allocate(latCellDegrees(nCells+1))
         allocate(lonVertexDegrees(nVertices+1))
         allocate(latVertexDegrees(nVertices+1))

         ! mrp: I thought I could change split in Mercator projection, 
         !  but it is not completely working yet.

         ! make value of split between 0 and 360
          split = mod(longitude_periodic_split,360.0)
          lonCellDegrees = mod(lonCell*180./pii + split,360.0) - split
          latCellDegrees = mod(latCell*180./pii + split,360.0) - split
          lonVertexDegrees = mod(lonVertex*180./pii + split,360.0) - split
          latVertexDegrees = mod(latVertex*180./pii + split,360.0) - split

!        lonCellDegrees = lonCell*180./pii
!        latCellDegrees = latCell*180./pii
!        lonVertexDegrees = lonVertex*180./pii
!        latVertexDegrees = latVertex*180./pii

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
                       lonCellDegrees, \
                       latCellDegrees, \
                       lonVertexDegrees, \
                       latVertexDegrees, \
                       nEdgesOnCell, \
                       cellsOnVertex, \
                       vertexMask, \
                       verticesOnCell, \
                       cellMask)

         deallocate(cellGhost)
         deallocate(cellHalo)
         deallocate(vertexGhost)
         deallocate(vertexHalo)

         deallocate(lonCellDegrees)
         deallocate(latCellDegrees)
         deallocate(lonVertexDegrees)
         deallocate(latVertexDegrees)

       block => block % next
    end do

   end subroutine mpas_init_mesh_paraview_catalyst!}}}

end module mpas_paraview_catalyst

! vim: foldmethod=marker
