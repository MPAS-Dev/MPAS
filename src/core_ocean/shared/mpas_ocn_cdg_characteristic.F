! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_ocn_cdg_characteristic
!
!> \brief  characteristic advectection of vertices and quadrature points
!          for the CDG advection scheme
!> \author Dave Lee
!> \date   December 2015
!> \details
!>  This module contains routines for advecting vertices and quadrature
!>  points on the sphere
!
!-----------------------------------------------------------------------

#define CDG_INT_ORDER_VERT 1
#define CDG_INT_ORDER_HORIZ 1

module ocn_cdg_characteristic

   use mpas_kind_types
   use mpas_dmpar
   use mpas_io_units

   use mpas_geometry_utils
   use ocn_particle_helpers

   use ocn_cdg_types
   use ocn_cdg_triangle

   implicit none
   private
   save 
   
   public :: ocn_cdg_char_integrate,                   &
             ocn_cdg_char_integrate_planar,            &
             ocn_cdg_char_departure_points,            &
             ocn_cdg_char_departure_points_planar,     &
             ocn_cdg_char_integrate_edge_planar,       &
             ocn_cdg_char_integrate_dep_pts_vertical,  &
             ocn_cdg_char_integrate_quad_pts_vertical

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_char_departure_points
!
!> \brief  Find departure points associated with each vertex
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This routine traces backward trajectories to locate the departure
!>  points associated with each vertex.
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_char_departure_points(meshPool, cdgPool, dt, departurePoints)!{{{

      type (mpas_pool_type), pointer, intent(in) :: meshPool
      type (mpas_pool_type), pointer, intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:), intent(out) :: departurePoints

      integer, pointer :: nVertices, nVertLevels, nCells, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, verticesOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, xCell, yCell, zCell

      integer :: iLevel, iVertex, iCell, jCell
      real (kind=RKIND), dimension(3) :: x1, x2, x3, x4, u1, u2, u3, u4, dx
      real (kind=RKIND), dimension(:,:), pointer :: uVertex, vVertex, wVertex, uVertexMid, vVertexMid, wVertexMid, areaBArray
    
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)

      call mpas_pool_get_array(cdgPool, 'uVertexVelocityNew', uVertex)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityNew', vVertex)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityNew', wVertex)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityMid', wVertexMid)

      do iVertex = 1, nVertices
         x1(1) = xVertex(iVertex)
         x1(2) = yVertex(iVertex)
         x1(3) = zVertex(iVertex)

         do iLevel = 1, maxLevelVertexBot(iVertex)
            u1(1) = uVertex(iLevel,iVertex)
            u1(2) = vVertex(iLevel,iVertex)
            u1(3) = wVertex(iLevel,iVertex)

            departurePoints(:,iLevel,iVertex) = x1(:)

            if (CDG_INT_ORDER_HORIZ .eq. 1) then
               !call particle_horizontal_movement(meshPool, departurePoints(:,iLevel,iVertex), -dt*u1)
               call mpas_particle_horizontal_movement(departurePoints(:,iLevel,iVertex), -dt*u1, .true.)
               cycle
            end if

            do iCell = 1, vertexDegree
               jCell = cellsOnVertex(iCell,iVertex)
               if (jCell .le. nCells .and. jCell .gt. 0) exit
            end do

            x2(:) = x1(:)
            !call particle_horizontal_movement(meshPool, x2, -0.5_RKIND*dt*u1)
            call mpas_particle_horizontal_movement(x2, -0.5_RKIND*dt*u1, .true.)
            !call find_and_interp_velocity(nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexMid, vVertexMid, wVertexMid, &
            !                              areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, radius, x2, jCell, u2)
            call find_and_interp_velocity(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexMid, vVertexMid, wVertexMid, &
                                          areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, x2, jCell, u2)

            dx(:) = -dt*u2(:)
            !call particle_horizontal_movement(meshPool, departurePoints(:,iLevel,iVertex), dx)
            call mpas_particle_horizontal_movement(departurePoints(:,iLevel,iVertex), dx, .true.)
         end do
      end do

   end subroutine ocn_cdg_char_departure_points!}}}

   subroutine ocn_cdg_char_departure_points_planar(meshPool, cdgPool, dt, departurePoints)!{{{

      type (mpas_pool_type), pointer, intent(in) :: meshPool
      type (mpas_pool_type), pointer, intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:), intent(out) :: departurePoints

      integer, pointer :: nVertices, nCells, vertexDegree
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelVertexBot
      integer, dimension(:,:), pointer :: cellsOnCell, cellsOnVertex, verticesOnCell, cellsOnEdge, edgesOnCell, edgesOnVertex
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, xCell, yCell, zCell

      logical :: done
      integer :: iLevel, iVertex, iCell, jCell, iEdge, jEdge
      real (kind=RKIND), dimension(3) :: x1, x2, u1, u2, dx
      real (kind=RKIND), dimension(:,:), pointer :: uVertexNew, uVertexMid, vVertexNew, vVertexMid, areaBArray

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnVertex', edgesOnVertex)
      call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityNew', uVertexNew)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityNew', vVertexNew)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)

      do iVertex = 1, nVertices
         do iLevel = 1, maxLevelVertexBot(iVertex)
            x1(1) = xVertex(iVertex)
            x1(2) = yVertex(iVertex)

            u1(1) = uVertexNew(iLevel,iVertex)
            u1(2) = vVertexNew(iLevel,iVertex)

            if (CDG_INT_ORDER_HORIZ .eq. 1) then
               departurePoints(:,iLevel,iVertex) = x1(:) - dt*u1(:)
               cycle
            end if

            x2(:) = x1(:) - 0.5_RKIND*dt*u1(:)

            done = .false.
            iCell = 0
            do while (.not. done)
               iCell = iCell + 1
               jCell = cellsOnVertex(iCell,iVertex)
               if (iCell .eq. vertexDegree) then
                  done = .true.
               end if
               if (jCell .gt. 0 .and. jCell .le. nCells) then
                  done = .true.
               end if
            end do

            call find_and_interp_velocity_planar(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, uVertexMid, vVertexMid, &
               areaBArray, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, x2, jCell, u2)

            departurePoints(:,iLevel,iVertex) = x1(:) - dt*u2(:)
         end do
      end do

   end subroutine ocn_cdg_char_departure_points_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_char_integrate
!
!> \brief  integrate a set of points along velocity characteristics at
!          a specified order of accuracy
!> \author Dave Lee
!> \date   May 2016
!> \details
!   Assumes characteristics follow great circle arcs on the sphere
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_char_integrate(meshPool, np, xp, yp, zp, uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, &!{{{
                                     nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, &
                                     verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, iCell, areaBArray, dt, xf)

      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in) :: np, iLevel, iCell, nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnCell
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:), intent(in) :: xp, yp, zp, xVertex, yVertex, zVertex, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), intent(in) :: uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, areaBArray
      real (kind=RKIND), dimension(:,:), intent(out) :: xf

      integer :: ip, jCell
      real (kind=RKIND), dimension(3) :: x1, x2, u1, u2, dx

      do ip = 1, np
         jCell = iCell

         x1(1) = xp(ip)
         x1(2) = yp(ip)
         x1(3) = zp(ip)
         xf(:,ip) = x1(:)

         !call find_and_interp_velocity(nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexOld, vVertexOld, wVertexOld, &
         !                              areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, radius, x1, jCell, u1)
         call find_and_interp_velocity(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexOld, vVertexOld, wVertexOld, &
                                       areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, x1, jCell, u1)

         if (CDG_INT_ORDER_HORIZ .eq. 1) then
            !call particle_horizontal_movement(meshPool, xf(:,ip), dt*u1(:))
            call mpas_particle_horizontal_movement(xf(:,ip), dt*u1(:), .true.)
            cycle
         end if

         x2(:) = x1(:)
         !call particle_horizontal_movement(meshPool, x2, 0.5_RKIND*dt*u1)
         call mpas_particle_horizontal_movement(x2, 0.5_RKIND*dt*u1, .true.)
         !call find_and_interp_velocity(nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexMid, vVertexMid, wVertexMid, &
         !                              areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, radius, x2, jCell, u2)
         call find_and_interp_velocity(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertexMid, vVertexMid, wVertexMid, &
                                       areaBArray, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, x2, jCell, u2)

         !call particle_horizontal_movement(meshPool, xf(:,ip), dt*u2(:))
         call mpas_particle_horizontal_movement(xf(:,ip), dt*u2(:), .true.)
      end do

   end subroutine ocn_cdg_char_integrate!}}}

   subroutine ocn_cdg_char_integrate_planar(meshPool, np, xp, yp, uVertex, vVertex, uVertexMid, vVertexMid, nCells, xCell, yCell, zCell, xVertex, yVertex, &!{{{
                                            verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, iCell, areaBArray, dt, xf)

      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(in) :: np, iLevel, iCell, nCells
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnCell, cellsOnVertex
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:), intent(in) :: xp, yp, xVertex, yVertex, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), intent(in) :: uVertex, vVertex, areaBArray, uVertexMid, vVertexMid
      real (kind=RKIND), dimension(:,:), intent(out) :: xf

      integer :: ip, jCell
      real (kind=RKIND), dimension(3) :: x1, x2, x3, x4, u1, u2, u3, u4

      do ip = 1, np
         x1(1) = xp(ip)
         x1(2) = yp(ip)

         jCell = iCell
         call find_and_interp_velocity_planar(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, uVertex, vVertex, &
            areaBArray, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, x1, jCell, u1)

         if (CDG_INT_ORDER_HORIZ .eq. 1) then
            xf(1,ip) = x1(1) + dt*u1(1)
            xf(2,ip) = x1(2) + dt*u1(2)
            cycle
         end if

         x2(:) = x1(:) + 0.5_RKIND*dt*u1(:)
         call find_and_interp_velocity_planar(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, uVertexMid, vVertexMid, &
            areaBArray, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, x2, jCell, u2)

         xf(1,ip) = x1(1) + dt*u2(1)
         xf(2,ip) = x1(2) + dt*u2(2)
      end do

   end subroutine ocn_cdg_char_integrate_planar!}}}

   subroutine find_and_interp_velocity(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, uVertex, vVertex, wVertex, areaBArray, &!{{{
                                       verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, xi, iCell, ui)

      type(mpas_pool_type), intent(in) :: meshPool
      integer, intent(in) :: nCells, iLevel
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), intent(in) :: xCell, yCell, zCell, xVertex, yVertex, zVertex, xi
      real (kind=RKIND), dimension(:,:), intent(in) :: uVertex, vVertex, wVertex, areaBArray
      integer, intent(inout) :: iCell
      real (kind=RKIND), dimension(3), intent(out) :: ui

      integer :: iVertex, jVertex, nEdges
      !real (kind=RKIND), dimension(6) :: uVertCell, vVertCell, wVertCell, xVertCell, yVertCell, zVertCell, lambda
      real (kind=RKIND), dimension(12) :: uVertCell, vVertCell, wVertCell, xVertCell, yVertCell, zVertCell, lambda
      real (kind=RKIND), dimension(3,12) :: vertCell

      !call mpas_get_nearby_cell_index_sphere(nCells, xCell, yCell, zCell, xi(1), xi(2), xi(3), iCell, cellsOnCell, nEdgesOnCell, radius)
      call mpas_get_nearby_cell_index(nCells, xCell, yCell, zCell, xi(1), xi(2), xi(3), meshPool, iCell, cellsOnCell, nEdgesOnCell)

      nEdges = nEdgesOnCell(iCell)
      do iVertex = 1, nEdges
         jVertex = verticesOnCell(iVertex,iCell)
         vertCell(1,iVertex) = xVertex(jVertex)
         vertCell(2,iVertex) = yVertex(jVertex)
         vertCell(3,iVertex) = zVertex(jVertex)
         uVertCell(iVertex) = uVertex(iLevel,jVertex)
         vVertCell(iVertex) = vVertex(iLevel,jVertex)
         wVertCell(iVertex) = wVertex(iLevel,jVertex)
      end do
      !lambda(1:nEdges) = mpas_wachspress_coordinates_spherical(nEdges, xVertCell, yVertCell, zVertCell, xi(1), xi(2), xi(3), radius, areaBArray(:,iCell))
      lambda(1:nEdges) = mpas_wachspress_coordinates(nEdges, vertCell, xi, meshPool, areaBArray(:,iCell))
      ui(1) = mpas_wachspress_interpolate(lambda(1:nEdges), uVertCell(1:nEdges))
      ui(2) = mpas_wachspress_interpolate(lambda(1:nEdges), vVertCell(1:nEdges))
      ui(3) = mpas_wachspress_interpolate(lambda(1:nEdges), wVertCell(1:nEdges))

   end subroutine find_and_interp_velocity!}}}

   subroutine find_and_interp_velocity_planar(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, uVertexVelocity, vVertexVelocity, areaBArray, &!{{{
                                              verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, xi, iCell, ui)

      type(mpas_pool_type), intent(in) :: meshPool
      integer, intent(in) :: nCells, iLevel
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), intent(in) :: xCell, yCell, zCell, xVertex, yVertex, xi
      real (kind=RKIND), dimension(:,:), intent(in) :: uVertexVelocity, vVertexVelocity, areaBArray
      integer, intent(inout) :: iCell
      real (kind=RKIND), dimension(2), intent(out) :: ui

      integer :: iVertex, jVertex, nEdges
      !real (kind=RKIND), dimension(6) :: uVertCell, vVertCell, xVertCell, yVertCell, lambda
      real (kind=RKIND), dimension(12) :: uVertCell, vVertCell, wVertCell, xVertCell, yVertCell, lambda
      real (kind=RKIND), dimension(3,12) :: vertCell

      !call mpas_get_nearby_cell_index_planar(xi(1), xi(2), xCell, yCell, nIncCells, incCells, iCell)
      call mpas_get_nearby_cell_index(nCells, xCell, yCell, zCell, xi(1), xi(2), xi(3), meshPool, iCell, cellsOnCell, nEdgesOnCell)

      nEdges = nEdgesOnCell(iCell)
      do iVertex = 1, nEdges
         jVertex = verticesOnCell(iVertex,iCell)
         vertCell(1,iVertex) = xVertex(jVertex)
         vertCell(2,iVertex) = yVertex(jVertex)
         vertCell(3,iVertex) = 0.0_RKIND
         uVertCell(iVertex) = uVertexVelocity(iLevel,jVertex)
         vVertCell(iVertex) = vVertexVelocity(iLevel,jVertex)
         wVertCell(iVertex) = 0.0_RKIND
      end do

      !lambda(1:nEdges) = mpas_wachspress_coordinates_planar(nEdges, xVertCell, yVertCell, xi(1), xi(2), areaBArray(:,iCell))
      lambda(1:nEdges) = mpas_wachspress_coordinates(nEdges, vertCell, xi, meshPool, areaBArray(:,iCell))

      ui(1) = mpas_wachspress_interpolate(lambda(1:nEdges), uVertCell(1:nEdges))
      ui(2) = mpas_wachspress_interpolate(lambda(1:nEdges), vVertCell(1:nEdges))

   end subroutine find_and_interp_velocity_planar!}}}

   subroutine ocn_cdg_char_integrate_edge_planar(meshPool, cdgPool, dt, departureEdges) !{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:), intent(out) :: departureEdges

      integer, pointer :: nEdges, nCells
      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell, cellsOnVertex, cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, xEdge, yEdge, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: uVertexMid, vVertexMid, edgeNormVecs, normVelocitySum, areaBArray

      integer :: iLevel, iEdge, iCell, cell1, cell2
      real (kind=RKIND) :: eNorm, uNorm
      real (kind=RKIND), dimension(3) :: xh, ui, uh

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocitySum)

      do iEdge = 1, nEdges
         do iLevel = 1, maxLevelEdgeTop(iEdge)
            ui(1) = normVelocitySum(iLevel,iEdge)*edgeNormVecs(1,iEdge)
            ui(2) = normVelocitySum(iLevel,iEdge)*edgeNormVecs(2,iEdge)

            if (CDG_INT_ORDER_HORIZ .eq. 1) then
               departureEdges(:,iLevel,iEdge) = dt*ui(:)
               cycle
            end if

            xh(1) = xEdge(iEdge) - 0.5_RKIND*dt*ui(1)
            xh(2) = yEdge(iEdge) - 0.5_RKIND*dt*ui(2)

            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)
            if (cell1 .gt. nCells .or. cell2 .gt. nCells) then
               departureEdges(:,iLevel,iEdge) = 0.0_RKIND
               cycle
            else
               if (normVelocitySum(iLevel,iEdge) .gt. 0.0_RKIND) then
                  iCell = cell1
               else
                  iCell = cell2
               end if
            end if
            call find_and_interp_velocity_planar(meshPool, nCells, xCell, yCell, zCell, xVertex, yVertex, uVertexMid, vVertexMid, &
               areaBArray, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, xh, iCell, uh)

            departureEdges(:,iLevel,iEdge) = dt*uh(:)*edgeNormVecs(:,iEdge)

            ! If 2nd order velocity is in opposite direction then fall back to 1st order
            uNorm = ui(1)*edgeNormVecs(1,iEdge) + ui(2)*edgeNormVecs(2,iEdge)
            eNorm = uh(1)*edgeNormVecs(1,iEdge) + uh(2)*edgeNormVecs(2,iEdge)
            if (eNorm*uNorm .lt. 0.0_RKIND) then
               departureEdges(:,iLevel,iEdge) = dt*ui(:)
            end if
         end do
      end do

   end subroutine ocn_cdg_char_integrate_edge_planar!}}}

   subroutine ocn_cdg_char_integrate_quad_pts_vertical(quad, nVertLevels, iLevel, zFullOld, zFullMid, wTopOld, wTopMid, dt, zb, dz, zf)!{{{

      type (cdgQuad3D), intent(in) :: quad
      integer, intent(in) :: nVertLevels, iLevel
      real (kind=RKIND), intent(in) :: dt, zb, dz
      real (kind=RKIND), dimension(:), intent(in) :: zFullOld, zFullMid, wTopOld, wTopMid
      real (kind=RKIND), dimension(:), intent(out) :: zf

      logical found
      integer :: jLevel, iQuad
      real (kind=RKIND), dimension(4) :: wi, zi, zh, zl, wh

      zi(1:quad % nv) = zb + dz*quad % qzv(:)
      zl(1:quad % nv) = (zi(1:quad % nv) - zFullOld(iLevel+1))/(zFullOld(iLevel) - zFullOld(iLevel+1))
      wi(1:quad % nv) = (1.0_RKIND - zl(1:quad % nv))*wTopOld(iLevel+1) + zl(1:quad % nv)*wTopOld(iLevel)
      zh(1:quad % nv) = zi(1:quad % nv) + 0.5_RKIND*dt*wi(1:quad % nv)

      if (CDG_INT_ORDER_VERT .eq. 1) then
         zf(1:quad % nv) = zi(1:quad % nv) + dt*wi(1:quad % nv)
         return
      end if

      do iQuad = 1, quad % nv
         found = .false.
         jLevel = 1
         do while (.not. found)
            if (jLevel .gt. nVertLevels) then
               write(0,*) 'cdg: ERROR! vertical quadrature point search'
               stop
            end if

            if (zFullMid(jLevel+1) .lt. zh(iQuad)) then
               found = .true.

               zl(iQuad) = (zh(iQuad) - zFullMid(jLevel+1))/(zFullMid(jLevel) - zFullMid(jLevel+1))
               ! TODO using mid time level vertical velocities breaks consistency
               !wh(iQuad) = (1.0_RKIND - zl(iQuad))*wTopMid(jLevel+1) + zl(iQuad)*wTopMid(jLevel)
               wh(iQuad) = (1.0_RKIND - zl(iQuad))*wTopOld(jLevel+1) + zl(iQuad)*wTopOld(jLevel)
               zf(iQuad) = zi(iQuad) + dt*wh(iQuad)
            end if

            jLevel = jLevel + 1
         end do
      end do

   end subroutine ocn_cdg_char_integrate_quad_pts_vertical!}}}

   subroutine ocn_cdg_char_integrate_dep_pts_vertical(nVertLevels, layerThicknessNew, layerThicknessMid, wTopNew, wTopMid, botDepth, dt, depPts)!{{{

      integer, intent(in) :: nVertLevels
      real (kind=RKIND), intent(in) :: botDepth, dt
      real (kind=RKIND), dimension(:), intent(in) :: layerThicknessNew, layerThicknessMid, wTopNew, wTopMid
      real (kind=RKIND), dimension(:), intent(out) :: depPts

      logical :: found
      integer :: iLevel, jLevel
      real (kind=RKIND) :: zh, zl, wh
      real (kind=RKIND), dimension(256) :: zFullNew, zFullMid

      depPts(1:nVertLevels) = zFullNew(1:nVertLevels)

      zFullNew(nVertLevels+1) = -botDepth
      zFullMid(nVertLevels+1) = -botDepth
      do iLevel = nVertLevels, 1, -1
         zFullNew(iLevel) = zFullNew(iLevel+1) + layerThicknessNew(iLevel)
         zFullMid(iLevel) = zFullMid(iLevel+1) + layerThicknessMid(iLevel)
      end do

      do iLevel = 2, nVertLevels
         if (CDG_INT_ORDER_VERT .eq. 1) then
            depPts(iLevel) = zFullNew(iLevel) - dt*wTopNew(iLevel)
         else
            zh = zFullNew(iLevel) - 0.5_RKIND*dt*wTopNew(iLevel)
            if (zh .gt. zFullNew(1)) then
               zh = zFullNew(1)
               write(0,*) 'cdg: ERROR! vertical departure point above sea surface', iLevel, nVertLevels
               stop
            else if (zh .lt. zFullNew(nVertLevels+1)) then
               zh = zFullNew(nVertLevels+1)
               write(0,*) 'cdg: ERROR! vertical departure point below sea bottom', iLevel, nVertLevels
               stop
            end if

            found = .false.
            jLevel = 1
            do while (.not. found)
               if (jLevel .gt. nVertLevels) then
                  write(0,*) 'cdg: ERROR! vertical departure point search', iLevel, nVertLevels, zFullNew(iLevel), zFullNew(nVertLevels+1), dt*wTopNew(iLevel), zh, dt, wTopNew(iLevel)
                  stop
               end if

               !if (zFullMid(jLevel+1) .lt. zh) then
               if (zFullNew(jLevel+1) .lt. zh) then
                  found = .true.

                  !zl = (zh - zFullMid(jLevel+1))/layerThicknessMid(jLevel)
                  zl = (zh - zFullNew(jLevel+1))/layerThicknessNew(jLevel)
                  if (zl .lt. -1.0E-12 .or. zl .gt. 1.0+1.0E-12) then
                     write(0,*) 'cdg: ERROR! vertical departure point local coord', zl
                     write(0,*) 'cdg: ERROR!', iLevel, jLevel, zFullNew(iLevel), zh, zFullNew(jLevel+1)
                     stop
                  end if
                  !wh = (1.0_RKIND - zl)*wTopMid(jLevel+1) + zl*wTopMid(jLevel)
                  wh = (1.0_RKIND - zl)*wTopNew(jLevel+1) + zl*wTopNew(jLevel)
                  depPts(iLevel) = zFullNew(iLevel) - dt*wh
                  if (depPts(iLevel) .gt. zFullNew(1)) then
                     depPts(iLevel) = zFullNew(1)
                     write(0,*) 'cdg: ERROR! vertical departure point above sea surface'
                     stop
                  else if (depPts(iLevel) .lt. zFullNew(nVertLevels+1)) then
                     depPts(iLevel) = zFullNew(nVertLevels+1)
                     write(0,*) 'cdg: ERROR! vertical departure point below sea bottom'
                     stop
                  end if
               end if

               jLevel = jLevel + 1
            end do
         end if
      end do

   end subroutine ocn_cdg_char_integrate_dep_pts_vertical!}}}

end module ocn_cdg_characteristic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
