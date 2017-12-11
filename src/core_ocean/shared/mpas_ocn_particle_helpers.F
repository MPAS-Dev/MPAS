! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  mpas_ocn_particle_helpers
!
!> \brief  MPAS ocean CDG advection utility routines
!> \author Dave Lee and Phillip Wolfram
!> \date   November 2017
!> \details
!   Contains utility routines to for lagrangian particle integration
!   and 
!
!-----------------------------------------------------------------------

module ocn_particle_helpers

   use mpas_kind_types
   use mpas_dmpar
   use mpas_io_units

   use mpas_sort
   use mpas_hash

   use mpas_geometry_utils
   use mpas_vector_reconstruction

   use ocn_constants

   implicit none
   private
   save 
   
   public :: mpas_initialize_wachspress_coefficients, &
             mpas_particle_horizontal_movement,       &
             ocn_vertex_reconstruction,               &
             old_ocn_spherical_area

contains

!***********************************************************************
!
!  routine initalize_fields
!
!> \brief   Initialize Wachspress coefficients
!> \author  Phillip Wolfram
!> \date    01/26/2015
!> \details     
!>  This routine inializes the B_i Wachspress coefficients which are 
!>  static in time
!         
!-----------------------------------------------------------------------

   subroutine mpas_initialize_wachspress_coefficients(domain, err) !{{{

      implicit none

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(in) :: domain

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: meshPool
      integer  :: nVertices, iCell, i, im1, i0, ip1, iVertex
      integer, pointer :: nCells
      real (kind=RKIND), dimension(:), allocatable :: xv,yv,zv
      real (kind=RKIND), dimension(3) :: v1, v2, v3
      integer, dimension(:,:), pointer :: verticesOnCell
      integer, dimension(:), pointer :: nCellVerticesArray
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray
      logical, pointer :: on_a_sphere, is_periodic
      real(kind=RKIND), pointer :: x_period, y_period, radius

      err = 0

      block => domain % blocklist
      do while (associated(block))
        ! setup pointers / get block
        call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
        call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
        call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
        call mpas_pool_get_config(meshPool, 'x_period', x_period)
        call mpas_pool_get_config(meshPool, 'y_period', y_period)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nCellVerticesArray)
        call mpas_pool_get_config(meshPool, 'sphere_radius', radius)
        call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
        call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
        call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
        call mpas_pool_get_array(meshPool, 'xCell', xCell)
        call mpas_pool_get_array(meshPool, 'yCell', yCell)
        call mpas_pool_get_array(meshPool, 'zCell', zCell)
        call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)

        ! compute B_i coefficients
        areaBArray = 0.0_RKIND
        do iCell = 1, nCells
          nVertices = nCellVerticesArray(iCell)
          allocate(xv(nVertices), yv(nVertices), zv(nVertices))
          if (on_a_sphere .or. .not. is_periodic) then
            xv = xVertex(verticesOnCell(:,iCell))
            yv = yVertex(verticesOnCell(:,iCell))
            zv = zVertex(verticesOnCell(:,iCell))
          else
            do iVertex=1,nVertices
              xv(iVertex) = mpas_fix_periodicity(xVertex(verticesOnCell(iVertex,iCell)), &
                xCell(iCell), x_period)
              yv(iVertex) = mpas_fix_periodicity(yVertex(verticesOnCell(iVertex,iCell)), &
                yCell(iCell), y_period)
              zv(iVertex) = zVertex(verticesOnCell(iVertex,iCell))
            end do
          end if
          do i = 1, nVertices
            ! compute first area B_i
            ! get vertex indices
            im1 = mod(nVertices + i - 2, nVertices) + 1
            i0  = mod(nVertices + i - 1, nVertices) + 1
            ip1 = mod(nVertices + i    , nVertices) + 1

            ! precompute B_i areas
            ! always the same because B_i independent of xp,yp,zp
            v1(1) = xv(im1)
            v1(2) = yv(im1)
            v1(3) = zv(im1)
            v2(1) = xv(i0)
            v2(2) = yv(i0)
            v2(3) = zv(i0)
            v3(1) = xv(ip1)
            v3(2) = yv(ip1)
            v3(3) = zv(ip1)
            if (.not. on_a_sphere) then
              areaBArray(i, iCell) = mpas_triangle_signed_area(v1, v2, v3, meshPool)
            else
              areaBArray(i, iCell) = old_ocn_spherical_area(v1(1), v1(2), v1(3), v2(1), v2(2), v2(3), v3(1), v3(2), v3(3), radius)
            end if
          end do
          deallocate(xv, yv, zv)

        end do

        block => block %  next
      end do

   end subroutine mpas_initialize_wachspress_coefficients !}}}

!***********************************************************************
!
!  routine particle_horizontal_movement
!
!> \brief   Compute horizontal movement for particle so particle stays
!>          in spherical shell
!> \author  Phillip Wolfram
!> \date    05/20/2014
!> \details
!>  This routine returns the particle position pParticle corresponding
!>  to an initial particle position pParticle for a Cartesian movemnt
!>  dpParticle.  If the calculation is onSphere, then the distance
!>  |dpParticle| must be along the great circle route of pParticle
!>  and the projection of pParticle + dpParticle on the spherical
!>  shell corresponding to pParticle.
!
!-----------------------------------------------------------------------
  subroutine mpas_particle_horizontal_movement(pParticle, dpParticle, onSphere, is_periodic, x_period, y_period) !{{{

     implicit none

     !-----------------------------------------------------------------
     ! input variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(in) :: dpParticle
     logical, intent(in) :: onSphere
     logical, intent(in), optional :: is_periodic
     real (kind=RKIND), intent(in), optional :: x_period, y_period

     !-----------------------------------------------------------------
     ! input / output variables
     !-----------------------------------------------------------------
     real (kind=RKIND), dimension(:), intent(inout) :: pParticle

     !-----------------------------------------------------------------
     ! local variables
     !-----------------------------------------------------------------
     real (kind=RKIND) :: lenPath, arcLen
     real (kind=RKIND) :: radiusShell
     real (kind=RKIND), dimension(size(pParticle)) :: pParticleTemp
     real (kind=RKIND), dimension(size(pParticle)) :: pParticleInterp
     real (kind=RKIND) :: alpha
     real (kind=RKIND), parameter :: eps=1e-10_RKIND
     character (len=StrKIND) :: message
     ! choosen based on the parameters, note that we loose about 6 - 7 units of precision because R is so large!
     ! therefore, eps = 1e-10 is conservative, if not too high!  this just helps with numerical stability
     !dpParticle    =   -4.2428037617887103E-011   4.3076544298828060E-011   5.0760704444480953E-011
     !pParticle     =    4444887.2990309987       -891565.00525021972        4476665.3916420965
     !pParticleTemp =    4444887.2990309987       -891565.00525021972        4476665.3916420965
     !mpas_arc_length =    0.0000000000000000      lenPath =    7.8945399869363434E-011

     ! may need a condition to determine if we need to project back to the sphere
     if(onSphere) then
       ! need to make sure new point is on the spherical shell

       ! get path length
       lenPath = sqrt(sum(dpParticle*dpParticle))

       ! consider case of particle not moving (need to have this code here in general)
       !if (lenPath < eps) then
       ! this is ok because this is only the case if the points are the same.  If there is a
       ! numerical instability it probably should be handled differently.
       if (lenPath < eps) then
         return
       end if

       ! get radius of particle's horizontal shell
       radiusShell = sqrt(sum(pParticle*pParticle))

       ! project endpoint to spherical shell containing pParticle
       pParticleTemp = pParticle + dpParticle
       pParticleTemp = ( radiusShell / sqrt(sum(pParticleTemp*pParticleTemp)) ) * pParticleTemp

       ! compute alpha parameter for spherical interpolant / extrapolant
       write(message, *) 'mpas_arc_length = ',  mpas_arc_length( pParticle(1), pParticle(2), pParticle(3), &
                         pParticleTemp(1), pParticleTemp(2), pParticleTemp(3)), 'lenPath = ', lenPath
       arcLen = mpas_arc_length(pParticle(1),pParticle(2),pParticle(3) , pParticleTemp(1),pParticleTemp(2),pParticleTemp(3))
       if (arcLen > eps) then
         alpha = lenPath / arcLen
       else
         return
       endif

       ! compute final position based on spherical interpolant
       call mpas_spherical_linear_interp(pParticleInterp, pParticle, pParticleTemp, alpha)
       pParticle = pParticleInterp
     else
       ! we are just on a plane so there is no need for spherical interpolation to keep
       ! the new particle location on a spherical shell
       pParticle = pParticle + dpParticle

       ! periodic fix to make sure particle advection stays in domain
       if (present(is_periodic) .and. present(x_period) .and. present(y_period)) then
         if (is_periodic) then
           pParticle(1) = mpas_fix_periodicity(pParticle(1), x_period/2.0_RKIND, x_period)
           pParticle(2) = mpas_fix_periodicity(pParticle(2), y_period/2.0_RKIND, y_period)
           !pParticle(3) = pParticle(3)
         end if
       end if
     endif

   end subroutine mpas_particle_horizontal_movement!}}}

   real (kind=RKIND) function old_ocn_spherical_area(ax, ay, az, bx, by, bz, cx, cy, cz, radius)!{{{

      real (kind=RKIND), intent(in) :: ax, ay, az, bx, by, bz, cx, cy, cz, radius

      real (kind=RKIND) :: AB, BC, CA       ! Lenth of triangle sides
      real (kind=RKIND) :: S                ! Length of semi-perimeter of TRI(ABC)
      real (kind=RKIND), dimension(3) :: aa, bb, cc

      !AB = old_ocn_arc_length(ax,ay,az , bx,by,bz, radius) / radius
      !BC = old_ocn_arc_length(bx,by,bz , cx,cy,cz, radius) / radius
      !CA = old_ocn_arc_length(cx,cy,cz , ax,ay,az, radius) / radius
      AB = mpas_arc_length(ax,ay,az , bx,by,bz) / radius
      BC = mpas_arc_length(bx,by,bz , cx,cy,cz) / radius
      CA = mpas_arc_length(cx,cy,cz , ax,ay,az) / radius
      S  = 0.5_RKIND * (AB + BC + CA)

      ! need max for degenerate case where one side length is zero (which makes arg to sqrt go negative)
      old_ocn_spherical_area = 4.0_RKIND*radius*radius*   &
        atan( sqrt(max(0.0_RKIND, tan(0.5_RKIND*S) * tan(0.5_RKIND*(S-AB)) * tan(0.5_RKIND*(S-BC)) * tan(0.5_RKIND*(S-CA))) ) )

      old_ocn_spherical_area = abs(old_ocn_spherical_area)

   end function old_ocn_spherical_area!}}}

!***********************************************************************
!
!  routine ocn_vertex_reconstruction
!
!> \brief   Reconstruct vertex velocity driver / interface
!> \author  Phillip Wolfram
!> \date    03/27/2014
!> \details
!>  Purpose: reconstruct vector field at vertex locations based on
!>           particular choice of reconstruction method
!>  Input: mesh meta data and vector component data residing at cell edges
!>         initialize_weights logical is to determine if weights should be initialized
!>  Output: reconstructed vector field (measured in X,Y,Z) located at vertices
!-----------------------------------------------------------------------
  subroutine ocn_vertex_reconstruction(filterNum, meshPool, scratchPool, particlesPool, layerThickness, u, & !{{{
                                       uvReconstructX, uvReconstructY, uvReconstructZ )

    implicit none

    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    type (mpas_pool_type), pointer, intent(in) :: scratchPool !< Input: Scratch variables
    type (mpas_pool_type), pointer, intent(in) :: particlesPool !< Input: particlefield variables
    integer, intent(in) :: filterNum  ! filtering strength employed
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness !< Input: layerThickness on cells
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: u !< Input: Velocity field on edges (normalVelocity)
    type (field2DReal), pointer, intent(out) :: uvReconstructX !< Output: X Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(out) :: uvReconstructY !< Output: Y Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(out) :: uvReconstructZ !< Output: Z Component of velocity reconstructed to vertices

    ! could add additional reconstruction techniques here with switch if desired

    ! assumption is made that mpas_init_reconstruct was previously called
    call ocn_RBFvertex(meshPool, filterNum, layerThickness, u, uvReconstructX, uvReconstructY, uvReconstructZ, .false., &
                       scratchPool, particlesPool)

  end subroutine ocn_vertex_reconstruction!}}}

!***********************************************************************
!
!  routine ocn_RBFvertex
!
!> \brief   Reconstruct vertex velocity using linear interpolation of
!>          RBFs reconstruction at cell centers
!> \author  Phillip Wolfram, Todd Ringler
!> \date    03/26/2014
!> \details
!>  Purpose: reconstruct vector field at vertex locations based on radial basis functions
!>  Input: mesh meta data and vector component data residing at cell edges
!>         initialize_weights logical is to determine if weights should be initialized
!>  Output: reconstructed vector field (measured in X,Y,Z) located at vertices
!-----------------------------------------------------------------------
  subroutine ocn_RBFvertex(meshPool, filterNum, layerThickness, u, uvReconstructX, uvReconstructY, uvReconstructZ, & !{{{
                           initialize_weights, scratchPool, particlesPool)

    implicit none

    ! inputs
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    type (mpas_pool_type), pointer, intent(in) :: scratchPool
    type (mpas_pool_type), pointer, intent(in) :: particlesPool !< Input: particlefield variables
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: u !< Input: Velocity field on edges
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness !< Input: layerThickness on cells
    integer, intent(in) :: filterNum  !< number of times to filter
    logical, intent(in) :: initialize_weights !< Input: Determine if weights for RBF should be pre-computed

    ! outputs
    type (field2DReal), pointer, intent(out) :: uvReconstructX !< Output: X Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(out) :: uvReconstructY !< Output: Y Component of velocity reconstructed to vertices
    type (field2DReal), pointer, intent(out) :: uvReconstructZ !< Output: Z Component of velocity reconstructed to vertices

    ! local / temporary arrays needed in the compute procedure
    type (field2DReal), pointer :: &
      ucReconstructX, ucReconstructY, ucReconstructZ, ucReconstructZonal, ucReconstructMeridional ! cell center values
    type (field2DReal), pointer :: ucStore, vcStore, wcStore
    type (field2DInteger), pointer :: boundaryVertex, boundaryCell, boundaryVertexGlobal

    ! get pointers
    call mpas_pool_get_field(scratchPool, 'ucReconstructX', ucReconstructX)
    call mpas_pool_get_field(scratchPool, 'ucReconstructY', ucReconstructY)
    call mpas_pool_get_field(scratchPool, 'ucReconstructZ', ucReconstructZ)
    call mpas_pool_get_field(scratchPool, 'ucReconstructZonal', ucReconstructZonal)
    call mpas_pool_get_field(scratchPool, 'ucReconstructMeridional', ucReconstructMeridional)
    call mpas_pool_get_field(scratchPool, 'boundaryVertexGlobal', boundaryVertexGlobal)

    ! allocate memory
    call mpas_allocate_scratch_field(ucReconstructX, .True.)
    call mpas_allocate_scratch_field(ucReconstructY, .True.)
    call mpas_allocate_scratch_field(ucReconstructZ, .True.)
    call mpas_allocate_scratch_field(ucReconstructZonal, .True.)
    call mpas_allocate_scratch_field(ucReconstructMeridional, .True.)
    call mpas_allocate_scratch_field(boundaryVertexGlobal, .True.)

    ucReconstructX % array = 0.0_RKIND
    ucReconstructY % array = 0.0_RKIND
    ucReconstructZ % array = 0.0_RKIND

    ! initialize weights (should be pre-initialized)
    if (initialize_weights) then
      call mpas_init_reconstruct(meshPool)
    end if

    ! get cell center reconstructed RBF values
    call mpas_reconstruct(meshPool, u, ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
      ucReconstructZonal % array, ucReconstructMeridional % array)

    ! need to do exchange for uc components (we don't use Zonal / Meridional for this calculation)
    call mpas_dmpar_exch_halo_field(ucReconstructX)
    call mpas_dmpar_exch_halo_field(ucReconstructY)
    call mpas_dmpar_exch_halo_field(ucReconstructZ)

    ! get boundaries
    call mpas_pool_get_field(meshPool,'boundaryVertex', boundaryVertex)
    call mpas_pool_get_field(meshPool,'boundaryCell', boundaryCell)


    if (filternum > 0) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! filter the cell velocity field
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ocn_second_order_shapiro_filter_ops(filterNum, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! store filter data &
      ! write data to file for output
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_pool_get_field(particlesPool, 'filteredVelocityU', ucStore)
      call mpas_pool_get_field(particlesPool, 'filteredVelocityV', vcStore)
      call mpas_pool_get_field(particlesPool, 'filteredVelocityW', wcStore)
      ucStore % array = ucReconstructX % array
      vcStore % array = ucReconstructY % array
      wcStore % array = ucReconstructZ % array

    end if

    ! interpolate to vertex locations for use in Wachspress
    call ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex % array, boundaryCell % array, &
      ucReconstructX % array, ucReconstructY % array, ucReconstructZ % array, &
      uvReconstructX % array, uvReconstructY % array, uvReconstructZ % array)

    ! handle boundary vertices (should be zero).  Can potentially remove if mpas_init_block initializes to 0 vs -1e34
    boundaryVertexGlobal % array = boundaryVertex % array
    call mpas_dmpar_exch_halo_field(boundaryVertexGlobal)
    ! definite change between these fields!
    uvReconstructX % array = uvReconstructX % array * (1.0_RKIND - boundaryVertexGlobal % array)
    uvReconstructY % array = uvReconstructY % array * (1.0_RKIND - boundaryVertexGlobal % array)
    uvReconstructZ % array = uvReconstructZ % array * (1.0_RKIND - boundaryVertexGlobal % array)

    ! do halo exchanges
    call mpas_dmpar_exch_halo_field(uvReconstructX)
    call mpas_dmpar_exch_halo_field(uvReconstructY)
    call mpas_dmpar_exch_halo_field(uvReconstructZ)

    ! deallocate memory
    call mpas_deallocate_scratch_field(ucReconstructX, .True.)
    call mpas_deallocate_scratch_field(ucReconstructY, .True.)
    call mpas_deallocate_scratch_field(ucReconstructZ, .True.)
    call mpas_deallocate_scratch_field(ucReconstructZonal, .True.)
    call mpas_deallocate_scratch_field(ucReconstructMeridional, .True.)
    call mpas_deallocate_scratch_field(boundaryVertexGlobal, .True.)

  end subroutine ocn_RBFvertex!}}}

!***********************************************************************
!
!  routine ocn_vector_cell_center_to_vertex
!
!> \brief   Interpolate cell center values to vertex values
!> \author  Phillip Wolfram
!> \date    05/27/2014
!> \details
!>  Purpose: interpolate vector field at vertex locations from cell center values
!>           using Barycentric (via Wachspress) interpolation
!>  Input: cell center data and mesh information
!>  Output: interpolated vertex values
!-----------------------------------------------------------------------
  subroutine ocn_vector_cell_center_to_vertex(meshPool, boundaryVertex, boundaryCell, & !{{{
      ucReconstructX, ucReconstructY, ucReconstructZ, &
      uvReconstructX, uvReconstructY, uvReconstructZ)

    implicit none

    ! input variables
    type (mpas_pool_type), pointer, intent(in) :: meshPool !< Input: Mesh information
    real (kind=RKIND), dimension(:,:), pointer, intent(in) :: ucReconstructX, & !< Input: X Cell center values
                                                              ucReconstructY, & !< Input: Y Cell center values
                                                              ucReconstructZ    !< Input: z Cell center values
    integer, dimension(:,:), pointer, intent(in) :: boundaryVertex, boundaryCell !< Input: Boundary flags

    ! output variables
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: uvReconstructX !< Output: Vertex Reconstructed X Velocity Component
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: uvReconstructY !< Output: Vertex Reconstructed Y Velocity Component
    real (kind=RKIND), dimension(:,:), pointer, intent(out) :: uvReconstructZ !< Output: Vertex Reconstructed Z Velocity Component

    ! local variables
    integer, pointer :: nVerticesSolve, nCells, vertexDegree, nVertLevels
    integer :: aVertex, aCell, aLevel
    integer, dimension(:,:), pointer :: cellsOnVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
    real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
    ! temporary arrays needed in the (to be constructed) init procedure
    ! note that lambda is going to be constant for this and could be cached
    real (kind=RKIND), dimension(:), allocatable :: lambda
    real (kind=RKIND), dimension(:,:), allocatable :: pointVertex
    real (kind=RKIND), dimension(3) :: pointInterp
    real (kind=RKIND) :: xp,yp,zp , sumArea, kiteArea
    logical, pointer :: is_periodic
    real(kind=RKIND), pointer :: x_period, y_period

    uvReconstructX = 0.0_RKIND
    uvReconstructY = 0.0_RKIND
    uvReconstructZ = 0.0_RKIND

    call mpas_pool_get_dimension(meshPool, 'vertexDegree', vertexDegree)

    allocate(lambda(vertexDegree), pointVertex(3,vertexDegree))

    ! setup pointers
    call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
    call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
    call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
    call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)

    call mpas_pool_get_array(meshPool, 'xCell', xCell)
    call mpas_pool_get_array(meshPool, 'yCell', yCell)
    call mpas_pool_get_array(meshPool, 'zCell', zCell)

    call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
    call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
    call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

    call mpas_pool_get_array(meshPool, 'kiteAreasOnVertex', kiteAreasOnVertex)

    call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
    call mpas_pool_get_config(meshPool, 'x_period', x_period)
    call mpas_pool_get_config(meshPool, 'y_period', y_period)

    ! loop over all vertices
    do aVertex = 1, nVerticesSolve
      ! could precompute the list as an optimization
      ! really, condition is any boundaryVertex in column greater than 0
      if(any(boundaryVertex(:,aVertex) < 1)) then
        ! get vertex location and cell center locations
        do aCell = 1, vertexDegree
          ! logical could be moved outside of code block as an optimization
          ! (then essentially would have two nearly identical code blocks...)
          if (is_periodic) then
            ! fix periodicity with respect to pointInterp (xVertex)
            pointVertex(1,aCell) = mpas_fix_periodicity(xCell(cellsOnVertex(aCell, aVertex)), xVertex(aVertex), x_period)
            pointVertex(2,aCell) = mpas_fix_periodicity(yCell(cellsOnVertex(aCell, aVertex)), yVertex(aVertex), y_period)
            pointVertex(3,aCell) = zCell(cellsOnVertex(aCell, aVertex))
          else
            pointVertex(1,aCell) = xCell(cellsOnVertex(aCell, aVertex))
            pointVertex(2,aCell) = yCell(cellsOnVertex(aCell, aVertex))
            pointVertex(3,aCell) = zCell(cellsOnVertex(aCell, aVertex))
          end if
        end do
        ! vertex point for reconstruction
        pointInterp(1) = xVertex(aVertex)
        pointInterp(2) = yVertex(aVertex)
        pointInterp(3) = zVertex(aVertex)
        ! get interpolation constants (could be cached / optimized with areaBin)
        lambda = mpas_wachspress_coordinates(vertexDegree, pointVertex , pointInterp, meshPool)
      else
        lambda = 0.0_RKIND
      end if

      do aLevel = 1, nVertLevels
        if(boundaryVertex(aLevel,aVertex) < 1) then
          ! perform interpolation
          uvReconstructX(aLevel,aVertex) = sum(ucReconstructX(aLevel,cellsOnVertex(:,aVertex)) * lambda)
          uvReconstructY(aLevel,aVertex) = sum(ucReconstructY(aLevel,cellsOnVertex(:,aVertex)) * lambda)
          uvReconstructZ(aLevel,aVertex) = sum(ucReconstructZ(aLevel,cellsOnVertex(:,aVertex)) * lambda)
        end if
      end do

      ! need to specify boundary conditions for the vertexes (outside this subroutine)

    end do

    deallocate(lambda, pointVertex)

  end subroutine ocn_vector_cell_center_to_vertex!}}}

!***********************************************************************
!
!  routine ocn_second_order_shapiro_filter_ops
!
!> \brief   Do Ntimes simple shapiro filtering operations, but make
!>          higher order
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: multiple applications of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
  subroutine ocn_second_order_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values
      integer, intent(in) :: Ntimes  ! number of filter applications
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness

      type (field2DReal), pointer :: ucStore, vcStore, wcStore

      call mpas_pool_get_field(scratchPool,'ucX',ucStore)
      call mpas_pool_get_field(scratchPool,'ucY',vcStore)
      call mpas_pool_get_field(scratchPool,'ucZ',wcStore)
      call mpas_allocate_scratch_field(ucStore,.True.)
      call mpas_allocate_scratch_field(vcStore,.True.)
      call mpas_allocate_scratch_field(wcStore,.True.)


      call ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)
      ucStore % array = 2.0_RKIND*ucReconstructX % array
      vcStore % array = 2.0_RKIND*ucReconstructY % array
      wcStore % array = 2.0_RKIND*ucReconstructZ % array
      call ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ)
      ucStore % array = ucStore % array - ucReconstructX % array
      vcStore % array = vcStore % array - ucReconstructY % array
      wcStore % array = wcStore % array - ucReconstructZ % array

      ! move temporary storage into final storage
      ucReconstructX % array = ucStore % array
      ucReconstructY % array = vcStore % array
      ucReconstructZ % array = wcStore % array

      ! deallocate temporary memory
      call mpas_deallocate_scratch_field(ucStore,.True.)
      call mpas_deallocate_scratch_field(vcStore,.True.)
      call mpas_deallocate_scratch_field(wcStore,.True.)

    end subroutine ocn_second_order_shapiro_filter_ops !}}}

!***********************************************************************
!
!  routine ocn_multiple_vector_shapiro_filter_ops
!
!> \brief   Do Ntimes simple shapiro filtering operations
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: multiple applications of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
  subroutine ocn_multiple_vector_shapiro_filter_ops(Ntimes, meshPool, scratchPool, boundaryVertex, boundaryCell, &
        layerThickness, ucReconstructX, ucReconstructY, ucReconstructZ) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      type (field2DInteger), pointer, intent(in) :: boundaryVertex, boundaryCell
      type (field2DReal), pointer, intent(inout) :: ucReconstructX, ucReconstructY, ucReconstructZ ! cell center values
      real (kind=RKIND), dimension(:,:), pointer, intent(in) :: layerThickness
      integer, intent(in) :: Ntimes  ! number of filter applications

      ! local variables
      integer atime

      do atime = 1,Ntimes
        !call ocn_simple_vector_shapiro_filter(meshPool, scratchPool, boundaryVertex, boundaryCell, &
        !  ucReconstructX, ucReconstructY, ucReconstructZ)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructX % array)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructY % array)
        call ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell % array, layerThickness, &
                                                ucReconstructZ % array)

      end do

    end subroutine ocn_multiple_vector_shapiro_filter_ops !}}}

!***********************************************************************
!
!  routine ocn_simple_vector_laplacian_filter
!
!> \brief   Do 1 pass of simple laplacian filter
!> \author  Phillip Wolfram
!> \date    08/01/2014
!> \details
!>  Purpose: one pass of digital shapiro filter (discrete Laplacian)
!>  Input: cell centered data and mesh information
!>  Output: filtered cell values
!-----------------------------------------------------------------------
    subroutine ocn_simple_vector_laplacian_filter(meshPool, scratchPool, boundaryCell, layerThickness, ucReconstruct) !{{{
      implicit none

      type (mpas_pool_type), pointer, intent(in) :: meshPool, scratchPool
      integer, dimension(:,:), pointer, intent(in) :: boundaryCell
      real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: ucReconstruct
      real (kind=RKIND), dimension(:,:), pointer, intent(in) ::  layerThickness

      ! local variables
      type (field2DReal), pointer :: ucTemp
      integer :: aCell, aNeigh, aLevel
      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND) :: volSum, cellVol

      ! allocate scratch memory
      call mpas_pool_get_field(scratchPool, 'ucTemp', ucTemp)
      call mpas_allocate_scratch_field(ucTemp,.True.)

      ! get values from pools
      call mpas_pool_get_dimension(meshPool,'nCellsSolve',nCellsSolve)
      call mpas_pool_get_dimension(meshPool,'nVertLevels',nVertLevels)
      call mpas_pool_get_array(meshPool,'nEdgesOnCell',nEdgesOnCell)
      call mpas_pool_get_array(meshPool,'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool,'areaCell',areaCell)

      ucTemp % array = 0.0_RKIND

      ! perform laplacian filtering
      do aCell = 1, nCellsSolve
        do aLevel = 1, nVertLevels
          volSum = nEdgesOnCell(aCell) * layerThickness(aLevel,aCell) * areaCell(aCell) * (1-boundaryCell(aLevel,aCell))
          ucTemp % array(aLevel, aCell) = ucReconstruct(aLevel, aCell) * volSum
          if (volSum /= 0 ) then
            ! loop over all neighbors
            do aNeigh = 1, nEdgesOnCell(aCell)
              cellVol = layerThickness(aLevel,cellsOnCell(aNeigh,aCell)) * areaCell(cellsOnCell(aNeigh,aCell)) &
                * (1-boundaryCell(aLevel, cellsOnCell(aNeigh,aCell)))
              volSum = volSum + cellVol
              ucTemp % array(aLevel, aCell) = ucTemp % array(aLevel, aCell) + ucReconstruct(aLevel,cellsOnCell(aNeigh,aCell)) &
                                            * cellVol
            end do
            ucTemp % array(aLevel, aCell) = ucTemp % array(aLevel, aCell) / volSum
          end if
        end do
      end do

      ! exchange halo values
      call mpas_dmpar_exch_halo_field(ucTemp)

      ! replace input values with filtered values
      ucReconstruct = ucTemp % array

      ! deallocate scratch memory
      call mpas_deallocate_scratch_field(ucTemp , .True.)

    end subroutine ocn_simple_vector_laplacian_filter !}}}

end module ocn_particle_helpers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
