! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_cdg_advection
!
!> \brief MPAS ocean characteristic discontinuous Galerkin tracer
!>        advection assembly routines
!> \author Dave Lee
!> \date   June 15
!> \details
!>  This module contains the routines to construct swept region fluxes
!>  and cell pre-volumes for the CDG scheme for fields of given 
!>  quadrature and basis orders
!
!-----------------------------------------------------------------------

#define CDG_EDGE_SCALE
#define CDG_QUAD_PT_HORIZ_ADV
#define CDG_QUAD_PT_VERT_ADV

module ocn_cdg_advection

   use mpas_kind_types
   use mpas_sort
   use mpas_hash
   use mpas_timer

   use ocn_constants
   
   use mpas_matrix_operations  
   use mpas_spline_interpolation
   use mpas_planar_projection

   use mpas_geometry_utils

   use ocn_particle_helpers

   use ocn_cdg_types
   use ocn_cdg_triangle
   use ocn_cdg_basis
   use ocn_cdg_characteristic
   use ocn_cdg_limiter

   implicit none
   private
   save

   public :: ocn_cdg_advection_integrate,            &
             ocn_cdg_advection_init_c_kj,            &
             ocn_cdg_advection_swept_regions,        &
             ocn_cdg_advection_swept_regions_planar, &
             ocn_cdg_advection_prev_flux,            &
             ocn_cdg_advection_prev_flux_planar,     &
             ocn_cdg_advection_mass_matrix,          &
             ocn_cdg_advection_mass_matrix_lumped,   &
             ocn_cdg_advection_error,                &
             ocn_cdg_advection_intersection,         &
             ocn_cdg_advection_intersection_plane,   &
             ocn_cdg_advection_edge_scaling

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_advection_init_c_kj
!
!> \brief MPAS ocean CDG advection initialization routine
!> \author Dave Lee
!> \date   October 2015
!> \details
!>  This routine initializes the basis function coefficients
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_advection_init_c_kj(meshPool, statePool, cdgPool, quad, basis, c_kj)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real(kind=RKIND), dimension(:,:,:,:), intent(out) :: c_kj

      integer, pointer :: nCells, nVertLevels, nTracers
      integer, dimension(:), pointer :: nEdgesOnCell
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness

      integer :: iCell, iEdge, iLevel, iQuad, jQuad, iBasis, jBasis, iTracer
      real (kind=RKIND) :: ansA, phi_i, phi_j, xl, yl, zl, wt

      integer, dimension(:), allocatable :: indx
      real (kind=RKIND), dimension(:), allocatable :: rhs
      real (kind=RKIND), dimension(:,:), allocatable :: mat, matInv

      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(indx(basis % nf))
      allocate(rhs(basis % nf))
      allocate(mat(basis % nf, basis % nf), matInv(basis % nf, basis % nf))

      do iCell = 1, nCells
         do iLevel = 1, nVertLevels
!           mat(:,:) = 0.0_RKIND
!           matInv(:,:) = 0.0_RKIND
!           do iEdge = 1, nEdgesOnCell(iCell)
!              do iQuad = 1, quad % nh
!                 xl = quad % xQuadL(iQuad,iEdge,iCell)
!                 yl = quad % yQuadL(iQuad,iEdge,iCell)

!                 do jQuad = 1, quad % nv
!                    zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
!                    wt = triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)

!                    do iBasis = 1, basis % nf
!                       phi_i = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(iBasis,iLevel,iCell), basis % mean_cur(iBasis,iLevel,iCell), iBasis)
!                       do jBasis = 1, basis % nf
!                          phi_j = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(jBasis,iLevel,iCell), basis % mean_cur(jBasis,iLevel,iCell), jBasis)
!                          mat(iBasis,jBasis) = mat(iBasis,jBasis) + wt*phi_i*phi_j
!                       end do ! jBasis
!                    end do ! iBasis
!                 end do ! jQuad
!              end do ! iQuad
!           end do ! iEdge
!           call mpas_migs(mat(:,:), basis % nf, matInv(:,:), indx)

           do iTracer = 1, nTracers
!              rhs = 0.0_RKIND
!              do iEdge = 1, nEdgesOnCell(iCell)
!                 do iQuad = 1, quad % nh
!                    xl = quad % xQuadL(iQuad,iEdge,iCell)
!                    yl = quad % yQuadL(iQuad,iEdge,iCell)

!                    do jQuad = 1, quad % nv
!                       wt = triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)
!                       zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)

!                       !ansA = 5.0_RKIND + 12.5_RKIND*(1.0_RKIND + tanh(0.0002_RKIND*(yl + yCell(iCell) - 0.5_RKIND*(577.350269189626+65673.5931203199))))
!                       !ansA = 5.0_RKIND + 12.5_RKIND*(1.0_RKIND + tanh(0.0002_RKIND*(yl + yCell(iCell) - 0.5_RKIND*(577.350269189626+77797.9487733021))))
!                       !ansA = 5.0_RKIND + 12.5_RKIND*(1.0_RKIND + tanh(0.0002_RKIND*(yl + yCell(iCell) - 38600.0_RKIND)))
                       ansA = 10.0_RKIND + 5.0_RKIND*(1.0_RKIND + tanh(0.0004_RKIND*(yCell(iCell) - 21155.0_RKIND)))

!                       do iBasis = 1, basis % nf
!                          phi_i = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(iBasis,iLevel,iCell), basis % mean_cur(iBasis,iLevel,iCell), iBasis)
!                          rhs(iBasis) = rhs(iBasis) + wt*ansA*phi_i
!                       end do ! iBasis
!                    end do ! jQuad
!                 end do ! iQuad
!              end do ! iEdge

!              ! Initialize the basis functions for this cell
!              c_kj(:,iTracer,iLevel,iCell) = matmul(matInv(:,:),rhs(:))
c_kj(:,iTracer,iLevel,iCell) = 0.0_RKIND
c_kj(1,iTracer,iLevel,iCell) = ansA
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      deallocate(indx)
      deallocate(rhs)
      deallocate(mat, matInv)

   end subroutine ocn_cdg_advection_init_c_kj!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_advection_swept_regions
!
!> \brief MPAS ocean CDG advection swept regions
!> \author Dave Lee
!> \date   September 2015
!> \details
!   Calculate the fluxes across edges
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_advection_swept_regions(meshPool, cdgPool, quad, basis, dt, fluxes, c_kj)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: c_kj
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: fluxes

      integer, pointer :: nCells, nEdges, nVertices, nVertLevels, nTracers
      integer, dimension(:), pointer :: nEdgesOnCell, nIncCells, maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, verticesOnCell, edgesOnCell, cellsOnCell, incCells
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xEdge, yEdge, zEdge, xVertex, yVertex, zVertex, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray, normalVelocity, edgeNormVecs, wTopOld, wTopMid
      real (kind=RKIND), dimension(:,:), pointer :: uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, &
                                                    layerThicknessOld, layerThicknessMid, edgeScale
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC, departurePoints, mean

      logical :: isBowtie
      integer :: iVertex, jVertex, iEdge, jEdge, iCell, jCell, kCell, iQuad, jQuad, iLevel, iBasis, nIntEdges, nTris, iTri, fCell, tCell, iTracer
      real (kind=RKIND) :: basisInto, basisFrom, triArea, zQuad, wt
      real (kind=RKIND), dimension(3) :: dx, eDep
      real (kind=RKIND), dimension(2) :: qInto, qFrom
      real (kind=RKIND), dimension(3,3,20) :: tris
      real (kind=RKIND), dimension(3,3) :: triLocal
      real (kind=RKIND), dimension(3,4) :: preImage
      real (kind=RKIND), dimension(3,20) :: intCell
      logical, dimension(:,:), allocatable :: edgeDone
      real (kind=RKIND), dimension(:), allocatable :: xQuadL, yQuadL, xQuadG, yQuadG, zQuadG, hq, qzf
      real (kind=RKIND), dimension(:,:), allocatable :: qFinalG, zFullOld, zFullMid

      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)

      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
      call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
      call mpas_pool_get_array(cdgPool, 'normVelocitySum', normalVelocity)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', wTopOld)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', wTopMid)

      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityMid', wVertexMid)
      call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)
      call mpas_pool_get_array(cdgPool, 'departurePoints', departurePoints)
      call mpas_pool_get_array(cdgPool, 'nIncCells', nIncCells)
      call mpas_pool_get_array(cdgPool, 'incCells', incCells)
      call mpas_pool_get_array(cdgPool, 'edgeScale', edgeScale)

      allocate(qFinalG(3,quad % nh))
      allocate(xQuadL(quad % nh), yQuadL(quad % nh), xQuadG(quad % nh), yQuadG(quad % nh), zQuadG(quad % nh))
      allocate(edgeDone(nVertLevels,nEdges))
      allocate(hq(nTracers))
      allocate(qzf(quad % nv))
      allocate(zFullOld(nVertLevels+1,nCells), zFullMid(nVertLevels+1,nCells))

      edgeDone(:,:) = .false.
      mean => basis % mean_cur

#ifdef CDG_QUAD_PT_VERT_ADV
      do iCell = 1, nCells
         zFullOld(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         zFullMid(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell), 1, -1
            zFullOld(iLevel,iCell) = zFullOld(iLevel+1,iCell) + layerThicknessOld(iLevel,iCell)
            zFullMid(iLevel,iCell) = zFullMid(iLevel+1,iCell) + layerThicknessMid(iLevel,iCell)
         end do
      end do
#endif

      do iCell = 1, nCells
         do iEdge = 1, nEdgesOnCell(iCell)
            jEdge = edgesOnCell(iEdge,iCell)

            ! Skip over boundary edges
            if (cellsOnEdge(1,jEdge) .gt. nCells .or. cellsOnEdge(1,jEdge) .lt. 1) cycle
            if (cellsOnEdge(2,jEdge) .gt. nCells .or. cellsOnEdge(1,jEdge) .lt. 1) cycle

            iVertex = verticesOnEdge(1,jEdge)
            jVertex = verticesOnEdge(2,jEdge)

            do iLevel = 1, maxLevelEdgeTop(jEdge)
               ! Determine the direction of the flux across the edge
               if (edgeDone(iLevel,jEdge)) cycle
               edgeDone(iLevel,jEdge) = .true.

               if (normalVelocity(iLevel,jEdge) .gt. 0.0_RKIND) then
                  tCell = cellsOnEdge(2,jEdge)
                  fCell = cellsOnEdge(1,jEdge)
               else
                  tCell = cellsOnEdge(1,jEdge)
                  fCell = cellsOnEdge(2,jEdge)
               end if

               ! Construct the edge preimage polygon in global coordinates
               call ocn_cdg_advection_pre_image(verticesOnEdge(:,jEdge), xVertex, yVertex, zVertex, departurePoints, iLevel, preImage)
               isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), preImage(:,1))
               ! Fall back to edge flux if CDG and FCT fluxes come form opposite sides
               if (.not. isBowtie) then
                  eDep(1) = xEdge(jEdge)
                  eDep(2) = yEdge(jEdge)
                  eDep(3) = zEdge(jEdge)
                  call mpas_particle_horizontal_movement(eDep, -dt*normalVelocity(iLevel,jEdge)*edgeNormVecs(:,jEdge), .true.)
                  isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), eDep)
               end if
               if (isBowtie) then
                  preImage(1,1) = xVertex(jVertex)
                  preImage(2,1) = yVertex(jVertex)
                  preImage(3,1) = zVertex(jVertex)
                  preImage(1,2) = xVertex(iVertex)
                  preImage(2,2) = yVertex(iVertex)
                  preImage(3,2) = zVertex(iVertex)
                  dx(:) = -dt*normalVelocity(iLevel,jEdge)*edgeNormVecs(:,jEdge)
                  call mpas_particle_horizontal_movement(preImage(:,1), dx, .true.)
                  call mpas_particle_horizontal_movement(preImage(:,2), dx, .true.)
               end if

               ! Loop over cells incident on this edge 
               do jCell = 1, nIncCells(jEdge)
                  kCell = incCells(jCell,jEdge)
                  if (kCell .gt. nCells) cycle

                  ! Find the intersecting polygons for each cell incident on the edge
                  call ocn_cdg_advection_intersection(radius, nEdgesOnCell(kCell), verticesOnCell(:,kCell), xVertex, yVertex, zVertex, 4, preImage, nIntEdges, intCell)

                  call ocn_cdg_advection_gen_tris(nIntEdges, intCell, nTris, tris)

                  ! Perform quadrature integration (on the plane)
                  do iTri = 1, nTris
                     ! Calculate the quadrature weight
                     triArea = old_ocn_spherical_area(tris(1,1,iTri), tris(2,1,iTri), tris(3,1,iTri), &
                        tris(1,2,iTri), tris(2,2,iTri), tris(3,2,iTri), tris(1,3,iTri), tris(2,3,iTri), tris(3,3,iTri), radius)
                     !triArea = mpas_triangle_signed_area_sphere(tris(:,1,iTri), tris(:,2,iTri), tris(:,3,iTri), radius)
#ifdef CDG_EDGE_SCALE
                     triArea = triArea*edgeScale(iLevel,jEdge)
#endif

                     call ocn_cdg_tri_trans_g2l(transGtoC(:,:,kCell), xCell(kCell), yCell(kCell), zCell(kCell), tris(1,1,iTri), tris(2,1,iTri), tris(3,1,iTri), triLocal(:,1))
                     call ocn_cdg_tri_trans_g2l(transGtoC(:,:,kCell), xCell(kCell), yCell(kCell), zCell(kCell), tris(1,2,iTri), tris(2,2,iTri), tris(3,2,iTri), triLocal(:,2))
                     call ocn_cdg_tri_trans_g2l(transGtoC(:,:,kCell), xCell(kCell), yCell(kCell), zCell(kCell), tris(1,3,iTri), tris(2,3,iTri), tris(3,3,iTri), triLocal(:,3))

                     ! Calculate the quadrature points for the intersecting polygon triangle on the plane
                     call ocn_cdg_tri_quad_pts_3d(quad, triLocal(:,1), triLocal(:,2), triLocal(:,3), xQuadL, yQuadL)

                     ! Map onto sphere using length preserving projection
                     call ocn_cdg_tri_map_lp(xCell(kCell), yCell(kCell), zCell(kCell), quad % nh, xQuadL, yQuadL, transGtoC(:,:,kCell), xQuadG, yQuadG, zQuadG)

                     ! Advect quadrature points forward (in the horiztonal)
                     call ocn_cdg_char_integrate(meshPool, quad % nh, xQuadG, yQuadG, zQuadG, uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, &
                        nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, kCell, areaBArray, dt, qFinalG)

                     ! And the vertical quadrature points also
#ifdef CDG_QUAD_PT_VERT_ADV
                     call ocn_cdg_char_integrate_quad_pts_vertical(quad, maxLevelCell(kCell), iLevel, zFullOld(:,kCell), zFullMid(:,kCell), &
                        wTopOld(:,kCell), wTopMid(:,kCell), dt, zFullOld(iLevel+1,kCell), layerThicknessOld(iLevel,kCell), qzf)
                     ! Transform vertical quadrature point to local coordinate
                     qzf = qzf - 0.5_RKIND*(zFullOld(iLevel,kCell) + zFullOld(iLevel+1,kCell))
#endif

                     do iQuad = 1, quad % nh
                        do jQuad = 1, quad % nv
                           zQuad = layerThicknessOld(iLevel,kCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                           wt = layerThicknessOld(iLevel,kCell)*triArea*quad % qw(iQuad)*quad % qwv(jQuad)
                           ! Evaluate basis function at quadrature point
                           do iTracer = 1, nTracers
                              hq(iTracer) = ocn_cdg_basis_eval_3d(basis, xQuadL(iQuad), yQuadL(iQuad), zQuad, basis % factor(:,iLevel,kCell), mean(:,iLevel,kCell), c_kj(:,iTracer,iLevel,kCell))
                           end do

                           ! Transform global quadrature coord to local planar coords of the cells either side of the edge
                           call ocn_cdg_tri_trans_g2l(transGtoC(:,:,tCell), xCell(tCell), yCell(tCell), zCell(tCell), qFinalG(1,iQuad), qFinalG(2,iQuad), qFinalG(3,iQuad), qInto)
                           call ocn_cdg_tri_trans_g2l(transGtoC(:,:,fCell), xCell(fCell), yCell(fCell), zCell(fCell), qFinalG(1,iQuad), qFinalG(2,iQuad), qFinalG(3,iQuad), qFrom)

                           do iBasis = 1, basis % nf
#ifdef CDG_QUAD_PT_VERT_ADV
                              basisInto = ocn_cdg_basis_eval_3d_i(basis, qInto(1), qInto(2), qzf(jQuad), basis % factor(iBasis,iLevel,tCell), mean(iBasis,iLevel,tCell), iBasis)
                              basisFrom = ocn_cdg_basis_eval_3d_i(basis, qFrom(1), qFrom(2), qzf(jQuad), basis % factor(iBasis,iLevel,fCell), mean(iBasis,iLevel,fCell), iBasis)
#else
                              basisInto = ocn_cdg_basis_eval_3d_i(basis, qInto(1), qInto(2), zQuad, basis % factor(iBasis,iLevel,tCell), mean(iBasis,iLevel,tCell), iBasis)
                              basisFrom = ocn_cdg_basis_eval_3d_i(basis, qFrom(1), qFrom(2), zQuad, basis % factor(iBasis,iLevel,fCell), mean(iBasis,iLevel,fCell), iBasis)
#endif
                              do iTracer = 1, nTracers
                                 fluxes(iBasis,iTracer,iLevel,tCell) = fluxes(iBasis,iTracer,iLevel,tCell) + wt*hq(iTracer)*basisInto
                                 fluxes(iBasis,iTracer,iLevel,fCell) = fluxes(iBasis,iTracer,iLevel,fCell) - wt*hq(iTracer)*basisFrom
                              end do ! iTracer
                           end do ! iBasis
                        end do ! jQuad
                     end do ! iQuad
                  end do ! iTri
               end do ! jCell
            end do ! iLevel
         end do ! iEdge
      end do ! iCell

      deallocate(qFinalG)
      deallocate(xQuadL, yQuadL, xQuadG, yQuadG, zQuadG)
      deallocate(edgeDone)
      deallocate(hq, qzf)
      deallocate(zFullOld, zFullMid)

   end subroutine ocn_cdg_advection_swept_regions!}}}

   subroutine ocn_cdg_advection_swept_regions_planar(meshPool, cdgPool, quad, basis, dt, fluxes, c_kj)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: c_kj
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: fluxes

      integer, pointer :: nCells, nEdges, nVertices, nTracers, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, nIncCells, maxLevelEdgeTop, maxLevelCell
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, verticesOnCell, edgesOnCell, cellsOnCell, cellsOnVertex, incCells
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, xEdge, yEdge, dvEdge, areaCell, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: areaBArray, normVelocityNew, edgeScale, triAreas, &
                                                    uVertexOld, vVertexOld, uVertexMid, vVertexMid
      real (kind=RKIND), dimension(:,:), pointer :: layerThicknessOld, layerThicknessMid, wTopOld, wTopMid, edgeNormVecs
      real (kind=RKIND), dimension(:,:,:), pointer :: departurePoints, mean, departureEdges

      logical :: isBowtie
      integer :: iEdge, jEdge, kEdge, iCell, jCell, kCell, iQuad, jQuad, iLevel, iBasis, nIntEdges, nTris, iTri, fCell, tCell, iTracer, iVertex, jVertex, kVertex
      real (kind=RKIND) :: basisInto, basisFrom, triArea, xDist, xOrig, zQuad, wt
      real (kind=RKIND), dimension(3) :: qInto, qFrom, eDep
      real (kind=RKIND), dimension(3,3,20) :: tris
      real (kind=RKIND), dimension(3,4) :: preImage
      real (kind=RKIND), dimension(3,20) :: intCell, cellCoords
      logical, dimension(:,:), allocatable :: edgeDone
      real (kind=RKIND), dimension(:), allocatable :: xQuad, yQuad, hq, qzf
      real (kind=RKIND), dimension(:,:), allocatable :: zFullOld, zFullMid, qFinal

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
      call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
      call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocityNew)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', wTopOld)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', wTopMid)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'departurePoints', departurePoints)
      call mpas_pool_get_array(cdgPool, 'nIncCells', nIncCells)
      call mpas_pool_get_array(cdgPool, 'incCells', incCells)
      call mpas_pool_get_array(cdgPool, 'edgeScale', edgeScale)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)
      call mpas_pool_get_array(cdgPool, 'departureEdges', departureEdges)

      allocate(qFinal(2,quad % nh))
      allocate(xQuad(quad % nh), yQuad(quad % nh))
      allocate(edgeDone(nVertLevels,nEdges))
      allocate(hq(nTracers))
      allocate(qzf(quad % nv))
      allocate(zFullOld(nVertLevels+1,nCells), zFullMid(nVertLevels+1,nCells))

      edgeDone(:,:) = .false.
      mean => basis % mean_cur

#ifdef CDG_QUAD_PT_VERT_ADV
      do iCell = 1, nCells
         zFullOld(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         zFullMid(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell), 1, -1
            zFullOld(iLevel,iCell) = zFullOld(iLevel+1,iCell) + layerThicknessOld(iLevel,iCell)
            zFullMid(iLevel,iCell) = zFullMid(iLevel+1,iCell) + layerThicknessMid(iLevel,iCell)
         end do
      end do
#endif

      do iVertex = 1, 4
         preImage(3,iVertex) = 0.0_RKIND
      end do 

      do iCell = 1, nCells
         do iEdge = 1, nEdgesOnCell(iCell)
            jEdge = edgesOnCell(iEdge,iCell)

            ! Skip over boundary edges
            if (cellsOnEdge(1,jEdge) .gt. nCells .or. cellsOnEdge(1,jEdge) .lt. 1) cycle
            if (cellsOnEdge(2,jEdge) .gt. nCells .or. cellsOnEdge(1,jEdge) .lt. 1) cycle

            iVertex = verticesOnEdge(1,jEdge)
            jVertex = verticesOnEdge(2,jEdge)

            do iLevel = 1, maxLevelEdgeTop(jEdge)
               if (edgeDone(iLevel,jEdge)) cycle
               edgeDone(iLevel,jEdge) = .true.

               ! Determine the direction of the flux across the edge
               if (normVelocityNew(iLevel,jEdge) .gt. 0.0_RKIND) then
                  tCell = cellsOnEdge(2,jEdge)
                  fCell = cellsOnEdge(1,jEdge)
               else
                  tCell = cellsOnEdge(1,jEdge)
                  fCell = cellsOnEdge(2,jEdge)
               end if

               ! Construct the edge pre-image polygon in global coordinates
               preImage(1,1) = departurePoints(1,iLevel,jVertex)
               preImage(2,1) = departurePoints(2,iLevel,jVertex)
               preImage(1,2) = departurePoints(1,iLevel,iVertex)
               preImage(2,2) = departurePoints(2,iLevel,iVertex)
               preImage(1,3) = xVertex(iVertex)
               preImage(2,3) = yVertex(iVertex)
               preImage(1,4) = xVertex(jVertex)
               preImage(2,4) = yVertex(jVertex)
               isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), preImage(:,1))
               ! Make sure that the CDG flux and the FCT flux are coming from the same side
               if (.not. isBowtie) then
                  eDep(1) = xEdge(jEdge) - dt*normVelocityNew(iLevel,jEdge)*edgeNormVecs(1,jEdge)
                  eDep(2) = yEdge(jEdge) - dt*normVelocityNew(iLevel,jEdge)*edgeNormVecs(2,jEdge)
                  eDep(3) = 0.0_RKIND
                  isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), eDep(:))
               end if
               if (isBowtie) then
                  preImage(1,1) = xVertex(jVertex) - departureEdges(1,iLevel,jEdge)
                  preImage(2,1) = yVertex(jVertex) - departureEdges(2,iLevel,jEdge)
                  preImage(1,2) = xVertex(iVertex) - departureEdges(1,iLevel,jEdge)
                  preImage(2,2) = yVertex(iVertex) - departureEdges(2,iLevel,jEdge)
               end if

               ! Loop over cells incident on this edge
               do jCell = 1, nIncCells(jEdge)
                  kCell = incCells(jCell,jEdge)
                  if (kCell .gt. nCells .or. kCell .lt. 1) cycle
                  if (iLevel .gt. maxLevelCell(kCell)) cycle

                  do kEdge = 1, nEdgesOnCell(kCell)
                     kVertex = verticesOnCell(kEdge,kCell)
                     cellCoords(1,kEdge) = xVertex(kVertex)
                     cellCoords(2,kEdge) = yVertex(kVertex)
                  end do

                  ! Find the intersecting polygons for each cell incident on the edge
                  call ocn_cdg_advection_intersection_plane(nEdgesOnCell(kCell), cellCoords, 4, preImage, nIntEdges, intCell)
                  call ocn_cdg_advection_gen_tris(nIntEdges, intCell, nTris, tris)

                  ! Perform quadrature integration (on the plane)
                  do iTri = 1, nTris
                     ! Calculate the quadrature weight
                     triArea = 0.5_RKIND*abs( (tris(1,1,iTri) - tris(1,3,iTri))*(tris(2,2,iTri) - tris(2,3,iTri)) - &
                                              (tris(2,1,iTri) - tris(2,3,iTri))*(tris(1,2,iTri) - tris(1,3,iTri)) )
#ifdef CDG_EDGE_SCALE
                     triArea = edgeScale(iLevel,jEdge)*triArea
#endif
                     ! Calculate the quadrature points for the intersecting polygon triangle on the plane
                     call ocn_cdg_tri_quad_pts_3d(quad, tris(:,1,iTri), tris(:,2,iTri), tris(:,3,iTri), xQuad, yQuad)

                     ! Advect quadrature points forward
                     call ocn_cdg_char_integrate_planar(meshPool, quad % nh, xQuad, yQuad, uVertexOld, vVertexOld, uVertexMid, vVertexMid, &
                        nCells, xCell, yCell, zCell, xVertex, yVertex, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, &
                        iLevel, kCell, areaBArray, dt, qFinal)

                     ! And the vertical quadrature points also
#ifdef CDG_QUAD_PT_VERT_ADV
                     call ocn_cdg_char_integrate_quad_pts_vertical(quad, maxLevelCell(kCell), iLevel, zFullOld(:,kCell), zFullMid(:,kCell), &
                        wTopOld(:,kCell), wTopMid(:,kCell), dt, zFullOld(iLevel+1,kCell), layerThicknessOld(iLevel,kCell), qzf)
                     ! Transform vertical quadrature point to local coordinate
                     qzf = qzf - 0.5_RKIND*(zFullOld(iLevel,kCell) + zFullOld(iLevel+1,kCell))
#endif

                     do iQuad = 1, quad % nh
                        do jQuad = 1, quad % nv
                           zQuad = layerThicknessOld(iLevel,kCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                           wt = layerThicknessOld(iLevel,kCell)*triArea*quad % qw(iQuad)*quad % qwv(jQuad)

                           ! Evaluate basis function at quadrature point
                           do iTracer = 1, nTracers
                              hq(iTracer) = ocn_cdg_basis_eval_3d(basis, xQuad(iQuad) - xCell(kCell), yQuad(iQuad) - yCell(kCell), zQuad, &
                                 basis % factor(:,iLevel,kCell), mean(:,iLevel,kCell), c_kj(:,iTracer,iLevel,kCell))
                           end do

                           qInto(1) = qFinal(1,iQuad) - xCell(tCell)
                           qInto(2) = qFinal(2,iQuad) - yCell(tCell)
                           qFrom(1) = qFinal(1,iQuad) - xCell(fCell)
                           qFrom(2) = qFinal(2,iQuad) - yCell(fCell)

                           do iBasis = 1, basis % nf
#ifdef CDG_QUAD_PT_VERT_ADV
                              basisInto = ocn_cdg_basis_eval_3d_i(basis, qInto(1), qInto(2), qzf(jQuad), basis % factor(iBasis,iLevel,tCell), mean(iBasis,iLevel,tCell), iBasis)
                              basisFrom = ocn_cdg_basis_eval_3d_i(basis, qFrom(1), qFrom(2), qzf(jQuad), basis % factor(iBasis,iLevel,fCell), mean(iBasis,iLevel,fCell), iBasis)
#else
                              basisInto = ocn_cdg_basis_eval_3d_i(basis, qInto(1), qInto(2), zQuad, basis % factor(iBasis,iLevel,tCell), mean(iBasis,iLevel,tCell), iBasis)
                              basisFrom = ocn_cdg_basis_eval_3d_i(basis, qFrom(1), qFrom(2), zQuad, basis % factor(iBasis,iLevel,fCell), mean(iBasis,iLevel,fCell), iBasis)
#endif
                              do iTracer = 1, nTracers
                                 fluxes(iBasis,iTracer,iLevel,tCell) = fluxes(iBasis,iTracer,iLevel,tCell) + wt*hq(iTracer)*basisInto
                                 fluxes(iBasis,iTracer,iLevel,fCell) = fluxes(iBasis,iTracer,iLevel,fCell) - wt*hq(iTracer)*basisFrom
                              end do ! iTracer
                           end do ! iBasis
                        end do ! jQuad
                     end do ! iQuad
                  end do ! iTri
               end do ! jCell
            end do ! iLevel
         end do ! iEdge
      end do ! iCell

      deallocate(qFinal)
      deallocate(xQuad, yQuad)
      deallocate(edgeDone)
      deallocate(hq, qzf)
      deallocate(zFullOld, zFullMid)

   end subroutine ocn_cdg_advection_swept_regions_planar!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_advection_pre_volume
!
!> \brief MPAS ocean CDG advection cell volumes at the previous time 
!         level
!> \author Dave Lee
!> \date   September 2015
!> \details
!   Calculate the volume of the tracer in the cells at the previous time
!   levels
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_advection_prev_flux(meshPool, statePool, cdgPool, quad, basis, dt, fluxes, c_kj_cur)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: fluxes
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: c_kj_cur

      integer, pointer :: nCells, nTracers, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, xCell, yCell, zCell, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, &
                                                    triAreas, areaBArray, layerThicknessOld, layerThicknessMid, wTopOld, wTopMid
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC, mean

      integer :: iEdge, iCell, iQuad, jQuad, iLevel, iBasis, jBasis, iTracer
      real (kind=RKIND) :: basisInto, xl, yl, zl, wt, phi_i, phi_j
      real (kind=RKIND), dimension(2) :: qFinalL
      real (kind=RKIND), dimension(:), allocatable :: hq, qzf
      real (kind=RKIND), dimension(:,:), allocatable :: qFinalG, zFullOld, zFullMid

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
      call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', wTopOld)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', wTopMid)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityOld', wVertexOld)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'wVertexVelocityMid', wVertexMid)
      call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(qFinalG(3,quad % nh))
      allocate(hq(nTracers))
      allocate(qzf(quad % nv))
      allocate(zFullOld(nVertLevels+1,nCells), zFullMid(nVertLevels+1,nCells))

      mean => basis % mean_cur

#ifdef CDG_QUAD_PT_VERT_ADV
      do iCell = 1, nCells
         zFullOld(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         zFullMid(maxLevelCell(iCell)+1,iCell) = -bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell), 1, -1
            zFullOld(iLevel,iCell) = zFullOld(iLevel+1,iCell) + layerThicknessOld(iLevel,iCell)
            zFullMid(iLevel,iCell) = zFullMid(iLevel+1,iCell) + layerThicknessMid(iLevel,iCell)
         end do
      end do
#endif

      ! Calcualte the tracer volume in each cell at the previous time level
      do iCell = 1, nCells
         do iLevel = 1, maxLevelCell(iCell)
            do iEdge = 1, nEdgesOnCell(iCell)
               call ocn_cdg_char_integrate(meshPool, quad % nh, quad % xQuadG(:,iEdge,iCell), quad % yQuadG(:,iEdge,iCell), quad % zQuadG(:,iEdge,iCell), &
                  uVertexOld, vVertexOld, wVertexOld, uVertexMid, vVertexMid, wVertexMid, nCells, xCell, yCell, zCell, xVertex, yVertex, zVertex, &
                  verticesOnCell, cellsOnCell, nEdgesOnCell, iLevel, iCell, areaBArray, dt, qFinalG)

               ! And the vertical quadrature points also
#ifdef CDG_QUAD_PT_VERT_ADV
               call ocn_cdg_char_integrate_quad_pts_vertical(quad, maxLevelCell(iCell), iLevel, zFullOld(:,iCell), zFullMid(:,iCell), &
                  wTopOld(:,iCell), wTopMid(:,iCell), dt, zFullOld(iLevel+1,iCell), layerThicknessOld(iLevel,iCell), qzf)
               ! Transform vertical quadrature point to local coordinate
               qzf = qzf - 0.5_RKIND*(zFullOld(iLevel,iCell) + zFullOld(iLevel+1,iCell))
#endif

               do iQuad = 1, quad % nh
                  do jQuad = 1, quad % nv
                     zl = layerThicknessOld(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     wt = layerThicknessOld(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)
                     do iTracer = 1, nTracers
                        hq(iTracer) = ocn_cdg_basis_eval_3d(basis, quad % xQuadL(iQuad,iEdge,iCell), quad % yQuadL(iQuad,iEdge,iCell), zl, &
                           basis % factor(:,iLevel,iCell), mean(:,iLevel,iCell), c_kj_cur(:,iTracer,iLevel,iCell))
                     end do
            
                     call ocn_cdg_tri_trans_g2l(transGtoC(:,:,iCell), xCell(iCell), yCell(iCell), zCell(iCell), qFinalG(1,iQuad), qFinalG(2,iQuad), qFinalG(3,iQuad), qFinalL)

                     do iBasis = 1, basis % nf
#ifdef CDG_QUAD_PT_VERT_ADV
                        basisInto = ocn_cdg_basis_eval_3d_i(basis, qFinalL(1), qFinalL(2), qzf(jQuad), basis % factor(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
#else
                        basisInto = ocn_cdg_basis_eval_3d_i(basis, qFinalL(1), qFinalL(2), zl, basis % factor(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
#endif
                        do iTracer = 1, nTracers
                           fluxes(iBasis,iTracer,iLevel,iCell) = fluxes(iBasis,iTracer,iLevel,iCell) + wt*hq(iTracer)*basisInto
                        end do ! iTracer
                     end do ! iBasis
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge
         end do ! iLevel
      end do ! iCell

      deallocate(qFinalG)
      deallocate(hq, qzf)
      deallocate(zFullOld, zFullMid)

   end subroutine ocn_cdg_advection_prev_flux!}}}

   subroutine ocn_cdg_advection_prev_flux_planar(meshPool, cdgPool, quad, basis, dt, fluxes, c_kj_cur)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: fluxes
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: c_kj_cur

      integer, pointer :: nCells, nTracers, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      integer, dimension(:,:), pointer :: verticesOnCell, cellsOnCell, edgesOnCell, cellsOnVertex
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, xCell, yCell, zCell, bottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, areaBArray, layerThicknessOld, layerThicknessMid, wTopOld, wTopMid, &
                                                    uVertexOld, uVertexMid, vVertexOld, vVertexMid

      integer :: iEdge, jEdge, iCell, iQuad, jQuad, iLevel, iBasis, jBasis, iTracer
      real (kind=RKIND), dimension(3) :: xf
      real (kind=RKIND) :: basisInto, weight, zm, zb, zl
      real (kind=RKIND), dimension(:), allocatable :: hq, xQuadInit, yQuadInit, zgf, zFullOld, zFullMid
      real (kind=RKIND), dimension(:,:), allocatable :: qFinal

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
      call mpas_pool_get_array(meshPool, 'wachspressAreaB', areaBArray)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThicknessOld)
      call mpas_pool_get_array(cdgPool, 'layerThicknessMid', layerThicknessMid)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityOld', wTopOld)
      call mpas_pool_get_array(cdgPool, 'vertAleVelocityMid', wTopMid)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityOld', uVertexOld)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityOld', vVertexOld)
      call mpas_pool_get_array(cdgPool, 'uVertexVelocityMid', uVertexMid)
      call mpas_pool_get_array(cdgPool, 'vVertexVelocityMid', vVertexMid)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(qFinal(2,quad % nh))
      allocate(xQuadInit(quad % nh), yQuadInit(quad % nh))
      allocate(hq(nTracers))
      allocate(zgf(quad % nv))
      allocate(zFullOld(nVertLevels+1), zFullMid(nVertLevels+1))

      ! Calcualte the tracer volume in each cell at the previous time level
      do iCell = 1, nCells
#ifdef CDG_QUAD_PT_VERT_ADV
         zFullOld(maxLevelCell(iCell)+1) = -bottomDepth(iCell)
         zFullMid(maxLevelCell(iCell)+1) = -bottomDepth(iCell)
         do iLevel = maxLevelCell(iCell), 1, -1
            zFullOld(iLevel) = zFullOld(iLevel+1) + layerThicknessOld(iLevel,iCell)
            zFullMid(iLevel) = zFullMid(iLevel+1) + layerThicknessMid(iLevel,iCell)
         end do
#endif
         do iLevel = 1, maxLevelCell(iCell)
#ifdef CDG_QUAD_PT_VERT_ADV
            call ocn_cdg_char_integrate_quad_pts_vertical(quad, maxLevelCell(iCell), iLevel, zFullOld, zFullMid, &
               wTopOld(:,iCell), wTopMid(:,iCell), dt, zFullOld(iLevel+1), layerThicknessOld(iLevel,iCell), zgf)
#endif

            ! Assemble the tracer concentration volume integral rhs at the previous time level for each tracer
            do iEdge = 1, nEdgesOnCell(iCell)
               jEdge = edgesOnCell(iEdge,iCell)
               xQuadInit(:) = quad % xQuadL(:,iEdge,iCell) + xCell(iCell)
               yQuadInit(:) = quad % yQuadL(:,iEdge,iCell) + yCell(iCell)

#ifdef CDG_QUAD_PT_HORIZ_ADV
               call ocn_cdg_char_integrate_planar(meshPool, quad % nh, xQuadInit, yQuadInit, uVertexOld, vVertexOld, uVertexMid, vVertexMid, &
                  nCells, xCell, yCell, zCell, xVertex, yVertex, verticesOnCell, cellsOnCell, cellsOnVertex, nEdgesOnCell, iLevel, iCell, &
                  areaBArray, dt, qFinal)
#endif

               do iQuad = 1, quad % nh
                  do jQuad = 1, quad % nv
                     zl = layerThicknessOld(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     ! Normalize by thickness at the appropriate time level
                     weight = layerThicknessOld(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)

                     do iTracer = 1, nTracers
                        hq(iTracer) = ocn_cdg_basis_eval_3d(basis, quad % xQuadL(iQuad,iEdge,iCell), quad % yQuadL(iQuad,iEdge,iCell), zl, &
                           basis % factor(:,iLevel,iCell), basis % mean_cur(:,iLevel,iCell), c_kj_cur(:,iTracer,iLevel,iCell))
                     end do
            
#ifdef CDG_QUAD_PT_HORIZ_ADV
                     xf(1) = qFinal(1,iQuad) - xCell(iCell)
                     xf(2) = qFinal(2,iQuad) - yCell(iCell)
#else
                     xf(1) = quad % xQuadL(iQuad,iEdge,iCell)
                     xf(2) = quad % yQuadL(iQuad,iEdge,iCell)
#endif
#ifdef CDG_QUAD_PT_VERT_ADV
                     xf(3) = zgf(jQuad) - 0.5_RKIND*(zFullOld(iLevel) + zFullOld(iLevel+1))
#else
                     xf(3) = zl
#endif
                     do iBasis = 1, basis % nf
                        basisInto = ocn_cdg_basis_eval_3d_i(basis, xf(1), xf(2), xf(3), basis % factor(iBasis,iLevel,iCell), basis % mean_cur(iBasis,iLevel,iCell), iBasis)
                        do iTracer = 1, nTracers
                           fluxes(iBasis,iTracer,iLevel,iCell) = fluxes(iBasis,iTracer,iLevel,iCell) + weight*hq(iTracer)*basisInto
                        end do
                     end do ! iBasis
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge
         end do ! iLevel
      end do ! iCell

      deallocate(qFinal)
      deallocate(xQuadInit, yQuadInit)
      deallocate(hq)
      deallocate(zgf)
      deallocate(zFullOld, zFullMid)

   end subroutine ocn_cdg_advection_prev_flux_planar!}}}

   subroutine ocn_cdg_advection_mass_matrix(meshPool, statePool, cdgPool, quad, basis, fluxes, c_kj_new)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: fluxes
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: c_kj_new

      integer, pointer :: nCells, nTracers
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: mean

      integer :: iEdge, iCell, iQuad, jQuad, iLevel, iBasis, jBasis, iTracer
      integer, dimension(:), allocatable :: indx
      real (kind=RKIND) :: xl, yl, zl, wt, phi_i, phi_j
      real (kind=RKIND), dimension(:,:), allocatable :: mat, matInv

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(mat(basis % nf,basis % nf), matInv(basis % nf, basis % nf))
      allocate(indx(basis % nf))

      mean => basis % mean_new

      ! Calcualte the tracer volume in each cell at the previous time level
      do iCell = 1, nCells
         do iLevel = 1, maxLevelCell(iCell)
            ! Assemble the mass matrix...
            mat(:,:) = 0.0_RKIND
            matInv(:,:) = 0.0_RKIND

            do iEdge = 1, nEdgesOnCell(iCell)
               do iQuad = 1, quad % nh
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)

                  do jQuad = 1, quad % nv
                     zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     wt = layerThickness(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)

                     do iBasis = 1, basis % nf
                        phi_i = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(iBasis,iLevel,iCell), mean(iBasis,iLevel,iCell), iBasis)
                        do jBasis = 1, basis % nf
                           phi_j = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(jBasis,iLevel,iCell), mean(jBasis,iLevel,iCell), jBasis)
                           mat(iBasis,jBasis) = mat(iBasis,jBasis) + wt*phi_i*phi_j
                        end do ! jBasis
                     end do ! iBasis
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge

            ! ...and take the inverse
            call mpas_migs(mat(:,:), basis % nf, matInv(:,:), indx)

            do iTracer = 1, nTracers
               c_kj_new(:,iTracer,iLevel,iCell) = matmul(matInv(:,:),fluxes(:,iTracer,iLevel,iCell))
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      deallocate(mat, matInv)
      deallocate(indx)

   end subroutine ocn_cdg_advection_mass_matrix!}}}

   subroutine ocn_cdg_advection_mass_matrix_lumped(meshPool, statePool, cdgPool, quad, basis, fluxes, c_kj_new)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(inout) :: fluxes
      real (kind=RKIND), dimension(:,:,:,:), intent(out) :: c_kj_new

      integer, pointer :: nCells, nTracers
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness

      integer :: iEdge, iCell, iQuad, jQuad, iLevel, iBasis, iTracer
      real (kind=RKIND) :: weight, xl, yl, zl, phi_i
      real (kind=RKIND), dimension(:), allocatable :: vec

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_dimension(cdgPool, 'nTracers', nTracers)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      allocate(vec(basis % nf))

      ! Calcualte the tracer volume in each cell at the previous time level
      do iCell = 1, nCells
         do iLevel = 1, maxLevelCell(iCell)
            ! Assemble the mass matrix...
            vec(:) = 0.0_RKIND

            do iEdge = 1, nEdgesOnCell(iCell)
               do iQuad = 1, quad % nh
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)

                  do jQuad = 1, quad % nv
                     zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     ! Normalize by thickness at the appropriate time level
                     weight = layerThickness(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)

                     do iBasis = 1, basis % nf
                        phi_i = ocn_cdg_basis_eval_3d_i(basis, xl, yl, zl, basis % factor(iBasis,iLevel,iCell), basis % mean_new(iBasis,iLevel,iCell), iBasis)
                        vec(iBasis) = vec(iBasis) + weight*phi_i*phi_i
                     end do ! iBasis
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge

            do iTracer = 1, nTracers
               do iBasis = 1, basis % nf
                  c_kj_new(iBasis,iTracer,iLevel,iCell) = fluxes(iBasis,iTracer,iLevel,iCell)/vec(iBasis)
               end do ! iBasis
            end do ! iTracer
         end do ! iLevel
      end do ! iCell

      deallocate(vec)

   end subroutine ocn_cdg_advection_mass_matrix_lumped!}}}

   subroutine ocn_cdg_advection_pre_image(verticesOnEdge, xVertex, yVertex, zVertex, departurePoints, iLevel, preImage)!{{{

      integer, intent(in) :: iLevel
      integer, dimension(:), intent(in) :: verticesOnEdge
      real (kind=RKIND), dimension(:), intent(in) :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:,:), intent(in) :: departurePoints
      real (kind=RKIND), dimension(:,:), intent(out) :: preImage

      integer :: aVertex, bVertex

      aVertex = verticesOnEdge(1)
      bVertex = verticesOnEdge(2)

      preImage(:,1) = departurePoints(:,iLevel,bVertex)
      preImage(:,2) = departurePoints(:,iLevel,aVertex)
      preImage(1,3) = xVertex(aVertex)
      preImage(2,3) = yVertex(aVertex)
      preImage(3,3) = zVertex(aVertex)
      preImage(1,4) = xVertex(bVertex)
      preImage(2,4) = yVertex(bVertex)
      preImage(3,4) = zVertex(bVertex)

   end subroutine ocn_cdg_advection_pre_image!}}}

   subroutine ocn_cdg_advection_intersection(radius, nEdgesOnCell, verticesOnCell, xVertex, yVertex, zVertex, nPreImage, preImage, nIntCell, intCell)!{{{

      integer, intent(in) :: nEdgesOnCell, nPreImage
      integer, dimension(:), intent(in) :: verticesOnCell
      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), dimension(:), intent(in) :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:), intent(in) :: preImage
      integer, intent(out) :: nIntCell
      real (kind=RKIND), dimension(:,:), intent(out) :: intCell

      integer :: iVertex, jVertex, dummy, err
      real (kind=RKIND) :: radInv
      real (kind=RKIND), dimension(3,20) :: targetCell, targetCellNorm, preImageR1
      real (kind=RKIND), dimension(3,20) :: work

      radInv = 1.0_RKIND/radius

      do iVertex = 1, nEdgesOnCell
         jVertex = verticesOnCell(iVertex)
         targetCell(1,iVertex) = xVertex(jVertex)*radInv
         targetCell(2,iVertex) = yVertex(jVertex)*radInv
         targetCell(3,iVertex) = zVertex(jVertex)*radInv
      end do

      ! Normalize the pre image
      do iVertex = 1, nPreImage
         preImageR1(:,iVertex) = preImage(:,iVertex)*radInv
      end do

      ! TODO: pre-compute edge normals and add to registry
      do iVertex = 1, nEdgesOnCell
         jVertex = mod(iVertex,nEdgesOnCell) + 1
         call edge_norm(targetCell(:,iVertex), targetCell(:,jVertex), targetCellNorm(:,iVertex))
      end do

      call clipagainstpolysphere(targetCell, nEdgesOnCell, targetCellNorm, preImageR1, nPreImage, intCell, nIntCell, work, 20, err)

      work = intCell
      dummy = nIntCell
      call remove_duplicates(work, dummy, intCell, nIntCell)

      ! Expand out to physical sphere
      do iVertex = 1, nIntCell
         intCell(:,iVertex) = intCell(:,iVertex)*radius
      end do

   end subroutine ocn_cdg_advection_intersection!}}}

   subroutine ocn_cdg_advection_arc_intersection(radius, p1, p2, p3, p4, intersect, isIntersection)!{{{

      real (kind=RKIND), intent(in) :: radius
      real (kind=RKIND), dimension(:), intent(in) :: p1, p2, p3, p4
      real (kind=RKIND), dimension(:), intent(out) :: intersect
      logical, intent(out) :: isIntersection

      real (kind=RKIND) :: arc1, arc2, arc3, radInv
      real (kind=RKIND), dimension(3) :: p1Norm, p2Norm, p3Norm, p4Norm, e2 

      ! Normalize to unit sphere
      radInv = 1.0_RKIND/radius
      p1Norm(:) = p1(:)*radInv
      p2Norm(:) = p2(:)*radInv
      p3Norm(:) = p3(:)*radInv
      p4Norm(:) = p4(:)*radInv

      ! Compute the normal vector for the second edge
      call edge_norm(p3Norm, p4Norm, e2)

      ! Call the fortran wrapper to Andrew's routine
      call intersectsphere(p1Norm, p2Norm, p3Norm, e2, intersect)

      ! Check that the intersection lies between the end points of the second edge
      !arc1 = ocn_arc_length(p3Norm(1), p3Norm(2), p3Norm(3), p4Norm(1), p4Norm(2), p4Norm(3), radius)
      !arc2 = ocn_arc_length(p3Norm(1), p3Norm(2), p3Norm(3), intersect(1), intersect(2), intersect(3), radius)
      !arc3 = ocn_arc_length(p4Norm(1), p4Norm(2), p4Norm(3), intersect(1), intersect(2), intersect(3), radius)
      arc1 = mpas_arc_length(p3Norm(1), p3Norm(2), p3Norm(3), p4Norm(1), p4Norm(2), p4Norm(3))
      arc2 = mpas_arc_length(p3Norm(1), p3Norm(2), p3Norm(3), intersect(1), intersect(2), intersect(3))
      arc3 = mpas_arc_length(p4Norm(1), p4Norm(2), p4Norm(3), intersect(1), intersect(2), intersect(3))

      if (abs(arc1 - arc2 - arc3) .lt. 1.0E-8_RKIND) then
         ! Expand to surface of the sphere
         intersect(:) = radius*intersect(:)
         isIntersection = .true.
      else
         isIntersection = .false.
      end if

   end subroutine ocn_cdg_advection_arc_intersection!}}}

   subroutine ocn_cdg_advection_intersection_plane(nEdgesOnCell, cellLocal, nPreImage, preImage, nIntCell, intCell)!{{{

      integer, intent(in) :: nEdgesOnCell, nPreImage
      real (kind=RKIND), dimension(:,:), intent(in) :: cellLocal, preImage
      integer, intent(out) :: nIntCell
      real (kind=RKIND), dimension(:,:), intent(out) :: intCell

      integer :: iVertex, jVertex, kVertex, dummy, err
      real (kind=RKIND), dimension(3,20) :: cellNorm, rwrk

      do iVertex = 1, nEdgesOnCell
         jVertex = mod(iVertex,nEdgesOnCell) + 1
         cellNorm(1,iVertex) = cellLocal(2,iVertex) - cellLocal(2,jVertex)
         cellNorm(2,iVertex) = cellLocal(1,jVertex) - cellLocal(1,iVertex)
         cellNorm(3,iVertex) = 0.0_RKIND
      end do

      ! All sides are linear
      call clipagainstpolyplane(cellLocal, nEdgesOnCell, cellNorm, preImage, nPreImage, intCell, nIntCell, rwrk, 20, err)

      rwrk = intCell
      dummy = nIntCell
      call remove_duplicates(rwrk, dummy, intCell, nIntCell)

   end subroutine ocn_cdg_advection_intersection_plane!}}}

   subroutine ocn_cdg_advection_gen_tris(nIntCell, intCell, nTris, tris)!{{{

      integer, intent(in) :: nIntCell
      real (kind=RKIND), dimension(:,:), intent(in) :: intCell
      integer, intent(out) :: nTris
      real (kind=RKIND), dimension(:,:,:), intent(out) :: tris

      integer :: iTri

       if (nIntCell .lt. 3) then
          nTris = 0
          return
       end if

       nTris = nIntCell - 2
       do iTri = 1, nTris
          tris(:,1,iTri) = intCell(:,1)
          tris(:,2,iTri) = intCell(:,iTri+1)
          tris(:,3,iTri) = intCell(:,iTri+2)
       end do

   end subroutine ocn_cdg_advection_gen_tris!}}}

   subroutine ocn_cdg_advection_gen_bowties(iEdge, fCell, tCell, nEdgesOnCell, verticesOnCell, edgesOnCell, xVertex, yVertex, zVertex, iLevel, bInt, depPts, fTri, tTri)!{{{

      integer, intent(in) :: iEdge, fCell, tCell, iLevel
      integer, dimension(:), intent(in) :: nEdgesOnCell
      integer, dimension(:,:), intent(in) :: verticesOnCell, edgesOnCell
      real (kind=RKIND), dimension(:), intent(in) :: xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:), intent(in) :: bInt
      real (kind=RKIND), dimension(:,:,:), intent(in) :: depPts
      real (kind=RKIND), dimension(:,:), intent(out) :: fTri, tTri

      integer :: iVertex, jVertex, aVertex, bVertex
      real (kind=RKIND) :: nsq, absq, ndotab, theta, pi
      real (kind=RKIND), dimension(3) :: va, vb, n, ab

      pi = 3.1415926535897932384626433832795028841971

      do iVertex = 1, nEdgesOnCell(fCell)
         if (edgesOnCell(iVertex,fCell) .eq. iEdge) then
            jVertex = mod(iVertex,nEdgesOnCell(fCell)) + 1
            aVertex = verticesOnCell(iVertex,fCell)
            bVertex = verticesOnCell(jVertex,fCell)
            va(1) = xVertex(aVertex)
            va(2) = yVertex(aVertex)
            va(3) = zVertex(aVertex)
            vb(1) = xVertex(bVertex)
            vb(2) = yVertex(bVertex)
            vb(3) = zVertex(bVertex)
            ! normal vector points into fCell
            call edge_norm(va, vb, n)
            ab = depPts(:,iLevel,bVertex) - depPts(:,iLevel,aVertex)
         end if
      end do

      nsq = n(1)*n(1) + n(2)*n(2) + n(3)*n(3)
      absq = ab(1)*ab(1) + ab(2)*ab(2) + ab(3)*ab(3)
      n = n/sqrt(nsq)
      ab = ab/sqrt(absq)

      ndotab = n(1)*ab(1) + n(2)*ab(2) + n(3)*ab(3)
      theta = acos(ndotab)

      if (abs(theta) .lt. 0.5_RKIND*pi) then
         ! bVertex originates from fCell
         fTri(:,1) = bInt
         fTri(:,2) = vb
         fTri(:,3) = depPts(:,iLevel,bVertex)
         tTri(:,1) = bInt
         tTri(:,2) = va
         tTri(:,3) = depPts(:,iLevel,aVertex)
      else
         ! bVertex originates from tCell
         fTri(:,1) = bInt
         fTri(:,2) = va
         fTri(:,3) = depPts(:,iLevel,aVertex)
         tTri(:,1) = bInt
         tTri(:,2) = vb
         tTri(:,3) = depPts(:,iLevel,bVertex)
      end if

   end subroutine ocn_cdg_advection_gen_bowties!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_advection_integrate
!
!> \brief MPAS ocean CDG advection volume integration
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This routine performs a surface integration by evaluating the 
!>  distribution at the quadrature points as a series expansion of the
!>  basis coefficients
!>  Note that the integration of thickness is implicit in the 
!>  integration of the thickness weighted tracer concentration
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_advection_integrate(meshPool, statePool, cdgPool, step, quad, basis, coeffs, volumeLocal)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      integer, intent(in) :: step
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs
      real (kind=RKIND), intent(out) :: volumeLocal
         
      integer :: iLevel, iCell, iEdge, iQuad, jQuad
      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: nEdgesOnCell, maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: areaCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: mean

      real (kind=RKIND) :: ans, xl, yl, zl, wt
      real (kind=RKIND), dimension(:), allocatable :: weight

      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, step)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      volumeLocal = 0.0_RKIND

      if (step .eq. 1) then
         mean => basis % mean_cur
      else if (step .eq. 2) then
         mean => basis % mean_new
      end if

      do iCell = 1, nCellsSolve
         do iLevel = 1, maxLevelCell(iCell)
            do iEdge = 1, nEdgesOnCell(iCell)
               do iQuad = 1, quad % nh
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)

                  do jQuad = 1, quad % nv
                     zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     wt = layerThickness(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)

                     ans = ocn_cdg_basis_eval_3d(basis, xl, yl, zl, basis % factor(:,iLevel,iCell), mean(:,iLevel,iCell), coeffs(:,1,iLevel,iCell))
                     volumeLocal = volumeLocal + wt*ans
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge
         end do ! iCell
      end do ! iLevel

   end subroutine ocn_cdg_advection_integrate!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_advection_error
!
!> \brief MPAS ocean CDG advection test error
!> \author Dave Lee
!> \date   July 2015
!> \details
!>  This routine integrates the L1 and L2 errors in the difference 
!>  betweem the final solution after a single rotation and the initial
!>  analytic solution
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_advection_error(meshPool, statePool, cdgPool, quad, basis, coeffs, errLocal, normLocal)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(in) :: cdgPool
      type (cdgQuad3D), intent(in) :: quad
      type (cdgBasis3D), intent(in) :: basis
      real (kind=RKIND), dimension(:,:,:,:), intent(in) :: coeffs
      real (kind=RKIND), intent(out) :: errLocal, normLocal

      integer :: iLevel, iCell, iEdge, iQuad, jQuad
      real (kind=RKIND) :: areaTri, ansA, ansN, xl, yl, zl, xg, yg, zg, dz, wt
      real (kind=RKIND) :: pi = 3.1415926535897932384626433832795028841971

      integer, pointer :: nCellsSolve, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: edgesOnCell
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, refBottomDepth
      real (kind=RKIND), dimension(:,:), pointer :: triAreas, layerThickness

      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
      call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

      errLocal = 0.0_RKIND
      normLocal = 0.0_RKIND

      do iCell = 1, nCellsSolve
         do iLevel = 1, nVertLevels
            do iEdge = 1, nEdgesOnCell(iCell)
               do iQuad = 1, quad % nh
                  xl = quad % xQuadL(iQuad,iEdge,iCell)
                  yl = quad % yQuadL(iQuad,iEdge,iCell)
                  xg = xl + xCell(iCell)
                  yg = yl + yCell(iCell)

                  do jQuad = 1, quad % nv
                     zl = layerThickness(iLevel,iCell)*(quad % qzv(jQuad) - 0.5_RKIND)
                     ! TODO: make general for variable layer depths
                     zg = -refBottomDepth(iLevel) + layerThickness(iLevel,iCell)*quad % qzv(jQuad)
                     wt = layerThickness(iLevel,iCell)*triAreas(iEdge,iCell)*quad % qw(iQuad)*quad % qwv(jQuad)
                     !ansA = ocn_cdg_advection_test_func(xg, yg, zg, radius)
                     ansA = 0.0_RKIND
                     ansN = ocn_cdg_basis_eval_3d(basis, xl, yl, zl, basis % factor(:,iLevel,iCell), basis % mean_new(:,iLevel,iCell), coeffs(:,1,iLevel,iCell))

                     errLocal = errLocal + wt*(ansA - ansN)*(ansA - ansN)
                     normLocal = normLocal + wt*ansA*ansA
                  end do ! jQuad
               end do ! iQuad
            end do ! iEdge
         end do ! iLevel
      end do ! iCell

   end subroutine ocn_cdg_advection_error!}}}

   subroutine ocn_cdg_advection_flux_scaling(meshPool, diagnosticsPool, cdgPool, dt, fluxes)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: diagnosticsPool
      type (mpas_pool_type), intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:), intent(out) :: fluxes

      integer, pointer :: nCells, nEdges, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, nIncCells, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, verticesOnCell, incCells
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, xEdge, yEdge, zEdge, dvEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, edgeNormVecs, layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: departurePoints

      logical :: isBowtie
      integer :: iVertex, jVertex, iEdge, jEdge, iCell, jCell, iLevel, nIntEdges, nPreImage, nTris, iTri, kVertex, fCell, tCell
      real (kind=RKIND) :: triArea, dist
      real (kind=RKIND), dimension(3) :: dx, eDep
      real (kind=RKIND), dimension(3,3,20) :: tris
      real (kind=RKIND), dimension(3,4) :: preImage
      real (kind=RKIND), dimension(3,20) :: intCell, cellCoords

      call mpas_pool_get_config(meshPool, 'sphere_radius', radius)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'zEdge', zEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThickness)
      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalVelocity)
      call mpas_pool_get_array(cdgPool, 'departurePoints', departurePoints)
      call mpas_pool_get_array(cdgPool, 'nIncCells', nIncCells)
      call mpas_pool_get_array(cdgPool, 'incCells', incCells)

      do iEdge = 1, nEdges
         ! Skip over boundary edges
         if (cellsOnEdge(1,iEdge) .gt. nCells .or. cellsOnEdge(1,iEdge) .lt. 1) cycle
         if (cellsOnEdge(2,iEdge) .gt. nCells .or. cellsOnEdge(2,iEdge) .lt. 1) cycle

         iVertex = verticesOnEdge(1,iEdge)
         jVertex = verticesOnEdge(2,iEdge)

         do iLevel = 1, maxLevelEdgeTop(iEdge)
            if (normalVelocity(iLevel,iEdge) .gt. 0.0_RKIND) then
               tCell = cellsOnEdge(2,iEdge)
               fCell = cellsOnEdge(1,iEdge)
            else
               tCell = cellsOnEdge(1,iEdge)
               fCell = cellsOnEdge(2,iEdge)
            end if

            ! Construct the edge pre-image polygon in global coordinates
            call ocn_cdg_advection_pre_image(verticesOnEdge(:,iEdge), xVertex, yVertex, zVertex, departurePoints, iLevel, preImage)
            isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), preImage(:,1))
            ! Fall back to edge flux if CDG and FCT fluxes come form opposite sides
            if (.not. isBowtie) then
               eDep(1) = xEdge(iEdge)
               eDep(2) = yEdge(iEdge)
               eDep(3) = zEdge(iEdge)
               call mpas_particle_horizontal_movement(eDep, -dt*normalVelocity(iLevel,iEdge)*edgeNormVecs(:,iEdge), .true.)
               isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), eDep)
            end if
            if (isBowtie) then
               preImage(1,1) = xVertex(jVertex)
               preImage(2,1) = yVertex(jVertex)
               preImage(3,1) = zVertex(jVertex)
               preImage(1,2) = xVertex(iVertex)
               preImage(2,2) = yVertex(iVertex)
               preImage(3,2) = zVertex(iVertex)
               dx(:) = -dt*normalVelocity(iLevel,iEdge)*edgeNormVecs(:,iEdge)
               call mpas_particle_horizontal_movement(preImage(:,1), dx, .true.)
               call mpas_particle_horizontal_movement(preImage(:,2), dx, .true.)
            end if

            ! Loop over cells incident on this edge 
            do iCell = 1, nIncCells(iEdge)
               jCell = incCells(iCell,iEdge)
               if (jCell .gt. nCells .or. jCell .lt. 1) cycle

               ! Find the intersecting polygons for each cell incident on the edge
               do jEdge = 1, nEdgesOnCell(jCell)
                  kVertex = verticesOnCell(jEdge,jCell)
                  cellCoords(1,jEdge) = xVertex(kVertex)
                  cellCoords(2,jEdge) = yVertex(kVertex)
                  cellCoords(3,jEdge) = zVertex(kVertex)
               end do
               call ocn_cdg_advection_intersection(radius, nEdgesOnCell(jCell), verticesOnCell(:,jCell), xVertex, yVertex, zVertex, 4, preImage, nIntEdges, intCell)
               call ocn_cdg_advection_gen_tris(nIntEdges, intCell, nTris, tris)

               do iTri = 1, nTris
                  triArea = old_ocn_spherical_area(tris(1,1,iTri), tris(2,1,iTri), tris(3,1,iTri), &
                     tris(1,2,iTri), tris(2,2,iTri), tris(3,2,iTri), tris(1,3,iTri), tris(2,3,iTri), tris(3,3,iTri), radius)
                  !triArea = mpas_triangle_signed_area_sphere(tris(:,1,iTri), tris(:,2,iTri), tris(:,3,iTri), radius)

                  fluxes(iLevel,iEdge) = fluxes(iLevel,iEdge) + triArea*layerThickness(iLevel,jCell)
               end do ! iTri
             end do ! iCell
          end do ! iLevel
       end do ! iEdge

   end subroutine ocn_cdg_advection_flux_scaling!}}}

   subroutine ocn_cdg_advection_flux_scaling_planar(meshPool, cdgPool, dt, fluxes)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt
      real (kind=RKIND), dimension(:,:), intent(out) :: fluxes

      integer, pointer :: nCells, nEdges, nVertLevels
      integer, dimension(:), pointer :: nEdgesOnCell, nIncCells, maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge, verticesOnEdge, verticesOnCell, edgesOnCell, cellsOnCell, cellsOnVertex, incCells
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, xVertex, yVertex, xEdge, yEdge
      real (kind=RKIND), dimension(:,:), pointer :: normVelocityNew, edgeNormVecs, layerThickness
      real (kind=RKIND), dimension(:,:,:), pointer :: depPts, depEdges

      logical :: isBowtie
      integer :: iEdge, jEdge, iCell, jCell, iLevel, nIntEdges, nTris, iTri, iVertex, jVertex, kVertex, fCell, tCell
      real (kind=RKIND) :: triArea
      real (kind=RKIND), dimension(3) :: eDep
      real (kind=RKIND), dimension(3,3,20) :: tris
      real (kind=RKIND), dimension(3,20) :: preImage, intCell, cellCoords

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'xEdge', xEdge)
      call mpas_pool_get_array(meshPool, 'yEdge', yEdge)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(cdgPool, 'departurePoints', depPts)
      call mpas_pool_get_array(cdgPool, 'nIncCells', nIncCells)
      call mpas_pool_get_array(cdgPool, 'incCells', incCells)
      call mpas_pool_get_array(cdgPool, 'layerThicknessOld', layerThickness)
      call mpas_pool_get_array(cdgPool, 'normVelocitySum', normVelocityNew)
      call mpas_pool_get_array(cdgPool, 'departureEdges', depEdges)

      do iVertex = 1, 4
         preImage(3,iVertex) = 0.0_RKIND
      end do 

      do iEdge = 1, nEdges
         ! Skip over boundary edges
         if (cellsOnEdge(1,iEdge) .gt. nCells .or. cellsOnEdge(1,iEdge) .lt. 1) cycle
         if (cellsOnEdge(2,iEdge) .gt. nCells .or. cellsOnEdge(2,iEdge) .lt. 1) cycle

         iVertex = verticesOnEdge(1,iEdge)
         jVertex = verticesOnEdge(2,iEdge)

         do iLevel = 1, maxLevelEdgeTop(iEdge)
            if (normVelocityNew(iLevel,iEdge) .gt. 0.0_RKIND) then
            !deDotN = depEdges(1,iLevel,iEdge)*edgeNormVecs(1,iEdge) + depEdges(2,iLevel,iEdge)*edgeNormVecs(2,iEdge)
            !if (deDotN .gt. 0.0_RKIND) then
               tCell = cellsOnEdge(2,iEdge)
               fCell = cellsOnEdge(1,iEdge)
            else
               tCell = cellsOnEdge(1,iEdge)
               fCell = cellsOnEdge(2,iEdge)
            end if

            ! Construct the edge pre-image polygon in global coordinates
            preImage(1,1) = depPts(1,iLevel,jVertex)
            preImage(2,1) = depPts(2,iLevel,jVertex)
            preImage(1,2) = depPts(1,iLevel,iVertex)
            preImage(2,2) = depPts(2,iLevel,iVertex)
            preImage(1,3) = xVertex(iVertex)
            preImage(2,3) = yVertex(iVertex)
            preImage(1,4) = xVertex(jVertex)
            preImage(2,4) = yVertex(jVertex)
            isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), preImage(:,1))
            ! Make sure that the CDG flux and the FCT flux are coming from the same side
            if (.not. isBowtie) then
               eDep(1) = xEdge(iEdge) - dt*normVelocityNew(iLevel,iEdge)*edgeNormVecs(1,iEdge)
               eDep(2) = yEdge(iEdge) - dt*normVelocityNew(iLevel,iEdge)*edgeNormVecs(2,iEdge)
               eDep(3) = 0.0_RKIND
               isBowtie = ocn_cdg_advection_bowtie_test(preImage(:,3), preImage(:,4), preImage(:,2), eDep(:))
            end if
            if (isBowtie) then
               preImage(1,1) = xVertex(jVertex) - depEdges(1,iLevel,iEdge)
               preImage(2,1) = yVertex(jVertex) - depEdges(2,iLevel,iEdge)
               preImage(1,2) = xVertex(iVertex) - depEdges(1,iLevel,iEdge)
               preImage(2,2) = yVertex(iVertex) - depEdges(2,iLevel,iEdge)
            end if

            ! Loop over cells incident on this edge 
            do iCell = 1, nIncCells(iEdge)
               jCell = incCells(iCell,iEdge)
               if (jCell .gt. nCells .or. jCell .lt. 1) cycle

               ! Find the intersecting polygons for each cell incident on the edge
               do jEdge = 1, nEdgesOnCell(jCell)
                  kVertex = verticesOnCell(jEdge,jCell)
                  cellCoords(1,jEdge) = xVertex(kVertex)
                  cellCoords(2,jEdge) = yVertex(kVertex)
               end do
               call ocn_cdg_advection_intersection_plane(nEdgesOnCell(jCell), cellCoords, 4, preImage, nIntEdges, intCell)
               call ocn_cdg_advection_gen_tris(nIntEdges, intCell, nTris, tris)

               do iTri = 1, nTris
                  triArea = 0.5_RKIND*abs( (tris(1,1,iTri) - tris(1,3,iTri))*(tris(2,2,iTri) - tris(2,3,iTri)) - &
                                           (tris(2,1,iTri) - tris(2,3,iTri))*(tris(1,2,iTri) - tris(1,3,iTri)) )

                  fluxes(iLevel,iEdge) = fluxes(iLevel,iEdge) + triArea*layerThickness(iLevel,jCell)
               end do ! iTri
             end do ! iCell
          end do ! iLevel
       end do ! iEdge

   end subroutine ocn_cdg_advection_flux_scaling_planar!}}}

! Rescale the edge fluxes for consistency with the dynamics
! This is done by evaluating the thickness fluxes using both the 
! standard and CDG schemes and then rescaling the CDG thickness
! weighted tracer concentrations by h_{std}/h_{cdg}

   subroutine ocn_cdg_advection_edge_scaling(meshPool, diagnosticsPool, cdgPool, dt)!{{{

      type(mpas_pool_type), intent(in) :: meshPool
      type(mpas_pool_type), intent(in) :: diagnosticsPool
      type(mpas_pool_type), intent(in) :: cdgPool
      real (kind=RKIND), intent(in) :: dt

      logical, pointer :: onSphere
      integer, pointer :: nCells, nEdges, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop
      integer, dimension(:,:), pointer :: cellsOnEdge
      real (kind=RKIND), dimension(:), pointer :: dvEdge
      real (kind=RKIND), dimension(:,:), pointer :: normalVelocity, edgeScale, contFluxEdge, edgeNormVecs
      real (kind=RKIND), dimension(:,:,:), pointer :: depEdges

      integer :: iEdge, iLevel, iCell, jCell
      real (kind=RKIND), dimension(:,:), allocatable :: flux_cdg

#ifndef CDG_EDGE_SCALE
      return
#endif

      call mpas_pool_get_config(meshPool, 'on_a_sphere', onSphere)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'edgeNormalVectors', edgeNormVecs)
      call mpas_pool_get_array(cdgPool, 'edgeScale', edgeScale)
      call mpas_pool_get_array(cdgPool, 'contFluxEdge', contFluxEdge)
      call mpas_pool_get_array(cdgPool, 'normVelocitySum', normalVelocity)
      call mpas_pool_get_array(cdgPool, 'departureEdges', depEdges)
      
      allocate(flux_cdg(nVertLevels,nEdges))
      flux_cdg(:,:) = 0.0_RKIND

      if (onSphere) then
         call ocn_cdg_advection_flux_scaling(meshPool, diagnosticsPool, cdgPool, dt, flux_cdg)
      else
         call ocn_cdg_advection_flux_scaling_planar(meshPool, cdgPool, dt, flux_cdg)
      end if
     
      do iEdge = 1, nEdges
         do iLevel = 1, maxLevelEdgeTop(iEdge)
            edgeScale(iLevel,iEdge) = 0.0_RKIND
            iCell = cellsOnEdge(1,iEdge)
            jCell = cellsOnEdge(2,iEdge)
            if (iCell .gt. nCells .or. iCell .lt. 1) cycle
            if (jCell .gt. nCells .or. jCell .lt. 1) cycle

            !deDotN = depEdges(1,iLevel,iEdge)*edgeNormVecs(1,iEdge) + depEdges(2,iLevel,iEdge)*edgeNormVecs(2,iEdge)

            ! continuity and cdg fluxes in opposite directions
!            if (normalVelocity(iLevel,iEdge)*contFluxEdge(iLevel,iEdge) .lt. 0.0_RKIND) then
!               edgeScale(iLevel,iEdge) = 0.0_RKIND
            !if (deDotN*contFluxEdge(iLevel,iEdge) .lt. 0.0_RKIND) then
            !   edgeScale(iLevel,iEdge) = -1.0_RKIND*abs(contFluxEdge(iLevel,iEdge)/flux_cdg(iLevel,iEdge))
!            else if (flux_cdg(iLevel,iEdge) .gt. 1.0E-12_RKIND) then
!               edgeScale(iLevel,iEdge) = abs(contFluxEdge(iLevel,iEdge)/flux_cdg(iLevel,iEdge))
!            end if
            
            if (flux_cdg(iLevel,iEdge) .gt. 1.0E-12_RKIND) then
               edgeScale(iLevel,iEdge) = abs(contFluxEdge(iLevel,iEdge)/flux_cdg(iLevel,iEdge))
            end if
         end do
      end do 

      deallocate(flux_cdg)

   end subroutine ocn_cdg_advection_edge_scaling!}}}

   function ocn_cdg_advection_bowtie_test(ei, ej, pi, pj) result (isBowtie)!{{{

      real (kind=RKIND), dimension(:), intent(in) :: ei, ej, pi, pj
      logical :: isBowtie

      real (kind=RKIND), dimension(3) :: a, b, c, ab, ac

      ! edge vector
      a(:) = ej(:) - ei(:)
      ! departure point i vector
      b(:) = pi(:) - ei(:)
      ! departure point j vector
      c(:) = pj(:) - ei(:)

      ab(1) = a(2)*b(3) - a(3)*b(2)
      ab(2) = a(3)*b(1) - a(1)*b(3)
      ab(3) = a(1)*b(2) - a(2)*b(1)

      ac(1) = a(2)*c(3) - a(3)*c(2)
      ac(2) = a(3)*c(1) - a(1)*c(3)
      ac(3) = a(1)*c(2) - a(2)*c(1)

      if (ab(1)*ac(1) + ab(2)*ac(2) + ab(3)*ac(3) .gt. 0.0_RKIND) then
         isBowtie = .false.
      else
         isBowtie = .true.
      end if

   end function ocn_cdg_advection_bowtie_test!}}}

   subroutine edge_norm( vi, vf, n )!{{{

      real (kind=RKIND), dimension(3), intent(in) :: vi, vf
      real (kind=RKIND), dimension(3), intent(out) :: n

      real (kind=RKIND), dimension(3) :: t, s

      t = vf - vi
      ! Note: s vector not normalized to sphere radius
      s = 0.5_RKIND*(vi + vf)

      n(1) = s(2)*t(3) - s(3)*t(2)
      n(2) = s(3)*t(1) - s(1)*t(3)
      n(3) = s(1)*t(2) - s(2)*t(1)

   end subroutine edge_norm!}}}

   subroutine remove_duplicates( vi, ni, vf, nf )!{{{

      real (kind=RKIND), dimension(:,:), intent(in) :: vi
      integer, intent(in) :: ni
      real (kind=RKIND), dimension(:,:), intent(out) :: vf
      integer, intent(out) :: nf

      logical :: add
      integer :: i, j
      real (kind=RKIND) :: rsq
      real (kind=RKIND), dimension(3) :: r

      nf = 0

      do i = 1, ni
         add = .true.
         do j = i + 1, ni
            r = vi(:,i) - vi(:,j)
            rsq = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
            if (rsq .lt. 1.0E-12_RKIND) add = .false.
            !if (rsq .lt. 1.0E-16_RKIND) add = .false.
         end do
         if (add) then
            nf = nf + 1
            vf(:,nf) = vi(:,i)
         end if
      end do

   end subroutine remove_duplicates!}}}

end module ocn_cdg_advection

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
