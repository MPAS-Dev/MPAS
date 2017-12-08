! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_cdg_triangle
!
!> \brief MPAS ocean triangle quadrature routines
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This module contains the routines for instantiating quadrature coords
!>  and quad_ws in local cell space and transforming from global to cell
!>  space.
!
!-----------------------------------------------------------------------

module ocn_cdg_triangle

   use mpas_kind_types
   use mpas_sort
   use mpas_hash
   use mpas_dmpar

   use ocn_constants

   use mpas_planar_projection

   use mpas_geometry_utils

   use ocn_particle_helpers
   use ocn_cdg_types

   implicit none
   private
   save

   public :: ocn_cdg_init_areas,          &
             ocn_cdg_tri_init,            &
             ocn_cdg_tri_init_3d,         &
             ocn_cdg_tri_init_quad_pts,   &
             ocn_cdg_tri_trans_g2l,       &
             ocn_cdg_tri_quad_pts,        &
             ocn_cdg_tri_quad_pts_3d,      &
             ocn_cdg_tri_integrate_full,  &
             ocn_cdg_tri_test,            &
             ocn_cdg_tri_map_lp

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_init_areas
!
!> \brief MPAS ocean CDG triangle quadrature initialization
!> \author Dave Lee
!> \date   January 2016
!> \details
!>  This routine initiates the triangle areas for each cell and the 
!   global to cell transformation matrices
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_init_areas(domain, err)!{{{

      type (domain_type), intent(inout) :: domain !<Input/Output: domain information
      integer, intent(out) :: err !< Output: the error flag

      type (block_type), pointer :: block
      type (dm_info), pointer :: dminfo

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: cdgPool

      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex, areaCell
      real (kind=RKIND), dimension(:,:), pointer :: triAreas
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC, transCtoG

      integer :: iCell, iEdge, jEdge, iVertex, jVertex

      err = 0

      block => domain % blockList
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'cdg', cdgPool)

         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'zCell', zCell)
         call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

         call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

         if (domain % on_a_sphere) then
            call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

            call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)
            call mpas_pool_get_array(cdgPool, 'transCellToGlobal', transCtoG)

            call define_local_to_global_transformations(xCell, yCell, zCell, nCells, transCtoG, transGtoC)
         end if

         do iCell = 1, nCells
            do iEdge = 1, nEdgesOnCell(iCell)
               jEdge = mod(iEdge,nEdgesOnCell(iCell))+1
               iVertex = verticesOnCell(iEdge,iCell)
               jVertex = verticesOnCell(jEdge,iCell)
               if (domain % on_a_sphere) then
                  triAreas(iEdge,iCell) = old_ocn_spherical_area(xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), &
                     xVertex(jVertex), yVertex(jVertex), zVertex(jVertex), xCell(iCell), yCell(iCell), zCell(iCell), radius)
               else
                  triAreas(iEdge,iCell) = 0.5_RKIND*abs( (xVertex(iVertex) - xCell(iCell))*(yVertex(jVertex) - yCell(iCell)) - &
                                                         (yVertex(iVertex) - yCell(iCell))*(xVertex(jVertex) - xCell(iCell)) )
               end if
            end do
         end do

         block => block % next
      end do

   end subroutine ocn_cdg_init_areas!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_init
!
!> \brief MPAS ocean CDG triangle quadrature initialization
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This routine initiates the triangle quadrature of a given order
!
!-----------------------------------------------------------------------

  subroutine ocn_cdg_tri_init(domain, quad, err)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain !<Input/Output: domain information
      type (cdgQuadrature), intent(inout) :: quad
      integer, intent(out) :: err !<Output: error flag

      type (mpas_pool_type), pointer :: cdgPool

      integer :: iQuad

      err = 0

      call mpas_pool_get_subpool(domain % blocklist % structs, 'cdg', cdgPool)

      call mpas_pool_get_config(domain % configs, 'config_cdg_quad_order_horiz', quad % order)
      call mpas_pool_get_config(domain % configs, 'config_cdg_num_quad_points_horiz', quad % nq)

      allocate(quad % qw(quad % nq))
      allocate(quad % qx(quad % nq))
      allocate(quad % qy(quad % nq))
      allocate(quad % qz(quad % nq))

      select case (quad % order)
         case (1)
            quad % qw = (/ 1.0_RKIND /)
            quad % qx = (/ 3.33333333333333333e-01_RKIND /)
            quad % qy = (/ 3.33333333333333333e-01_RKIND /)
         case (2)
            quad % qw = (/ &
               3.33333333333333333333333333333333333333333e-01_RKIND, &
               3.33333333333333333333333333333333333333333e-01_RKIND, &
               3.33333333333333333333333333333333333333333e-01_RKIND /)
            quad % qx = (/ &
               1.66666666666666666666666666666666666666666e-01_RKIND, &
               6.66666666666666666666666666666666666666666e-01_RKIND, &
               1.66666666666666666666666666666666666666666e-01_RKIND /)
            quad % qy = (/ &
               6.66666666666666666666666666666666666666666e-01_RKIND, &
               1.66666666666666666666666666666666666666666e-01_RKIND, &
               1.66666666666666666666666666666666666666666e-01_RKIND /)
         case (3:4)
            quad % qw = (/ &
               1.09951743655321843e-01_RKIND, &
               1.09951743655321857e-01_RKIND, &
               1.09951743655321885e-01_RKIND, &
               2.23381589678011389e-01_RKIND, &
               2.23381589678011527e-01_RKIND, &
               2.23381589678011527e-01_RKIND /)
            quad % qx = (/ &
               9.15762135097704655e-02_RKIND, &
               8.16847572980458514e-01_RKIND, &
               9.15762135097710761e-02_RKIND, &
               1.08103018168070275e-01_RKIND, &
               4.45948490915965612e-01_RKIND, &
               4.45948490915964113e-01_RKIND /)
            quad % qy = (/ &
               8.16847572980458514e-01_RKIND, &
               9.15762135097710761e-02_RKIND, &
               9.15762135097704655e-02_RKIND, &
               4.45948490915964113e-01_RKIND, &
               1.08103018168070275e-01_RKIND, &
               4.45948490915965612e-01_RKIND /)
         case (5)
            quad % qw = (/ &
               1.31356049751916795e-02_RKIND, &
               1.31358306034076201e-02_RKIND, &
               1.37081973800151392e-02_RKIND, &
               1.17419193291163376e-01_RKIND, &
               1.17420611913379477e-01_RKIND, &
               1.24012589655715613e-01_RKIND, &
               1.24015246126072495e-01_RKIND, &
               1.25930230276426303e-01_RKIND, &
               1.25933026682913923e-01_RKIND, &
               2.25289469095714456e-01_RKIND /)
            quad % qx = (/ &
               0.00000000000000000e+00_RKIND, &
               1.00000000000000000e+00_RKIND, &
               0.00000000000000000e+00_RKIND, &
               2.67327353118498978e-01_RKIND, &
               6.72817552946136210e-01_RKIND, &
               6.49236350054349654e-02_RKIND, &
               6.71649853904175198e-01_RKIND, &
               6.54032456800035522e-02_RKIND, &
               2.69376706913982855e-01_RKIND, &
               3.38673850389605513e-01_RKIND /)
            quad % qy = (/ &
               1.00000000000000000e+00_RKIND, &
               0.00000000000000000e+00_RKIND, &
               0.00000000000000000e+00_RKIND, &
               6.72819921871012694e-01_RKIND, &
               2.67328859948191944e-01_RKIND, &
               6.71653011149382917e-01_RKIND, &
               6.49251690028951334e-02_RKIND, &
               2.69378936645285116e-01_RKIND, &
               6.54054874919145490e-02_RKIND, &
               3.38679989302702156e-01_RKIND /)
         case (6)
            quad % qw = (/ &
               3.80680718529555623e-02_RKIND, &
               3.83793553077528410e-02_RKIND, &
               4.62004567445618367e-02_RKIND, &
               5.34675894441989999e-02_RKIND, &
               8.37558269657456833e-02_RKIND, &
               1.01644833025517037e-01_RKIND, &
               1.01861524461366940e-01_RKIND, &
               1.11421831660001677e-01_RKIND, &
               1.12009450262946106e-01_RKIND, &
               1.24787571437558295e-01_RKIND, &
               1.88403488837394911e-01_RKIND /)
            quad % qx = (/ &
               5.72549866774768601e-02_RKIND, &
               8.95362640024579104e-01_RKIND, &
               6.84475748456514044e-01_RKIND, &
               6.87462559150295305e-02_RKIND, &
               6.15676205575839575e-01_RKIND, &
               6.27946141197789465e-01_RKIND, &
               6.29091383418635686e-02_RKIND, &
               6.83782119205099126e-02_RKIND, &
               2.87529458374392255e-01_RKIND, &
               3.28783556413134614e-01_RKIND, &
               3.12290405013644801e-01_RKIND /)
            quad % qy = (/ &
               8.95498146789879490e-01_RKIND, &
               6.18282212503219533e-02_RKIND, &
               2.33437384976827311e-02_RKIND, &
               6.00302757472630025e-02_RKIND, &
               3.33461808341377175e-01_RKIND, &
               1.59189185992151483e-01_RKIND, &
               6.55295093705452469e-01_RKIND, &
               3.09117685428267230e-01_RKIND, &
               6.36426509179620181e-01_RKIND, &
               7.70240056424634223e-02_RKIND, &
               3.52344786445899505e-01_RKIND /)
         case (8)
            quad % qw = (/ &
               2.72303141744348991e-02_RKIND, &
               2.72303141744349199e-02_RKIND, &
               2.72303141744349199e-02_RKIND, &
               2.72303141744349789e-02_RKIND, &
               2.72303141744349789e-02_RKIND, &
               2.72303141744349997e-02_RKIND, &
               3.24584976231980793e-02_RKIND, &
               3.24584976231980793e-02_RKIND, &
               3.24584976231981001e-02_RKIND, &
               9.50916342672845638e-02_RKIND, &
               9.50916342672846193e-02_RKIND, &
               9.50916342672846193e-02_RKIND, &
               1.03217370534718286e-01_RKIND, &
               1.03217370534718314e-01_RKIND, &
               1.03217370534718314e-01_RKIND, &
               1.44315607677787283e-01_RKIND /)
            quad % qx = (/ &
               7.28492392955404355e-01_RKIND, &
               8.39477740995753056e-03_RKIND, &
               2.63112829634638112e-01_RKIND, &
               8.39477740995753056e-03_RKIND, &
               7.28492392955404355e-01_RKIND, &
               2.63112829634638112e-01_RKIND, &
               5.05472283170310122e-02_RKIND, &
               5.05472283170309566e-02_RKIND, &
               8.98905543365938087e-01_RKIND, &
               4.59292588292723236e-01_RKIND, &
               8.14148234145536387e-02_RKIND, &
               4.59292588292723125e-01_RKIND, &
               1.70569307751760324e-01_RKIND, &
               1.70569307751760046e-01_RKIND, &
               6.58861384496479685e-01_RKIND, &
               3.33333333333333370e-01_RKIND /)
            quad % qy = (/ &
               8.39477740995753056e-03_RKIND, &
               2.63112829634638112e-01_RKIND, &
               7.28492392955404355e-01_RKIND, &
               7.28492392955404355e-01_RKIND, &
               2.63112829634638112e-01_RKIND, &
               8.39477740995753056e-03_RKIND, &
               5.05472283170309566e-02_RKIND, &
               8.98905543365938087e-01_RKIND, &
               5.05472283170310122e-02_RKIND, &
               8.14148234145536387e-02_RKIND, &
               4.59292588292723125e-01_RKIND, &
               4.59292588292723236e-01_RKIND, &
               6.58861384496479685e-01_RKIND, &
               1.70569307751760324e-01_RKIND, &
               1.70569307751760046e-01_RKIND, &
               3.33333333333333370e-01_RKIND /)
         case (9)
            quad % qw = (/ &
               1.58676858667487208e-02_RKIND, &
               2.19524732703951786e-02_RKIND, &
               2.40354401213296598e-02_RKIND, &
               2.58522468388647786e-02_RKIND, &
               2.71951393759608216e-02_RKIND, &
               3.02097786027936584e-02_RKIND, &
               3.70093240446550606e-02_RKIND, &
               4.11482921825866571e-02_RKIND, &
               4.26331605467379776e-02_RKIND, &
               4.71413336863812371e-02_RKIND, &
               5.45129844125978591e-02_RKIND, &
               6.26632599630084636e-02_RKIND, &
               6.31379657675310846e-02_RKIND, &
               7.14623133641135444e-02_RKIND, &
               7.51048615652924606e-02_RKIND, &
               7.98259878444318866e-02_RKIND, &
               8.16607475819435963e-02_RKIND, &
               9.37481686311500140e-02_RKIND, &
               1.04838836333477403e-01_RKIND /)
            quad % qx = (/ &
               2.26739052759332704e-01_RKIND, &
               4.77345862087794129e-02_RKIND, &
               2.26577168977105115e-02_RKIND, &
               9.10074385862343016e-01_RKIND, &
               4.41452661673673585e-02_RKIND, &
               4.79944340675050984e-01_RKIND, &
               7.42657808541620557e-01_RKIND, &
               7.43369623518591927e-01_RKIND, &
               2.79454959355581213e-02_RKIND, &
               3.71861932583309532e-02_RKIND, &
               2.22639561442096401e-01_RKIND, &
               1.16082059855864395e-01_RKIND, &
               4.73822270420208358e-01_RKIND, &
               4.77758170054016440e-01_RKIND, &
               6.46387881792721997e-01_RKIND, &
               2.85357695207302253e-01_RKIND, &
               2.04236860041029755e-01_RKIND, &
               1.59370884213907937e-01_RKIND, &
               3.95698265017060125e-01_RKIND /)
            quad % qy = (/ &
               0.00000000000000000e+00_RKIND, &
               9.16183156802148568e-01_RKIND, &
               7.97193825386026345e-01_RKIND, &
               4.44666861644595901e-02_RKIND, &
               4.81588383854628099e-02_RKIND, &
               5.01294615157430568e-01_RKIND, &
               3.03405081749971196e-02_RKIND, &
               2.22245578824042445e-01_RKIND, &
               5.25527023486726308e-01_RKIND, &
               2.39263537482135413e-01_RKIND, &
               7.29063709376736702e-01_RKIND, &
               6.62507673462198188e-01_RKIND, &
               4.60334709656892230e-02_RKIND, &
               4.01038691325781238e-01_RKIND, &
               1.65342747538830548e-01_RKIND, &
               4.92973630851354261e-01_RKIND, &
               1.19056565447230756e-01_RKIND, &
               3.66261159763432431e-01_RKIND, &
               2.27511600022304139e-01_RKIND /)
         end select

      do iQuad = 1, quad % nq
         quad % qz(iQuad) = 1.0_RKIND - quad % qx(iQuad) - quad % qy(iQuad)
      end do

   end subroutine ocn_cdg_tri_init!}}}

   subroutine ocn_cdg_tri_init_3d(domain, quad, err)!{{{

      implicit none

      type (domain_type), intent(inout) :: domain !<Input/Output: domain information
      type (cdgQuad3D), intent(inout) :: quad
      integer, intent(out) :: err !<Output: error flag

      type (mpas_pool_type), pointer :: cdgPool
      type (cdgQuadrature) :: quad2d

      call mpas_pool_get_subpool(domain % blocklist % structs, 'cdg', cdgPool)
      call mpas_pool_get_config(domain % configs, 'config_cdg_quad_order_horiz', quad % order_horiz)
      call mpas_pool_get_config(domain % configs, 'config_cdg_num_quad_points_horiz', quad % nh)
      call mpas_pool_get_config(domain % configs, 'config_cdg_quad_order_vert', quad % order_vert)
      call mpas_pool_get_config(domain % configs, 'config_cdg_quad_order_vert', quad % nv)

      ! Initialize the 2D quadrature points
      quad2d % fieldName = 'hq'
      call ocn_cdg_tri_init(domain, quad2d, err)
      allocate(quad % qw(quad % nh), quad % qx(quad % nh), quad % qy(quad % nh), quad % qz(quad % nh))
      quad % qw = quad2d % qw
      quad % qx = quad2d % qx
      quad % qy = quad2d % qy
      quad % qz = quad2d % qz
      deallocate(quad2d % qw, quad2d % qx, quad2d % qy, quad2d % qz)

   end subroutine ocn_cdg_tri_init_3d!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_init_quad_pts
!
!> \brief MPAS ocean CDG initialization of local and global quadratre 
!>        points
!> \author Dave Lee
!> \date   October 2015
!> \details
!>  This routine initializes the quadrature points in local and global
!   coordinates for each cell
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_tri_init_quad_pts(meshPool, cdgPool, onSphere, quad)!{{{

      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: cdgPool
      logical, intent(in) :: onSphere
      type (cdgQuadrature), intent(inout) :: quad

      integer, pointer :: nCells
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell
      real (kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex, xCell, yCell, zCell
      real (kind=RKIND), dimension(:,:,:), pointer :: transGtoC

      integer :: iVertex, jVertex, iEdge, jEdge, iCell
      real (kind=RKIND), dimension(2) :: p1, p2, p3

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'zCell', zCell)
      call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
      call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
      call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
      call mpas_pool_get_array(cdgPool, 'transGlobalToCell', transGtoC)

      call mpas_log_write( 'cdg: initializing quadrature points $i $i', intArgs=(/ quad % order, quad % nq /) )

      do iCell = 1, nCells
         do iEdge = 1, nEdgesOnCell(iCell)
            jEdge = mod(iEdge,nEdgesOnCell(iCell))+1
            iVertex = verticesOnCell(iEdge,iCell)
            jVertex = verticesOnCell(jEdge,iCell)

            if (onSphere) then
               ! Transform global coords to local cell coords and calculate the quadrature points
               call ocn_cdg_tri_trans_g2l_quad(quad, transGtoC(:,:,iCell), xCell(iCell), yCell(iCell), zCell(iCell), &
                                               xVertex(iVertex), yVertex(iVertex), zVertex(iVertex), xVertex(jVertex), yVertex(jVertex), zVertex(jVertex), &
                                               xCell(iCell), yCell(iCell), zCell(iCell), quad % xQuadL(:,iEdge,iCell), quad % yQuadL(:,iEdge,iCell))

               ! Transform the local cell coords of the quadrature points to global coords
               call ocn_cdg_tri_map_lp(xCell(iCell), yCell(iCell), zCell(iCell), quad % nq, quad % xQuadL(:,iEdge,iCell), quad % yQuadL(:,iEdge,iCell), &
                                       transGtoC(:,:,iCell), quad % xQuadG(:,iEdge,iCell), quad % yQuadG(:,iEdge,iCell), quad % zQuadG(:,iEdge,iCell))
            else
               p1(1) = xVertex(iVertex) - xCell(iCell)
               p1(2) = yVertex(iVertex) - yCell(iCell)
               p2(1) = xVertex(jVertex) - xCell(iCell)
               p2(2) = yVertex(jVertex) - yCell(iCell)
               p3(:) = 0.0_RKIND

               call ocn_cdg_tri_quad_pts(quad, p1, p2, p3, quad % xQuadL(:,iEdge,iCell), quad % yQuadL(:,iEdge,iCell))

               quad % xQuadG(:,iEdge,iCell) = quad % xQuadL(:,iEdge,iCell) + xCell(iCell)
               quad % yQuadG(:,iEdge,iCell) = quad % yQuadL(:,iEdge,iCell) + yCell(iCell)
               quad % zQuadG(:,iEdge,iCell) = 0.0_RKIND
            end if
         end do
      end do

   end subroutine ocn_cdg_tri_init_quad_pts!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_quad_pts
!
!> \brief MPAS ocean CDG triangle quadrature points in local coords
!> \author Dave Lee
!> \date   July 2015
!> \details
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

   subroutine ocn_cdg_tri_quad_pts(quad, p1, p2, p3, qx, qy)!{{{

      type (cdgQuadrature), intent(in) :: quad
      real(kind=RKIND), dimension(2), intent(in) :: p1, p2, p3
      real(kind=RKIND), dimension(:), intent(out) :: qx, qy

      integer :: iQuad

      do iQuad = 1, quad % nq
         qx(iQuad) = p1(1)*quad % qx(iQuad) + p2(1)*quad % qy(iQuad) + p3(1)*quad % qz(iQuad)
         qy(iQuad) = p1(2)*quad % qx(iQuad) + p2(2)*quad % qy(iQuad) + p3(2)*quad % qz(iQuad)
      end do

   end subroutine ocn_cdg_tri_quad_pts!}}}

   subroutine ocn_cdg_tri_quad_pts_3d(quad, p1, p2, p3, qx, qy)!{{{

      type (cdgQuad3D), intent(in) :: quad
      real(kind=RKIND), dimension(3), intent(in) :: p1, p2, p3
      real(kind=RKIND), dimension(:), intent(out) :: qx, qy

      integer :: iQuad

      do iQuad = 1, quad % nh
         qx(iQuad) = p1(1)*quad % qx(iQuad) + p2(1)*quad % qy(iQuad) + p3(1)*quad % qz(iQuad)
         qy(iQuad) = p1(2)*quad % qx(iQuad) + p2(2)*quad % qy(iQuad) + p3(2)*quad % qz(iQuad)
      end do

   end subroutine ocn_cdg_tri_quad_pts_3d!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_trans_g2l
!
!> \brief MPAS ocean CDG triangle global to local coordinate transform
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This routine transforms a global to local planar coordinates
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_tri_trans_g2l(transGtoL, xOrig, yOrig, zOrig, xGlobal, yGlobal, zGlobal, local)!{{{

      real (kind=RKIND), dimension(3,3), intent(in) :: transGtoL
      real (kind=RKIND), intent(in) :: xOrig, yOrig, zOrig, xGlobal, yGlobal, zGlobal
      real (kind=RKIND), dimension(2), intent(out) :: local

      real (kind=RKIND) :: dist
      real (kind=RKIND), dimension(3) :: gVec

      ! Global vectors from cell center to edge vertices
      gVec(1) = xGlobal - xOrig
      gVec(2) = yGlobal - yOrig
      gVec(3) = zGlobal - zOrig

      local(1) = transGtoL(1,1)*gVec(1) + transGtoL(1,2)*gVec(2) + transGtoL(1,3)*gVec(3)
      local(2) = transGtoL(2,1)*gVec(1) + transGtoL(2,2)*gVec(2) + transGtoL(2,3)*gVec(3)

      ! Ensure transformation is distance preserving
      !dist = ocn_arc_length(xGlobal, yGlobal, zGlobal, xOrig, yOrig, zOrig, radius)
      dist = mpas_arc_length(xGlobal, yGlobal, zGlobal, xOrig, yOrig, zOrig)
      dist = dist/sqrt(local(1)*local(1) + local(2)*local(2))
      local(:) = dist*local(:)

   end subroutine ocn_cdg_tri_trans_g2l!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_trans_g2l_quad
!
!> \brief MPAS ocean CDG triangle global to local coordinate transform
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  This routine transforms quadrature coords to local cell coordinates
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_tri_trans_g2l_quad(quad, transGtoL, xOrig, yOrig, zOrig, &!{{{
                                         xVert1, yVert1, zVert1, xVert2, yVert2, zVert2, xVert3, yVert3, zVert3, xQuad, yQuad)

      type (cdgQuadrature), intent(in) :: quad
      real (kind=RKIND), dimension(3,3), intent(in) :: transGtoL
      real (kind=RKIND), intent(in) :: xOrig, yOrig, zOrig, xVert1, yVert1, zVert1, xVert2, yVert2, zVert2, xVert3, yVert3, zVert3
      real (kind=RKIND), dimension(:), intent(out) :: xQuad, yQuad

      integer :: iQuad
      real (kind=RKIND) :: dist
      real (kind=RKIND), dimension(2) :: vert1Local, vert2Local, vert3Local
      real (kind=RKIND), dimension(3) :: vert1Vec, vert2Vec, vert3Vec

      ! Global vectors from cell center to edge vertices
      vert1Vec(1) = xVert1 - xOrig
      vert1Vec(2) = yVert1 - yOrig
      vert1Vec(3) = zVert1 - zOrig

      vert2Vec(1) = xVert2 - xOrig
      vert2Vec(2) = yVert2 - yOrig
      vert2Vec(3) = zVert2 - zOrig

      vert3Vec(1) = xVert3 - xOrig
      vert3Vec(2) = yVert3 - yOrig
      vert3Vec(3) = zVert3 - zOrig

      ! Project global vectors to local coordinates (preserving lengths on the plane
      vert1Local(1) = transGtoL(1,1)*vert1Vec(1) + transGtoL(1,2)*vert1Vec(2) + transGtoL(1,3)*vert1Vec(3)
      vert1Local(2) = transGtoL(2,1)*vert1Vec(1) + transGtoL(2,2)*vert1Vec(2) + transGtoL(2,3)*vert1Vec(3)
      !dist = ocn_arc_length(xOrig, yOrig, zOrig, xVert1, yVert1, zVert1, radius)
      dist = mpas_arc_length(xOrig, yOrig, zOrig, xVert1, yVert1, zVert1)
      if (dist .gt. 1.0E-12_RKIND) then
         dist = dist/sqrt(vert1Local(1)*vert1Local(1) + vert1Local(2)*vert1Local(2))
         vert1Local(:) = dist*vert1Local(:)
      end if

      vert2Local(1) = transGtoL(1,1)*vert2Vec(1) + transGtoL(1,2)*vert2Vec(2) + transGtoL(1,3)*vert2Vec(3)
      vert2Local(2) = transGtoL(2,1)*vert2Vec(1) + transGtoL(2,2)*vert2Vec(2) + transGtoL(2,3)*vert2Vec(3)
      !dist = ocn_arc_length(xOrig, yOrig, zOrig, xVert2, yVert2, zVert2, radius)
      dist = mpas_arc_length(xOrig, yOrig, zOrig, xVert2, yVert2, zVert2)
      if (dist .gt. 1.0E-12_RKIND) then
         dist = dist/sqrt(vert2Local(1)*vert2Local(1) + vert2Local(2)*vert2Local(2))
         vert2Local(:) = dist*vert2Local(:)
      end if

      vert3Local(1) = transGtoL(1,1)*vert3Vec(1) + transGtoL(1,2)*vert3Vec(2) + transGtoL(1,3)*vert3Vec(3)
      vert3Local(2) = transGtoL(2,1)*vert3Vec(1) + transGtoL(2,2)*vert3Vec(2) + transGtoL(2,3)*vert3Vec(3)
      !dist = ocn_arc_length(xOrig, yOrig, zOrig, xVert3, yVert3, zVert3, radius)
      dist = mpas_arc_length(xOrig, yOrig, zOrig, xVert3, yVert3, zVert3)
      if (dist .gt. 1.0E-12_RKIND) then
         dist = dist/sqrt(vert3Local(1)*vert3Local(1) + vert3Local(2)*vert3Local(2))
         vert3Local(:) = dist*vert3Local(:)
      end if

      do iQuad = 1, quad % nq
         xQuad(iQuad) = vert1Local(1)*quad % qx(iQuad) + vert2Local(1)*quad % qy(iQuad) + vert3Local(1)*quad % qz(iQuad)
         yQuad(iQuad) = vert1Local(2)*quad % qx(iQuad) + vert2Local(2)*quad % qy(iQuad) + vert3Local(2)*quad % qz(iQuad)
      end do

   end subroutine ocn_cdg_tri_trans_g2l_quad!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_integrate_full
!
!> \brief MPAS ocean cdg triangle surface integral
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  Integrate the volume of a distribution defined by an analytic 
!>  function and return the total volume
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_tri_integrate_full(domain, quad, volume)!{{{

      type (domain_type), intent(in) :: domain
      type (cdgQuadrature), intent(in) :: quad
      real (kind=RKIND), intent(out) :: volume

      type (dm_info), pointer :: dminfo

      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: cdgPool

      integer :: iCell, iEdge, jEdge, iQuad, iVertex, jVertex
      integer, pointer :: nCellsSolve
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: verticesOnCell

      real (kind=RKIND) :: resultLocal, radsq, detJ, xl, yl, xg, yg, zg
      real (kind=RKIND), pointer :: radius
      real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell, xVertex, yVertex, zVertex
      real (kind=RKIND), dimension(:,:), pointer :: triAreas

      real (kind=RKIND), dimension(2) :: pl0, pl1, pl2
      real (kind=RKIND), dimension(3) :: pg0, pg1, pg2
      real (kind=RKIND), dimension(3,3) :: A
      real (kind=RKIND), dimension(2,16) :: pq
      real (kind=RKIND), dimension(:), allocatable :: weight

      dminfo => domain % dminfo

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
         call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block_ptr % structs, 'cdg', cdgPool)

         call mpas_pool_get_config(meshPool, 'sphere_radius', radius)

         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)

         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
         call mpas_pool_get_array(meshPool, 'xCell', xCell)
         call mpas_pool_get_array(meshPool, 'yCell', yCell)
         call mpas_pool_get_array(meshPool, 'zCell', zCell)
         call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
         call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
         call mpas_pool_get_array(meshPool, 'zVertex', zVertex)

         call mpas_pool_get_array(cdgPool, 'triAreas', triAreas)

         allocate(weight(quad % nq))

         resultLocal = 0.0_RKIND

         do iCell = 1, nCellsSolve
            ! Test the cell rotation matrix
            call ocn_cdg_tri_cell_to_pole(xCell(iCell), yCell(iCell), zCell(iCell), radius, A)

            do iEdge = 1, nEdgesOnCell(iCell)
               jEdge = (mod(iEdge,nEdgesOnCell(iCell)))+1

               iVertex = verticesOnCell(iEdge,iCell)
               jVertex = verticesOnCell(jEdge,iCell)

               pg0(1) = 0.0_RKIND
               pg0(2) = 0.0_RKIND
               pg0(3) = 1.0_RKIND

               pg1(1) = A(1,1)*xVertex(iVertex) + A(1,2)*yVertex(iVertex) + A(1,3)*zVertex(iVertex)
               pg1(2) = A(2,1)*xVertex(iVertex) + A(2,2)*yVertex(iVertex) + A(2,3)*zVertex(iVertex)
               pg1(3) = A(3,1)*xVertex(iVertex) + A(3,2)*yVertex(iVertex) + A(3,3)*zVertex(iVertex)

               pg2(1) = A(1,1)*xVertex(jVertex) + A(1,2)*yVertex(jVertex) + A(1,3)*zVertex(jVertex)
               pg2(2) = A(2,1)*xVertex(jVertex) + A(2,2)*yVertex(jVertex) + A(2,3)*zVertex(jVertex)
               pg2(3) = A(3,1)*xVertex(jVertex) + A(3,2)*yVertex(jVertex) + A(3,3)*zVertex(jVertex)

               pl0(:) = 0.0_RKIND

               pl1(1) = radius*pg1(1)/pg1(3)
               pl1(2) = radius*pg1(2)/pg1(3)

               pl2(1) = radius*pg2(1)/pg2(3)
               pl2(2) = radius*pg2(2)/pg2(3)

               weight(:) = triAreas(iEdge,iCell)*quad % qw(:)

               do iQuad = 1, quad % nq
                  pq(1,iQuad) = pl1(1)*quad % qx(iQuad) + pl2(1)*quad % qy(iQuad) + pl0(1)*quad % qz(iQuad)
                  pq(2,iQuad) = pl1(2)*quad % qx(iQuad) + pl2(2)*quad % qy(iQuad) + pl0(2)*quad % qz(iQuad)

                  detJ = radius*radius*radius/((radius*radius + pq(1,iQuad)*pq(1,iQuad) + pq(2,iQuad)*pq(2,iQuad))**1.5_RKIND)

                  weight(iQuad) = weight(iQuad)*detJ

                  resultLocal = resultLocal + weight(iQuad)
               end do
            end do
         end do

         deallocate(weight)

         block_ptr => block_ptr % next
      end do

      ! Sum over all blocks
      call mpas_dmpar_sum_real(dminfo, resultLocal, volume)

   end subroutine ocn_cdg_tri_integrate_full!}}}

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  routine ocn_cdg_tri_test
!
!> \brief MPAS ocean cdg triangle unit test
!> \author Dave Lee
!> \date   June 2015
!> \details
!>  Integrate the volume of a distribution at a given resolution and 
!>  compare against an analytic solution to test convergence
!
!-----------------------------------------------------------------------

   subroutine ocn_cdg_tri_test(domain, quad, iErr)!{{{

      type (domain_type), intent(in) :: domain
      type (cdgQuadrature), intent(in) :: quad
      integer, intent(out) :: iErr

      real (kind=RKIND) :: pi = 3.1415926535897932384626433832795028841971_RKIND
      real (kind=RKIND) :: volA, volN

      iErr = 0

      if (domain % on_a_sphere .eqv. .false.) return

      call ocn_cdg_tri_integrate_full(domain, quad, volN)

      !volA = 2.0_RKIND*sqrt(2.0_RKIND)*pi/15.0_RKIND
      volA = 4.0_RKIND*pi

      ! Write result to the log file
      write(0,*) 'cdg: triangle test volume', volN, volA, 1.0_RKIND - volN/volA

   end subroutine ocn_cdg_tri_test!}}}

   ! Length preserving mapping from global to planar coordinates
   subroutine ocn_cdg_tri_map_lp(xo, yo, zo, np, xl, yl, g2l, xg, yg, zg)!{{{

      integer, intent(in) :: np
      real (kind=RKIND), intent(in) :: xo, yo, zo
      real (kind=RKIND), dimension(:), intent(in) :: xl, yl
      real (kind=RKIND), dimension(:,:), intent(in) :: g2l
      real (kind=RKIND), dimension(:), intent(out) :: xg, yg, zg

      integer :: i
      real (kind=RKIND), dimension(3) :: orig, da

      do i = 1, np
         orig(1) = xo
         orig(2) = yo
         orig(3) = zo

         ! Local to global mapping is the transpose of the global to local mapping
         da(1) = g2l(1,1)*xl(i) + g2l(2,1)*yl(i)
         da(2) = g2l(1,2)*xl(i) + g2l(2,2)*yl(i)
         da(3) = g2l(1,3)*xl(i) + g2l(2,3)*yl(i)

         !call particle_horizontal_movement(meshPool, orig, da)
         call mpas_particle_horizontal_movement(orig, da, .true.)

         xg(i) = orig(1)
         yg(i) = orig(2)
         zg(i) = orig(3)
      end do

   end subroutine ocn_cdg_tri_map_lp!}}}

   subroutine ocn_cdg_tri_cell_to_pole(x, y, z, R, A)!{{{

      real (kind=RKIND), intent(in) :: x, y, z, R
      real (kind=RKIND), dimension(3,3), intent(out) :: A

      integer :: i, j, k
      real (kind=RKIND) :: theta, phi
      real (kind=RKIND), dimension(3,3) :: Rz, Ry

      theta = atan2(y,x)
      phi = acos(z/R)

      ! Rotate -theta radians about the z-axis
      Rz = 0.0_RKIND
      Rz(3,3) = 1.0_RKIND
      Rz(1,1) = +cos(-theta)
      Rz(1,2) = -sin(-theta)
      Rz(2,1) = +sin(-theta)
      Rz(2,2) = +cos(-theta)

      ! Rotate -phi radians about the y-axis
      Ry = 0.0_RKIND
      Ry(2,2) = 1.0_RKIND
      Ry(1,1) = +cos(-phi)
      Ry(1,3) = +sin(-phi)
      Ry(3,1) = -sin(-phi)
      Ry(3,3) = +cos(-phi)

      A = 0.0_RKIND
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               A(i,j) = A(i,j) + Ry(i,k)*Rz(k,j)
            end do
         end do
      end do

   end subroutine ocn_cdg_tri_cell_to_pole!}}}

end module ocn_cdg_triangle

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
