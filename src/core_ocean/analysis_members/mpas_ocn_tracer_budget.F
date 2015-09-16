! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_budget
!
!> \brief MPAS ocean analysis mode member: tracer_budget
!> \author Milena Veneziani
!> \date   Sept 2015
!> \details
!>  MPAS ocean analysis mode member: tracer_budget
!>  
!>  5. In src/core_ocean/analysis_members/Makefile, add your
!>     new analysis member to the list of members. See another analysis member
!>     in that file for an example. 
!>     NOTE: If your analysis member depends on other files, add a dependency
!>           line for the member and list them there. See okubo weiss for an example.
!>
!-----------------------------------------------------------------------

module ocn_tracer_budget

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_timekeeping
   use mpas_stream_manager

   use ocn_constants
   use ocn_diagnostics_routines

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

   public :: ocn_init_tracer_budget, &
             ocn_compute_tracer_budget, &
             ocn_restart_tracer_budget, &
             ocn_finalize_tracer_budget

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_tracer_budget
!
!> \brief   Initialize MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all initializations required for the
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_tracer_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_init_tracer_budget!}}}

!***********************************************************************
!
!  routine ocn_compute_tracer_budget
!
!> \brief   Compute MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all computation required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_compute_tracer_budget(domain, timeLevel, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      integer, intent(in) :: timeLevel

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracerBudgetAMPool
      type (dm_info) :: dminfo
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: statePool
      type (mpas_pool_type), pointer :: meshPool
      type (mpas_pool_type), pointer :: scratchPool
      type (mpas_pool_type), pointer :: diagnosticsPool
      type (mpas_pool_type), pointer :: tracerBudgetAM

      ! Here are some example variables which may be needed for your analysis member
      integer, pointer :: nVertLevels, nCellsSolve, nEdgesSolve, nVerticesSolve, num_tracers
      integer :: iTracer, k, iCell
      integer, dimension(:), pointer :: maxLevelCell, maxLevelEdgeTop, maxLevelVertexBot

      real (kind=RKIND), dimension(:), pointer ::  areaCell, dcEdge, dvEdge

      err = 0

      dminfo = domain % dminfo

      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'state', statePool)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
         call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
         call mpas_pool_get_subpool(block % structs, 'tracerBudgetAM', tracerBudgetAMPool)

         ! Here are some example variables which may be needed for your analysis member
         call mpas_pool_get_dimension(statePool, 'num_tracers', num_tracers)

         call mpas_pool_get_dimension(block % dimensions, 'nVertLevels', nVertLevels)
         call mpas_pool_get_dimension(block % dimensions, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nEdgesSolve', nEdgesSolve)
         call mpas_pool_get_dimension(block % dimensions, 'nVerticesSolve', nVerticesSolve)

         call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
         call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
         call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
         call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
         call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
         call mpas_pool_get_array(meshPool, 'maxLevelVertexBot', maxLevelVertexBot)

         ! Computations which are functions of nCells, nEdges, or nVertices
         ! must be placed within this block loop
         ! Here are some example loops
         do iCell = 1,nCellsSolve
            do k = 1, maxLevelCell(iCell)
               do iTracer = 1, num_tracers
               ! computations on tracers(iTracer,k, iCell)
               end do
            end do
         end do

         block => block % next
      end do

      ! mpi gather/scatter calls may be placed here.
      ! Here are some examples.  See mpas_oac_global_stats.F for further details.
!      call mpas_dmpar_sum_real_array(dminfo, nVariables, sumSquares(1:nVariables), reductions(1:nVariables))
!      call mpas_dmpar_min_real_array(dminfo, nMins, mins(1:nMins), reductions(1:nMins))
!      call mpas_dmpar_max_real_array(dminfo, nMaxes, maxes(1:nMaxes), reductions(1:nMaxes))

      ! Even though some variables do not include an index that is decomposed amongst
      ! domain partitions, we assign them within a block loop so that all blocks have the
      ! correct values for writing output.
      block => domain % blocklist
      do while (associated(block))
         call mpas_pool_get_subpool(block % structs, 'tracerBudgetAM', tracerBudgetAMPool)

         ! assignment of final tracerBudgetAM variables could occur here.

         block => block % next
      end do

      write(stdoutUnit,*) ' **** THIS IS A TEST ***** '

   end subroutine ocn_compute_tracer_budget!}}}

!***********************************************************************
!
!  routine ocn_restart_tracer_budget
!
!> \brief   Save restart for MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts computation required to save a restart state
!>  for the MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_restart_tracer_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_restart_tracer_budget!}}}

!***********************************************************************
!
!  routine ocn_finalize_tracer_budget
!
!> \brief   Finalize MPAS-Ocean analysis member
!> \author  Milena Veneziani
!> \date    Sept 2015
!> \details
!>  This routine conducts all finalizations required for this
!>  MPAS-Ocean analysis member.
!
!-----------------------------------------------------------------------

   subroutine ocn_finalize_tracer_budget(domain, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (domain_type), intent(inout) :: domain

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      err = 0

   end subroutine ocn_finalize_tracer_budget!}}}

end module ocn_tracer_budget

! vim: foldmethod=marker
