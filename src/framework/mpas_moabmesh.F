

module mpas_moabmesh
!  use, intrinsic :: ISO_C_BINDING

  use mpas_log     
  use mpas_derived_types, only: dm_info, domain_type
  use mpas_field_routines
  use mpas_sort
  use mpas_stream_manager
  use mpas_pool_routines
  !use mpas_vector_operations

  implicit none
#include "moab/MOABConfig.h"

  integer, public :: MPOID !   app id on moab side, for MPAS ocean

  contains

  SUBROUTINE errorout(ierr, message)
  integer ierr
  character*(*) message
  if (ierr.ne.0) then
    print *, message
    call exit (1)
  end if
  return
  end subroutine
  
  subroutine mpas_moab_instance(domain)
    
    type (domain_type), intent(inout) :: domain

    type (block_type), pointer :: block
    type (mpas_pool_type), pointer :: meshPool
    integer, pointer :: nCells, nVertices, maxEdges
    integer  :: pid, nblocks
    real (kind=RKIND), dimension(:), allocatable :: xv,yv,zv
    integer, dimension(:,:), pointer :: verticesOnCell
    integer, dimension(:), pointer :: nEdgesOnCell
    integer, dimension(:), pointer :: indexToVertexID, indexToCellID
    real(kind=RKIND), dimension(:), pointer :: xVertex, yVertex, zVertex
    real (kind=RKIND), dimension(:), pointer :: xCell, yCell, zCell
    logical, pointer :: on_a_sphere, is_periodic
    real(kind=RKIND), pointer :: x_period, y_period
    integer, pointer :: nCellsSolve, nEdgesSolve, nVerticesSolve
    
   

    integer , external :: iMOAB_RegisterFortranApplication, &
        iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_DetermineGhostEntities, &
        iMOAB_DefineTagStorage,  iMOAB_SetIntTagStorage
    integer :: c_comm, i1, j1, ic, lastvertex
    character*12 appname
    integer :: ierr, num_verts_in_cells, ext_comp_id
    real(kind=RKIND), allocatable, target :: moab_vert_coords(:)
    integer dimcoord, dimen, mbtype, block_ID, proc_id
    integer  ,allocatable , target :: all_connects(:)
    character*100 outfile, wopts, localmeshfile, lnum, tagname
    integer tagtype, numco, tag_sto_len, ent_type, tagindex

    c_comm = domain % dminfo % comm
    appname = 'MPAS_MESH'//CHAR(0)
    MPOID = -1 ! initialized negative , so we know if this register passes
    ext_comp_id = 17  ! some number less than 33 ! maybe this should be input
    ierr = iMOAB_RegisterFortranApplication(appname, c_comm, ext_comp_id, pid)
    MPOID = pid ! this is exported, need for send to work
    call errorout(ierr, 'fail to register MPAS_MOAB mesh')
    proc_id = domain % dminfo % my_proc_id
    call mpas_log_write('MOAB MPAS app pid:  $i  task  $i ', intArgs=(/pid, proc_id/) )

!   blocks should be merged if there is more than one block per task
    nblocks = 0
    block => domain % blocklist
    do while (associated(block)) !{{{
       nblocks = nblocks + 1
        ! allocate scratch memory 
       call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
       call mpas_pool_get_array(meshPool, 'xVertex', xVertex)
       call mpas_pool_get_array(meshPool, 'yVertex', yVertex)
       call mpas_pool_get_array(meshPool, 'zVertex', zVertex)
       call mpas_pool_get_dimension(meshPool, 'nVertices', nVertices)
       call mpas_pool_get_dimension(meshPool, 'maxEdges', maxEdges)
       call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
       call mpas_pool_get_config(meshPool, 'is_periodic', is_periodic)
       call mpas_pool_get_config(meshPool, 'x_period', x_period)
       call mpas_pool_get_config(meshPool, 'y_period', y_period)
       call mpas_pool_get_array(meshPool, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(meshPool, 'indexToVertexID', indexToVertexID)
       call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
       call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
       call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve', nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve', nVerticesSolve)
!       call mpas_pool_get_array(meshPool, 'xCell', xCell)
!       call mpas_pool_get_array(meshPool, 'yCell', yCell)
!       call mpas_pool_get_array(meshPool, 'zCell', zCell)

#ifdef MPAS_DEBUG
       call mpas_log_write(' number of vertices:: $i number of cells:: $i solve: v:$i c:$i', intArgs=(/nVertices, nCells, nVerticesSolve, nCellsSolve/) )
#endif
        ! first create local vertices, nVerticesSolve
        allocate(moab_vert_coords(3*nVertices))
        do i1 =1, nVertices
          moab_vert_coords(3*i1-2) = xVertex(i1)
          moab_vert_coords(3*i1-1) = yVertex(i1)
          moab_vert_coords(3*i1  ) = zVertex(i1)
          ! call mpas_log_write('i:: $i coords:: $r $r $r $r', intArgs=(/i1/), realArgs=(/moab_vert_coords(3*i1-2),moab_vert_coords(3*i1-1), moab_vert_coords(3*i1)/) )
        enddo
        dimcoord = 3*nVertices
        dimen = 3
        ierr = iMOAB_CreateVertices(pid, dimcoord, dimen, moab_vert_coords)
        call errorout(ierr, 'fail to create vertices')
    ! fill now connectivity array, nCellsSolve; fist pad to max nc
        num_verts_in_cells = nCellsSolve * maxEdges
        allocate(all_connects(num_verts_in_cells))
!    collect all vertices, and also pad 
        j1 = 0
        do ic=1, nCellsSolve
          do i1 = 1, nEdgesOnCell(ic)
            j1 = j1 + 1
            all_connects(j1) = verticesOnCell( i1, ic)
          enddo
          lastvertex = verticesOnCell( nEdgesOnCell (ic), ic)
          ! pad the rest with the last vertex
          do i1 =   nEdgesOnCell (ic) + 1, maxEdges
            j1 = j1 + 1
            all_connects(j1) = lastvertex ! repeat the last vertex (pad)
          enddo
          ! call mpas_log_write('cell: $i  v:: $i $i $i $i $i $i', intArgs=(/ic, all_connects(j1-5), all_connects(j1-4), all_connects(j1-3), all_connects(j1-2), all_connects(j1-1), all_connects(j1)/) )
       enddo

       mbtype = 4 !  polygon

       block_ID = 100 + nblocks ! this will be for coarse mesh

       ierr = iMOAB_CreateElements( pid, nCellsSolve, mbtype, maxEdges, all_connects, block_ID );
       call errorout(ierr, 'fail to create polygons')
!      set the global id for vertices
!   first, retrieve the tag
       tagname='GLOBAL_ID'//CHAR(0)
       tagtype = 0  ! dense, integer
       numco = 1
       ierr = iMOAB_DefineTagStorage(pid, tagname, tagtype, numco,  tagindex )
       call errorout(ierr, 'fail to get global id tag')
! now set the values
       ent_type = 0 ! vertex type
       ierr = iMOAB_SetIntTagStorage ( pid, tagname, nVertices , ent_type, indexToVertexID )
       call errorout(ierr, 'fail to set global id tag for vertices')
       ! set global id tag for elements
       ent_type = 1 ! now set the global id tag on elements
       ierr = iMOAB_SetIntTagStorage ( pid, tagname, nCellsSolve, ent_type, indexToCellID)
       call errorout(ierr, 'fail to set global id tag for polygons')
       ! get next block
#ifdef MPAS_DEBUG
       if (proc_id.lt. 5) then
         write(lnum,"(I0.2)")proc_id
         localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
         wopts = CHAR(0)
         ierr = iMOAB_WriteMesh(pid, localmeshfile, wopts)
         call errorout(ierr, 'fail to write local mesh file')
       endif
#endif
       ierr = iMOAB_ResolveSharedEntities( pid, nVertices, indexToVertexID );
       call errorout(ierr, 'fail to resolve shared entities')

       deallocate (moab_vert_coords)
       deallocate (all_connects)
       block => block % next

    end do !}}}
    

    if (nblocks .ne. 1) then
      call errorout(1, 'more than one block per task')
    endif 
    !     write out the mesh file to disk, in parallel
    outfile = 'whole.h5m'//CHAR(0)
    wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
    ierr = iMOAB_WriteMesh(pid, outfile, wopts)
    call errorout(ierr, 'fail to write the mesh file')
    
    
  end subroutine mpas_moab_instance
end module mpas_moabmesh


