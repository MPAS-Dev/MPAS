

module mpas_moabmesh
!  use, intrinsic :: ISO_C_BINDING

  use mpas_log     
  use mpas_derived_types, only: dm_info, domain_type
  use mpas_field_routines
  use mpas_sort

  implicit none
#include "moab/MOABConfig.h"

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

  subroutine create_mpas_moabmesh(domain, nEdgesOnCell, verticesOnCell, &
     xVertex, yVertex, zVertex )

    type (field2DInteger), pointer :: verticesOnCell
    type (field1DInteger), pointer :: nEdgesOnCell
    type (field1dReal), pointer :: xVertex, yVertex, zVertex
    type (domain_type), pointer :: domain
    integer :: c_comm, pid, i1, j1, ic, lastvertex
    character*12 appname

    integer :: num_verts_cells  !   total number of vertices in the cells
    integer  ,allocatable , target :: all_connects(:)
    integer  ,allocatable , target :: moabvh(:)
    integer nv, currentval
    ! RKIND must be double
    real(kind=RKIND), allocatable, target :: moab_vert_coords(:)
    integer dimcoord, dimen, nLocalCells, maxv, nedg, mbtype, block_ID, proc_id

    character*100 outfile, wopts, localmeshfile, lnum
    integer num_layers, dimgh, bridge !  for ghosting


! do we really need this?
    integer , external :: iMOAB_InitializeFortran, iMOAB_RegisterFortranApplication, &
        iMOAB_CreateVertices, iMOAB_WriteMesh, iMOAB_CreateElements, &
        iMOAB_ResolveSharedEntities, iMOAB_DetermineGhostEntities

    integer :: ierr, num_verts_in_cells, ext_comp_id

    c_comm = domain % dminfo % comm
    appname = 'MPAS_MESH'//CHAR(0)
    ext_comp_id = 17  ! some number less than 33
    ierr = iMOAB_RegisterFortranApplication(appname, c_comm, ext_comp_id, pid)
    call errorout(ierr, 'fail to register MPAS_MOAB mesh')
    proc_id = domain % dminfo % my_proc_id
    call mpas_log_write('MOAB MPAS app pid:  $i  task  $i ', intArgs=(/pid, proc_id/) )

    nLocalCells=verticesOnCell%dimSizes(2)
    maxv = verticesOnCell%dimSizes(1)
    call mpas_log_write(' vertsOnCell:  $i  $i ', intArgs=(/maxv, nLocalCells/) )
    call mpas_log_write(' nEdgesOnCell: $i ', intArgs=(/nEdgesOnCell%dimSizes(1)/) )
    call mpas_log_write(' xVertex: $i ', intArgs=(/xVertex%dimSizes(1)/) )

    num_verts_in_cells = nLocalCells * maxv
    allocate(all_connects(num_verts_in_cells))
!    collect all vertices
    j1 = 0
    do ic=1, nLocalCells
       do i1 = 1, nEdgesOnCell % array(ic)
         j1 = j1 + 1
         all_connects(j1) = verticesOnCell % array( i1, ic)
       enddo
       lastvertex = verticesOnCell % array( nEdgesOnCell % array(ic), ic)
       ! pad the rest with the last vertex
       do i1 =  nEdgesOnCell % array(ic) + 1, maxv
         j1 = j1 + 1
         all_connects(j1) = lastvertex ! repeat the last vertex (pad)
       enddo
    enddo

    call mpas_quicksort(num_verts_in_cells, all_connects)
    ! now compress the duplicates, keep only one vertex; then do binary search

    ! allocate(moabvh(num_verts_in_cells)) ! used to compress

   !  allocate(moabconn(moab_dim_cquads))
    nv=1
    currentval = all_connects(1)
    do i1=2,num_verts_in_cells
      if (all_connects(i1) .ne. currentval ) then
        nv=nv+1
        currentval = all_connects(i1)
        all_connects(nv) = currentval
      endif
    enddo
    call mpas_log_write(' number of vertices:: $i ', intArgs=(/nv/) )
!     write (stderrUnit,100) (all_connects(j1), j1=1,nv)
! 100 format(10I6)
    ! num local vertices = nv; moabnv will be the global ids, that will be used for resolve sharing
    allocate(moabvh(nv))
    !do i1=1,nv
    moabvh(:) = all_connects(1:nv) ! copy the first
    !enddo
    ! first create local vertices, accumulate xv, yv, zv from xVertex, ...
    allocate(moab_vert_coords(3*nv))
    do i1 =1, nv
      j1 = moabvh(i1)
      moab_vert_coords(3*i1-2) = xVertex % array(j1)
      moab_vert_coords(3*i1-1) = yVertex % array(j1)
      moab_vert_coords(3*i1  ) = zVertex % array(j1)
    enddo
    dimcoord = 3*nv
    dimen = 3
    ierr = iMOAB_CreateVertices(pid, dimcoord, dimen, moab_vert_coords)
    call errorout(ierr, 'fail to create vertices')
    ! build now connectivity arrays using binary search
    j1 = 0
    do ic=1, nLocalCells
       nedg = nEdgesOnCell % array(ic)
       do i1 = 1, nedg
         j1=j1+1
         all_connects(j1) = mpas_binary_search(moabvh, 1, 1, nv, verticesOnCell % array( i1, ic) )

       enddo
       lastvertex = all_connects(j1)
       ! pad the rest with the last vertex
       do i1 = nedg + 1, maxv
         j1 = j1 + 1
         all_connects(j1) = lastvertex ! repeat the last vertex (pad)
       enddo
    enddo
    mbtype = 4 !  polygon

    block_ID = 100 ! this will be for coarse mesh

    ierr = iMOAB_CreateElements( pid, nLocalCells, mbtype, maxv, all_connects, block_ID );
    call errorout(ierr, 'fail to create polygons')

    ierr = iMOAB_ResolveSharedEntities( pid, nv, moabvh );
    call errorout(ierr, 'fail to resolve shared entities')

! write in serial, on each task, before ghosting
    if (proc_id.lt. 10) then
        write(lnum,"(I0.2)")proc_id
        localmeshfile = 'owned_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(pid, localmeshfile, wopts)
        call errorout(ierr, 'fail to write local mesh file')
      endif

!     write out the mesh file to disk, in parallel
      outfile = 'whole.h5m'//CHAR(0)
      wopts   = 'PARALLEL=WRITE_PART'//CHAR(0)
      ierr = iMOAB_WriteMesh(pid, outfile, wopts)
      call errorout(ierr, 'fail to write the mesh file')

      ! (iMOAB_AppID pid, int * ghost_dim, int *num_ghost_layers, int * bridge_dim )
      dimgh = 2 ! will ghost polygons, topological dim 2
      bridge = 0 ! use vertex as bridge
      num_layers = 1 ! so far, one layer only
      ierr = iMOAB_DetermineGhostEntities( pid, dimgh, num_layers, bridge)
      call errorout(ierr, 'fail to determine ghosts')

      ! write in serial, on each task
      if (proc_id .lt. 10) then
        write(lnum,"(I0.2)")proc_id
        localmeshfile = 'localmesh_'//trim(lnum)// '.h5m' // CHAR(0)
        wopts = CHAR(0)
        ierr = iMOAB_WriteMesh(pid, localmeshfile, wopts)
        call errorout(ierr, 'fail to write local mesh file')
      endif

    deallocate (moab_vert_coords)
    deallocate (moabvh)
    deallocate (all_connects)
  end subroutine create_mpas_moabmesh
  

end module mpas_moabmesh


