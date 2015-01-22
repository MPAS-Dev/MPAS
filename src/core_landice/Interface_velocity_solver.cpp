// ===================================================
//! Includes
// ===================================================

#include <algorithm>
#include "Interface_velocity_solver.hpp"
//#include <lifev/ice_sheet/interface_with_mpas/Interface.hpp>
//#include <lifev/ice_sheet/solver/BuildMeshFromBareData.hpp>
//#include <lifev/ice_sheet/solver/Extrude3DMesh.hpp>

// ===================================================
//! Namespaces
// ===================================================

//typedef std::list<exchange> exchangeList_Type;

// ice_problem pointer

int Ordering = 0; //ordering ==0 means that the mesh is extruded layerwise, whereas ordering==1 means that the mesh is extruded columnwise.
MPI_Comm comm, reducedComm;
bool isDomainEmpty = true;
bool initialize_velocity = true;
bool first_time_step = true;
int nCells_F, nEdges_F, nVertices_F;
int nCellsSolve_F, nEdgesSolve_F, nVerticesSolve_F;
int nVertices, nEdges, nTriangles, nGlobalVertices, nGlobalEdges,
    nGlobalTriangles;
int maxNEdgesOnCell_F;
int const *cellsOnEdge_F, *cellsOnVertex_F, *verticesOnCell_F,
    *verticesOnEdge_F, *edgesOnCell_F, *indexToCellID_F, *nEdgesOnCells_F;
std::vector<double> layersRatio, levelsNormalizedThickness;
int nLayers;
double const *xCell_F, *yCell_F, *zCell_F, *areaTriangle_F;
std::vector<double> xCellProjected, yCellProjected, zCellProjected;
const double unit_length = 1000;
const double T0 = 273.15;
const double minThick = 1e-2; //10m
const double minBeta = 1e-5;
//void *phgGrid = 0;
std::vector<int> edgesToReceive, fCellsToReceive, indexToTriangleID,
    verticesOnTria, trianglesOnEdge, trianglesPositionsOnEdge, verticesOnEdge;
std::vector<int> indexToVertexID, vertexToFCell, indexToEdgeID, edgeToFEdge,
    mask, fVertexToTriangleID, fCellToVertex;
std::vector<double> temperatureOnTetra, velocityOnVertices, velocityOnCells,
    elevationData, thicknessData, betaData, smb_F, thicknessOnCells;
std::vector<bool> isVertexBoundary, isBoundaryEdge;
;
int numBoundaryEdges;
double radius;

exchangeList_Type const *sendCellsList_F = 0, *recvCellsList_F = 0;
exchangeList_Type const *sendEdgesList_F = 0, *recvEdgesList_F = 0;
exchangeList_Type const *sendVerticesList_F = 0, *recvVerticesList_F = 0;
exchangeList_Type sendCellsListReversed, recvCellsListReversed,
    sendEdgesListReversed, recvEdgesListReversed;

exchange::exchange(int _procID, int const* vec_first, int const* vec_last,
    int fieldDim) :
    procID(_procID), vec(vec_first, vec_last), buffer(
        fieldDim * (vec_last - vec_first)), doubleBuffer(
        fieldDim * (vec_last - vec_first)) {
}

extern "C" {

// ===================================================
//! Interface functions
// ===================================================
int velocity_solver_init_mpi(int* fComm) {
  // get MPI_Comm from Fortran
  comm = MPI_Comm_f2c(*fComm);

  return 0;
}

void velocity_solver_export_2d_data(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F) {
  if (isDomainEmpty)
    return;

  import2DFields(lowerSurface_F, thickness_F, beta_F, minThick);

  velocity_solver_export_2d_data__(reducedComm, elevationData, thicknessData,
      betaData, indexToVertexID);
}

void velocity_solver_set_grid_data(int const* _nCells_F, int const* _nEdges_F,
    int const* _nVertices_F, int const* _nLayers, int const* _nCellsSolve_F,
    int const* _nEdgesSolve_F, int const* _nVerticesSolve_F,
    int const* _maxNEdgesOnCell_F, double const* radius_F,
    int const* _cellsOnEdge_F, int const* _cellsOnVertex_F,
    int const* _verticesOnCell_F, int const* _verticesOnEdge_F,
    int const* _edgesOnCell_F, int const* _nEdgesOnCells_F,
    int const* _indexToCellID_F, double const* _xCell_F, double const* _yCell_F,
    double const* _zCell_F, double const* _areaTriangle_F,
    int const* sendCellsArray_F, int const* recvCellsArray_F,
    int const* sendEdgesArray_F, int const* recvEdgesArray_F,
    int const* sendVerticesArray_F, int const* recvVerticesArray_F) {

  nCells_F = *_nCells_F;
  nEdges_F = *_nEdges_F;
  nVertices_F = *_nVertices_F;
  nLayers = *_nLayers;
  nCellsSolve_F = *_nCellsSolve_F;
  nEdgesSolve_F = *_nEdgesSolve_F;
  nVerticesSolve_F = *_nVerticesSolve_F;
  maxNEdgesOnCell_F = *_maxNEdgesOnCell_F;
  radius = *radius_F;
  cellsOnEdge_F = _cellsOnEdge_F;
  cellsOnVertex_F = _cellsOnVertex_F;
  verticesOnCell_F = _verticesOnCell_F;
  verticesOnEdge_F = _verticesOnEdge_F;
  edgesOnCell_F = _edgesOnCell_F;
  nEdgesOnCells_F = _nEdgesOnCells_F;
  indexToCellID_F = _indexToCellID_F;
  xCell_F = _xCell_F;
  yCell_F = _yCell_F;
  zCell_F = _zCell_F;
  areaTriangle_F = _areaTriangle_F;
  mask.resize(nVertices_F);

  thicknessOnCells.resize(nCellsSolve_F);

  sendCellsList_F = new exchangeList_Type(unpackMpiArray(sendCellsArray_F));
  recvCellsList_F = new exchangeList_Type(unpackMpiArray(recvCellsArray_F));
  sendEdgesList_F = new exchangeList_Type(unpackMpiArray(sendEdgesArray_F));
  recvEdgesList_F = new exchangeList_Type(unpackMpiArray(recvEdgesArray_F));
  sendVerticesList_F = new exchangeList_Type(
      unpackMpiArray(sendVerticesArray_F));
  recvVerticesList_F = new exchangeList_Type(
      unpackMpiArray(recvVerticesArray_F));

  if (radius > 10) {
    xCellProjected.resize(nCells_F);
    yCellProjected.resize(nCells_F);
    zCellProjected.assign(nCells_F, 0.);
    for (int i = 0; i < nCells_F; i++) {
      double r = std::sqrt(
          xCell_F[i] * xCell_F[i] + yCell_F[i] * yCell_F[i]
              + zCell_F[i] * zCell_F[i]);
      xCellProjected[i] = radius * std::asin(xCell_F[i] / r);
      yCellProjected[i] = radius * std::asin(yCell_F[i] / r);
    }
    xCell_F = &xCellProjected[0];
    yCell_F = &yCellProjected[0];
    zCell_F = &zCellProjected[0];
  }
}

void velocity_solver_init_l1l2(double const* levelsRatio_F) {
#ifdef LIFEV
  velocityOnVertices.resize(2 * nVertices * (nLayers + 1), 0.);
  velocityOnCells.resize(2 * nCells_F * (nLayers + 1), 0.);

  if (isDomainEmpty)
    return;

  layersRatio.resize(nLayers);
  // !!Indexing of layers is reversed
  for (int i = 0; i < nLayers; i++)
    layersRatio[i] = levelsRatio_F[nLayers - 1 - i];
  //std::copy(levelsRatio_F, levelsRatio_F+nLayers, layersRatio.begin());
  mapCellsToVertices(velocityOnCells, velocityOnVertices, 2, nLayers, Ordering);

  velocity_solver_init_l1l2__(layersRatio, velocityOnVertices, initialize_velocity);
  initialize_velocity = false;
#endif
}




void velocity_solver_solve_l1l2(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F,
    double const* temperature_F, double* u_normal_F, double* heatIntegral_F,
    double* viscosity_F) {

#ifdef LIFEV

  std::fill(u_normal_F, u_normal_F + nEdges_F * nLayers, 0.);

  double localSum(0), sum(0);

  for (int i = 0; i < nCellsSolve_F; i++) {
    localSum = std::max(localSum,
        std::fabs(thickness_F[i] - thicknessOnCells[i]));
  }

  MPI_Allreduce(&localSum, &sum, 1, MPI_DOUBLE, MPI_MAX, comm);

  std::cout << "Thickness change: " << sum << std::endl;
  std::copy(thickness_F, thickness_F + nCellsSolve_F, &thicknessOnCells[0]);

  if (!isDomainEmpty) {
    std::vector<double> temperatureData(nLayers * nVertices);

    import2DFields(lowerSurface_F, thickness_F, beta_F, minThick);

    for (int index = 0; index < nVertices; index++) {
      int iCell = vertexToFCell[index];
      for (int il = 0; il < nLayers; il++) {
        temperatureData[index + il * nVertices] = temperature_F[iCell * nLayers
            + (nLayers - il - 1)] + T0;
      }
    }


    velocity_solver_solve_l1l2__(elevationData, thicknessData, betaData,
        temperatureData, indexToVertexID, velocityOnVertices);



    get_prism_velocity_on_FEdges(u_normal_F, velocityOnVertices, edgeToFEdge);
  }


   mapVerticesToCells (velocityOnVertices, &velocityOnCells[0], 2, nLayers, Ordering);

   allToAll (u_normal_F,  &sendEdgesListReversed, &recvEdgesListReversed, nLayers);
   allToAll (u_normal_F,  sendEdgesList_F, recvEdgesList_F, nLayers);


#endif

}

void velocity_solver_export_l1l2_velocity() {

  if (isDomainEmpty)
     return;

  std::vector<double> regulThk(thicknessData);
      for (int index = 0; index < nVertices; index++)
        regulThk[index] = std::max(1e-4, thicknessData[index]);

      std::vector<int> mpasIndexToVertexID(nVertices);
      for (int i = 0; i < nVertices; i++) {
        mpasIndexToVertexID[i] = indexToCellID_F[vertexToFCell[i]];
      }
#ifdef LIFEV
  velocity_solver_export_l1l2_velocity__(layersRatio, elevationData, regulThk, mpasIndexToVertexID, reducedComm);
#endif
}

void velocity_solver_init_fo(double const *levelsRatio_F) {

  velocityOnVertices.resize(2 * nVertices * (nLayers + 1), 0.);
  velocityOnCells.resize(2 * nCells_F * (nLayers + 1), 0.);

  if (isDomainEmpty)
    return;

  layersRatio.resize(nLayers);
  // !!Indexing of layers is reversed
  for (int i = 0; i < nLayers; i++)
    layersRatio[i] = levelsRatio_F[nLayers - 1 - i];
  //std::copy(levelsRatio_F, levelsRatio_F+nLayers, layersRatio.begin());

  mapCellsToVertices(velocityOnCells, velocityOnVertices, 2, nLayers, Ordering);

#ifdef LIFEV
  velocity_solver_init_fo__(layersRatio, velocityOnVertices, indexToVertexID, initialize_velocity);
#endif
  //    iceProblemPtr->initializeSolverFO(layersRatio, velocityOnVertices, thicknessData, elevationData, indexToVertexID, initialize_velocity);
  initialize_velocity = false;
}

void velocity_solver_solve_fo(double const* lowerSurface_F,
    double const* thickness_F, double const* beta_F,
    double const* temperature_F, double* u_normal_F, double* heatIntegral_F,
    double* viscosity_F) {

  std::fill(u_normal_F, u_normal_F + nEdges_F * nLayers, 0.);

  if (!isDomainEmpty) {

#ifdef LIFEV
   double localSum(0), sum(0);

   for (int i = 0; i < nCellsSolve_F; i++) {
   localSum = std::max(localSum,
   std::fabs(thickness_F[i] - thicknessOnCells[i]));
   }

   MPI_Allreduce(&localSum, &sum, 1, MPI_DOUBLE, MPI_MAX, comm);

   std::cout << "Thickness change: " << sum << std::endl;
   std::copy(thickness_F, thickness_F + nCellsSolve_F, &thicknessOnCells[0]);
#endif



    import2DFields(lowerSurface_F, thickness_F, beta_F, minThick);

    std::vector<double> regulThk(thicknessData);
    for (int index = 0; index < nVertices; index++)
      regulThk[index] = std::max(1e-4, thicknessData[index]);

    importP0Temperature(temperature_F);

    velocity_solver_solve_fo__(nLayers, nGlobalVertices, nGlobalTriangles,
        Ordering, first_time_step, indexToVertexID, indexToTriangleID, minBeta,
        regulThk, levelsNormalizedThickness, elevationData, thicknessData,
        betaData, temperatureOnTetra, velocityOnVertices);

    std::vector<int> mpasIndexToVertexID(nVertices);
    for (int i = 0; i < nVertices; i++) {
      mpasIndexToVertexID[i] = indexToCellID_F[vertexToFCell[i]];
    }
    get_tetraP1_velocity_on_FEdges(u_normal_F, velocityOnVertices, edgeToFEdge,
        mpasIndexToVertexID);
  }

  mapVerticesToCells(velocityOnVertices, &velocityOnCells[0], 2, nLayers,
      Ordering);

  std::vector<double> velOnEdges(nEdges * nLayers);
  for (int i = 0; i < nEdges; i++) {
    for (int il = 0; il < nLayers; il++) {
      velOnEdges[i * nLayers + il] = u_normal_F[edgeToFEdge[i] * nLayers + il];
    }
  }

  allToAll(u_normal_F, &sendEdgesListReversed, &recvEdgesListReversed, nLayers);

  allToAll(u_normal_F, sendEdgesList_F, recvEdgesList_F, nLayers);

  first_time_step = false;


#ifdef LIFEV

  std::vector<int> edgesProcId(nEdges_F), trianglesProcIds(nVertices_F);
  getProcIds(edgesProcId, recvEdgesList_F);
  getProcIds(trianglesProcIds, recvVerticesList_F);

  int localSumInt(0), sumInt(0);

  for (int i = 0; i < nEdges; i++) {
    for (int il = 0; il < 1; il++) {
      if (std::fabs(
          velOnEdges[i * nLayers + il]
              - u_normal_F[edgeToFEdge[i] * nLayers + il]) > 1e-9)
      //  if(edgeToFEdge[i]>nEdgesSolve_F)
          {
        localSumInt++;
        int edge = edgeToFEdge[i];
        int gEdge = indexToEdgeID[i];
        ID fVertex0 = verticesOnEdge_F[2 * edge] - 1;
        ID fVertex1 = verticesOnEdge_F[2 * edge + 1] - 1;
        ID triaId0 = fVertexToTriangleID[fVertex0];
        ID triaId1 = fVertexToTriangleID[fVertex1];
        ID procTria0 = trianglesProcIds[fVertex0];
        ID procTria1 = trianglesProcIds[fVertex1];
        std::cout << "vs( " << velOnEdges[i * nLayers + il] << ", "
            << u_normal_F[edgeToFEdge[i] * nLayers + il] << ")  ";
        std::cout << "edge: " << edge << ", gEdge: " << gEdge << ", on proc: "
            << edgesProcId[edgeToFEdge[i]];
        if (triaId0 != NotAnId) {
          std::cout << ". first tria0: " << triaId0 << " on proc: "
              << procTria0;
        }
        if (triaId1 != NotAnId) {
          std::cout << ".. second tria0:" << std::endl;
        }
        if ((triaId0 == NotAnId) || (triaId1 == NotAnId)) {
          std::cout << ". and to Tria: " << triaId1 << " on proc: " << procTria1
              << std::endl;
        }

      }

      //localSum = std::max(localSum, std::fabs(velOnEdges[i*nLayers+il] - u_normal_F[edgeToFEdge[i]*nLayers+il]));
    }
  }

  MPI_Allreduce(&localSumInt, &sumInt, 1, MPI_INT, MPI_SUM, comm);

  int localNum(sendEdgesListReversed.size()), num(0);

  MPI_Allreduce(&localNum, &num, 1, MPI_INT, MPI_SUM, comm);

  std::cout << "Edges change: " << sumInt << " " << num << std::endl;

#endif



}


void velocity_solver_export_fo_velocity() {

  if (isDomainEmpty)
    return;

  velocity_solver_export_fo_velocity__(reducedComm);
}

void velocity_solver_finalize() {
  velocity_solver_finalize__();
  delete sendCellsList_F;
  delete recvCellsList_F;
  delete sendEdgesList_F;
  delete recvEdgesList_F;
  delete sendVerticesList_F;
  delete recvVerticesList_F;
}

/*duality:
 *
 *   mpas(F) |  lifev
 *  ---------|---------
 *   cell    |  vertex
 *   vertex  |  triangle
 *   edge    |  edge
 *
 */

void velocity_solver_compute_2d_grid(int const* verticesMask_F) {
  int numProcs, me;

  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &me);
  std::vector<int> partialOffset(numProcs + 1), globalOffsetTriangles(
      numProcs + 1), globalOffsetVertices(numProcs + 1), globalOffsetEdge(
      numProcs + 1);

  std::vector<int> triangleToFVertex;
  triangleToFVertex.reserve(nVertices_F);
  std::vector<int> fVertexToTriangle(nVertices_F, NotAnId);
  bool changed = false;
  for (int i(0); i < nVerticesSolve_F; i++) {
    if ((verticesMask_F[i] & 0x02) && !isGhostTriangle(i)) {
      fVertexToTriangle[i] = triangleToFVertex.size();
      triangleToFVertex.push_back(i);
    }
    changed = changed || (verticesMask_F[i] != mask[i]);
  }

  for (int i(0); i < nVertices_F; i++)
    mask[i] = verticesMask_F[i];

  if (changed)
    std::cout << "mask changed!!" << std::endl;

  if ((me == 0) && (triangleToFVertex.size() == 0))
    for (int i(0); i < nVerticesSolve_F; i++) {
      if (!isGhostTriangle(i)) {
        fVertexToTriangle[i] = triangleToFVertex.size();
        triangleToFVertex.push_back(i);
        break;
      }
    }

  nTriangles = triangleToFVertex.size();

  initialize_iceProblem(nTriangles);

  //Compute the global number of triangles, and the localOffset on the local processor, such that a globalID = localOffset + index
  int localOffset(0);
  nGlobalTriangles = 0;
  computeLocalOffset(nTriangles, localOffset, nGlobalTriangles);

  //Communicate the globalIDs, computed locally, to the other processors.
  indexToTriangleID.resize(nTriangles);

  //To make local, not used
  fVertexToTriangleID.assign(nVertices_F, NotAnId);
  //   std::vector<int> fVertexToTriangleID(nVertices_F, NotAnId);
  for (int index(0); index < nTriangles; index++)
    fVertexToTriangleID[triangleToFVertex[index]] = index + localOffset;

  allToAll(fVertexToTriangleID, sendVerticesList_F, recvVerticesList_F);

  for (int index(0); index < nTriangles; index++)
    indexToTriangleID[index] = fVertexToTriangleID[triangleToFVertex[index]];

  //Compute triangle edges
  std::vector<int> fEdgeToEdge(nEdges_F), edgesToSend, trianglesProcIds(
      nVertices_F);
  getProcIds(trianglesProcIds, recvVerticesList_F);

  int interfaceSize(0);

  std::vector<int> fEdgeToEdgeID(nEdges_F, NotAnId);
  edgesToReceive.clear();
  edgeToFEdge.clear();
  isBoundaryEdge.clear();
  trianglesOnEdge.clear();

  edgesToReceive.reserve(nEdges_F - nEdgesSolve_F);
  edgeToFEdge.reserve(nEdges_F);
  trianglesOnEdge.reserve(nEdges_F * 2);
  edgesToSend.reserve(nEdgesSolve_F);
  isBoundaryEdge.reserve(nEdges_F);

  //first, we compute boundary edges (boundary edges must be the first edges)
  for (int i = 0; i < nEdges_F; i++) {
    ID fVertex1(verticesOnEdge_F[2 * i] - 1), fVertex2(
        verticesOnEdge_F[2 * i + 1] - 1);
    ID triaId_1 = fVertexToTriangleID[fVertex1];
    ID triaId_2 = fVertexToTriangleID[fVertex2];
    bool isboundary = (triaId_1 == NotAnId) || (triaId_2 == NotAnId);

    ID iTria1 = fVertexToTriangle[fVertex1];
    ID iTria2 = fVertexToTriangle[fVertex2];
    if (iTria1 == NotAnId)
      std::swap(iTria1, iTria2);
    bool belongsToLocalTriangle = (iTria1 != NotAnId) || (iTria2 != NotAnId);

    if (belongsToLocalTriangle) {
      if (isboundary) {
        fEdgeToEdge[i] = edgeToFEdge.size();
        edgeToFEdge.push_back(i);
        trianglesOnEdge.push_back(iTria1);
        trianglesOnEdge.push_back(iTria2);
        isBoundaryEdge.push_back(true);
      } else
        interfaceSize += (iTria2 == NotAnId);
    }
  }

  numBoundaryEdges = edgeToFEdge.size();

  //procOnInterfaceEdge contains the pairs <id of interface edge, rank of processor that owns the non local element adjacent to the edge>.
  std::vector < std::pair<int, int> > procOnInterfaceEdge;
  procOnInterfaceEdge.reserve(interfaceSize);

  //then, we compute the other edges
  for (int i = 0; i < nEdges_F; i++) {

    ID fVertex1(verticesOnEdge_F[2 * i] - 1), fVertex2(
        verticesOnEdge_F[2 * i + 1] - 1);
    ID iTria1 = fVertexToTriangle[fVertex1];
    ID iTria2 = fVertexToTriangle[fVertex2];

    ID triaId_1 = fVertexToTriangleID[fVertex1]; //global Triangle
    ID triaId_2 = fVertexToTriangleID[fVertex2]; //global Triangle

    if (iTria1 == NotAnId) {
      std::swap(iTria1, iTria2);
      std::swap(fVertex1, fVertex2);
    }

    bool belongsToAnyTriangle = (triaId_1 != NotAnId) || (triaId_2 != NotAnId);
    bool isboundary = (triaId_1 == NotAnId) || (triaId_2 == NotAnId);
    bool belongsToLocalTriangle = (iTria1 != NotAnId);
    bool isMine = i < nEdgesSolve_F;

    if (belongsToLocalTriangle && !isboundary) {
      fEdgeToEdge[i] = edgeToFEdge.size();
      edgeToFEdge.push_back(i);
      trianglesOnEdge.push_back(iTria1);
      trianglesOnEdge.push_back(iTria2);
      isBoundaryEdge.push_back(false);
      if (iTria2 == NotAnId)
        procOnInterfaceEdge.push_back(
            std::make_pair(fEdgeToEdge[i], trianglesProcIds[fVertex2]));
    }

    if (belongsToAnyTriangle && isMine) {
      edgesToSend.push_back(i);
      if (!belongsToLocalTriangle)
        edgesToReceive.push_back(i);
    }

  }

  //Compute the global number of edges, and the localOffset on the local processor, such that a globalID = localOffset + index
  computeLocalOffset(edgesToSend.size(), localOffset, nGlobalEdges);

  //Communicate the globalIDs, computed locally, to the other processors.
  for (ID index = 0; index < edgesToSend.size(); index++)
    fEdgeToEdgeID[edgesToSend[index]] = index + localOffset;

  allToAll(fEdgeToEdgeID, sendEdgesList_F, recvEdgesList_F);

  nEdges = edgeToFEdge.size();
  indexToEdgeID.resize(nEdges);
  for (int index = 0; index < nEdges; index++)
    indexToEdgeID[index] = fEdgeToEdgeID[edgeToFEdge[index]];

  //Compute vertices:
  std::vector<int> fCellsToSend;
  fCellsToSend.reserve(nCellsSolve_F);

  vertexToFCell.clear();
  vertexToFCell.reserve(nCells_F);

  fCellToVertex.assign(nCells_F, NotAnId);
  std::vector<int> fCellToVertexID(nCells_F, NotAnId);

  fCellsToReceive.clear();

  // if(! isDomainEmpty)
  // {
  fCellsToReceive.reserve(nCells_F - nCellsSolve_F);
  for (int i = 0; i < nCells_F; i++) {
    bool isMine = i < nCellsSolve_F;
    bool belongsToLocalTriangle = false;
    bool belongsToAnyTriangle = false;
    int nEdg = nEdgesOnCells_F[i];
    for (int j = 0; j < nEdg; j++) {
      ID fVertex(verticesOnCell_F[maxNEdgesOnCell_F * i + j] - 1);
      ID iTria = fVertexToTriangle[fVertex];
      ID triaId = fVertexToTriangleID[fVertex];
      belongsToLocalTriangle = belongsToLocalTriangle || (iTria != NotAnId);
      belongsToAnyTriangle = belongsToAnyTriangle || (triaId != NotAnId);
    }

    if (belongsToAnyTriangle && isMine) {
      fCellsToSend.push_back(i);
      if (!belongsToLocalTriangle)
        fCellsToReceive.push_back(i);
    }

    if (belongsToLocalTriangle) {
      fCellToVertex[i] = vertexToFCell.size();
      vertexToFCell.push_back(i);
    }
  }
  //        }

  //Compute the global number of vertices, and the localOffset on the local processor, such that a globalID = localOffset + index
  computeLocalOffset(fCellsToSend.size(), localOffset, nGlobalVertices);

  //Communicate the globalIDs, computed locally, to the other processors.
  for (int index = 0; index < int(fCellsToSend.size()); index++)
    fCellToVertexID[fCellsToSend[index]] = index + localOffset;

  allToAll(fCellToVertexID, sendCellsList_F, recvCellsList_F);

  nVertices = vertexToFCell.size();
  std::cout << "\n nvertices: " << nVertices << " " << nGlobalVertices << "\n"
      << std::endl;
  indexToVertexID.resize(nVertices);
  for (int index = 0; index < nVertices; index++)
    indexToVertexID[index] = fCellToVertexID[vertexToFCell[index]];

  createReverseCellsExchangeLists(sendCellsListReversed, recvCellsListReversed,
      fVertexToTriangleID, fCellToVertexID);

  //construct the local vector vertices on triangles making sure the area is positive
  verticesOnTria.resize(nTriangles * 3);
  double x[3], y[3], z[3];
  for (int index = 0; index < nTriangles; index++) {
    int iTria = triangleToFVertex[index];

    for (int j = 0; j < 3; j++) {
      int iCell = cellsOnVertex_F[3 * iTria + j] - 1;
      verticesOnTria[3 * index + j] = fCellToVertex[iCell];
      x[j] = xCell_F[iCell];
      y[j] = yCell_F[iCell];
      //  z[j] = zCell_F[iCell];
    }
    if (signedTriangleArea(x, y) < 0)
      std::swap(verticesOnTria[3 * index + 1], verticesOnTria[3 * index + 2]);
  }

  //construct the local vector vertices on edges
  trianglesPositionsOnEdge.resize(2 * nEdges);
  isVertexBoundary.assign(nVertices, false);

  verticesOnEdge.resize(2 * nEdges);

  //contains the local id of a triangle and the global id of the edges of the triangle.
  //dataForGhostTria[4*i] contains the triangle id
  //dataForGhostTria[4*i+1+k] contains the global id of the edge (at position k = 0,1,2) of the triangle.
  //Possible Optimization: for our purposes it would be enough to store two of the three edges of a triangle.
  std::vector<int> dataForGhostTria(nVertices_F * 4, NotAnId);

  //*
  for (int iV = 0; iV < nVertices; iV++) {
    int fCell = vertexToFCell[iV];
    int nEdg = nEdgesOnCells_F[fCell];
    int j = 0;
    bool isBoundary;
    do {
      int fVertex = verticesOnCell_F[maxNEdgesOnCell_F * fCell + j++] - 1;
      isBoundary = !(verticesMask_F[fVertex] & 0x02);
    } while ((j < nEdg) && (!isBoundary));
    isVertexBoundary[iV] = isBoundary;
  }
  /*/
   for(int index=0; index<numBoundaryEdges; index++)
   {
   int fEdge = edgeToFEdge[index];
   int v1 = fCellToVertex[cellsOnEdge_F[2*fEdge]-1];
   int v2 = fCellToVertex[cellsOnEdge_F[2*fEdge+1]-1];
   isVertexBoundary[ v1 ] = isVertexBoundary[ v2 ] = true;
   }
   //*/
  //computing the position and local ids of the triangles adjacent to edges.
  //in the case an adjacent triangle is not local, the position and id will be NotAnId.
  //We will get the needed information about the non local edge leter on,
  //for this purpose we fill the vector dataForGhostTria.
  for (int index = 0; index < nEdges; index++) {
    int p1, p2, v1, v2, v3, t, position;
    int fEdge = edgeToFEdge[index];
    int fCell1 = cellsOnEdge_F[2 * fEdge] - 1;
    int fCell2 = cellsOnEdge_F[2 * fEdge + 1] - 1;
    verticesOnEdge[2 * index] = p1 = fCellToVertex[fCell1];
    verticesOnEdge[2 * index + 1] = p2 = fCellToVertex[fCell2];

    for (int k = 0; k < 2 && (t = trianglesOnEdge[2 * index + k]) != NotAnId;
        k++) {
      v1 = verticesOnTria[3 * t];
      v2 = verticesOnTria[3 * t + 1];
      v3 = verticesOnTria[3 * t + 2];
      position = (((p1 == v2) && (p2 == v3)) || ((p1 == v3) && (p2 == v2)))
          + 2 * (((p1 == v1) && (p2 == v3)) || ((p1 == v3) && (p2 == v1)));
      trianglesPositionsOnEdge[2 * index + k] = position;
      int dataIndex = 4 * triangleToFVertex[t];

      dataForGhostTria[dataIndex] = t;
      dataForGhostTria[dataIndex + position + 1] = indexToEdgeID[index];
    }
  }

  createReverseEdgesExchangeLists(sendEdgesListReversed, recvEdgesListReversed,
      fVertexToTriangleID, fEdgeToEdgeID);

  //send the information about ghost elements.
  allToAll(dataForGhostTria, sendVerticesList_F, recvVerticesList_F, 4);

  //retrieving the position, local id and owner processor of non local triangles adjacent to interface edges.
  for (int index = 0; index < (int) procOnInterfaceEdge.size(); index++) {
    int iEdge = procOnInterfaceEdge[index].first;
    int fEdge = edgeToFEdge[iEdge];
    ID fVertex1(verticesOnEdge_F[2 * fEdge] - 1), fVertex2(
        verticesOnEdge_F[2 * fEdge + 1] - 1);
    if ((ID) fVertexToTriangle[fVertex1] == NotAnId)
      fVertex2 = fVertex1;
    int dataIndex = 4 * fVertex2;
    int edgeID = indexToEdgeID[iEdge];
    int position = (dataForGhostTria[dataIndex + 2] == edgeID)
        + 2 * (dataForGhostTria[dataIndex + 3] == edgeID);
    trianglesOnEdge[2 * iEdge + 1] = dataForGhostTria[dataIndex];
    trianglesPositionsOnEdge[2 * iEdge + 1] = position;
  }

  if (isDomainEmpty)
    return;

#ifdef LIFEV
    std::vector<double> verticesCoords(3 * nVertices);

    for (int index = 0; index < nVertices; index++) {
      int iCell = vertexToFCell[index];
      verticesCoords[index * 3] = xCell_F[iCell] / unit_length;
      verticesCoords[index * 3 + 1] = yCell_F[iCell] / unit_length;
      verticesCoords[index * 3 + 2] = zCell_F[iCell] / unit_length;
    }

    velocity_solver_compute_2d_grid__(nGlobalTriangles,
        nGlobalVertices, nGlobalEdges, indexToVertexID, verticesCoords,
        isVertexBoundary, verticesOnTria,isBoundaryEdge, trianglesOnEdge,
        trianglesPositionsOnEdge, verticesOnEdge, indexToEdgeID,
        indexToTriangleID, procOnInterfaceEdge );
#else
    velocity_solver_compute_2d_grid__(reducedComm);
#endif

  /*

   //initialize the mesh
   iceProblemPtr->mesh2DPtr.reset (new RegionMesh<LinearTriangle>() );

   //construct the mesh nodes
   constructNodes ( * (iceProblemPtr->mesh2DPtr), indexToVertexID, verticesCoords, isVertexBoundary, nGlobalVertices, 3);

   //construct the mesh elements
   constructElements ( * (iceProblemPtr->mesh2DPtr), indexToTriangleID, verticesOnTria, nGlobalTriangles);

   //construct the mesh facets
   constructFacets ( * (iceProblemPtr->mesh2DPtr), isBoundaryEdge, trianglesOnEdge, trianglesPositionsOnEdge, verticesOnEdge, indexToEdgeID, procOnInterfaceEdge, nGlobalEdges, 3);

   Switch sw;
   std::vector<bool> elSign;
   checkVolumes ( * (iceProblemPtr->mesh2DPtr), elSign, sw );
   */
}

void velocity_solver_extrude_3d_grid(double const* levelsRatio_F,
    double const* lowerSurface_F, double const* thickness_F) {

  if (isDomainEmpty)
    return;

  layersRatio.resize(nLayers);
  // !!Indexing of layers is reversed
  for (int i = 0; i < nLayers; i++)
    layersRatio[i] = levelsRatio_F[nLayers - 1 - i];
  //std::copy(levelsRatio_F, levelsRatio_F+nLayers, layersRatio.begin());

  levelsNormalizedThickness.resize(nLayers + 1);

  levelsNormalizedThickness[0] = 0;
  for (int i = 0; i < nLayers; i++)
    levelsNormalizedThickness[i + 1] = levelsNormalizedThickness[i]
        + layersRatio[i];

  std::vector<int> mpasIndexToVertexID(nVertices);
  for (int i = 0; i < nVertices; i++)
    mpasIndexToVertexID[i] = indexToCellID_F[vertexToFCell[i]];

  //construct the local vector of coordinates
  std::vector<double> verticesCoords(3 * nVertices);

  for (int index = 0; index < nVertices; index++) {
    int iCell = vertexToFCell[index];
    verticesCoords[index * 3] = xCell_F[iCell] / unit_length;
    verticesCoords[index * 3 + 1] = yCell_F[iCell] / unit_length;
    verticesCoords[index * 3 + 2] = zCell_F[iCell] / unit_length;
  }

  velocity_solver_extrude_3d_grid__(nLayers, nGlobalTriangles, nGlobalVertices,
      nGlobalEdges, Ordering, reducedComm, indexToVertexID, mpasIndexToVertexID,
      verticesCoords, isVertexBoundary, verticesOnTria, isBoundaryEdge,
      trianglesOnEdge, trianglesPositionsOnEdge, verticesOnEdge, indexToEdgeID,
      indexToTriangleID);
}
}

void get_tetraP1_velocity_on_FEdges(double * uNormal,
    const std::vector<double>& velocityOnVertices,
    const std::vector<int>& edgeToFEdge,
    const std::vector<int>& mpasIndexToVertexID) {

  int columnShift = (Ordering == 1) ? 1 : nVertices;
  int layerShift = (Ordering == 0) ? 1 : nLayers + 1;

  UInt nPoints3D = nVertices * (nLayers + 1);

  //the velocity on boundary edges is set to zero by construction
  for (int i = numBoundaryEdges; i < nEdges; i++) {
    ID lId0 = verticesOnEdge[2 * i];
    ID lId1 = verticesOnEdge[2 * i + 1];
    int iCell0 = vertexToFCell[lId0];
    int iCell1 = vertexToFCell[lId1];

    double nx = xCell_F[iCell1] - xCell_F[iCell0];
    double ny = yCell_F[iCell1] - yCell_F[iCell0];
    double n = sqrt(nx * nx + ny * ny);
    nx /= n;
    ny /= n;

    int iEdge = edgeToFEdge[i];
    //prism lateral face is splitted with the diagonal that start from p0 (by construction).
    double coeff0 =
        (mpasIndexToVertexID[lId0] > mpasIndexToVertexID[lId1]) ?
            1. / 3. : 1. / 6.;
    double coeff1 = 0.5 - coeff0;
    for (int il = 0; il < nLayers; il++) {
      int ilReversed = nLayers - il - 1;
      ID lId3D1 = il * columnShift + layerShift * lId0;
      ID lId3D2 = il * columnShift + layerShift * lId1;
      //not accurate
      uNormal[iEdge * nLayers + ilReversed] = nx
          * (coeff1 * velocityOnVertices[lId3D1]
              + coeff0 * velocityOnVertices[lId3D2]
              + coeff0 * velocityOnVertices[lId3D1 + columnShift]
              + coeff1 * velocityOnVertices[lId3D2 + columnShift]);

      lId3D1 += nPoints3D;
      lId3D2 += nPoints3D;
      uNormal[iEdge * nLayers + ilReversed] += ny
          * (coeff1 * velocityOnVertices[lId3D1]
              + coeff0 * velocityOnVertices[lId3D2]
              + coeff0 * velocityOnVertices[lId3D1 + columnShift]
              + coeff1 * velocityOnVertices[lId3D2 + columnShift]);
    }
  }

}

void get_prism_velocity_on_FEdges(double * uNormal,
    const std::vector<double>& velocityOnVertices,
    const std::vector<int>& edgeToFEdge) {

  //fix this
  int columnShift = (Ordering == 1) ? 1 : nVertices;
  int layerShift = (Ordering == 0) ? 1 : nLayers + 1;

  UInt nPoints3D = nVertices * (nLayers + 1);

  for (int i = numBoundaryEdges; i < nEdges; i++) {
    ID lId0 = verticesOnEdge[2 * i];
    ID lId1 = verticesOnEdge[2 * i + 1];
    int iCell0 = vertexToFCell[lId0];
    int iCell1 = vertexToFCell[lId1];

    double nx = xCell_F[iCell1] - xCell_F[iCell0];
    double ny = yCell_F[iCell1] - yCell_F[iCell0];
    double n = sqrt(nx * nx + ny * ny);
    nx /= n;
    ny /= n;

    ID iEdge = edgeToFEdge[i];

    for (int il = 0; il < nLayers; il++) {
      int ilReversed = nLayers - il - 1;
      ID lId3D1 = il * columnShift + layerShift * lId0;
      ID lId3D2 = il * columnShift + layerShift * lId1;
      //not accurate
      uNormal[iEdge * nLayers + ilReversed] = 0.25 * nx
          * (velocityOnVertices[lId3D1] + velocityOnVertices[lId3D2]
              + velocityOnVertices[lId3D1 + columnShift]
              + velocityOnVertices[lId3D2 + columnShift]);

      lId3D1 += nPoints3D;
      lId3D2 += nPoints3D;
      uNormal[iEdge * nLayers + ilReversed] += 0.25 * ny
          * (velocityOnVertices[lId3D1] + velocityOnVertices[lId3D2]
              + velocityOnVertices[lId3D1 + columnShift]
              + velocityOnVertices[lId3D2 + columnShift]);
    }
  }

}

void mapVerticesToCells(const std::vector<double>& velocityOnVertices,
    double* velocityOnCells, int fieldDim, int numLayers, int ordering) {
  int lVertexColumnShift = (ordering == 1) ? 1 : nVertices;
  int vertexLayerShift = (ordering == 0) ? 1 : numLayers + 1;

  int nVertices3D = nVertices * (numLayers + 1);
  for (UInt j = 0; j < nVertices3D; ++j) {
    int ib = (ordering == 0) * (j % lVertexColumnShift)
        + (ordering == 1) * (j / vertexLayerShift);
    int il = (ordering == 0) * (j / lVertexColumnShift)
        + (ordering == 1) * (j % vertexLayerShift);

    int iCell = vertexToFCell[ib];
    int cellIndex = iCell * (numLayers + 1) + il;
    int vertexIndex = j;
    for (int dim = 0; dim < fieldDim; dim++) {
      velocityOnCells[cellIndex] = velocityOnVertices[vertexIndex];
      cellIndex += nCells_F * (numLayers + 1);
      vertexIndex += nVertices3D;
    }
  }

  for (int dim = 0; dim < fieldDim; dim++) {
    allToAll(&velocityOnCells[dim * nCells_F * (numLayers + 1)],
        &sendCellsListReversed, &recvCellsListReversed, (numLayers + 1));
    allToAll(&velocityOnCells[dim * nCells_F * (numLayers + 1)],
        sendCellsList_F, recvCellsList_F, (numLayers + 1));
  }
}

void createReverseCellsExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& fVertexToTriangleID,
    const std::vector<int>& fCellToVertexID) {
  sendListReverse_F.clear();
  receiveListReverse_F.clear();
  //std::map<int, std::vector<int> > sendMap;
  std::map<int, std::map<int, int> > sendMap, receiveMap;
  std::vector<int> cellsProcId(nCells_F), trianglesProcIds(nVertices_F);
  getProcIds(cellsProcId, recvCellsList_F);
  getProcIds(trianglesProcIds, recvVerticesList_F);

  //std::cout << "SendList " ;
  for (int i = 0; i < nVertices; i++) {
    int iCell = vertexToFCell[i];
    if (iCell < nCellsSolve_F)
      continue;
    bool belongToTriaOnSameProc = false;
    int j(0);
    int nEdg = nEdgesOnCells_F[iCell];
    do {
      ID fVertex(verticesOnCell_F[maxNEdgesOnCell_F * iCell + j] - 1);
      ID triaId = fVertexToTriangleID[fVertex];
      belongToTriaOnSameProc = (triaId != NotAnId)
          && (trianglesProcIds[fVertex] == cellsProcId[iCell]);
    } while ((belongToTriaOnSameProc == false) && (++j < nEdg));
    if (!belongToTriaOnSameProc) {
      sendMap[cellsProcId[iCell]].insert(
          std::make_pair(fCellToVertexID[iCell], iCell));
      // std::cout<< "(" << cellsProcId[iCell] << "," << iCell << ") ";
    }

  }
  //std::cout <<std::endl;

  for (std::map<int, std::map<int, int> >::const_iterator it = sendMap.begin();
      it != sendMap.end(); it++) {
    std::vector<int> sendVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      sendVec[i++] = iter->second;
    sendListReverse_F.push_back(
        exchange(it->first, &sendVec[0], &sendVec[0] + sendVec.size()));
  }

  //std::cout << "ReceiveList " ;
  for (UInt i = 0; i < fCellsToReceive.size(); i++) {
    int iCell = fCellsToReceive[i];
    int nEdg = nEdgesOnCells_F[iCell];
    for (int j = 0; j < nEdg; j++) {
      ID fVertex(verticesOnCell_F[maxNEdgesOnCell_F * iCell + j] - 1);
      ID triaId = fVertexToTriangleID[fVertex];
      if (triaId != NotAnId) {
        receiveMap[trianglesProcIds[fVertex]].insert(
            std::make_pair(fCellToVertexID[iCell], iCell));
        // std::cout<< "(" << trianglesProcIds[fVertex] << "," << iCell << ") ";
      }
    }
  }
  //std::cout <<std::endl;

  for (std::map<int, std::map<int, int> >::const_iterator it =
      receiveMap.begin(); it != receiveMap.end(); it++) {
    std::vector<int> receiveVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      receiveVec[i++] = iter->second;
    receiveListReverse_F.push_back(
        exchange(it->first, &receiveVec[0],
            &receiveVec[0] + receiveVec.size()));
  }
}

void createReverseEdgesExchangeLists(exchangeList_Type& sendListReverse_F,
    exchangeList_Type& receiveListReverse_F,
    const std::vector<int>& fVertexToTriangleID,
    const std::vector<int>& fEdgeToEdgeID) {
  sendListReverse_F.clear();
  receiveListReverse_F.clear();
  //std::map<int, std::vector<int> > sendMap;
  std::map<int, std::map<int, int> > sendMap, receiveMap;
  std::vector<int> edgesProcId(nEdges_F), trianglesProcIds(nVertices_F);
  getProcIds(edgesProcId, recvEdgesList_F);
  getProcIds(trianglesProcIds, recvVerticesList_F);

  //std::cout << "EdgesSendList " ;
  for (int i = 0; i < nEdges; i++) {
    int iEdge = edgeToFEdge[i];
    if (iEdge < nEdgesSolve_F)
      continue;
    bool belongToTriaOnSameProc = false;
    int j(0);
    do {
      ID fVertex(verticesOnEdge_F[2 * iEdge + j] - 1);
      ID triaId = fVertexToTriangleID[fVertex];
      belongToTriaOnSameProc = (triaId != NotAnId)
          && (trianglesProcIds[fVertex] == edgesProcId[iEdge]);
    } while ((belongToTriaOnSameProc == false) && (++j < 2));
    if (!belongToTriaOnSameProc) {
      sendMap[edgesProcId[iEdge]].insert(
          std::make_pair(fEdgeToEdgeID[iEdge], iEdge));
      //std::cout<< "(" << edgesProcId[iEdge] << "," << iEdge << ") ";
    }
  }
  //std::cout <<std::endl;

  for (std::map<int, std::map<int, int> >::const_iterator it = sendMap.begin();
      it != sendMap.end(); it++) {
    std::vector<int> sendVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      sendVec[i++] = iter->second;
    sendListReverse_F.push_back(
        exchange(it->first, &sendVec[0], &sendVec[0] + sendVec.size()));
  }

  //std::cout << "EdgesReceiveList " ;
  for (UInt i = 0; i < edgesToReceive.size(); i++) {
    int iEdge = edgesToReceive[i];
    for (int j = 0; j < 2; j++) {
      ID fVertex(verticesOnEdge_F[2 * iEdge + j] - 1);
      ID triaId = fVertexToTriangleID[fVertex];
      if (triaId != NotAnId) {
        receiveMap[trianglesProcIds[fVertex]].insert(
            std::make_pair(fEdgeToEdgeID[iEdge], iEdge));
        //  std::cout<< "(" << trianglesProcIds[fVertex] << "," << iEdge << ") ";
      }
    }
  }
  // std::cout <<std::endl;

  for (std::map<int, std::map<int, int> >::const_iterator it =
      receiveMap.begin(); it != receiveMap.end(); it++) {
    std::vector<int> receiveVec(it->second.size());
    int i = 0;
    for (std::map<int, int>::const_iterator iter = it->second.begin();
        iter != it->second.end(); iter++)
      receiveVec[i++] = iter->second;
    receiveListReverse_F.push_back(
        exchange(it->first, &receiveVec[0],
            &receiveVec[0] + receiveVec.size()));
  }
}

void mapCellsToVertices(const std::vector<double>& velocityOnCells,
    std::vector<double>& velocityOnVertices, int fieldDim, int numLayers,
    int ordering) {
  int lVertexColumnShift = (ordering == 1) ? 1 : nVertices;
  int vertexLayerShift = (ordering == 0) ? 1 : numLayers + 1;

  int nVertices3D = nVertices * (numLayers + 1);
  for (UInt j = 0; j < nVertices3D; ++j) {
    int ib = (ordering == 0) * (j % lVertexColumnShift)
        + (ordering == 1) * (j / vertexLayerShift);
    int il = (ordering == 0) * (j / lVertexColumnShift)
        + (ordering == 1) * (j % vertexLayerShift);

    int iCell = vertexToFCell[ib];
    int cellIndex = iCell * (numLayers + 1) + il;
    int vertexIndex = j;
    for (int dim = 0; dim < fieldDim; dim++) {
      velocityOnVertices[vertexIndex] = velocityOnCells[cellIndex];
      cellIndex += nCells_F * (numLayers + 1);
      vertexIndex += nVertices3D;
    }
  }
}

bool isGhostTriangle(int i, double relTol) {
  double x[3], y[3], area;

  for (int j = 0; j < 3; j++) {
    int iCell = cellsOnVertex_F[3 * i + j] - 1;
    x[j] = xCell_F[iCell];
    y[j] = yCell_F[iCell];
  }

  area = std::fabs(signedTriangleArea(x, y));
  return false; //(std::fabs(areaTriangle_F[i]-area)/areaTriangle_F[i] > relTol);
}

double signedTriangleArea(const double* x, const double* y) {
  double u[2] = { x[1] - x[0], y[1] - y[0] };
  double v[2] = { x[2] - x[0], y[2] - y[0] };

  return 0.5 * (u[0] * v[1] - u[1] * v[0]);
}

double signedTriangleAreaOnSphere(const double* x, const double* y,
    const double *z) {
  double u[3] = { x[1] - x[0], y[1] - y[0], z[1] - z[0] };
  double v[3] = { x[2] - x[0], y[2] - y[0], z[2] - z[0] };

  double crossProduct[3] = { u[1] * v[2] - u[2] * v[1], u[2] * v[0]
      - u[0] * v[2], u[0] * v[1] - u[1] * v[0] };
  double area = 0.5
      * std::sqrt(
          crossProduct[0] * crossProduct[0] + crossProduct[1] * crossProduct[1]
              + crossProduct[2] * crossProduct[2]);
  return
      (crossProduct[0] * x[0] + crossProduct[1] * y[0] + crossProduct[2] * z[0]
          > 0) ? area : -area;
}

//TO BE FIXED, Access To verticesOnCell_F is not correct
void extendMaskByOneLayer(int const* verticesMask_F,
    std::vector<int>& extendedFVerticesMask) {
  extendedFVerticesMask.resize(nVertices_F);
  extendedFVerticesMask.assign(&verticesMask_F[0],
      &verticesMask_F[0] + nVertices_F);
  for (int i = 0; i < nCells_F; i++) {
    bool belongsToMarkedTriangle = false;
    int nEdg = nEdgesOnCells_F[i];
    for (UInt k = 0; k < nEdg && !belongsToMarkedTriangle; k++)
      belongsToMarkedTriangle = belongsToMarkedTriangle
          || verticesMask_F[verticesOnCell_F[maxNEdgesOnCell_F * i + k] - 1];
    if (belongsToMarkedTriangle)
      for (UInt k = 0; k < nEdg; k++) {
        ID fVertex(verticesOnCell_F[maxNEdgesOnCell_F * i + k] - 1);
        extendedFVerticesMask[fVertex] = !isGhostTriangle(fVertex);
      }
  }
}

void import2DFields(double const * lowerSurface_F, double const * thickness_F,
    double const * beta_F, double eps) {
  elevationData.assign(nVertices, 1e10);
  thicknessData.assign(nVertices, 1e10);
  std::map<int, int> bdExtensionMap;
  if (beta_F != 0)
    betaData.assign(nVertices, 1e10);

  for (int iV = 0; iV < nVertices; iV++) {
    if (isVertexBoundary[iV]) {
      int c;
      int fCell = vertexToFCell[iV];
      int nEdg = nEdgesOnCells_F[fCell];
      for (int j = 0; j < nEdg; j++) {
        int fEdge = edgesOnCell_F[maxNEdgesOnCell_F * fCell + j] - 1;
        bool keep = (mask[verticesOnEdge_F[2 * fEdge] - 1] & 0x02)
            && (mask[verticesOnEdge_F[2 * fEdge + 1] - 1] & 0x02);
        if (!keep)
          continue;

        int c0 = cellsOnEdge_F[2 * fEdge] - 1;
        int c1 = cellsOnEdge_F[2 * fEdge + 1] - 1;
        c = (fCellToVertex[c0] == iV) ? c1 : c0;
        double elev = thickness_F[c] + lowerSurface_F[c]; // - 1e-8*std::sqrt(pow(xCell_F[c0],2)+std::pow(yCell_F[c0],2));
        if (elevationData[iV] > elev) {
          elevationData[iV] = elev;
          bdExtensionMap[iV] = c;
        }
      }
    }
  }

  for (std::map<int, int>::iterator it = bdExtensionMap.begin();
      it != bdExtensionMap.end(); ++it) {
    int iv = it->first;
    int ic = it->second;
    thicknessData[iv] = thickness_F[ic] / unit_length + eps; //- 1e-8*std::sqrt(pow(xCell_F[ic],2)+std::pow(yCell_F[ic],2));
    elevationData[iv] = std::max(
        thicknessData[iv] + lowerSurface_F[ic] / unit_length,
        118. / 1028. * thicknessData[iv]);
    if (beta_F != 0)
      //  betaData[ iv ] = beta_F[ic]/unit_length;
      betaData[iv] =
          (lowerSurface_F[ic] > -910. / 1028. * thickness_F[ic]) ?
              beta_F[ic] / unit_length : 0;
  }

  for (int index = 0; index < nVertices; index++) {
    int iCell = vertexToFCell[index];

    if (!isVertexBoundary[index]) {
      thicknessData[index] = thickness_F[iCell] / unit_length + eps;
      elevationData[index] = std::max(
          (lowerSurface_F[iCell] / unit_length) + thicknessData[index],
          118. / 1028. * thicknessData[index]);
    }
  }

  if (beta_F != 0) {
    for (int index = 0; index < nVertices; index++) {
      int iCell = vertexToFCell[index];

      if (!isVertexBoundary[index])
        //betaData[ index ] = beta_F[iCell]/unit_length;
        betaData[index] =
            (lowerSurface_F[iCell] > -910. / 1028. * thickness_F[iCell]) ?
                beta_F[iCell] / unit_length : 0;
    }
  }

}

void importP0Temperature(double const * temperature_F) {
  int lElemColumnShift = (Ordering == 1) ? 3 : 3 * indexToTriangleID.size();
  int elemLayerShift = (Ordering == 0) ? 3 : 3 * nLayers;
  temperatureOnTetra.resize(3 * nLayers * indexToTriangleID.size());
  for (int index = 0; index < nTriangles; index++) {
    for (int il = 0; il < nLayers; il++) {
      double temperature = 0;
      int ilReversed = nLayers - il - 1;
      int nPoints = 0;
      for (int iVertex = 0; iVertex < 3; iVertex++) {
        int v = verticesOnTria[iVertex + 3 * index];
        if (!isVertexBoundary[v]) {
          int iCell = vertexToFCell[v];
          temperature += temperature_F[iCell * nLayers + ilReversed];
          nPoints++;
        }
      }
      if (nPoints == 0)
        temperature = T0;
      else
        temperature = temperature / nPoints + T0;
      for (int k = 0; k < 3; k++)
        temperatureOnTetra[index * elemLayerShift + il * lElemColumnShift + k] =
            temperature;
    }
  }

}

void createReducedMPI(int nLocalEntities, MPI_Comm& reduced_comm_id) {
  int numProcs, me;
  MPI_Group world_group_id, reduced_group_id;
  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &me);
  std::vector<int> haveElements(numProcs);
  int nonEmpty = int(nLocalEntities > 0);
  MPI_Allgather(&nonEmpty, 1, MPI_INT, &haveElements[0], 1, MPI_INT, comm);
  std::vector<int> ranks;
  for (int i = 0; i < numProcs; i++) {
    if (haveElements[i])
      ranks.push_back(i);
  }

  MPI_Comm_group(comm, &world_group_id);
  MPI_Group_incl(world_group_id, ranks.size(), &ranks[0], &reduced_group_id);
  MPI_Comm_create(comm, reduced_group_id, &reduced_comm_id);
}

void computeLocalOffset(int nLocalEntities, int& localOffset,
    int& nGlobalEntities) {
  int numProcs, me;
  MPI_Comm_size(comm, &numProcs);
  MPI_Comm_rank(comm, &me);
  std::vector<int> offsetVec(numProcs);

  MPI_Allgather(&nLocalEntities, 1, MPI_INT, &offsetVec[0], 1, MPI_INT, comm);

  localOffset = 0;
  for (int i = 0; i < me; i++)
    localOffset += offsetVec[i];

  nGlobalEntities = localOffset;
  for (int i = me; i < numProcs; i++)
    nGlobalEntities += offsetVec[i];
}

void getProcIds(std::vector<int>& field, int const * recvArray) {
  int me;
  MPI_Comm_rank(comm, &me);
  field.assign(field.size(), me);

  //unpack recvArray and set the proc rank into filed
  for (int i(1), procID, size; i < recvArray[0]; i += size) {
    procID = recvArray[i++];
    size = recvArray[i++];
    if (procID == me)
      continue;
    for (int k = i; k < i + size; k++)
      field[recvArray[k]] = procID;
  }
}

void getProcIds(std::vector<int>& field, exchangeList_Type const * recvList) {
  int me;
  MPI_Comm_rank(comm, &me);
  field.assign(field.size(), me);
  exchangeList_Type::const_iterator it;

  for (it = recvList->begin(); it != recvList->end(); ++it) {
    if (it->procID == me)
      continue;
    for (int k = 0; k < (int) it->vec.size(); k++)
      field[it->vec[k]] = it->procID;
  }
}

exchangeList_Type unpackMpiArray(int const * array) {
  exchangeList_Type list;
  for (int i(1), procID, size; i < array[0]; i += size) {
    procID = array[i++];
    size = array[i++];
    list.push_back(exchange(procID, &array[i], &array[i + size]));
  }
  return list;
}

void allToAll(std::vector<int>& field, int const * sendArray,
    int const * recvArray, int fieldDim) {
  exchangeList_Type sendList, recvList;

  //unpack sendArray and build the sendList class
  for (int i(1), procID, size; i < sendArray[0]; i += size) {
    procID = sendArray[i++];
    size = sendArray[i++];
    sendList.push_back(
        exchange(procID, &sendArray[i], &sendArray[i + size], fieldDim));
  }

  //unpack recvArray and build the recvList class
  for (int i(1), procID, size; i < recvArray[0]; i += size) {
    procID = recvArray[i++];
    size = recvArray[i++];
    recvList.push_back(
        exchange(procID, &recvArray[i], &recvArray[i + size], fieldDim));
  }

  int me;
  MPI_Comm_rank(comm, &me);

  exchangeList_Type::iterator it;
  for (it = recvList.begin(); it != recvList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Irecv(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID,
        it->procID, comm, &it->reqID);
  }

  for (it = sendList.begin(); it != sendList.end(); ++it) {
    if (it->procID == me)
      continue;
    for (ID i = 0; i < it->vec.size(); i++)
      for (int iComp = 0; iComp < fieldDim; iComp++)
        it->buffer[fieldDim * i + iComp] = field[fieldDim * it->vec[i] + iComp];

    MPI_Isend(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID, me,
        comm, &it->reqID);
  }

  for (it = recvList.begin(); it != recvList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

    for (int i = 0; i < int(it->vec.size()); i++)
      for (int iComp = 0; iComp < fieldDim; iComp++)
        field[fieldDim * it->vec[i] + iComp] = it->buffer[fieldDim * i + iComp];
  }

  for (it = sendList.begin(); it != sendList.end(); ++it) {
    if (it->procID == me)
      continue;
    MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
  }
}

void allToAll(std::vector<int>& field, exchangeList_Type const * sendList,
    exchangeList_Type const * recvList, int fieldDim) {
  int me;
  MPI_Comm_rank(comm, &me);

  for (int iComp = 0; iComp < fieldDim; iComp++) {
    exchangeList_Type::const_iterator it;
    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Irecv(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID,
          it->procID, comm, &it->reqID);
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      for (ID i = 0; i < it->vec.size(); i++)
        it->buffer[i] = field[fieldDim * it->vec[i] + iComp];

      MPI_Isend(&(it->buffer[0]), it->buffer.size(), MPI_INT, it->procID, me,
          comm, &it->reqID);
    }

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

      for (int i = 0; i < int(it->vec.size()); i++)
        field[fieldDim * it->vec[i] + iComp] = it->buffer[i];
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
    }
  }
}

void allToAll(double* field, exchangeList_Type const * sendList,
    exchangeList_Type const * recvList, int fieldDim) {
  int me;
  MPI_Comm_rank(comm, &me);

  for (int iComp = 0; iComp < fieldDim; iComp++) {
    exchangeList_Type::const_iterator it;

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Irecv(&(it->doubleBuffer[0]), it->doubleBuffer.size(), MPI_DOUBLE,
          it->procID, it->procID, comm, &it->reqID);
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      for (ID i = 0; i < it->vec.size(); i++)
        it->doubleBuffer[i] = field[fieldDim * it->vec[i] + iComp];

      MPI_Isend(&(it->doubleBuffer[0]), it->doubleBuffer.size(), MPI_DOUBLE,
          it->procID, me, comm, &it->reqID);
    }

    for (it = recvList->begin(); it != recvList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);

      for (int i = 0; i < int(it->vec.size()); i++)
        field[fieldDim * it->vec[i] + iComp] = it->doubleBuffer[i];
    }

    for (it = sendList->begin(); it != sendList->end(); ++it) {
      if (it->procID == me)
        continue;
      MPI_Wait(&it->reqID, MPI_STATUS_IGNORE);
    }
  }
}

int initialize_iceProblem(int nTriangles) {
  bool keep_proc = nTriangles > 0;

  createReducedMPI(keep_proc, reducedComm);

  isDomainEmpty = !keep_proc;


#ifdef LIFEV
if(!isDomainEmpty) {
  velocity_solver_initialize_iceProblem__(keep_proc, reducedComm);
}
#endif

  // initialize ice problem pointer
  if (keep_proc) {
    std::cout << nTriangles
        << " elements of the triangular grid are stored on this processor"
        << std::endl;
  } else {
    std::cout
        << "No elements of the triangular grid are stored on this processor"
        << std::endl;
  }

  return 0;
}


int prismType(long long int const* prismVertexMpasIds, int& minIndex)
{
  int PrismVerticesMap[6][6] = {{0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1}, {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0}};
  minIndex = std::min_element (prismVertexMpasIds, prismVertexMpasIds + 3) - prismVertexMpasIds;

  int v1 (prismVertexMpasIds[PrismVerticesMap[minIndex][1]]);
  int v2 (prismVertexMpasIds[PrismVerticesMap[minIndex][2]]);

  return v1  > v2;
}

 void tetrasFromPrismStructured (long long int const* prismVertexMpasIds, long long int const* prismVertexGIds, long long int tetrasIdsOnPrism[][4])
  {
      int PrismVerticesMap[6][6] = {{0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1}, {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0}};

      int tetraOfPrism[2][3][4] = {{{0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3}}, {{0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3}}};

      int tetraAdjacentToPrismLateralFace[2][3][2] = {{{1, 2}, {0, 1}, {0, 2}}, {{0, 2}, {0, 1}, {1, 2}}};
      int tetraFaceIdOnPrismLateralFace[2][3][2] = {{{0, 0}, {1, 1}, {2, 2}}, {{0, 0}, {1, 1}, {2, 2}}};
      int tetraAdjacentToBottomFace = 0; //does not depend on type;
      int tetraAdjacentToUpperFace = 2; //does not depend on type;
      int tetraFaceIdOnBottomFace = 3; //does not depend on type;
      int tetraFaceIdOnUpperFace = 0; //does not depend on type;

      int minIndex;
      int prismT = prismType(prismVertexMpasIds, minIndex);

      long long int reorderedPrismLIds[6];

      for (int ii = 0; ii < 6; ii++)
      {
          reorderedPrismLIds[ii] = prismVertexGIds[PrismVerticesMap[minIndex][ii]];
      }

      for (int iTetra = 0; iTetra < 3; iTetra++)
          for (int iVertex = 0; iVertex < 4; iVertex++)
          {
              tetrasIdsOnPrism[iTetra][iVertex] = reorderedPrismLIds[tetraOfPrism[prismT][iTetra][iVertex]];
          }
  }


  void computeMap()
  {
    int PrismVerticesMap[6][6] = {{0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1}, {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0}};

    int tetraOfPrism[2][3][4] = {{{0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3}}, {{0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3}}};

    int TetraFaces[4][3] = {{0 , 1 , 3}, {1 , 2 , 3}, {0 , 3 , 2}, {0 , 2 , 1}};

    int PrismFaces[5][4] = {{0 , 1 , 4 , 3}, {1 , 2 , 5 , 4}, {0 , 3 , 5 , 2}, {0 , 2 , 1 , -1}, {3 , 4 , 5, -1}};


    for(int pType=0; pType<2; ++pType){
      std::cout<< "pType: " << pType <<std::endl;
      for(int minIndex = 0; minIndex<6; ++minIndex){
        std::cout<< "mIndex: " << minIndex <<std::endl;
        for(int pFace = 0; pFace<5; ++pFace){
          for(int tFace =0; tFace<4; ++tFace){
            for(int iTetra =0; iTetra<3; ++iTetra){
              int count=0;
              for(int in =0; in<3; ++in){
                int node=PrismVerticesMap[minIndex][tetraOfPrism[pType][iTetra][TetraFaces[tFace][in]]];
                for(int i=0; i<4; ++i)
                  count += (node == PrismFaces[pFace][i]);
              }
              if(count == 3)
                std::cout << pFace << " " << tFace << " " << iTetra << std::endl;
              }
            }
          }
        }
      std::cout<<std::endl;
      }
    std::cout<<std::endl;
    }


  void tetrasFromPrismStructured (int const* prismVertexMpasIds, int const* prismVertexGIds, int tetrasIdsOnPrism[][4])
  {
      int PrismVerticesMap[6][6] = {{0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3}, {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1}, {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0}};

      int tetraOfPrism[2][3][4] = {{{0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3}}, {{0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3}}};

      int minIndex = std::min_element (prismVertexMpasIds, prismVertexMpasIds + 3) - prismVertexMpasIds;

      int v1 (prismVertexMpasIds[PrismVerticesMap[minIndex][1]]);
      int v2 (prismVertexMpasIds[PrismVerticesMap[minIndex][2]]);

      int prismType = v1  > v2;

      for (int iTetra = 0; iTetra < 3; iTetra++)
          for (int iVertex = 0; iVertex < 4; iVertex++)
          {
              tetrasIdsOnPrism[iTetra][iVertex] = prismVertexGIds[tetraOfPrism[prismType][iTetra][iVertex]];
          }

      // return;

      int reorderedPrismLIds[6];

      for (int ii = 0; ii < 6; ii++)
      {
          reorderedPrismLIds[ii] = prismVertexGIds[PrismVerticesMap[minIndex][ii]];
      }

      for (int iTetra = 0; iTetra < 3; iTetra++)
          for (int iVertex = 0; iVertex < 4; iVertex++)
          {
              tetrasIdsOnPrism[iTetra][iVertex] = reorderedPrismLIds[tetraOfPrism[prismType][iTetra][iVertex]];
          }
  }




  void setBdFacesOnPrism (const std::vector<std::vector<std::vector<int> > >& prismStruct, const std::vector<int>& prismFaceIds, std::vector<int>& tetraPos, std::vector<int>& facePos)
  {
    int numTriaFaces = prismFaceIds.size() - 2;
    tetraPos.assign(numTriaFaces,-1);
    facePos.assign(numTriaFaces,-1);


    for (int iTetra (0), k (0); (iTetra < 3 && k < numTriaFaces); iTetra++)
    {
      bool found;
      for (int jFaceLocalId = 0; jFaceLocalId < 4; jFaceLocalId++ )
      {
        found = true;
        for (int ip (0); ip < 3 && found; ip++)
        {
          int localId = prismStruct[iTetra][jFaceLocalId][ip];
          int j = 0;
          found = false;
          while ( (j < prismFaceIds.size()) && !found )
          {
            found = (localId == prismFaceIds[j]);
            j++;
          }
        }
        if (found)
        {
          tetraPos[k] = iTetra;
          facePos[k] = jFaceLocalId;
          k += found;
          break;
        }
      }
    }
  }

