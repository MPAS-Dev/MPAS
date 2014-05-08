#include "MPASAdaptor.h"

#include "MPASAdaptorAPIMangling.h"

#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkUnsignedCharArray.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkMPI.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIController.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

#include "vtkTrivialProducer.h"
#include "vtkXMLPUnstructuredGridWriter.h"

#include <float.h>
#include <sstream>
#include <string>

#include "vtkSmartPointer.h"
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

using namespace std;

namespace MPAS
{
  int rank;
  int totalRank;

  // Parameters from MPAS renamed for primal mesh use
  int nPrimalCells;		// Number of primal cells
  int nPrimalVerts;		// Number of primal vertices
  int nPrimalVertsPerCell;	// Maximum number of vertices per primal cell
  int nPrimalGhosts;		// Number of ghost cells
  int* primalGhostCell;		// Ghost cell indices
  int* primalGhostHalo;		// Ghost cell levels
  int totalPrimalCells;
  int* nEdgesOnCell;		// Primal cell number of sides
  int* verticesOnCell;		// Point indices for dual cells
  int* vertexMask;              // Valid values for primal cells

  // Parameters from MPAS renamed for dual mesh use
  int nDualCells;		// Number of dual cells
  int nDualVerts;		// Number of dual vertices
  int nDualVertsPerCell;	// Maximum number of vertices per dual cell
  int nDualGhosts;		// Number of dual ghost cells
  int* dualGhostCell;		// Ghost cell indices
  int* dualGhostHalo;		// Ghost cell levels
  int totalDualCells;
  int* cellsOnVertex;		// Point indices for dual cells
  int* cellMask;                // Valid values for dual cells

  int nVertLevels;		// Number of vertex depth levels

  const int PRIMAL = 0;
  const int DUAL   = 1;

  string rankName[] = {"primalRank", "dualRank"};
  string maskName[] = {"primalMask", "dualMask"};

  map<string, vector<bool> > usedCells;

  vtkCPProcessor* Processor = NULL;
  vtkCPDataDescription* Data = NULL;
}

using namespace MPAS;


//////////////////////////////////////////////////////////////////////////
//
// Create the coprocessing grid one time
//
//////////////////////////////////////////////////////////////////////////

#define grid_is_necessary(dataDesc, name)            \
  (dataDesc->GetIfGridIsNecessary(name "-primal") || \
   dataDesc->GetIfGridIsNecessary(name "-dual"))

extern "C" void coprocessor_create_grid_(
                           int* nCells_,
                           int* maxEdges_,
                           int* nGhostCell_,
                           int* cellGhosts_,
                           int* cellHalos_,

                           int* nVertices_,
                           int* vertexDegree_,
                           int* nGhostVertex_,
                           int* vertexGhosts_,
                           int* vertexHalos_,

                           int* nVertLevels_,

                           double* xCell_,
                           double* yCell_,
                           double* zCell_,

                           double* xVertex_,
                           double* yVertex_,
                           double* zVertex_,

                           double* lonCell_,
                           double* latCell_,

                           double* lonVertex_,
                           double* latVertex_,

                           int* nEdgesOnCell_,
                           int* cellsOnVertex_,
                           int* vertexMask_,

                           int* verticesOnCell_,
                           int* cellMask_)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &totalRank);

  nEdgesOnCell = nEdgesOnCell_;
  cellsOnVertex = cellsOnVertex_;
  vertexMask = vertexMask_;
  verticesOnCell = verticesOnCell_;
  cellMask = cellMask_;
  nVertLevels = *nVertLevels_;

  // Set number of cells in primal and dual mesh
  nPrimalCells = *nCells_;
  nPrimalVerts = *nVertices_;
  nPrimalVertsPerCell = *maxEdges_;
  nPrimalGhosts = *nGhostCell_;
  primalGhostCell = cellGhosts_;
  primalGhostHalo = cellHalos_;
  totalPrimalCells = nPrimalCells * nVertLevels;

  nDualCells = *nVertices_;
  nDualVerts = *nCells_;
  nDualVertsPerCell = *vertexDegree_;
  nDualGhosts = *nGhostVertex_;
  dualGhostCell = vertexGhosts_;
  dualGhostHalo = vertexHalos_;
  totalDualCells = nDualCells * nVertLevels;

  // We set the time data to 0 assuming that all pipelines will
  // want to update with these values.
  Data->SetTimeData(0, 0);
  // Processor->RequestDataDescription() only checks which grids
  // will be needed but doesn't do any substantial computation.
  Processor->RequestDataDescription(Data);
  if (grid_is_necessary(Data, "X_Y_NLAYER"))
    // For X,Y cartesian with all layers create 3D cells in 3D space
    create_xy3D_grids(Data, xCell_, yCell_,
                      xVertex_, yVertex_,
                      -10000.0);
  if (grid_is_necessary(Data, "X_Y_Z_1LAYER"))
    // For X,Y,Z spherical with any one layer create 2D cells in 3D space
    create_xyz2D_grids(Data, xCell_, yCell_, zCell_,
                       xVertex_, yVertex_, zVertex_);
  if (grid_is_necessary(Data, "X_Y_Z_NLAYER"))
    // For X,Y,Z spherical with all layers create 3D cells in 3D space
    create_xyz3D_grids(Data, xCell_, yCell_, zCell_,
                       xVertex_, yVertex_, zVertex_,
                       -50000.0);
  if (grid_is_necessary(Data, "LON_LAT_1LAYER"))
    // For lon/lat cartesian with any one layer create 2D cells in 3D space
    create_lonlat2D_grids(Data, lonCell_, latCell_,
                          lonVertex_, latVertex_);
  if (grid_is_necessary(Data, "LON_LAT_NLAYER"))
    // For lon/lat cartesian with all layers create #D cells in 3D space
    create_lonlat3D_grids(Data, lonCell_, latCell_,
                          lonVertex_, latVertex_,
                          -0.1);
}

#define check_datasets(call)        \
  call("X_Y_NLAYER-primal", 3);     \
  call("X_Y_NLAYER-dual", 3);       \
  call("X_Y_Z_1LAYER-primal", 2);   \
  call("X_Y_Z_1LAYER-dual", 2);     \
  call("X_Y_Z_NLAYER-primal", 3);   \
  call("X_Y_Z_NLAYER-dual", 3);     \
  call("LON_LAT_1LAYER-primal", 2); \
  call("LON_LAT_1LAYER-dual", 2);   \
  call("LON_LAT_NLAYER-primal", 3); \
  call("LON_LAT_NLAYER-dual", 3);

//////////////////////////////////////////////////////////////////////////
//
// Register data for the mesh
//
//////////////////////////////////////////////////////////////////////////

static void register_data(
                      char* fname,
                      int* dim0,
                      int* dim1,
                      double* data,
                      const char *datasetName,
                      int cellDim)
{
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(
      Data->GetInputDescriptionByName (datasetName)->GetGrid ());

  int numLevels = *dim0;
  int numCells = *dim1;

  vector<bool>& useCell = usedCells[datasetName];
  if (useCell.size() != numCells) {
    return;
  }

  vtkStdString name(fname);
  string::size_type pos = name.find(' ');
  vtkStdString varName = name.substr(0, pos);
  vtkFloatArray* arr = vtkFloatArray::New();
  arr->SetName(varName.c_str());
  grid->GetCellData()->AddArray(arr);

  // 2D cells have a component for every layer
  if (cellDim == 2) {
    arr->SetNumberOfComponents(numLevels);
    arr->Allocate(numCells);
    float* value = new float[numLevels];
    for (int j = 0; j < numCells; j++) {
      if (useCell[j] == 1) {
        int findx = j * numLevels;
        for (int lev = 0; lev < numLevels; lev++)
          value[lev] = (float) data[findx + lev];
        arr->InsertNextTuple(value);
      }
    }
    delete [] value;
  }

  // 3D cells have a component for every cell
  else {
    arr->SetNumberOfComponents(1);
    arr->Allocate(numCells * numLevels);
    for (int j = 0; j < numCells; j++) {
      for (int lev = 0; lev < numLevels; lev++) {
        if (useCell[j] == 1) {
          int findx = j * numLevels + lev;
          arr->InsertNextValue((float) data[findx]);
        }
      }
    }
  }
  arr->Delete();
}

extern "C" void coprocessor_register_data(
                                     char* fname,
                                     int* dim0,
                                     int* dim1,
                                     double* data)
{
#define coprocessor_register_dataset(name, dim) \
  do                                            \
  {                                             \
    if (Data->GetIfGridIsNecessary(name))       \
      register_data(fname, dim0, dim1, data,    \
        name, dim);                             \
  } while (0)

  check_datasets(coprocessor_register_dataset)

#undef coprocessor_register_dataset
}

//////////////////////////////////////////////////////////////////////////
//
// Tracer data is 3D where the first dimension is the elements of the group
// and so the index of each variable must be given
// Dim 2 is the number of vertex levels and Dim 3 is the number of cells
// So for each cell, each depth, all variable values
//
//////////////////////////////////////////////////////////////////////////

static void register_tracer_data(
                      int* tindex,
                      char* fname,
                      int* dim0,
                      int* dim1,
                      int* dim2,
                      double* data,
                      const char *datasetName,
                      int cellDim)
{
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(
      Data->GetInputDescriptionByName (datasetName)->GetGrid ());

  int varIndx = *tindex - 1;
  int numTracers = *dim0;
  int numLevels = *dim1;
  int numCells = *dim2;

  int perCell = numLevels * numTracers;
  int perLevel = numTracers;

  vector<bool>& useCell = usedCells[datasetName];
  if (useCell.size() != numCells) {
    return;
  }

  vtkStdString name(fname);
  string::size_type pos = name.find(' ');
  vtkStdString varName = name.substr(0, pos);
  vtkFloatArray* arr = vtkFloatArray::New();
  arr->SetName(varName.c_str());
  grid->GetCellData()->AddArray(arr);

  // 2D cells have a component for every layer
  if (cellDim == 2) {
    arr->SetNumberOfComponents(numLevels);
    arr->Allocate(numCells);
    float* value = new float[numLevels];
    for (int cell = 0; cell < numCells; cell++) {
      if (useCell[cell] == 1) {
        for (int lev = 0; lev < numLevels; lev++) {
          int findx = (cell * perCell) + (lev * perLevel) + varIndx;
          value[lev] = (float) data[findx];
        }
        arr->InsertNextTuple(value);
      }
    }
    delete [] value;
  }

  // 3D cells have a component for every cell
  else {
    arr->SetNumberOfComponents(1);
    arr->Allocate(numCells * numLevels);
    for (int cell = 0; cell < numCells; cell++) {
      for (int lev = 0; lev < numLevels; lev++) {
        if (useCell[cell] == 1) {
          int findx = (cell * perCell) + (lev * perLevel) + varIndx;
          arr->InsertNextValue((float) data[findx]);
        }
      }
    }
  }
  arr->Delete();
}

extern "C" void coprocessor_register_tracer_data(
                                     int* tindex,
                                     char* fname,
                                     int* dim0,
                                     int* dim1,
                                     int* dim2,
                                     double* data)
{
#define coprocessor_register_tracer_dataset(name, dim)            \
  do                                                              \
  {                                                               \
    if (Data->GetIfGridIsNecessary(name))                         \
      register_tracer_data(tindex, fname, dim0, dim1, dim2, data, \
        name, dim);                                               \
  } while (0)

  check_datasets(coprocessor_register_tracer_dataset)

#undef coprocessor_register_tracer_dataset
}

//////////////////////////////////////////////////////////////////////////
//
// Load data into the cartesian mesh
//
//////////////////////////////////////////////////////////////////////////

static void add_data(
                      int* itime,
                      char* fname,
                      int* dim0,
                      int* dim1,
                      double* data,
                      const char *datasetName,
                      int cellDim)
{
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(
      Data->GetInputDescriptionByName (datasetName)->GetGrid ());

  int numLevels = *dim0;
  int numCells = *dim1;

  vector<bool>& useCell = usedCells[datasetName];
  if (useCell.size() != numCells) {
    return;
  }

  vtkStdString name(fname);
  string::size_type pos = name.find(' ');
  vtkStdString varName = name.substr(0, pos);
  vtkFloatArray* arr =  vtkFloatArray::SafeDownCast(
        grid->GetCellData()->GetArray(varName.c_str()));
  float* ptr = arr->GetPointer(0);

  // 2D cells have a component for every layer
  if (cellDim == 2) {
    int cindx = 0;
    for (int j = 0; j < numCells; j++) {
      if (useCell[j] == 1) {
        int findx = j * numLevels;
        for (int lev = 0; lev < numLevels; lev++)
          ptr[cindx++] = (float) data[findx + lev];
      }
    }
  }

  // 3D cells have a component for every cell
  else {
    int cindx = 0;
    for (int j = 0; j < numCells; j++) {
      for (int lev = 0; lev < numLevels; lev++) {
        if (useCell[j] == 1) {
          int findx = j * numLevels + lev;
          ptr[cindx++] = (float) data[findx];
        }
      }
    }
  }
}

extern "C" void coprocessor_add_data(int* itime,
                                     char* fname,
                                     int* dim0,
                                     int* dim1,
                                     double* data)
{
#define coprocessor_add_dataset(name, dim)     \
  do                                           \
  {                                            \
    if (Data->GetIfGridIsNecessary(name))      \
      add_data(itime, fname, dim0, dim1, data, \
        name, dim);                            \
  } while (0)

  check_datasets(coprocessor_add_dataset)

#undef coprocessor_add_dataset
}

//////////////////////////////////////////////////////////////////////////
//
// Load data into the cartesian mesh
//
//////////////////////////////////////////////////////////////////////////

static void add_tracer_data(
                      int* itime,
                      int* tindex,
                      char* fname,
                      int* dim0,
                      int* dim1,
                      int* dim2,
                      double* data,
                      const char *datasetName,
                      int cellDim)
{
  vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(
      Data->GetInputDescriptionByName (datasetName)->GetGrid ());

  int varIndx = *tindex - 1;
  int numTracers = *dim0;
  int numLevels = *dim1;
  int numCells = *dim2;

  int perCell = numLevels * numTracers;
  int perLevel = numTracers;

  vector<bool>& useCell = usedCells[datasetName];
  if (useCell.size() != numCells) {
    return;
  }

  vtkStdString name(fname);
  string::size_type pos = name.find(' ');
  vtkStdString varName = name.substr(0, pos);
  vtkFloatArray* arr =  vtkFloatArray::SafeDownCast(
        grid->GetCellData()->GetArray(varName.c_str()));
  float* ptr = arr->GetPointer(0);

  // 2D cells have a component for every layer
  if (cellDim == 2) {
    int cindx = 0;
    for (int cell = 0; cell < numCells; cell++) {
      if (useCell[cell] == 1) {
        for (int lev = 0; lev < numLevels; lev++) {
          int findx = (cell * perCell) + (lev * perLevel) + varIndx;
          ptr[cindx++] = (float) data[findx];
        }
      }
    }
  }

  // 3D cells have a component for every cell
  else {
    int cindx = 0;
    for (int cell = 0; cell < numCells; cell++) {
      for (int lev = 0; lev < numLevels; lev++) {
        if (useCell[cell] == 1) {
          int findx = (cell * perCell) + (lev * perLevel) + varIndx;
          ptr[cindx++] = (float) data[findx];
        }
      }
    }
  }
}

extern "C" void coprocessor_add_tracer_data(int* itime,
                                            int* tindex,
                                            char* fname,
                                            int* dim0,
                                            int* dim1,
                                            int* dim2,
                                            double* data)
{
#define coprocessor_add_tracer_dataset(name, dim)                   \
  do                                                                \
  {                                                                 \
    if (Data->GetIfGridIsNecessary(name))                           \
      add_tracer_data(itime, tindex, fname, dim0, dim1, dim2, data, \
        name, dim);                                                 \
  } while (0)

  check_datasets(coprocessor_add_tracer_dataset)

#undef coprocessor_add_tracer_dataset
}

extern "C" void mpas_init_paraview_catalyst()
{
  // Initialize the adaptor.
  if(!Processor) {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
    VTK_CREATE(vtkCPPythonScriptPipeline, pipeline);
    pipeline->Initialize("mpas.py");
    Processor->AddPipeline(pipeline);

    Data = vtkCPDataDescription::New();
    Data->AddInput("X_Y_NLAYER-primal");
    Data->AddInput("X_Y_NLAYER-dual");
    Data->AddInput("X_Y_Z_1LAYER-primal");
    Data->AddInput("X_Y_Z_1LAYER-dual");
    Data->AddInput("X_Y_Z_NLAYER-primal");
    Data->AddInput("X_Y_Z_NLAYER-dual");
    Data->AddInput("LON_LAT_1LAYER-primal");
    Data->AddInput("LON_LAT_1LAYER-dual");
    Data->AddInput("LON_LAT_NLAYER-primal");
    Data->AddInput("LON_LAT_NLAYER-dual");
  }
  else {
    vtkGenericWarningMacro("Already initialized Catalyst");
  }
}

// doWork will be 0 if no pipelines need to be run
extern "C" void mpas_check_coprocess(int* step, int* doWork)
{
  int rstep = *step;
  double time = static_cast<double>(rstep);
  Data->SetTimeData(time, rstep);
  *doWork = Processor->RequestDataDescription(Data);
}

extern "C" void mpas_coprocess()
{
  Processor->CoProcess(Data);
}

extern "C" void mpas_finalize_paraview_catalyst()
{
  if(Data) {
    Data->Delete();
    Data = NULL;
  }
  else {
    vtkGenericWarningMacro("Nothing to finalize");
  }
  if(Processor) {
    Processor->Finalize();
    Processor->Delete();
    Processor = NULL;
  }
  else {
    vtkGenericWarningMacro("Nothing to finalize");
  }
}
