#include <ttkSeparatrixGeometricStability.h>

#include <vtkInformation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkDataSet.h>
#include <vtkAppendFilter.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSignedCharArray.h>
#include <vtkSmartPointer.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <algorithm>
#include <cmath>

vtkStandardNewMacro(ttkSeparatrixGeometricStability);

ttkSeparatrixGeometricStability::ttkSeparatrixGeometricStability() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

int ttkSeparatrixGeometricStability::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    return 1;
  }
  if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkSeparatrixGeometricStability::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

void ttkSeparatrixGeometricStability::computePointId(const int &cellId_1,
                                                     const int &cellId_2,
                                                     const int &sourceGlobalId,
                                                     vtkDataSet *separatrix_1,
                                                     vtkIdType &sourcePointId) {


  vtkIdList *pointIds = separatrix_1->GetCell(cellId_1)->GetPointIds();
  int id_0 = pointIds->GetId(0);
  int id_1 = pointIds->GetId(1);
  pointIds = separatrix_1->GetCell(cellId_2)->GetPointIds();
  int id_2 = pointIds->GetId(0);
  int id_3 = pointIds->GetId(1);
  ttkSimplexIdTypeArray *cellId = ttkSimplexIdTypeArray::SafeDownCast(
    separatrix_1->GetPointData()->GetArray(ttk::MorseSmaleCellIdName));

  if((int)cellId->GetValue(id_0) == sourceGlobalId)
    sourcePointId = id_0;
  else if((int)cellId->GetValue(id_1) == sourceGlobalId)
    sourcePointId = id_1;
  else if((int)cellId->GetValue(id_2) == sourceGlobalId)
    sourcePointId = id_2;
  else if((int)cellId->GetValue(id_3) == sourceGlobalId)
    sourcePointId = id_3;
}

float ttkSeparatrixGeometricStability::computeDistance(
  const double (&coords_1)[3], const double (&coords_2)[3]) {
  return std::sqrt(std::pow(coords_1[0] - coords_2[0], 2)
                   + pow(coords_1[1] - coords_2[1], 2)
                   + pow(coords_1[2] - coords_2[2], 2));
}

int ttkSeparatrixGeometricStability::execute3D(
  vtkDataSet *separatrix_1,
  vtkDataSet *separatrix_2,
  vtkFloatArray *distanceToSeparatrix_2,
  vtkIntArray *closestSeparatrixId_2) {
    
  vtkCellData *cellData_1 = separatrix_1->GetCellData();
  ttkSimplexIdTypeArray *separatrixIds_1 = ttkSimplexIdTypeArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSeparatrixIdName));
  ttkSimplexIdTypeArray *sourceIds = ttkSimplexIdTypeArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSourceIdName));
  vtkSignedCharArray *separatrixTypeArray_1 = vtkSignedCharArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSeparatrixTypeName));

  vtkPoints *separatrixPoints_1 = separatrix_1->GetPoints();

  int cellId_1 = 0;
  int n_cells = separatrix_1->GetNumberOfCells();
  std::vector<float> distances;
  std::vector<int> argMinIds;
  while(cellId_1 < n_cells) {

    int separatrixId = separatrixIds_1->GetValue(cellId_1);
    int cellId_2 = cellId_1 + 1;
    int nextSeparatrixId = separatrixIds_1->GetValue(cellId_2);
    while(nextSeparatrixId == separatrixId && cellId_2 < n_cells) {
      nextSeparatrixId = separatrixIds_1->GetValue(++cellId_2);
    }
    cellId_2--;

    int separatrixType_1 = (int)separatrixTypeArray_1->GetValue(cellId_1);
    if(separatrixType_1 != 0 && separatrixType_1 != 2) {
      distances.push_back(-1);
      argMinIds.push_back(-1);
      cellId_1 = cellId_2 + 1;
      continue;
    }
    int separatrixType_2 = separatrixType_1 == 2 ? 2 : 1;
    vtkSmartPointer<vtkThreshold> thresholdFilterSeparatrix
      = vtkSmartPointer<vtkThreshold>::New();
    thresholdFilterSeparatrix->SetInputData(separatrix_2);
    thresholdFilterSeparatrix->SetInputArrayToProcess(0, 0, 0,
                                      vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                      ttk::MorseSmaleSeparatrixTypeName);
    thresholdFilterSeparatrix->SetLowerThreshold(separatrixType_2);
    thresholdFilterSeparatrix->SetUpperThreshold(separatrixType_2);
    thresholdFilterSeparatrix->Update();
    vtkDataSet *filteredSeparatrix_2 = thresholdFilterSeparatrix->GetOutput();

    if(!filteredSeparatrix_2->GetPoints()){
      const std::string mess = "2_separatrices of type "
                                + std::to_string(separatrixType_2)
                                + " missing for 1_separatrices of type "
                                + std::to_string(separatrixType_1) + ".";
      this->printErr(mess);
      return 0;
    }
    int sourceGlobalId = sourceIds->GetValue(cellId_1);
    vtkIdType sourcePointId;

    computePointId(
      cellId_1, cellId_2, sourceGlobalId, separatrix_1, sourcePointId);

    double sourceCoords[3];
    separatrixPoints_1->GetPoint(sourcePointId, sourceCoords);

    vtkSmartPointer<vtkThreshold> thresholdOtherSeparatrix = vtkSmartPointer<vtkThreshold>::New();
    thresholdOtherSeparatrix->SetInputData(filteredSeparatrix_2);
    thresholdOtherSeparatrix->SetInputArrayToProcess(0, 0, 0, vtkDataObject::CELL, ttk::MorseSmaleSourceIdName);
    thresholdOtherSeparatrix->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    thresholdOtherSeparatrix->SetInvert(true);
    thresholdOtherSeparatrix->SetLowerThreshold(sourceGlobalId);
    thresholdOtherSeparatrix->SetUpperThreshold(sourceGlobalId);
    thresholdOtherSeparatrix->Update();


    vtkDataSet* otherSeparatrix_2 = thresholdOtherSeparatrix->GetOutput();
    double currentPointCoords[3];
    vtkPoints *otherSeparatrixPoints_2 = otherSeparatrix_2->GetPoints();
    otherSeparatrixPoints_2->GetPoint(0, currentPointCoords);
    int threadNumber = omp_get_num_threads();
    std::vector<float> distanceForEachThread(threadNumber);
    distanceForEachThread[0] = computeDistance(sourceCoords, currentPointCoords);

    std::vector<int> argMinIdForEachThread(threadNumber);
    argMinIdForEachThread[0] = 0;

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber)
    #endif // TTK_ENABLE_OPENMP
    for(int i = 1; i < otherSeparatrix_2->GetNumberOfPoints(); i++) {
        otherSeparatrixPoints_2->GetPoint(i, currentPointCoords);
        float currentDistance
          = computeDistance(sourceCoords, currentPointCoords);
        int threadId = omp_get_thread_num();
        if(currentDistance < distanceForEachThread[threadId]) {
          argMinIdForEachThread[threadId] = i;
          distanceForEachThread[threadId] = currentDistance;
      }
    }
    float distanceMin=distanceForEachThread[0];
    int argMinId = argMinIdForEachThread[0];
    for (int i = 1; i< threadNumber; i++){
      if(distanceForEachThread[i]<distanceMin){
        argMinId = argMinIdForEachThread[i];
        distanceMin = distanceForEachThread[i];
      }
    }
    distances.push_back(distanceMin);
    ttkSimplexIdTypeArray *otherSeparatrixIds_2
      = ttkSimplexIdTypeArray::SafeDownCast(
        otherSeparatrix_2->GetCellData()->GetArray(
          ttk::MorseSmaleSeparatrixIdName));
    bool endLoop = false;
    for(unsigned int i = 0; i < otherSeparatrix_2->GetNumberOfCells(); i++) {
      vtkIdList *pointsIds = otherSeparatrix_2->GetCell(i)->GetPointIds();
      for(int j = 0; j < pointsIds->GetNumberOfIds(); j++) {
        if(pointsIds->GetId(j) == argMinId) {
          argMinIds.push_back(otherSeparatrixIds_2->GetValue(i));
          endLoop = true;
          break;
        }
      }
      if(endLoop)
        break;
    }
    cellId_1 = cellId_2 + 1;
  }
  vtkIdType currentCellId = 0;
  int currentSeparatrixId = separatrixIds_1->GetValue(currentCellId);
  int separatrixCount = 0;
  float newValueDistance = distances[separatrixCount];
  int newValueArgMinSeparatrix = argMinIds[separatrixCount];
  distanceToSeparatrix_2->InsertNextValue(newValueDistance);
  closestSeparatrixId_2->InsertNextValue(newValueArgMinSeparatrix);
  for(int j = 1; j < n_cells; j++) {
    if(separatrixIds_1->GetValue(j) != currentSeparatrixId) {
      currentSeparatrixId = separatrixIds_1->GetValue(j);
      separatrixCount++;
      newValueDistance = distances[separatrixCount];
      newValueArgMinSeparatrix = argMinIds[separatrixCount];
    }
    distanceToSeparatrix_2->InsertNextValue(newValueDistance);
    closestSeparatrixId_2->InsertNextValue(newValueArgMinSeparatrix);
  }
  return 1;
}

int ttkSeparatrixGeometricStability::execute2D(
  vtkDataSet *separatrix_1,
  vtkFloatArray *distanceToSeparatrix_1,
  vtkIntArray *closestSeparatrixId_1) {
    
  vtkCellData *cellData_1 = separatrix_1->GetCellData();
  ttkSimplexIdTypeArray *separatrixIds_1 = ttkSimplexIdTypeArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSeparatrixIdName));
  ttkSimplexIdTypeArray *sourceIds = ttkSimplexIdTypeArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSourceIdName));
  vtkSignedCharArray *separatrixTypeArray_1 = vtkSignedCharArray::SafeDownCast(
    cellData_1->GetArray(ttk::MorseSmaleSeparatrixTypeName));

  vtkPoints *separatrixPoints_1 = separatrix_1->GetPoints();

  int cellId_1 = 0;
  int n_cells = separatrix_1->GetNumberOfCells();
  std::vector<float> distances;
  std::vector<int> argMinIds;
  while(cellId_1 < n_cells) {

    int separatrixId = separatrixIds_1->GetValue(cellId_1);
    int cellId_2 = cellId_1 + 1;
    int nextSeparatrixId = separatrixIds_1->GetValue(cellId_2);
    while(nextSeparatrixId == separatrixId && cellId_2 < n_cells) {
      nextSeparatrixId = separatrixIds_1->GetValue(++cellId_2);
    }
    cellId_2--;

    int separatrixType_1 = (int)separatrixTypeArray_1->GetValue(cellId_1);
    if(separatrixType_1 != 0 && separatrixType_1 != 2) {
      distances.push_back(-1);
      argMinIds.push_back(-1);
      cellId_1 = cellId_2 + 1;
      continue;
    }
    vtkSmartPointer<vtkThreshold> thresholdFilterSeparatrix
      = vtkSmartPointer<vtkThreshold>::New();
    thresholdFilterSeparatrix->SetInputData(separatrix_1);
    thresholdFilterSeparatrix->SetInputArrayToProcess(0, 0, 0,
                                      vtkDataObject::FIELD_ASSOCIATION_CELLS,
                                      ttk::MorseSmaleSeparatrixTypeName);
    thresholdFilterSeparatrix->SetLowerThreshold(1-separatrixType_1);
    thresholdFilterSeparatrix->SetUpperThreshold(1-separatrixType_1);
    thresholdFilterSeparatrix->Update();
    vtkDataSet *filteredSeparatrix_1 = thresholdFilterSeparatrix->GetOutput();
    if(!filteredSeparatrix_1->GetPoints()){
      const std::string mess = "1_separatrices of type "
                                + std::to_string(1-separatrixType_1)
                                + " missing for 1_separatrices of type "
                                + std::to_string(separatrixType_1) + ".";
      this->printErr(mess);
      return 0;
    }
    int sourceGlobalId = sourceIds->GetValue(cellId_1);
    vtkIdType sourcePointId;

    vtkSmartPointer<vtkThreshold> thresholdOtherSeparatrix = vtkSmartPointer<vtkThreshold>::New();
    thresholdOtherSeparatrix->SetInputData(filteredSeparatrix_1);
    thresholdOtherSeparatrix->SetInputArrayToProcess(0, 0, 0, vtkDataObject::CELL, ttk::MorseSmaleSourceIdName);
    thresholdOtherSeparatrix->SetThresholdFunction(vtkThreshold::THRESHOLD_BETWEEN);
    thresholdOtherSeparatrix->SetInvert(true);
    thresholdOtherSeparatrix->SetLowerThreshold(sourceGlobalId);
    thresholdOtherSeparatrix->SetUpperThreshold(sourceGlobalId);
    thresholdOtherSeparatrix->Update();

    vtkDataSet* otherSeparatrix_1 = thresholdOtherSeparatrix->GetOutput();
    computePointId(
      cellId_1, cellId_2, sourceGlobalId, separatrix_1, sourcePointId);

    double sourceCoords[3];
    separatrixPoints_1->GetPoint(sourcePointId, sourceCoords);
    double currentPointCoords[3];
    vtkPoints *filteredSeparatrixPoints_1 = otherSeparatrix_1->GetPoints();
    otherSeparatrix_1->GetPoint(0, currentPointCoords);
    int threadNumber = omp_get_num_threads();

    std::vector<float> distanceForEachThread(threadNumber, FLT_MAX);
    std::vector<int> argMinIdForEachThread(threadNumber);

    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(threadNumber)
    #endif // TTK_ENABLE_OPENMP
    for(int i = 0; i < otherSeparatrix_1->GetNumberOfPoints(); i++) {
        filteredSeparatrixPoints_1->GetPoint(i, currentPointCoords);
        float currentDistance
          = computeDistance(sourceCoords, currentPointCoords);
        int threadId = omp_get_thread_num();
        if(currentDistance < distanceForEachThread[threadId]) {
          argMinIdForEachThread[threadId] = i;
          distanceForEachThread[threadId] = currentDistance;
      }
    }
    float distanceMin=distanceForEachThread[0];
    int argMinId = argMinIdForEachThread[0];
    for (int i = 1; i< threadNumber; i++){
      if(distanceForEachThread[i]<distanceMin){
        argMinId = argMinIdForEachThread[i];
        distanceMin = distanceForEachThread[i];
      }
    }
    distances.push_back(distanceMin);
    ttkSimplexIdTypeArray *filteredSeparatrixIds_1
      = ttkSimplexIdTypeArray::SafeDownCast(
        otherSeparatrix_1->GetCellData()->GetArray(
          ttk::MorseSmaleSeparatrixIdName));
    bool endLoop = false;
    for(unsigned int i = 0; i < otherSeparatrix_1->GetNumberOfCells(); i++) {
      vtkIdList *pointsIds = otherSeparatrix_1->GetCell(i)->GetPointIds();
      for(int j = 0; j < pointsIds->GetNumberOfIds(); j++) {
        if(pointsIds->GetId(j) == argMinId) {
          argMinIds.push_back(filteredSeparatrixIds_1->GetValue(i));
          endLoop = true;
          break;
        }
      }
      if(endLoop)
        break;
    }
    cellId_1 = cellId_2 + 1;
  }
  vtkIdType currentCellId = 0;
  int currentSeparatrixId = separatrixIds_1->GetValue(currentCellId);
  int separatrixCount = 0;
  float newValueDistance = distances[separatrixCount];
  int newValueArgMinSeparatrix = argMinIds[separatrixCount];
  distanceToSeparatrix_1->InsertNextValue(newValueDistance);
  closestSeparatrixId_1->InsertNextValue(newValueArgMinSeparatrix);
  for(int j = 1; j < n_cells; j++) {
    if(separatrixIds_1->GetValue(j) != currentSeparatrixId) {
      currentSeparatrixId = separatrixIds_1->GetValue(j);
      separatrixCount++;
      newValueDistance = distances[separatrixCount];
      newValueArgMinSeparatrix = argMinIds[separatrixCount];
    }
    distanceToSeparatrix_1->InsertNextValue(newValueDistance);
    closestSeparatrixId_1->InsertNextValue(newValueArgMinSeparatrix);
  }
  return 1;
}

int ttkSeparatrixGeometricStability::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  vtkDataSet *separatrix_1 = vtkDataSet::GetData(inputVector[0]);
  vtkDataSet *separatrix_2 = vtkDataSet::GetData(inputVector[1]);
  if(!separatrix_1 && !separatrix_2) {
    vtkMultiBlockDataSet* mbSeparatrix_1 = vtkMultiBlockDataSet::GetData(inputVector[0]);
    vtkMultiBlockDataSet* mbSeparatrix_2 = vtkMultiBlockDataSet::GetData(inputVector[1]);

    if(!mbSeparatrix_1 && !mbSeparatrix_2){
      this->printErr("Invalid inputs.");
      return 0;
    }

    if(mbSeparatrix_2){
      if(mbSeparatrix_1->GetNumberOfBlocks() != mbSeparatrix_2->GetNumberOfBlocks()){
        this->printErr("Both inputs should have the same number of vtkDataSet.");
        return 0;
      }
    }
      int n_block = mbSeparatrix_1->GetNumberOfBlocks();
      vtkMultiBlockDataSet* mbOutputSeparatrix_1 = vtkMultiBlockDataSet::GetData(outputVector, 0);
      mbOutputSeparatrix_1->SetNumberOfBlocks(n_block);
      for (int  i = 0 ; i < n_block; i++){
        vtkDataSet* separatrix_1_i = vtkDataSet::SafeDownCast(mbSeparatrix_1->GetBlock(i));
        vtkDataSet* separatrix_2_i = vtkDataSet::SafeDownCast(mbSeparatrix_2->GetBlock(i));
        vtkNew<vtkFloatArray> distanceToSeparatrix_i;
        distanceToSeparatrix_i->SetNumberOfComponents(1);
        vtkNew<vtkIntArray> closestSeparatrixId_i;
        closestSeparatrixId_i->SetNumberOfComponents(1);
        int status = 0; // this integer checks if the base code returns an error
        if(separatrix_2_i){
          distanceToSeparatrix_i->SetName(ttk::DistanceToSeparatrix2Name);
          closestSeparatrixId_i->SetName(ttk::ClosestSeparatrix2IdName);
          status = this->execute3D(
            separatrix_1_i, separatrix_2_i, distanceToSeparatrix_i, closestSeparatrixId_i);
        }else{
          distanceToSeparatrix_i->SetName(ttk::DistanceToSeparatrix1Name);
          closestSeparatrixId_i->SetName(ttk::ClosestSeparatrix1IdName);
          status = this->execute2D(separatrix_1_i, distanceToSeparatrix_i, closestSeparatrixId_i);
        }
        if(status != 1)
        return 0;

        vtkSmartPointer<vtkDataSet> outputSeparatrix_1_i = vtkSmartPointer<vtkDataSet>::Take(separatrix_1_i->NewInstance());
        outputSeparatrix_1_i->ShallowCopy(separatrix_1_i);

        outputSeparatrix_1_i->GetCellData()->AddArray(distanceToSeparatrix_i);
        outputSeparatrix_1_i->GetCellData()->AddArray(closestSeparatrixId_i);
        mbOutputSeparatrix_1->SetBlock(i , outputSeparatrix_1_i);
      }
  }
  else{
  vtkNew<vtkFloatArray> distanceToSeparatrix;
  distanceToSeparatrix->SetNumberOfComponents(1);

  vtkNew<vtkIntArray> closestSeparatrixId;
  closestSeparatrixId->SetNumberOfComponents(1);
  int status = 0; // this integer checks if the base code returns an error

  if(separatrix_2){
    distanceToSeparatrix->SetName(ttk::DistanceToSeparatrix2Name);
    closestSeparatrixId->SetName(ttk::ClosestSeparatrix2IdName); 
    status = this->execute3D(
    separatrix_1, separatrix_2, distanceToSeparatrix, closestSeparatrixId);

  }else{
    distanceToSeparatrix->SetName(ttk::DistanceToSeparatrix1Name);
    closestSeparatrixId->SetName(ttk::ClosestSeparatrix1IdName); 
    status = this->execute2D(separatrix_1, distanceToSeparatrix, closestSeparatrixId);
  }

  if(status != 1)
    return 0;

  vtkDataSet *outputSeparatrix_1 = vtkDataSet::GetData(outputVector, 0);
  outputSeparatrix_1->ShallowCopy(separatrix_1);

  outputSeparatrix_1->GetCellData()->AddArray(distanceToSeparatrix);
  outputSeparatrix_1->GetCellData()->AddArray(closestSeparatrixId);

  }

  return 1;
}
