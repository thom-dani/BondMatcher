
/// \ingroup vtk
/// \class ttkSeparatrixGeometricStability
/// \author Thomas Daniel <thomas.daniel@lip6.fr>
/// \date January 2025
///
/// \brief TTK VTK-filter that wraps the ttk::SeparatrixGeometricStability
/// module.
///
/// \param Input vtkDataSet.
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// \sa ttk::SeparatrixGeometricStability
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkSeparatrixGeometricStabilityModule.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>

// VTK Includes
#include <ttkAlgorithm.h>

#include <SeparatrixGeometricStability.h>

class TTKSEPARATRIXGEOMETRICSTABILITY_EXPORT ttkSeparatrixGeometricStability
  : public ttkAlgorithm,
    protected ttk::SeparatrixGeometricStability {
private:
  std::string OutputArrayName{"AveragedScalarField"};

public:
  vtkSetMacro(OutputArrayName, const std::string &);
  vtkGetMacro(OutputArrayName, std::string);

  static ttkSeparatrixGeometricStability *New();
  vtkTypeMacro(ttkSeparatrixGeometricStability, ttkAlgorithm);

protected:
  ttkSeparatrixGeometricStability();
  ~ttkSeparatrixGeometricStability() override = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  float computeDistance(const double (&coords_1)[3],
                        const double (&coords_2)[3]);

  void computePointId(const int &cellId_1,
                      const int &cellId_2,
                      const int &sourceGlobalId,
                      vtkDataSet *separatrix_1,
                      vtkIdType &sourcePointId);

  int execute3D(vtkDataSet *separatrix_1,
              vtkDataSet *separatrix_2,
              vtkFloatArray *distanceToWall,
              vtkIntArray *closestSeparatrixId_2);

  int execute2D(vtkDataSet *separatrix_1,
              vtkFloatArray *distanceToWall,
              vtkIntArray *closestSeparatrixId_2);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
