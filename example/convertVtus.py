import os
import vtk
import sys

def convert_vtkid_to_int(input_file, target_arrays):
    # Read the .vtu file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(input_file)
    reader.Update()
    
    unstructured_grid = reader.GetOutput()
    cell_data = unstructured_grid.GetCellData()
    point_data = unstructured_grid.GetPointData()
    
    # Convert specified cell data arrays from vtkIdTypeArray to vtkIntArray
    for array_name in target_arrays['cell_data']:
        vtk_array = cell_data.GetArray(array_name)
        if vtk_array and isinstance(vtk_array, vtk.vtkIdTypeArray):
            print(f"Converting {array_name} (Cell Data) from vtkIdTypeArray to vtkIntArray")
            new_array = vtk.vtkIntArray()
            new_array.SetName(array_name)
            new_array.SetNumberOfComponents(vtk_array.GetNumberOfComponents())
            new_array.SetNumberOfTuples(vtk_array.GetNumberOfTuples())
            
            for i in range(vtk_array.GetNumberOfTuples()):
                new_array.SetTuple(i, vtk_array.GetTuple(i))
            
            cell_data.RemoveArray(array_name)
            cell_data.AddArray(new_array)
        else:
            print(f"Skipping {array_name} (Cell Data): Not found or not of type vtkIdTypeArray")
    
    # Convert specified point data arrays from vtkIdTypeArray to vtkIntArray
    for array_name in target_arrays['point_data']:
        vtk_array = point_data.GetArray(array_name)
        if vtk_array and isinstance(vtk_array, vtk.vtkIdTypeArray):
            print(f"Converting {array_name} (Point Data) from vtkIdTypeArray to vtkIntArray")
            new_array = vtk.vtkIntArray()
            new_array.SetName(array_name)
            new_array.SetNumberOfComponents(vtk_array.GetNumberOfComponents())
            new_array.SetNumberOfTuples(vtk_array.GetNumberOfTuples())
            
            for i in range(vtk_array.GetNumberOfTuples()):
                new_array.SetTuple(i, vtk_array.GetTuple(i))
            
            point_data.RemoveArray(array_name)
            point_data.AddArray(new_array)
        else:
            print(f"Skipping {array_name} (Point Data): Not found or not of type vtkIdTypeArray")
    
    # Overwrite the original file with the modified data
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(input_file)
    writer.SetInputData(unstructured_grid)
    writer.Write()
    print(f"Original .vtu file overwritten: {input_file}")

def process_vtu_files_in_directory(input_directory, target_arrays):
    # Walk through all the directories and files
    for root, dirs, files in os.walk(input_directory):
        for file in files:
            if file.endswith(".vtu"):
                input_vtu = os.path.join(root, file)
                print(f"Processing file: {input_vtu}")
                convert_vtkid_to_int(input_vtu, target_arrays)

# Main script execution
if __name__ == "__main__":
    # Use the current directory as the root
    input_directory = os.getcwd()
    
    # Define arrays for cell data and point data
    target_arrays = {
        'cell_data': ["ArrayName_1", "ArrayName_2", "ArrayName_3"],  # Cell data arrays to convert
        'point_data': ["CellId"]  # Point data array to convert
    }
    
    process_vtu_files_in_directory(input_directory, target_arrays)

