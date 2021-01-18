#ifndef vtkVoxelizer.h
#define vtkVoxelizer .h

#include "vtkImageAlgorithm.h"

class vtkVoxelizer : public vtkImageAlgorithm
{
public:
    vtkTypeMacro(vtkVoxelizer, vtkImageAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent) override;

    static vtkVoxelizer *New();

    vtkSetMacro(Resolution, double);

protected:
    vtkVoxelizer();
    ~vtkVoxelizer() override = default;
    int RequestData(vtkInformation *request, vtkInformationVector **inputVector,
                    vtkInformationVector *outputVector) override;
    int FillInputPortInformation(int port, vtkInformation *info) override;

private:
    vtkVoxelizer(const vtkVoxelizer &) = delete;
    void operator=(const vtkVoxelizer &) = delete;
};

#endif