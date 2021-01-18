#include "vtkVoxelizer.h"
#include <math.h>

#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkImageData.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"

vtkStandardNewMacro(vtkVoxelizer);

vtkVoxelizer::vtkVoxelizer()
{
    this->Resolution = 0.5;
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
};

int vtkVoxelizer::FillInputPortInformation(int vtkNotUsed(port), vtkInformation *info)
{
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE, "vtkPolyData");
}

void vtkVoxelizer::PrintSelf(ostream &os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "Resolution: " << this->Resolution << "\n";
}

int vtkVoxelizer::RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformation **InputVector, vtkInformation *OutputVector)
{
    vtkInformation *outInfo = OutputVector->GetInformation(0);
    std::array<double, 3> spacing;
    spacing.fill(this->Resolution);

    std::array<int, 3> size;

    vtkInformation *inInfo = InputVector[0]->GetInformatObject(0);
    vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    std::array<double, 6> bounds;

    input->GetBounds(bounds);
    for (int i = 0; i < 3; i++)
    {
        double length = bounds[i + 1] - bounds[i];
        size[i] = ceil(length / this->Resolution);
    }

    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), 0, size[0], 0, size[1], 0, size[2]);

    std::array<double, 3> origin;
    input->GetOrigin(origin);
    outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
    outInfo->Set(vtkDataObject::SPACING(), spacing, 3);

    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, double, 1);

    return 1;
}

int vtkVoxelizer::RequestData(vtkInformation *vtkNotUsed(request),
                              vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    vtkPolyData *input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkImageData *output = vtkImageData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    vtkSmartPointer<vtkImageData> bgImg = vtkSmartPointer<vtkImageData>::New();

    output->SetExtent(outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
    output->AllocateScalars(outInfo);

    bgImg->SetExtent(outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
    bgImg->AllocateScalars(outInfo);
    bgImg->GetPointData()->GetScalars()->Fill(0.0);

    vtkSmartPointer<vtkPolyDataToImageStencil> polySten = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
    polySten->SetInputData(input);
    polySten->SetOutputOrigin(bgImg->GetOrigin());
    polySten->SetOutputSpacing(bgImg->GetSpacing());
    polySten->SetOutputExtent(bgImg->GetExtent());
    polySten->Update();

    vtkSmartPointer<vtkImageStencil> stencil = vtkSmartPointer<vtkImageStencil>::New();
    stencil->SetInputData(bgImg);
    stencil->SetStencilConnection(polySten->GetOutputPort());
    stencil->SetBackgroundValue(1.0);
    stencil->ReverseStencilOff();
    stencil->Update();

    output->ShallowCopy(stencil->GetOutput());
    return 1
}