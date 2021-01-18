#include "gtest/gtest.h"
#include "vtkVoxelizer.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkSTLReader.h"

char *path;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    path = argv[1];
    return RUN_ALL_TESTS();
}

class StlFixture : public testing::Test
{
public:
    vtkPolyData *StlFile;

protected:
    void SetUp() override
    {

        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(path);
        reader->Update();
        this->StlFile = reader.GetOutput();
    }
};

TEST(StlFixture, vtkVoxelizer)
{
    stlData = StlFixture->StlFile;

    vtkSmartPointer<vtkVoxelizer> vox = vtkSmartPointer<vtkVoxelizer>::New();
    vox->SetInputData(stlData);
    vox->SetResolution(0.5);

    vox->Update();
}