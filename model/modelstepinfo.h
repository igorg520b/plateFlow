#ifndef MODELSTEPINFO_H
#define MODELSTEPINFO_H

#include <H5Cpp.h>
#include <vector>

namespace icy { struct StepInfo;}

struct icy::StepInfo
{
    int stepNumber = 0;
    double time = 0;
    double currentStepFactor = 1;
    double tentativeStepFactor = 1;

    void Reset()
    {
        stepNumber = 0;
        time = 0;
        currentStepFactor = 1;
        tentativeStepFactor = 1;
    }

    void HDF5Save(H5::H5File *file);
    static void HDF5Read(H5::H5File *file, std::vector<icy::StepInfo> &vec);
    static void HDF5CreateDS(H5::H5File *file);
    static void HDF5Trim(int count, H5::H5File *file);
};

#endif // MODELSTEPINFO_H
