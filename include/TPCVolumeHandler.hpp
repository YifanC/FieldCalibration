#include <array>
#include <vector>
#include <cmath>
#include "ThreeVector.hpp"

#ifndef TPCVOLUMEHANDLER_H
#define TPCVOLUMEHANDLER_H

class TPCVolumeHandler {
public:
    TPCVolumeHandler();

    TPCVolumeHandler(const ThreeVector<float> &, const ThreeVector<float> &, const ThreeVector<unsigned long> &);

    TPCVolumeHandler(const std::array<float, 3> &, const std::array<float, 3> &, const std::array<unsigned long, 3> &);

    TPCVolumeHandler(float[3], float[3], unsigned long[3]);

    std::vector<ThreeVector<float>> GetNormalVectors() const;

    std::vector<ThreeVector<float>> GetNormVectorOffset() const;

    ThreeVector<float> GetDetectorSize() const;

    ThreeVector<float> GetDetectorOffset() const;

    ThreeVector<unsigned long> GetDetectorResolution() const;

    ThreeVector<float> GetMapMinimum() const;

    ThreeVector<float> GetMapMaximum() const;

private:
    ThreeVector<float> DetectorSize;
    ThreeVector<float> DetectorOffset;
    ThreeVector<unsigned long> DetectorResolution;
    ThreeVector<float> MapCoordMinimum;
    ThreeVector<float> MapCoordMaximum;

    void CalcNormal(ThreeVector<float> &, ThreeVector<float> &);

    void CalcMapCoordExtreme();

    std::vector<ThreeVector<float>> NormalVector;
    std::vector<ThreeVector<float>> NormalVectorOffSet;

};

#endif