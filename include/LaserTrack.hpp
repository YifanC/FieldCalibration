#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <tuple>
#include <array>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <string>
#include <iomanip>

#include <TFile.h>
#include <TH3.h>

#include "ThreeVector.hpp"
#include "TPCVolumeHandler.hpp"
#include <TVector3.h>

#ifndef LASERTRACK_H
#define LASERTRACK_H

class LaserTrack {
private:
    unsigned NumberOfTracksegments;
    std::array<float, 2> TrackAngles;
    ThreeVector<float> LaserPosition;
    ThreeVector<float> EntryPoint;
    ThreeVector<float> ExitPoint;

//   std::vector<ThreeVector<float>> LaserTrue;
    std::vector<ThreeVector<float>> LaserReco;
    std::vector<ThreeVector<float>> LaserDisplacement;

    void FillTrack();

    void FindBoundaries(const TPCVolumeHandler &);

    // Displacement algorithms
    void DerivativeDisplAlgo();

    void ClosestNStretchDisplAlgo(int Nstep);

    void ClosestPointDisplAlgo(int Nstep);

    void LinearStretchDisplAlgo(int Nstep);

    bool CorrMapFlag;//For ClosestPointDisplAlgo

public:
    LaserTrack();

    // Constructor with reco data input with entry point, exit point, and reco track
    LaserTrack(const TVector3 &InEntryPoint, const TVector3 &InExitPoint, const std::vector<TVector3> &RecoTrack);

    LaserTrack(const ThreeVector<float> &InEntryPoint, const ThreeVector<float> &InExitPoint,
               const std::vector<ThreeVector<float>> &RecoTrack);

    LaserTrack(const std::vector<ThreeVector<float>> &RecoTrack, const std::vector<ThreeVector<float>> &LaserDisp);

    LaserTrack(std::array<float, 2> &, ThreeVector<float> &, const TPCVolumeHandler &);

    LaserTrack(const unsigned int, std::array<float, 2> &, ThreeVector<float> &, const TPCVolumeHandler &);

    // Displacement Algorithm names. Add new algorithm name if new algo is introduced
    enum DisplacementAlgo {
        TrackDerivative,
        ClosestNStretch,
        ClosestPoint,
        LinearStretch
    };

    ThreeVector<float> GetPoyntingVector();

    void DistortTrack(std::string, const TPCVolumeHandler &);

    void CalcDisplacement(const DisplacementAlgo &Algo, int Nstep);

    void AddCorrectionToReco();

    void AddCorrectionToRecoPart(std::vector<ThreeVector<float>> CorrPart);

    void Displacement(LaserTrack LaserTrackEnd);

    void AddToDisplacement(ThreeVector<float> &, unsigned long);

    ThreeVector<float> GetLaserPosition();

    std::array<float, 2> GetAngles();

    unsigned long GetNumberOfSamples() const;

    ThreeVector<float> GetSamplePosition(const unsigned int &) const;

    ThreeVector<float> GetDisplacement(const unsigned int &) const;

    ThreeVector<float> GetEntryPoint();

    ThreeVector<float> GetExitPoint();

    ThreeVector<float> GetFront() const;

    ThreeVector<float> GetBack() const;

    void AppendSample(ThreeVector<float> &SamplePosition, ThreeVector<float> &SampleDisplacement);

    void AppendSample(ThreeVector<float> &);

    void AppendSample(float SamplePos_x, float SamplePos_y, float SamplePos_z, float SampleCorr_x, float SampleCorr_y,
                      float SampleCorr_z);

    void AppendSample(float, float, float);

    std::vector<ThreeVector<float>> GetReco();

    std::vector<ThreeVector<float>> GetTrackDisp();

    static void DistortTracks(std::vector<LaserTrack> &, const std::string &, const TPCVolumeHandler &);
};

#endif