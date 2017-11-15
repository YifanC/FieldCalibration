#include <array>
#include <vector>
#include <iostream>
#include <string>
#include <iterator>

#include <TPolyLine3D.h>
#include <TMultiGraph.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TApplication.h>
// #include <TH3.h>
// #include <TFile.h>

#include "LaserTrack.hpp"
#include "Interpolation3D.hpp"

#ifndef LASER_H
#define LASER_H

class Laser {
private:
    std::vector<LaserTrack> LaserTrackSet;

public:
    Laser();

    Laser(std::vector<LaserTrack>);

    std::vector<LaserTrack>::iterator begin();

    std::vector<LaserTrack>::iterator end();

    void AppendTrack(const LaserTrack &InputTrack);

    void Mergewith(const Laser &LaserB);

//    static Laser Merge(std::vector<Laser> &LaserVec);

    LaserTrack GetTrack(const unsigned long int &);

    LaserTrack GetLastTrack();

    LaserTrack GetFirstTrack();

    void SortTracks();

    unsigned long int GetNumberOfTracks() const;

    std::vector<LaserTrack> GetTrackSet() const;

    void DistortTrackSet(std::string, TPCVolumeHandler &);

    void CalcDisplacement(const LaserTrack::DisplacementAlgo &Algo, int Nstep);

//    void AddCorrectionToReco(bool plus);
    void AddCorrectionToReco();

    void SetDisplacement(Laser LaserReco, bool Corr);

    void InterpolateTrackSet(const std::vector<LaserTrack> &, const Delaunay &);

    void InterpolateTrackSet(const Laser &, const Delaunay &);

    void DrawTrack(const unsigned long int &);
//   void Interpolate

};

#endif