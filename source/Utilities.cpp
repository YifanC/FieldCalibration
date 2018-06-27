//
// Created by matthias on 28.04.17.
//

#include "../include/Utilities.hpp"
#include "../include/ThreeVector.hpp"

std::vector<Laser> ReachedExitPoint(const Laser &LaserSet, float ExitBoundary) {
    /*
     * Function that splits the input data set into two sets. The first set contains the
     * tracks that have a distance between the last reconstructed point and the expected
     * exit point smaller than the specified distance (ExitBoundary). All other LaserTracks
     * will be put in the second set.
     */

    std::vector<Laser> Selection;
    Selection.resize(2);

    for (auto &Track : LaserSet.GetTrackSet()) {
        auto ExitToLast = Track.GetExitPoint() - Track.GetBack();
        auto d = ExitToLast.GetNorm();

        if (d < ExitBoundary) {
            Selection.front().AppendTrack(Track);
        } else {
            Selection.back().AppendTrack(Track);
        }
    }
    return Selection;
}

std::vector<Laser> SplitTrackSet(const Laser &LaserSet, unsigned int Downsample) {
    /*
     * This function separates the input laser data set into multiple smaller data sets,
     * the number of produced set is specified by the Downsample value.
     */
    std::vector<Laser> Sets;
    Sets.resize(Downsample);

    for (auto &Track : LaserSet.GetTrackSet()) {

        auto SourceTrack = Track.GetReco();

        for (unsigned long offset = 0; offset < Downsample; offset++) {
            std::vector<ThreeVector<float>> SampledRecoTrack;

            for (unsigned long idx = offset; idx < SourceTrack.size(); idx += Downsample) {
                SampledRecoTrack.push_back(SourceTrack[idx]);
            }
            LaserTrack SampledTrack(Track.GetEntryPoint(), Track.GetExitPoint(), SampledRecoTrack);
            Sets[offset].AppendTrack(SampledTrack);
            SampledRecoTrack.clear();

        }
    }
    return Sets;
}

std::vector<Laser> SplitTrackSetDisp(const Laser &LaserSet, unsigned int Downsample) {
    /*
     * This function separates the laser track set into multiple smaller data sets,
     * the number of produced set is specified by the Downsample value.
     * The laser set can be corrected or raw reco track sets.
     */
    std::vector<Laser> Sets;
    Sets.resize(Downsample);

    for (auto &Track : LaserSet.GetTrackSet()) {

        auto SourceTrack = Track.GetReco();
        auto SourceDisp = Track.GetTrackDisp();

        for (unsigned long offset = 0; offset < Downsample; offset++) {
            std::vector<ThreeVector<float>> TrackSample;
            std::vector<ThreeVector<float>> SampleDisp;

            for (unsigned long idx = offset; idx < SourceTrack.size(); idx += Downsample) {
                TrackSample.push_back(SourceTrack[idx]);
                SampleDisp.push_back(SourceDisp[idx]);
            }
            LaserTrack SampledTrack(TrackSample, SampleDisp);
            Sets[offset].AppendTrack(SampledTrack);

            TrackSample.clear();
            SampleDisp.clear();

        }
    }
    return Sets;
}

std::vector<Laser> InterlacedIterTrackSamples(const Laser &LaserSet) {
    // This function separates the input laser data set into 2 data subsamples for iterated correction.
    std::vector<Laser> Sets;
    Sets.resize(2);

    unsigned int number = 0;

    for (auto &Track : LaserSet.GetTrackSet()) {

        if(number%2 == 0){
            Sets[0].AppendTrack(Track);
        }
        if(number%2 == 1){
            Sets[1].AppendTrack(Track);
        }

        number++;
    }
    return Sets;
}


// Merge two "Laser"s into one "Laser"
// the merged Laser will be independent to LaserA and LaserB
Laser MergeLaser(const Laser &LaserA, const Laser &LaserB) {

    // Initialize output of type Laser
    Laser MergedLaser;

    MergedLaser.Mergewith(LaserA);
    MergedLaser.Mergewith(LaserB);

    return MergedLaser;
}

// Acknowledge the zero distortion at anode into mesh by forming the information into "LaserTrack"
LaserTrack Anode(TPCVolumeHandler &TPCVolume){

    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    int AnodeSize = Resolution[1]*Resolution[2];
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    std::vector<ThreeVector<float>> AnodePoints;
    std::vector<ThreeVector<float>> AnodeDisp(AnodeSize,ThreeVector<float>(0.,0.,0.));

    for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
        for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {

            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
            ThreeVector<float> grid = {0., ybin * Unit[1] + TPCVolume.GetDetectorOffset()[1], zbin * Unit[2] + TPCVolume.GetDetectorOffset()[2]};
            AnodePoints.push_back(grid);

        }
    }

    return LaserTrack(AnodePoints,AnodeDisp);
}

//TODO
// Acknowledge the zero distortion at anode into mesh by forming the information into "LaserTrack"
LaserTrack BoundaryCondition(TPCVolumeHandler &TPCVolume){

    ThreeVector<unsigned long> Resolution = TPCVolume.GetDetectorResolution();
    int AnodeSize = Resolution[1]*Resolution[2];
    ThreeVector<float> Unit = {TPCVolume.GetDetectorSize()[0] / (Resolution[0] - 1),
                               TPCVolume.GetDetectorSize()[1] / (Resolution[1] - 1),
                               TPCVolume.GetDetectorSize()[2] / (Resolution[2] - 1)};

    std::vector<ThreeVector<float>> AnodePoints;
    std::vector<ThreeVector<float>> AnodeDisp(AnodeSize,ThreeVector<float>(0.,0.,0.));

    for (unsigned ybin = 0; ybin < Resolution[1]; ybin++) {
        for (unsigned zbin = 0; zbin < Resolution[2]; zbin++) {

            //Push back the location (x,y,z coord) of Anode Points. Anode sits at x=0
            ThreeVector<float> grid = {0., ybin * Unit[1] + TPCVolume.GetDetectorOffset()[1], zbin * Unit[2] + TPCVolume.GetDetectorOffset()[2]};
            AnodePoints.push_back(grid);

        }
    }

    return LaserTrack(AnodePoints,AnodeDisp);
}