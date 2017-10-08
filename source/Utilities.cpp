//
// Created by matthias on 28.04.17.
//

#include "../include/Utilities.hpp"

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